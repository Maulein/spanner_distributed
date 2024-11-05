//
//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
//
#ifndef BOOST_GRAPH_BREADTH_FIRST_SEARCH_HPP
#define BOOST_GRAPH_BREADTH_FIRST_SEARCH_HPP

/*
  Breadth First Search Algorithm (Cormen, Leiserson, and Rivest p. 470)
*/
#include <boost/config.hpp>
#include <vector>
#include <boost/pending/queue.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include "visitors.hpp"  // After adding record_clusters()
//#include <boost/graph/visitors.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/overloading.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/two_bit_color_map.hpp>
#include <boost/graph/detail/mpi_include.hpp>
#include <boost/concept/assert.hpp>
#include <boost/mpi/collectives.hpp> // For Broadcast()
#include <boost/range/algorithm.hpp> // For Random Shuffle
#include <cstdlib>   // For srand()
#include <ctime>     // Fir time()

#include <boost/property_map/property_map.hpp>
#include <boost/property_map/parallel/parallel_property_maps.hpp>
//#include <boost/graph/distributed/detail/queue.ipp>
#include <boost/graph/distributed/distributed_graph_utility.hpp>
#include BOOST_GRAPH_MPI_INCLUDE(<boost/graph/distributed/concepts.hpp>)

namespace boost
{

template < class Visitor, class Graph > struct BFSVisitorConcept
{
    void constraints()
    {
        BOOST_CONCEPT_ASSERT((CopyConstructibleConcept< Visitor >));
        vis.initialize_vertex(u, g);
        vis.discover_vertex(u, g);
        vis.examine_vertex(u, g);
        vis.examine_edge(e, g);
        vis.tree_edge(e, g);
        vis.non_tree_edge(e, g);
        vis.gray_target(e, g);
        vis.black_target(e, g);
        vis.finish_vertex(u, g);
    }
    Visitor vis;
    Graph g;
    typename graph_traits< Graph >::vertex_descriptor u;
    typename graph_traits< Graph >::edge_descriptor e;
};

// Set up the vertex names
enum vertex_id_t { one, two, three, four, five, six, seven };
char vertex_names_b[] = { '1', '2', '3', '4', '5', '6', '7' };



// Multiple-source version
template < class IncidenceGraph, class Buffer, class BFSVisitor, class ColorMap,
    class SourceIterator >
void breadth_first_visit(const IncidenceGraph& g, SourceIterator sources_begin,
    SourceIterator sources_end, Buffer& Q,Buffer& Q1, BFSVisitor vis, ColorMap color)
{
    BOOST_CONCEPT_ASSERT((IncidenceGraphConcept< IncidenceGraph >));
    typedef graph_traits< IncidenceGraph > GTraits;
    typedef typename GTraits::vertex_descriptor Vertex;
    BOOST_CONCEPT_ASSERT((BFSVisitorConcept< BFSVisitor, IncidenceGraph >));
    BOOST_CONCEPT_ASSERT((ReadWritePropertyMapConcept< ColorMap, Vertex >));
    typedef typename property_traits< ColorMap >::value_type ColorValue;
    typedef color_traits< ColorValue > Color;
    typename GTraits::out_edge_iterator ei, ei_end;
    

    auto pg =  process_group(g);
    size_t n = num_vertices(g); 
    int np = num_processes(pg);
    auto pro_id = process_id(pg);
  auto owner = get(boost::vertex_owner, g);
  auto local = get(boost::vertex_local, g);
    Vertex v_first = (vertices(g).first[0]);
    int global_index_begin = g.distribution().global(get(owner,v_first), get(local,v_first));
    Vertex v_first_temp = vertex(0, g);
    // No of rounds
    int k = 3; // Bring as input using sources_begin
    int k_find = 0 ; 
    for ( ; *sources_begin != v_first_temp; )
    { 
	    k_find++;
    	    v_first_temp = vertex(k_find, g);
    }
    k = k_find;
    int last_round = k-1;
    int  r = k-1;
 
    // Generate shifts
    std::vector<int> shifts(last_round+2);
    size_t isolated = 0;
    BGL_FORALL_VERTICES_T(v, g, IncidenceGraph) { 
	if(out_degree(v,g) == 0)
		isolated++;
    }

    size_t global_isolated = boost::parallel::all_reduce(pg, isolated, std::plus<size_t>());
    size_t n_global = boost::parallel::all_reduce(pg, n, std::plus<size_t>());
    size_t n_global_reduced = n_global - global_isolated; 

    double p = (1 - (1 / pow(n_global_reduced, (double)1/k))) ; 
    for(size_t i = 2; i < last_round+2; i++)
 	shifts[i] = ceil(n_global_reduced * p * pow((1-p), last_round - (i-1)));
    shifts[0] = 0;
    shifts[1] = ceil(n_global_reduced* pow((1-p), r));
  
    for(size_t i = 1 ; i < last_round+2 ; i++)
   	shifts[i] += shifts[i-1];
    //Permute vertices
    bool permute = true;
    std::vector<size_t> vertex_permute(n);
    for (size_t i = 0; i < n ; i++)
	    vertex_permute[i] = i;
    srand(time(0));
    
    boost::range::random_shuffle(vertex_permute);
    
    // Common seed to decide which process to get center incase no of processors < vertices to be picked as center in a round 
    int seed = 0; //broadcast to oll
    // BFS
    size_t num_to_add, num_to_add_total, num_added = 0;
    size_t index = 0;
    int ldd_round = 0; // LDD baseline round, begin after radius number of HCIS rounds for growing the clusters
    Buffer *current_frontier = &Q;
    Buffer *next_frontier = &Q1;
    size_t num_visited = 0;
    size_t processor_begin = 0;    
    size_t processor_end = 0;
    while(ldd_round <= last_round) 
    {	
	int start = shifts[ldd_round];
    	int end = shifts[ldd_round + 1];
    	num_to_add_total = end - start; 
 	// divide over processors
	if ( num_to_add_total >= np ) { //  Number of Processors is more than the number of vertices to be picked as centers
	     //processors get extra vertex in round robin fashion
	     num_to_add = num_to_add_total/np;
	     size_t extra = num_to_add_total % np;
	     if(extra !=0) {
             processor_end = (processor_begin + extra) % np;
	     if (processor_begin < processor_end) {
     	     	if ((pro_id >= processor_begin) && (pro_id < processor_end ))  
			num_to_add++; 
	     	else   
	 		num_to_add;
	     }
	     else {

     	     	if ((pro_id >= processor_begin) || (pro_id < processor_end ))  
			num_to_add++; 
	     	else   
	 		num_to_add;
	     }
	     processor_begin = processor_end;
	     }
	}
	else { // Number of Processors is less than the number of vertices to be picked as centers
	       // Randomly pick as many processors to add 1 vertex as center
    	     std::vector<size_t> use(np);
	     int i = 0;
             for (; i < num_to_add_total ; i++)
		    use[i] = 1;	       
             for (; i < np ; i++)
		    use[i] = 0;

    	      if (process_id(pg) == 0)
	    	seed = rand();
    	      broadcast (pg, seed, 0);
	      srand(seed);
		boost::range::random_shuffle(use);
  	    	 if (use[process_id(pg)] == 0 ) // dont useprocessor
			num_to_add = 0;
	     	 else
			num_to_add = 1;
	}
        if( ldd_round == last_round )
		num_to_add = n - num_added;
        if(num_visited >= n)
		num_to_add = 0; 
        if (num_to_add > 0) {
	     int i;
	     for(i = 0; ((i < num_to_add) && (num_added + i) < n) ; i++ ) { 
		   if (permute)
          		index =  vertex_permute[num_added + i];
        	   else
          		index =  num_added + i;
                   int global_index = global_index_begin + index;;
		   Vertex v = vertex(global_index,g);
            	   ColorValue v_color = get(color, v);
		   if (v_color == Color::white()) {
        		   put(color, v, Color::gray());
			   vis.initialize_vertex(v,g);  // To set distance to 0
        		   vis.discover_vertex(v, g);
        		   current_frontier->push(v);
		   }
	     }
       	     num_added += i ;
	}
	synchronize(color);
// FGV - Tie Breaking: Set the lowest ID neighbor as parent for each of the vertices in next frontier - visitor.hpp
    	while (!(current_frontier->empty())) {
        	Vertex u = current_frontier->top();    //Pop from current frontier
		current_frontier->pop();
        	vis.examine_vertex(u, g);		
		if(ldd_round != last_round) {      		
        	for (boost::tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
            		Vertex v = target(*ei, g);
          		vis.examine_edge(*ei, g);
            		ColorValue v_color = get(color, v);
            		if (v_color == Color::white()) {
                		vis.discover_vertex(v, g);
                		vis.tree_edge(*ei, g);
            	    		put(color, v, Color::gray());
				next_frontier->push(v);  // Push in next frontier
 		       	}
            		else  {
                		vis.non_tree_edge(*ei, g);
                		if (v_color == Color::gray())
                    			vis.gray_target(*ei, g);
                		else
                    			vis.black_target(*ei, g);
            		}
        	} // end for
		}
        	put(color, u, Color::black());
        	vis.finish_vertex(u, g);
		num_visited++;
		
    	} // end while
	while (!(next_frontier->empty()))
	{
		current_frontier->push(next_frontier->top());
		next_frontier->pop();
	}
	ldd_round++;
    }

} //breadth_first_visit

// Single-source version
template < class IncidenceGraph, class Buffer, class BFSVisitor,
    class ColorMap >
void breadth_first_visit(const IncidenceGraph& g,
    typename graph_traits< IncidenceGraph >::vertex_descriptor s, Buffer& Q,
Buffer& Q1, BFSVisitor vis, ColorMap color)
{
    typename graph_traits< IncidenceGraph >::vertex_descriptor sources[1]
        = { s };
    breadth_first_visit(g, sources, sources + 1, Q, Q1, vis, color);
}

template < class VertexListGraph, class SourceIterator, class Buffer,
    class BFSVisitor, class ColorMap >
void breadth_first_search(const VertexListGraph& g,
    SourceIterator sources_begin, SourceIterator sources_end, Buffer& Q,Buffer& Q1,
    BFSVisitor vis, ColorMap color)
{
    // Initialization
    typedef typename property_traits< ColorMap >::value_type ColorValue;
    typedef color_traits< ColorValue > Color;
    typename boost::graph_traits< VertexListGraph >::vertex_iterator i, i_end;
    for (boost::tie(i, i_end) = vertices(g); i != i_end; ++i)
    {
//        vis.initialize_vertex(*i, g);
        put(color, *i, Color::white());
    }
    breadth_first_visit(g, sources_begin, sources_end, Q, Q1,vis, color);
}

template < class VertexListGraph, class Buffer, class BFSVisitor,
    class ColorMap >
void breadth_first_search(const VertexListGraph& g,
    typename graph_traits< VertexListGraph >::vertex_descriptor s, Buffer& Q,Buffer& Q1,
    BFSVisitor vis, ColorMap color)
{
    typename graph_traits< VertexListGraph >::vertex_descriptor sources[1]
        = { s };
    breadth_first_search(g, sources, sources + 1, Q, Q1, vis, color);
}

namespace graph
{
    struct bfs_visitor_event_not_overridden
    {
    };
}

template < class Visitors = null_visitor > class bfs_visitor
{
public:
    bfs_visitor() {}
    bfs_visitor(Visitors vis) : m_vis(vis) {}

    template < class Vertex, class Graph >
    graph::bfs_visitor_event_not_overridden initialize_vertex(
        Vertex u, Graph& g)
    {
        invoke_visitors(m_vis, u, g, ::boost::on_initialize_vertex());
        return graph::bfs_visitor_event_not_overridden();
    }

    template < class Vertex, class Graph >
    graph::bfs_visitor_event_not_overridden discover_vertex(Vertex u, Graph& g)
    {
        invoke_visitors(m_vis, u, g, ::boost::on_discover_vertex());
        return graph::bfs_visitor_event_not_overridden();
    }

    template < class Vertex, class Graph >
    graph::bfs_visitor_event_not_overridden examine_vertex(Vertex u, Graph& g)
    {
        invoke_visitors(m_vis, u, g, ::boost::on_examine_vertex());
        return graph::bfs_visitor_event_not_overridden();
    }

    template < class Edge, class Graph >
    graph::bfs_visitor_event_not_overridden examine_edge(Edge e, Graph& g)
    {
        invoke_visitors(m_vis, e, g, ::boost::on_examine_edge());
        return graph::bfs_visitor_event_not_overridden();
    }

    template < class Edge, class Graph >
    graph::bfs_visitor_event_not_overridden tree_edge(Edge e, Graph& g)
    {
        invoke_visitors(m_vis, e, g, ::boost::on_tree_edge());
        return graph::bfs_visitor_event_not_overridden();
    }

    template < class Edge, class Graph >
    graph::bfs_visitor_event_not_overridden non_tree_edge(Edge e, Graph& g)
    {
        invoke_visitors(m_vis, e, g, ::boost::on_non_tree_edge());
        return graph::bfs_visitor_event_not_overridden();
    }

    template < class Edge, class Graph >
    graph::bfs_visitor_event_not_overridden gray_target(Edge e, Graph& g)
    {
        invoke_visitors(m_vis, e, g, ::boost::on_gray_target());
        return graph::bfs_visitor_event_not_overridden();
    }

    template < class Edge, class Graph >
    graph::bfs_visitor_event_not_overridden black_target(Edge e, Graph& g)
    {
        invoke_visitors(m_vis, e, g, ::boost::on_black_target());
        return graph::bfs_visitor_event_not_overridden();
    }

    template < class Vertex, class Graph >
    graph::bfs_visitor_event_not_overridden finish_vertex(Vertex u, Graph& g)
    {
        invoke_visitors(m_vis, u, g, ::boost::on_finish_vertex());
        return graph::bfs_visitor_event_not_overridden();
    }

    BOOST_GRAPH_EVENT_STUB(on_initialize_vertex, bfs)
    BOOST_GRAPH_EVENT_STUB(on_discover_vertex, bfs)
    BOOST_GRAPH_EVENT_STUB(on_examine_vertex, bfs)
    BOOST_GRAPH_EVENT_STUB(on_examine_edge, bfs)
    BOOST_GRAPH_EVENT_STUB(on_tree_edge, bfs)
    BOOST_GRAPH_EVENT_STUB(on_non_tree_edge, bfs)
    BOOST_GRAPH_EVENT_STUB(on_gray_target, bfs)
    BOOST_GRAPH_EVENT_STUB(on_black_target, bfs)
    BOOST_GRAPH_EVENT_STUB(on_finish_vertex, bfs)

protected:
    Visitors m_vis;
};
template < class Visitors >
bfs_visitor< Visitors > make_bfs_visitor(Visitors vis)
{
    return bfs_visitor< Visitors >(vis);
}
typedef bfs_visitor<> default_bfs_visitor;

namespace detail
{

    template < class VertexListGraph, class ColorMap, class BFSVisitor, class P,
        class T, class R >
    void bfs_helper(VertexListGraph& g,
        typename graph_traits< VertexListGraph >::vertex_descriptor s,
        ColorMap color, BFSVisitor vis,
        const bgl_named_params< P, T, R >& params, boost::mpl::false_)
    {
        typedef graph_traits< VertexListGraph > Traits;
        // Buffer default
        typedef typename Traits::vertex_descriptor Vertex;
        typedef boost::queue< Vertex > queue_t;
        queue_t Q,Q1;
        breadth_first_search(g, s,
            choose_param(get_param(params, buffer_param_t()), boost::ref(Q), boost::ref(Q1))
                .get(),
            vis, color);
    }

#ifdef BOOST_GRAPH_USE_MPI
    template < class DistributedGraph, class ColorMap, class BFSVisitor,
        class P, class T, class R >
    void bfs_helper(DistributedGraph& g,
        typename graph_traits< DistributedGraph >::vertex_descriptor s,
        ColorMap color, BFSVisitor vis,
        const bgl_named_params< P, T, R >& params, boost::mpl::true_);
#endif // BOOST_GRAPH_USE_MPI

    //-------------------------------------------------------------------------
    // Choose between default color and color parameters. Using
    // function dispatching so that we don't require vertex index if
    // the color default is not being used.

    template < class ColorMap > struct bfs_dispatch
    {
        template < class VertexListGraph, class P, class T, class R >
        static void apply(VertexListGraph& g,
            typename graph_traits< VertexListGraph >::vertex_descriptor s,
            const bgl_named_params< P, T, R >& params, ColorMap color)
        {
            bfs_helper(g, s, color,
                choose_param(get_param(params, graph_visitor),
                    make_bfs_visitor(null_visitor())),
                params,
                boost::mpl::bool_<
                    boost::is_base_and_derived< distributed_graph_tag,
                        typename graph_traits<
                            VertexListGraph >::traversal_category >::value >());
        }
    };

    template <> struct bfs_dispatch< param_not_found >
    {
        template < class VertexListGraph, class P, class T, class R >
        static void apply(VertexListGraph& g,
            typename graph_traits< VertexListGraph >::vertex_descriptor s,
            const bgl_named_params< P, T, R >& params, param_not_found)
        {
            null_visitor null_vis;

            bfs_helper(g, s,
                make_two_bit_color_map(num_vertices(g),
                    choose_const_pmap(
                        get_param(params, vertex_index), g, vertex_index)),
                choose_param(get_param(params, graph_visitor),
                    make_bfs_visitor(null_vis)),
                params,
                boost::mpl::bool_<
                    boost::is_base_and_derived< distributed_graph_tag,
                        typename graph_traits<
                            VertexListGraph >::traversal_category >::value >());
        }
    };

} // namespace detail

#if 1
// Named Parameter Variant
template < class VertexListGraph, class P, class T, class R >
void breadth_first_search(const VertexListGraph& g,
    typename graph_traits< VertexListGraph >::vertex_descriptor s,
    const bgl_named_params< P, T, R >& params)
{
    // The graph is passed by *const* reference so that graph adaptors
    // (temporaries) can be passed into this function. However, the
    // graph is not really const since we may write to property maps
    // of the graph.
    VertexListGraph& ng = const_cast< VertexListGraph& >(g);
    typedef typename get_param_type< vertex_color_t,
        bgl_named_params< P, T, R > >::type C;
    detail::bfs_dispatch< C >::apply(
        ng, s, params, get_param(params, vertex_color));
}
#endif

// This version does not initialize colors, user has to.

template < class IncidenceGraph, class P, class T, class R >
void breadth_first_visit(const IncidenceGraph& g,
    typename graph_traits< IncidenceGraph >::vertex_descriptor s,
    const bgl_named_params< P, T, R >& params)
{
    // The graph is passed by *const* reference so that graph adaptors
    // (temporaries) can be passed into this function. However, the
    // graph is not really const since we may write to property maps
    // of the graph.
    IncidenceGraph& ng = const_cast< IncidenceGraph& >(g);

    typedef graph_traits< IncidenceGraph > Traits;
    // Buffer default
    typedef typename Traits::vertex_descriptor vertex_descriptor;
    typedef boost::queue< vertex_descriptor > queue_t;
    queue_t Q,Q1;

    breadth_first_visit(ng, s,
        choose_param(get_param(params, buffer_param_t()), boost::ref(Q), boost::ref(Q)).get(),
        choose_param(
            get_param(params, graph_visitor), make_bfs_visitor(null_visitor())),
        choose_pmap(get_param(params, vertex_color), ng, vertex_color));
}

namespace graph
{
    namespace detail
    {
        template < typename Graph, typename Source >
        struct breadth_first_search_impl
        {
            typedef void result_type;
            template < typename ArgPack >
            void operator()(
                const Graph& g, const Source& source, const ArgPack& arg_pack)
            {
                using namespace boost::graph::keywords;
                typename boost::graph_traits< Graph >::vertex_descriptor
                    sources[1]
                    = { source };
                boost::queue<
                    typename boost::graph_traits< Graph >::vertex_descriptor >
                    Q;
                boost::breadth_first_search(g, &sources[0], &sources[1],
                    boost::unwrap_ref(arg_pack[_buffer | boost::ref(Q)]),
                    arg_pack[_visitor | make_bfs_visitor(null_visitor())],
                    boost::detail::make_color_map_from_arg_pack(g, arg_pack));
            }
        };
    }
    BOOST_GRAPH_MAKE_FORWARDING_FUNCTION(breadth_first_search, 2, 4)
}

#if 0
  // Named Parameter Variant
  BOOST_GRAPH_MAKE_OLD_STYLE_PARAMETER_FUNCTION(breadth_first_search, 2)
#endif

} // namespace boost

//#include BOOST_GRAPH_MPI_INCLUDE(<boost/graph/distributed/breadth_first_search.hpp>)
#include BOOST_GRAPH_MPI_INCLUDE("bfs1hop_distributed.hpp")
#endif // BOOST_GRAPH_BREADTH_FIRST_SEARCH_HPP
