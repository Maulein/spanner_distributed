#include <boost/graph/use_mpi.hpp>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/distributed/graphviz.hpp>
// Standard Library includes
#include <fstream>
#include <string>

// For choose_min_reducer
#include <boost/graph/distributed/distributed_graph_utility.hpp>

#include "bfs1hop_distributed.hpp"

#include <boost/unordered/unordered_flat_map.hpp>
#include <iostream>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp> 

#include <boost/graph/distributed/delta_stepping_shortest_paths.hpp>
#include <cstdlib>   // For srand() random()
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>

#define RMAT
#define graph500
//#define cahepph
//#define webnotre

#include  "rmat_graph_generator.hpp"
#include <boost/random/linear_congruential.hpp>

#include "analyze.hpp"

#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
    	std::cout << ex.what() << std::endl;
    	abort();
}
#endif

using namespace boost;
using boost::graph::distributed::mpi_process_group;
namespace boost {

//Custom Vertex Properties
struct inMIS_t {
  typedef vertex_property_tag kind;
};

struct frontier_t {
  typedef vertex_property_tag kind;
};
struct active_t {
  typedef vertex_property_tag kind;
};
struct more_t {
  typedef vertex_property_tag kind;
};
struct parent_HCIS_t {
  typedef vertex_property_tag kind;
};
struct propagate_t {
  typedef vertex_property_tag kind;
};

struct cluster_t {
  typedef vertex_property_tag kind;
};

struct parent_t {
  typedef vertex_property_tag kind;
};

struct pred_numeric_t {
  typedef vertex_property_tag kind;
};

struct map_t {
  typedef vertex_property_tag kind;
};

struct out_t {
  typedef vertex_property_tag kind;
};

struct spanner_inter_t {
  typedef edge_property_tag kind;
};

struct spanner_intra_t {
  typedef edge_property_tag kind;
};
struct picked_edge_t {
  typedef edge_property_tag kind;
};

  enum vertex_inMIS_t { vertex_inMIS };
  enum vertex_frontier_t { vertex_frontier };
  enum vertex_active_t { vertex_active };
  enum vertex_more_t { vertex_more };
  enum vertex_parent_HCIS_t { vertex_parent_HCIS };
  enum vertex_propagate_t { vertex_propagate };
  enum vertex_cluster_t { vertex_cluster };
  enum vertex_parent_t { vertex_parent};
  enum vertex_pred_numeric_t { vertex_pred_numeric};
  enum vertex_map_t { vertex_map};
  enum vertex_out_t { vertex_out};
  enum edge_spanner_inter_t { edge_spanner_inter};
  enum edge_spanner_intra_t { edge_spanner_intra};
  enum edge_picked_edge_t { edge_picked_edge};

  BOOST_INSTALL_PROPERTY(vertex, inMIS);
  BOOST_INSTALL_PROPERTY(vertex, frontier);
  BOOST_INSTALL_PROPERTY(vertex, active);
  BOOST_INSTALL_PROPERTY(vertex, more);
  BOOST_INSTALL_PROPERTY(vertex, parent_HCIS);
  BOOST_INSTALL_PROPERTY(vertex, propagate);
  BOOST_INSTALL_PROPERTY(vertex, cluster);
  BOOST_INSTALL_PROPERTY(vertex, parent);
  BOOST_INSTALL_PROPERTY(vertex, pred_numeric);
  BOOST_INSTALL_PROPERTY(vertex, map);
  BOOST_INSTALL_PROPERTY(edge, spanner_inter);
  BOOST_INSTALL_PROPERTY(edge, spanner_intra);
  BOOST_INSTALL_PROPERTY(edge, picked_edge);

  namespace graph { namespace distributed { namespace detail {
      template<typename T>
	struct rank_accumulate_reducer {
 	   	BOOST_STATIC_CONSTANT(bool, non_default_resolver = true);

      		template<typename K>
    		  T operator()(const K&) const { return T(0); }

    	 	 template<typename K>
    		  T operator()(const K&, const T& x, const T& y) const { return x + y; }
  	};


      template<typename T>
	struct deg_reducer {
 	   	BOOST_STATIC_CONSTANT(bool, non_default_resolver = true);

      		template<typename K>
    		  T operator()(const K&) const { return T(0); }
    	 	 template<typename K>
    		  T operator()(const K&, const T& x, const T& y) const { return T(0); }
  	};
  } } } // end namespace detail distributed graph
}
  typedef compressed_sparse_row_graph <directedS, 
					distributedS<mpi_process_group>>::vertex_descriptor vertex_descriptor;

  typedef compressed_sparse_row_graph<directedS, 
				 	property<vertex_out_degree_t, size_t,  		//Degree 
				 		property<vertex_degree_t, size_t,  	//Coverage
					 	property<vertex_root_t,bool,  			//Center or not
						property<vertex_distance_t, std::size_t,	//BFS Level
						property<propagate_t,size_t,
						property<cluster_t,vertex_descriptor,
						property<frontier_t, bool,
						property<vertex_predecessor_t, vertex_descriptor,
					        property<parent_t,vertex_descriptor,
						property<map_t, size_t>>>>> >>>>>, //Numeric map
					// Edge Properties
					property<spanner_intra_t, bool>,
                                      	no_property, distributedS<mpi_process_group> >
    				Graph;

typedef adjacency_list<vecS,
                         distributedS<mpi_process_group, vecS>,
                         undirectedS>::vertex_descriptor vertex_descriptor_sg;
typedef adjacency_list<vecS,
                         distributedS<mpi_process_group, vecS>,
                         undirectedS,
			 // Vertex Properties
			 	property<vertex_out_degree_t, size_t,  		//Degree 
			 	property<vertex_degree_t, size_t,  	//Coverage
			 	property<vertex_root_t,bool,  			//Center or not
				property<inMIS_t, bool, 
				property<propagate_t,size_t,
				property<cluster_t,vertex_descriptor_sg,
				property<frontier_t, bool,
				property<active_t, bool,
				property<more_t, bool,
				property<map_t, size_t>>>>> >>>>>, //Numeric map
			// Edge Properties
				no_property, no_property, listS>
    SubGraph;

    
typedef adjacency_list<listS,
                         distributedS<mpi_process_group, vecS>,
                         undirectedS,
			 // Vertex Properties
			 property<vertex_distance_t, size_t>,
			 //Edge Properties
			 property<edge_weight_t,int>>
   SpannerGraph;
 
typedef boost::graph_traits< Graph >::vertex_iterator Viter;

typedef double time_type;

inline time_type get_time()
{
  return MPI_Wtime();
}

std::string print_time(time_type t)
{
  std::ostringstream out;
  out << std::setiosflags(std::ios::fixed) << std::setprecision(2) << t;
  return out.str();
}

void 
calc_coverage(SubGraph &g, int k, int round){		 
  	mpi_process_group pg;
	property_map<SubGraph, map_t>::type
  		 map = get(map_t(), g);
        // Internal property map to store degree of each vertex
	property_map<SubGraph,  vertex_out_degree_t>::type
  		  deg = get(vertex_out_degree, g);

        // Internal property map to store 2-hop coverage of each vertex
	property_map<SubGraph,  vertex_degree_t>::type
  		 coverage = get(vertex_degree, g);
        // Internal property map to store whether vertex in frontier or not
	property_map<SubGraph,  frontier_t>::type
  		 frontier = get(frontier_t(), g);
	typedef typename property_map< SubGraph, vertex_degree_t >::type DegMap;
    	typedef typename property_traits<DegMap >::value_type deg_type;
  	deg.set_reduce(boost::graph::distributed::detail::deg_reducer<deg_type>());
	
	typedef typename property_map< SubGraph, frontier_t >::type FrontierMap;
    	typedef typename property_traits<FrontierMap >::value_type front_type;
  	frontier.set_reduce(boost::graph::distributed::detail::deg_reducer<front_type>());
	size_t self_loop = 0;
	// Store Degree
	if( round <= 1) {  // Computing first time
       		BGL_FORALL_VERTICES(v, g, SubGraph) {
			size_t degree = out_degree(v,g);
			if(edge(v,v,g).second == true) {
				degree -= 2;    //??
				self_loop++;
			}
			put(deg, v, degree);	
  		}
	}
	
	size_t global_self_loop =  boost::parallel::all_reduce(process_group(g), self_loop, std::plus<size_t>());
	// Add up till k hops to get coverages
	if (k == 2){
       		BGL_FORALL_VERTICES(v, g, SubGraph) {
    			put(coverage, v, get(deg, v));
  		}
	} 
	else if (k == 3){
	   // Store Coverage
  	   BGL_FORALL_VERTICES_T(v, g, SubGraph) {
	 	BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){ 
         		request(deg,target(e,g)); 
         		request(frontier,target(e,g)); 
		}
  	   }
	   synchronize(deg);
	   synchronize(frontier);


  	   BGL_FORALL_VERTICES_T(v, g, SubGraph) {
	     if(round == 0 || get(frontier, v)) {
	 	size_t sum = 0;
	 	BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){ 
				sum += get(deg,target(e,g)); 
		}
		if(sum > 0 && edge(v,v,g).second ==  true)
			sum = sum - 2 * get(deg, v);      

		put(coverage,v,sum);
  	     }
	     else
		put(coverage, v, 0);
	  }	  
	}
}
void find_subgraph(Graph &g, float length, float percentage, std::ofstream& f_log, float extend, size_t& size_sub)
{
	// Internal property map to store visit count of a vetrex (during random walk)
	property_map<Graph,  vertex_degree_t>::type
  		 picked = get(vertex_degree, g);

        // Internal property map to store whether vertex is randomly selected seed or not
	property_map<Graph,  map_t>::type
  		 map = get(map_t(), g);
	

        // Internal property map to store whether vertex is to be pushed in queue or not
	property_map<Graph,  frontier_t>::type
  		 put_queue = get(frontier_t(), g);

  	property_map<Graph,  spanner_intra_t>::type
    		picked_edge = get(spanner_intra_t(), g);
		
	property_map<Graph,  vertex_out_degree_t>::type
  		  deg = get(vertex_out_degree, g);
  	
	picked.set_consistency_model(boost::parallel::cm_bidirectional);
  	put_queue.set_consistency_model(boost::parallel::cm_bidirectional);
	typedef typename property_map< Graph, vertex_degree_t >::type DegMap;
    	typedef typename property_traits<DegMap >::value_type deg_type;
  	deg.set_reduce(boost::graph::distributed::detail::deg_reducer<deg_type>());

	auto pro_id = process_id(process_group(g));
    	typedef graph_traits< Graph > GTraits;
    	typedef typename GTraits::vertex_descriptor Vertex;
	
    	typename graph_traits< Graph >::out_edge_iterator edge_it2, edge_it2_end;
    	Vertex v_first = (vertices(g).first[0]);
    	int np = num_processes(process_group(g));
	size_t n = num_vertices(g);
	size_t m = num_edges(g);
  	auto v_invalid = vertex(n+1,g);
	size_t isolated = 0;	
	size_t no_seeds = m  * percentage / 100; // Number of random starting edges for Random Walk
	size_t global_n =  boost::parallel::all_reduce(process_group(g), n, std::plus<size_t>());
	size_t global_m =  boost::parallel::all_reduce(process_group(g), m, std::plus<size_t>());
	size_t global_seeds =  boost::parallel::all_reduce(process_group(g), no_seeds, std::plus<size_t>());
  	typedef property_map<Graph, vertex_owner_t>::const_type OwnerMap;
  	typedef property_map<Graph, vertex_local_t>::const_type LocalMap;

  	OwnerMap owner = get(vertex_owner, g);
  	LocalMap local = get(vertex_local, g);
	BGL_FORALL_VERTICES(v, g, Graph) { 
		put(picked, v, 0);
		put(deg,v,out_degree(v,g));
		put(put_queue,v,false);
    		put(map, v, g.distribution().global(get(owner,v), get(local,v))); 
    		size_t out_deg = out_degree(v,g);
		if(out_deg  == 0)
			isolated++;
	}
	
	BGL_FORALL_VERTICES(v,g,Graph) {
		BGL_FORALL_OUTEDGES(v,e,g,Graph) {
			put(picked_edge,e,false);
			request(deg, target(e,g));
		}
	}
	synchronize(deg);
	size_t random_vertex = 0, random_edge = 0;

	  	srand(time(0));
		typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;
		std::pair<edge_iterator, edge_iterator> ei = edges(g);
		std::vector<edge_iterator> ei_vector(n);
		size_t chunk_size = ceil((float)m / n);
		size_t chunk_no = 0;
		size_t edge_id = 0; 
		for(  edge_iterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter){
			if(edge_id % chunk_size == 0) {
				ei_vector[chunk_no] = edge_iter;
				chunk_no++ ;	
			}
			edge_id++;
		}
		for (size_t i = 0; i < no_seeds;i++ ) {

			random_edge = (int) (rand() % m);
			size_t specific_chunk = (int)(random_edge/chunk_size);	
			size_t j = specific_chunk * chunk_size;
			Vertex random_v;
			Vertex random_u;
			for(  edge_iterator edge_iter = ei_vector[specific_chunk]; edge_iter != ei.second; ++edge_iter)
			{
				if(j==random_edge) {
					random_u = target(*edge_iter, g);
					random_v = source(*edge_iter, g);
					put(picked_edge, *edge_iter, true);
					put(picked, random_v, 1);
					put(picked, random_u, 1);
				if(extend == 0) {		
					if(out_degree(random_v,g) >= get(deg,random_u))
							put(put_queue,random_v,true);
					else
							put(put_queue,random_u,true);
				}
				else if(extend == 1) {		
					if(out_degree(random_v,g) < get(deg,random_u))
							put(put_queue,random_v,true);
					else
							put(put_queue,random_u,true);
				}
				else if(extend == 2) {		
					int choose = rand()%2;
					if(choose == 0)
							put(put_queue,random_v,true);
					else
							put(put_queue,random_u,true);
				}
				break;
				}
				j++;
			}	
		}
	synchronize(put_queue);
      		typedef typename boost::graph::parallel::process_group_type<Graph>::type 
        		process_group_type;
      		typedef typename property_map<Graph, vertex_owner_t>
        		::const_type vertex_owner_map;
		typedef boost::graph::distributed::distributed_queue<
        	        process_group_type, vertex_owner_map, queue<Vertex> > queue_t;
      		queue_t Q(process_group(g),
                	get(vertex_owner, g));

      		queue_t Q1(process_group(g),
                	get(vertex_owner, g));
		if(length > 0) {
		BGL_FORALL_VERTICES(v, g, Graph){
			if(get(put_queue, v) == true) { 
				Q.push(v);
			}
		}
		}

	   for (size_t step = 0; step < length; step++) {
	     if((step%2) == 0) {
    		while (!(Q.empty())) {
        		Vertex u = Q.top();    //Pop from current frontier
			Q.pop();
			//Randomly choose neighbour
			Vertex random_neigh;
            		boost::tie(edge_it2, edge_it2_end) = out_edges(u, g);
			size_t degree = out_degree(u, g);
			if(degree != 0) {
				do {	
					random_edge = rand() % degree;
					random_neigh = target(*(edge_it2 + (int)random_edge), g);
					if(degree <= 3 && random_neigh == u) /////??
						break;
				}  while(random_neigh == u);   // Dont pick self loop
				if(random_neigh != u) {	
					Q1.push(random_neigh);
					put(picked_edge, *(edge_it2 + (int)random_edge) , true);
					put(picked, random_neigh, 1);
					}
			}
		}
	      }	     
	    else {
    		while (!(Q1.empty())) {
        		Vertex u = Q1.top();    //Pop from current frontier
			Q1.pop();
			//Randomly choose neighbour
			Vertex random_neigh;
            		boost::tie(edge_it2, edge_it2_end) = out_edges(u, g);
			size_t degree = out_degree(u, g);
			if(degree != 0) {
				do {	
					random_edge = rand() % degree;
					random_neigh = target(*(edge_it2 + (int)random_edge), g);
					if(degree <= 3 && random_neigh == u) /////??
						break;
				}  while(random_neigh == u);   // Dont pick self loop
				if(random_neigh != u) {	
					Q.push(random_neigh);
					put(picked_edge, *(edge_it2 + (int)random_edge) , true);
					put(picked, random_neigh, 1);
					}
			}
		}
	      }         
	}
	size_t count = 0;
	size_t count_edge = 0;
		BGL_FORALL_VERTICES(v, g, Graph){
			if(get(picked, v) == 1) { 
				count++;
			}
			BGL_FORALL_OUTEDGES(v,e,g,Graph){
				if(get(picked_edge, e) == true) { 
					count_edge++;
				}
			}
		}
	size_t global_count =  boost::parallel::all_reduce(process_group(g), count, std::plus<size_t>());
	size_t global_count_edge =  boost::parallel::all_reduce(process_group(g), count_edge, std::plus<size_t>());
	if(pro_id == 0)
		f_log<<global_seeds<<","<<global_count<<","<<global_count_edge<<",";
	size_sub = count;
	put_queue.clear();
}

void create_subgraph(Graph &g, SubGraph &sub_g, int k) {
	auto pro_id = process_id(process_group(g));
 	typedef graph_traits< Graph > GTraits;
    	typedef typename GTraits::vertex_descriptor Vertex;

 	typedef graph_traits< SubGraph > SubGTraits;
    	typedef typename SubGTraits::vertex_descriptor SubVertex;
  	Vertex v_first_g= (vertices(g).first[0]);
  	SubVertex v_first = (vertices(sub_g).first[0]);
	property_map<Graph,  vertex_degree_t>::type
  		 picked = get(vertex_degree, g);
	
  	property_map<Graph,  spanner_intra_t>::type
    		picked_edge = get(spanner_intra_t(), g);
	
	property_map<Graph,  map_t>::type
  		 map_g = get(map_t(), g);
	
	property_map<Graph,  vertex_out_degree_t>::type
  		  deg_g = get(vertex_out_degree, g);

	property_map<SubGraph,  map_t>::type
  		 map = get(map_t(), sub_g);

	property_map<SubGraph,  vertex_out_degree_t>::type
  		  deg = get(vertex_out_degree, sub_g);

	typedef typename property_map< Graph, frontier_t >::type MapMap;
    	typedef typename property_traits<MapMap >::value_type map_type;
  	map_g.set_reduce(boost::graph::distributed::detail::deg_reducer<map_type>());

  	typedef property_map<Graph, vertex_owner_t>::const_type OwnerMap_g;
	typedef property_map<Graph, vertex_local_t>::const_type LocalMap_g;

  	OwnerMap_g owner_g = get(vertex_owner, g);
  	LocalMap_g local_g = get(vertex_local, g);
  	
	size_t n = num_vertices(sub_g);
	size_t global_n =  boost::parallel::all_reduce(process_group(sub_g), n, std::plus<size_t>());
    	int np = num_processes(process_group(sub_g));
	BGL_FORALL_VERTICES(v, sub_g, SubGraph) { 
    		put(map, v, sub_g.distribution().global(owner(v), local(v)));
	}
	BGL_FORALL_VERTICES(v,g,Graph) {
		BGL_FORALL_OUTEDGES(v,e,g,Graph) {
		   if(get(picked_edge, e) ==  true)	{
			Vertex src_g = source(e,g);
			Vertex target_g = target(e,g);
			size_t src_g_num = get(map_g, src_g);
			size_t target_g_num = g.distribution().global(get(owner_g,target_g), get(local_g,target_g));
			SubVertex src = vertex(src_g_num, sub_g);
			SubVertex target =  vertex( target_g_num, sub_g);
			add_edge(src,target,sub_g);
		  }   
		}	
	}
	size_t m = num_edges(sub_g);
	size_t global_m =  boost::parallel::all_reduce(process_group(sub_g), m, std::plus<size_t>());
  	if (process_id(process_group(g)) == 0) { 
		std::cout<<process_id(process_group(g))<<"Process Done\n";
		std::cout<<pro_id<<"Process: Number of vertices "<<n<<" and edges "<<m<<"\n";
	}
	picked.clear();
	picked_edge.clear();
}
void map_back(Graph &g, SubGraph &sub_g) {

	auto pro_id = process_id(process_group(g));
	// Internal property map to store whether vertex is center or not
	property_map<Graph,  vertex_root_t>::type
  		 root_g = get(vertex_root, g);
        
	property_map<Graph,  map_t>::type
  		 map_g = get(map_t(), g);
	
	// Internal property map to store whether vertex is center or not
	property_map<SubGraph,  vertex_root_t>::type
  		 root = get(vertex_root, sub_g);
        
	property_map<SubGraph,  map_t>::type
  		 map = get(map_t(), sub_g);
	

  	typedef property_map<Graph, vertex_owner_t>::const_type OwnerMap;
  	typedef property_map<Graph, vertex_local_t>::const_type LocalMap;

  	OwnerMap owner = get(vertex_owner, g);
  	LocalMap local = get(vertex_local, g);

	BGL_FORALL_VERTICES(v, g, Graph) {
		put(root_g, v, false);
    		put(map_g, v, g.distribution().global(get(owner,v), get(local,v))); 
	}
	BGL_FORALL_VERTICES(v, sub_g, SubGraph) {
		if(get(root, v) == true) {
			put(root_g,vertex(get(map, v),g), true);	
		}
	}
}

void 
cleanup(Graph &g) {
	// Internal property map to store degree of a vetrex
	property_map<Graph,  vertex_out_degree_t>::type
 		 deg = get(vertex_out_degree, g);

	property_map<Graph,  vertex_degree_t>::type
 		 coverage = get(vertex_degree, g);
	
	property_map<Graph,  vertex_root_t>::type
  		 root = get(vertex_root, g);
	
	property_map<Graph,  vertex_distance_t>::type
  		 distance = get(vertex_distance, g);

	property_map<Graph,  propagate_t>::type
  		 propagate = get(propagate_t(), g);
	property_map<Graph,  cluster_t>::type
  		 cluster = get(cluster_t(), g);
	property_map<Graph,  frontier_t>::type
  		 frontier = get(frontier_t(), g);
 	property_map<Graph,  vertex_predecessor_t>::type
  		 pred = get(vertex_predecessor, g);  // Deposit	

	property_map<Graph,  parent_t>::type
  		 parent = get(parent_t(), g);
	property_map<Graph,  map_t>::type
  		 map = get(map_t(), g);
  	property_map<Graph,  spanner_intra_t>::type
    		spanner_intra = get(spanner_intra_t(), g);
  	BGL_FORALL_VERTICES(v, g, Graph) {
  		BGL_FORALL_OUTEDGES(v, e, g, Graph) {
 	   		put(spanner_intra, e, false);
  		}
	}
	
  	auto v_invalid = vertex(num_vertices(g)+1,g);
  	BGL_FORALL_VERTICES(v, g, Graph) {
		put(frontier, v, false);
    		put(parent, v, v_invalid);
		put(deg, v, 0);  
		put(coverage, v, 0);  
		put(propagate, v, 0);  
		put(root, v, false);  
    		put(distance, v, (std::numeric_limits<std::size_t>::max)());
    		put(cluster, v, v_invalid);
    		put(pred, v, v_invalid);
  	}
}

void 
Subcleanup(SubGraph &g) {
	property_map<SubGraph,  vertex_out_degree_t>::type
 		 deg = get(vertex_out_degree, g);
	property_map<SubGraph,  vertex_degree_t>::type
 		 coverage = get(vertex_degree, g);
	property_map<SubGraph,  vertex_root_t>::type
  		 root = get(vertex_root, g);
	property_map<SubGraph,  inMIS_t>::type
  		 inMIS = get(inMIS_t(), g);
	property_map<SubGraph,  propagate_t>::type
  		 propagate = get(propagate_t(), g);
	property_map<SubGraph,  cluster_t>::type
  		 cluster = get(cluster_t(), g);
	property_map<SubGraph,  active_t>::type
  		 active = get(active_t(), g);
	property_map<SubGraph,  frontier_t>::type
  		 frontier = get(frontier_t(), g);
	property_map<SubGraph,  more_t>::type
  		 more = get(more_t(), g);
	property_map<SubGraph,  map_t>::type
  		 map = get(map_t(), g);
	
  	auto v_invalid = vertex(num_vertices(g)+1,g);
  	BGL_FORALL_VERTICES(v, g, SubGraph) {
		put(frontier, v, false);
    		put(active, v, false);
		put(inMIS, v, false);
		put(deg, v, 0);  
		put(coverage, v, 0);  
		put(propagate, v, 0);  
		put(root, v, false);  
    		put(more, v, false);
    		put(cluster, v, v_invalid);
  	}
}

void 
HCIS(SubGraph &g, int k,  std::ofstream& f_log, int iterations ) {
	property_map<SubGraph,  vertex_out_degree_t>::type
 		 deg = get(vertex_out_degree, g);
	property_map<SubGraph,  vertex_degree_t>::type
 		 coverage = get(vertex_degree, g);
	property_map<SubGraph,  active_t>::type
  		 active = get(active_t(), g);
	property_map<SubGraph,  inMIS_t>::type
  		 inMIS = get(inMIS_t(), g);
	property_map<SubGraph,  vertex_root_t>::type
  		 root = get(vertex_root, g);
	property_map<SubGraph,  cluster_t>::type
  		 parent = get(cluster_t(), g);
	property_map<SubGraph,  frontier_t>::type
  		 frontier = get(frontier_t(), g);
	property_map<SubGraph,  more_t>::type
  		 more = get(more_t(), g);
	property_map<SubGraph,  map_t>::type
  		 map = get(map_t(), g);
	property_map<SubGraph,  propagate_t>::type
  		 propagate = get(propagate_t(), g);
	
  	typedef property_map<SubGraph, vertex_owner_t>::const_type OwnerMap;
  	typedef property_map<SubGraph, vertex_local_t>::const_type LocalMap;

  	OwnerMap owner = get(vertex_owner, g);
  	LocalMap local = get(vertex_local, g);

  	auto v_invalid = vertex(num_vertices(g)+1,g);
  	BGL_FORALL_VERTICES(v, g, SubGraph) {
		put(frontier, v, true);
    		put(active, v, true);
		put(inMIS, v, false);
    		put(parent, v, v_invalid);
		put(deg, v, 0);  
		put(root, v, false);  
    		put(more, v, false);
    		put(map, v, g.distribution().global(get(owner,v), get(local,v))); 
  	}
	size_t n = num_vertices(g);

	size_t centers_count  = 0;	
	size_t total_centers  = 0;	
  	auto pg =  process_group(g);
	size_t global_n =  boost::parallel::all_reduce(pg, n, std::plus<size_t>());
	frontier.set_max_ghost_cells(0);
	
	int iter = 0 ;
	size_t global_frontier_size =  boost::parallel::all_reduce(pg, n, std::plus<size_t>());
	size_t global_frontier_size_old =  0;
	auto pro_id = process_id(process_group(g));
	iterations = 1;

	while (global_frontier_size > 0 && global_frontier_size_old != global_frontier_size && iter < iterations) {
		iter++;
		global_frontier_size_old = global_frontier_size;
	       	calc_coverage(g, k, iter); // calculate coverage only if its active
  		if (process_id(process_group(g)) == 0) 
			std::cout<<"calc coverage done iteration"<<iter<<"\n";
		size_t useless= 0 ;
		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			if(get(coverage, v) == 0)
			{
				useless++;
				put(active, v, false);
			}
		}	
		size_t global_useless =  boost::parallel::all_reduce(pg, useless, std::plus<size_t>());
		if(pro_id == 0)      
			f_log<<global_useless<<",";	
  		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
	    		BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
				request(active, target(e,g));
	    		}
		}
		synchronize(active);

		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			if (get(active,v))
				put(inMIS, v, true);
			put(parent,v, v_invalid);
			put(propagate , v, 0 );
			put(more,v,false);
			if(get(frontier,v)) {
		 	   BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
			     if (get(active, target(e,g)) == true) {
			   	request(coverage, target(e,g));
  		  		request(map, target(e,g));
			     } 
			   }
			}
		}
		synchronize(coverage);
		synchronize(map);
		//Find out if not in MIS due to immediate neighbour
  		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			if( get(frontier, v)) {
			   size_t max_cov = get(coverage,v);	
		 	   BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
			   	if(  target(e,g) != v &&  get(active, target(e,g))) {
					size_t target_cov =  get(coverage,target(e,g)) ;
					if(get(parent,v) == v_invalid) {  //Set parent the first time
         		 			if (get(coverage,v) < target_cov || ((get(coverage,v) == target_cov) && (get(map,v) > get(map, target(e,g))))){	
							put(inMIS, v, false);
							put(parent, v, target(e,g));
							max_cov = target_cov;
							put(more, v, true);
						}
					}
					else
					{
						if((max_cov < target_cov) || ((max_cov == target_cov )&& (get(map,get(parent,v)) > get(map, target(e,g))))) {
							put(inMIS, v, false);
							put(parent, v, target(e,g));
							max_cov = target_cov;
							put(more, v, true);
						}
			 		}
			   	}
		       	   }
			   put(propagate, v, max_cov); 
			}
		}
	 
		synchronize(g);
		int remove_hops = k - 2;
		while(remove_hops >= 1) {
		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
		 	   BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
			   	request(propagate, target(e,g));
			   	request(parent, target(e,g));
			   	request(more, target(e,g));
			   }
		}
		synchronize (propagate);
		synchronize (more);
		synchronize (parent);

  		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
		 	   BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
			   	request(map, get(parent, target(e,g)));
			   	request(map, get(parent, v));
			   }
		}
		synchronize (map);

		//Find out if not in MIS due to 'more'
  		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			   size_t max_cov = get(propagate,v);	
		 	   BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
			   	if(target(e,g) != v && get(more, target(e,g))) {
					size_t target_cov =  get(propagate,target(e,g)) ;
					if((get(frontier,target(e,g)) == false) && (get(more,target(e,g)) == true) && (get(parent,target(e,g))!=v))
							put(inMIS, v, false);
					    
         		 		if (get(parent, target(e,g)) != v && 
			 		   ((max_cov < target_cov ) || 
					     	((max_cov == target_cov) && (get(map, get(parent,v)) > get(map, get(parent, target(e,g))))))){	
							put(inMIS, v, false);
							put(parent, v, get(parent, target(e,g)));
							max_cov = target_cov;
		 		} 
		       	   }
			   put(propagate, v, max_cov);
			}
		    }
		    remove_hops--;
		    synchronize(g);
		}
		// Mark neighbours of in as out due to immediate meighbour
  		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			if(get(active,v)) {
		 	   BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
			   	request(inMIS, target(e,g));
			   	request(frontier, target(e,g));
			   	request(coverage, target(e,g));
			   	request(parent, target(e,g));
			   }
			}
			put(more, v, false);
		}
	       synchronize(inMIS);
	       synchronize(frontier);
	       synchronize(coverage);
	       synchronize(parent);
  		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			if(get(active,v)) {
	 	   		BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
					auto ngh = target(e,g);
						if( target(e,g) != v &&  get(frontier, ngh) && get(inMIS, ngh)) {
		       					put(active, v, false);
							put(inMIS, v, false);
							put(more, v, true);
					  	}
		  		 }
			}
			else  {
				size_t cov = 0;
	 	   		BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
					auto ngh = target(e,g);
						if( target(e,g) != v &&  get(frontier, ngh) && get(inMIS, ngh) && (cov < get(coverage,ngh))) {
							put(more, v, true);
							cov = get(coverage,target(e,g));
							put(parent , v , target(e,g));
					  	}
		  		 }
			}
		}		
		synchronize(g);
		// Mark neighbours of in as out due to more
		remove_hops = k - 2;
		while(remove_hops >= 1) {
  		    BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			if(get(active,v)) {
		 	   BGL_FORALL_OUTEDGES_T(v,e,g,SubGraph){
				request(more,target(e,g)); 
			   	request(parent, target(e,g));
			   }
			}
		   }
		   synchronize (parent);
		   synchronize (more);
  		   BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			if(get(active,v)) {
		 		bool flag = true;
	 	   		BGL_FORALL_OUTEDGES_T(v, e, g, SubGraph){
					auto ngh = target(e,g);
					if( target(e,g) != v &&  get(parent, ngh) != v && get(more, ngh)) {
		       				put(active,v, false);
						put(inMIS,v, false);
						flag = false;
		  			}
		   		}
			}
		   }		
		   remove_hops--;
		   synchronize(g);
		}
		centers_count  = 0;	
  		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			if(get(frontier,v) && get(inMIS, v)){
				put(root, v, true);
				centers_count++;
			}
			put(frontier, v, false);
  		}
		size_t new_frontier_size = 0;
	  	BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			if(get(active,v) && !get(inMIS, v)) {
				put(frontier, v, true);
				new_frontier_size++;
			}
  		}
  		BGL_FORALL_VERTICES_T(v, g, SubGraph) {
	 	   	BGL_FORALL_OUTEDGES_T(v, e, g, SubGraph){
				auto ngh = target(e,g);
				request(frontier, ngh);
			}
		}
		synchronize(frontier);
		if(iterations > 1 ) {
  		   BGL_FORALL_VERTICES_T(v, g, SubGraph) {
			size_t reduced = 0;
	 	   	BGL_FORALL_OUTEDGES_T(v, e, g, SubGraph){
				auto ngh = target(e,g);
				if(ngh != v && get(frontier, ngh)) 
		       			reduced +=1;	
			} 
			put(deg, v , reduced);
		  }
		}	
		global_frontier_size =  boost::parallel::all_reduce(pg, new_frontier_size, std::plus<size_t>());
		size_t global_centers_count =  boost::parallel::all_reduce(pg, centers_count, std::plus<size_t>());
		if(process_id(process_group(g)) == 0){
			std::cout<<process_id(pg)<<"  global_centers_count = "<<global_centers_count<<"\n";
			total_centers +=global_centers_count;
		}
	}
	if(process_id(process_group(g)) == 0) {
		std::cout<<process_id(pg)<<"  Total Centers Count = "<<total_centers<<" in "<<iter<<" iterations\n";
		f_log<< total_centers<<","<<iter<<",";
	}
// clear all ghost cells
active.clear();
inMIS.clear();
frontier.clear();
more.clear();
parent.clear();
coverage.clear();
propagate.clear();
}

void 
ldd(Graph &g, int k, SpannerGraph &sg, std::ofstream& f_log ){

  // Compute BFS levels from root
  property_map<Graph, vertex_distance_t>::type 
    distance = get(vertex_distance, g);
 //  Internal property map to store predecessor
  property_map<Graph,  vertex_predecessor_t>::type
    pred = get(vertex_predecessor, g);
  // Internal property map to store cluster id
  property_map<Graph,  cluster_t>::type
    cluster = get(cluster_t(), g);
  // Internal property map to store final parent
  property_map<Graph,  parent_t>::type
    parent = get(parent_t(), g);
  // Internal property map to store numeric equivalent of vertex descriptor
  property_map<Graph,  map_t>::type
    map = get(map_t(), g); 
  property_map<Graph,  propagate_t>::type
    pred_numeric = get(propagate_t(), g);
  // Internal property map to store degree of each vertex
  property_map<Graph,  vertex_out_degree_t>::type
   deg = get(vertex_out_degree, g);

  // Internal property map to store whether edge is part of spanner or not (Tree)
  property_map<Graph,  spanner_intra_t>::type
    spanner_intra = get(spanner_intra_t(), g);
  spanner_intra.set_consistency_model(boost::parallel::cm_bidirectional);
  
  BGL_FORALL_VERTICES(v, g, Graph) {
  // Initialize distances to infinity and set reduction operation to 'min'
    put(distance, v, (std::numeric_limits<std::size_t>::max)());
    put(pred_numeric, v, (std::numeric_limits<std::size_t>::max)());
  }
  distance.set_reduce(boost::graph::distributed::choose_min_reducer<std::size_t>());
  distance.set_consistency_model(boost::parallel::cm_bidirectional);

  pred_numeric.set_reduce(boost::graph::distributed::choose_min_reducer<std::size_t>());
 
  size_t max_degree = 0;
  size_t isolated = 0;
  // Initialize parents to invalid
  auto v_invalid = vertex(num_vertices(g)+1,g);
  BGL_FORALL_VERTICES(v, g, Graph) {
    put(pred, v, v_invalid);
    put(cluster, v, v_invalid);
    size_t degree_v = get(deg,v);
    if(degree_v == 0)
	isolated++;
    if(max_degree < degree_v)
	max_degree = degree_v;
  }
  size_t global_isolated = boost::parallel::all_reduce(process_group(g), isolated, std::plus<size_t>());
  if(process_id(process_group(g)) == 0)
  {
	std::cout<<"Isolated vertices "<<global_isolated<<"\n";
  	std::cout<<process_id(process_group(g))<<"Process Maximum Degree "<<max_degree<<"\n"; 
  }
  // Initialize edges to non spanner edges
  BGL_FORALL_VERTICES(v, g, Graph) {
  	BGL_FORALL_OUTEDGES(v, e, g, Graph) {
 	   	put(spanner_intra, e, false);
  	}
  }
  //Start spanner construction
  graph_traits<Graph>::vertex_descriptor start = vertex(k,g);
	time_type gen_start, gen_end;
	gen_start = get_time();
  breadth_first_search
    (g, start,
     visitor(make_bfs_visitor(std::make_pair(put_property(distance, 0, on_initialize_vertex()), //Initialize distance to 0
			     	std::make_pair(initialize_pred_nums(pred_numeric, map, on_initialize_vertex()), //
		     		std::make_pair(record_distances(distance, on_tree_edge()),
				std::make_pair(record_pred_nums(pred_numeric, map, on_tree_edge()), 	
				std::make_pair(record_parents(parent, pred_numeric, spanner_intra, distance, map, on_examine_vertex()),  
		     		record_clusters(cluster, parent, on_examine_vertex())))))))));

  	gen_end = get_time();
  auto pg =  process_group(g);
  auto pid = process_id(pg);
        if(pid == 0){
		f_log<<print_time(gen_end - gen_start) <<",";
		std::cout << "INFO: bfs time " << print_time(gen_end - gen_start) << std::endl;
	}
  size_t n = num_vertices(g);  
  //Add all vetices to graph spanner 'sg'
  // Add Intercluster Edges using unordered flat map
  size_t spanner_edges = 0;
  size_t spanner_inter_edges = 0;
  typedef boost::unordered_flat_map<size_t, size_t, 
	                            boost::hash<size_t>, std::equal_to< size_t>,
				    std::allocator<std::pair<const  size_t, size_t>>> MyMap;
  MyMap hashmap;
  
	gen_start = get_time();
  BGL_FORALL_VERTICES(v, g, Graph){
    BGL_FORALL_OUTEDGES(v, e, g, Graph){
      	vertex_descriptor ngh = target(e, g);
	request(cluster,ngh);
	request(distance,ngh);
	}
  }
  synchronize(cluster);
  synchronize(distance);
  
  	gen_end = get_time();
        if(pid == 0) {
		f_log<<print_time(gen_end - gen_start) <<",";
		std::cout << "INFO: Request cluster and level time " << print_time(gen_end - gen_start) << std::endl;
	}
	gen_start = get_time();
  	typedef property_map<Graph, vertex_owner_t>::const_type OwnerMap;
  	typedef property_map<Graph, vertex_local_t>::const_type LocalMap;

  	OwnerMap owner = get(vertex_owner, g);
  	LocalMap local = get(vertex_local, g);
  BGL_FORALL_VERTICES(v, g, Graph){
    BGL_FORALL_OUTEDGES(v, e, g, Graph){
      vertex_descriptor src = source(e, g);
      vertex_descriptor ngh = target(e, g);
      vertex_descriptor c_src = get(cluster,src);
      vertex_descriptor c_ngh = get(cluster,ngh);
      size_t src_num =  g.distribution().global(get(owner,src), get(local,src));
      size_t ngh_num =  g.distribution().global(get(owner,ngh), get(local,ngh)); 
      size_t c_src_num = g.distribution().global(get(owner,c_src), get(local,c_src));      
      size_t c_ngh_num = g.distribution().global(get(owner,c_ngh), get(local,c_ngh));     
      size_t l_src = get(distance,src);
      size_t l_ngh = get(distance,ngh);
 
      if (c_src != c_ngh){ // Edge between 2 clusters	
	  size_t key_part1 = src_num;     
	  size_t key_part2 = c_ngh_num; 

      	  if(src_num > c_ngh_num) { 
	  	key_part1 = c_ngh_num;     
	  	key_part2 = src_num; 
	  }

          if( (l_src > l_ngh) || ((l_src == l_ngh) && (c_src > c_ngh))) 
	  {
	 	  //Key for Hashing      
		  size_t x = (key_part1 << 32) + key_part2;
		  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  		  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  		  x = x ^ (x >> 31);
	          
		  if( !hashmap.contains(x) ) {
      			hashmap.emplace(x, x);
  			spanner_edges++;
  			spanner_inter_edges++;
		}	
       	  }
    }
   } 
  }
  	gen_end = get_time();
        if(pid == 0) {
		f_log<<print_time(gen_end - gen_start) <<",";
		std::cout << "INFO: Intercluster edges added time " << print_time(gen_end - gen_start) << std::endl;
	}
  BGL_FORALL_VERTICES(v, g, Graph) {
	BGL_FORALL_OUTEDGES( v, e, g, Graph) {
	  request(spanner_intra, e);
	}
  }
		  
  synchronize(spanner_intra); 
	gen_start = get_time();
  BGL_FORALL_VERTICES(v, g, Graph) {
      BGL_FORALL_OUTEDGES(v, e, g, Graph) {
        vertex_descriptor src = source(e, g);
        vertex_descriptor ngh = target(e, g);
        size_t src_num =  g.distribution().global(get(owner,src), get(local,src));
        size_t ngh_num =  g.distribution().global(get(owner,ngh), get(local,ngh)); 
        if(get(spanner_intra,e) ==  true)
  		spanner_edges++;

      }
  }
  
  	gen_end = get_time();
        if(pid == 0) {
		std::cout << "INFO: count Intracluster edges time " << print_time(gen_end - gen_start) << std::endl;
		f_log<<print_time(gen_end - gen_start) <<","; }
  auto spg =  process_group(sg);
  synchronize(pg);
  //  all reduce
  size_t global_spanner_edges = boost::parallel::all_reduce(pg, spanner_edges, std::plus<size_t>());
  size_t total;
  total =  boost::parallel::all_reduce(pg, num_vertices(g), std::plus<size_t>());
  if (process_id(process_group(g)) == 0) {
  	std::cout<<"Num of Vertices in Original Graph is "<<total<<"\n";
  	f_log<<total<<",";
  }

  size_t M = num_edges(g);
  total =  boost::parallel::all_reduce(pg, M, std::plus<size_t>());
  if (process_id(process_group(g)) == 0) {
	std::cout<<"Num of Edges in Original Graph is "<<total<<"\n";
  	std::cout<<"Spanner Edges from original graph  "<<global_spanner_edges<<"\n";
  	f_log<<total<<",";
  	f_log<<global_spanner_edges<<",";
  }
  if (process_id(process_group(sg)) == 0) {
        std::cout.flush();
  }	
synchronize(pg);
}

// Graph from File
// argv[1] = file name
// argv[2] = k
// argv[3] = percentage

// RMAT Graph
// argv[1] = n
// argv[2] = m
// argv[3] = k
// argv[4] = percentage
int
main(int argc, char* argv[]){
	boost::mpi::environment env(argc, argv);
	mpi_process_group pg;
  	mpi_process_group::process_id_type id = process_id(pg);	
	std::ofstream f_log ;
	if(id == 0) 
		f_log.open ("Spanner_log.txt", std::ios_base::app);
	time_type gen_start, gen_end;
	time_type gen_start_total, gen_end_total;
	time_type gen_start_sub, gen_end_sub;
	gen_start_total = get_time();
	gen_start = get_time();
	int dig = 4;
	float percentage = 1;
	int i = 0;
	int flag = 0;
	float  num = 0;
	if(argc > dig ) {
	while(argv[dig][i]!='\0') {
		 if(argv[dig][i] == '.')
			flag = 1;
		 else {
		   if(!flag)
			num = (num *10) + (argv[dig][i] - '0');
		   else {
			num = num + (float)((argv[dig][i] - '0') /(float)pow(10,flag));
			flag++;
		   }
		 }
		   i++;
	}}
	if(num != 0)
		percentage = num; 

	srand(time(0));
	int n_RMAT = 10; 
	long m = 20;
#ifdef graph500 
	double a = .57;
	double b = .19;
	double c = .19;
#elif defined cahepph
	double a = .42;
	double b = .19;
	double c = .19;
#elif defined webnotre
	double a = .48;
	double b = .20;
	double c = .21;
#endif
	double d = 1-(a+b+c);	

	
	//m and n input
	dig = 1;
	if(argc > 1)
	{ 
		int num = 0;
		int i = 0;
		while(argv[dig][i]!='\0') {
		   num = (num *10) + (argv[dig][i] - '0');
		   i++;
		}
	        n_RMAT = num;
	}
	dig = 2;
	if(argc > 2)
	{ 
		long num = 0;
		int i = 0;
		while(argv[dig][i]!='\0') {
		   num = (num *10) + (argv[dig][i] - '0');
		   i++;
		}
		m = num;
	}
	
	typedef parallel::variant_distribution<mpi_process_group> Distribution;	
    	Distribution distrib = parallel::block(pg, n_RMAT);


  	typedef keep_local_edges<parallel::variant_distribution<mpi_process_group>,
                           mpi_process_group::process_id_type>
    			   EdgeFilter; 

	typedef scalable_rmat_iterator<mpi_process_group, Distribution, rand48, Graph>
    			   RMATIter;

  	if (id == 0) printf("INFO: Generating graph.\n");

  	rand48 gen;
  	Graph g(RMATIter(pg, distrib, gen, n_RMAT, m, a, b, c, d, true),
          	RMATIter(), n_RMAT, pg, distrib);
	synchronize(g);

  	gen_end = get_time();
        if(id == 0) {
		std::cout << "INFO: Graph Gen time: " << print_time(gen_end - gen_start) << std::endl;
		f_log<<print_time(gen_end - gen_start) <<"\n";
	}
	//parameter k input
	int k = 3; //default - FGV5
	int k_arg = 3; 
	if(argc > k_arg)
		k = argv[k_arg][0] - '0';


	size_t n = num_vertices(g);
	size_t global_n =  boost::parallel::all_reduce(process_group(g), n, std::plus<size_t>());

	const int variations_1 = 1;
	float perc[variations_1]= {1}; //{.5, .25, .1}; //{1, 5};
	const int variations_2 = 1;
	float extend_array[variations_2] = {1};//0, 1, 2}; // Choose (0-HighDegree 1-LowDegree 2-arbitrary) end point
	const int variations_3 = 1;
	float length_array[variations_3] = {4};//2,4,8};//,32};//150, 128, 64, 32, 16, 8};
	const int variations_4 = 1;
	float iterations_array[variations_4] = {1};
	for(int times = 0 ; times < variations_1 ; times++) { 
		for(int times3 = 0 ; times3 < variations_3 ; times3++) { 
			for(int times2 = 0 ; times2 < variations_2 ; times2++) { 
				for(int times4 = 0 ; times4 < variations_4 ; times4++) { 
					percentage = perc[times];
					float extend = extend_array[times2];
					float length = length_array[times3];
					int iterations = iterations_array[times4];

					f_log<<percentage<<","<<extend<<","<<length<<","<<iterations<<"\n";
	const int no_of_runs = 3;
	for(int run_id = 0 ; run_id <no_of_runs ; run_id++) { 

  	gen_start = get_time();
	size_t size_subg;
  	gen_start_sub = get_time();
	find_subgraph(g, length, percentage, f_log, extend, size_subg);

	size_t global_size_subg =  boost::parallel::all_reduce(process_group(g), size_subg, std::plus<size_t>());
  	gen_end_sub = get_time();
        if(id == 0) {
		std::cout << "INFO: Finding Subgraph " << print_time(gen_end_sub - gen_start_sub) << std::endl;
		f_log<<print_time(gen_end_sub - gen_start_sub) <<",";
	}
  	gen_start_sub = get_time();
{	
	SubGraph sub_g(global_n);
	create_subgraph(g, sub_g, k);
  	gen_end_sub = get_time();
        if(id == 0) {
		std::cout << "INFO: Creating Subgraph " << print_time(gen_end_sub - gen_start_sub) << std::endl;
		f_log<<print_time(gen_end_sub - gen_start_sub) <<",";
	}
  	gen_start_sub = get_time();
	HCIS(sub_g, k, f_log, iterations);
  	gen_end_sub = get_time();
        if(id == 0) {
		std::cout << "INFO: HCIS on Subgraph " << print_time(gen_end_sub - gen_start_sub) << std::endl;
		f_log<<print_time(gen_end_sub - gen_start_sub) <<",";
	}
  	gen_start_sub = get_time();
	map_back(g, sub_g);
	Subcleanup(sub_g);
} 
 	gen_end_sub = get_time();
        if(id == 0) {
		std::cout << "INFO: Mapping Subgraph to Original Graph" << print_time(gen_end_sub - gen_start_sub) << std::endl;
		f_log<<print_time(gen_end_sub - gen_start_sub) <<",";
	}
  	gen_end = get_time();
        if(id == 0) {
		std::cout << "INFO: Finding centers " << print_time(gen_end - gen_start) << std::endl;
		f_log<<print_time(gen_end - gen_start) <<",";
	}
  	gen_start = get_time();

  	gen_end = get_time();
        if(id == 0){
		std::cout << "INFO: Center Picking time: " << print_time(gen_end - gen_start) << std::endl;
		std::cout<<"\n\n\n    CONSTRUCT GRAPH SPANNER \n\n\n";	
		f_log<<","<<print_time(gen_end - gen_start) <<",";
	}
	SpannerGraph sg;
  	gen_start = get_time();

  	ldd(g, k, sg, f_log);
  	gen_end = get_time();
        if(id == 0) {
		f_log<<print_time(gen_end - gen_start) <<",";
		std::cout << "INFO: Spanner Construct time: " << print_time(gen_end - gen_start) << std::endl;
	}
	gen_end_total = get_time();
        if(id == 0) {
		f_log<<print_time(gen_end_total - gen_start_total) <<"\n";
		std::cout << "INFO: Total Time: " << print_time(gen_end_total - gen_start_total) << std::endl;
	}
	f_log.flush();
	cleanup(g);
	}
      	if(id == 0)
		f_log<<"\n";
	}}}}
	return 0;
}
