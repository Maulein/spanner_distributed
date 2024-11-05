#include <boost/graph/use_mpi.hpp>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/distributed/compressed_sparse_row_graph.hpp>
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
#include <cstdlib>   // For srand()

#define graph500
//#define cahepph
//#define webnotre

#include  "rmat_graph_generator.hpp"
#include <boost/random/linear_congruential.hpp>

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

struct spanner_intra_t {
  typedef edge_property_tag kind;
};

  enum vertex_cluster_t { vertex_cluster };
  enum vertex_parent_t { vertex_parent};
  enum vertex_pred_numeric_t { vertex_pred_numeric};
  enum vertex_map_t { vertex_map};
  enum edge_spanner_intra_t { edge_spanner_intra};

  BOOST_INSTALL_PROPERTY(vertex, cluster);
  BOOST_INSTALL_PROPERTY(vertex, parent);
  BOOST_INSTALL_PROPERTY(vertex, pred_numeric);
  BOOST_INSTALL_PROPERTY(vertex, map);
  BOOST_INSTALL_PROPERTY(edge, spanner_intra);
}
  typedef compressed_sparse_row_graph <directedS, 
					distributedS<mpi_process_group>>::vertex_descriptor vertex_descriptor;

  typedef compressed_sparse_row_graph<directedS, 
				 	property<vertex_out_degree_t, size_t,  		//Degree 
				 		property<vertex_degree_t, size_t,  	//Coverage
						property<vertex_distance_t, std::size_t,	//BFS Level
						property<cluster_t,vertex_descriptor,
						property<vertex_predecessor_t, vertex_descriptor,
					        property<parent_t,vertex_descriptor,
					        property<pred_numeric_t,size_t,
						property<map_t, size_t>>>>> >>>, //Numeric map
					// Edge Properties
					property<spanner_intra_t, bool>,
                                      	no_property, distributedS<mpi_process_group> >
    				Graph;
    
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
coverage(Graph &g){		 
  	mpi_process_group pg;
  	synchronize(pg);
	
        // Internal property map to store degree of each vertex
	property_map<Graph,  vertex_out_degree_t>::type
  		  deg = get(vertex_out_degree, g);

        // Internal property map to store 2-hop coverage of each vertex
	property_map<Graph,  vertex_degree_t>::type
  		 coverage = get(vertex_degree, g);

	// Store Degree
       	BGL_FORALL_VERTICES(v, g, Graph) {
    		put(deg, v, out_degree(v,g));
  	}
	
	// Store Coverage
  	BGL_FORALL_VERTICES_T(v, g, Graph) {
	 	size_t sum = 0;
	 	BGL_FORALL_OUTEDGES_T(v,e,g,Graph){ 
         		sum += get(deg,target(e,g)); 
		}
		put(coverage,v,sum);
  	}
}

void 
ldd(Graph &g, int k, SpannerGraph &sg,std::ofstream& f_log ){

	time_type gen_start, gen_end;
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

  // Internal property map to store numeric equivalent of predeccesor
  property_map<Graph,  pred_numeric_t>::type
    pred_numeric = get(pred_numeric_t(), g);
  
  // Internal property map to store whether edge is part of spanner or not (Tree)
  property_map<Graph,  spanner_intra_t>::type
    spanner_intra = get(spanner_intra_t(), g);

  spanner_intra.set_consistency_model(boost::parallel::cm_bidirectional);
  
  	typedef property_map<Graph, vertex_owner_t>::const_type OwnerMap;
  	typedef property_map<Graph, vertex_local_t>::const_type LocalMap;

  	OwnerMap owner = get(vertex_owner, g);
  	LocalMap local = get(vertex_local, g);
  BGL_FORALL_VERTICES(v, g, Graph) {
  // Initialize distances to infinity and set reduction operation to 'min'
    put(distance, v, (std::numeric_limits<std::size_t>::max)());
  // Initialize distances to infinity and set reduction operation to 'min'
    put(pred_numeric, v, (std::numeric_limits<std::size_t>::max)());
  // Initialize map  
    put(map, v, g.distribution().global(get(owner,v), get(local,v))); 

  }
  distance.set_reduce(boost::graph::distributed::choose_min_reducer<std::size_t>());
  distance.set_consistency_model(boost::parallel::cm_bidirectional);

  pred_numeric.set_reduce(boost::graph::distributed::choose_min_reducer<std::size_t>());
 
  // Initialize parents to invalid
  auto pg = process_group(g);
  size_t n = num_vertices(g);
  size_t global_n = boost::parallel::all_reduce(pg, n, std::plus<size_t>());
  auto v_invalid = vertex(global_n -1,g);
  BGL_FORALL_VERTICES(v, g, Graph) {
    put(pred, v, v_invalid);
    put(cluster, v, v_invalid);
    put(parent, v, v_invalid);
  } 
  // Initialize edges to non spanner edges
  BGL_FORALL_VERTICES(v, g, Graph) {
	  BGL_FORALL_OUTEDGES(v, e, g, Graph) {
 		 put(spanner_intra, e, false);
  	  }
  }
  //Start spanner construction
  graph_traits<Graph>::vertex_descriptor start = vertex(k,g);
	gen_start = get_time();
  breadth_first_search
    (g, start,
     visitor(make_bfs_visitor(std::make_pair(put_property(distance, 0, on_initialize_vertex()), //Initialize distance to 0
			     	std::make_pair(initialize_pred_nums(pred_numeric, map, on_initialize_vertex()), //
		     		std::make_pair(record_distances(distance, on_tree_edge()),
				std::make_pair(record_pred_nums(pred_numeric, map, on_tree_edge()), 	
				std::make_pair(record_parents(parent, pred_numeric, spanner_intra, distance, map, on_examine_vertex()),  
		     		record_clusters(cluster, parent, on_examine_vertex())))))))));

  auto pid = process_id(pg);
  	gen_end = get_time();
        if(pid == 0) {
		std::cout << "INFO: bfs time: " << print_time(gen_end - gen_start) << std::endl;
		f_log<<print_time(gen_end - gen_start) <<",";
	}
	gen_start = get_time();
  // Add Intercluster Edges using unordered flat map

  BGL_FORALL_VERTICES(v, g, Graph) {
	  BGL_FORALL_OUTEDGES(v, e, g, Graph) {
   		 request(cluster, target(e,g));
 		 request(distance, target(e,g)); 
  	  }
  }
  synchronize(cluster);
  synchronize(distance);

  	gen_end = get_time();
        if(pid == 0) {
		std::cout << "INFO: request time: " << print_time(gen_end - gen_start) << std::endl;
		f_log<<print_time(gen_end - gen_start) <<",";
	}
	gen_start = get_time();
  size_t spanner_edges = 0;
  size_t spanner_edges_intra = 0;
  size_t spanner_edges_inter = 0;
  typedef boost::unordered_flat_map<size_t, size_t, 
	                            boost::hash<size_t>, std::equal_to< size_t>,
				    std::allocator<std::pair<const  size_t, size_t>>> MyMap;

  MyMap hashmap;
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
  			spanner_edges_inter++;
      	  	}	
       	  }
    }
   } 
  }

  	gen_end = get_time();
        if(pid == 0) {
		std::cout << "INFO:Intercluster add time: " << print_time(gen_end - gen_start) << std::endl;
		f_log<<print_time(gen_end - gen_start) <<",";
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
        if(get(spanner_intra,e) ==  true)
	 { 
  		spanner_edges++;
  		spanner_edges_intra++;
	 }

      }
  }
  	gen_end = get_time();
        if(pid == 0) {
		std::cout << "INFO: count intracluster time: " << print_time(gen_end - gen_start) << std::endl;
		f_log<<print_time(gen_end - gen_start) <<",";
	}
	gen_start = get_time();
  auto spg =  process_group(sg);
  synchronize(pg);
  synchronize(spg);

  //  all reduce
  size_t global_spanner_edges = boost::parallel::all_reduce(pg, spanner_edges, std::plus<size_t>());
  size_t global_spanner_edges_intra = boost::parallel::all_reduce(pg, spanner_edges_intra, std::plus<size_t>());
  size_t global_spanner_edges_inter = boost::parallel::all_reduce(pg, spanner_edges_inter, std::plus<size_t>());
  size_t total;
  total =  boost::parallel::all_reduce(pg, num_vertices(g), std::plus<size_t>());
  if (process_id(process_group(g)) == 0) {
  	std::cout<<"Num of Vertices in Original Graph is "<<total<<"\n";
	f_log<<total<<",";
   }	
  total =  boost::parallel::all_reduce(pg, num_edges(g), std::plus<size_t>());
  if (process_id(process_group(g)) == 0) {
   	std::cout<<"Num of Edges in Original Graph is "<<total<<"\n";
	f_log<<total<<",";
   }	


  if (process_id(process_group(sg)) == 0) {
  	std::cout<<"Spanner Edges from original graph "<<global_spanner_edges<<"\n";
  	f_log<<",,,"<<global_spanner_edges<<","<<global_spanner_edges_intra<<","<<global_spanner_edges_inter<<",";
  }
  total =  boost::parallel::all_reduce(pg, num_vertices(sg), std::plus<size_t>());

  if (process_id(process_group(g)) == 0)  
	std::cout<<"Num of Vertices in Spanner Graph created copy is "<<total<<"\n";

  total =  boost::parallel::all_reduce(pg, num_edges(sg), std::plus<size_t>());
  if (process_id(process_group(sg)) == 0) {
  	std::cout<<"Num of edges in Spanner Graph created copy is "<<total<<"\n";
        std::cout.flush();
  }	
  synchronize(pg);
}

void 
cleanup(Graph &g) {
	// Internal property map to store whether vertex is center or not
	property_map<Graph,  vertex_distance_t>::type
  		 distance = get(vertex_distance, g);

        // Internal property map to store cluster of vertex
	property_map<Graph,  cluster_t>::type
  		 cluster = get(cluster_t(), g);
 	property_map<Graph,  vertex_predecessor_t>::type
  		 pred = get(vertex_predecessor, g);  // Deposit	

	property_map<Graph,  parent_t>::type
  		 parent = get(parent_t(), g);

	property_map<Graph,  pred_numeric_t>::type
  		 pred_numeric = get(pred_numeric_t(), g);


	// Internal property map to store more
	property_map<Graph,  map_t>::type
  		 map = get(map_t(), g);
  // Internal property map to store whether edge is part of spanner or not (Tree)
  property_map<Graph,  spanner_intra_t>::type
    spanner_intra = get(spanner_intra_t(), g);

  // Initialize edges to non spanner edges
  BGL_FORALL_VERTICES(v, g, Graph) {
	  BGL_FORALL_OUTEDGES(v, e, g, Graph) {
 		 put(spanner_intra, e, false);
  	  }
  } 
  auto pg = process_group(g);
  size_t n = num_vertices(g);
  size_t global_n = boost::parallel::all_reduce(pg, n, std::plus<size_t>());
  auto v_invalid = vertex(global_n -1,g);
  	BGL_FORALL_VERTICES(v, g, Graph) {
    		put(parent, v, v_invalid);
    		put(distance, v, (std::numeric_limits<std::size_t>::max)());
   		put(pred_numeric, v, (std::numeric_limits<std::size_t>::max)());
    		put(cluster, v, v_invalid);
    		put(pred, v, v_invalid);

  	}
}
int
main(int argc, char* argv[]){
	boost::mpi::environment env(argc, argv);
	mpi_process_group pg;
  	mpi_process_group::process_id_type id = process_id(pg);	
	std::ofstream f_log ;
	if(id == 0) 
		f_log.open ("All_log.txt", std::ios_base::app);
	time_type gen_start, gen_end;
	time_type gen_start_total, gen_end_total;
	gen_start_total = get_time();

	gen_start = get_time();

	srand(time(0));
	int n = 1000; 
	long m = 4000;
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
	int dig = 1;
	if(argc > 1)
	{ 
		int num = 0;
		int i = 0;
		while(argv[dig][i]!='\0') {
		   num = (num *10) + (argv[dig][i] - '0');
		   i++;
		}
	        n = num;
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
    	Distribution distrib = parallel::block(pg, n);
  	typedef keep_local_edges<parallel::variant_distribution<mpi_process_group>,
                           mpi_process_group::process_id_type>
    			   EdgeFilter; 

	typedef scalable_rmat_iterator<mpi_process_group, Distribution, rand48, Graph>
    			   RMATIter;

  	if (id == 0) printf("INFO: Generating graph.\n");

  	rand48 gen;
  	Graph g(RMATIter(pg, distrib, gen, n, m, a, b, c, d, true),
          	RMATIter(), n, pg, distrib);

	synchronize(g);
  	gen_end = get_time();
        if(id == 0){
		std::cout << "INFO: Graph Gen time: " << print_time(gen_end - gen_start) << std::endl;
		f_log<<print_time(gen_end - gen_start) <<"\n";
	}
	
	//parameter k input
	int k = 3; //Default - FGV5
	int k_arg = 3; //3rd argument is k
	if(argc > k_arg)
		k = argv[k_arg][0] - '0';
	const int no_of_runs = 3;
	for(int run_id = 0 ; run_id < no_of_runs ; run_id++) {
 
	SpannerGraph sg;
        if(id == 0)
		std::cout<<"\n\n\n    CONSTRUCT GRAPH SPANNER \n\n\n";
  	gen_start = get_time();
  	ldd(g, k, sg, f_log);
  	gen_end = get_time();
        if(id == 0) {
		std::cout << "INFO: Spanner Construct time: " << print_time(gen_end - gen_start) << std::endl;
		f_log<<print_time(gen_end - gen_start) <<",";
	}
	gen_end_total = get_time();
        if(id == 0) {
		std::cout << "INFO: Total Time: " << print_time(gen_end_total - gen_start_total) << std::endl;
		f_log<<print_time(gen_end_total - gen_start_total) <<"\n";
	}
	f_log.flush();
	cleanup(g);
	}
	return 0;
	//return boost::report_errors();
}
