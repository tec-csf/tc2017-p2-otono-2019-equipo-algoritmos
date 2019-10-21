#include <iostream>
#include <stack>
#include <map>
#include <queue>
#include <vector>
#include <iostream>
#include <chrono>
// Librerias de Boost
#include <boost/graph/adjacency_list.hpp> 
#include <boost/graph/graph_traits.hpp> 
#include <boost/graph/graphviz.hpp> 
#include <boost/tuple/tuple.hpp> 
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>

#include <stdio.h>
#include<stdlib.h>


#include<time.h>
#include <ctime> 





using namespace boost;
using namespace std;

typedef property<edge_weight_t, int> EdgeWeightProperty; 
typedef adjacency_list<listS, vecS, directedS, no_property, EdgeWeightProperty> grafoInicial; 
typedef boost::adjacency_list<boost::vecS, boost::hash_setS, boost::undirectedS, uint32_t, uint32_t, boost::no_property> graph_t;
typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;

typedef boost::graph_traits<grafoInicial>::out_edge_iterator oe_it;
typedef boost::graph_traits<grafoInicial>::edge_iterator e_it;
typedef boost::graph_traits<grafoInicial>::vertex_iterator v_it;
typedef boost::graph_traits<grafoInicial>::vertex_descriptor Vertex;
typedef boost::graph_traits<grafoInicial>::edge_descriptor Edge;


grafoInicial g(13);
void AddEdge(int s, int d, int w, grafoInicial &a);
void EraseN(int cual, grafoInicial &a);
void AddN(grafoInicial &grafo);
void EraseEdge(int s, int d, grafoInicial &graf);
void AddEdgeDiff(int s, int d, int w, grafoInicial &a);
void dfs(grafoInicial &a, unsigned long startingNode);
void bfs(grafoInicial &a, unsigned long startingNode);
void kruskal();
void prim();
void dijkstra();
vector<vector<int>> floydWarshall(grafoInicial g);

struct oper{
    bool operator()(const std::tuple<int,Vertex, Vertex>& one, const std::tuple<int,Vertex, Vertex>& two) { return std::get<0>(one) > std::get<0>(two); }
};

struct oper2{
    bool operator()(const Edge& one, const Edge& two){ return boost::get(boost::edge_weight_t(), g, one) > boost::get(boost::edge_weight_t(), g, two); }
};

int main()
{
    AddEdge(1 ,3, 8, g);
    AddEdge(1, 4, 8, g);
	AddEdge(2, 5, 7, g);
    AddEdge(3, 2, 7, g);
    AddEdge(3, 5, 8, g);
    AddEdge(3, 10, 4, g);
    AddEdge(4, 7, 3, g);
    AddEdge(4, 5, 1, g);
    AddEdge(4, 8, 2, g);
    AddEdge(5, 6, 9, g);
    AddEdge(6, 13, 4, g);
    AddEdge(7, 4, 6, g);
    AddEdge(8, 7, 3, g);
    AddEdge(8, 9, 3, g);
    AddEdge(9, 12, 4, g);
    AddEdge(9, 10, 2, g);
    AddEdge(10, 6, 6, g);
    AddEdge(10, 3, 10, g);
    AddEdge(11, 12, 6, g);
    AddEdge(12, 9, 2, g);
    AddEdge(12, 14, 9, g);
    AddEdge(12, 11, 8, g);
    AddEdge(13, 14, 6, g);
    
	unsigned t0_insertar, t1_insertar;
	t0_insertar=clock();    //Corre tiempo

	AddEdge(14, 13, 2, g);
    graph_traits<grafoInicial>::vertex_iterator vi, vi_fin, sig; 
    
    t1_insertar = clock();  //Termina de ejecutarse
	
	double time_insertar = (double(t1_insertar-t0_insertar)/CLOCKS_PER_SEC); //Se calcula el tiempo
	
    unsigned t0_insertarV, t1_insertarV;
	t0_insertarV=clock();    //Corre tiempo
	
    tie(vi, vi_fin) = vertices(g); 
    t1_insertarV = clock();  //Termina de ejecutarse
	
	double time_insertarV = (double(t1_insertarV-t0_insertarV)/CLOCKS_PER_SEC); //Se calcula el tiempo
	
    write_graphviz(cout, g);
    unsigned t0_erase, t1_erase;
	t0_erase=clock();    //Corre tiempo
	EraseEdge(13,14,g);
    t1_erase = clock();  //Termina de ejecutarse
	
	double time_erase = (double(t1_erase-t0_erase)/CLOCKS_PER_SEC); //Se calcula el tiempo
   
    
    AddEdge(13,14,6,g);/*Lo retornamos*/
    
    
    
    cout<<"\nRecorrido de amplitud:\n";
    
	

	
	unsigned t0_bfs, t1_bfs;
	t0_bfs=clock();    //Corre tiempo
	bfs(g, 1); 
	t1_bfs = clock();  //Termina de ejecutarse
	
	double time_bfs = (double(t1_bfs-t0_bfs)/CLOCKS_PER_SEC); //Se calcula el tiempo
	
	
    cout<<"\nRecorrido de profundidad:\n";
    unsigned t0_dfs, t1_dfs;
	t0_dfs=clock();    //Corre tiempo
	dfs(g, 1); 
	t1_dfs = clock();  //Termina de ejecutarse
	
	double time_dfs = (double(t1_dfs-t0_dfs)/CLOCKS_PER_SEC); //Se calcula el tiempo
    
    unsigned t0_eraseN, t1_eraseN;
	t0_eraseN=clock();    //Corre tiempo
	EraseN(13,g);
    t1_eraseN = clock();  //Termina de ejecutarse
	
	double time_eraseN = (double(t1_eraseN-t0_eraseN)/CLOCKS_PER_SEC); //Se calcula el tiempo
    
    cout<<"\nAlgoritmo de Prim:\n";
    prim();
    
	cout<<"\nAlgoritmo de Kruskal:\n";
    kruskal(); 
    
	cout<<"\nAlgoritmo de Dijkstra:\n";
    dijkstra(); 

    cout<<"\nAlgoritmo de Floyd Warshall:\n";
    vector<vector<int>> FloydWarhsall = floydWarshall(g); 
    
    
    
    std::cout << "\nExecution Time Insert: " << time_insertar << std::endl; //Imprimimos tiempo de ejecucion
    std::cout << "\nExecution Time Insert V: " << time_insertarV << std::endl; //Imprimimos tiempo de ejecucion
    std::cout << "\nExecution Time Erase Edge: " << time_erase << std::endl; //Imprimimos tiempo de ejecucion
    std::cout << "\nExecution Time Erase Node: " << time_eraseN << std::endl; //Imprimimos tiempo de ejecucion
    std::cout << "\nExecution Time BFS: " << time_bfs << std::endl; //Imprimimos tiempo de ejecucion
    std::cout << "\nExecution Time DFS: " << time_dfs << std::endl; //Imprimimos tiempo de ejecucion
    return 0;
}

void AddEdge(int s, int d, int w, grafoInicial &a) { 
add_edge(s, d, w, a); 
}

void AddEdgeDiff(int s, int d, int w, grafoInicial &a)
{
    add_edge(s,d,w,a);
    add_edge(d,s,w,a);
}

void EraseEdge(int s, int d, grafoInicial &graf) { 
remove_edge(s,d,graf); 

}

void EraseN(int cual, grafoInicial &a) { 
remove_vertex(cual, a); 

}

void AddN(grafoInicial &grafo) { 
add_vertex(grafo); 
}

void dfs(grafoInicial &a, unsigned long startingNode)
{
    stack<Vertex> visiting;
    map<Vertex,bool> travelled;
    Vertex checking = a.vertex_set()[startingNode];
    visiting.push(checking);
    travelled[checking] = true;
    while(!visiting.empty())
    {
        checking = visiting.top();
        visiting.pop();
        pair<oe_it,oe_it> iterators = out_edges(checking,a);
        for(oe_it it = iterators.first; it != iterators.second; ++it)
        {
            Vertex target = boost::target(*it,a);
            if(travelled[target]) continue;
            travelled[target] = true;
            visiting.push(target);
        }
        cout << checking << "-";

    }
    cout<<"\n";

}

void bfs(grafoInicial &a, unsigned long startingNode){
    std::queue<Vertex> visiting;
    map<Vertex,bool> travelled;
    Vertex checking = a.vertex_set()[startingNode];
    visiting.push(checking);
    travelled[checking] = true;
    while(!visiting.empty()){
        checking = visiting.front();
        visiting.pop();
        pair<oe_it,oe_it> iterators = out_edges(checking,a);
        for(oe_it it = iterators.first; it != iterators.second; ++it){
            Vertex target = boost::target(*it,a);
            if(travelled[target]) continue;
            travelled[target] = true;
            visiting.push(target);
        }
        cout << checking << "-";
    }
    cout<<"\n";
}

void prim()
{
    using namespace boost;
  typedef adjacency_list < vecS, vecS, undirectedS,
    property<vertex_distance_t, int>, property < edge_weight_t, int > > Graph;
  typedef std::pair < int, int >E;
  const int num_nodes = 14;
  E edges[] = { 

E(0, 3),
E(0, 2),
E(1, 2),
E(1, 4),
E(2, 9),
E(2, 4),
E(3, 6),
E(3, 4),
E(3, 7),
E(4, 5),
E(5, 9),
E(5, 12),
E(6, 7),
E(7, 8),
E(8, 9),
E(8, 11	),
E(10, 11),
E(11, 13),
E(12, 13)

  };

  int weights[] = { 8,8,7,7,14,8,9,1,2,9,6,4,3,3,2,4,14,9,8 };


	unsigned t0, t1;
	t0=clock();    //Corre tiempo

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
  Graph g(num_nodes);

  for (std::size_t j = 0; j < sizeof(edges) / sizeof(E); ++j) {
    graph_traits<Graph>::edge_descriptor e; bool inserted;
    boost::tie(e, inserted) = add_edge(edges[j].first, edges[j].second, g);
    weightmap[e] = weights[j];
  }
#else
  Graph g(edges, edges + sizeof(edges) / sizeof(E), weights, num_nodes);
  property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
#endif
  std::vector < graph_traits < Graph >::vertex_descriptor >
    p(num_vertices(g));

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
  property_map<Graph, vertex_distance_t>::type distance = get(vertex_distance, g);
  property_map<Graph, vertex_index_t>::type indexmap = get(vertex_index, g);
  

  prim_minimum_spanning_tree
    (g, *vertices(g).first, &p[0], distance, weightmap, indexmap, 
     default_dijkstra_visitor());
#else
  prim_minimum_spanning_tree(g, &p[0]);
#endif

	  std::cout<<"\n";



	
	
  for (std::size_t i = 0; i != p.size(); ++i)
    if (p[i] != i)
      std::cout << "parent[" << i << "] = " << p[i] << std::endl;
    else
      std::cout << "parent[" << i << "] = no parent" << std::endl;


	t1 = clock();  //Termina de ejecutarse
	
	double time = (double(t1-t0)/CLOCKS_PER_SEC); //Se calcula el tiempo
	std::cout << "\nExecution Time: " << time << std::endl; //Imprimimos tiempo de ejecucion
	
}

void kruskal()
{
    using namespace boost;
  typedef adjacency_list < vecS, vecS, undirectedS,
    no_property, property < edge_weight_t, int > > Graph;
  typedef graph_traits < Graph >::edge_descriptor Edge;
  typedef std::pair<int, int> E;

  const int num_nodes = 14;
  E edge_array[] = { 

E(0, 3),
E(0, 2),
E(1, 2),
E(1, 4),
E(2, 9),
E(2, 4),
E(3, 6),
E(3, 4),
E(3, 7),
E(4, 5),
E(5, 9),
E(5, 12),
E(6, 7),
E(7, 8),
E(8, 9),
E(8, 11	),
E(10, 11),
E(11, 13),
E(12, 13)

  };

  int weights[] = { 8,8,7,7,14,8,9,1,2,9,6,4,3,3,2,4,14,9,8 };

  std::size_t num_edges = sizeof(edge_array) / sizeof(E);
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
  Graph g(num_nodes);
  property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
  for (std::size_t j = 0; j < num_edges; ++j) {
    Edge e; bool inserted;
    boost::tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
    weightmap[e] = weights[j];
  }
#else
  Graph g(edge_array, edge_array + num_edges, weights, num_nodes);
#endif
  property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
  std::vector < Edge > spanning_tree;





unsigned t0_kruskal, t1_kruskal;
	t0_kruskal=clock();    //Corre tiempo
	
  kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
  


  std::cout << "Print the edges in the MST:" << std::endl;
  for (std::vector < Edge >::iterator ei = spanning_tree.begin();
       ei != spanning_tree.end(); ++ei) {
    std::cout << source(*ei, g) << " <--> " << target(*ei, g)
      << " with weight of " << weight[*ei]
      << std::endl;
  }
  
    t1_kruskal = clock();  //Termina de ejecutarse
	
	double time_kruskal = (double(t1_kruskal-t0_kruskal)/CLOCKS_PER_SEC); //Se calcula el tiempo
	std::cout << "\nExecution Time: " << time_kruskal << std::endl; //Imprimimos tiempo de ejecucion

  std::ofstream fout("figs/kruskal-eg.dot");
  fout << "graph A {\n"
    << " rankdir=LR\n"
    << " size=\"3,3\"\n"
    << " ratio=\"filled\"\n"
    << " edge[style=\"bold\"]\n" << " node[shape=\"circle\"]\n";
  graph_traits<Graph>::edge_iterator eiter, eiter_end;
  for (boost::tie(eiter, eiter_end) = edges(g); eiter != eiter_end; ++eiter) {
    fout << source(*eiter, g) << " -- " << target(*eiter, g);
    if (std::find(spanning_tree.begin(), spanning_tree.end(), *eiter)
        != spanning_tree.end())
      fout << "[color=\"black\", label=\"" << get(edge_weight, g, *eiter)
           << "\"];\n";
    else
      fout << "[color=\"gray\", label=\"" << get(edge_weight, g, *eiter)
           << "\"];\n";
  }
  fout << "}\n";

}

void dijkstra()
{
 typedef adjacency_list < listS, vecS, directedS,
    no_property, property < edge_weight_t, int > > graph_t;
  typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
  typedef std::pair<int, int> Edge;

  
  //ASI SE AGREGAN NODOS
  const int num_nodes = 14;
 
  enum nodes {A,B,C,D,E,F,G,H,I,J,K,L,M,N};
  char name[] = "ABCDEFGHIJKLMN";
  Edge edge_array[] = { 

Edge(A, D),
Edge(A, C),
Edge(B, C),
Edge(B, E),
Edge(C, J),
Edge(C, E),
Edge(D, G),
Edge(D, E),
Edge(D, H),
Edge(E, F),
Edge(F, J),
Edge(F, M),
Edge(G, H),
Edge(H, I),
Edge(I, J),
Edge(I, L),
Edge(K, L),
Edge(L, N),
Edge(M, N)
  };
    int weights[] = {8,8,7,7,14,8,9,1,2,9,6,4,3,3,2,4,14,9,8};

 
  
  
  
  int num_arcs = sizeof(edge_array) / sizeof(Edge);
  graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
  property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
  std::vector<vertex_descriptor> p(num_vertices(g));
  std::vector<int> d(num_vertices(g));
  vertex_descriptor s = vertex(A, g);
  
  
  
  
  
  
  unsigned t0_dj, t1_dj;
	t0_dj=clock();    //Corre tiempo
	
	
	
  
  dijkstra_shortest_paths_no_color_map(g, s,
                                       predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
                                       distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));

	
	

  std::cout << "distances and parents:" << std::endl;
  graph_traits < graph_t >::vertex_iterator vi, vend;
  for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) {
    std::cout << "distance(" << name[*vi] << ") = " << d[*vi] << ", ";
    std::cout << "parent(" << name[*vi] << ") = " << name[p[*vi]] << std::
      endl;
                                             
    
  }
  
  t1_dj = clock();  //Termina de ejecutarse
	
	double time_dj = (double(t1_dj-t0_dj)/CLOCKS_PER_SEC); //Se calcula el tiempo
	std::cout << "\nExecution Time: " << time_dj << std::endl; //Imprimimos tiempo de ejecucion
  std::cout << std::endl;

}

vector<vector<int>> floydWarshall(grafoInicial g)
{
	
	
	unsigned t0_floyd, t1_floyd;
	t0_floyd=clock();    //Corre tiempo
	
	
    vector<vector<int>> distances;
    unsigned long vNum = num_vertices(g);
    distances.resize(vNum);
    pair<v_it,v_it> iterators = vertices(g);
    for(v_it it = iterators.first; it != iterators.second; ++it)
    {
        for(v_it it2 = iterators.first; it2 != iterators.second; ++it2)
        {
            if(it == it2)
            {
                distances[*it].push_back(0);
                continue;
            }
            else distances[*it].push_back(100);
        }
    }
    pair<e_it,e_it> edIt = edges(g);
    for(e_it it = edIt.first; it != edIt.second; ++it)
    {
        distances[boost::source(*it,g)][boost::target(*it,g)] = boost::get(boost::edge_weight_t(),g, *it);
    }

    for(int k = 0; k < vNum; ++k){
        for(int i = 0; i < vNum; ++i){
            for(int j = 0; j < vNum; ++j){
                if(distances[i][j] > distances[i][k] + distances[k][j])
                    distances[i][j] = distances[i][k] + distances[k][j];
            }
        }
    }
    
    
    t1_floyd = clock();  //Termina de ejecutarse
	
	double time_floyd = (double(t1_floyd-t0_floyd)/CLOCKS_PER_SEC); //Se calcula el tiempo
	std::cout << "\nExecution Time: " << time_floyd << std::endl; //Imprimimos tiempo de ejecucion
	
	
	
    
    
    for (int r = 0; r<vNum;++r){
    	for (int s = 0; s<vNum; ++s){
    		cout<<" || " << distances[r][s];
		}
		cout<<"\n";
	}
	
	cout<<"\n\n\n";
    return distances;
}
