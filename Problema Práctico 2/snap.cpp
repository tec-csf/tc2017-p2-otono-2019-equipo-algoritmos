/***************************************
 ** Daniela Vignau |   Allan Sánchez  **
 **    A01021698   |      A01379951   **
 ***************************************/

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <iomanip>
#include <bits/stdc++.h>
#include <chrono>

#define INF 1000000
using namespace std;

const int TOTAL_NODES = 15;
PNGraph graph;
int matValue[TOTAL_NODES][TOTAL_NODES];

void insertNode(int id, PNGraph graph) {
  graph->AddNode(id);
}

void insertVertex(int srcId, int desId, int weight, PNGraph graph) {
  graph->AddEdge(srcId, desId);
  matValue[srcId][desId] = weight;
}

void deleteNode(int id, PNGraph graph) {
  graph->DelNode(id);
}

void deleteVertex(int srcId, int desId, PNGraph graph) {
  graph->DelEdge(srcId, desId);
  matValue[srcId][desId] = 0;
}

void sortChilds(vector<pair <int,int> >& adyacents) {
    sort(adyacents.begin(), adyacents.end()/*,
        [](const pair<int,int> & a, const pair<int,int> & b) -> bool
    {
    return a.second > b.second;
    }*/);
}

vector<pair <int,int> > getAdyacentChildsID(int id) {
   vector<pair <int,int> > adyacents;
    TNGraph::TNodeI nodeI = graph->GetNI(id);
    for (int i = 1; i < TOTAL_NODES; ++i){
        if (nodeI.IsOutNId(i)){
            adyacents.push_back(make_pair(i,matValue[id][i]));
        }
    }
    sortChilds(adyacents);
    return adyacents;
}


struct DisjointSets {
    int *parent, *rnk;
    int n;

    DisjointSets(int n)   {
        this->n = n;
        parent = new int[n+1];
        rnk = new int[n+1];
        for (int i = 0; i <= n; ++i)  {
            rnk[i] = 0;
            parent[i] = i;
        }
    }

    int findSet(int x)  {
        if (x != parent[x])
            parent[x] = findSet(parent[x]);
        return parent[x];
    }

    void unionSet(int x, int y)   {
        x = findSet(x), y = findSet(y);

        if (rnk[x] > rnk[y]) parent[y] = x;
        else parent[x] = y;

        if (rnk[x] == rnk[y]) rnk[y]++;
    }
};

struct oper{
    bool operator()(const std::tuple<int,int,int>& one, const std::tuple<int,int,int>& two)
    {
        return std::get<0>(one) > std::get<0>(two);
    }
};

// ALGORITHMS

void bfs(int startId) {
   //COMPLEXITY: O(|y| + |E|)
    queue<pair<int,int> > q;
    bool visited[TOTAL_NODES];
    fill(visited, visited+TOTAL_NODES, false);
    q.push(make_pair(startId,0));

    while (!q.empty()){
        int actId = q.front().first;
        q.pop();
        if (!visited[actId]){
            cout << actId << " - ";
            vector<pair<int,int> > adyacents = getAdyacentChildsID(actId);
            for (int i = 0; i < adyacents.size(); ++i){
                q.push(adyacents.at(i));
            }
        }
        visited[actId] = true;
    }
}

void dfs(bool visited[], int actId) {
   //COMPLEXITY: O(|y| + |E|)
    if (!visited[actId]){
        cout << actId << " - ";
        visited[actId] = true;
        vector<pair <int,int> > adyacents = getAdyacentChildsID(actId);

        for (int i = 0; i < adyacents.size(); ++i)
            dfs(visited, adyacents.at(i).first);
    }
}


void dijkstra(int startId) {
   //COMPLEXITY: O(|y|2)
    queue<pair<int,int> > q;
    queue<pair<int,int> > ps;
    bool visited[TOTAL_NODES];
    int distances[TOTAL_NODES];
    fill(visited, visited+TOTAL_NODES, false);
    fill(distances, distances+TOTAL_NODES,INF);
    int parentId, distanceToNode;

    q.push(make_pair(startId,0));
    distances[startId] = 0;

cout << "Moving to    Parent\tT. weight" << endl;
    while (!q.empty()) {
        int actId = q.front().first;
        q.pop();
        if (!visited[actId]) {
            vector<pair<int,int> > adyacents = getAdyacentChildsID(actId);
            for (int i = 0; i < adyacents.size(); ++i) {
            parentId = actId;
                if (distances[actId] + adyacents.at(i).second < distances[adyacents.at(i).first]) {
                    distances[adyacents.at(i).first] = distances[actId] + adyacents.at(i).second;
                    q.push(adyacents.at(i));
                    ps.push(make_pair(distances[i], parentId));
                    distanceToNode = distances[adyacents.at(i).first];
                }
                cout << "     " << adyacents.at(i).first << "\t\t" << parentId << "\t   " << distanceToNode << endl;
            }
        }
        visited[actId] = true;
    }

  cout << "\nOrigin\tWeight" << endl;
  int i = 1;
  while(!ps.empty()) {
   if(i == startId) {
     cout << "   " << i << "\t  inf" << endl;
     ++i;
   } else {
     cout << "   " << i << "\t  "  << distances[i] << endl;
     ps.pop();
     ++i;
   }
   if(i == TOTAL_NODES) break;
  }
}

void prim(int startId) {
   //COMPLEXITY: O(|y|3)
    queue<pair<int,int> > q;
    bool visited[TOTAL_NODES];
    int distances[TOTAL_NODES];
    queue<pair<int,int> > ps;
    fill(visited, visited+TOTAL_NODES, false);
    fill(distances, distances+TOTAL_NODES,INF);
    int distanceToNode, parentId;

    q.push(make_pair(startId,0));
    distances[startId] = 0;

    cout << "Parent \tMoving to    \tT. weight" << endl;
    while (!q.empty()) {
      int actId = q.front().first;
      q.pop();
      vector<pair<int,int> > adyacents = getAdyacentChildsID(actId); // getting adyacent children from actId
      parentId = actId;
      for (int i = 0; i < adyacents.size(); ++i) {
        //cout << "At node: " << parentId << endl;
        vector<pair<int,int> > childAdy = getAdyacentChildsID(adyacents.at(i).first);
        for(int j = 0; j < childAdy.size(); ++j) {
          //cout << "Checking adyacents from " << parentId << ": " << adyacents.at(i).first << endl;
          //cout << "Checking adyecents from " << adyacents.at(i).first << " is: " << childAdy.at(j).first << endl;
          if(childAdy.at(j).first != parentId) {
          //cout << "There are no adyacents from " << childAdy.at(j).first << " to " << parentId << endl;
            if (adyacents.at(i).second < distances[adyacents.at(i).first]) {
              //cout << "Distance: " << adyacents.at(i).second << endl;
              distances[adyacents.at(i).first] = adyacents.at(i).second;
              q.push(adyacents.at(i));
              ps.push(make_pair(distances[i],adyacents.at(i).second));
              distanceToNode = adyacents.at(i).second;
            }
          } else {
            distanceToNode = adyacents.at(i).second;
            //cout << "Distance: " << adyacents.at(i).second << endl;
          }
        }
        cout << "     " << parentId << "\t\t" << adyacents.at(i).first << "\t   " << distanceToNode << endl;
      }
    }

// cout << "Origin\tWeight" << endl;
//    int i = 1;
//    while(!ps.empty()) {
//     if(i == startId) {
//   cout << "   " << i << "\t   -" << endl;
//   ++i;
//     } else {
//  cout << "   " << i << "\t   "  << distances[i] << endl;
//     ps.pop();
//     ++i;
//     }
//     if(i == TOTAL_NODES) break;
//    }

}

vector<std::tuple<int, int, int>> kruskal(PNGraph graph) {
   //COMPLEXITY: O(|E| log |E|)
  priority_queue<std::tuple<int,int,int>,vector<std::tuple<int,int,int>>, oper> next;
  DisjointSets ds(graph->GetNodes());
  vector<std::tuple<int, int, int>> pathFollowed;
  std::tuple<int,int,int> now;
  int x, y, set_o, set_t;

  for (TNEANet::TEdgeI EI = graph->BegEI(); EI < graph->EndEI(); E++i) {
    x = EI.GetSrcNId();
    y = EI.GetDstNId();
      next.push(graph->GetIntAttrDatE(EI, "weight"), x, y);
      int set_u = ds.findSet(x);
      int set_v = ds.findSet(y);
      cout << "SET x "  << set_u << endl;
      cout << "SET y "  << set_u << endl;
      if (set_u != set_v) {
          cout << x << " - " << y << endl;
          ds.unionSet(set_u, set_v);
      }
  }

  while(!next.empty()){
    now = next.top();
    next.pop();
    x = std::get<1>(now);
    y = std::get<2>(now);
    set_o = ds.findSet(x);
    set_t = ds.findSet(x);
    if (set_o != set_t){
        pathFollowed.push_back(make_tuple(std::get<0>(now), x, y));
        ds.unionSet(set_o, set_t);
    }
  }
  return pathFollowed;

}

void floydWarshall() {
  //COMPLEXITY: O(|y|3)
  int matrix[TOTAL_NODES][TOTAL_NODES];
  for (int j = 1; j < TOTAL_NODES; ++j) {
    for (int i = 1; i < TOTAL_NODES; ++i) {
      matrix[i][j] = matValue[i][j];
      if (matrix[i][j] == 0) matrix[i][j] = INF;
      if (i == j) matrix[i][j] = 0;
    }
  }

  for (int k = 1; k < TOTAL_NODES; ++k) {
    for (int j = 1; j < TOTAL_NODES; ++j) {
      for (int i = 1; i < TOTAL_NODES; ++i) {
        if (matrix[k][j] + matrix[i][k] < matrix[i][j]) {
          matrix[i][j] = matrix[k][j] + matrix[i][k];
        }
      }
    }
  }

  for (int j = 1; j < TOTAL_NODES; ++j) {
    for (int i = 1; i < TOTAL_NODES; ++i) {
      if (matrix[i][j] == INF) cout << "oo";
      else {
        if (matrix[i][j] == 0) cout << "  ";
        else cout << setfill('0') << setw(2) << matrix[i][j];
      }
      cout << " ";
    }
    cout << endl;
  }

}

int main(int argc, char * argv[]) {
  typedef PNGraph PGraph; //  creating directed graph
  graph = TNGraph::New();
  double elapsedTime;

  //Inserting 14 total nodes from given graph
  for (int n = 0; n < 15; n++)
    insertNode(n, graph);

  cout << "\t\tExecution times (ms)" << endl;

 //These are simply used to calculate how long it takes to add or delete a node or a vertex
  chrono::high_resolution_clock::time_point startIN = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  insertNode(16, graph);
  chrono::high_resolution_clock::time_point endIN = chrono::high_resolution_clock::now();
  cout << "Insert Node:\t        " << chrono::duration_cast<chrono::milliseconds>(endIN - startIN).count() << endl;

  chrono::high_resolution_clock::time_point startDN = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  deleteNode(16, graph);
  chrono::high_resolution_clock::time_point endDN = chrono::high_resolution_clock::now();
  cout << "Delete Node:\t        " << chrono::duration_cast<chrono::milliseconds>(endDN - startDN).count() << endl;

  insertNode(16, graph);
  insertNode(17, graph);

  chrono::high_resolution_clock::time_point startIV = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  insertVertex(16, 17, 9, graph);
  chrono::high_resolution_clock::time_point endIV = chrono::high_resolution_clock::now();
  cout << "Insert Edge:\t        " << chrono::duration_cast<chrono::milliseconds>(endIV - startIV).count() << endl;

  chrono::high_resolution_clock::time_point startDV = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  deleteVertex(16, 17, graph);
  chrono::high_resolution_clock::time_point endDV = chrono::high_resolution_clock::now();
  cout << "Delete Edge:\t        " << chrono::duration_cast<chrono::milliseconds>(endDV - startDV).count() << endl;

  //Inserting vertices for given graph
  insertVertex(1, 3, 8, graph);
  insertVertex(1, 4, 8, graph);
  insertVertex(2, 5, 7, graph);
  insertVertex(3, 2, 7, graph);
  insertVertex(3, 5, 8, graph);
  insertVertex(3, 10, 4, graph);
  insertVertex(4, 5, 1, graph);
  insertVertex(4, 7, 3, graph);
  insertVertex(4, 8, 2, graph);
  insertVertex(5, 6, 9, graph);
  insertVertex(6, 13, 4, graph);
  insertVertex(7, 4, 6, graph);
  insertVertex(8, 7, 3, graph);
  insertVertex(8, 9, 3, graph);
  insertVertex(9, 10, 2, graph);
  insertVertex(9, 12, 4, graph);
  insertVertex(10, 3, 10, graph);
  insertVertex(10, 6, 6, graph);
  insertVertex(11, 12, 6, graph);
  insertVertex(12, 9, 2, graph);
  insertVertex(12, 11, 8, graph);
  insertVertex(12, 14, 9, graph);
  insertVertex(13, 14, 6, graph);
  insertVertex(14, 13, 2, graph);

  bool visited[TOTAL_NODES];
  fill(visited, visited+TOTAL_NODES, false);

  cout << "\nAlgorithms\n\nDFS" << endl;
  chrono::high_resolution_clock::time_point startDFS = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  dfs(visited, 1);
  chrono::high_resolution_clock::time_point endDFS = chrono::high_resolution_clock::now();
  elapsedTime = chrono::duration_cast<chrono::milliseconds>(endDFS - startDFS).count();
  cout << "\nElapsed Time: " << elapsedTime << setprecision(9) << endl;

  cout << "\n\nBFS" << endl;
  chrono::high_resolution_clock::time_point startBFS = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  bfs(1);
  chrono::high_resolution_clock::time_point endBFS = chrono::high_resolution_clock::now();
  elapsedTime = chrono::duration_cast<chrono::milliseconds>(endBFS - startBFS).count();
  cout << "\nElapsed Time: " << elapsedTime << setprecision(9) << endl;

  cout << "\n\nFloyd Warshall" << endl;
  chrono::high_resolution_clock::time_point startFW = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  floydWarshall();
  chrono::high_resolution_clock::time_point endFW  = chrono::high_resolution_clock::now();
  elapsedTime = chrono::duration_cast<chrono::milliseconds>(endFW - startFW).count();
  cout << "Elapsed Time: " << elapsedTime << setprecision(9) << endl;

  cout << "\n\nKruskal" << endl;
  chrono::high_resolution_clock::time_point startK = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  //vector<std::tuple<int, int, int>> Kruskal = kruskal(graph);
  chrono::high_resolution_clock::time_point endK = chrono::high_resolution_clock::now();
  elapsedTime = chrono::duration_cast<chrono::milliseconds>(endK - startK).count();
  cout << "Elapsed Time: " << elapsedTime << setprecision(9) << endl;

  cout << "\n\nPrim" << endl;
  chrono::high_resolution_clock::time_point startPrim = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  prim(1);
  chrono::high_resolution_clock::time_point endPrim = chrono::high_resolution_clock::now();
  elapsedTime = chrono::duration_cast<chrono::milliseconds>(endPrim - startPrim).count();
  cout << "Elapsed Time: " << elapsedTime << setprecision(9) << endl;

  cout << "\n\nDijkstra" << endl;
  chrono::high_resolution_clock::time_point startDij = chrono::high_resolution_clock::now();
  ios_base::sync_with_stdio(false);
  dijkstra(1);
  chrono::high_resolution_clock::time_point endDij = chrono::high_resolution_clock::now();
  elapsedTime = chrono::duration_cast<chrono::milliseconds>(endDij - startDij).count();
  cout << "Elapsed Time: " << elapsedTime << setprecision(9) << endl;

}

/*
FOR COMPILING --> IT MUST BE IN C++11 IN ORDER TO USE high_resolution_clock
g++ -std=c++11 -Wall -O3 -DNDEBUG -fopenmp -o testgraph testgraph.cpp  ../../snap-core/Snap.o -I../../snap-core -I../../snap-adv -I../../glib-core -I../../snap-exp  -lrt
*/
