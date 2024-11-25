#include <iostream>
#include <vector>
#include <queue>
#include <map>
#include <limits>

#include <algorithm>
#ifndef ROAD_NETWORKS_H
#define ROAD_NETWORKS_H
#include "road_network.h"
#endif // ROAD_NETWORKS_H
using road_network::Neighbor;
using road_network::NodeID;
using road_network::distance_t;

using std::cout;
using std::endl;
using std::vector;
using std::map;
using std::pair;
using std::numeric_limits;
using std::priority_queue;


// typedef uint32_t NodeID;
// typedef uint32_t distance_t;

// struct Neighbor {
//     NodeID node;
//     distance_t distance;
//     int edge_spc;
//     Neighbor(NodeID node, distance_t distance, int edge_spc);
// };

// Constructor implementation
// Neighbor::Neighbor(NodeID node, distance_t distance, int edge_spc) : 
//     node(node), distance(distance), edge_spc(edge_spc) {}

const int INF = numeric_limits<int>::max();

// Structure to represent an edge
struct NodeInfo {
    int nodeID;
    int distance;
    int spc;
    std::vector<int> pred;
    NodeInfo(int id, int dist, int spc) : nodeID(id), distance(dist), spc(spc) {}
    NodeInfo(){}
};


// Graph class
class SimpleGraph {
    map<int, vector<NodeInfo>> nodes;
    map<int,NodeInfo> node_data; // tmp data for dijkstra
    

public:
    void clear() {
        nodes.clear();
        node_data.clear();
    }

    // Add a node to the graph
    void addNode(int id);

    // Add a directed edge between two nodes
    void addEdge(int src, int dest, int distance, int spc, bool add_reverse);

    void updateEdge(int src, int dest, int distance, int spc);

    map<int, NodeInfo> dijkstra(int source);

    map<int, NodeInfo> dijkstra_exclude_given_Node(int source, int exclude_node);
    /*
        for two borders, sd can achieve via
        1. the current cut 
        2. via partition (since via cut level is dijk from cut nodes)
    */   
    map<pair<NodeID, NodeID>, pair<distance_t,int>> create_min_connect_graph (NodeID cut, std::vector<Neighbor> &rm_node_nbrs,
                                map<pair<NodeID, NodeID>, distance_t> &dist_via_partition,
                                map<pair<NodeID, NodeID>, distance_t> &dist_via_cut,
                                std::set<NodeID> &cut_sets,std::map<std::pair<NodeID, NodeID>, std::pair<distance_t, int>> &tmp_edges,
                                bool debug_flag);

    map<pair<NodeID, NodeID>, pair<distance_t,int>> create_min_connect_graph (NodeID cut, std::vector<Neighbor> &rm_node_nbrs,
                                map<pair<NodeID, NodeID>, distance_t> &dist_via_partition,
                                map<pair<NodeID, NodeID>, distance_t> &dist_via_cut,
                                bool debug_flag);
    

};
