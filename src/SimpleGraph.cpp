#include "SimpleGraph.h"


// Add a node to the graph
void SimpleGraph::addNode(int id) {
    if (nodes.find(id) == nodes.end()) {
        nodes[id] = vector<NodeInfo>();
    }
}

// Add a directed edge between two nodes
void SimpleGraph::addEdge(int src, int dest, int distance, int spc, bool add_reverse) {
    addNode(src);
    vector<NodeInfo>& nbrs = nodes[src];
    bool updated = false;
    for (int i = 0; i < nbrs.size(); ++i) {
        if (nbrs[i].nodeID == dest) {

            if ( nbrs[i].distance > distance) {
                nbrs[i].distance = distance;
                nbrs[i].spc = spc;
            } else if (nbrs[i].distance == distance) {
                nbrs[i].spc += spc;
            }
            
            updated = true;
            break;
        }
    }

    if (!updated)
        nodes[src].push_back(NodeInfo(dest,distance,spc));
    if (add_reverse)
        addEdge(dest, src, distance, spc, false);
}

void SimpleGraph::updateEdge(int src, int dest, int distance, int spc) {
    auto it = nodes.find(src);
    if (it == nodes.end()) return;
    vector<NodeInfo>& nbrs = it->second;
    for (int i = 0; i < nbrs.size();++i) {
        if (nbrs[i].nodeID == dest) {
            nbrs[i].distance = distance;
            nbrs[i].spc = spc;
            break;
        }
    }
}

map<int, NodeInfo> SimpleGraph::dijkstra(int source) {
    for (auto& pair : nodes) {
        int node = pair.first;
        NodeInfo info(node, INF, 0); // Provide initial values for nodeID, distance, and spc
        node_data[node] =  info; // Initialize node data
    }
    node_data[source].distance = 0; // Set distance of source node to 0
    node_data[source].spc = 1;

    priority_queue<pair<int, int>, vector<pair<int, int>>, std::greater<pair<int, int>>> pq;
    pq.push({0, source});
    while (!pq.empty()) {
        int cur_node = pq.top().second;
        int d = pq.top().first;
        pq.pop();
        // printf("exploreing %d\n", cur_node);
        // printf("\t%d %d\n", d, node_data[cur_node].distance);
        if (d > node_data[cur_node].distance) {
            continue;
        }
        for (const auto& edge : nodes[cur_node]) {
            int v = edge.nodeID;
            int dist = edge.distance;
            int spc = edge.spc;
            int cur_dist = node_data[cur_node].distance + dist;

            if (cur_dist < node_data[v].distance) {
                node_data[v].distance = cur_dist;
                node_data[v].spc = node_data[cur_node].spc * spc;
                pq.push({node_data[v].distance, v});
                // printf("\t to %d dist %d spc %d\n", v,node_data[v].distance,node_data[v].spc );
            } else if (cur_dist == node_data[v].distance){
                // printf("\t add spc %d dist %d spc %d\n", v,node_data[v].distance,node_data[cur_node].spc * spc);
                node_data[v].spc += node_data[cur_node].spc * spc;
            }
        }
    }

    return node_data;
}


map<int, NodeInfo> SimpleGraph::dijkstra_exclude_given_Node(int source, int exclude_node) {
    for (auto& pair : nodes) {
        int node = pair.first;
        NodeInfo info(node, INF, 0); // Provide initial values for nodeID, distance, and spc
        node_data[node] =  info; // Initialize node data
        node_data[node].pred = vector<int>();
    }
    node_data[source].distance = 0; // Set distance of source node to 0
    node_data[source].spc = 1;

    priority_queue<pair<int, int>, vector<pair<int, int>>, std::greater<pair<int, int>>> pq;
    pq.push({0, source});
    while (!pq.empty()) {
        int cur_node = pq.top().second;
        int d = pq.top().first;
        pq.pop();
        // printf("exploreing %d\n", cur_node);
        // printf("\t%d %d\n", d, node_data[cur_node].distance);
        if (d > node_data[cur_node].distance) {
            continue;
        }
        for (const auto& edge : nodes[cur_node]) {
            int v = edge.nodeID;
            if (v == exclude_node) continue; // skip exploring from a given node

            int dist = edge.distance;
            int spc = edge.spc;
            int cur_dist = node_data[cur_node].distance + dist;

            if (cur_dist < node_data[v].distance) {
                node_data[v].distance = cur_dist;
                node_data[v].spc = node_data[cur_node].spc * spc;
                
                pq.push({cur_dist, v});
                node_data[v].pred = node_data[cur_node].pred;
                node_data[v].pred.push_back(cur_node);

                // printf("\t to %d dist %d spc %d\n", v,node_data[v].distance,node_data[v].spc );
            } else if (cur_dist == node_data[v].distance){
                // printf("\t add spc %d dist %d spc %d\n", v,node_data[v].distance,node_data[cur_node].spc * spc);
                node_data[v].spc += node_data[cur_node].spc * spc;
                node_data[v].pred.push_back(cur_node);
            }
        }
    }

    return node_data;
}

/*
    for two borders, sd can achieve via
    1. the current cut 
    2. via partition (since via cut level is dijk from cut nodes)
*/   
map<pair<NodeID, NodeID>, pair<distance_t,int>> SimpleGraph::create_min_connect_graph (NodeID cut, std::vector<Neighbor> &rm_node_nbrs,
                            map<pair<NodeID, NodeID>, distance_t> &dist_via_partition,
                            map<pair<NodeID, NodeID>, distance_t> &dist_via_cut,
                            std::set<NodeID> &cut_sets,std::map<std::pair<NodeID, NodeID>, std::pair<distance_t, int>> &tmp_edges,
                            bool debug_flag = false) {
    clear();
    // printf("nbr size %d dis_p size %d dist_c size %d \n", rm_node_nbrs.size(), dist_via_partition.size(), dist_via_cut.size());

    for (auto it: rm_node_nbrs) {
        addEdge(cut, it.node, it.distance, it.edge_spc, true);
    }
    map<pair<NodeID, NodeID>, pair<distance_t,int>> new_edges;
    for (int i = 0; i < rm_node_nbrs.size();++i) {
        int n1 = rm_node_nbrs[i].node;
        if (debug_flag) {
            printf("Simple: nbr %d --\n", n1);
        }

        for (int j = i+1; j < rm_node_nbrs.size(); ++j) {
            int n2 = rm_node_nbrs[j].node;
            
            int b1 = std::min(n1,n2);
            int b2 = std::max(n1,n2);

            // sd not achieve via cut node
            if (debug_flag) {
                printf("------- inner nbr %d --\n", n2);
                printf("---------- (%d,%d) dist via cut %d  via partition %d\n",b1,b2, dist_via_cut[{b1,b2}],dist_via_partition[{b1,b2}]);
            }

            if (dist_via_cut[{b1,b2}] > dist_via_partition[{b1,b2}]) continue;

            

            int dist_direct = rm_node_nbrs[i].distance + rm_node_nbrs[j].distance;
            int spc_direct = rm_node_nbrs[i].edge_spc * rm_node_nbrs[j].edge_spc;


            //sd is achieved by cut node but not this cut node
            if (dist_direct != dist_via_cut[{b1,b2}]) continue;
            if (debug_flag) {
                printf("+++++ adding achieve via this node append to new edge   %d-%d dist %d %d spc %d\n",b1,b2,dist_direct,spc_direct);
            }
            // sd achieved by this cut node, add back 
            addEdge(b1,b2,dist_direct,spc_direct,true);
            new_edges[{b1,b2}]={dist_direct,spc_direct};

            if (cut_sets.find(b1) != cut_sets.end() 
            || cut_sets.find(b2) != cut_sets.end() ) {
                tmp_edges[{b1,b2}] = {dist_direct,spc_direct};
            }
        }
    }
    return new_edges;

}



/*
    for two borders, sd can achieve via
    1. the current cut 
    2. via partition (since via cut level is dijk from cut nodes)
*/   
map<pair<NodeID, NodeID>, pair<distance_t,int>> SimpleGraph::create_min_connect_graph (NodeID cut, std::vector<Neighbor> &rm_node_nbrs,
                            map<pair<NodeID, NodeID>, distance_t> &dist_via_partition,
                            map<pair<NodeID, NodeID>, distance_t> &dist_via_cut,
                            bool debug_flag = false) {
    clear();
    // printf("nbr size %d dis_p size %d dist_c size %d \n", rm_node_nbrs.size(), dist_via_partition.size(), dist_via_cut.size());

    for (auto it: rm_node_nbrs) {
        addEdge(cut, it.node, it.distance, it.edge_spc, true);
    }
    map<pair<NodeID, NodeID>, pair<distance_t,int>> new_edges;
    for (int i = 0; i < rm_node_nbrs.size();++i) {
        int n1 = rm_node_nbrs[i].node;
        if (debug_flag) {
            printf("Simple: nbr %d --\n", n1);
        }

        for (int j = 0; j < rm_node_nbrs.size(); ++j) {
            int n2 = rm_node_nbrs[j].node;
            
            int b1 = std::min(n1,n2);
            int b2 = std::max(n1,n2);

            // sd not achieve via cut node
            if (debug_flag) {
                printf("------- inner nbr %d --\n", n2);
                printf("---------- (%d,%d) dist via cut %d  via partition %d\n",b1,b2, dist_via_cut[{b1,b2}],dist_via_partition[{b1,b2}]);
            }

            if (dist_via_cut[{b1,b2}] > dist_via_partition[{b1,b2}]) continue;

            

            int dist_direct = rm_node_nbrs[i].distance + rm_node_nbrs[j].distance;
            int spc_direct = rm_node_nbrs[i].edge_spc * rm_node_nbrs[j].edge_spc;


            //sd is achieved by cut node but not this cut node
            if (dist_direct != dist_via_cut[{b1,b2}]) continue;
            if (debug_flag) {
                printf("----------achieve via this node append to new edge   %d-%d dist %d %d spc %d\n",b1,b2,dist_direct,spc_direct);
            }
            // sd achieved by this cut node, add back 
            addEdge(b1,b2,dist_direct,spc_direct,true);
            new_edges[{b1,b2}]={dist_direct,spc_direct};
            // printf("adding edges %d %d dist %d  spc %d\n",b1,b2,dist_direct,spc_direct);
            

        }
    }
    return new_edges;

}



// // input: the cut node that be removed, and the sd for any border pairs
// // output: the edges that need to be added to the original graph to maintain the sd and spc
// map<pair<uint32_t,uint32_t>, pair<uint32_t,uint32_t>> validate(int excluded_node,std::map<std::pair<uint32_t, uint32_t>,  uint32_t> dist_via_cuts) {
//     map<pair<uint32_t,uint32_t>, pair<uint32_t,uint32_t>> added_edges; 

//     // Iterating through the map
//     for (auto it = nodes.begin(); it != nodes.end(); ++it) {
//         int b1 = it->first;
//         if (b1 == excluded_node) continue;
//         dijkstra_exclude_given_Node(b1, excluded_node);
//         // Inner loop starting from i+1
//         for (auto it2 = std::next(it); it2 != nodes.end(); ++it2) {
//             int b2 = it2->first;
//             if (b2 == excluded_node) continue;
//             if (node_data[b2].distance == dist_via_cuts[{std::min(b1,b2), std::max(b1,b2)}]) {

//             }
            
            

//         }
//     }


// }


// int main() {
//     // Example usage
//     SimpleGraph graph;

//     /*
//         2. rm cut and create minimun graph test 
//     */
    
//     // // test 1 
//     // NodeID cut = 0;
//     // Neighbor b1(1,1,1);
//     // Neighbor b2(2,2,2);
//     // Neighbor b3(3,3,3);
//     // std::vector<Neighbor> rm_nbrs;
//     // rm_nbrs.push_back(b1);
//     // rm_nbrs.push_back(b2);
//     // rm_nbrs.push_back(b3);

//     // map<pair<NodeID, NodeID>, distance_t> dist_via_cut;
//     // dist_via_cut[{1,2}] = 3;
//     // dist_via_cut[{1,3}] = 4;
//     // dist_via_cut[{2,3}] = 5;

//     // map<pair<NodeID, NodeID>, distance_t> dist_via_partition;
//     // dist_via_partition[{1,2}] = 3;
//     // dist_via_partition[{1,3}] = 4;
//     // dist_via_partition[{2,3}] = 5;
    
//     // graph.create_min_connect_graph(cut, rm_nbrs, dist_via_partition,dist_via_cut);

//     // // Dijkstra's algorithm from node 0
//     // map<int, NodeInfo> node_data = graph.dijkstra_exclude_given_Node(3,0);

//     // // Output shortest distances from node 0 to other nodes
//     // cout << "NodeID\tDistance\tSPC\n";
//     // for (const auto& pair : node_data) {
//     //     cout << pair.second.nodeID << "\t" << pair.second.distance << "\t\t" << pair.second.spc << endl;
//     //     // printf("pred is ");
//     //     // for (auto &it_pred : pair.second.pred) {
//     //     //     printf(" %d ", it_pred);
//     //     // }
//     //     printf("\n");
//     // }

//     // test 2 
//     NodeID cut = 0;
//     Neighbor b1(1,1,1);
//     Neighbor b2(2,2,2);
//     Neighbor b3(3,3,3);
//     std::vector<Neighbor> rm_nbrs;
//     rm_nbrs.push_back(b1);
//     rm_nbrs.push_back(b2);
//     rm_nbrs.push_back(b3);

//     map<pair<NodeID, NodeID>, distance_t> dist_via_cut;
//     dist_via_cut[{1,2}] = 3;
//     dist_via_cut[{1,3}] = 4;
//     dist_via_cut[{2,3}] = 5;

//     map<pair<NodeID, NodeID>, distance_t> dist_via_partition;
//     dist_via_partition[{1,2}] = 2;
//     dist_via_partition[{1,3}] = 1;
//     dist_via_partition[{2,3}] = 5;
    
//     graph.create_min_connect_graph(cut, rm_nbrs, dist_via_partition,dist_via_cut);

//     // Dijkstra's algorithm from node 0
//     map<int, NodeInfo> node_data = graph.dijkstra_exclude_given_Node(1,0);

//     // Output shortest distances from node 0 to other nodes
//     cout << "NodeID\tDistance\tSPC\n";
//     for (const auto& pair : node_data) {
//         cout << pair.second.nodeID << "\t" << pair.second.distance << "\t\t" << pair.second.spc << endl;
//         // printf("pred is ");
//         // for (auto &it_pred : pair.second.pred) {
//         //     printf(" %d ", it_pred);
//         // }
//         printf("\n");
//     }

    




//     // /*
//     //     1. basic add edge and dijk test
//     // */
//     // // Add edges
//     // graph.addEdge(1,3,1,2,true);
//     // graph.addEdge(1,5,3,4,true);
//     // graph.addEdge(1,5,1,2,true);
//     // graph.addEdge(1,5,1,1,true);
//     // graph.addEdge(1,7,1,2,true);
//     // graph.addEdge(5,9,1,1,true);
//     // graph.addEdge(5,9,1,2,true);
//     // graph.addEdge(7,9,1,1,true);
//     // graph.addEdge(9,10,2,2,true);

//     // // Dijkstra's algorithm from node 0
//     // map<int, NodeInfo> node_data = graph.dijkstra_exclude_given_Node(1,0);

//     // // Output shortest distances from node 0 to other nodes
//     //  cout << "NodeID\tDistance\tSPC\n";
//     // for (const auto& pair : node_data) {
//     //     cout << pair.second.nodeID << "\t" << pair.second.distance << "\t\t" << pair.second.spc << endl;
//     //     printf("pred is ");
//     //     for (auto &it_pred : pair.second.pred) {
//     //         printf(" %d ", it_pred);
//     //     }
//     //     printf("\n");
//     // }

//     return 0;
// }
