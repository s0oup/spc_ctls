#include "road_network.h"
#include "util.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <thread>
#include <sys/resource.h>

#include <string>

using namespace std;
using namespace road_network;

// #define REMOVE_REDUNDANT
#define CONTRACT

const size_t repeats = 5;
const size_t nr_query_tests = 10;
const size_t nr_buckets = 10;
const size_t bucket_size = 10000;
const distance_t bucket_min = 1000;

const size_t MB = 1024 * 1024;


// generate distance-based query: 
// ./generate_bucket_query dataset/NY/d.NY.gr dataset/NY/idx_NY dataset/NY/NY_bucket_query 5

// diff dataset/NY/ans.NY.bucket_0 dataset/NY/NY_bucket_query_0
// diff dataset/NY/ans.NY.bucket_1 dataset/NY/NY_bucket_query_1
// diff dataset/NY/ans.NY.bucket_2 dataset/NY/NY_bucket_query_2
// diff dataset/NY/ans.NY.bucket_3 dataset/NY/NY_bucket_query_3
// diff dataset/NY/ans.NY.bucket_4 dataset/NY/NY_bucket_query_4

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        cout << "syntax: " << argv[0] << " <graph_dir> <index_dir> <query_file_format> <num_of_files>" << endl;
        return 0;
    }
    const char* graph_name = argv[1];
    std::stringstream graph_stream;
    graph_stream << graph_name;
    std::string graph_name_dir = graph_stream.str();
    cout << endl << "reading graph from " << graph_name_dir << endl;
    fstream fs(graph_name_dir);
    Graph g;
    read_graph(g, fs);
    fs.close();
    cout << "read " << g.node_count() << " vertices and " << g.edge_count() << " edges" << flush;
    distance_t diameter = g.diameter(true);
    cout << " (diameter=" << g.diameter(false) << "|" << diameter << ")" << endl;

    const char* index_dir = argv[2];
    std::string index_dir_str(index_dir);
    fstream index_stream(index_dir_str);
    
    std::string queryFileFormat = argv[3];
    size_t num_of_files = std::stoi(argv[4]);
    queryFileFormat += "_";


    printf("start load index ");
    cout<< index_dir_str << std::endl;

    // read index
    util::start_timer();
    ContractionIndex con_index(index_stream);
    double read_index_time = util::stop_timer();
    cout << "read index in " << read_index_time << "s (" << con_index.size() / MB << " MB)" << endl;
    

    std::stringstream query_dir; 
    for (size_t i = 0; i < num_of_files; ++ i) {
        query_dir.str("");
        query_dir << queryFileFormat << i;

        std::cerr << "output file name is " << query_dir.str() << " ";
        std::ofstream outputFile(query_dir.str());

        cout << "generating queries by distance: " << flush;
        vector<vector<pair<NodeID,NodeID>>> query_buckets(nr_buckets);
        util::start_timer();
        g.random_pairs(query_buckets, bucket_min, bucket_size, con_index);
        cout << " in " << util::stop_timer() << "s" << endl;
        spc_distance_t result_arr[nr_buckets][bucket_size] = {{0}};
        for (size_t bucket = 0; bucket < query_buckets.size(); bucket++)
        {
            util::start_timer();
            int tmp_c=0;
            for (pair<NodeID,NodeID> q : query_buckets[bucket]){
                // con_index.get_distance(q.first, q.second);
                result_arr[bucket][tmp_c++] = con_index.get_distance(q.first, q.second);
            }

            // wirte out for this bucket 
            if (outputFile.is_open()) {          
                for (int count =0; count < bucket_size; ++count) {
                    outputFile << query_buckets[bucket][count].first << " " << query_buckets[bucket][count].second 
                        << " " << static_cast<distance_t>(result_arr[bucket][count]) << " " << static_cast<spc_t>(result_arr[bucket][count]>> DIST_SHIFT) << "\n";
                }

            }else {
                std::cerr << "==========Unable to open file for writing.\n";
            }    
        }


        outputFile.close();
        
        }



        

   

    return 0;
}
