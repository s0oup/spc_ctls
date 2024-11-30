#include "road_network.h"
#include "util.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <chrono>
#include <thread>
using namespace std;
using namespace road_network;

const size_t nr_buckets = 10;
const size_t bucket_size = 100000;
const size_t MB = 1024 * 1024;

void load_query(vector<pair<NodeID,NodeID>>& queries, std::ifstream& queryFile)
{
    NodeID u, v;
    distance_t d;
    spc_t spc;
    while (queryFile >> u >> v >> d >> spc) {
        queries.push_back(make_pair(u, v));
    }
  
}

int main(int argc, char* argv[])
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <index_path> <queries_format> <result_format> n" << std::endl;
        return 1; // Return an error code
    }

    const char* index_dir = argv[1];
    std::string index_dir_str(index_dir);
    fstream index_stream(index_dir_str);
    
    std::string queryFileFormat = argv[2];
    std::string resultFileFormat = argv[3];
    size_t num_of_files = std::stoi(argv[4]);
    queryFileFormat += "_";
    resultFileFormat += "_";

    printf("start load index ");
    cout<< index_dir_str << std::endl;

    // read index
    util::start_timer();
    ContractionIndex con_index(index_stream);
    double read_index_time = util::stop_timer();
    cout << "read index in " << read_index_time << "s (" << con_index.size() / MB << " MB)" << endl;
    cout << "index tree height " << con_index.height() << endl;
    cout << "ave cut size " << con_index.avg_cut_size() << endl;
    cout << "max cut size " << con_index.max_cut_size() << endl;
    cout << "max CA(largest label size) " << con_index.max_label_count() << endl;
    cout << endl;




    for (size_t i = 0; i < num_of_files; ++i) {
        std::string queryFilePath =  queryFileFormat + std::to_string(i);
        std::ifstream queryFile(queryFilePath);
        if (!queryFile.is_open()) {
            std::cerr << "Error opening queries file: " << queryFilePath << std::endl;
            return 1; // Return an error code
        }
        cout << "processing query file " << queryFilePath << endl;

    
        std::string resultFilePath = resultFileFormat + std::to_string(i);
        std::ofstream resultFile(resultFilePath);
        if (!resultFile.is_open()) {
            std::cerr << "Error opening result file: " << queryFilePath << std::endl;
            return 1; // Return an error code
        }

        vector<pair<NodeID,NodeID>> queries;
        load_query(queries,queryFile);
        queryFile.close();


        vector<double> bucket_time_recorder(10, 0.0);
        std::vector<spc_distance_t> spc_distances(queries.size()); 

        auto overall_start = std::chrono::high_resolution_clock::now();
        for (size_t bucket_c = 0; bucket_c < nr_buckets; ++bucket_c ) {
            util::start_timer();
            for (size_t query_c = 0; query_c < bucket_size; ++query_c) {
                spc_distances[bucket_c*bucket_size+query_c] = con_index.get_distance(queries[bucket_c*bucket_size+query_c].first, queries[bucket_c*bucket_size+query_c].second);
            }
            bucket_time_recorder[bucket_c] = util::stop_timer();
        }
        auto overall_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = overall_end - overall_start;
        double total_hoplinks = 0;

        cout<< "run time for each bucket (in sec). Bucket size: " << bucket_size << endl;
        for (size_t bucket_c = 0; bucket_c < nr_buckets; ++bucket_c) {
            cout << bucket_time_recorder[bucket_c] << " ";
        }
        cout << endl;

        cout << "ave hoplinks for each bucket: \n";
        for (size_t bucket_c = 0; bucket_c < nr_buckets; ++bucket_c) {
            std::vector<std::pair<NodeID, NodeID>> subset(queries.begin() + bucket_c*bucket_size, queries.begin() + (bucket_c+1)*bucket_size);
            double random_hoplinks = con_index.avg_hoplinks(subset);
            total_hoplinks += random_hoplinks;
            cout<< random_hoplinks << " ";
        }
        cout << endl;

        cout << "Overall ave query time: " << elapsed.count() << " s \n";
        cout<< "Ave query time (in micro sec): " << (std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count())/( (double)(nr_buckets*bucket_size))<< " mu s\n"; 
        cout << "Ave hoplinks: " << total_hoplinks/nr_buckets << "\n\n";
        

        for (size_t j = 0; j < nr_buckets*bucket_size; ++j) {
            resultFile << queries[j].first << " " << queries[j].second << " " << static_cast<distance_t> (spc_distances[j]) << " " <<static_cast<spc_t>(spc_distances[j]>>32)<< "\n";
        }    
    }
    
    return 0;
}
