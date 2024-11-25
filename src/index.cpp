#include "road_network.h"
#include "util.h"
#include <sstream>
#include <iostream>
#include <cstring>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace road_network;

int main(int argc, char** argv)
{
    const char* graph_dir = argv[1];
    std::string graph_dir_str(graph_dir);
    std::string dataset_dir = graph_dir_str.substr(0, graph_dir_str.find("d."));

    const char* index_dir = argv[2];
    std::string index_dir_str(index_dir);
    std::ofstream outputFile(index_dir_str);

    std::stringstream filenameStream;
    filenameStream << graph_dir;
    std::string graphFileName = filenameStream.str();
    fstream fs(graphFileName);
    Graph g;
    read_graph(g, fs);
    vector<Neighbor> closest;
    g.contract(closest);

    // construct index
    vector<CutIndex> ci;

    util::start_timer();
    g.create_cut_index(ci, 0.2);

    for (auto it : external_test_set) {
        printf("ci[%d] level %u\n", it, ci[it].cut_level);
    }
    
    ContractionIndex con_index(ci, closest);
    double indexing_time = util::stop_timer();
    cout << "build " << graph_dir_str << " index  in " << indexing_time << endl;
    // write index

    util::start_timer();
    if (outputFile.is_open()) {
        con_index.write(outputFile);
        outputFile.close();
    } else {
        std::cerr << "Error: Unable to open file for writing.\n";
    }
    double io_time = util::stop_timer();
    cout << "finish write out " << graph_dir_str << " index  in " << io_time << endl;


    return 0;
}
