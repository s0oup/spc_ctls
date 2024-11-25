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
const size_t nr_queries = 1000000;
const size_t nr_query_tests = 10;
const size_t nr_buckets = 10;
const size_t bucket_size = 10000;
const distance_t bucket_min = 1000;

const size_t MB = 1024 * 1024;

struct ResultData
{
    size_t label_count;
    size_t max_label_count;
    size_t index_size;
    size_t index_height;
    double index_time;
    double avg_cut_size;
    size_t max_cut_size;
    size_t pruning_2hop;
    size_t pruning_3hop;
    size_t pruning_tail;
    double random_query_time;
    double random_hoplinks;
    vector<double> bucket_query_times;
    vector<double> bucket_hoplinks;
};

struct FileResults
{
    string filename;
    vector<ResultData> results;
    FileResults(string filename, vector<ResultData> results) : filename(filename), results(results) {}
};

ostream& operator<<(ostream& os, const FileResults& fr)
{
    if (fr.results.empty())
        return os;
    os << endl << "Summary for " << fr.filename << ":" << endl << setprecision(5);
    os << "Index size (MB): " << util::summarize(fr.results, [](ResultData r) -> double { return r.index_size; }) << endl;
    os << "Index time (s): " << util::summarize(fr.results, [](ResultData r) -> double { return r.index_time; }) << endl;
    os << "Index height: " << util::summarize(fr.results, [](ResultData r) -> double { return r.index_height; }) << endl;
    os << "Avg cut size: " << util::summarize(fr.results, [](ResultData r) -> double { return r.avg_cut_size; }) << endl;
    os << "Max cut size: " << util::summarize(fr.results, [](ResultData r) -> double { return r.max_cut_size; }) << endl;
    os << "Query time (s): " << util::summarize(fr.results, [](ResultData r) -> double { return r.random_query_time; }) << endl;
    os << "Avg Hoplinks: " << util::summarize(fr.results, [](ResultData r) -> double { return r.random_hoplinks; }) << endl;
    if (!fr.results[0].bucket_query_times.empty())
        for (size_t bucket = 0; bucket < nr_buckets; bucket++)
        {
            os << "Bucket " << bucket << ": time = " << util::summarize(fr.results, [bucket](ResultData r) -> double { return r.bucket_query_times[bucket]; }) * (nr_queries / bucket_size);
            os << ", hoplinks = " << util::summarize(fr.results, [bucket](ResultData r) -> double { return r.bucket_hoplinks[bucket]; }) << endl;
        }
    return os;
}

#ifdef PRUNING
size_t get_2hop_pruning(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (NodeID node = 1; node < ci.size(); node++)
        total += ci[node].pruning_2hop;
    return total;
}

size_t get_3hop_pruning(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (NodeID node = 1; node < ci.size(); node++)
        total += ci[node].pruning_3hop;
    return total;
}

size_t get_tail_pruning(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (NodeID node = 1; node < ci.size(); node++)
        total += ci[node].pruning_tail;
    return total;
}
#endif

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "syntax: " << argv[0] << " [balance] <filename> ... <filename>" << endl;
        return 0;
    }
    // check for balance parameter
    double balance = atof(argv[1]);
    int file_start = 2;
    if (balance == 0.0)
    {
        balance = 0.2;
        file_start = 1;
    }

#ifdef NO_SHORTCUTS
    cout << "shortcuts disabled" << endl;
#else
    cout << "shortcuts enabled" << endl;
#endif

#ifdef PRUNING
    cout << "pruning enabled" << endl;
#else
    cout << "pruning disabled" << endl;
#endif

#ifdef CONTRACT2D
    cout << "path contraction enabled" << endl;
#else
    cout << "path contraction disabled" << endl;
#endif

#ifdef MULTI_THREAD
    cout << "multi-threading enabled" << endl;
    cout << "threads supported by hardware: " << thread::hardware_concurrency() << endl;
#else
    cout << "multi-threading disabled" << endl;
#endif

#ifdef NDEBUG
    srand(time(nullptr));
#endif
std::string directory = "/home/soup/spc-shortes_path_counting/road-networks-main/";

    vector<FileResults> file_results;
    for (int f = file_start; f < argc; f++)
    {
        const char* filename = argv[f];
        vector<ResultData> results;
        bool use_buckets;
/////////////////////////////////////////////
        /*
            write n times and write all result out
            now only build graph once to generate queries 

        */
        std::stringstream filenameStream;
        filenameStream << directory << filename;
        std::string graphFileName = filenameStream.str();
        cout << endl << "reading graph from " << graphFileName << endl;
        fstream fs(graphFileName);
        Graph g;
        read_graph(g, fs);
        fs.close();
        cout << "read " << g.node_count() << " vertices and " << g.edge_count() << " edges" << flush;
        distance_t diameter = g.diameter(true);
        cout << " (diameter=" << g.diameter(false) << "|" << diameter << ")" << endl;
        // check for redundant edges
        vector<Edge> redundant_edges;
        util::start_timer();
        
#ifdef REMOVE_REDUNDANT
        g.get_redundant_edges(redundant_edges);
        for (Edge e : redundant_edges)
            g.remove_edge(e.a, e.b);
        cout << "removed " << redundant_edges.size() << " redundant edges in " << util::stop_timer() << "s" << endl;
            cout << "found " << redundant_edges.size() << " redundant edges in " << util::stop_timer() << "s" << endl;
#else
        
#endif
#ifdef CONTRACT
        util::start_timer();
        size_t old_size = g.node_count();
        vector<Neighbor> closest;
        g.contract(closest);
        cout << "contracted to " << g.node_count() << " vertices (" << g.node_count() * 100 / max<size_t>(1, old_size) << "%) and "
            << g.edge_count() << " edges in " << util::stop_timer() << "s" << endl;
#endif
#ifdef CONTRACT2D
        size_t deg2nodes = 0;
        for (NodeID node : g.get_nodes())
            if (g.degree(node) == 2)
                deg2nodes++;
        cout << deg2nodes << " of these vertices (" << deg2nodes * 100 / max<size_t>(1, g.node_count()) << "%) have degree 2" << endl;
#endif
#ifdef NDEBUG
        g.randomize();
#endif
        ResultData result = {};
        // construct index
        Graph::show_progress(true);
        vector<CutIndex> ci;
        util::start_timer();
        size_t shortcuts = g.create_cut_index(ci, balance);
#ifdef PRUNING
        result.pruning_2hop = get_2hop_pruning(ci);
        result.pruning_3hop = get_3hop_pruning(ci);
        result.pruning_tail = get_tail_pruning(ci);
#endif
#ifdef CONTRACT
        ContractionIndex con_index(ci, closest);
#else
        ContractionIndex con_index(ci);
#endif
        result.index_time = util::stop_timer();
        result.index_size = con_index.size() / MB;
        result.label_count = con_index.label_count();
        result.max_label_count = con_index.max_label_count();
        result.index_height = con_index.height();
        result.avg_cut_size = con_index.avg_cut_size();
        result.max_cut_size = con_index.max_cut_size();
        cout << "created index of size " << result.index_size << " MB in " << result.index_time << "s using " << shortcuts << " shortcuts" << endl;
        cout << "#labels=" << result.label_count << " (max " << result.max_label_count << ")" << ", avg/max cut size=" << setprecision(3) << result.avg_cut_size << "/" << result.max_cut_size << ", height=" << result.index_height << endl;
        cout << "partition tree contains " << con_index.non_empty_cuts() << " non-empty cuts (" << 100 * con_index.non_empty_cuts() / con_index.uncontracted_count() << "% of uncontracted vertices)" << endl;
#ifdef PRUNING
        size_t unpruned_labels = max<size_t>(1, result.label_count + result.pruning_tail);
        cout << "3-HOP pruning could remove " << result.pruning_3hop << " labels (" << result.pruning_3hop * 100 / unpruned_labels << "%)" << endl;
        cout << "2-HOP pruning could remove " << result.pruning_2hop << " labels (" << result.pruning_2hop * 100 / unpruned_labels << "%)" << endl;
        cout << "tail pruning *has* removed " << result.pruning_tail << " labels (" << result.pruning_tail * 100 / unpruned_labels << "%)" << endl;
#endif
        g.reset(); // needed for distance testing

        // show memory consumption
        rusage usage;
        if (getrusage(RUSAGE_SELF, &usage) != -1)
            cout << "maximum memory used: " << usage.ru_maxrss / 1024 << " MB" << endl;
                    
/////////////////////////////////////////////////////////


        for (size_t i = 0; i < repeats; i++)
        {
            filenameStream.str("");

            filenameStream << graphFileName << "_query_" << i;
            std::string resultFileName = filenameStream.str();
            // size_t pos = resultFileName.find(directory);
            // if (pos != std::string::npos) {
            //     // Erase the prefix from the fullPath
            //     resultFileName.erase(pos, directory.length());
            // }
            std::cerr << "output file name is " << resultFileName;
            std::ofstream outputFile(resultFileName);

            use_buckets = diameter >= bucket_min * nr_buckets;
distance_t result_arr[nr_buckets][bucket_size] = {{0}};
            if (use_buckets)
            {
                cout << "generating queries by distance: " << flush;
                vector<vector<pair<NodeID,NodeID>>> query_buckets(nr_buckets);
                util::start_timer();
                g.random_pairs(query_buckets, bucket_min, bucket_size, con_index);
                cout << " in " << util::stop_timer() << "s" << endl;
                for (size_t bucket = 0; bucket < query_buckets.size(); bucket++)
                {
                    util::start_timer();
int tmp_c=0;
                    for (pair<NodeID,NodeID> q : query_buckets[bucket]){
                        // con_index.get_distance(q.first, q.second);
result_arr[bucket][tmp_c++]=con_index.get_distance(q.first, q.second);
                    }
                    result.bucket_query_times.push_back(util::stop_timer());
                    result.bucket_hoplinks.push_back(con_index.avg_hoplinks(query_buckets[bucket]));
                    cout << "ran " << query_buckets[bucket].size() << " queries (bucket " << bucket << ") in " << result.bucket_query_times.back() << "s (hoplinks=" << result.bucket_hoplinks.back() << ")" << endl;


if (outputFile.is_open()) {          
    for (int count =0; count < bucket_size; ++count) {
        outputFile << query_buckets[bucket][count].first << " " << query_buckets[bucket][count].second << " " << result_arr[bucket][count] << "\n";
    }

}else {
    std::cerr << "==========Unable to open file for writing.\n";
}                
                }
outputFile.close();
            }
            results.push_back(result);
        }
        if (repeats > 1)
            file_results.push_back(FileResults(filename, results));
    }
    for (FileResults& fr : file_results)
        cout << fr;
    return 0;
}
