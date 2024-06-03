#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <chrono>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <memory>

extern "C" {
    #include "nauty.h"
}

//// Global Variables
bool TRIANGLE_FREE = true;
bool DENSE_SET_FREE = true;

int k_dense = 2;
int dense_set_size = 6;
int k_sparse = 2;
int sparse_set_size = 8;
bool is_cout_extremals = false;

static optionblk options;
static statsblk stats;

int calculate_degree_sum(std::vector<std::vector<bool>> &graph){
    
    std::vector<int> degree_seq;
    for(const auto& row: graph){
        int temp_degree = 0;
        for(const auto& entry: row){
            if(entry){
                temp_degree++;
            }
        }
        degree_seq.push_back(temp_degree);
    }

    int sum = 0;
    for(const auto& degree: degree_seq){
        sum += degree;
    }
    return sum;
}

void print_graph(std::vector<std::vector<bool>>& graph){
    for(const auto& row: graph){
        for(const auto& entry: row){
            std::cout << entry << ",";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    return;
}

graph* create_nauty_graph_from_adj_matrix(const std::vector<std::vector<bool>>& adj_matrix) {
    int n = adj_matrix.size();

    graph* g = new graph[n * MAXM];
    EMPTYGRAPH(g, MAXM, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (adj_matrix[i][j]) {
                ADDONEEDGE(g, i, j, MAXM);
            }
        }
    }

    return g;
}

std::vector<std::vector<int>> combinations(int n, int r) {
    if (n < r) {
        return {};
    }

    std::vector<std::vector<int>> comb_list;
    std::vector<bool> bitmask(n, 0);
    for (int i = n-r; i < n; ++i) {
        bitmask[i] = 1;
    }

    do {
        std::vector<int> current_combination;
        for (int i = 0; i < n; ++i) {
            if (bitmask[i]) {
                current_combination.push_back(i);
            }
        }
        comb_list.push_back(current_combination);
    } while (std::next_permutation(bitmask.begin(), bitmask.end()));

    return comb_list;
}

bool include_k_sparse_x_set(std::vector<std::vector<bool>> &graph, std::vector<std::vector<int>> &combs){

    if(combs.size() == 0) return false;
    for(const auto& comb : combs){
        bool is_set_k_sparse = true;
        for(int i = sparse_set_size - 1; i >= 0; i--){
            int sub_degree_of_vert_i = 0;
            for(int j = 0; j < sparse_set_size; j++){
                if(i == j) continue;
                if(graph[comb[i]][comb[j]]){
                    sub_degree_of_vert_i++;
                }
            }
            if(sub_degree_of_vert_i >= (k_sparse + 1)){
                is_set_k_sparse = false;
                break;
            }
        }
        if(is_set_k_sparse) return true;
    }
    return false;
}

bool include_j_dense_x_set(std::vector<std::vector<bool>> &graph, std::vector<std::vector<int>> &combs){

    if(combs.size() == 0) return false;
    for(const auto& comb : combs){
        bool is_set_j_dense = true;
        for(int i = dense_set_size - 1; i >= 0; i--){
            int sub_degree_of_vert_i = 0;
            for(int j = 0; j < dense_set_size; j++){
                if(i == j) continue;
                if(graph[comb[i]][comb[j]]){
                    sub_degree_of_vert_i++;
                }
            }
            if(sub_degree_of_vert_i <= (dense_set_size - 2 - k_dense)){
                is_set_j_dense = false;
                break;
            }
        }
        if(is_set_j_dense) return true;
    }
    return false;
}

std::vector<std::vector<std::vector<bool>>> add_vertex_to_graph(std::vector<std::vector<bool>> prev_graph, 
    std::vector<std::vector<int>>& combination_sets_sparse, std::vector<std::vector<int>>& combination_sets_dense){
    
    std::vector<std::vector<std::vector<bool>>> new_graphs;
    std::vector<graph*> canonical_labelings;
    int n = prev_graph.size();
    
    // Loop over 2^n possible combinations (subsets)
    for (int subset_int = 0; subset_int < (1 << n); subset_int++) {
        std::vector<int> subset;

        // Loop over bits of i
        for (int j = 0; j < n; j++) {
            // Check if the jth bit in i is set
            if (subset_int & (1 << j)) {
                subset.push_back(j);
            }
        }

        // Eliminate Triangles. Subsets must be independent sets
        if(TRIANGLE_FREE){
            bool is_ind_set = true;
            for(int i = 0; i < subset.size(); i++){
                for(int j = i + 1; j < subset.size(); j++){
                    if(prev_graph[subset[i]][subset[j]]){
                        is_ind_set = false;
                        break;
                    }
                }
                if(!is_ind_set) break;
            }

            // If not ind. set -> new graph will have triangles -> pass this subset
            if(!is_ind_set){
                continue;
            }
        }

        // Create the new graph
        std::vector<std::vector<bool>> new_graph;
        for(int i = 0; i < n; i++){
            std::vector<bool> new_row = prev_graph[i];
            // if the row (vertex i) is in subset new row will have 1 added to end. Otherwise 0
            if(std::find(subset.begin(), subset.end(), i) != subset.end()){
                new_row.push_back(1);
            }
            else{
                new_row.push_back(0);
            }
            new_graph.push_back(new_row);
        }

        // Add last row which the new vertex.
        std::vector<bool> new_row (n + 1, false);
        for(const auto& i : subset){
            new_row[i] = 1;
        }
        new_graph.push_back(new_row);

        // Check if there is any k-sparse _ set.
        // We assumed the input graph had not such sets.
        // So if there is in the new graph, such a set will have to include the new vertex.
        if(include_k_sparse_x_set(new_graph, combination_sets_sparse)) continue;
        
        if(DENSE_SET_FREE){
            //or j-dense set -> pass this subset
            bool is_there_j_dense_x_set = include_j_dense_x_set(new_graph, combination_sets_dense);
            if(is_there_j_dense_x_set) continue;
        }

        // calculate canonical labeling
        graph* g1 = create_nauty_graph_from_adj_matrix(new_graph);
        graph* new_graph_canonical_labeling = new graph[new_graph.size() * MAXM];
        int lab1[MAXN], ptn[MAXN], orbits[MAXN];
        densenauty(g1, lab1, ptn, orbits, &options, &stats, MAXM, new_graph.size(), new_graph_canonical_labeling);
        delete[] g1;

        bool is_isomorphically_new = true;
        for(const auto& existing_graph_canonical_label: canonical_labelings){
            bool isomorph = true;
            for (int i = 0; i < new_graph.size() * MAXM; ++i) {
                if (new_graph_canonical_labeling[i] != existing_graph_canonical_label[i]) {
                    isomorph = false;
                    break;
                }
            }
            if(isomorph){
                is_isomorphically_new = false;
                break;
            }
        }

        if(is_isomorphically_new){
            new_graphs.push_back(new_graph);
            canonical_labelings.push_back(new_graph_canonical_labeling);
        }
    }

    for (int i = 0; i < canonical_labelings.size(); ++i) {
        delete[] canonical_labelings[i];
    }

    return new_graphs;
}

static std::vector<int> degree_sums;

void isomorphism_checks(int thread_id, std::mutex &mutex, int& shared_index, int& end_index,
    std::vector<std::vector<std::vector<bool>>> &new_graphs, std::vector<bool> &labels, std::vector<std::shared_mutex> &add_label_mutexes, std::mutex &mutex2, std::vector<graph*>& canonical_labelings){
    
    while(true){
        
        int local_index;
        {
            std::lock_guard<std::mutex> lock(mutex);
            if (shared_index >= end_index) {
                return;
            }
            local_index = shared_index;
            shared_index++;
        }
        // iso check
        bool is_isomorphically_new = true;
        for(int j = local_index - 1; j >= 0; j--){
            if(degree_sums[local_index] != degree_sums[j]){
                break;
            }
            {
                add_label_mutexes[j].lock_shared();
                if(!labels[j]) continue;
                add_label_mutexes[j].unlock_shared();
            }
            // canon label check
            bool isomorph = true;
            for (int i = 0; i < new_graphs[local_index].size() * MAXM; ++i) {
                if (canonical_labelings[local_index][i] != canonical_labelings[j][i]) {
                    isomorph = false;
                    break;
                }
            }
            if(isomorph){
                is_isomorphically_new = false;
                break;
            }
        }
        if(!is_isomorphically_new) {
            std::lock_guard<std::shared_mutex> lock1(add_label_mutexes[local_index]);
            std::lock_guard<std::mutex> lock2(mutex2);
            labels[local_index] = false;
        }
        if(local_index % 10000 == 0 && local_index > 0){
            std::cout << "Isomorphism Phase: Graphs Processed / All To Be Processed: " << local_index << " / " << end_index << std::endl;
        }
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }
}

std::vector<std::vector<std::vector<bool>>> find_graphs_of_one_more_cardinality(std::vector<std::vector<std::vector<bool>>> prev_graphs){

    int thread_count = std::thread::hardware_concurrency();
    
    //////////////////////////////////
    //       Generation Phase       //
    //////////////////////////////////

    std::vector<std::vector<std::vector<bool>>> new_graphs;
    int no_of_prev_graphs = prev_graphs.size();
    new_graphs.reserve(no_of_prev_graphs * 1000); // RESERVE, HARDCODED

    int n = prev_graphs[0].size(); // prev graph size

    // Combinations for check k-sparse _ sets.
    // If there is any new k-sparse set, it has to include the new vertex.
    std::vector<std::vector<int>> combination_sets_sparse = combinations(n, sparse_set_size - 1);
    for(int i = 0; i < combination_sets_sparse.size(); i++){
        combination_sets_sparse[i].push_back(n);
    }
    std::vector<std::vector<int>> combination_sets_dense = {};
    if(DENSE_SET_FREE){
        combination_sets_dense = combinations(n, dense_set_size - 1);
        for(int i = 0; i < combination_sets_dense.size(); i++){
            combination_sets_dense[i].push_back(n);
        }
    }

    std::vector<std::thread> worker_threads;
    std::mutex m1; // for getting the shared index
    std::mutex m2; // for adding new graphs to the new_graphs vector

    int shared_index = 0;

    for(int i = 0; i < thread_count; i++){
        worker_threads.push_back(std::thread([](std::mutex &mutex1, std::mutex &mutex2, int& shared_index, int& end_index, std::vector<std::vector<std::vector<bool>>> &prev_graphs,
         std::vector<std::vector<std::vector<bool>>> &new_graphs, std::vector<std::vector<int>>& combination_sets_sparse, std::vector<std::vector<int>>& combination_sets_dense){
            while(true){
                int local_index;
                {
                    std::lock_guard<std::mutex> lock(mutex1);
                    if (shared_index >= end_index) {
                        return;
                    }
                    local_index = shared_index;
                    shared_index++;
                }
                const auto& temp_results = add_vertex_to_graph(std::ref(prev_graphs[local_index]), std::ref(combination_sets_sparse), std::ref(combination_sets_dense));
                {
                    std::lock_guard<std::mutex> lock(mutex2);
                    for(const auto& graph: temp_results){
                        new_graphs.push_back(graph);
                    }
                }
                if(local_index % 1000 == 0 && local_index > 0){
                    std::cout << "Generation Phase: Graphs Processed / All To Be Processed: " << local_index << " / " << end_index << std::endl;
                }
            }
        }, std::ref(m1), std::ref(m2), std::ref(shared_index), std::ref(no_of_prev_graphs), std::ref(prev_graphs), std::ref(new_graphs), std::ref(combination_sets_sparse), std::ref(combination_sets_dense)));
    }

    for(auto& thread: worker_threads){
        thread.join();
    }
    worker_threads.clear(); // clear threads for isomorphism phase

    //////////////////////////////////
    //      Isomorphism Checks      //
    //////////////////////////////////

    std::vector<std::vector<std::vector<bool>>> non_iso_new_graphs;

    // sort graphs by their degree sums
    std::sort(new_graphs.begin(), new_graphs.end(), [](std::vector<std::vector<bool>> &graph1, std::vector<std::vector<bool>> &graph2)
                                                    {
                                                        return (calculate_degree_sum(graph1) > calculate_degree_sum(graph2));
                                                    });

    int new_graphs_size = new_graphs.size();

    // Label vector initialized as true. If false, graph is isomorphic to another graph and will be labeled false as to be removed
    std::vector<bool> add_label(new_graphs_size, true);

    degree_sums.clear();
    for(int i = 0; i < new_graphs_size; i++){
        degree_sums.push_back(calculate_degree_sum(new_graphs[i]));
    }
    std::vector<std::shared_mutex> add_label_mutexes(new_graphs_size);
    
    // Create canonical labelings for each graph
    std::vector<graph*> canonical_labelings(new_graphs_size); // A Vector to hold nauty canonical labelings, so that they are calculated once
    for(int i = 0; i < new_graphs_size; i++){

        graph* g1 = create_nauty_graph_from_adj_matrix(new_graphs[i]);
        canonical_labelings[i] = new graph[new_graphs[i].size() * MAXM];
        int lab1[MAXN], ptn[MAXN], orbits[MAXN];
        densenauty(g1, lab1, ptn, orbits, &options, &stats, MAXM, new_graphs[i].size(), canonical_labelings[i]);
        delete[] g1;
    }

    shared_index = 0;
    int ending_index = new_graphs_size;

    std::mutex m3; // for getting the shared index
    std::mutex m4; // for updating the add_label vector

    for(int i = 0; i < thread_count; i++){
        worker_threads.push_back(std::thread(isomorphism_checks, i, std::ref(m3), std::ref(shared_index), std::ref(ending_index),
        std::ref(new_graphs), std::ref(add_label), std::ref(add_label_mutexes), std::ref(m4), std::ref(canonical_labelings)));
    }

    for(auto& thread: worker_threads){
        thread.join();
    }

    for (int i = 0; i < new_graphs_size; ++i) {
        delete[] canonical_labelings[i];
    }

    int no_of_uniques = 0;
    for(int i = 0; i < add_label.size(); i++){
        if(add_label[i]) no_of_uniques++;
    }
    non_iso_new_graphs.reserve(no_of_uniques);
    for(int i = 0; i < add_label.size(); i++){
        if(add_label[i]){
            non_iso_new_graphs.push_back(new_graphs[i]);
        }
    }

    return non_iso_new_graphs;
}

int main(){

    auto time_start = std::chrono::high_resolution_clock::now();

    // nauty objects
    options = {0,FALSE,FALSE,FALSE,TRUE,FALSE,CONSOLWIDTH, NULL,NULL,NULL,NULL,NULL,NULL,NULL,100,0,1,0,&dispatch_graph,FALSE,NULL};
    options.getcanon = TRUE;  // Important to get canonical labeling

    std::vector<std::vector<std::vector<bool>>> one_vertex_graph_list = {{{0}}};
    std::vector<std::vector<std::vector<std::vector<bool>>>> graphs_vector;
    graphs_vector.push_back(one_vertex_graph_list);
    int last_index = 0;

    int most_cardinality = 64; // Max vertex number program will go for

    for(int i = 1; i < most_cardinality; i++){
        graphs_vector.push_back(find_graphs_of_one_more_cardinality(graphs_vector[i-1]));
        std::cout << std::endl;
        std::cout << "Number of size " << i + 1 << " graphs: " << graphs_vector[i].size() << std::endl;

        if(graphs_vector[i].size() == 0){
            last_index = i - 1;
            break;
        }
    }

    if(is_cout_extremals){
        for(auto graph: graphs_vector[last_index]){
            print_graph(graph);
        }
    }

    auto time_end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    std::cout << "Time: " << elapsed << " milliseconds" << std::endl;

    return 0;
}