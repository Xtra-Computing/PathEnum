//
// Created by Shixuan Sun on 2020/6/30.
//

#include <queue>
#include <unordered_set>
#include <random>
#include <iterator>
#include <fstream>
#include <chrono>
#include <iostream>

#include "util/log/log.h"
#include "util/graph/directed_graph.h"
#include "util/graph/graph_operation.h"

typedef std::pair<uint32_t, uint32_t> pi;

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

void generate_pairs(uint32_t num_pairs, std::vector<uint32_t>& src_array, std::vector<uint32_t>& dst_array,
        std::vector<pi>& pair_array) {
    auto start = std::chrono::high_resolution_clock::now();
    uint32_t step_length = num_pairs * 0.2;

    for (uint32_t i = 0; i < num_pairs; ++i) {
        uint32_t src, dst;
        do {
            src = *select_randomly(src_array.begin(), src_array.end());
            dst = *select_randomly(dst_array.begin(), dst_array.end());
        } while (src == dst);

        pair_array.emplace_back(std::make_pair(src, dst));

        if ((i + 1) % step_length == 0) {
            log_info("%u pairs have been generated", i + 1);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("generate pairs time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void generate_pairs_within_constrain(DirectedGraph *digraph, uint32_t num_pairs, std::vector<uint32_t>& src_array, std::vector<uint32_t>& dst_array,
        std::vector<pi>& pair_array, uint32_t length_constrain) {
    auto start = std::chrono::high_resolution_clock::now();
    uint32_t step_length = num_pairs * 0.2;
    uint32_t *distance;
    distance = (uint32_t*)malloc(sizeof(uint32_t)*(digraph->num_vertices()));

    log_info("Generator begin\n");

    for (uint32_t i = 0; i < num_pairs; ) {
        uint32_t src, dst;
        src = *select_randomly(src_array.begin(), src_array.end());
        GraphOperation::bfs_distance_from_target(digraph, src, distance, length_constrain, length_constrain+1);
        std::vector<uint32_t> dst2;
        for(auto iter = dst_array.begin(); iter!= dst_array.end();iter++){
            if((distance[*iter] > length_constrain) || ((*iter) == src)){
		continue;
            }
            else
            {
		dst2.push_back(*iter);
            }
        }
        if(dst2.empty())
            continue;
        i++;
        dst = *select_randomly(dst2.begin(), dst2.end());

        pair_array.emplace_back(std::make_pair(src, dst));

        if ((i + 1) % step_length == 0) {
            log_info("%u pairs have been generated", i + 1);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    log_info("generate pairs time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

void output_to_file(const std::string& file_path, std::vector<pi>& pair_array) {
    auto start = std::chrono::high_resolution_clock::now();
    std::ofstream ofs(file_path, std::ios::binary);
    uint32_t pair_size = sizeof(pi);
    uint32_t length = pair_array.size();
    uint64_t size = (uint64_t)length * pair_size;

    ofs.write(reinterpret_cast<const char *>(&pair_size), 4);
    ofs.write(reinterpret_cast<const char *>(&length), 4);
    ofs.write(reinterpret_cast<const char *>(&pair_array.front()), size);

    auto end = std::chrono::high_resolution_clock::now();
    log_info("store pairs file time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0);
}

int main(int argc, char *argv[]) {
    std::string input_graph_folder(argv[1]);
    std::string output_folder(argv[2]);
    uint32_t num_pairs = 1000;
    double hot_vertex_percentage = 0.1;
    uint32_t length_constrain = 5;

    if (argc > 3) {
        num_pairs = atoi(argv[3]);
    }

    if (argc > 4) {
        hot_vertex_percentage = atof(argv[4]);
    }

    DirectedGraph digraph;
    digraph.load_csr(input_graph_folder);
    digraph.print_metadata();

    // Separate hot vertices and unhot vertices.
    log_info("num of pairs = %u, percentage of hot vertices = %.2lf", num_pairs, hot_vertex_percentage);
    log_info("separate hot vertices with unhot vertices...");
    auto start = std::chrono::high_resolution_clock::now();

    uint32_t num_vertices = digraph.num_vertices();
    uint32_t num_hot_vertices = (uint32_t)(num_vertices * hot_vertex_percentage);

    std::priority_queue<pi, std::vector<pi>, std::greater<pi>> pq;
    for (uint32_t u = 0; u < num_hot_vertices; ++u) {
        pq.emplace(std::make_pair(digraph.num_neighbors(u), u));
    }

    for (uint32_t u = num_hot_vertices; u < num_vertices; ++u) {
        auto top = pq.top();
        uint32_t cur_degree = digraph.num_neighbors(u);

        if (cur_degree > top.first) {
            pq.pop();
            pq.emplace(std::make_pair(cur_degree, u));
        }
    }

    std::vector<uint32_t> hot_vertices;
    std::unordered_set<uint32_t> hot_vertices_flag;
    while (!pq.empty()) {
        auto top = pq.top();
        hot_vertices.emplace_back(top.second);
        hot_vertices_flag.emplace(top.second);
        pq.pop();
    }

    auto end1 = std::chrono::high_resolution_clock::now();
    log_info("obtain hot vertices time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start).count() / 1000.0);

    std::vector<uint32_t> unhot_vertices;
    for (uint32_t u = 0; u < num_vertices; ++u) {
        if (hot_vertices_flag.find(u) == hot_vertices_flag.end()) {
            unhot_vertices.emplace_back(u);
        }
    }

    auto end2 = std::chrono::high_resolution_clock::now();
    log_info("obtain unhot vertices time: %.3lf seconds",
             std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start).count() / 1000.0);

    // Generate pairs.
    std::vector<pi> pairs;

    // 1. hot to hot
    log_info("generate hot to hot pairs...");
    std::string output_file_path = output_folder + std::string("/hot2hot_pairs.bin");
    //generate_pairs(num_pairs, hot_vertices, hot_vertices, pairs);

    //use bfs to generate vertex pairs whose distance is within the length constrain
    generate_pairs_within_constrain(&digraph, num_pairs, hot_vertices, hot_vertices, pairs, length_constrain);
    output_to_file(output_file_path, pairs);
    pairs.clear();

    // 2. hot to unhot
    log_info("generate hot to unhot pairs...");
    output_file_path = output_folder + std::string("/hot2unhot_pairs.bin");
    //generate_pairs(num_pairs, hot_vertices, unhot_vertices, pairs);
    generate_pairs_within_constrain(&digraph, num_pairs, hot_vertices, unhot_vertices, pairs, length_constrain);
    output_to_file(output_file_path, pairs);
    pairs.clear();

    // 3. unhot to hot
    log_info("generate unhot to hot pairs...");
    output_file_path = output_folder + std::string("/unhot2hot_pairs.bin");
    //generate_pairs(num_pairs, unhot_vertices, hot_vertices, pairs);
    generate_pairs_within_constrain(&digraph, num_pairs, unhot_vertices, hot_vertices, pairs, length_constrain);
    output_to_file(output_file_path, pairs);
    pairs.clear();

    // 4. unhot to unhot
    log_info("generate unhot to unhot pairs...");
    output_file_path = output_folder + std::string("/unhot2unhot_pairs.bin");
    //generate_pairs(num_pairs, unhot_vertices, unhot_vertices, pairs);
    generate_pairs_within_constrain(&digraph, num_pairs, unhot_vertices, unhot_vertices, pairs, length_constrain);
    output_to_file(output_file_path, pairs);
    pairs.clear();

    log_info("done.");
    return 0;
}
