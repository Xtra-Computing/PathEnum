//
// Created by Shixuan Sun on 2020/6/30.
//
#include <chrono>
#include <future>
#include <numeric>
#include <iostream>

#include "util/log/log.h"
#include "util/graph/directed_graph.h"
#include "util/io/io.h"

#include "cycle_enumerator.h"

bool execute_within_time_limit(CycleEnumerator* enumerator, uint32_t src, uint32_t dst,
        CycleEnumerator::query_method method_type, uint64_t time_limit) {

    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [enumerator, src, dst, method_type](){
        return enumerator->execute(src, dst, method_type);
    });

    std::future_status status;
    do {
        status = future.wait_for(std::chrono::seconds(time_limit));
        if (status == std::future_status::deferred) {
            log_error("Deferred.");
            exit(-1);
        } else if (status == std::future_status::timeout) {
            g_exit = true;
        }
    } while (status != std::future_status::ready);

    return !g_exit;
}

int main(int argc, char *argv[]) {
    /**
     * Parse command.
     */
    std::string input_graph_folder(argv[1]);
    std::string input_query_file(argv[2]);
    std::string input_result_folder(argv[3]);
    std::string input_method_type(argv[4]);
    std::string input_length_constraint(argv[5]);
    std::string input_per_query_time_limit(argv[6]);

    CycleEnumerator::query_method method_type = CycleEnumerator::query_method::IDX_DFS;

    if (input_method_type == "IDX_DFS") {
        method_type = CycleEnumerator::query_method::IDX_DFS;
    }
    else if (input_method_type == "IDX_JOIN") {
        method_type = CycleEnumerator::query_method::IDX_JOIN;
    }
    else if (input_method_type == "PATH_ENUM") {
        method_type = CycleEnumerator::query_method::PATH_ENUM;
    }
    else if (input_method_type == "SPECTRUM") {
        method_type = CycleEnumerator::query_method::SPECTRUM;
    }
    else {
        std::cerr << "The input method type does not exist." << std::endl;
        exit(-1);
    }

    uint32_t length_constraint = std::stoul(input_length_constraint);
    uint32_t per_query_time_limit = std::stoul(input_per_query_time_limit);
    uint64_t target_num_results = std::numeric_limits<uint64_t>::max();

    if (argc > 7) {
        std::string input_target_number_results(argv[7]);
        target_num_results = std::stoull(input_target_number_results);
    }

    /**
     * Print command.
     */
    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << "Command:\n";
    std::cout << "Input graph folder: " << input_graph_folder << '\n';
    std::cout << "Input query file: " << input_query_file << '\n';
    std::cout << "Input result folder: " << input_result_folder << '\n';
    std::cout << "Input method: " << input_method_type << '\n';
    std::cout << "Input length constraint: " << length_constraint << '\n';
    std::cout << "Input per query time limit: " << per_query_time_limit << " seconds\n";
    std::cout << "Input target number of results: " << target_num_results << '\n';
    std::cout << "---------------------------------------------------------------------------" << std::endl;

    /**
     * Initialize graph, queries and enumerator.
     */
    DirectedGraph digraph;
    digraph.load_csr(input_graph_folder);
    digraph.print_metadata();

    std::vector<std::pair<uint32_t, uint32_t>> queries;
    IO::read(input_query_file, queries);

    CycleEnumerator enumerator;
    enumerator.initialize(&digraph, length_constraint);
    enumerator.target_number_results_ = target_num_results;
    /**
     * Execute queries.
     */
    uint32_t terminated_query_count = 0;
    uint32_t num_queries = queries.size();
    printf("num of queries: %u", num_queries);

    for (uint32_t i = 0; i < num_queries; ++i) {
        auto query = queries[i];

        std::string query_status = "Complete";
        if (!execute_within_time_limit(&enumerator, query.first, query.second, method_type, per_query_time_limit)) {
            terminated_query_count += 1;
            query_status = "Time Out";
        }

        printf("ID %u, src %u, dst %u, status %s, result count %lu, preprocessing time %.6lf seconds, query time %.6lf seconds,\n",
               i, query.first, query.second, query_status.c_str(), enumerator.result_count_arr_.back(),
               enumerator.preprocess_time_arr_.back() / (double)1000000000, enumerator.query_time_arr_.back() / (double)1000000000);
    }

    /**
     * Dump brief results to console.
     */
    uint64_t total_query_time = std::accumulate(enumerator.query_time_arr_.begin(), enumerator.query_time_arr_.end(), 0ull);
    uint64_t total_preprocess_time = std::accumulate(enumerator.preprocess_time_arr_.begin(), enumerator.preprocess_time_arr_.end(), 0ull);
    uint64_t total_result_count = std::accumulate(enumerator.result_count_arr_.begin(), enumerator.result_count_arr_.end(), 0ull);
    uint64_t total_partial_result_count = std::accumulate(enumerator.partial_result_count_arr_.begin(), enumerator.partial_result_count_arr_.end(), 0ull);
    uint64_t total_invalid_partial_result_count = std::accumulate(enumerator.invalid_partial_result_count_arr_.begin(), enumerator.invalid_partial_result_count_arr_.end(), 0ull);
    uint64_t total_neighbors_access_count = std::accumulate(enumerator.neighbors_access_count_arr_.begin(), enumerator.neighbors_access_count_arr_.end(), 0ull);
    uint64_t total_conflict_count = std::accumulate(enumerator.conflict_count_arr_.begin(), enumerator.conflict_count_arr_.end(), 0ull);
    uint64_t total_forward_bfs_time = std::accumulate(enumerator.forward_bfs_time_arr_.begin(), enumerator.forward_bfs_time_arr_.end(), 0ull);
    uint64_t total_backward_bfs_time = std::accumulate(enumerator.backward_bfs_time_arr_.begin(), enumerator.backward_bfs_time_arr_.end(), 0ull);
    uint64_t total_construct_bigraph_time = std::accumulate(enumerator.construct_bigraph_time_arr_.begin(), enumerator.construct_bigraph_time_arr_.end(), 0ull);
    uint64_t total_full_fledged_estimation_time = std::accumulate(enumerator.full_fledged_estimation_time_arr_.begin(), enumerator.full_fledged_estimation_time_arr_.end(), 0ull);
    uint64_t total_left_dfs_time = std::accumulate(enumerator.left_dfs_time_arr_.begin(), enumerator.left_dfs_time_arr_.end(), 0ull);
    uint64_t total_right_dfs_time = std::accumulate(enumerator.right_dfs_time_arr_.begin(), enumerator.right_dfs_time_arr_.end(), 0ull);
    uint64_t total_join_time = std::accumulate(enumerator.join_time_arr_.begin(), enumerator.join_time_arr_.end(), 0ull);


    double average_estimate_accuracy = std::accumulate(enumerator.estimate_accuracy_arr_.begin(), enumerator.estimate_accuracy_arr_.end(), 0.0);

    log_info("num of out of time queries: %u", terminated_query_count);
    log_info("average query time: %.6lf seconds", (double)total_query_time / num_queries / 1000000000.0);
    log_info("average preprocess time: %.6lf seconds", (double)total_preprocess_time / num_queries / 1000000000.0);
    log_info("average forward bfs time: %.6lf seconds", (double)total_forward_bfs_time / num_queries / 1000000000.0);
    log_info("average backward bfs time: %.6lf seconds", (double)total_backward_bfs_time / num_queries / 1000000000.0);
    log_info("average construct bigraph time: %.6lf seconds", (double)total_construct_bigraph_time / num_queries / 1000000000.0);
    log_info("average full fledged estimation time: %.6lf seconds", (double)total_full_fledged_estimation_time / num_queries / 1000000000.0);
    log_info("average left dfs time: %.6lf seconds", (double)total_left_dfs_time / num_queries / 1000000000.0);
    log_info("average right dfs time time: %.6lf seconds", (double)total_right_dfs_time / num_queries / 1000000000.0);
    log_info("average join time: %.6lf seconds", (double)total_join_time / num_queries / 1000000000.0);
    log_info("average result count: %.2lf", (double)total_result_count / num_queries);
    log_info("average partial result count: %.2lf", (double)total_partial_result_count / num_queries);
    log_info("average invalid partial result count: %.2lf", (double)total_invalid_partial_result_count / num_queries);
    log_info("average neighbors access count: %.2lf", (double)total_neighbors_access_count / num_queries);
    log_info("average conflict count: %.2lf", (double)total_conflict_count / num_queries);
    log_info("average estimate accuracy: %.4lf", average_estimate_accuracy);

    /**
     * Dump detailed results to file.
     */
    log_info("dump results into file...");
    std::string file_path_prefix = input_result_folder + "/" + input_method_type + "-" + input_length_constraint
                                   + "-" + std::to_string(num_queries);

    std::string file_path = file_path_prefix + std::string("-query_time.bin");
    IO::write(file_path, enumerator.query_time_arr_);

    file_path = file_path_prefix + std::string("-preprocess_time.bin");
    IO::write(file_path, enumerator.preprocess_time_arr_);

    file_path = file_path_prefix + std::string("-forward_bfs_time.bin");
    IO::write(file_path, enumerator.forward_bfs_time_arr_);

    file_path = file_path_prefix + std::string("-backward_bfs_time.bin");
    IO::write(file_path, enumerator.backward_bfs_time_arr_);

    file_path = file_path_prefix + std::string("-construct_bigraph_time.bin");
    IO::write(file_path, enumerator.construct_bigraph_time_arr_);

    file_path = file_path_prefix + std::string("-result_count.bin");
    IO::write(file_path, enumerator.result_count_arr_);

    file_path = file_path_prefix + std::string("-partial_result_count.bin");
    IO::write(file_path, enumerator.partial_result_count_arr_);

    file_path = file_path_prefix + std::string("-invalid_partial_result_count.bin");
    IO::write(file_path, enumerator.invalid_partial_result_count_arr_);

    file_path = file_path_prefix + std::string("-neighbors_access_count.bin");
    IO::write(file_path, enumerator.neighbors_access_count_arr_);

    file_path = file_path_prefix + std::string("-conflict_count.bin");
    IO::write(file_path, enumerator.conflict_count_arr_);

    file_path = file_path_prefix + std::string("-preliminary_estimated_result_count.bin");
    IO::write(file_path, enumerator.preliminary_estimated_result_count_arr);

    file_path = file_path_prefix + std::string("-full_fledged_estimated_result_count.bin");
    IO::write(file_path, enumerator.full_fledged_estimated_result_count_arr);

    file_path = file_path_prefix + std::string("-full_fledged_estimation_time.bin");
    IO::write(file_path, enumerator.full_fledged_estimation_time_arr_);

    file_path = file_path_prefix + std::string("-left_dfs_time.bin");
    IO::write(file_path, enumerator.left_dfs_time_arr_);

    file_path = file_path_prefix + std::string("-right_dfs_time.bin");
    IO::write(file_path, enumerator.right_dfs_time_arr_);

    file_path = file_path_prefix + std::string("-join_time.bin");
    IO::write(file_path, enumerator.join_time_arr_);

    file_path = file_path_prefix + std::string("-join_cost.bin");
    IO::write(file_path, enumerator.estimated_join_cost_arr_);

    file_path = file_path_prefix + std::string("-dfs_cost.bin");
    IO::write(file_path, enumerator.estimated_dfs_cost_arr_);

    file_path = file_path_prefix + std::string("-left_relation_size.bin");
    IO::write(file_path, enumerator.estimated_left_relation_size_arr_);

    file_path = file_path_prefix + std::string("-right_relation_size.bin");
    IO::write(file_path, enumerator.estimated_right_relation_size_arr_);

    file_path = file_path_prefix + std::string("-cut_position.bin");
    IO::write(file_path, enumerator.cut_position_arr_);

    file_path = file_path_prefix + std::string("-preliminary_selection.bin");
    IO::write(file_path, enumerator.preliminary_selection_arr_);

    file_path = file_path_prefix + std::string("-full_fledged_selection.bin");
    IO::write(file_path, enumerator.full_fledged_selection_arr_);

    file_path = file_path_prefix + std::string("-memory_cost.bin");
    IO::write(file_path, enumerator.calculated_memory_cost_arr_);

    file_path = file_path_prefix + std::string("-index_edge_count.bin");
    IO::write(file_path, enumerator.index_edge_count_arr_);

    file_path = file_path_prefix + std::string("-index_vertex_count.bin");
    IO::write(file_path, enumerator.index_vertex_count_arr_);

    file_path = file_path_prefix + std::string("-algorithm_selection.bin");
    IO::write(file_path, enumerator.algorithm_selected_arr_);

    file_path = file_path_prefix + std::string("-path_enum_exclude_index_construction.bin");
    IO::write(file_path, enumerator.path_enum_time_arr_);

    file_path = file_path_prefix + std::string("-dfs_spectrum.bin");
    IO::write(file_path, enumerator.dfs_spectrum_time_arr_);

    file_path = file_path_prefix + std::string("-join_spectrum.bin");
    IO::write(file_path, enumerator.join_spectrum_time_arr_);
    log_info("done.");
    return 0;
}
