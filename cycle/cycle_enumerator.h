//
// Created by Shixuan Sun on 2020/6/30.
//

#ifndef XTRAGRAPHCOMPUTING_CYCLE_ENUMERATOR_H
#define XTRAGRAPHCOMPUTING_CYCLE_ENUMERATOR_H

#include <vector>
#include "util/sparsepp/spp.h"
#include "util/graph/directed_graph.h"

extern bool g_exit;

#define LISTING_RESULTS
#define OUTPUT_BUFFER_RATIO (128 * 1024 * 1024)

class CycleEnumerator {
public:
    enum query_method {
        IDX_DFS,            // Hop constraint dfs search method on the bigraph.
        IDX_JOIN,           // Hop constraint join method on the bigraph.
        PATH_ENUM,          // Complete algorithm.
        SPECTRUM            // Execute spectrum analysis.
    };
public:
    /**
     * Performance counter.
     */
    std::vector<uint64_t> result_count_arr_;
    std::vector<uint64_t> query_time_arr_;  // In nanoseconds.
    std::vector<uint64_t> preprocess_time_arr_; // In nanoseconds.
    std::vector<uint64_t> partial_result_count_arr_; // Exclude result.
    std::vector<uint64_t> invalid_partial_result_count_arr_;
    std::vector<uint64_t> neighbors_access_count_arr_;
    std::vector<uint64_t> conflict_count_arr_;
    std::vector<uint64_t> preliminary_estimated_result_count_arr;
    std::vector<uint64_t> full_fledged_estimated_result_count_arr;
    std::vector<double> estimate_accuracy_arr_;

    uint64_t query_time_;
    uint64_t preprocess_time_;
    uint64_t partial_result_count_;
    uint64_t invalid_partial_result_count_;
    uint64_t neighbors_access_count_;
    uint64_t conflict_count_;
    uint64_t preliminary_estimated_result_count_;
    uint64_t full_fledged_estimated_result_count_;
    double estimate_accuracy_;

    /**
     * Detailed performance counter for bigraph index.
     */
    std::vector<uint64_t> forward_bfs_time_arr_;
    std::vector<uint64_t> backward_bfs_time_arr_;
    std::vector<uint64_t> construct_bigraph_time_arr_;
    std::vector<uint64_t> index_edge_count_arr_;
    std::vector<uint64_t> index_vertex_count_arr_;

    uint64_t forward_bfs_time_;
    uint64_t backward_bfs_time_;
    uint64_t construct_bigraph_time_;
    uint64_t index_edge_count_;
    uint64_t index_vertex_count_;

    /**
     * Detailed performance counter for join.
     */
    std::vector<uint64_t> full_fledged_estimation_time_arr_;
    std::vector<uint64_t> left_dfs_time_arr_;
    std::vector<uint64_t> right_dfs_time_arr_;
    std::vector<uint64_t> join_time_arr_;

    uint64_t full_fledged_estimation_time_;
    uint64_t left_dfs_time_;
    uint64_t right_dfs_time_;
    uint64_t join_time_;

    /**
     * Tune query optimizer.
     */
    std::vector<uint64_t> estimated_dfs_cost_arr_;
    std::vector<uint64_t> estimated_join_cost_arr_;
    std::vector<uint64_t> estimated_left_relation_size_arr_;
    std::vector<uint64_t> estimated_right_relation_size_arr_;
    std::vector<uint64_t> cut_position_arr_;
    std::vector<uint64_t> preliminary_selection_arr_;           // 0: DFS, 1: JOIN
    std::vector<uint64_t> full_fledged_selection_arr_;          // 0: DFS, 1: JOIN

    uint64_t estimated_dfs_cost_;
    uint64_t estimated_join_cost_;
    uint64_t preliminary_selection_;
    uint64_t full_fledged_selection_;

    /**
     * Estimate the memory consumption of the index or relation.
     */
    std::vector<uint64_t> calculated_memory_cost_arr_;
    uint64_t calculated_memory_cost_;

    /**
     * Spectrum analysis.
     */
    uint32_t algorithm_selected_; // 0: DFS, 1: JOIN
    uint64_t preliminary_estimation_time_;
    uint64_t path_enum_exclude_index_construction_time_;

    std::vector<uint64_t> dfs_spectrum_time_;
    std::vector<uint64_t> join_spectrum_time_;

    std::vector<uint32_t> algorithm_selected_arr_;
    std::vector<uint64_t> preliminary_estimation_time_arr_;
    std::vector<uint64_t> path_enum_time_arr_;
    std::vector<uint64_t> dfs_spectrum_time_arr_;
    std::vector<uint64_t> join_spectrum_time_arr_;
public:
    uint64_t target_number_results_;

private:
    /**
     * Input parameter.
     */
    DirectedGraph *digraph_;
    uint32_t length_constraint_;
    uint32_t src_;
    uint32_t dst_;
    query_method method_;

    /**
     * Helper data structures.
     */
    uint64_t count_;
    uint32_t *stack_;
    bool *visited_;

    /**
     * Bipartite graph index based method.
     */
    std::pair<uint8_t, uint8_t>* distance_;
    uint32_t* updated_values_;
    uint32_t* bucket_degree_sum_;

    uint32_t *buckets_offset_;
    uint32_t *buckets_;

    spp::sparse_hash_map<uint32_t, uint32_t> single_bigraph_;
    uint32_t *single_bigraph_offset_;
    uint32_t *single_bigraph_adj_;

    spp::sparse_hash_map<uint32_t, uint32_t> single_reverse_bigraph_;
    uint32_t *single_reverse_bigraph_offset_;
    uint32_t *single_reverse_bigraph_adj_;

    uint32_t min_cut_position_;
    uint64_t estimated_left_relation_size_;
    uint64_t estimated_right_relation_size_;
    uint64_t estimated_result_count_;

    uint32_t *left_relation_;
    uint64_t left_relation_size_;
    uint32_t *left_cursor_;
    uint32_t *left_partial_begin_;
    uint32_t *left_partial_end_;
    uint32_t left_part_length_;
    uint32_t *right_relation_;
    uint64_t right_relation_size_;
    uint32_t *right_cursor_;
    uint32_t right_part_length_;
    uint32_t *right_partial_begin_;
    uint32_t *right_partial_end_;
    spp::sparse_hash_map<uint32_t, std::pair<uint32_t*, uint64_t>> index_table_;

    /**
     * For listing test purpose.
     */
#ifdef LISTING_RESULTS
    uint32_t *output_buffer_;
#endif

public:
    CycleEnumerator();
    ~CycleEnumerator();

    void initialize(DirectedGraph* digraph, uint32_t length_constraint);
    uint64_t execute(uint32_t src, uint32_t dst, query_method method);
    void clear();

private:
    /**
     * DFS on the bipartite graph index.
     */
    void fast_build_bigraph();
    void dfs_on_bigraph(uint32_t u, uint32_t k);
    void clear_bigraph();

    /**
     * Spectrum analysis.
     */
    void dfs_spectrum();
    void dfs_on_bigraph(uint32_t depth, const uint32_t *order, uint32_t start_position);

    void join_spectrum();

    void generate_connected_order(std::vector<std::vector<uint32_t>>& orders);
    void generate_connected_order(std::vector<std::vector<uint32_t >> &orders, std::vector<uint32_t> &order,
                                  uint32_t cur_depth, int cur_lower, int cur_upper);
    void order_correctness(std::vector<uint32_t>& order);

    /**
     * Path count estimator and join operator.
     */
    bool preliminary_cardinality_estimator();
    void generate_single_join_plan();
    void single_join_on_bigraph();
    void left_dfs(uint32_t u, uint32_t k);
    void right_dfs(uint32_t u, uint32_t k);
    void single_join();

private:
    void initialize_performance_counter();
    void clear_performance_counter();
    void update_performance_counter();
    void reset_performance_counter();
};


#endif //XTRAGRAPHCOMPUTING_CYCLE_ENUMERATOR_H
