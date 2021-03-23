//
// Created by Shixuan Sun on 2020/6/30.
//

#include <cassert>
#include <chrono>
#include <queue>
#include <iostream>
#include "cycle_enumerator.h"
#include "util/graph/graph_operation.h"
#include "util/log/log.h"
#include "util/sparsepp/spp.h"

bool g_exit = false;

#define BUCKET_ID(i, j, l) ((i)*(l) + (j))

CycleEnumerator::CycleEnumerator() {
    stack_ = nullptr;
    visited_ = nullptr;
    length_constraint_ = 0;
    count_ = 0;

    buckets_offset_ = nullptr;
    buckets_ = nullptr;
    single_bigraph_offset_ = nullptr;
    single_bigraph_adj_ = nullptr;

    single_reverse_bigraph_adj_ = nullptr;
    single_reverse_bigraph_offset_ = nullptr;

    updated_values_ = nullptr;
    distance_ = nullptr;
}

CycleEnumerator::~CycleEnumerator() {
    clear();
}

void CycleEnumerator::initialize(DirectedGraph *digraph, uint32_t length_constraint) {
    assert(length_constraint > 0 && digraph != nullptr);

    digraph_ = digraph;
    length_constraint_ = length_constraint;

    uint64_t size = sizeof(uint32_t) * (length_constraint + 1);
    stack_ = (uint32_t*)malloc(size);
    memset(stack_, 0, size);

    size = sizeof(bool) * digraph->num_vertices();
    visited_ = (bool*)malloc(size);
    memset(visited_, 0, size);

    size = sizeof(uint32_t) * digraph->num_vertices();
    updated_values_ = (uint32_t*) malloc(size);
    memset(updated_values_, 0, size);

    size = sizeof(std::pair<uint8_t, uint8_t>) * digraph->num_vertices();
    distance_ = (std::pair<uint8_t, uint8_t>*) malloc(size);
    memset((uint8_t*)distance_, static_cast<uint8_t>(length_constraint) + 1, size);

    // Initialize the bipartite graph index.
    buckets_offset_ = (uint32_t*)malloc(sizeof(uint32_t) * ((length_constraint_ + 1) * (length_constraint_ + 1) + 1));
    memset(buckets_offset_, 0, sizeof(uint32_t) * ((length_constraint_ + 1) * (length_constraint_ + 1) + 1));
    bucket_degree_sum_ = (uint32_t*)malloc(sizeof(uint32_t) * (length_constraint_ + 1) * (length_constraint_ + 1) * (length_constraint_ + 1));
    memset(bucket_degree_sum_, 0, sizeof(uint32_t) * (length_constraint_ + 1) * (length_constraint_ + 1) * (length_constraint_ + 1));

    count_ = 0;
    initialize_performance_counter();

#ifdef LISTING_RESULTS
    output_buffer_ = (uint32_t*)malloc(sizeof(uint32_t) * length_constraint_ * (128 * 1024 * 1024));
#endif
}

void CycleEnumerator::clear() {
    free(stack_);
    stack_ = nullptr;
    free(visited_);
    visited_ = nullptr;

    free(bucket_degree_sum_);
    bucket_degree_sum_ = nullptr;
    free(buckets_offset_);
    buckets_offset_ = nullptr;
    free(buckets_);
    buckets_ = nullptr;

    length_constraint_ = 0;
    count_ = 0;

    free(updated_values_);
    updated_values_ = nullptr;
    free(distance_);
    distance_ = nullptr;

    clear_performance_counter();
}

uint64_t CycleEnumerator::execute(uint32_t src, uint32_t dst, query_method method) {
    assert(length_constraint_ > 0);
    auto start = std::chrono::high_resolution_clock::now();

    src_ = src;
    dst_ = dst;
    method_ = method;

    if (method_ == query_method::IDX_DFS) {
        fast_build_bigraph();
        preprocess_time_ = forward_bfs_time_ + backward_bfs_time_ + construct_bigraph_time_;
        dfs_on_bigraph(src_, 0);
        clear_bigraph();
    }
    else if (method_ == query_method::IDX_JOIN) {
        fast_build_bigraph();
        preprocess_time_ = forward_bfs_time_ + backward_bfs_time_ + construct_bigraph_time_;
        preliminary_cardinality_estimator();

        auto estimation_start = std::chrono::high_resolution_clock::now();

        generate_single_join_plan();

        auto estimation_end = std::chrono::high_resolution_clock::now();
        full_fledged_estimation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(estimation_end - estimation_start).count();

        single_join_on_bigraph();
        clear_bigraph();
    }
    else if (method_ == query_method::PATH_ENUM) {
        fast_build_bigraph();
        preprocess_time_ = forward_bfs_time_ + backward_bfs_time_ + construct_bigraph_time_;
        if (preliminary_cardinality_estimator()) {
            preliminary_selection_ = 1;
            auto estimation_start = std::chrono::high_resolution_clock::now();

            generate_single_join_plan();

            auto estimation_end = std::chrono::high_resolution_clock::now();
            full_fledged_estimation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(estimation_end - estimation_start).count();

            if (full_fledged_selection_ == 0) {
                dfs_on_bigraph(src_, 0);
            }
            else {
                single_join_on_bigraph();
            }
        }
        else {
            preliminary_selection_ = 0;
            dfs_on_bigraph(src_, 0);
        }
        clear_bigraph();
    }
    else if (method_ == query_method::SPECTRUM) {
        fast_build_bigraph();
        preprocess_time_ = forward_bfs_time_ + backward_bfs_time_ + construct_bigraph_time_;
        /***
         * Optimization time.
         */
        if (preliminary_cardinality_estimator()) {
            preliminary_selection_ = 1;

            generate_single_join_plan();

            if (full_fledged_selection_ == 0) {
                algorithm_selected_ = 0;
            }
            else {
                algorithm_selected_ = 1;
            }
        }
        else {
            preliminary_selection_ = 0;
            algorithm_selected_ = 0;
        }

        dfs_spectrum();

        join_spectrum();

        path_enum_exclude_index_construction_time_ = preliminary_estimation_time_;
        path_enum_exclude_index_construction_time_ += preliminary_selection_ == 0 ? 0 : full_fledged_estimation_time_;

        path_enum_exclude_index_construction_time_ += algorithm_selected_ == 0 ? dfs_spectrum_time_[0] : join_spectrum_time_[0];

        printf("Algorithm selected %d, Optimization %.6lf seconds, DFS %.6lf seconds, Join %.6lf seconds, PathEnum %.6lf seconds\n",
               algorithm_selected_, full_fledged_estimation_time_ / 1000000000.0,
               dfs_spectrum_time_[0] / 1000000000.0, join_spectrum_time_[0] / 1000000000.0,
               path_enum_exclude_index_construction_time_ / 1000000000.0);

        clear_bigraph();
    }
    auto end = std::chrono::high_resolution_clock::now();
    query_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    /**
     * The time spent on updating performance counter is excluded from the query time.
     */
    uint64_t temp_count = count_;

    update_performance_counter();
    reset_performance_counter();
    return temp_count;
}


void CycleEnumerator::initialize_performance_counter() {
    reset_performance_counter();
    preprocess_time_arr_.reserve(256);
    result_count_arr_.reserve(256);
    query_time_arr_.reserve(256);
    partial_result_count_arr_.reserve(256);
    invalid_partial_result_count_arr_.reserve(256);
    neighbors_access_count_arr_.reserve(256);
    conflict_count_arr_.reserve(256);
    preliminary_estimated_result_count_arr.reserve(256);
    full_fledged_estimated_result_count_arr.reserve(256);

    forward_bfs_time_arr_.reserve(256);
    backward_bfs_time_arr_.reserve(256);
    construct_bigraph_time_arr_.reserve(256);
    estimate_accuracy_arr_.reserve(256);

    full_fledged_estimation_time_arr_.reserve(256);
    left_dfs_time_arr_.reserve(256);
    right_dfs_time_arr_.reserve(256);
    join_time_arr_.reserve(256);

    estimated_join_cost_arr_.reserve(256);
    estimated_dfs_cost_arr_.reserve(256);
    estimated_left_relation_size_arr_.reserve(256);
    estimated_right_relation_size_arr_.reserve(256);
    cut_position_arr_.reserve(256);
    preliminary_selection_arr_.reserve(256);
    full_fledged_selection_arr_.reserve(256);

    calculated_memory_cost_arr_.reserve(256);
    index_vertex_count_arr_.reserve(256);
    index_edge_count_arr_.reserve(256);

    join_spectrum_time_arr_.reserve(256);
    dfs_spectrum_time_arr_.reserve(256);
    algorithm_selected_arr_.reserve(256);
    preliminary_estimation_time_arr_.reserve(256);
    path_enum_time_arr_.reserve(256);
}

void CycleEnumerator::clear_performance_counter() {
    reset_performance_counter();
    preprocess_time_arr_.clear();
    result_count_arr_.clear();
    query_time_arr_.clear();
    partial_result_count_arr_.clear();
    invalid_partial_result_count_arr_.clear();
    neighbors_access_count_arr_.clear();
    conflict_count_arr_.clear();

    preliminary_estimated_result_count_arr.clear();
    full_fledged_estimated_result_count_arr.clear();
    backward_bfs_time_arr_.clear();
    construct_bigraph_time_arr_.clear();
    estimate_accuracy_arr_.clear();

    full_fledged_estimation_time_arr_.clear();
    left_dfs_time_arr_.clear();
    right_dfs_time_arr_.clear();
    join_time_arr_.clear();

    estimated_join_cost_arr_.clear();
    estimated_dfs_cost_arr_.clear();
    estimated_left_relation_size_arr_.clear();
    estimated_right_relation_size_arr_.clear();
    cut_position_arr_.clear();
    preliminary_selection_arr_.clear();
    full_fledged_selection_arr_.clear();

    calculated_memory_cost_arr_.clear();
    index_edge_count_arr_.clear();
    index_vertex_count_arr_.clear();

    join_spectrum_time_arr_.clear();
    dfs_spectrum_time_arr_.clear();
    algorithm_selected_arr_.clear();
    preliminary_estimation_time_arr_.clear();
    path_enum_time_arr_.clear();
}

void CycleEnumerator::reset_performance_counter() {
    count_ = 0;
    query_time_ = 0;
    preprocess_time_ = 0;
    partial_result_count_ = 0;
    invalid_partial_result_count_ = 0;
    neighbors_access_count_ = 0;
    conflict_count_ = 0;
    preliminary_estimated_result_count_ = 0;
    full_fledged_estimated_result_count_ = 0;

    forward_bfs_time_ = 0;
    backward_bfs_time_ = 0;
    construct_bigraph_time_ = 0;
    estimate_accuracy_ = 0;

    full_fledged_estimation_time_ = 0;
    left_dfs_time_ = 0;
    right_dfs_time_ = 0;
    join_time_ = 0;

    estimated_join_cost_ = 0;
    estimated_dfs_cost_ = 0;
    preliminary_selection_ = 0;
    full_fledged_selection_ = 0;

    calculated_memory_cost_ = 0;
    index_edge_count_ = 0;
    index_vertex_count_ = 0;

    algorithm_selected_ = 0; // 0: DFS, 1: JOIN
    preliminary_estimation_time_ = 0;
    path_enum_exclude_index_construction_time_ = 0;

    dfs_spectrum_time_.clear();
    join_spectrum_time_.clear();
}

void CycleEnumerator::update_performance_counter() {
    result_count_arr_.emplace_back(count_);
    query_time_arr_.emplace_back(query_time_);
    preprocess_time_arr_.emplace_back(preprocess_time_);
    partial_result_count_arr_.emplace_back(partial_result_count_);
    invalid_partial_result_count_arr_.emplace_back(invalid_partial_result_count_);
    neighbors_access_count_arr_.emplace_back(neighbors_access_count_);
    conflict_count_arr_.emplace_back(conflict_count_);
    full_fledged_estimated_result_count_arr.emplace_back(full_fledged_estimated_result_count_);
    preliminary_estimated_result_count_arr.emplace_back(preliminary_estimated_result_count_);

    forward_bfs_time_arr_.emplace_back(forward_bfs_time_);
    backward_bfs_time_arr_.emplace_back(backward_bfs_time_);
    construct_bigraph_time_arr_.emplace_back(construct_bigraph_time_);
    estimate_accuracy_arr_.emplace_back(estimate_accuracy_);

    full_fledged_estimation_time_arr_.emplace_back(full_fledged_estimation_time_);
    right_dfs_time_arr_.emplace_back(right_dfs_time_);
    left_dfs_time_arr_.emplace_back(left_dfs_time_);
    join_time_arr_.emplace_back(join_time_);

    estimated_join_cost_arr_.emplace_back(estimated_join_cost_);
    estimated_dfs_cost_arr_.emplace_back(estimated_dfs_cost_);
    estimated_left_relation_size_arr_.emplace_back(estimated_left_relation_size_);
    estimated_right_relation_size_arr_.emplace_back(estimated_right_relation_size_);
    cut_position_arr_.emplace_back(min_cut_position_);
    preliminary_selection_arr_.emplace_back(preliminary_selection_);
    full_fledged_selection_arr_.emplace_back(full_fledged_selection_);

    calculated_memory_cost_arr_.emplace_back(calculated_memory_cost_);

    index_edge_count_arr_.emplace_back(index_edge_count_);
    index_vertex_count_arr_.emplace_back(index_vertex_count_);

    algorithm_selected_arr_.emplace_back(algorithm_selected_);
    path_enum_time_arr_.emplace_back(path_enum_exclude_index_construction_time_);
    preliminary_estimation_time_arr_.emplace_back(preliminary_estimation_time_);
    dfs_spectrum_time_arr_.insert(dfs_spectrum_time_arr_.end(), dfs_spectrum_time_.begin(), dfs_spectrum_time_.end());
    join_spectrum_time_arr_.insert(join_spectrum_time_arr_.end(), join_spectrum_time_.begin(), join_spectrum_time_.end());
}

void CycleEnumerator::fast_build_bigraph() {
    uint64_t num_edges = 0;
    auto forward_bfs_start = std::chrono::high_resolution_clock::now();

    uint32_t num_vertices = digraph_->num_vertices();
    uint32_t updated_values_count = 0;

    // Forward breadth-first search.
    if (length_constraint_ > std::numeric_limits<uint8_t>::max()) {
        std::cerr << "The length constraint is greater than uint8_t." << std::endl;
        exit(-1);
    }

    auto k = static_cast<uint8_t>(length_constraint_);
    std::queue<uint32_t> q;

    visited_[src_] = true;
    visited_[dst_] = true;

    updated_values_[updated_values_count++] = src_;
    updated_values_[updated_values_count++] = dst_;

    q.push(src_);
    distance_[src_].first = 0;

    while (!q.empty()) {
        uint32_t v = q.front();
        q.pop();

        if (distance_[v].first < k - 1) {
            uint8_t next_distance = distance_[v].first + 1;
            auto out_neighbors = digraph_->out_neighbors(v);
            for (uint32_t i = 0; i < out_neighbors.second; ++i) {
                uint32_t vv = out_neighbors.first[i];

                if (!visited_[vv]) {
                    visited_[vv] = true;
                    distance_[vv].first = next_distance;
                    updated_values_[updated_values_count++] = vv;
                    q.push(vv);
                }
            }
        }
    }

    auto forward_bfs_end = std::chrono::high_resolution_clock::now();
    forward_bfs_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(forward_bfs_end - forward_bfs_start).count();

    // Backward breadth-first search.
    std::vector<std::vector<uint32_t>> temp_buckets((length_constraint_ + 1) * (length_constraint_ + 1));

    // One cache line contains 64 boolean values and each AVX256 instruction can manipulate 32 boolean values.
    if (num_vertices / updated_values_count > 64 * 32) {
        for (uint32_t i = 0; i < updated_values_count; ++i) {
            visited_[updated_values_[i]] = false;
        }
    }
    else {
        memset(visited_, 0, sizeof(bool) * num_vertices);
    }

    q.push(dst_);
    visited_[dst_] = true;
    visited_[src_] = true;
    distance_[dst_].second = 0;

    uint32_t active_vertices_count = 0;

    while (!q.empty()) {
        uint32_t v = q.front();
        q.pop();

        if (distance_[v].second < k - 1) {
            uint32_t next_distance = distance_[v].second + 1;
            auto in_neighbors = digraph_->in_neighbors(v);
            for (uint32_t i = 0; i < in_neighbors.second; ++i) {
                uint32_t vv = in_neighbors.first[i];

                if (!visited_[vv] && distance_[vv].first + next_distance <= k) {
                    visited_[vv] = true;
                    distance_[vv].second = next_distance;
                    q.push(vv);

                    // Put vertices into temp buckets based on the distance_ to src and dst respectively.
                    active_vertices_count += 1;
                    uint32_t bucket_id = BUCKET_ID(distance_[vv].first, distance_[vv].second, k + 1);
                    temp_buckets[bucket_id].push_back(vv);
                }
            }
        }
    }

    auto backward_bfs_end = std::chrono::high_resolution_clock::now();
    backward_bfs_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(backward_bfs_end - forward_bfs_end).count();

    // Put vertices into different buckets. Include the src and dst vertex.
    active_vertices_count += 1;
    buckets_ = (uint32_t*)malloc(sizeof(uint32_t) * active_vertices_count);
    buckets_[0] = src_;

    uint32_t offset = 1;
    for (uint32_t i = 0; i < length_constraint_ + 1; ++i) {
        for (uint32_t j = 0; j < length_constraint_ + 1; ++j) {
            uint32_t bucket_id = BUCKET_ID(i, j, length_constraint_ + 1);
            buckets_offset_[bucket_id] = offset;
            memcpy(buckets_ + offset, temp_buckets[bucket_id].data(), sizeof(uint32_t) * temp_buckets[bucket_id].size());
            offset += temp_buckets[bucket_id].size();
            temp_buckets[bucket_id].clear();
        }
    }
    buckets_offset_[(length_constraint_ + 1) * (length_constraint_ + 1)] = offset;

    // Construct the forward bipartite graph.
    std::vector<uint32_t> temp_bigraph_adj;
    temp_bigraph_adj.reserve(1024);

    std::vector<std::vector<uint32_t>> temp_adj(length_constraint_);
    single_bigraph_offset_ = (uint32_t*)malloc(sizeof(uint32_t) * length_constraint_ * (active_vertices_count + 1));
    memset(single_bigraph_offset_, 0, sizeof(uint32_t) * length_constraint_ * (active_vertices_count + 1));

    uint32_t cur_bucket_id = 0;
    uint32_t cur_bucket_offset = buckets_offset_[cur_bucket_id + 1];
    uint32_t cur_bucket_max_degree_offset = cur_bucket_id * (length_constraint_ + 1);

    for (uint32_t i = 0; i < active_vertices_count; ++i) {
        uint32_t v = buckets_[i];

        // Find the bucket id.
        if (i != 0) {
            while (i >= cur_bucket_offset) {
                cur_bucket_id += 1;
                cur_bucket_max_degree_offset = cur_bucket_id * (length_constraint_ + 1);
                cur_bucket_offset = buckets_offset_[cur_bucket_id + 1];
            }
        }

        auto out_neighbors = digraph_->out_neighbors(v);

        for (uint32_t j = 0; j < out_neighbors.second; ++j) {
            uint32_t vv = out_neighbors.first[j];

            if (vv == dst_) {
                temp_adj[0].push_back(vv);
            }
            else if (visited_[vv] && distance_[vv].second < k) {
                temp_adj[distance_[vv].second].push_back(vv);
            }
        }

        uint32_t temp_offset = i * length_constraint_;

        uint32_t local_degree = 0;
        for (uint32_t j = 0; j < length_constraint_; ++j) {
            single_bigraph_offset_[temp_offset + j] = temp_bigraph_adj.size();
            temp_bigraph_adj.insert(temp_bigraph_adj.end(), temp_adj[j].begin(), temp_adj[j].end());
            local_degree += temp_adj[j].size();
            temp_adj[j].clear();

            // Collect statistics for the preliminary estimator.
            if (i != 0) {
                bucket_degree_sum_[cur_bucket_max_degree_offset + j] += local_degree;
            }
            else {
                bucket_degree_sum_[j] = local_degree;
            }
        }

        single_bigraph_offset_[temp_offset + length_constraint_] = temp_bigraph_adj.size();
        single_bigraph_[v] = temp_offset;
    }

    uint32_t temp_res_for_test1 = temp_bigraph_adj.size();
    num_edges += temp_res_for_test1;
    single_bigraph_adj_ = (uint32_t*)malloc(sizeof(uint32_t) * temp_bigraph_adj.size());
    memcpy(single_bigraph_adj_, temp_bigraph_adj.data(), sizeof(uint32_t) * temp_bigraph_adj.size());
    temp_bigraph_adj.clear();

    // Construct the backward bipartite graph.
    single_reverse_bigraph_offset_ = (uint32_t*)malloc(sizeof(uint32_t) * length_constraint_ * (active_vertices_count + 1));
    memset(single_reverse_bigraph_offset_, 0, sizeof(uint32_t) * length_constraint_ * (active_vertices_count + 1));

    buckets_[0] = dst_;
    for (uint32_t i = 0; i < active_vertices_count; ++i) {
        uint32_t v = buckets_[i];

        auto in_neighbors = digraph_->in_neighbors(v);

        for (uint32_t j = 0; j < in_neighbors.second; ++j) {
            uint32_t vv = in_neighbors.first[j];

            if (vv == src_) {
                temp_adj[0].push_back(vv);
            }
            else if (visited_[vv] && distance_[vv].first < length_constraint_) {
                temp_adj[distance_[vv].first].push_back(vv);
            }
        }

        uint32_t temp_offset = i * length_constraint_;
        for (uint32_t j = 0; j < length_constraint_; ++j) {
            single_reverse_bigraph_offset_[temp_offset + j] = temp_bigraph_adj.size();
            temp_bigraph_adj.insert(temp_bigraph_adj.end(), temp_adj[j].begin(), temp_adj[j].end());
            temp_adj[j].clear();
        }
        single_reverse_bigraph_offset_[temp_offset + length_constraint_] = temp_bigraph_adj.size();
        single_reverse_bigraph_[v] = temp_offset;
    }

    uint32_t temp_res_for_test2 = temp_bigraph_adj.size();
    //assert(temp_res_for_test1 == temp_res_for_test2);
    num_edges += temp_res_for_test2;
    single_reverse_bigraph_adj_ = (uint32_t*)malloc(sizeof(uint32_t) * temp_bigraph_adj.size());
    memcpy(single_reverse_bigraph_adj_, temp_bigraph_adj.data(), sizeof(uint32_t) * temp_bigraph_adj.size());

    if (num_vertices / updated_values_count > 16 * 8) {
        for (uint32_t i = 0; i < updated_values_count; ++i) {
            visited_[updated_values_[i]] = false;
            distance_[updated_values_[i]] = {k + 1, k + 1};
        }
    }
    else {
        memset(visited_, 0, sizeof(bool) * num_vertices);
        memset((uint8_t*)distance_, k + 1, sizeof(std::pair<uint8_t, uint8_t>) * num_vertices);
    }

    auto construct_bigraph_end = std::chrono::high_resolution_clock::now();
    construct_bigraph_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(construct_bigraph_end - backward_bfs_end).count();

    /**
     * Compute the memory cost of the index.
     */
    // Edges
    calculated_memory_cost_ = num_edges;
    // Offset
    calculated_memory_cost_ += (length_constraint_ + 1) * active_vertices_count * 2;
    // Hash table
    calculated_memory_cost_ += active_vertices_count * 2 * 2;
    calculated_memory_cost_ *= sizeof(uint32_t);

    index_edge_count_ = num_edges;
    index_vertex_count_ = num_vertices;
}

void
CycleEnumerator::dfs_on_bigraph(uint32_t depth, const uint32_t *order, uint32_t start_position) {
    uint32_t position = order[depth];
    uint32_t start, end;
    uint32_t *neighbors = nullptr;

    if (position > start_position) {
        // Go forward.
        uint32_t u = stack_[position - 1];
        if (u == dst_) {
            stack_[position] = dst_;
            if (depth == length_constraint_) {
                count_ += 1;
            }
            else {
                dfs_on_bigraph(depth + 1, order, start_position);
            }
            return;
        }

        uint32_t budget = length_constraint_ - (position - 1) - 1;

        uint32_t neighbor_offset = single_bigraph_[u];
        start = single_bigraph_offset_[neighbor_offset];
        end = single_bigraph_offset_[neighbor_offset + budget + 1];
        neighbors = single_bigraph_adj_;
    }
    else {
        // Go backward.
        uint32_t u = stack_[position + 1];

        if (u == dst_) {
            if (position > 0) {
               stack_[position] = dst_;
               dfs_on_bigraph(depth + 1, order, start_position);
            }
        }

        uint32_t budget = (position + 1) - 1;
        uint32_t neighbor_offset = single_reverse_bigraph_[u];
        start = single_reverse_bigraph_offset_[neighbor_offset];
        end = single_reverse_bigraph_offset_[neighbor_offset + budget + 1];
        neighbors = single_reverse_bigraph_adj_;
    }

    for (uint32_t i = start; i < end; ++i) {
        uint32_t v = neighbors[i];

        if (!visited_[v]) {
            if ((v == src_ && position != 0) || (position == 0 && v != src_))
                continue;
            stack_[position] = v;
            if (depth >= length_constraint_ - 1) {
                stack_[0] = src_;
                stack_[length_constraint_]= dst_;
                count_ += 1;
            }
            else {
                visited_[v] = true;
                dfs_on_bigraph(depth + 1, order, start_position);
                visited_[v] = false;
            }
        }
    }
}

void CycleEnumerator::dfs_on_bigraph(uint32_t u, uint32_t k) {
    stack_[k] = u;
    visited_[u] = true;

    uint64_t temp_count = count_;

    // k is cost; length_constraint_ - k is the remaining budget; minus 1 is the cost of moving to a out neighbor.
    uint32_t budget = length_constraint_ - k - 1;
    uint32_t neighbor_offset = single_bigraph_[u];
    uint32_t start = single_bigraph_offset_[neighbor_offset];
    uint32_t end = single_bigraph_offset_[neighbor_offset + budget + 1];

    neighbors_access_count_ += (end - start);

    for (uint32_t i = start; i < end; ++i) {
        if (g_exit || count_ >= target_number_results_) goto EXIT;

        uint32_t v = single_bigraph_adj_[i];
        if (v == dst_) {
            // Emit the result.
            stack_[k + 1] = dst_;
            count_ += 1;
        }
        else if (k == length_constraint_ - 2 && !visited_[v]) {
            // Emit the result.
            stack_[k + 1] = v;
            stack_[k + 2] = dst_;
            count_ += 1;
            partial_result_count_ += 1;
        }
        else if (!visited_[v]) {
            dfs_on_bigraph(v, k + 1);
        }
        else {
            conflict_count_ += 1;
        }
    }

    EXIT:
    visited_[u] = false;
    partial_result_count_ += 1;
    if (temp_count == count_) {
        invalid_partial_result_count_ += 1;
    }
}

void CycleEnumerator::clear_bigraph() {
    memset(buckets_offset_, 0, sizeof(uint32_t) * ((length_constraint_ + 1) * (length_constraint_ + 1) + 1));
    memset(bucket_degree_sum_, 0, sizeof(uint32_t) * (length_constraint_ + 1) * (length_constraint_ + 1) * (length_constraint_ + 1));
    free(single_bigraph_adj_);
    single_bigraph_adj_ = nullptr;
    free(single_bigraph_offset_);
    single_bigraph_offset_ = nullptr;
    single_bigraph_.clear();

    free(single_reverse_bigraph_adj_);
    single_reverse_bigraph_adj_ = nullptr;
    free(single_reverse_bigraph_offset_);
    single_reverse_bigraph_offset_ = nullptr;
    single_reverse_bigraph_.clear();
}

void CycleEnumerator::generate_single_join_plan() {
    uint32_t num_vertices = digraph_->num_vertices();
    std::vector<uint64_t> current_value(num_vertices, 0);
    std::vector<uint64_t> previous_value(num_vertices, 0);

    // Estimate the number of paths to dst.
    std::vector<uint64_t> num_path_to_dst(length_constraint_ + 1, 0);
    for (uint32_t budget = 1; budget < length_constraint_; ++budget) {
        // i is the distance to dst and j is the distance from src.
        uint64_t global_sum = 0;
        uint32_t remaining_budget = budget - 1;

        for (uint32_t i = 1; i <= budget; ++i) {
            for (uint32_t j = 1; j <= length_constraint_ - budget; ++j) {
                uint32_t bucket_id = BUCKET_ID(j, i, length_constraint_ + 1);
                for (uint32_t k = buckets_offset_[bucket_id]; k < buckets_offset_[bucket_id + 1]; ++k) {
                    uint32_t u = buckets_[k];
                    // budget is i - 1.
                    uint32_t neighbor_offset = single_bigraph_[u];
                    uint32_t start = single_bigraph_offset_[neighbor_offset];
                    uint32_t end = single_bigraph_offset_[neighbor_offset + remaining_budget + 1];

                    uint64_t local_sum = 0;
                    for (uint32_t l = start; l < end; ++l) {
                        uint32_t v = single_bigraph_adj_[l];
                        if (v == dst_) {
                            local_sum += 1;
                        } else {
                            local_sum += previous_value[v];
                        }
                    }
                    current_value[u] = local_sum;
                    global_sum += local_sum;
                }
            }
        }
        num_path_to_dst[budget] = global_sum;
        previous_value.swap(current_value);
    }
    {
        estimated_result_count_ = 0;
        uint32_t neighbor_offset = single_bigraph_[src_];
        uint32_t start = single_bigraph_offset_[neighbor_offset];
        uint32_t end = single_bigraph_offset_[neighbor_offset + length_constraint_];
        for (uint32_t i = start; i < end; ++i) {
            uint32_t u = single_bigraph_adj_[i];
            if (u == dst_) {
                estimated_result_count_ += 1;
            } else {
                estimated_result_count_ += previous_value[u];
            }
        }
    }
    // TODO: After testing, remove the reset.
    std::fill(previous_value.begin(), previous_value.end(), 0);
    std::fill(current_value.begin(), current_value.end(), 0);

    // Estimate the number of paths
    std::vector<uint64_t> num_path_from_src(length_constraint_ + 1, 0);
    for (uint32_t budget = 1; budget < length_constraint_; ++budget) {
        // i is the distance from src and j is the distance to dst.
        uint64_t global_sum = 0;
        uint32_t remaining_budget = budget - 1;

        for (uint32_t i = 1; i <= budget; ++i) {
            for (uint32_t j = 1; j <= length_constraint_ - budget; ++j) {
                uint32_t bucket_id = BUCKET_ID(i, j, length_constraint_ + 1);

                for (uint32_t k = buckets_offset_[bucket_id]; k < buckets_offset_[bucket_id + 1]; ++k) {
                    uint32_t u = buckets_[k];
                    // budget is i - 1.
                    uint32_t neighbor_offset = single_reverse_bigraph_[u];
                    uint32_t start = single_reverse_bigraph_offset_[neighbor_offset];
                    uint32_t end = single_reverse_bigraph_offset_[neighbor_offset + remaining_budget + 1];

                    uint64_t local_sum = 0;
                    for (uint32_t l = start; l < end; ++l) {
                        uint32_t v = single_reverse_bigraph_adj_[l];
                        if (v == src_) {
                            local_sum += 1;
                        } else {
                            local_sum += previous_value[v];
                        }
                    }

                    current_value[u] = local_sum;
                    global_sum += local_sum;
                }
            }
        }
        num_path_from_src[budget] = global_sum;
        previous_value.swap(current_value);
    }

    uint64_t min_sum = std::numeric_limits<uint64_t>::max();
    for (uint32_t i = 1; i < length_constraint_; ++i) {
        uint64_t cur_sum = num_path_to_dst[i] + num_path_from_src[length_constraint_ - i];
        if (cur_sum < min_sum) {
            min_sum = cur_sum;
            min_cut_position_ = length_constraint_ - i;
        }
    }

    estimated_left_relation_size_ = num_path_from_src[min_cut_position_];
    estimated_right_relation_size_ = num_path_to_dst[length_constraint_ - min_cut_position_];
    full_fledged_estimated_result_count_ = estimated_result_count_;

    // Estimated cost of DFS: the total number of intermediate results.
    estimated_dfs_cost_ = 0;
    for (uint32_t i = 1; i < length_constraint_; ++i) {
        estimated_dfs_cost_ += num_path_from_src[i];
    }

    // Materialization cost of the partial results + loop over the results. 1.05 is the penalty of checking the duplicate vertices.
    estimated_join_cost_ = estimated_left_relation_size_ + estimated_right_relation_size_ + (uint64_t)(estimated_result_count_ * 1.05);

    full_fledged_selection_ = estimated_join_cost_ < estimated_dfs_cost_ ? 1 : 0;
}

void CycleEnumerator::left_dfs(uint32_t u, uint32_t k) {
    stack_[k] = u;
    visited_[u] = true;

    // k is cost; length_constraint_ - k is the remaining budget; minus 1 is the cost of moving to a out neighbor.
    uint32_t budget = length_constraint_ - k - 1;
    uint32_t neighbor_offset = single_bigraph_[u];
    uint32_t start = single_bigraph_offset_[neighbor_offset];
    uint32_t end = single_bigraph_offset_[neighbor_offset + budget + 1];

    neighbors_access_count_ += (end - start);

    for (uint32_t i = start; i < end; ++i) {
        if (g_exit) goto EXIT;

        uint32_t v = single_bigraph_adj_[i];
        if (v == dst_) {
            // Emit the result.
            stack_[k + 1] = dst_;
            count_ += 1;
        }
        else if (k == min_cut_position_ - 1 && !visited_[v]) {
            stack_[k + 1] = v;
            // Copy the result to the buffer.
            std::copy(left_partial_begin_, left_partial_end_, left_cursor_);
            left_cursor_ += left_part_length_;
            left_relation_size_ += 1;
            partial_result_count_ += 1;
        }
        else if (!visited_[v]) {
            left_dfs(v, k + 1);
        }
        else {
            conflict_count_ += 1;
        }
    }

    EXIT:
    partial_result_count_ += 1;
    visited_[u] = false;
}

void CycleEnumerator::right_dfs(uint32_t u, uint32_t k) {
    stack_[k] = u;
    visited_[u] = true;

    // k is cost; length_constraint_ - k is the remaining budget; minus 1 is the cost of moving to a out neighbor.
    uint32_t budget = length_constraint_ - k - 1;
    uint32_t neighbor_offset = single_bigraph_[u];
    uint32_t start = single_bigraph_offset_[neighbor_offset];
    uint32_t end = single_bigraph_offset_[neighbor_offset + budget + 1];

    neighbors_access_count_ += (end - start);

    for (uint32_t i = start; i < end; ++i) {
        if (g_exit) goto EXIT;

        uint32_t v = single_bigraph_adj_[i];
        if (v == dst_) {
            // Emit the result.
            stack_[k + 1] = dst_;
            std::copy(right_partial_begin_, right_partial_end_, right_cursor_);
            right_cursor_ += right_part_length_;
            right_relation_size_ += 1;
        }
        else if (k == length_constraint_ - 2 && !visited_[v]) {
            stack_[k + 1] = v;
            stack_[k + 2] = dst_;
            // Copy the result to the buffer.
            std::copy(right_partial_begin_, right_partial_end_, right_cursor_);
            right_cursor_ += right_part_length_;
            right_relation_size_ += 1;
            partial_result_count_ += 1;
        }
        else if (!visited_[v]) {
            right_dfs(v, k + 1);
        }
        else {
            conflict_count_ += 1;
        }
    }

    EXIT:
    partial_result_count_ += 1;
    visited_[u] = false;
}

void CycleEnumerator::single_join() {
    left_cursor_ = left_relation_;
    uint32_t left_key_position = left_part_length_ - 1;
    for (uint64_t i = 0; i < left_relation_size_; ++i) {
        // Initialize visited table.
        for (uint32_t j = 0; j < left_part_length_; ++j) {
            uint32_t u = left_cursor_[j];
            visited_[u] = true;
        }

        // Join with the partitions.
        uint32_t key = left_cursor_[left_key_position];
        if (index_table_.contains(key)) {
            auto partitions = index_table_[key];
            right_cursor_ = partitions.first;
            for (uint64_t j = 0; j < partitions.second; ++j) {
                if(g_exit){
                    return;
                }
                for (uint32_t k = 1; k < right_part_length_; ++k) {
                    uint32_t u = right_cursor_[k];
                    if (u == dst_) {
                        count_ += 1;
                        break;
                    } else if (visited_[u]) {
                        break;
                    }
                }
                right_cursor_ += right_part_length_;
            }
        }

        // Clear visited table.
        for (uint32_t j = 0; j < left_part_length_; ++j) {
            uint32_t u = left_cursor_[j];
            visited_[u] = false;
        }

        left_cursor_ += left_part_length_;
    }
}

void CycleEnumerator::single_join_on_bigraph() {
    // Initialize.
    left_part_length_ = min_cut_position_ + 1;
    right_part_length_ = length_constraint_ - min_cut_position_ + 1;
    left_relation_size_ = 0;
    right_relation_size_ = 0;
    left_partial_begin_ = stack_;
    left_partial_end_ = left_partial_begin_ + left_part_length_;
    right_partial_begin_ = stack_ + min_cut_position_;
    right_partial_end_ = right_partial_begin_ + right_part_length_;

    left_relation_ = (uint32_t*)malloc(sizeof(uint32_t) * left_part_length_ * estimated_left_relation_size_);
    right_relation_ = (uint32_t*)malloc(sizeof(uint32_t) * right_part_length_ * estimated_right_relation_size_);

    auto left_dfs_start = std::chrono::high_resolution_clock::now();
    // Allocate the memory for the materialization.
    left_cursor_ = left_relation_;

    left_dfs(src_, 0);
    auto left_dfs_end = std::chrono::high_resolution_clock::now();
    left_dfs_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(left_dfs_end - left_dfs_start).count();

    // Allocate the memory for the materialization.
    right_cursor_ = right_relation_;

    for (uint32_t i = 1; i <= min_cut_position_; ++i) {
        for (uint32_t j = 1; j <= length_constraint_ - min_cut_position_; ++j) {
            uint32_t bucket_id = BUCKET_ID(i, j, length_constraint_ + 1);
            for (uint32_t k = buckets_offset_[bucket_id]; k < buckets_offset_[bucket_id + 1]; ++k) {
                uint32_t u = buckets_[k];
                uint32_t* cursor = right_cursor_;
                right_dfs(u, min_cut_position_);
                uint64_t temp_count = (right_cursor_ - cursor) / right_part_length_;
                index_table_[u] = std::make_pair(cursor, temp_count);
            }
        }
    }

    auto right_dfs_end = std::chrono::high_resolution_clock::now();
    right_dfs_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(right_dfs_end - left_dfs_end).count();

    single_join();

    estimate_accuracy_ = estimated_result_count_ == 0? 1 : (double)count_ / estimated_result_count_;

    // Release the memory.

    auto join_end = std::chrono::high_resolution_clock::now();

    join_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(join_end - right_dfs_end).count();

    free(left_relation_);
    free(right_relation_);
    index_table_.clear();
//    auto free_end = std::chrono::high_resolution_clock::now();
    
//    uint64_t free_time = std::chrono::duration_cast<std::chrono::nanoseconds>(free_end - join_end).count();
//    uint64_t malloc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(left_dfs_start - estimation_end).count();
//    printf("Left relation size: %u, Right relation size: %u\n", left_relation_size_, right_relation_size_);
//    printf("Free time: %.6lf, Malloc time: %.6lf\n", free_time/1000000.0, malloc_time/1000000.0);

}

bool CycleEnumerator::preliminary_cardinality_estimator() {
    auto estimation_start = std::chrono::high_resolution_clock::now();
    // the first bucket stores src.
    preliminary_estimated_result_count_ = bucket_degree_sum_[length_constraint_ - 1];
    for (uint32_t i = 1; i < length_constraint_; ++i) {
        uint32_t degree_sum = 1;
        uint32_t vertex_sum = 1;
        for (uint32_t j = 1; j <= i; ++j) {
            for (uint32_t k = 1; k <= length_constraint_ - i; ++k) {
                uint32_t budget = length_constraint_ - i - 1;
                uint32_t cur_degree = bucket_degree_sum_[BUCKET_ID(j, k, length_constraint_ + 1) * (length_constraint_ + 1) + budget];
                uint32_t cur_vertex = buckets_offset_[BUCKET_ID(j, k, length_constraint_ + 1) + 1] - buckets_offset_[BUCKET_ID(j, k, length_constraint_ + 1)];
                degree_sum += cur_degree;
                vertex_sum += cur_vertex;
            }
        }

        preliminary_estimated_result_count_ *= (degree_sum / vertex_sum) > 1 ? (degree_sum / vertex_sum) : 1;
    }

    auto estimation_end = std::chrono::high_resolution_clock::now();
    preliminary_estimation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(estimation_end - estimation_start).count();

    // If the number of estimated results > 100000, then invoke the fine grained optimizer. Our experiment results
    // show that this estimation generally underestimates the number of results.
    return preliminary_estimated_result_count_ > 100000;
}

void CycleEnumerator::generate_connected_order(std::vector<std::vector<uint32_t>> &orders) {
    std::vector<uint32_t> order(length_constraint_ + 1);
    for (uint32_t i = 0; i <= length_constraint_; ++i) {
        order[0] = i;
        generate_connected_order(orders, order, 1, i, i);
    }
}

void
CycleEnumerator::generate_connected_order(std::vector<std::vector<uint32_t >> &orders, std::vector<uint32_t> &order,
                                          uint32_t cur_depth, int cur_lower, int cur_upper) {
    if (cur_depth > length_constraint_) {
        order_correctness(order);
        orders.emplace_back(order);
        return;
    }

    if (cur_lower > 0) {
        order[cur_depth] = cur_lower - 1;
        generate_connected_order(orders, order, cur_depth + 1, cur_lower - 1, cur_upper);
    }

    if (cur_upper < length_constraint_) {
        order[cur_depth] = cur_upper + 1;
        generate_connected_order(orders, order, cur_depth + 1, cur_lower, cur_upper + 1);
    }
}

void CycleEnumerator::order_correctness(std::vector<uint32_t> &order) {
    std::vector<bool> visited(length_constraint_ + 1, false);
    assert(order.size() == visited.size());

    uint32_t start = order[0];
    assert(start >= 0 && start <= length_constraint_);
    visited[start] = true;

    for (uint32_t i = 1; i <= length_constraint_; ++i) {
        uint32_t cur = order[i];
        assert(cur >= 0 && cur <= length_constraint_);
        assert(!visited[cur]);
        visited[cur] = true;

        if (cur < start) {
            assert(visited[cur + 1]);
        }
        else {
            assert(visited[cur - 1]);
        }
    }
}

void CycleEnumerator::dfs_spectrum() {
    std::vector<std::vector<uint32_t>> orders;

    generate_connected_order(orders);
    printf("Source %d => Target %d, Left Deep Join Order Count %d\n", src_, dst_, orders.size());

    uint64_t correct_count = 0;
    int id = 0;
    for (auto& order : orders) {
        auto start_time = std::chrono::high_resolution_clock::now();
        count_ = 0;
        if (order[0] == 0) {
            dfs_on_bigraph(src_, 0);
            correct_count = count_;
        }
        else {
            uint32_t start = order[0];

            for (uint32_t i = 1; i <= start; ++i) {
                for (uint32_t j = 1; j <= length_constraint_ - start; ++j) {
                    uint32_t bucket_id = BUCKET_ID(i, j, length_constraint_ + 1);
                    for (uint32_t k = buckets_offset_[bucket_id]; k < buckets_offset_[bucket_id + 1]; ++k) {
                        uint32_t u = buckets_[k];
                        stack_[start] = u;
                        visited_[u] = true;
                        dfs_on_bigraph(1, order.data(), start);
                        visited_[u] = false;
                    }
                }
            }

            stack_[start] = dst_;
            visited_[dst_] = true;
            dfs_on_bigraph(1, order.data(), start);
            visited_[dst_] = false;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto enumerate_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
        dfs_spectrum_time_.emplace_back(enumerate_time);

        printf("%i: %zu, %.3lf ms\n", id++, count_, enumerate_time / 1000000.0);
        assert(correct_count == count_);
    }
}

void CycleEnumerator::join_spectrum() {
    auto estimation_start = std::chrono::high_resolution_clock::now();
    uint32_t num_vertices = digraph_->num_vertices();
    std::vector<uint64_t> current_value(num_vertices, 0);
    std::vector<uint64_t> previous_value(num_vertices, 0);

    // Estimate the number of paths to dst.
    std::vector<uint64_t> num_path_to_dst(length_constraint_ + 1, 0);
    for (uint32_t budget = 1; budget < length_constraint_; ++budget) {
        // i is the distance to dst and j is the distance from src.
        uint64_t global_sum = 0;
        uint32_t remaining_budget = budget - 1;

        for (uint32_t i = 1; i <= budget; ++i) {
            for (uint32_t j = 1; j <= length_constraint_ - budget; ++j) {
                uint32_t bucket_id = BUCKET_ID(j, i, length_constraint_ + 1);
                for (uint32_t k = buckets_offset_[bucket_id]; k < buckets_offset_[bucket_id + 1]; ++k) {
                    uint32_t u = buckets_[k];
                    // budget is i - 1.
                    uint32_t neighbor_offset = single_bigraph_[u];
                    uint32_t start = single_bigraph_offset_[neighbor_offset];
                    uint32_t end = single_bigraph_offset_[neighbor_offset + remaining_budget + 1];

                    uint64_t local_sum = 0;
                    for (uint32_t l = start; l < end; ++l) {
                        uint32_t v = single_bigraph_adj_[l];
                        if (v == dst_) {
                            local_sum += 1;
                        } else {
                            local_sum += previous_value[v];
                        }
                    }
                    current_value[u] = local_sum;
                    global_sum += local_sum;
                }
            }
        }
        num_path_to_dst[budget] = global_sum;
        previous_value.swap(current_value);
    }
    {
        estimated_result_count_ = 0;
        uint32_t neighbor_offset = single_bigraph_[src_];
        uint32_t start = single_bigraph_offset_[neighbor_offset];
        uint32_t end = single_bigraph_offset_[neighbor_offset + length_constraint_];
        for (uint32_t i = start; i < end; ++i) {
            uint32_t u = single_bigraph_adj_[i];
            if (u == dst_) {
                estimated_result_count_ += 1;
            } else {
                estimated_result_count_ += previous_value[u];
            }
        }
    }
    // TODO: After testing, remove the reset.
    std::fill(previous_value.begin(), previous_value.end(), 0);
    std::fill(current_value.begin(), current_value.end(), 0);

    // Estimate the number of paths
    std::vector<uint64_t> num_path_from_src(length_constraint_ + 1, 0);
    for (uint32_t budget = 1; budget < length_constraint_; ++budget) {
        // i is the distance from src and j is the distance to dst.
        uint64_t global_sum = 0;
        uint32_t remaining_budget = budget - 1;

        for (uint32_t i = 1; i <= budget; ++i) {
            for (uint32_t j = 1; j <= length_constraint_ - budget; ++j) {
                uint32_t bucket_id = BUCKET_ID(i, j, length_constraint_ + 1);

                for (uint32_t k = buckets_offset_[bucket_id]; k < buckets_offset_[bucket_id + 1]; ++k) {
                    uint32_t u = buckets_[k];
                    // budget is i - 1.
                    uint32_t neighbor_offset = single_reverse_bigraph_[u];
                    uint32_t start = single_reverse_bigraph_offset_[neighbor_offset];
                    uint32_t end = single_reverse_bigraph_offset_[neighbor_offset + remaining_budget + 1];

                    uint64_t local_sum = 0;
                    for (uint32_t l = start; l < end; ++l) {
                        uint32_t v = single_reverse_bigraph_adj_[l];
                        if (v == src_) {
                            local_sum += 1;
                        } else {
                            local_sum += previous_value[v];
                        }
                    }

                    current_value[u] = local_sum;
                    global_sum += local_sum;
                }
            }
        }
        num_path_from_src[budget] = global_sum;
        previous_value.swap(current_value);
    }

    uint64_t min_sum = std::numeric_limits<uint64_t>::max();
    for (uint32_t i = 1; i < length_constraint_; ++i) {
        uint64_t cur_sum = num_path_to_dst[i] + num_path_from_src[length_constraint_ - i];
        if (cur_sum < min_sum) {
            min_sum = cur_sum;
            min_cut_position_ = length_constraint_ - i;
        }
    }

    auto estimation_end = std::chrono::high_resolution_clock::now();
    full_fledged_estimation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(estimation_end - estimation_start).count();


    uint32_t id = 0;
    uint32_t selected_min_cut_position = min_cut_position_;

    estimated_left_relation_size_ = num_path_from_src[min_cut_position_];
    estimated_right_relation_size_ = num_path_to_dst[length_constraint_ - min_cut_position_];

    {
        count_ = 0;
        printf("Source %d => Target %d, Bushy Join Order Count %d\n", src_, dst_, length_constraint_ - 1);
        auto start_time = std::chrono::high_resolution_clock::now();
        single_join_on_bigraph();
        auto end_time = std::chrono::high_resolution_clock::now();
        auto enumerate_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
        join_spectrum_time_.emplace_back(enumerate_time);
        printf("%i: %zu, %.3lf ms, %d\n", id++, count_, enumerate_time / 1000000.0, min_cut_position_);
    }


    uint64_t correct_count = count_;

    for (uint32_t i = 1; i < length_constraint_; ++i) {
        if (i != selected_min_cut_position) {
            count_ = 0;
            min_cut_position_ = i;
            estimated_left_relation_size_ = num_path_from_src[min_cut_position_];
            estimated_right_relation_size_ = num_path_to_dst[length_constraint_ - min_cut_position_];

            uint64_t max_size = 42;
            max_size = max_size * 1024 * 1024 * 1024 / 4;
            assert(estimated_left_relation_size_ < max_size && estimated_right_relation_size_ < max_size);
            auto start_time = std::chrono::high_resolution_clock::now();

            single_join_on_bigraph();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto enumerate_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
            join_spectrum_time_.emplace_back(enumerate_time);
            printf("%i: %zu, %.3lf ms, %d\n", id++, count_, enumerate_time / 1000000.0, min_cut_position_);
            assert(correct_count == count_);
        }
    }

}
