//
// Created by Shixuan Sun on 2020/6/30.
//

#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

#include <string>
#include <vector>

class DirectedGraph {
public:
    /**
     * The meta information.
     */
    uint32_t num_vertices_;
    uint32_t num_edges_;
    uint32_t max_num_in_neighbors_;
    uint32_t max_num_out_neighbors_;

    std::vector<uint32_t> out_degree_distribution;
    std::vector<uint32_t> in_degree_distribution;

    /**
     * CSR representation.
     */
    uint32_t* out_offset_;
    uint32_t* out_adj_;
    uint32_t* out_degree_;

    uint32_t* in_offset_;
    uint32_t* in_adj_;
    uint32_t* in_degree_;
public:
    explicit DirectedGraph() : num_vertices_(0), num_edges_(0), max_num_in_neighbors_(0), max_num_out_neighbors_(0),
                               out_offset_(nullptr), out_adj_(nullptr), out_degree_(nullptr),
                               in_offset_(nullptr), in_adj_(nullptr), in_degree_(nullptr) {}

    ~DirectedGraph() {
        free(out_offset_);
        free(in_offset_);
        free(out_adj_);
        free(in_adj_);
        free(in_degree_);
        free(out_degree_);
    }

    inline uint32_t num_vertices() {
        return num_vertices_;
    }

    inline uint32_t num_edges() {
        return num_edges_;
    }

    inline uint32_t max_num_in_neighbors() {
        return max_num_in_neighbors_;
    }

    inline uint32_t max_num_out_neighbors() {
        return max_num_out_neighbors_;
    }

    inline std::pair<uint32_t*, uint32_t> out_neighbors(uint32_t u) {
        return std::make_pair(out_adj_ + out_offset_[u], out_degree_[u]);
    }

    inline std::pair<uint32_t*, uint32_t> in_neighbors(uint32_t u) {
        return std::make_pair(in_adj_ + in_offset_[u], in_degree_[u]);
    }

    inline uint32_t num_out_neighbors(uint32_t u) {
        return out_degree_[u];
    }

    inline uint32_t num_in_neighbors(uint32_t u) {
        return in_degree_[u];
    }

    inline uint32_t num_neighbors(uint32_t u) {
        return in_degree_[u] + out_degree_[u];
    }

public:
    /**
     * Graph I/O operations.
     */
    void load_csr(const std::string& graph_dir);

    void store_csr(const std::string& graph_dir);

    void load_edge_list(const std::string& graph_dir, char skip='#');

    void store_edge_list(const std::string& graph_dir);

    void print_metadata();

private:
    /**
     * Graph I/O helper functions.
     */
    void load_csr_degree_file(const std::string &degree_file_path, uint32_t *&degree);
    void load_csr_adj_file(const std::string &adj_file_path, const uint32_t *degree, uint32_t *&offset,
                           uint32_t *&adj);

    void load_edge_list(const std::string &edge_list_file_path, std::vector<std::pair<uint32_t, uint32_t>> &lines,
                        uint32_t &max_vertex_id, char skip);
    void load_edge_list(std::vector<std::pair<uint32_t, uint32_t> > &edge_list, uint32_t max_vertex_id,
                        uint32_t *&offset, uint32_t *&degree, uint32_t *&adj, uint32_t key_pos);

    /**
     * Graph manipulation functions.
     */
    void collect_metadata();
};


#endif // DIRECTED_GRAPH_H
