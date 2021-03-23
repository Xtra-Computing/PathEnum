//
// Created by Shixuan Sun on 2020/7/4.
//

#include <cstring>
#include <queue>
#include <algorithm>
#include "graph_operation.h"

void GraphOperation::bfs_distance_to_target(DirectedGraph *g, uint32_t u, uint32_t *distance, uint32_t invalid_value) {
    bool* visited = (bool*)malloc(sizeof(bool) * g->num_vertices());
    memset(visited, 0, sizeof(bool) * g->num_vertices());
    std::fill(distance, distance + g->num_vertices(), invalid_value);
    std::queue<uint32_t> q;

    q.push(u);
    visited[u] = true;
    distance[u] = 0;

    while (!q.empty()) {
        uint32_t v = q.front();
        q.pop();

        auto in_neighbors = g->in_neighbors(v);
        for (uint32_t i = 0; i < in_neighbors.second; ++i) {
            uint32_t vv = in_neighbors.first[i];

            if (!visited[vv]) {
                visited[vv] = true;
                distance[vv] = distance[v] + 1;
                q.push(vv);
            }
        }
    }

    free(visited);
}

void
GraphOperation::bfs_distance_to_target(DirectedGraph *g, uint32_t u, uint32_t *distance, uint32_t distance_constraint,
                                       uint32_t invalid_value) {
    bool* visited = (bool*)malloc(sizeof(bool) * g->num_vertices());
    memset(visited, 0, sizeof(bool) * g->num_vertices());
    std::fill(distance, distance + g->num_vertices(), invalid_value);
    std::queue<uint32_t> q;

    q.push(u);
    visited[u] = true;
    distance[u] = 0;

    while (!q.empty()) {
        uint32_t v = q.front();
        q.pop();

        if (distance[v] < distance_constraint) {
            auto in_neighbors = g->in_neighbors(v);
            for (uint32_t i = 0; i < in_neighbors.second; ++i) {
                uint32_t vv = in_neighbors.first[i];

                if (!visited[vv]) {
                    visited[vv] = true;
                    distance[vv] = distance[v] + 1;
                    q.push(vv);
                }
            }
        }
    }

    free(visited);
}

void
GraphOperation::bfs_distance_from_target(DirectedGraph *g, uint32_t u, uint32_t *distance, uint32_t distance_constraint,
                                       uint32_t invalid_value) {
    bool* visited = (bool*)malloc(sizeof(bool) * g->num_vertices());
    memset(visited, 0, sizeof(bool) * g->num_vertices());
    std::fill(distance, distance + g->num_vertices(), invalid_value);
    std::queue<uint32_t> q;

    q.push(u);
    visited[u] = true;
    distance[u] = 0;

    while (!q.empty()) {
        uint32_t v = q.front();
        q.pop();

        if (distance[v] < distance_constraint) {
            auto out_neighbors = g->out_neighbors(v);
            for (uint32_t i = 0; i < out_neighbors.second; ++i) {
                uint32_t vv = out_neighbors.first[i];

                if (!visited[vv]) {
                    visited[vv] = true;
                    distance[vv] = distance[v] + 1;
                    q.push(vv);
                }
            }
        }
    }

    free(visited);
}
