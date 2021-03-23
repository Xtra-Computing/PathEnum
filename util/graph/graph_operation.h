//
// Created by Shixuan Sun on 2020/7/4.
//

#ifndef XTRAGRAPHCOMPUTING_GRAPH_OPERATION_H
#define XTRAGRAPHCOMPUTING_GRAPH_OPERATION_H

#include "directed_graph.h"
class GraphOperation {
public:
    static void bfs_distance_to_target(DirectedGraph *g, uint32_t u, uint32_t *distance, uint32_t invalid_value);
    static void bfs_distance_to_target(DirectedGraph *g, uint32_t u, uint32_t *distance,
                                       uint32_t distance_constraint, uint32_t invalid_value);
    static void bfs_distance_from_target(DirectedGraph *g, uint32_t u, uint32_t *distance,
                                       uint32_t distance_constraint, uint32_t invalid_value);
};


#endif //XTRAGRAPHCOMPUTING_GRAPH_OPERATION_H
