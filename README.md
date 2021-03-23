# PathEnum
## Introduction

We study the hop-constrained s-t path enumeration (HcPE) problem, which takes a graph G, two distinct vertices s,t and a hop constraint k as input, and outputs all paths from s to t whose length is at most k. The state-of-the-art algorithms suffer from severe performance issues caused by the costly pruning operations during enumeration for the workloads with the large search space. Consequently, these algorithms hardly meet the real-time constraints of many online applications. In this paper, we propose PathEnum, an efficient index-based algorithm towards real-time HcPE. For an input query, PathEnum first builds a light-weight index aiming to reduce the number of edges involved in the enumeration, and develops efficient index-based approaches for enumeration, one based on depth-first search and the other based on joins. We further develop a query optimizer based on a join-based cost model to optimize the search order. We conduct experiments with 15 real-world graphs. Our experiment results show that PathEnum outperforms the state-of-the-art approaches by orders of magnitude in terms of the query time, throughput and response time.

For the details, please refer to our SIGMOD'2021 paper
"PathEnum: Towards Real-Time Hop-Constrained s-t Path Enumeration"
by [Dr. Shixuan Sun](https://shixuansun.github.io/), [Yuhang Chen](https://alexcyh7.github.io/),
[Dr. Bingsheng He](https://www.comp.nus.edu.sg/~hebs/), and [Dr. Bryan Hooi](https://bhooi.github.io/).
If you have any further questions, please feel free to contact us.


## Compile
Under the root directory of the project, execute the following commands to compile the source code 

```zsh
mkdir build
cd build
cmake ..
make
```

## Test
Under the root directory of the project, execute the following commands to test the correctness of the binary file.

```zsh
cd cycle/script
python test_cycle_enumerator.py ../../build/cycle/CycleEnumerator.out
```

## Preprocessing
After compiling the source code, you can find the binary file `GenerateVertexPairs.out` and `EdgeList2DirectedCSR.out` under the `build/toolset` directory. Before executing PathEnum, you need to covert the graph file into the format (CSR) that can be consumed by PathEnum, and generate the file containing queries.

### Input Graph Format

PathEnum works on directed and unlabeled graphs.
In the graph file, each line is a directed edge, which consits of the source vertex ID and the destination vertex ID.
    
Example:
    
```zsh
0 1
0 2
1 2
```

### Convert Graph
First, you need to convert the input graph to the format (CSR) that can be consumed by PathEnum with the command './EdgeList2DirectedCSR.out input_graph_folder output_graph_folder [skip_character]' where 'input_graph_folder' contains the input file and 'output_graph_folder' sets
the folder storing the converted graph. If you set 'skip_character', our tool will skip the line begin with that character. Note that before running the tool, you must change your graph file name to 'b_edge_list.bin'.

Example (Convert the input graph in '../../dataset/web-gooble/' folder to the CSR format, and store the output graph into '../../dataset/web-google/for_demo'. The skip character is '#'.):

```zsh
mkdir dataset/web-google/for_demo
cd build/toolset
./EdgeList2DirectedCSR.out ../../dataset/web-google/ ../../dataset/web-google/for_demo/ #
```

### Generate Queries

Next, you need to generate query pairs based on the graph with the command
'./GeneratedVertexPairs.out graph_folder query_file [pairs_number] [hot_veterx_percentage]' where 'graph_folder' sets the folder stored the converted graph and 'query_file' sets the file containing a number of queries. 'pairs_number' sets the number of queries to be generated. 'hot_vertex_percentage' sets the percentage of vertexes are hot based on the vertex degree. The tool will generate four kinds of query sets:
hot vertex to hot vertex, hot vretex to unhot vertex, unhot vertex to hot vertex, and unhot vertex to unhot vertex.

Example (Generate 1000 queries for the web-google graph. Set the top 10% vertexes in terms of the vertex degree as the hub vertex):

```zsh
./GenerateVertexPairs.out ../../dataset/web-google/for_demo/ ../../dataset/web-google/for_demo/ 1000 0.1
```
## Execute
After compiling the source code, you can find the binary file 'CycleEnumerator.out'
under the 'build/cycle' directory. Execute the binary with the following
command './CycleEnumerator graph_folder query_file result_folder method k per_query_time_limit [target_number_of_results]',
in which 'graph_folder' specifies the input folder containing the graph and 'query_file' specifies the
input of the queries. The 'result_folder' parameter sets the folder storing the logs. 'method' denotes the method evaluating queries, which can be 'IDX_DFS', 'IDX_JOIN' and 'PATH_ENUM'. 'k' is the length constraint. The 'per_query_time_limit' parameter configures the
time budget for the query. If the query cannot be completed within the time limit,
then the program will terminate the query and return the number of results found. The 'target_number_of_results' parameter sets the maximum number of
embeddings that you would like to find. If the number of embeddings enumerated
reaches the limit or all the results have been found, then the program will terminate. If you want to find all results, then omit the 'target_number_of_results' example. 

Example (Execute the query with the IDX_DFS method, and find all results.
The length constraint is 5 and the time limit for each query is 120 seconds. The logs are stored in the current directory):

```zsh
cd ../cycle
./CycleEnumerator.out ../../dataset/web-google/for_demo/ ../../dataset/web-google/for_demo/hot2hot_pairs.bin ./ "IDX_DFS" 5 120
```

## Output
We enumerate results without listing them explicitly. The metrics are listed as follows.

|Keyword                               | Description                                                 | Algorithms                     |
| :----------------------------------: | :---------------------------------------------------------: | :----------------------------: |
| Query time                           | the elapsed time on each query                              | All                            |
| Preprocess time                      | the elapsed time on the index construction and optimization | All                            |
| Full fledged estimation time         | the elapsed time on the full fledged estimation             | IDX_JOIN, PATH_ENUM            |
| Partial result count                 | the number of partial results generated                     | All                            |
| Result count                         | the number of results reported                              | All                            |
| Preliminary estimation result count  | the number of results estimated by preliminary estimation   | IDX_JOIN, PATH_ENUM            |
| Full fledged estimation result count | the number of results estimated by fledged estimation       | IDX_JOIN, PATH_ENUM            |


## Settings

In our paper, we report all results by default, while set the target number of results as 1000 when examining the response time.

## Datasets

You can download the real-world dataset used in our paper from [SNAP](http://snap.stanford.edu/data/) and [Network Repository](http://networkrepository.com/networks.php).
