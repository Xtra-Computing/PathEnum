import numpy as np
import struct
import os
import sys
import shutil
import glob
from subprocess import Popen, PIPE


def read_array(file_path):
    with open(file_path, 'rb') as f:
        (element_size, array_length) = struct.unpack('2I', f.read(8))

        if element_size == 4:
            dtype = np.int32
        elif element_size == 8:
            dtype = np.int64
        else:
            print('The element size {0} cannot be supported'.format(element_size))
            exit(-1)
        res_arr = np.fromfile(f, dtype=dtype)

        if array_length != len(res_arr):
            print('The array length cannot match {0}:{1}'.format(array_length, len(res_arr)))
            exit(-1)
    return res_arr


def generate_args(binary, *params):
    arguments = [binary]
    arguments.extend(list(params))
    return arguments


def execute_binary(args):
    process = Popen(' '.join(args), shell=True, stdout=PIPE, stderr=PIPE)
    (std_output, std_error) = process.communicate()
    process.wait()
    rc = process.returncode

    return rc, std_output, std_error


def check_correctness(binary_path, graph_folder, query_file, result_folder, query_method, length_constraint,
                      correct_result_file):
    # Execute the binary.
    execution_args = generate_args(binary_path, graph_folder, query_file, result_folder, query_method, length_constraint, '120')
    (rc, std_output, std_error) = execute_binary(execution_args)
    if rc == 0:
        # Load correct result.
        expected_results = read_array(correct_result_file)

        # Load output result.
        result_file_path_pattern = '{0}/{1}-{2}-*-result_count.bin'.format(result_folder, query_method, length_constraint)
        output_results = read_array(glob.glob(result_file_path_pattern)[0])

        for query_id, query_result in enumerate(output_results):
            if query_result != expected_results[query_id]:
                print('Query method {0} with length constraint {1} on query {2} is error. Output {3},' \
                      'Expected {4}'.format(query_method, length_constraint, query_id,
                                            query_result, expected_results[query_id]))
                exit(-1)
    else:
        print('Query method {0} with length constraint {1} is error: {2} {3}'.format(query_method,
                                                                                     length_constraint, std_output, std_error))
        exit(-1)

    print('Query method {0} with length constraint {1} passes the test'.format(query_method, length_constraint))


if __name__ == '__main__':
    input_binary_path = sys.argv[1]
    if not os.path.isfile(input_binary_path):
        print('The binary {0} does not exist.'.format(input_binary_path))
        exit(-1)

    default_graph_folder = 'soc-Epinions1'
    default_test_queries_file = 'test_queries.bin'
    default_test_queries_correct_results_file = 'test_queries_length_4_correct_results.bin'

    dir_path = os.path.dirname(os.path.realpath(__file__))
    input_graph_folder = '{0}/test_dataset/{1}'.format(dir_path, default_graph_folder)
    input_query_file = '{0}/test_dataset/{1}/{2}'.format(dir_path, default_graph_folder, default_test_queries_file)
    input_result_folder = '{0}/test_dataset/{1}/output'.format(dir_path, default_graph_folder)
    input_correct_results_file = '{0}/test_dataset/{1}/{2}'.format(dir_path, default_graph_folder,
                                                                   default_test_queries_correct_results_file)

    # Create output folder.
    if not os.path.exists(input_result_folder):
        os.makedirs(input_result_folder)

    # Check correctness.
    query_method_list = ['IDX_DFS', 'IDX_JOIN', 'PATH_ENUM']
    for query_method in query_method_list:
        check_correctness(input_binary_path, input_graph_folder, input_query_file, input_result_folder, query_method, '4',
                          input_correct_results_file)

    # Remove output folder.
    if os.path.exists(input_result_folder):
        shutil.rmtree(input_result_folder)

    print('All test cases successfully complete.')
