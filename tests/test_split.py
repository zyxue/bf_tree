import numpy as np

import utils as U


def test_calc_total_num_nodes_above():
    assert U.calc_total_num_nodes_above(0) == 0
    assert U.calc_total_num_nodes_above(1) == 1
    assert U.calc_total_num_nodes_above(2) == 3
    assert U.calc_total_num_nodes_above(3) == 7


def test_calc_total_num_nodes_above_quaternary():
    assert U.calc_total_num_nodes_above(0, 4) == 0
    assert U.calc_total_num_nodes_above(1, 4) == 1
    assert U.calc_total_num_nodes_above(2, 4) == 5
    assert U.calc_total_num_nodes_above(3, 4) == 21


def test_calc_level():
    assert U.calc_level(0) == 0
    assert U.calc_level(1) == 1
    assert U.calc_level(2) == 1
    assert U.calc_level(3) == 2
    assert U.calc_level(4) == 2
    assert U.calc_level(5) == 2
    assert U.calc_level(6) == 2
    assert U.calc_level(7) == 3
    assert U.calc_level(8) == 3
    assert U.calc_level(9) == 3
    assert U.calc_level(10) == 3
    assert U.calc_level(11) == 3
    assert U.calc_level(12) == 3
    assert U.calc_level(13) == 3
    assert U.calc_level(14) == 3
    assert U.calc_level(15) == 4


def test_calc_level_quaternary():
    assert U.calc_level(0, 4) == 0
    assert U.calc_level(1, 4) == 1
    assert U.calc_level(2, 4) == 1
    assert U.calc_level(3, 4) == 1
    assert U.calc_level(4, 4) == 1
    assert U.calc_level(5, 4) == 2
    assert U.calc_level(6, 4) == 2
    assert U.calc_level(7, 4) == 2
    assert U.calc_level(8, 4) == 2
    assert U.calc_level(9, 4) == 2
    assert U.calc_level(10, 4) == 2
    assert U.calc_level(11, 4) == 2
    assert U.calc_level(12, 4) == 2
    assert U.calc_level(13, 4) == 2
    assert U.calc_level(14, 4) == 2
    assert U.calc_level(15, 4) == 2
    assert U.calc_level(21, 4) == 3


def test_ith_node():
    assert U.calc_ith_node(0, 0) == 1
    assert U.calc_ith_node(1, 1) == 1
    assert U.calc_ith_node(2, 1) == 2
    assert U.calc_ith_node(3, 2) == 1
    assert U.calc_ith_node(4, 2) == 2
    assert U.calc_ith_node(5, 2) == 3
    assert U.calc_ith_node(6, 2) == 4
    assert U.calc_ith_node(7, 3) == 1
    assert U.calc_ith_node(13, 3) == 7


def test_ith_node_quaternary():
    assert U.calc_ith_node(0, 0, 4) == 1
    assert U.calc_ith_node(1, 1, 4) == 1
    assert U.calc_ith_node(2, 1, 4) == 2
    assert U.calc_ith_node(3, 1, 4) == 3
    assert U.calc_ith_node(4, 1, 4) == 4
    assert U.calc_ith_node(5, 2, 4) == 1
    assert U.calc_ith_node(6, 2, 4) == 2
    assert U.calc_ith_node(7, 2, 4) == 3
    assert U.calc_ith_node(13, 2, 4) == 9

    # with no level information
    assert U.calc_ith_node(0, nbr=4) == 1
    assert U.calc_ith_node(1, nbr=4) == 1
    assert U.calc_ith_node(2, nbr=4) == 2
    assert U.calc_ith_node(3, nbr=4) == 3
    assert U.calc_ith_node(4, nbr=4) == 4
    assert U.calc_ith_node(5, nbr=4) == 1
    assert U.calc_ith_node(6, nbr=4) == 2
    assert U.calc_ith_node(7, nbr=4) == 3
    assert U.calc_ith_node(13, nbr=4) == 9


def test_calc_children_id():
    # binary
    assert U.calc_children_id(0, 0) == (1, 2)
    assert U.calc_children_id(2, 1) == (5, 6)
    assert U.calc_children_id(3, 2) == (7, 8)
    assert U.calc_children_id(5, 2) == (11, 12)
    assert U.calc_children_id(6, 2) == (13, 14)
    assert U.calc_children_id(7, 3) == (15, 16)

    # quarternary
    assert U.calc_children_id(0, 0, 4) == (1, 2, 3, 4)
    assert U.calc_children_id(1, 1, 4) == (5, 6, 7, 8)
    assert U.calc_children_id(2, 1, 4) == (9, 10, 11, 12)
    assert U.calc_children_id(5, 2, 4) == (21, 22, 23, 24)
    assert U.calc_children_id(6, 2, 4) == (25, 26, 27, 28)


def test_split():
    subarrays = list(U.split(range(8)))
    assert (1, [0, 1, 2, 3]) in subarrays
    assert (2, [4, 5, 6, 7]) in subarrays
    assert (7, [0]) in subarrays
    assert (8, [1]) in subarrays
    assert (9, [2]) in subarrays
    assert (10, [3]) in subarrays
    assert (11, [4]) in subarrays
    assert (12, [5]) in subarrays
    assert (13, [6]) in subarrays
    assert (14, [7]) in subarrays


def test_split_in_memoryview():
    """memoryview won't work with array of strings"""
    memview = memoryview(np.array(range(8)))
    f = lambda x: (np.array(x[0]).tolist(), x[1])
    subarrays = list(map(f, U.split(memview)))
    assert (7, [0]) in subarrays
    assert (8, [1]) in subarrays
    assert (9, [2]) in subarrays
    assert (10, [3]) in subarrays
    assert (11, [4]) in subarrays
    assert (12, [5]) in subarrays
    assert (13, [6]) in subarrays
    assert (14, [7]) in subarrays


def test_split_indexes():
    """memoryview won't work with array of strings"""
    # e.g. splitting an array with 10 elements
    subarrays = list(U.split_indexes(0, 10))
    # for __ in sorted(subarrays):
    #     print __
    assert (7, 0, 1) in subarrays
    assert (8, 1, 2) in subarrays
    assert (9, 2, 3) in subarrays
    assert (10, 3, 5) in subarrays
    assert (11, 5, 6) in subarrays
    assert (12, 6, 7) in subarrays
    assert (13, 7, 8) in subarrays
    assert (14, 8, 10) in subarrays

    # quaternary
    subarrays = list(U.split_indexes(0, 15, 4))
    for __ in sorted(subarrays):
        print __
    assert (1, 0, 3) in subarrays
    assert (2, 3, 6) in subarrays
    assert (3, 6, 9) in subarrays
    assert (4, 9, 15) in subarrays
    assert (5, 0, 1) in subarrays
    assert (6, 1, 2) in subarrays
    assert (7, 2, 3) in subarrays
    assert (9, 3, 4) in subarrays
    assert (10, 4, 5) in subarrays
    assert (11, 5, 6) in subarrays
    assert (13, 6, 7) in subarrays
    assert (14, 7, 8) in subarrays
    assert (15, 8, 9) in subarrays
    assert (17, 9, 10) in subarrays
    assert (18, 10, 11) in subarrays
    assert (19, 11, 12) in subarrays
    assert (20, 12, 15) in subarrays
    assert (81, 12, 13) in subarrays
    assert (82, 13, 14) in subarrays
    assert (83, 14, 15) in subarrays


test_calc_total_num_nodes_above()
test_calc_total_num_nodes_above_quaternary()
test_calc_level()
test_calc_level_quaternary()
test_ith_node()
test_ith_node_quaternary()
test_calc_children_id()
test_split_indexes()
test_split()
