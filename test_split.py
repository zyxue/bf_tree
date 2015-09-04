
from utils import calc_total_num_nodes_above, calc_ith_node, calc_children_id, split

def test_calc_total_num_nodes_above():
    assert calc_total_num_nodes_above(0) == 0
    assert calc_total_num_nodes_above(1) == 1
    assert calc_total_num_nodes_above(2) == 3
    assert calc_total_num_nodes_above(3) == 7


def test_ith_node():
    assert calc_ith_node(0, 0) == 1
    assert calc_ith_node(1, 1) == 1
    assert calc_ith_node(2, 1) == 2
    assert calc_ith_node(3, 2) == 1
    assert calc_ith_node(4, 2) == 2
    assert calc_ith_node(5, 2) == 3
    assert calc_ith_node(6, 2) == 4
    assert calc_ith_node(7, 3) == 1
    assert calc_ith_node(13, 3) == 7


def test_calc_children_id():
    assert calc_children_id(0, 0) == (1, 2)
    assert calc_children_id(2, 1) == (5, 6)
    assert calc_children_id(3, 2) == (7, 8)
    assert calc_children_id(5, 2) == (11, 12)
    assert calc_children_id(6, 2) == (13, 14)
    assert calc_children_id(7, 3) == (15, 16)


def test_split():
    subarrays = []
    split(range(8), subarrays)
    assert ([0], 7) in subarrays
    assert ([1], 8) in subarrays
    assert ([2], 9) in subarrays
    assert ([3], 10) in subarrays
    assert ([4], 11) in subarrays
    assert ([5], 12) in subarrays
    assert ([6], 13) in subarrays
    assert ([7], 14) in subarrays


test_calc_total_num_nodes_above()
test_ith_node()
test_calc_children_id()
test_split()


subarrays = []
split(range(10), subarrays)
for k, i in enumerate(sorted(subarrays, key=lambda x: x[1])):
    print k+ 1, i

