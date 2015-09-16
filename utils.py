import time
import math                     # numpy log doesn't take arbitray base
import numpy as np
from functools import update_wrapper


def decorator(d):
    "Make function d a decorator: d wraps a function fn."
    def _d(fn):
        return update_wrapper(d(fn), fn)
    update_wrapper(_d, d)
    return _d


@decorator
def memo(f):
    """Decorator that caches the return value for each call to f(args).
    Then when called again with same args, we can just look it up."""
    cache = {}
    def _f(*args):
        try:
            return cache[args]
        except KeyError:
            cache[args] = result = f(*args)
            return result
        except TypeError:
            # some element of args can't be a dict key
            return f(args)
    return _f


@decorator
def timeit(f):
    def new_f(*args, **kwargs):
        bt = time.time()
        r = f(*args, **kwargs)
        et = time.time()        
        print "time spent on {0}: {1:.6f}s".format(f.func_name, et - bt)
        return r
    return new_f


def pretty_usage(num):
    """convert file size to a pretty human readable format"""
    # http://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0 and num > -1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0
    return "%3.1f %s" % (num, 'PB')


def split(array, current_id=0, current_level=0):
    """split array recursively, keeping track of subarray ids, starting from 0"""
    n = len(array)
    if n == 1:
        return
    else:
        c1_id, c2_id = calc_children_id(current_id, current_level)
        current_level += 1

        i = n // 2
        subarray1 = array[:i]
        yield (c1_id, subarray1)

        subarray2 = array[i:]
        yield (c2_id, subarray2)

        for a in split(subarray1, c1_id, current_level):
            yield a
        for a in split(subarray2, c2_id, current_level):
            yield a
        return



def split_indexes(beg, end, current_id=0, current_level=0, n_branches=2):
    """only return indexes to save memory usage

    split array recursively, keeping track of subarray ids, starting from 0

    :param b: when it's 2, it's a binary tree
    """
    delta = end - beg
    if delta == 1:
        return
    else:
        c1_id, c2_id = calc_children_id(current_id, current_level)
        current_level += 1

        mid = delta // 2 + beg
        yield (c1_id, beg, mid)
        yield (c2_id, mid, end)

        for a in split_indexes(beg, mid, c1_id, current_level):
            yield a
        for a in split_indexes(mid, end, c2_id, current_level):
            yield a
        return 


# def split(array, subarrays):
#     n = len(array)
#     if n == 1:
#         return
#     else:
#         i = n / 2
#         subarray1 = array[:i]
#         subarrays.append(subarray1)
#         subarray2 = array[i:]
#         subarrays.append(subarray2)
#         split(subarray1, subarrays)
#         split(subarray2, subarrays)
#         return 


def calc_total_num_nodes_above(level, nbr=2):
    """
    :param nbr: number of branches
    """
    if level < 0:
        raise ValueError('level cannot be below 0: {0}'.format(level))

    if level == 0:              # there is no node above level 0
        return 0

    # # because it's above level
    # level -= 1
    # res = (nbr ** (level + 1) - 1) // (nbr - 1)

    # combining the above two commented lines since level doesn't change
    # effectively
    res = (nbr ** (level) - 1) // (nbr - 1)

    # print level, res
    return res


def calc_level(cid, nbr=2):
    """given an current id (cid), calculate its level"""
    b, n  = nbr, cid
    if n == 0:
        return 0
    l =  math.log((n + 1) * (b - 1) + 1, b) - 1
    if l == int(l):
        return l
    elif l > int(l):
        return int(l + 1)
    else:
        raise ValueError('l < int(l) with l equals {0}'.format(r))


def calc_ith_node(cid, level=None, nbr=2):
    """calculate this is the ith node at this level"""
    if level is None:
        level = calc_level(cid, nbr)
    res = cid - calc_total_num_nodes_above(level, nbr) + 1
    # print res
    return res


def calc_children_id(cid, current_level=None, nbr=2):
    if current_level is None:
        current_level = calc_level(cid, nbr)

    num_nodes_above = calc_total_num_nodes_above(current_level, nbr)
    # i: ith node at the current_level, +1 because of 0-based index
    ith = cid - num_nodes_above + 1

    num_nodes_at_current_level = nbr ** current_level
    nth_children = nbr * ith + (num_nodes_above + num_nodes_at_current_level - 1)
    children = [nth_children]
    for i in xrange(1, nbr):
        children.append(nth_children - i)
    children.reverse()
    return tuple(children)


def kmerize(seq, kmer_size):
    res = []
    for i in xrange(len(seq) + 1 - kmer_size):
        res.append(seq[i:i + kmer_size])
    return res
