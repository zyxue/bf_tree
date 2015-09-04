import time
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


def calc_total_num_nodes_above(level):
    if level < 0:
        raise ValueError('level cannot be below 0: {0}'.format(level))

    if level == 0:
        return 0

    return 2 ** level - 1


def calc_ith_node(current_id, level):
    return current_id - calc_total_num_nodes_above(level) + 1


def calc_children_id(current_id, current_level):
    if current_level == 0 and current_id == 0:
        return (1, 2)

    num_nodes_above = calc_total_num_nodes_above(current_level)
    # i: ith node at the current_level, +1 because of 0-based index
    i = current_id - num_nodes_above + 1

    num_nodes_at_current_level = 2 ** current_level
    c2 = 2 * i + (num_nodes_above + num_nodes_at_current_level - 1)
    c1 = c2 - 1
    return c1, c2
