from bitarray import bitarray
import mmh3

from utils import timeit


class BloomFilterBuilder(object):
    """adapted from http://www.maxburstein.com/blog/creating-a-simple-bloom-filter/"""

    def __init__(self, size, hash_count):
        self.size = size
        self.hash_count = hash_count
        self.bit_array = bitarray(size)
        self.bit_array.setall(0)

    def _hash(self, string, seed):
        return mmh3.hash(string, seed)
        # return hash(string + str(seed))

    @timeit                     # on the scale of 0.000043s
    def add(self, string):
        for seed in xrange(self.hash_count):
            result = self._hash(string, seed) % self.size
            self.bit_array[result] = 1

    # def __contains__(self, string):
    #     """lookup for string to see if it's in the bloomfilter"""
    #     for seed in xrange(self.hash_count):
    #         result = self._hash(string, seed) % self.size
    #         if self.bit_array[result] == 0:
    #             return False
    #     return True

    # def to_file(self, filename):
    #     with open(filename, 'wb') as opf:
    #         self.bit_array.tofile(opf)
        

    # @classmethod
    # def from_file(cls, filename):
    #     bf = cls()
    #     with open(filename, 'rb') as inf:
    #         bf.bit_array = bitarray.bitarray().fromfile(inf)
