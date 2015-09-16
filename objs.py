from bitarray import bitarray
import mmh3

from utils import timeit, pretty_usage


class BloomFilterBuilder(object):
    """adapted from http://www.maxburstein.com/blog/creating-a-simple-bloom-filter/"""

    def __init__(self, size, hash_count):
        self.size = size
        self.hash_count = hash_count
        self.bit_array = bitarray(size)
        self.bit_array.setall(0)
        self.num_uniq_elems = 0
        self.fpr = None

    def _hash(self, string, seed):
        return mmh3.hash(string, seed)
        # return hash(string + str(seed))

    def add_uniq(self, string):
        """
        only add the string when it's new, compared to using contains, then add,
        add_uniq doesn't calculate hash function more than hash_count times
        """
        is_new = False
        for seed in xrange(self.hash_count):
            result = self._hash(string, seed) % self.size
            if is_new:
                self.bit_array[result] = 1
            else:
                if self.bit_array[result] == 0:
                    # meaning it's new
                    is_new = True
                    self.bit_array[result] = 1
                # else: need to continue check for all hash_count
        if is_new:
            self.num_uniq_elems += 1

    def calc_fpr(self):
        """
        Calculate false positive rate

        ref: http://www.maxburstein.com/blog/creating-a-simple-bloom-filter/
        """
        self.fpr = (1 - (1 - 1 / float(self.size)) ** (self.hash_count * self.num_uniq_elems)) ** self.hash_count
        return self.fpr


class BloomFilter(object):
    def __init__(self, size, hash_count, buffer_):
        """
        :param buffer_: a blob from sqlite3
        """
        self.size = size
        self.hash_count = hash_count
        self.bit_array = bitarray()
        self.bit_array.frombytes(str(buffer_))


    def _hash(self, string, seed):
        return mmh3.hash(string, seed)
        # return hash(string + str(seed))

    # def add(self, string):
    #     for seed in xrange(self.hash_count):
    #         result = self._hash(string, seed) % self.size
    #         self.bit_array[result] = 1

    def __contains__(self, string):
        """lookup for string to see if it's in the bloomfilter"""
        for seed in xrange(self.hash_count):
            result = self._hash(string, seed) % self.size
            if self.bit_array[result] == 0:
                return False
        return True

    def __str__(self):
        # self.size / 8.: bit => bytes
        return '<bf size: {0} ({1}), hash_count: {2}>'.format(
            self.size, pretty_usage(self.size / 8.), self.hash_count)

    def __unicode__(self):
        return str(self)

    def __repr__(self):
        return str(self)

    # def to_file(self, filename):
    #     with open(filename, 'wb') as opf:
    #         self.bit_array.tofile(opf)
        

    # @classmethod
    # def from_file(cls, filename):
    #     bf = cls()
    #     with open(filename, 'rb') as inf:
    #         bf.bit_array = bitarray.bitarray().fromfile(inf)
