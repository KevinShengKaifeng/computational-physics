import matplotlib.pyplot as plt
import time


class Ran(object):

    def __init__(self, j=time.time()):
        j = int(j) % 2**64
        self.v = 4101842887655102017
        self.w = 1
        self.u = j ^ self.v
        self.int_rand()
        self.v = self.u
        self.int_rand()
        self.w = self.v
        self.int_rand()

    def int_rand(self):
        self.u = (self.u * 2862933555777941757 + 7046029254386353087) & 0xffffffffffffffff
        self.v ^= self.v >> 17
        self.v ^= (self.v << 31) & 0xffffffffffffffff
        self.v ^= self.v >> 8
        self.w = (4294957665*(self.w & 0xffffffff) + (self.w >> 32)) & 0xffffffffffffffff
        x = self.u ^ (self.u << 21) & 0xffffffffffffffff
        x ^= x >> 35
        x ^= (x << 4) & 0xffffffffffffffff
        return ((x + self.v) ^ self.w) & 0xffffffffffffffff

    def float_rand(self):
        return self.int_rand()/2**64


def statistical_test():
    ran = Ran()
    x = []
    y = []
    for i in range(1000):
        x.append(ran.float_rand())
        y.append(ran.float_rand())
    plt.plot(x, y, 'r.')
    plt.show()


if __name__ == "__main__":
    statistical_test()
