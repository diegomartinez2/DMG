#! /usr/local/bin/python3.6
"""
Random number generatrion with Box-Muller algorithm
"""
# source: https://gist.github.com/komasaru/
import math
import random
import sys
import traceback


class RndnumBoxMuller:
    M     = 10        # Average
    S     = 2.5       # Standard deviation
    N     = 10000     # Number to generate
    SCALE = N // 100  # Scale for histgram

    def __init__(self):
        self.hist = [0 for _ in range(self.M * 5)]

    def generate_rndnum(self):
        """ Generation of random numbers """
        try:
            for _ in range(self.N):
                res = self.__rnd()
                self.hist[res[0]] += 1
                self.hist[res[1]] += 1
        except Exception as e:
            raise

    def display(self):
        """ Display """
        try:
            for i in range(0, self.M * 2 + 1):
                print("{:>3}:{:>4} | ".format(i, self.hist[i]), end="")
                for j in range(1, self.hist[i] // self.SCALE + 1):
                    print("*", end="")
                print()
        except Exception as e:
            raise

    def __rnd(self):
        """ Generation of random integers """
        try:
            r_1 = random.random()
            r_2 = random.random()
            x = self.S \
              * math.sqrt(-2 * math.log(r_1)) \
              * math.cos(2 * math.pi * r_2) \
              + self.M
            y = self.S \
              * math.sqrt(-2 * math.log(r_1)) \
              * math.sin(2 * math.pi * r_2) \
              + self.M
            return [math.floor(x), math.floor(y)]
        except Exception as e:
            raise


if __name__ == '__main__':
    try:
        obj = RndnumBoxMuller()
        obj.generate_rndnum()
        obj.display()
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)
