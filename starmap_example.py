#!/usr/bin/env python3
from multiprocessing.dummy import Pool as ThreadPool 

def write(i, x):
    print(i, "---", x)

a = ["1","2","3"]
b = ["4","5","6"]

pool = ThreadPool(2)
pool.starmap(write, zip(a,b))
pool.close()
pool.join()
