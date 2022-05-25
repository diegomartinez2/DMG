import multiprocessing
from functools import partial

# Manager to create shared object.
manager = multiprocessing.Manager()

# Create a global variable.
dictionary = manager.dict()


# Each process will run this function.
def fetch(lock, item):
    # Do cool stuff
    lock.acquire()
    dictionary[item] = 0
    # Write to stdout or logfile, etc.
    lock.release()


if __name__ == '__main__':

    # Each item will map to function.
    my_list = range(5)

    # Create a lock.
    l = manager.Lock()

    # Multiprocess pool.
    pool = multiprocessing.Pool(processes=4)

    func = partial(fetch, l)

    # Map each item to function
    pool.map(func, my_list)

    # Close pool.
    pool.close()

    # Wait for all thread.
    pool.join()

    # Print meta_info
    print(dictionary)
