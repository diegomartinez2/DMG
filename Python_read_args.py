#!/usr/bin/env python
import argparse
import sys
import subprocess
print("this is script.py")

def parse_args():
    parser=argparse.ArgumentParser(description="a script to do stuff")
    parser.add_argument("R1_file")
    parser.add_argument("--threads", type=int, default=16)
    args=parser.parse_args()

    print("the inputs are:")
    for arg in vars(args):
        print("{} is {}".format(arg, getattr(args, arg)))

    return args

def stats(file):
    cmd=["seqkit", "stats", file]
    subprocess.run(cmd)

def process_files(things):
    print("this is the process files function")
    print(things.R1_file)
    # cmd=["seqkit", "stats", things.R1_file]
    # subprocess.run(cmd)
    stats(things.R1_file)

def main():
    print("this is the main function")
    inputs=parse_args()
    print(inputs.R1_file)
    print(inputs.threads)
    process_files(inputs)

if __name__ == '__main__':
    main()
