#!/usr/bin/env python
#another parser for imput from the terminal
import argparse

argParser = argparse.ArgumentParser()
argParser.add_argument("-pop", "--POPULATION", help="your name")
argParser.add_argument("-i", "--int", type=int, help="your numeric age ")
args = argParser.parse_args()
print("args=%s" % args)

print("args.POPULATION=%s" % args.POPULATION)
