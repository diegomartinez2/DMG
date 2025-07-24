#!/usr/bin/env python
# python program to check if a path exists
#if path doesn’t exist we create a new path
import os
try:
   os.makedirs("pythondirectory")
except FileExistsError:
   # directory already exists
   pass
