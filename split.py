import os
import argparse
import sys

######################## Args Parser ####################################
parser = argparse.ArgumentParser(description='3D plot/PCA/kmeans clustering')
parser.add_argument("-f", '--file', help="file_name.csv", action="store")
args = parser.parse_args()

from itertools import chain, islice

def chunks(iterable, n):
   "chunks(ABCDE,2) => AB CD E"
   iterable = iter(iterable)
   while True:
       # store one line in memory,
       # chain it to an iterator on the rest of the chunk
       yield chain([next(iterable)], islice(iterable, n-1))

l = 500
file_large = args.file
with open(file_large) as bigfile:
    for i, lines in enumerate(chunks(bigfile, l)):
        file_split = '{}.{}'.format(file_large, i)
        with open(file_split, 'w') as f:
            f.writelines(lines)
