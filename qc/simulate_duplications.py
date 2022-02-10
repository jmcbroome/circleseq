#script that takes a set of paired reads and simulates PCR duplication of the original reads by generating additional copies and rotating their content
#it saves the expanded fasta files and a logging file containing the names of duplicated reads and how much rotation was selected

import argparse
import numpy as np
import random

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threshold', type = int, help = 'Set a minimum number of times a base must be seen. default 2', default = 2)
    parser.add_argument('-i', '--input', help = 'prefixes of input fastqs. Assumes fastq file ends in "...R1.fq.gz".', required = True)
    parser.add_argument('-o', '--output', help = 'name of output fastqs. Default is the same as input with ".duplicated" inserted before the file extension', default = None)
    parser.add_argument('-s', '--stats', help = 'write names of reads duplicated and rotations used')
    args = parser.parse_args()
    return args

