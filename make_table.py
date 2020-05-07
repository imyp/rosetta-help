#!/usr/bin/env python
import pandas as pd
from argparse import ArgumentParser
import sys
import os
import textwrap

# Setup arguments
parser = ArgumentParser(description="Create mutation table from a textfile containing mutations for various mutants.")

parser.add_argument('prefix',
        help="Prefix of the mutant ID that will be prepended to a number.")
parser.add_argument('mut_list',
        help="""
        Path to text file containing mutations from various mutants. Each line contains 
        mutations from a single mutant and mutations are comma-separated and written in 
        the following format: A24L
        """)
parser.add_argument('mut_table',
        help="""
        Name of the table containing mutations. This table will be saved to the same
        directory as the mutation list.
        """)

args = parser.parse_args()

out_dir, in_name = os.path.split(args.mut_list)
out_path = os.path.join(out_dir, args.mut_table)

# Read the mut_list
df = pd.read_csv(args.mut_list, sep=';', names=['mutations'])

# Find the number of digits needed to enumerate all mutants.
lines = len(df)
digits = len(str(lines))

# Append a zero-filled enumeration to prefix to create an id.
def make_id(line):
    num = str(line.name)
    znum = num.zfill(digits)
    return args.prefix+znum

# Apply to all mutants to make id column.
df['id'] = df.apply(make_id, axis=1)

df.to_csv(out_path, sep=';', index=False)
