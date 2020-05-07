#!/usr/bin/env python
"""Creates mutfiles for the mutants in `mut_table` with the correct numbering
for the desired `pdbid` and `chains`. 

This is done by using the renumbering table in the `structure`
directory that corresponds to the chosen `pdbid`.

If no renumbering table exists for your `pdbid`, you can create one
with the `make_renum.py` script found in the `structure` directory.

Resfiles will be written to the following directory: 

resfiles/<pdbid>/<mut_table_dir>/<out_dir>/

Where <pdbid> is the chosen value for `pdbid`, <mut_table_dir> is the
name of the directory your `mut_table` resides in, and <out_dir> is
your chosen value for `out_dir`.
"""
import os
import re
from os.path import split, join, isdir
import pandas as pd
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def make_mutfile(mutstr, chains, renum):
    """
    Create mutfile from mutation string, chains and renumbering table.
    
    Parameters
    ----------  
    mutstr : str
        String containing comma-separated mutations for a mutants.
    chains : str
        Letters of the chains where mutations are taking place.
    renum : pandas.DataFrame
        The renumbering table with numbering for chains.
    """
    
    # Split string to individual mutations
    muts = re.split(',\s*', mutstr)
    
    # Make chain list
    chains = list(chains)
    
    # Create mutfile header
    n = len(chains)*len(muts)
    mutfile = [f'total {n}\n{n}']

    for chain in chains:
        for mut in muts:
            
            # Split mutation into old and new residue, as well as number
            mut_pattern = '^(?P<old>\D)(?P<num>\d+)(?P<new>\D)$' 
            m = re.match(mut_pattern, mut)
            
            # renumber according to chain
            num = renum[renum.pos==int(m.group('num'))][chain].values[0]
            
            # append mutline to mutfile
            mutline = f"{m.group('old')} {num} {m.group('new')}"
            mutfile.append(mutline)
    return '\n'.join(mutfile)



def write_mutfile(name, mutstr, chains, renum):
    """
    Write mutfile specified by the arguments given.

    Parameters
    ----------

    name : str
        Name of file that is going to be written.
    mutstr : str
        Mutant string with mutations separated by commas.
    chains : str
        A string containing the names of the chains where mutation occurs.
    renum : pandas.DataFrame
        Renumbering data frame with information of structure.
    """
    with open(name, 'w') as f:
        f.write(make_mutfile(mutstr, chains, renum))



parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

parser.add_argument('mut_table',
        help = """a CSV file generated using the `make_table.py`-script in the 
        `mut_tables`-directory """)

parser.add_argument('pdbid',
        help="""
        The PDBID of the structure that is going to be used in Rosetta 
        protocols. This argument makes the program choose the correct 
        renumbering table for translating mutations to numbering found in 
        in the cleaned PDB structures.
        """)

parser.add_argument('chains',
        help="""
        A string containing the chains where mutations are going to be placed.
        If you want mutations to happen in chains A and L of a PDB structure, 
        use the string "AL" as an argument.
        """)

parser.add_argument('out_dir',
        help="Name of output directory that mutants are going to be saved to.")

args = parser.parse_args()

renum_path = "structures/{pdbid}/{pdbid}_renum.csv".format(pdbid=args.pdbid)

print("Reading mutation table: {} ...\n".format(args.mut_table))
mut_df = pd.read_csv(args.mut_table, sep=';')

print("Reading renumbering table: {} ...\n".format(renum_path))
renum_df = pd.read_csv(renum_path)

table_dir_path = split(args.mut_table)[0]
table_dir_name = split(table_dir_path)[1]


out_path = join('mutfiles',join(args.pdbid,join(table_dir_name, args.out_dir)))
if not isdir(out_path):
    print("Creating output directory: {} ... \n".format(out_path))
    os.makedirs(out_path)
else:
    print("Found output directory: {} ...\n".format(out_path))

def apply_write_mut(line):
    "Write mutfile when given a line from a mutation table"
    filename = join(out_path, line.id+".mut")
    mutstr = line.mutations
    return write_mutfile(filename, mutstr, args.chains, renum_df)
    

print('Writing resfiles to output directory...\n')

mut_df.apply(apply_write_mut, axis=1)

print("Wrote files.")
