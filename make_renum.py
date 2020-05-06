#!/usr/bin/env python
import os
from argparse import ArgumentParser
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import pandas as pd
import requests
from os.path import exists, basename

pdbfile = ""

fastaseq = {}

def check_and_print_pdb(count, residue_buffer, residue_letter):
    global pdbfile
    global fastaseq
  # Check that CA, N and C are present!def check_and_print_pdb( outid, residue_buffer )
    hasCA, hasN, hasC = False, False, False

    for line in residue_buffer:
        atomname = line[12:16]
        # Only add bb atoms if they have occupancy!
        occupancy = float(line[55:60])
        if atomname == " CA " and occupancy > 0.0:
            hasCA = True
        if atomname == " N  " and occupancy > 0.0:
            hasN = True
        if atomname == " C  " and occupancy > 0.0:
            hasC = True

  # if all three backbone atoms are present withoccupancy proceed to print the residue
    if hasCA and hasN and hasC:
        for line in residue_buffer:
            # add linear residue count
            newnum = '%4d ' % count
            line_edit = line[0:22] + newnum + line[27:]
            # write the residue line
            pdbfile = pdbfile + line_edit

    # finally print residue letter into fasta strea
        chain = line[21]
        try:
            fastaseq[chain] += residue_letter
        except KeyError:
            fastaseq[chain] = residue_letter
    # count up residue number
        count = count + 1
        return True
    return False


def open_pdb( filename ):
    '''Open the PDB given in the filename (or equivalent).
    If the file is not found, then try downloading it from the internet.

    Returns: (lines, filename_stem)
    '''
    if os.path.exists(filename):
        print("Found existing PDB file at", filename)

    stem = os.path.basename(filename)
    stem = stem[:-4]

    lines = open(filename, 'r').readlines()

    return lines, stem

def clean_pdb(pdbid):
    '''Cleans PDBs by removing extraneous information, converting residues
    names and renumbering.

    Parameters
    ----------
    pbdid : str
        PDB file should be specified with the .pdb file handle.

    Returns
    -------
    None

    It writes both the cleaned PDB and a fasta file of the cleaned
    sequence to files in the cwd.
    '''
    longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                  'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                  'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                  'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                  'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
    shit_stat_insres, shit_stat_altpos, shit_stat_misdns = False, False, False
    global pdbfile
    global fastaseq
    pdbfile = ""
    fastaseq = {}
    
    files_to_unlink = []
    lines, filename_stem = open_pdb( pdbid )
    oldresnum , count, residue_buffer, residue_letter = '   ', 1, [], ''


    for line in lines:

        if line.startswith('ENDMDL'): break  # Only take the first NMR model
        if len(line) > 21:
            if line[0:4] != "ATOM" and line[0:6] != 'HETATM':
                continue

            line_edit = line
            resn = line[17:20]

            # Only process residues we know are valid.
            if resn not in longer_names:
                continue

            resnum = line_edit[22:27]

            # Is this a new residue
            if resnum != oldresnum:
                if residue_buffer != []:  # is there a residue in the buffer ?
                    # Are there any missing densities in the residue buffer?
                    if not check_and_print_pdb(count, residue_buffer, residue_letter):
                        # if unsuccessful
                        shit_stat_misdns = True
                    else:
                        count = count + 1

                    residue_buffer = []

                residue_letter = longer_names[resn]

            oldresnum = resnum

            insres = line[26]
            if insres != ' ':
                shit_stat_insres = True

            altpos = line[16]
            if altpos != ' ':
                shit_stat_altpos = True
                if altpos == 'A':
                    line_edit = line_edit[:16]+' '+line_edit[17:]
                else:
                    # Don't take the second and following alternate locations
                    continue

            residue_buffer.append(line_edit)


    if residue_buffer != []: # is there a residue in the buffer ?
        if not check_and_print_pdb(count, residue_buffer, residue_letter):
            # if unsuccessful
            shit_stat_misdns = True
        else:
            count = count + 1

    flag_altpos = "---"
    if shit_stat_altpos:
        flag_altpos = "ALT"
    flag_insres = "---"
    if shit_stat_insres:
        flag_insres = "INS"
    flag_misdns = "---"
    if shit_stat_misdns:
        flag_misdns = "DNS"

    nres = len("".join(fastaseq.values()))

    flag_successful = "OK"

    if nres <= 0:
        flag_successful = "BAD"

    print(filename_stem, "%5d" % nres, flag_altpos,  flag_insres,  flag_misdns, flag_successful)

    if nres > 0:
        outfile = filename_stem + "_" + "clean" + ".pdb"

        with open(outfile, 'w') as outid:
            outid.write(pdbfile)
            outid.write("TER\n")

        fastaseq = ["".join(fastaseq.values())]
        with open(filename_stem + ".fasta", 'w') as handle:
            handle.write('>'+filename_stem+'\n')
            handle.writelines(fastaseq)
            handle.write('\n')

    if len(files_to_unlink) > 0:
        for file in files_to_unlink:
            os.unlink(file)



def download_and_clean(pdbid):
    '''
    Download original and cleaned pdb from pdbid. 
    '''

    pdburl = 'https://files.rcsb.org/download/{}.pdb'.format(pdbid)
    pdbfile = pdbid+".pdb"
    
    print("Getting {} from {} ... \n".format(pdbid, pdburl))

    r = requests.get(pdburl)
    with open(pdbfile, 'w') as file:
        file.write(r.text)
    
    print("Cleaning {} ... \n".format(pdbfile))

    clean_pdb(pdbfile)

def create_renum_table(pdbid):
    """
    Create renumbered table from original and cleaned pdb.
    """
    # Parse pdbs
    parser = PDBParser()
    pdb = parser.get_structure(id='pdb', file=f'{pdbid}.pdb')
    clean = parser.get_structure(id='renumbered', file=f'{pdbid}_clean.pdb')

    # Create pdb df
    res_id = list(r.id[1] for r in pdb.get_residues())
    res_name = list(seq1(r.resname) for r in pdb.get_residues())
    df = pd.DataFrame(dict(pos=res_id, res=res_name))
    
    # Create clean df
    get_res = lambda c : list(r.id[1] for r in c.get_residues())
    data = dict((c.id, get_res(c)) for c in clean.get_chains())
    df2 = pd.DataFrame(data)
    
    # Merge to create renumbered table
    renum_table = df.merge(df2, left_index=True, right_index=True)
    return renum_table

def renum_from_pdbid(pdbid):
    download_and_clean(pdbid)
    
    print("Creating renumbering table ... \n")
    renum_table = create_renum_table(pdbid)
    return renum_table


parser = ArgumentParser(description="""
    Downloads and cleans PDB file using PDBID and creates renumbering table.
    The renumbering table is used to map positions in the sequence of the 
    monomer to corresponding positions in the chains of the cleaned PDB.
    A renumbering table is necessary to create mutfiles or resfiles.
    """)

parser.add_argument('pdbid', help='PDBID of the structure that should be downloaded, cleaned and renumbered')

args = parser.parse_args()

if not os.path.isdir(args.pdbid):
    print("Creating directory: {}".format(args.pdbid))
    os.mkdir(args.pdbid)

os.chdir(args.pdbid)

df = renum_from_pdbid(args.pdbid)

df.to_csv(f'{args.pdbid}_renum.csv',index=False)
print("Wrote renumbering table.")


os.chdir('..')
