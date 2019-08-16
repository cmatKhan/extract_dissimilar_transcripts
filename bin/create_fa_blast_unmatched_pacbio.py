# -*- coding: utf-8 -*-
"""
name: create_fa_blast_unmatched_pacbio.py
purpose: extract a .fa of sequences not matched in blast results FOR PACBIO
input: .fa used to create blast results in question and the unique tx names of the transcripts matched to some sequence by blast -- PACBIO
output: .fa of unmatched tx
written by: chase mateusiak, chase.mateusiak@gmail.com, chase.mateusiak@ucsf.edu
credit:
date: 20190813
Note:
"""
from Bio import SeqIO
import pandas as pd
import sys
import os

def dfFromFa (path_to_fa):
# purpose: create dataframe from .fasta
    with open(path_to_fa) as fasta_file:  # Will close handle cleanly
        identifiers = []
        seq = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            identifiers.append(seq_record.id)
            seq.append(str(seq_record.seq))

    my_dict = {'id':identifiers, 'seq': seq}

    return(pd.DataFrame.from_dict(my_dict))

def removeRowsInList(df, col, a_list):
    return(df[~df[col].isin(a_list)])

# open .fa and create dataframe with columns: txids, sequence
fa_path = sys.argv[1]

fa_df = dfFromFa(fa_path)

# get list of unique names in blast results, stripping \n from each item in list

similar_txids_path = sys.argv[2]

with open(similar_txids_path) as unique_blast_res:
    uniq_blast_res = list(map(lambda txid: txid.strip(), unique_blast_res))

# create dataframe of records not in the blast output

unmatched_records = removeRowsInList(fa_df, 'id',uniq_blast_res)

cols = ['id', 'seq']

unmatched_records = unmatched_records[cols]

# get basename of blast_results_input
basename = os.path.basename(fa_path)
base_noext = os.path.splitext(basename)[0]

# print .fa of unmatched sequences
fa = open(base_noext+'_blast_unmatched.fa', 'x')
[fa.write('>'+name+'\n'+seq+'\n') for name,seq in zip(unmatched_records['id'], unmatched_records['seq'])]
fa.close()
