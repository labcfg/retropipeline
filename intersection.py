import sys
import os
import re
import pandas as pd
import numpy as np
from Bio import pairwise2
from Bio.Seq import Seq
from collections import Counter, defaultdict
import tqdm
import pysam

refway = '../../reference/hg38.fa'
ref = pysam.Fastafile(refway)


def intersection(filename, inswindow):
    rep = pd.read_table(replib_inputdir + readsname + '_ematch.txt', '\t')
    replib_group = rep.groupby(['CHR', 'STRAND'])
    replib_tree = {}
    for name, group in tqdm_notebook(replib_group, desc='rep: '+readsname):
        """
        in original script (for ALU 5') was:
        if name[1] == '+':
        '+' - strand of ALU (5')

        In code below, '-' - strand of LINE (3')
        """
        if name[1] == '-':
            start_group = [pos-int((pos-r_site)*0.9)
                             for r_site, pos in zip(list(group['CTAG']), list(group['START']))]
            end_group = [pos+inswindow+1 for pos in list(group['START'])]
        else:
            end_group = [pos+int((r_site-pos)*0.9)+1
                             for r_site, pos in zip(list(group['CTAG']), list(group['END']))]
            start_group = [pos-inswindow for pos in list(group['END'])]
        replib_tree[name[0] + name[1]] = it.IntervalTree(it.Interval(start, end)
         for start, end in zip(start_group, end_group))



def create_db(x, seq, window):
    x.columns = ['CHR', 'START', 'END', 'STRAND', 'NAME']
    chromosome_name = ['chr{}'.format(i) for i in range(1, 22+1)] + ['chrX', 'chrY']
    x = x[x['CHR'].isin(chromosome_name)]
    seq_control = []
    seq_recompile = re.compile(seq)
    for i, row in tqdm.tqdm_notebook(x.iterrows(), total=x.shape[0]):
        try:
            row_info = [row['CHR'], row['START'], row['END'], row['STRAND'], row['NAME']]
            pos_left = int(row['START']) - window
            flank_left = ref.fetch(row['CHR'], pos_left-1, int(row['START'])-1).upper()
            fiter = list(seq_recompile.finditer(flank_left))
            if fiter:
                fiter_pos = fiter[-1].start()
                ctag_left_pos = int(row['START']) - (len(flank_left) - fiter_pos)
                row_info.append(ctag_left_pos)
            else:
                row_info.append(np.nan)
            pos_right = int(row['END']) + window
            flank_right = ref.fetch(row['CHR'], int(row['END']), pos_right).upper()
            fiter = list(seq_recompile.finditer(flank_right))
            if fiter:
                fiter_pos = fiter[0].start()
                ctag_right_pos = int(row['END']) + fiter_pos + 1
                row_info.append(ctag_right_pos)
            else:
                row_info.append(np.nan)
            seq_control.append(row_info)
        except:
            pass
    colnames = ['CHR', 'START', 'END', 'STRAND', 'NAME', 'CTAG_LEFT', 'CTAG_RIGHT']
    return pd.DataFrame(seq_control, columns=colnames)