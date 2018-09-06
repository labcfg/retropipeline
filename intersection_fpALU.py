import pandas as pd
import intervaltree as it
from collections import defaultdict
from collections import Counter
import sys, os, re
from os import listdir
from os.path import isfile, join
from tqdm import tqdm_notebook, tnrange
import distance
from operator import itemgetter
import numpy as np
from datetime import datetime
from joblib import Parallel, delayed



def intersection(filename, inputdir, outputdir, replib_inputdir, inswindow):

    readsname = os.path.splitext(filename)[0].split('_humanread')[0]

    rep = pd.read_table(replib_inputdir + readsname + '_ematch.txt', '\t')
    replib_group = rep.groupby(['CHR', 'STRAND'])
    replib_tree = {}
    for name, group in tqdm_notebook(replib_group, desc='rep: '+readsname):
        if name[1] == '+':
            start_group = [pos-int((pos-r_site)*0.9)
                             for r_site, pos in zip(list(group['CTAG']), list(group['START']))]
            end_group = [pos+inswindow+1 for pos in list(group['START'])]
        else:
            end_group = [pos+int((r_site-pos)*0.9)+1
                             for r_site, pos in zip(list(group['CTAG']), list(group['END']))]
            start_group = [pos-inswindow for pos in list(group['END'])]
        replib_tree[name[0] + name[1]] = it.IntervalTree(it.Interval(start, end)
         for start, end in zip(start_group, end_group))

    df_out = open(outputdir + readsname + '.txt', 'w')
    df_out.write('\t'.join(['CLUSTER_ID',
                            'READNAME',
                            'CHR',
                            'POS',
                            'INS_STRAND',
                            'RE',
                            'RE_AMOUNT',
                            'RE_HAMMING',
                            'R1',
                            'R2',
                            'TLEN',
                            'CIGAR_R1',
                            'MDFLAG_R1',
                            'MD_SUM',
                            'MISMATCH',
                            'INSERTION',
                            'DELETION',
                            'NUM_READS',
                            'NUM_BC']) + '\n')
    df = pd.read_table(inputdir + filename, '\t')
    for i in tqdm_notebook(range(np.shape(df)[0]), desc=readsname):
        row = df.iloc[i, ]
        if row['CHR'] + row['INS_STRAND'] in replib_tree:
            if len(replib_tree[row['CHR'] + row['INS_STRAND']][row['POS']]) == 0:
                pd.DataFrame(row).T.to_csv(df_out, sep='\t', mode='a', header=None, index=None)
            else:
                pass
        else:
            pd.DataFrame(row).T.to_csv(df_out, sep='\t', mode='a', header=None, index=None)
    df_out.close()


def main(inputdir, outputdir, replib_inputdir, inswindow, n_core):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'
    replib_inputdir = os.path.abspath(replib_inputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]
    onlyfiles = [f for f in onlyfiles if re.search('humanread', f)]

    if len(onlyfiles) == 1:
        filename = onlyfiles[0]
        stat_series = intersection(filename,
                              inputdir, outputdir, replib_inputdir, inswindow)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(intersection)(filename,
                                            inputdir, outputdir, replib_inputdir, inswindow)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()
