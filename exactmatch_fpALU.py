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


def intersection(filename, inputdir, outputdir, replibrary, min_read, window):

    rep = pd.read_table(replibrary, compression='bz2')
    rep.columns = ['CHR', 'START', 'END', 'CTAG', 'STRAND', 'NAME']
    replib_group = rep.groupby(['CHR', 'STRAND'])
    replib_tree = {}
    for name, group in tqdm_notebook(replib_group, desc='replib'):
        if name[1] == '+':
            start_group = [pos-window for pos in list(group['START'])]
            end_group = [pos+1+window for pos in list(group['START'])]
        else:
            end_group = [pos+2+window for pos in list(group['END'])]
            start_group = [pos+1-window for pos in list(group['END'])]
        group_list = map(list, group.values)
        group_list = ['\t'.join([x[0], str(x[1]), str(x[2]), str(x[3]), x[4], x[5]]) for x in group_list]
        group_list = [(x,y) for x,y in zip(group_list, list(group.index))]
        replib_tree[name[0] + name[1]] = it.IntervalTree(it.Interval(start, end, data)
         for start, end, data in zip(start_group, end_group, group_list))

    readsname = os.path.splitext(filename)[0].split('_humanread')[0]

    replib_out = open(outputdir + readsname + '_ematch.txt', 'w')
    replib_out.write('\t'.join(['CHR',
                                'START',
                                'END',
                                'CTAG',
                                'STRAND',
                                'NAME',
                                'READNAME',
                                'IDX']) + '\n')
    df = pd.read_table(inputdir + filename, '\t')
    for i in tqdm_notebook(range(np.shape(df)[0]), desc=readsname):
        row = df.iloc[i, ]
        if row['NUM_READS'] >= min_read:
            if row['CHR'] + row['INS_STRAND'] in replib_tree:
                finter = replib_tree[row['CHR'] + row['INS_STRAND']][row['POS']]
                if len(finter) > 0:
                    for x in finter:
                        replib_out.write(x.data[0] + '\t' + row['READNAME'] + '\t' + str(x.data[1]) + '\n')
                        replib_tree[row['CHR'] + row['INS_STRAND']].removei(x.begin, x.end, x.data)
    replib_out.close()


def main(inputdir, outputdir, replibrary, min_read, n_core, window):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]
    onlyfiles = [f for f in onlyfiles if re.search('humanread', f)]

    if len(onlyfiles) == 1:
        filename = onlyfiles[0]
        stat_series = intersection(filename,
                              inputdir, outputdir, replibrary, min_read, window)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(intersection)(filename,
                                            inputdir, outputdir, replibrary, min_read, window)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()
