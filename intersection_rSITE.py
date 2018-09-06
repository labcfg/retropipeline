import pandas as pd
import numpy as np
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


def intersection_rsite(filename, inputdir, outputdir, restrict_dict, restrictase, inswindow):

    readsname = os.path.splitext(filename)[0]

    df = pd.read_table(inputdir + filename, '\t')
    df[restrictase] = pd.Series(np.zeros(df.shape[0]), index = df.index)
    df_group = df.groupby(['CHR', 'INS_STRAND'])
    for name, group in tqdm_notebook(df_group, desc=readsname):
        if name[0] in restrict_dict:
            if name[1] == '+':
                insrange = zip(np.array(group['POS'])-np.array(group['TLEN'])+3,
                               np.array(group['POS'])+inswindow+1,
                               list(group.index))
            else:
                insrange = zip(np.array(group['POS'])-inswindow,
                               np.array(group['POS'])+np.array(group['TLEN']-2),
                               list(group.index))
            for start, end, idx in insrange:
                if len(restrict_dict[name[0]][start:end]) > 0:
                    df.set_value(idx, restrictase, 1)

    df.to_csv(outputdir + filename, sep='\t', index=None)


def main(inputdir, outputdir, restrictway, restrictase, inswindow, n_core):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]
    restrict_dict = {}
    restrict = pd.read_table(restrictway, compression='bz2')
    restrict.columns = ['CHR', 'POS']
    restrict_group = restrict.groupby(['CHR'])
    for name, group in  tqdm_notebook(restrict_group, desc='restrict'):
        start_group = np.array(group['POS'])
        end_group = start_group+1
        restrict_dict[name] = it.IntervalTree(it.Interval(start, end)
         for start, end in zip(start_group, end_group))

    if len(onlyfiles) == 1:
        filename = onlyfiles[0]
        stat_series = intersection_rsite(filename,
                              inputdir, outputdir, restrict_dict, restrictase, inswindow)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(intersection_rsite)(filename,
                                            inputdir, outputdir, restrict_dict, restrictase, inswindow)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()