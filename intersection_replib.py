import pandas as pd
import numpy as np
import intervaltree as it
from collections import defaultdict
from collections import Counter
import sys, os, re
from os import listdir
from os.path import isfile, join
from tqdm import tqdm_notebook, tnrange, tqdm
import distance
from operator import itemgetter
import numpy as np
from datetime import datetime
from joblib import Parallel, delayed


get_repfunc = lambda d,o: (1-d/1000)*o
get_overlap = lambda x,y: max(0,min(x[1],y[1])-max(x[0],y[0]))

def intersection_replib(filename, inputdir, outputdir, replib_dict):

    readsname = os.path.splitext(filename)[0]

    df = pd.read_table(inputdir+filename)
    new_colums = list(df.columns)
    new_colums.extend(['REPEAT_FUNC', 'REPEAT_NAME'])

    df['REPEAT_FUNC'] = pd.Series(np.zeros(df.shape[0]), index = df.index)
    df['REPEAT_NAME'] = pd.Series(['*' for x in range(df.shape[0])], index = df.index)
    df_group = df.groupby(['CHR', 'INS_STRAND'])
    for name, group in tqdm_notebook(df_group, desc=readsname):
        if name[0] in replib_dict:
            if name[1] == '+':
                for p, t, idx in zip(group['POS'], group['TLEN'], group.index):
                    start, end = p-t, p+1
                    finter = replib_dict[name[0]][start:end]
                    if len(finter) > 0:
                        repfunc = []
                        repname_list = [x.data[1] for x in finter]
                        for x in finter:
                            overlap = get_overlap([x.begin, x.end], [start, end])/(end - start)
                            repfunc.append(get_repfunc(x.data[0], overlap))
                        repname = repname_list[repfunc.index(max(repfunc))]
                        df.set_value(idx, 'REPEAT_FUNC', max(repfunc))
                        df.set_value(idx, 'REPEAT_NAME', repname)
            else:
                for p, t, idx in zip(group['POS'], group['TLEN'], group.index):
                    start, end = p, p+t+1
                    finter = replib_dict[name[0]][start:end]
                    if len(finter) > 0:
                        repfunc = []
                        repname_list = [x.data[1] for x in finter]
                        for x in finter:
                            overlap = get_overlap([x.begin, x.end], [start, end])/(end - start)
                            repfunc.append(get_repfunc(x.data[0], overlap))
                        repname = repname_list[repfunc.index(max(repfunc))]
                        df.set_value(idx, 'REPEAT_FUNC', max(repfunc))
                        df.set_value(idx, 'REPEAT_NAME', repname)

    df.to_csv(outputdir + filename, sep='\t', index=None)


def main(inputdir, outputdir, repeatway, n_core):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'
    replib_inputdir = os.path.abspath(repeatway) + '/'

    autosomeXY = list(range(1, 23))
    autosomeXY.append('X')
    autosomeXY.append('Y')
    autosomeXY = ['chr' + str(x) for x in autosomeXY]

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]

    replib = pd.read_table(repeatway)
    #replib = replib[replib['CHR'] in autosomeXY]
    replib_group = replib.groupby(['CHR'])
    replib_dict = {}
    for name, group in tqdm_notebook(replib_group, desc='replib'):
        if name in autosomeXY:
            start_group = np.array(group['START'])
            end_group = np.array(group['END'])+1
            data_group = [(x,y) for x,y in zip(list(group['DIV']), list(group['NAME']))]
            replib_dict[name] = it.IntervalTree(it.Interval(begin, end, data)
             for begin, end, data in zip(start_group, end_group, data_group))

    if len(onlyfiles) == 1:
        filename = onlyfiles[0]
        stat_series = intersection_replib(filename,
                              inputdir, outputdir, replib_dict)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(intersection_replib)(filename,
                                            inputdir, outputdir, replib_dict)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()
