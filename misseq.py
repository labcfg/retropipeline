import pandas as pd
import numpy as np
import pysam
from Bio import pairwise2
from Bio.Seq import Seq
import distance
from collections import defaultdict
from collections import Counter
import sys, os, re
from os import listdir
from os.path import isfile, join
from tqdm import tqdm_notebook, tnrange
from operator import itemgetter
from datetime import datetime
from joblib import Parallel, delayed


def misseq(filename, inputdir, outputdir, refway, mseq, mname, shift):

    df = pd.read_table(inputdir + filename)
    ref = pysam.Fastafile(refway)

    readsname = os.path.splitext(filename)[0]

    idx = list(df.index)
    if mseq:
        mseq_list = [mseq for x in range(df.shape[0])]
    else:
        mseq_list = list(df['RE'])
    df[mname] = pd.Series(np.repeat(len(mseq_list[0]),df.shape[0]),index = df.index)
    if mseq not in ['AGCT', 'CTAG']:
        for i in tqdm_notebook(range(df.shape[0]), desc=readsname):
            row = df.iloc[i, ]
            miss_seq = mseq_list[i]
            pos = int(row['POS'])
            if row['INS_STRAND'] == '+':
                seq = Seq(ref.fetch(row['CHR'], pos, pos+shift)).reverse_complement()
                seq = str(seq).upper()
            else:
                seq = ref.fetch(row['CHR'], pos-shift-1, pos-1)
                seq = seq.upper()
            score = []
            for j in range(len(seq) - len(mseq_list[i]) + 1):
                score.append(distance.hamming(seq[j:len(mseq_list[i])+j], mseq_list[i]))
            df.set_value(idx[i], mname, min(score))
    else:
        for i in tqdm_notebook(range(df.shape[0]), desc=readsname):
            row = df.iloc[i, ]
            miss_seq = mseq_list[i]
            pos = int(row['POS'])
            if row['INS_STRAND'] == '+':
                seq = Seq(ref.fetch(row['CHR'], pos-int(row['TLEN']), pos+shift))
            else:
                seq = Seq(ref.fetch(row['CHR'], pos-shift-1, pos+int(row['TLEN'])-1-2))
            if seq.find(miss_seq, 0) != -1:
                df.set_value(idx[i], mname, 1)
            else:
                df.set_value(idx[i], mname, 0)

    df.to_csv(outputdir + filename, sep='\t', index=None)


def main(inputdir, outputdir, refway, mseq, mname, shift, n_core):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]

    if len(onlyfiles) == 1:
        filename = onlyfiles[0]
        stat_series = misseq(filename,
                        inputdir, outputdir, refway, mseq, mname, shift)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(misseq)(filename,
                                            inputdir, outputdir, refway, mseq, mname, shift)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()
