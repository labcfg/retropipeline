import pandas as pd
import numpy as np
from collections import defaultdict
from collections import Counter
import sys, os, re
from os import listdir
from os.path import isfile, join
from tqdm import tqdm_notebook, tnrange
from operator import itemgetter
from datetime import datetime


def main(inputdir, outputdir,
         re_hamming,
         flank_errors,
         rsite_name,
         repeat,
         m_primer, primer_name,
         m_re, re_name):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]

    for filename in onlyfiles:
        df = pd.read_table(inputdir + filename, '\t')
        if re_hamming:
            df = df[df['RE_HAMMING'] < re_hamming]
        if flank_errors:
            df['TMP'] = df['MISMATCH'].values + df['INSERTION'].values + df['DELETION'].values
            df = df[df['TMP'] < flank_errors]
            df = df.drop('TMP', 1)
        if rsite_name:
            df = df[df[rsite_name] == 0]
        if repeat:
            df = df[df['REPEAT_FUNC'] < repeat]
        if m_primer:
            df = df[df['MISS_'+primer_name+'_HAMMING'] > m_primer]
        if m_re:
            df = df[df['MISS_'+re_name+'_HAMMING'] > m_re]
        df.to_csv(outputdir + filename, index=None, sep='\t')
