from Bio import SeqIO
from Bio.Seq import Seq
import sys, os, re
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join
from collections import namedtuple, defaultdict
from tqdm import tqdm, tnrange
from joblib import Parallel, delayed
import multiprocessing
import gzip
from IPython.display import display
from operator import itemgetter
from datetime import datetime
import pyximport
pyximport.install()
from cimple_func import hamming
import distance


'''
check: is r1 is real r1?
(find primer in start of seq with shift)
'''
def is_r1(record, primer, shift, mist):
    for i in range(shift):
        if hamming(primer, str(record.seq[i : len(primer)+i]), mist):
            return (True)
    return (False)

def trim_primers(record, primers, re_parts, shift, mist):
    len_primer = len(primer)
    for i in range(shift):
        if hamming(primer, str(record.seq[i : len_primer + i]), mist):
            for elem in [primer, str(Seq(ad1).reverse_complement())]:
                if record.seq[len_primer+len(re_part)+i :].find(elem, 0) != -1:
                    return (info(is_good=False, read=None,
                                 alu_barcode=None,
                                 errors=np.array([0, 0, 0, 1, 0])))
            if record.seq[len_primer+len(re_part)+i :].find(restrict_site, 0) != -1:
                return (info(is_good=False, read=None,
                             alu_barcode=None,
                             errors=np.array([0, 0, 1, 0, 0])))
            record.description = ''
            record.name = ''
            re_part = str(record.seq[len_primer+i : len_primer+len(re_part)+i])
            record = record[len_primer+len(re_part)+i :]
            record = trimm_short_reads(str(Seq(target).reverse_complement()), record,
                                       mid_mist=mid_mist_short_reads[0], 
                                       end_mist=end_mist_short_reads[0],
                                       initial_pos=place_of_search_tail[0],
                                       save_pos=len(r2_start))
            try:
                point = poly_n_point(str(record.seq), 'T', poly_n_win_size_r1, poly_n_th_r1, poly_n_shift_r1, kind='head')
                record = record[point:]
            except:
                return (info(is_good=False, read=None,
                             alu_barcode=None,
                             errors=np.array([0, 0, 0, 0, 0])))
                print('wow')
            if record is None:
                return (info(is_good=False, read=None,
                             alu_barcode=None,
                             errors=np.array([0, 0, 0, 0, 0])))                
            alu_bar = '__pr12bq:' + re_part + '__pr12bq:' + str(record.seq)
            return (info(is_good=True, read=record,
                         alu_barcode=alu_bar,
                         errors=np.array([0, 0, 0, 0, 0])))
    return (info(is_good=False, read=None,
                 alu_barcode=None,
                 errors=np.array([1, 0, 0, 0, 0])))
