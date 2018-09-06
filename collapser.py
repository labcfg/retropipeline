import pandas as pd
import re
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


'''
get best re (retroelement part)
1. find max_amount of re_part
2. if max_amount not one return min hamming with target re_part
return(seq, amount, hamming)
'''
def get_best_re(x, target_re):
    x_count = sorted(dict(Counter(x)).items(), key=itemgetter(1), reverse=True)
    x_count_max = [(re,a) for re,a in x_count if a == x_count[0][1]]
    if len(x_count_max) == 1:
        return((x_count_max[0][0],x_count_max[0][1],distance.hamming(x_count_max[0][0], target_re)))
    else:
        x_ham = [(re,a,distance.hamming(re, target_re)) for re,a in x_count_max]
        return(min(x_ham, key=itemgetter(2)))


def collapser(filename, inputdir, outputdir, target_re):

    readsname = os.path.splitext(filename)[0]
    print(readsname)
    df = pd.read_table(inputdir + filename, '\t')
    df = df.drop_duplicates('READNAME', keep=False)
    df['BARCODE'].fillna('', inplace=True)
    df['BARCODE_Q'].fillna('', inplace=True)
    humantable = open(outputdir + readsname + '_humanread.txt', 'w')
    humantable.write('\t'.join(['CLUSTER_ID',
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
                                'NUM_BC',
                                'NUM_TLEN']) + '\n')
    pctable = open(outputdir + readsname + '_pcread.txt', 'w')
    pctable.write('\t'.join(['CLUSTER_ID',
                             'ID_LIST',
                             'FILENAME',
                             'READNAME',
                             'CHR',
                             'POS',
                             'INS_STRAND',
                             'RE',
                             'RE_AMOUNT',
                             'RE_HAMMING',
                             'RE_LIST',
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
                             'NUM_BC',
                             'NUM_TLEN',
                             'BARCODE_LIST',
                             'BARCODE_Q_LIST',
                             'TLEN_LIST']) + '\n')

    df['MDR1_value'] = df['MDFLAG_R1'].apply(lambda x: sum([int(i) for i in re.findall(r'\d+', x)]))
    df_group = df.groupby(['CHR', 'INS_STRAND', 'POS'])

    cluster_id = 0
    for (chrom, strand, pos), group in tqdm(df_group, desc=readsname):
        cluster_id += 1
        best_row = group.loc[group['MDR1_value'] == max(list(group['MDR1_value']))].iloc[0]
        best_re = get_best_re(list(group['RE']), target_re)
        num_bc = len(set(list(group['BARCODE'])))
        num_tlen = len(set(list(group['TLEN'])))
        cigar_del = sum([int(x) for x in re.findall(r'(\d+)+D', best_row['CIGAR_R1'])])
        cigar_ins = sum([int(x) for x in re.findall(r'(\d+)+I', best_row['CIGAR_R1'])])
        md_mm = len(re.findall('[A-Z]', best_row['MDFLAG_R1'].split('MD:Z:')[1])) - cigar_del
        md_sum = best_row['MDR1_value'] + md_mm
        humantable.write('\t'.join([str(cluster_id),
                                    best_row['READNAME'],
                                    chrom,
                                    str(pos),
                                    strand,
                                    best_re[0],
                                    str(best_re[1]),
                                    str(best_re[2]),
                                    best_row['R1'],
                                    best_row['R2'],
                                    str(best_row['TLEN']),
                                    best_row['CIGAR_R1'],
                                    best_row['MDFLAG_R1'],
                                    str(md_sum),
                                    str(md_mm),
                                    str(cigar_ins),
                                    str(cigar_del),
                                    str(np.shape(group)[0]),
                                    str(num_bc),
                                    str(num_tlen)]) + '\n')
        pctable.write('\t'.join([str(cluster_id),
                                 ','.join([str(x) for x in list(group['ID'])]),
                                 best_row['FILENAME'],
                                 best_row['READNAME'],
                                 chrom,
                                 str(pos),
                                 strand,
                                 best_re[0],
                                 str(best_re[1]),
                                 str(best_re[2]),
                                 ','.join(list(group['RE'])),
                                 best_row['R1'],
                                 best_row['R2'],
                                 str(best_row['TLEN']),
                                 best_row['CIGAR_R1'],
                                 best_row['MDFLAG_R1'],
                                 str(md_sum),
                                 str(md_mm),
                                 str(cigar_ins),
                                 str(cigar_del),
                                 str(np.shape(group)[0]),
                                 str(num_bc),
                                 str(num_tlen),
                                 ''.join(list(group['BARCODE'])),
                                 ''.join(list(group['BARCODE_Q'])),
                                 ','.join([str(x) for x in list(group['TLEN'])])]) + '\n')


def main(inputdir, outputdir, target_re, n_core):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]

    if len(onlyfiles) == 1:
        filename = onlyfiles[0]
        stat_series = collapser(filename,
                              inputdir, outputdir, target_re)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(collapser)(filename,
                                            inputdir, outputdir, target_re)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()
