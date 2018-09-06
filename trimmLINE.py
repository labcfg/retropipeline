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

info = namedtuple('info', 'is_good read alu_barcode errors')

"""
Trimming for LINE. Difference between ALU and LINE:
- Multiple primers and multiples RE fragments
- Trimm polyA-rich fragments 
"""


'''
Some functions for trimm polyN
'''
# count rolling statistics (of N freq)
def count_in_rolling_window(seq, n, win_size):
    count = np.ones(len(seq))
    for i in range(len(seq)-win_size+1):
        sub_seq = seq[i:i+win_size]
        count[i+win_size-1] = sub_seq.count(n) / win_size
    return count
# recuresive search for end of polyN
def detect_poly_n_mod(seq, n, win_size, th, shift, rev=0):
    count_n = count_in_rolling_window(seq, n, win_size)
    cut_point = np.argwhere(count_n < th).reshape(1, -1)[0][rev]
    count_shift = 0
    if shift > 0 and len(seq) - (cut_point + shift + 1) > 0:
        for i in range(shift):
            if count_n[cut_point+i+1] > th:
                count_shift += 1
    if count_shift > 0:
        return detect_poly_n_mod(seq, n, win_size, th, shift, rev=rev+1)
    
    return cut_point
# usless function
# TODO: improve. set double win_size for additional correction
def poly_n_point(seq, n, win_size, th, shift, kind='head'):
    if kind == 'head':
        point = detect_poly_n_mod(seq, n, win_size, th, shift)
    elif kind == 'tail':
        point = detect_poly_n_mod(seq[::-1], n, win_size, th, shift)
    else:
        sys.exit("Variable kind must be 'head' or 'tail'")
    return point



'''
check: is r1 is real r1?
(find primer in start of seq with shift)
'''
def is_r1(record, primer_dict, shift, mist):
    for i in range(shift):
        for primer, _ in primer_dict.items():
            if hamming(primer, str(record.seq[i : len(primer)+i]), mist):
                return True
    return False


def trimm_short_reads(target, record, mid_mist, end_mist, initial_pos=None, save_pos=0, cut_pos=0):
    # initial pos - from tail 
    seq = str(record.seq)
    if len(seq) < len(target)+1:
        return record
    if initial_pos is None:
        initial_pos = len(seq)
    initial_pos = len(seq) - initial_pos
    ham_i = []
    for i in range(initial_pos, len(seq) - len(target) + 1):
        ham = distance.hamming(target, seq[i:i+len(target)])
        ham_i.append((ham, i))
        min_ham = sorted(ham_i, key=lambda x: (x[0], x[1]))[0]
    if min_ham[0] <= mid_mist:
        p = min_ham[1] + save_pos - cut_pos
        if p < 1:
            return None
        return record[:p]
    k = 1
    for i in range(len(seq) - len(target) + 1, len(seq) - len(target) + 1 + end_mist):
        if seq[i : (i + len(target[:-k]))] == target[:-k]:
            p = i + save_pos - cut_pos
            if p < 1:
                return None
            return record[:p]
        k += 1
    return record

'''
ONLY R1 (read1)
trimming primers

1. find primer with hamming in start of seq with shift
2. find primer and adapters in natural_zone* and remove if in
3. find restrict_site (seq like 'AGCT' or 'CTAG' - restriction sites) and remove if in
4. trimm primer and return seq

errors = [primer, adapters, restrict_site, natural_zone, r2_start**]

* natural_zone - zone without primer or adapters
** r2_start - nucleotids in start of trim seq (for example: 'CT') (only r2)
'''
def trim_primers (record, primer_dict, ad1, shift, mist, restrict_site, target, r2_start, 
                  mid_mist_short_reads, end_mist_short_reads, place_of_search_tail,
                  poly_n_win_size_r1, poly_n_th_r1, poly_n_shift_r1, re_mist):

    for i in range(shift):
        for primer, re_part in zip(primer_dict['primer'], primer_dict['re_part']):
            sys.exit(primer)
            len_primer = len(primer)
            if hamming(primer, str(record.seq[i : len_primer + i]), mist):
                #sys.exit('First step, Done!')
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
                re_i = 0
                for i_nucl in range(1, len(re_part)):
                    seq_part = str(record.seq[len_primer+i : len(primer)+i_nucl+i])
                    if hamming(seq_part, re_part[:i_nucl], re_mist):
                        re_i = i_nucl
                        
                re_part = str(record.seq[len_primer+i : len_primer+re_i+i])
                record = record[len_primer+len(re_part)+i :]
                record = trimm_short_reads(str(Seq(target).reverse_complement()), record,
                                        mid_mist=mid_mist_short_reads[0], 
                                        end_mist=end_mist_short_reads[0],
                                        initial_pos=place_of_search_tail[0],
                                        save_pos=len(r2_start))
                try:
                    #point = poly_n_point(str(record.seq), 'T', poly_n_win_size_r1, poly_n_th_r1, poly_n_shift_r1, kind='head')
                    #record = record[point:]
                    1 + 1
                except:
                    return (info(is_good=False, read=None,
                                alu_barcode=None,
                                errors=np.array([0, 0, 0, 0, 0])))
                    #print('wow')
                if record is None:
                    return (info(is_good=False, read=None,
                                alu_barcode=None,
                                errors=np.array([0, 0, 0, 0, 0])))
                alu_bar = primer + '__pr12bq:' + re_part + '__pr12bq:' + str(record.seq)
                return (info(is_good=True, read=record,
                            alu_barcode=alu_bar,
                            errors=np.array([0, 0, 0, 0, 0])))
    return (info(is_good=False, read=None,
                 alu_barcode=None,
                 errors=np.array([1, 0, 0, 0, 0])))


'''
ONLY R2 (read2)
trimming adapters

1. find adapters with hamming in start of seq with shift
    pos(ad1) - pos(ad2) = barlen (length of barcode - UMI)
2. find primer and adapters in natural_zone* and remove if in
3. find restrict_site (seq like 'AGCT' or 'CTAG' - restriction sites) and remove if in
4. find r2_start** and remove if not (r2_start can be empty: '')
5. trimm adapters (with UMI) and return seq (with UMI and its quality)

errors = [primer, adapters, restrict_site, natural_zone, r2_start]

* natural_zone - zone without primer or adapters
** r2_start - nucleotids in start of trim seq (for example: 'CT') (only r2)
'''
def trim_adapters (record, ad1, ad2, blen, shift, mist, restrict_site, r2_start, target, re_part,
                   mid_mist_short_reads, end_mist_short_reads, place_of_search_tail,
                   poly_n_win_size_r2, poly_n_th_r2, poly_n_shift_r2):

    len_ad1 = len(ad1)
    len_ad2 = len(ad2)
    trim_len = len_ad1 + blen + len_ad2
    for i in range(shift):
        seq1 = record.seq[i : len_ad1+i]
        seq2 = record.seq[len_ad1+blen+i : len_ad1+blen+len_ad2+i]
        if hamming(ad1, str(seq1), mist) and hamming(ad2, str(seq2), mist):
            if record.seq[trim_len+i : trim_len+i+len(r2_start)] != r2_start:
                return(info(is_good=False, read=None,
                            alu_barcode=None,
                            errors=np.array([0, 0, 0, 0, 1])))
            if record.seq[len_ad1+blen+len_ad2+i :].find(restrict_site, 0) != -1:
                return (info(is_good=False, read=None,
                             alu_barcode=None,
                             errors=np.array([0, 0, 1, 0, 0])))
            record.description = ''
            record.name = ''
            barcode = str(record.seq[len_ad1+i : len_ad1+blen+i])
            barcode_q = [chr(x + 33)
                         for x in record.letter_annotations['phred_quality']]
            barcode_q = ''.join(barcode_q)
            barcode_q = barcode_q[len_ad1+i : len_ad1+blen+i]
            record = record[trim_len+i :]
            record = trimm_short_reads(str(Seq(target).reverse_complement()), record,
                                       mid_mist=mid_mist_short_reads[1],
                                       end_mist=end_mist_short_reads[1],
                                       initial_pos=place_of_search_tail[1],
                                       cut_pos=len(re_part))
            try:
                #point = poly_n_point(str(record.seq), 'A', poly_n_win_size_r2, poly_n_th_r2, poly_n_shift_r2, kind='tail')
                #record = record[:(-point+1)]
                1 + 1
            except:
                #print('wow')
                return (info(is_good=False, read=None,
                             alu_barcode=None,
                             errors=np.array([0, 0, 0, 0, 0])))
            if record is None:
                return (info(is_good=False, read=None,
                             alu_barcode=None,
                             errors=np.array([0, 0, 0, 0, 0])))
            alu_bar = '__pr12bq:' + str(record.seq) + '__pr12bq:' + str(barcode) + '__pr12bq:' + str(barcode_q)
            return (info(is_good=True, read=record,
                         alu_barcode=alu_bar,
                         errors=np.array([0, 0, 0, 0, 0])))
    return(info(is_good=False, read=None,
                alu_barcode=None,
                errors=np.array([0, 1, 0, 0, 0])))


'''
Humanreadble for error statistics
'''
def concate_errors(type_err, amount_err):
    result = ','.join([str(t) + '-' + str(a)
                      for t, a in zip(type_err, amount_err)])
    return (result)


'''
Count lines if file
'''
def count_lines(filepath):
    with gzip.open(filepath) as f:
        return (sum(1 for _ in f))

'''
trim reads
1. SeqIO.parse - generator (read 'read' by 'read')
2. Chaos = T | F - if r1 and r2 mix in files (R1, R2)
3. Trim and if good write in goodfile else badfile
4. Print and return statistics
'''
def trim_reads(filename1, filename2, inputdir, outputdir,
               primer_dict, ad1, ad2, blen, shift, mist,
               restrict_site, r2_start, is_short_flank, chaos, trimm_n, skip_short_reads,
               mid_mist_short_reads, end_mist_short_reads, place_of_search_tail,
               min_seq_len_after_trimm, poly_n_win_size, poly_n_th, poly_n_shift, re_mist):

    fname = filename1.split('R1')[0].rsplit('.', 1)[0]
    fnameout = filename1.split('R1')[0]
    outputfile1, ext = os.path.splitext(filename1)
    outputfile2, ext = os.path.splitext(filename2)
    goodr1 = open(outputdir + fnameout + '_R1_good.fastq', 'w')
    goodr2 = open(outputdir + fnameout + '_R2_good.fastq', 'w')
    badr1 = open(outputdir + fnameout + '_R1_bad.fastq', 'w')
    badr2 = open(outputdir + fnameout + '_R2_bad.fastq', 'w')
    
    metafile = open(outputdir + fnameout + '_meta.txt', 'w')
    
    print('start: ' + fname)
    unzip_r1 = gzip.open(inputdir + filename1, 'rt')
    unzip_r2 = gzip.open(inputdir + filename2, 'rt')
    original_R1_reads = SeqIO.parse(unzip_r1, "fastq")
    original_R2_reads = SeqIO.parse(unzip_r2, "fastq")

    count = np.array([0, 0, 0, 0, 0])
    elem = ('primer', 'adapters', 'restrict_site', 'natural_zone', 'r2_start')
    count_stat = {'fname':fname, 'all':0, 'good':0, 'bad':0,
                   'primer':0, 'adapters':0,
                   'restrict_site':0, 'natural_zone':0, 'r2_start':0}
    count_stat_col = ['fname', 'all', 'good', 'bad',
                      'primer', 'adapters',
                      'restrict_site', 'natural_zone', 'r2_start']
    bar = tnrange(int(count_lines(inputdir+filename1)/4), desc=fname)
    original_R12 = zip(original_R1_reads, original_R2_reads)
    for i in bar:
        count_stat['all'] += 1
        r1, r2 = next(original_R12)
        if chaos:
            if is_r1(r1, primer_dict, shift, mist):
                rx1 = r1
                rx2 = r2
            else:
                rx1 = r2
                rx2 = r1
        else:
            rx1 = r1
            rx2 = r2
        if len(str(rx1.seq)) > skip_short_reads and len(str(rx2.seq)) > skip_short_reads:
            if trimm_n:
                rx1 = rx1[:-trimm_n]
                rx2 = rx2[:-trimm_n]
            if is_short_flank:
                fr1 = trim_primers(rx1, primer_dict, '$',
                                   shift, mist, restrict_site, ad2+r2_start, r2_start,
                                   mid_mist_short_reads, end_mist_short_reads,
                                   place_of_search_tail, poly_n_win_size[0], poly_n_th[0], poly_n_shift[0], re_mist)
            else:
                fr1 = trim_primers(rx1, primer_dict, ad1,
                                   shift, mist, restrict_site, ad2+r2_start, r2_start,
                                   mid_mist_short_reads, end_mist_short_reads,
                                   place_of_search_tail, poly_n_win_size[0], poly_n_th[0], poly_n_shift[0], re_mist)
            if fr1.is_good:
                target_primer = fr1.alu_barcode.split('__')[0]
                target_re_part = fr1.alu_barcode.split('__')[1]
                fr2 = trim_adapters(rx2, ad1, ad2, blen,
                                    shift, mist, restrict_site, r2_start, target_primer, target_re_part,
                                    mid_mist_short_reads, end_mist_short_reads,
                                    place_of_search_tail, poly_n_win_size[1], poly_n_th[1], poly_n_shift[1])
                if fr2.is_good and len(str(fr1.read.seq)) > min_seq_len_after_trimm[0] and len(str(fr2.read.seq)) > min_seq_len_after_trimm[1]:
                    count_stat['good'] += 1
                    meta_header = fr1.read.id + fr1.alu_barcode + fr2.alu_barcode
                    fr1.read.id += '__ct:' + str(count_stat['good']) + '__ct:' + 'r1'
                    fr2.read.id += '__ct:' + str(count_stat['good']) + '__ct:' + 'r2'
                    goodr1.write(fr1.read.format('fastq'))
                    goodr2.write(fr2.read.format('fastq'))
                    metafile.write(str(meta_header) + '\n')
                else:
                    rx2.description += (' reason:' +
                            concate_errors(elem, (np.char.mod('%d', fr2.errors))))
                    badr1.write(rx1.format('fastq'))
                    badr2.write(rx2.format('fastq'))
                    count = np.sum([count, fr2.errors], axis=0)
            else:
                rx2.description += (' reason:' +
                        concate_errors(elem, (np.char.mod('%d', fr1.errors))))
                badr1.write(rx1.format('fastq'))
                badr2.write(rx2.format('fastq'))
                count = np.sum([count, fr1.errors], axis=0)
        else:
            rx2.description += (' reason: skip_short_reads')
            badr1.write(rx1.format('fastq'))
            badr2.write(rx2.format('fastq')) 
    
    metafile.close()
    unzip_r1.close()
    unzip_r2.close()
    goodr1.close()
    goodr2.close()
    badr1.close()
    badr2.close()

    count_stat['primer'] = count[0]
    count_stat['adapters'] = count[1]
    count_stat['restrict_site'] = count[2]
    count_stat['natural_zone'] = count[3]
    count_stat['r2_start'] = count[4]

    count_stat['bad'] = round((count_stat['all'] - count_stat['good']) / count_stat['all'], 2)
    count_stat['good'] = round(1 - count_stat['bad'], 2)
    count_stat_pd = pd.Series(count_stat, index = count_stat_col)
    count_elem = concate_errors(elem, (np.char.mod('%d', count)))

    return(count_stat_pd)


'''
main
1. find r1,r2-pairs in input directory (*.fastq.gz)
2. run parallel trimming
3. display stat
'''
def main(inputdir, outputdir,
         primer_dict, ad1, ad2, blen,
         shift, mist,
         restrict_site, r2_start, is_short_flank,
         chaos, n_core, trimm_n, skip_short_reads,
         mid_mist_short_reads, end_mist_short_reads, place_of_search_tail,
         min_seq_len_after_trimm,
         poly_n_win_size, poly_n_th, poly_n_shift, re_mist):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'
    
    if isinstance(primer_dict, str):
        primer_dict = pd.read_csv(primer_dict, '\t')
        primer_dict.columns = ['name', 'primer', 're_part']
        primer_dict = primer_dict.to_dict(orient='list')

    # Read files in folder
    onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]

    r1_files = {}
    r2_files = {}

    for filename in onlyfiles:
        filename = filename.rstrip()
        if re.search('R1', filename):
            key_filename = filename.split('R1')[0]
            r1_files[key_filename] = filename
        elif re.search('R2', filename):
            key_filename = filename.split('R2')[0]
            r2_files[key_filename] = filename

    conform_files = []
    nonconform_files = []

    for key in r1_files:
        if key in r2_files:
            conform_files.append((r1_files[key], r2_files[key]))
            del r2_files[key]
        else: nonconform_files.append(r1_files[key])

    conform_files = sorted(conform_files, key=itemgetter(0))
    nonconform_files = nonconform_files + list(r2_files.values())

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    stat_col = ['fname', 'reads', 'good', 'bad',
                'primer', 'adapters',
                'restrict_site', 'natural_zone', 'r2_start']

    if len(conform_files) == 1:
        stat_series = trim_reads(conform_files[0][0], conform_files[0][1],
                              inputdir, outputdir,
                              primer_dict, ad1, ad2, blen,
                              shift, mist,
                              restrict_site, r2_start, is_short_flank, chaos, trimm_n, skip_short_reads,
                              mid_mist_short_reads, end_mist_short_reads, place_of_search_tail,
                              min_seq_len_after_trimm, poly_n_win_size, poly_n_th, poly_n_shift, re_mist)
        stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(trim_reads)(filename1, filename2,
                                            inputdir, outputdir,
                                            primer_dict, ad1, ad2, blen,
                                            shift, mist,
                                            restrict_site, r2_start, is_short_flank, chaos, trimm_n, skip_short_reads, mid_mist_short_reads, end_mist_short_reads, place_of_search_tail,
                              min_seq_len_after_trimm, poly_n_win_size, poly_n_th, poly_n_shift, re_mist)
                                                for filename1, filename2 in conform_files)
        stat_df = pd.concat(stat_series, axis=1).transpose()

    stat_df_total = stat_df.sum(axis = 0)
    stat_df_total.fname = 'total'
    total_bad = sum([stat_df_total.primer,
                     stat_df_total.adapters,
                     stat_df_total.restrict_site,
                     stat_df_total.natural_zone,
                     stat_df_total.r2_start])
    stat_df_total.bad = round(total_bad / stat_df_total['all'], 2)
    stat_df_total.good = round(1 - stat_df_total.bad, 2)
    stat_df.loc[np.shape(stat_df)[0]] = stat_df_total

    display(stat_df)

    stat_name = ''.join([str(before.year),
                         str(before.month),
                         str(before.day),
                         str(before.hour),
                         str(before.minute),
                         str(before.second)])
    stat_df.to_csv(outputdir + 'statistics_' + stat_name + '.csv', sep=' ', index=False)

    after = datetime.now()
    delta_time = after - before

    # Write logfile
    logfile = open(outputdir + 'logfile_' + stat_name + '.log', 'w')
    logfile.write('#DESC OF PARS:\n' +
                  '#MIST - max hamming (for primer and ads)\n' +
                  '#SHIFT - for search primer or ads\n' +
                  '#PRIMER - seq of primer\n' +
                  '#AD1 - seq of adapter1\n' +
                  '#AD2 - seq of adapter2\n' +
                  '#BLEN - len of UMI(barcode)\n' +
                  '#RESTRICT_SITE - \'AGCT\' or \'CTAG\'\n' +
                  '#R2_START - seq after adapters\n' +
                  '#CHAOS - mixed files\n' +
                  '#N_CORE - number of active cores\n')
    logfile.write('time_start = ' + str(before) + '\n')
    logfile.write('time_end = ' + str(after) + '\n')
    logfile.write('duration (in sec) = ' + str(round(delta_time.total_seconds(), 2)) + '\n')
    logfile.write('MIST = ' + str(mist) + '\n' +
                  'SHIFT = ' + str(shift) + '\n' +
                  'PRIMER = ' + 'MULTIPLE' + '\n' +
                  'AD1 = ' + ad1 + '\n' +
                  'AD2 = ' + ad2 + '\n' +
                  'BLEN = ' + str(blen) + '\n' +
                  'RESTRICT_SITE = ' + restrict_site + '\n' +
                  'R2_START = ' + r2_start + '\n' +
                  'CHAOS = ' + str(chaos) + '\n'
                  'N_CORE = ' + str(n_core) + '\n' +
                  'INPUTDIR = ' + inputdir + '\n' +
                  'OUTPUTDIR = ' + outputdir + '\n')
    logfile.close()
    stat_df.to_csv(outputdir + 'logfile_' + stat_name + '.log', index=False, sep=' ', mode='a')

    if len(nonconform_files) != 0:
        print ('I can\'t read this files' + str(nonconform_files))
