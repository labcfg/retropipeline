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

import time

import pandas as pd
import re
from collections import defaultdict
from collections import Counter
import sys, os, re
from os import listdir
from os.path import isfile, join
import numpy as nps

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
        return((x_count_max[0][0],
                x_count_max[0][1],
                distance.hamming(x_count_max[0][0], target_re)))
    else:
        x_ham = [(re,a,distance.hamming(re, target_re)) for re,a in x_count_max]
        return(min(x_ham, key=itemgetter(2)))


def merge_rows(list_rows,
               target_re,
               metacluster_id,
               chrom,
               strand,
               files_name,
               colnames_filters,
               blen):
    all_mdflag_int = [re.findall(r'\d+', row['MDFLAG_R1']) for row in list_rows]
    best_idx = all_mdflag_int.index(max(all_mdflag_int))
    best_row = list_rows[best_idx]
    pos = Counter([row['POS'] for row in list_rows]).most_common(1)[0][0]
    num_reads_by_files = defaultdict()
    num_barcodes_by_files = defaultdict()
    num_tlen_by_files = defaultdict()
    filecluster_list = defaultdict()
    barcode = defaultdict()
    barcode_q = defaultdict()
    tlen_dict = defaultdict()
    re_list = []
    for i in files_name:
        filecluster_list[i] = '0'
        barcode[i] = blen*'R'
        barcode_q[i] = blen*'R'
        num_reads_by_files[i] = 0
        tlen_dict[i] = []
    for row in list_rows:
        #print(row)
        num_reads_by_files[row['FILENAME']] += row['NUM_READS']
        barcode[row['FILENAME']] += row['BARCODE_LIST']
        barcode_q[row['FILENAME']] += row['BARCODE_Q_LIST']
        re_list.extend(row['RE_LIST'].split(','))
        if filecluster_list[row['FILENAME']] == '0':
            filecluster_list[row['FILENAME']] = str(int(row['CLUSTER_ID']))
        else:
            filecluster_list[row['FILENAME']] += ','+str(int(row['CLUSTER_ID']))
        tlen_dict[row['FILENAME']] += str(row['TLEN_LIST']).split(',')
    best_re = get_best_re(re_list, target_re)
    for i in files_name:
        num_barcodes_by_files[i] = len(set(re.findall(r'.'*blen, barcode[i]))) - 1
        num_tlen_by_files[i] = len(set(tlen_dict[i]))
    ht = [str(metacluster_id),
          best_row['FILENAME'],
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
          str(best_row['MD_SUM']),
          str(best_row['MISMATCH']),
          str(best_row['INSERTION']),
          str(best_row['DELETION']),
          '\t'.join([str(num_reads_by_files[x]) for x in files_name]),
          '\t'.join([str(num_barcodes_by_files[x]) for x in files_name]),
          '\t'.join([str(num_tlen_by_files[x]) for x in files_name])]
    ht[18:18] = [str(best_row[col]) for col in colnames_filters]
    pct = [str(metacluster_id),
           ';'.join([str(key) for key, value in filecluster_list.items()]),
           ';'.join([str(value) for key, value in filecluster_list.items()]),
           best_row['FILENAME'],
           best_row['READNAME'],
           chrom,
           str(pos),
           strand,
           best_re[0],
           str(best_re[1]),
           str(best_re[2]),
           ','.join(re_list),
           best_row['R1'],
           best_row['R2'],
           str(best_row['TLEN']),
           best_row['CIGAR_R1'],
           best_row['MDFLAG_R1'],
           str(best_row['MD_SUM']),
           str(best_row['MISMATCH']),
           str(best_row['INSERTION']),
           str(best_row['DELETION']),
           '\t'.join([str(num_reads_by_files[x]) for x in files_name]),
           '\t'.join([str(num_barcodes_by_files[x]) for x in files_name]),
           '\t'.join([str(num_tlen_by_files[x]) for x in files_name]),
           '\t'.join([str(barcode[x]) for x in files_name]),
           '\t'.join([str(barcode_q[x]) for x in files_name]),
           '\t'.join([','.join(tlen_dict[x]) for x in files_name])]
    pct[20:20] = [str(best_row[col]) for col in colnames_filters]
    ht = '\t'.join(ht)
    pct = '\t'.join(pct)
    table_row = {'ht': ht, 'pct': pct}
    return(table_row)


def metaclustering(df, window,
                   metacluster_id, chrom, strand, target_re,
                   table1, table2, files_name, colnames_filters, blen):
    is_cluster_open = False
    for index, row in df.iterrows():
        if not is_cluster_open:
            is_cluster_open = True
            metacluster_id += 1
            list_row = [row]
            pos = row['POS']
        else:
            if abs(int(row['POS']) - pos) <= window:
                list_row.append(row)
            else:
                is_cluster_open = False
                trow = merge_rows(list_row,
                                  target_re,
                                  metacluster_id,
                                  chrom,
                                  strand,
                                  files_name,
                                  colnames_filters,
                                  blen)
                table1.write(trow['ht'] + '\n')
                table2.write(trow['pct'] + '\n')

                is_cluster_open = True
                metacluster_id += 1
                list_row = [row]
                pos = row['POS']

    if is_cluster_open:
        is_cluster_open = False
        trow = merge_rows(list_row,
                          target_re,
                          metacluster_id,
                          chrom,
                          strand,
                          files_name,
                          colnames_filters,
                          blen)
        table1.write(trow['ht'] + '\n')
        table2.write(trow['pct'] + '\n')
    return(metacluster_id)


def main(inputdir, pcdir, outputdir, window, target_re, blen):

    inputdir = os.path.abspath(inputdir) + '/'
    pcdir = os.path.abspath(pcdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = sorted([f for f in listdir(inputdir) if isfile(join(inputdir, f)) and not re.search('pcread', f)])
    pcfiles = sorted([f for f in listdir(pcdir) if isfile(join(pcdir, f)) and re.search('pcread', f)])
    #pcfiles = [x for x in pcfiles if re.search('pcread', x)]
    #pcfiles = sorted([x for x in pcfiles if x.split('_pcread')[0] + '.txt' in onlyfiles])

    # dummy check
    paired_files = []
    unpaired_files = []
    for target_file in onlyfiles:
        is_pair = False
        for pcfile in pcfiles:
            if is_pair: break
            if target_file.split('__')[0] == pcfile.split('__')[0]:
                is_pair = True
                paired_files.append((target_file, pcfile))
        if not is_pair:
            unpaired_files.append(target_file)

    def bold_text(x):
        return "\033[1m" + x + "\033[0m"

    print(bold_text("Paired:"))
    for x1, x2 in paired_files:
        print(x1 + ' <----- ' + x2)
    print('\n')
    print(bold_text("Unpaired:"))
    if len(unpaired_files) > 0:
        for x in unpaired_files:
            print(x)
    else:
        print("NULL", end='\n'*2)

    files_name = sorted([x.split('___')[0] + '__' if re.search('humanread', x) else os.path.splitext(x)[0] for x in onlyfiles])

    table1 = open(outputdir + 'metatable_humanread.txt', 'w')
    table2 = open(outputdir + 'metatable_pcread.txt', 'w')

    autosomeXY = ['chr' + str(x) for x in range(1, 23)]
    autosomeXY.extend(['chrX', 'chrY'])

    time.sleep(0.5)

    is_meta_table = False
    
    for hfile, pcfile in tqdm(paired_files):
        htable = pd.read_table(inputdir + hfile, '\t')
        htable_cols = list(htable.columns)
        htable = htable.sort_values(by='CLUSTER_ID')
        pctable = pd.read_table(pcdir + pcfile, '\t')
        pctable['BARCODE_LIST'].fillna('', inplace=True)
        pctable['BARCODE_Q_LIST'].fillna('', inplace=True)
        pctable_cols = list(pctable.columns)
        pctable = pctable.sort_values(by='CLUSTER_ID')
        pccluster = list(pctable['CLUSTER_ID'])
        hcluster = list(htable['CLUSTER_ID'])
        if not is_meta_table:
            #pctable_cut = pctable[[x in hcluster for x in pccluster]]
            #htable.index = hcluster
            #pctable_cut.index = hcluster
            #return (htable_cols, pctable_cols)
            #merge_cols = list(set(htable_cols) & set(pctable_cols))
            meta_table = pd.merge(htable, pctable, how='left')
            is_meta_table = True

            colnames_ht = htable_cols[:]
            colnames_ht[0] = 'METACLUSTER_ID'
            colnames_ht.insert(1, 'FILENAME')
            colnames_ht.remove('NUM_READS')
            colnames_ht.remove('NUM_BC')
            colnames_ht.remove('NUM_TLEN')
            colnames_ht.extend([x + '_NUM_READS' for x in files_name])
            colnames_ht.extend([x + '_NUM_BC' for x in files_name])
            colnames_ht.extend([x + '_NUM_TLEN' for x in files_name])
            table1.write('\t'.join(colnames_ht) + '\n')
            colnames_pct = colnames_ht[:]
            colnames_pct[1:1] = ['FILES', 'FILES_ID']
            colnames_pct.insert(11, 'RE_LIST')
            colnames_pct.extend([x + '_BARCODE_LIST' for x in files_name])
            colnames_pct.extend([x + '_BARCODE_Q_LIST' for x in files_name])
            colnames_pct.extend([x + '_TLEN_LIST' for x in files_name])
            table2.write('\t'.join(colnames_pct) + '\n')
            colnames_filters = htable_cols[20:]
            #print(pd.isnull(meta_table['FILENAME']))
            if np.any(pd.isnull(meta_table['FILENAME'])):
                sys.exit('FUCKING_FILENAME    ' + hfile)

        else:
            #sys.exit('Success')
            #pctable_cut = pctable[[x in hcluster for x in pccluster]]
            #htable.index = hcluster
            #pctable_cut.index = hcluster
            #merge_cols = list(set(htable_cols) & set(pctable_cols))
            super_table = pd.merge(htable, pctable, how='left')
            if np.any(pd.isnull(super_table['FILENAME'])):
                sys.exit('FUCKING_FILENAME    ' + hfile)
            meta_table = meta_table.append(super_table)

    mgroup = meta_table.groupby(['CHR', 'INS_STRAND'])
    metacluster_id = 0
    for name, group in tqdm_notebook(mgroup, desc='chrom+strand'):
        group = group.sort_values(['POS'])
        metacluster_id = metaclustering(group, window, metacluster_id,
                                        name[0], name[1], target_re,
                                        table1, table2, files_name, colnames_filters, blen)
    table1.close()
    table2.close()
