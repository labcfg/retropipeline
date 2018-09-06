from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys, os, re
from os import listdir
from os.path import isfile, join
import pandas as pd
import subprocess
from simplesam import Reader, Writer
import pandas as pd
import numpy as np
from tqdm import tqdm_notebook, tnrange
from datetime import datetime

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

'''
Count lines if samfile
'''
def count_samlines(filepath):
    with open(filepath) as f:
        i = 0
        for _ in f:
            if _[0] != '@':
                i += 1
        return(i)


def tempate_switch(inputdir, outputdir,
                   primer, main_flank_len, template_switch_md, refway, bwaway):

    readsname = 'metatable'
    filename = 'metatable_humanread.txt'

    df = pd.read_table(inputdir+filename, '\t')
    handle = open(outputdir+readsname+'.fastq', 'wt')
    for i in tqdm_notebook(range(df.shape[0]), desc='create fastq'):
        row = df.iloc[i, ]
        seq = Seq(primer + str(row['RE']) + str(row['R1'])[0:main_flank_len])
        rec = SeqRecord(seq, id=''.join([row['READNAME'],
                                        '__',
                                        str(row['METACLUSTER_ID'])]))
        rec.description = ''
        rec.letter_annotations["phred_quality"] = [40] * len(seq)
        recs = [ rec ]
        SeqIO.write(recs, handle, 'fastq')
    handle.close()

    bwaaln = ' '.join([bwaway,
                       'aln',
                       refway,
                       outputdir+readsname+'.fastq',
                       '>',
                       outputdir+readsname+'.sai'])
    print(bwaaln)
    p_aln = subprocess.Popen(bwaaln, stderr=subprocess.PIPE, shell = True)
    logline = p_aln.stderr.read().decode()
    bwasamse = ' '.join([bwaway,
                        'samse',
                        refway,
                        outputdir+readsname+'.sai',
                        outputdir+readsname+'.fastq',
                        '>',
                        outputdir+readsname+'.sam'])
    print(bwasamse)
    p_samse = subprocess.Popen(bwasamse, stderr=subprocess.PIPE, shell = True)
    logline = p_samse.stderr.read().decode()

    samfile = open(outputdir+readsname+'.sam', 'r')
    sam = Reader(samfile)
    tablefile = open(outputdir+readsname+'_tmp.txt', 'w')
    tablefile.write('\t'.join(['METACLUSTER_ID',
                               'MDFLAG',
                               'MATCH',
                               'TAGS',
                               'TEMPLATE_SWITCH_STRAND']) + '\n')
    bar = tnrange(int(count_samlines(outputdir+readsname+'.sam')),
                  desc='read samfile')
    for i in bar:
        r = next(sam)
        metacluster_id = r.qname.split('__')[1]
        if r.mapped:
            for x in r._tags:
                if re.search('MD', x):
                    mdflag = x
                else:
                    mdflag = '*'
            match = re.findall(r'\d+', mdflag)
            match = sum([int(x) for x in match])
            tags = ','.join([str(x) for x in r._tags])
            if r.reverse:
                strand = '-'
            else:
                strand = '+'
        else:
            mdflag = '*'
            tags = '*'
            strand = '*'
            match = 0
        tablefile.write('\t'.join([str(metacluster_id),
                                   mdflag,
                                   str(match),
                                   tags,
                                   strand]) + '\n')
    samfile.close()
    tablefile.close()

    df_tmp = pd.read_table(outputdir+readsname+'_tmp.txt', '\t')
    mdcolumn = list(np.repeat('*', df.shape[0]))
    mdmatch_list = list(np.zeros(df.shape[0]))
    tags = list(np.repeat('*', df.shape[0]))
    strand = list(np.repeat('*', df.shape[0]))
    for i in tqdm_notebook(range(df.shape[0]), desc='ts'):
        row = df.iloc[i, ]
        meta_id = row['METACLUSTER_ID']
        minitable = df_tmp[df_tmp['METACLUSTER_ID'] == meta_id]
        if len(minitable) == 1:
            mdcolumn[i] = list(minitable['MDFLAG'])[0]
            mdmatch_list[i] = list(minitable['MATCH'])[0]
            tags[i] = list(minitable['TAGS'])[0]
            strand[i] = list(minitable['TEMPLATE_SWITCH_STRAND'])[0]
        elif len(minitable) > 1:
            mdmatch = list(minitable['MATCH'])
            mdcolumn[i] = list(minitable['MDFLAG'])[mdmatch.index(max(mdmatch))]
            mdmatch_list[i] = max(mdmatch)
            tags[i] = list(minitable['TAGS'])[mdmatch.index(max(mdmatch))]
            strand[i] = list(minitable['TEMPLATE_SWITCH_STRAND'])[mdmatch.index(max(mdmatch))]
        else:
            continue

    df['MDFLAG_TEMPLATE'] = mdcolumn
    df['MDMATCH'] = mdmatch_list
    #df = df[df['MDMATCH'] < template_switch_md]
    #df['TAGS'] = tags
    #df['TEMPLATE_SWITCH_STRAND'] = strand
    df.to_csv(outputdir+readsname+'_ts.txt', index=None, sep='\t')


def main(inputdir, outputdir,
         primer, main_flank_len, template_switch_md, refway, bwaway):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    tempate_switch(inputdir, outputdir,
                   primer, main_flank_len, template_switch_md, refway, bwaway)
