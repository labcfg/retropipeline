from simplesam import Reader, Writer
import sys, os, re
from os import listdir
from os.path import isfile, join
from tqdm import tqdm, tnrange
from datetime import datetime
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from collections import defaultdict
import gline
import time
import pysam

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

def get_needed_pars(r1, r2, sam):
    if r1.is_reverse:
        ins_strand = '+'
        pos = r2.pos + abs(r1.tlen) - 1
    else:
        ins_strand = '-'
        pos = r1.pos
    r1_mdflag = [x for x in r1.get_tags() if x[0] == 'MD']
    r2_mdflag = [x for x in r2.get_tags() if x[0] == 'MD']
    try:
        r1_mdflag = 'MD:Z:' + r1_mdflag[0][1]
        r2_mdflag = 'MD:Z:' + r2_mdflag[0][1]
    except:
        print(r1_mdflag, r2_mdflag)
    result = [sam.getrname(r1.rname), str(pos), ins_strand, str(abs(r1.tlen)),
              r1.cigarstring, r2.cigarstring, r1_mdflag, r2_mdflag]
    return(result)


def sam2table(filename, flags, inputdir, outputdir):

    autosomeXY = list(range(1, 23))
    autosomeXY.append('X')
    autosomeXY.append('Y')
    autosomeXY = ['chr' + str(x) for x in autosomeXY]

    readsname = os.path.splitext(filename)[0]
    
    preprocessing_dir = os.path.dirname(os.path.dirname(inputdir)) + '/preprocessing/'
    init_file_pathway = preprocessing_dir + readsname + 'meta.txt'
    temp_file = open(init_file_pathway, 'r')
    temp_file.close()
    #sys.exit(init_file_pathway)

    samfile = open(inputdir + filename, 'r')
    tablefile = open(outputdir + readsname + '.txt', 'w')
    #errorfile = open(outputdir + readsname + '_error.sam', 'w')
    #error = Writer(errorfile)
    sam = pysam.AlignmentFile(samfile)
    colnames = ['ID','FILENAME','READNAME',
                'CHR','POS','INS_STRAND','RE','R1','R2',
                'TLEN','CIGAR_R1','CIGAR_R2','MDFLAG_R1','MDFLAG_R2',
                'BARCODE','BARCODE_Q']
    tablefile.write('\t'.join(colnames) + '\n')
    print('start: ' + readsname)
    bar = tnrange(int(count_samlines(inputdir+filename)), desc=readsname)
    r1_bool = False
    r2_bool = False
    id_count = 1
    for i in bar:
        r = next(sam)
        try:
            if sam.getrname(r.rname) in autosomeXY and r.flag in flags:
                if r.qname.split('__ct:')[-1] == 'r1':
                    r1 = r
                    r1_bool = True
                    if r2_bool:
                        r1_name = r1.qname.split('__ct:')[0][:-1]
                        r2_name = r2.qname.split('__ct:')[0][:-1]
                        if r1_name == r2_name and abs(r1.tlen) == abs(r2.tlen):
                            c_nmbr = int(r1.qname.split('__ct:')[1])
                            real_qname = gline.getline(init_file_pathway, c_nmbr)
                            real_qname = real_qname.strip()
                            real_name_list = real_qname.split('__a12bq:')
                            line = get_needed_pars(r1, r2, sam)
                            super_line = [real_name_list[0], line[0], line[1], line[2],
                                        real_name_list[1], real_name_list[2], real_name_list[3], line[3],
                                        line[4], line[5], line[6], line[7],
                                        real_name_list[4], real_name_list[5]]
                            super_line = '\t'.join(super_line)
                            tablefile.write('\t'.join([str(id_count),
                                                    readsname,
                                                    super_line]) + '\n')
                            #sys.exit(line)
                            r1_bool, r2_bool = False, False
                            id_count += 1
                else:
                    r2 = r
                    r2_bool = True
                    if r1_bool:
                        r1_name = r1.qname.split('__ct:')[0][:-1]
                        r2_name = r2.qname.split('__ct:')[0][:-1]
                        if r1_name == r2_name and abs(r1.tlen) == abs(r2.tlen):
                            c_nmbr = int(r1.qname.split('__ct:')[1])
                            real_qname = gline.getline(init_file_pathway, c_nmbr)
                            real_qname = real_qname.strip()
                            real_name_list = real_qname.split('__a12bq:')
                            line = get_needed_pars(r1, r2, sam)
                            super_line = [real_name_list[0], line[0], line[1], line[2],
                                        real_name_list[1], real_name_list[2], real_name_list[3], line[3],
                                        line[4], line[5], line[6], line[7],
                                        real_name_list[4], real_name_list[5]]
                            super_line = '\t'.join(super_line)
                            tablefile.write('\t'.join([str(id_count),
                                                        readsname,
                                                        super_line]) + '\n')
                            #sys.exit(line)
                            r1_bool, r2_bool = False, False
                            id_count += 1
            else:
                pass
                #error.write(r)
        except:
            pass

    samfile.close()
    tablefile.close()
    #errorfile.close()
    return(0)

def main(inputdir, outputdir, flags, n_core):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.sam')]

    if len(onlyfiles) == 1:
        filename = onlyfiles[0]
        stat_series = sam2table(filename, flags,
                              inputdir, outputdir)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(sam2table)(filename, flags,
                                            inputdir, outputdir)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()
