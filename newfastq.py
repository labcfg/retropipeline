from Bio import SeqIO
import sys, os, re
import numpy as np
from os import listdir
from os.path import isfile, join
import pandas as pd

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
from IPython.display import display
from operator import itemgetter
from datetime import datetime


'''
Count lines if file
'''
def count_lines(filepath):
    with open(filepath) as f:
        return (sum(1 for _ in f))

def recovery(readname_list, inputdirfq, filename1, filename2, name, outputdir):
    R1_reads = SeqIO.parse(inputdirfq + filename1, "fastq")
    R2_reads = SeqIO.parse(inputdirfq + filename2, "fastq")
    goodr1 = open(outputdir + name + '_R1_good.fastq', 'w')
    goodr2 = open(outputdir + name + '_R2_good.fastq', 'w')
    R12_reads = zip(R1_reads, R2_reads)
    readname_list_set = set(readname_list)
    bar = tnrange(int(count_lines(inputdirfq+filename1)/4), desc=name)
    for i in bar:
        r1,r2 = next(R12_reads)
        if r1.id.split('__abq')[0] in readname_list_set:
            goodr1.write(r1.format('fastq'))
            goodr2.write(r2.format('fastq'))
            #readname_list.remove(r1.id.split('__abq')[0])
    goodr1.close()
    goodr2.close()
    return (0)


def main(inputdir, outputdir, inputdirfq):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    # Read files in folder
    onlyfiles = [f for f in listdir(inputdirfq) if isfile(join(inputdirfq, f))]

    r1_files = {}
    r2_files = {}

    for filename in onlyfiles:
        filename = filename.rstrip()
        if re.search('good', filename):
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
            conform_files.append([r1_files[key], r2_files[key]])
            del r2_files[key]
        else: nonconform_files.append(r1_files[key])

    nonconform_files = nonconform_files + list(r2_files.values())

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    if len(nonconform_files) != 0:
        print ('I can\'t read this files' + str(nonconform_files))

    filedict = {}
    for filename1, filename2 in conform_files:
        readsname = filename1.split('R1')[0]
        readsname = readsname.rsplit('.', 1)[0]
        filedict[readsname] = [filename1, filename2]

    inputfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]
    for f in inputfiles:
        name, ext = os.path.splitext(f)
        if name in list(filedict.keys()):
            data = pd.read_table(inputdir + f, '\t')
            data = list(data['READNAME'])
            file1 = list(filedict[name])[0]
            file2 = list(filedict[name])[1]
            recovery(data, inputdirfq, file1, file2, name, outputdir)
