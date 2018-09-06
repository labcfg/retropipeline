import sys, argparse, os, re
import pandas as pd
from os import listdir
from os.path import join

parser = argparse.ArgumentParser(description='Create testALU fastq')
parser.add_argument('-i', '--input', help='folder or file with mixcr results')
parser.add_argument('-o', '--output', help='file for processed data')
parser.add_argument('-s', '--sequence', type=str, help='base sequence')
parser.add_argument('-d', '--distance', type=int, help='max edit distance')
parser.add_argument('-ind', '--indels', const=True, default=False,
                    help='is indels allowed (default: prohibited)', nargs='?')

if len(sys.argv) == 0:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

input_smth = os.path.abspath(args.input)
outputfile = os.path.abspath(args.output)
base_seq = args.sequence
max_dist = args.distance
flag = args.indels


    for i in tqdm_notebook(range(df.shape[0]), desc='create fastq'):
        row = df.iloc[i, ]
        seq = Seq(primer + str(row['RE']) + str(row('READ1'))[0:main_flank_len])
        rec = SeqRecord(seq, id=''.join([row['READNAME'],
                                        '__',
                                        row['METACLUSTER_ID'])
        rec.letter_annotations["phred_quality"] = [40] * len(seq)
        recs = [ rec ]
        handle = file(outputdir+readsname+'.fastq', 'at')
        SeqIO.write(recs, handle, 'fastq')
