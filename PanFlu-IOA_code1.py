#!/bin/python3

import re
import os
import glob
import random
import argparse
import pandas as pd
from Bio.Align.Applications import MafftCommandline
from Bio import SeqIO
from Bio import AlignIO
from collections import defaultdict

from multiprocessing import Process
import multiprocessing

def findMaxFreqAA(rowdata):
	aa = rowdata.value_counts().index[0]
	return aa

def countMaxFreqAA(rowdata):
	aa = rowdata.value_counts().index[0]
	aa_count = rowdata.value_counts()[aa]
	return aa_count

def getConsenseSeq(seqs, seqlen, outdir, step, t):
	#getConsenseSeq(seqs, seq_len, args.outdir, step)
	data = pd.DataFrame.from_dict(seqs)
	data_cons = pd.DataFrame()
	data_cons['Cons_seq'] = data.apply(findMaxFreqAA, axis=1)
	data_cons['Cons_count'] = data.apply(countMaxFreqAA, axis=1)
	data_cons['Cons_ratio'] = data_cons["Cons_count"]*100/seqlen
	data_cons = data_cons.drop(data_cons[data_cons["Cons_seq"]=='-'].index)
	#if os.path.exists(f'{outdir}/{t}.step{step}.stat.xls'):
	#	os.remove(f'{outdir}/{t}.step{step}.stat.xls')
	#	os.remove(f'{outdir}/{t}.step{step}.cons.fasta')
	data_cons.to_csv(f'{outdir}/step{step}.{t}.stat.xls', sep='\t', index=False)
	with open(f'{outdir}/step{step}.{t}.cons.fasta', 'w') as outdata:
		seq = ''.join(data_cons['Cons_seq'].to_list())
		outdata.write(f'>step{step}.{t}.cons\n{seq}\n')

def seqAlign(infile, outfile, outdir, step, t, semaphore):
	semaphore.acquire()
	mafft_cline = MafftCommandline(input=infile)
	aln_res, aln_error = mafft_cline()
	with open(outfile, "w") as handle:
		handle.write(aln_res)
	align = AlignIO.read(outfile, "fasta")
	seqs = defaultdict(list)
	seq_len = 0
	for recode in align:
		seqs[recode.id] = list(recode.seq)
	seq_len = len(seqs)
	getConsenseSeq(seqs, seq_len, outdir, step, t)
	semaphore.release()
	#return seqs, seq_len

def getRandom(indir):
	min_seq_len = 100000000000000000000000
	seq_count = {}
	random_count = {}
	for file in glob.glob(f'{indir}/*fa*'):
		if re.search('[fa|fasta]$', file):
			seq_data = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
			seq_num = len(seq_data)
			if seq_num < min_seq_len:
				min_seq_len = seq_num
			seq_count[file] = seq_num
		else:
			continue
	for key in seq_count:
		random_count[key] = round(seq_count[key]/min_seq_len)
	return random_count

def getIterFile(indir, iternumber, outdir, step):
	#random_count = getRandom(indir)
	for t in range(iternumber):
		if os.path.exists(f'{outdir}/step{step}.{t}.extract.fasta'):
			os.remove(f'{outdir}/step{step}.{t}.extract.fasta')
		with open(f'{outdir}/step{step}.{t}.extract.fasta', 'a') as fafile:
			for file in glob.glob(f'{indir}/*fa*'):
				if re.search('[fa|fasta]$', file):
					seq_data = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
					seq_count = len(seq_data)
					if seq_count > 10000:
						seq_ratio = 0.0625
					elif seq_count > 5000:
						seq_ratio = 0.125
					elif seq_count > 1000:
						seq_ratio = 0.25
					elif seq_count > 250:
						seq_ratio = 0.5
					else:
						seq_ratio = 1
					#print (f'Extract {random_count[file]*seqcount} sequence from {file}')
					random_seq = random.sample(list(seq_data.keys()),round(seq_ratio*seq_count))
					for seqid in random_seq:
						fafile.write(f'>{seqid}\n{seq_data[seqid].seq}\n')
				else:
					continue

def combineSeq(outdir, step):
	combine_file = open(f'{outdir}/step{step+1}.combine.fasta', 'w')
	for file in glob.glob(f'{outdir}/step{step}.*.cons.fasta'):
		with open(file) as infile:
			for line in infile:
				combine_file.write(line)
	combine_file.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-i","--indir", help="input fasta data path")
	parser.add_argument("-o","--outdir", help="output data path, do not same as input data path")
	parser.add_argument("-p","--parallel", type=int, help="max tasks analysis")
	parser.add_argument("-t","--iternumber", type=int,help="iter number for consense output")
	args = parser.parse_args()
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)
	maxtask = args.parallel
	semaphore = multiprocessing.Semaphore(maxtask)
	step = 1
	getIterFile(args.indir, args.iternumber, args.outdir, step)
	for fa in glob.glob(f'{args.outdir}/step1.*.extract.fasta'):
		t = os.path.basename(fa).split('.')[1]
		processes = []
		p = Process(target=seqAlign, args=(fa, f'{args.outdir}/step{step}.{t}.aln.fasta', args.outdir, step, t, semaphore))
		p.start()
		processes.append(p)
	for i in processes:
		p.join()
		#getConsenseSeq(seqs, seq_len, args.outdir, step, t)
	combineSeq(args.outdir, step)
	step += 1
	seqAlign(f'{args.outdir}/step{step}.combine.fasta', f'{args.outdir}/step{step}.combine.aln.fasta', args.outdir, step, t, semaphore)
	#getConsenseSeq(seqs, seq_len, args.outdir, step, 0)