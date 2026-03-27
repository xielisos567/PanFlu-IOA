#!/bin/python3

import re
import os
import glob
import shutil
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

def getConsenseSeq(seqs, seqlen, outdir, faname, outfile, t):
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
	data_cons.to_csv(f'{outdir}/{faname}.{t}.stat.xls', sep='\t', index=False)
	with open(f'{outfile}', 'a') as outdata:
		seq = ''.join(data_cons['Cons_seq'].to_list())
		outdata.write(f'>{faname}.{t}.cons\n{seq}\n')

def seqAlign(infile, outfile, outdir, faname, t):
	mafft_cline = MafftCommandline(input=infile)
	aln_res, aln_error = mafft_cline()
	with open(f'{infile}.aln.fa', "w") as handle:
		handle.write(aln_res)
	align = AlignIO.read(f'{infile}.aln.fa', "fasta")
	seqs = defaultdict(list)
	seq_len = 0
	for recode in align:
		seqs[recode.id] = list(recode.seq)
	seq_len = len(seqs)
	#shutil.move(f'{infile}.aln.fa', f'{outdir}/tmp')
	data = pd.DataFrame.from_dict(seqs)
	data_cons = pd.DataFrame()
	data_cons['Cons_seq'] = data.apply(findMaxFreqAA, axis=1)
	data_cons['Cons_count'] = data.apply(countMaxFreqAA, axis=1)
	data_cons['Cons_ratio'] = data_cons["Cons_count"]*100/seq_len
	data_cons = data_cons.drop(data_cons[data_cons["Cons_seq"]=='-'].index)
	#if os.path.exists(f'{outdir}/{t}.step{step}.stat.xls'):
	#	os.remove(f'{outdir}/{t}.step{step}.stat.xls')
	#	os.remove(f'{outdir}/{t}.step{step}.cons.fasta')
	data_cons.to_csv(f'{outdir}/tmp/{faname}.{t}.stat.xls', sep='\t', index=False)
	#shutil.move(f'{outdir}/{faname}.{t}.stat.xls', f'{outdata}/tmp')
	with open(f'{outfile}', 'a') as outdata:
		seq = ''.join(data_cons['Cons_seq'].to_list())
		outdata.write(f'>{faname}.{t}.cons\n{seq}\n')
	#getConsenseSeq(seqs, seq_len, outdir, faname, t)
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

def getSeqAln(fa, outdir, finaldir, max_seq_count, semaphore):
	semaphore.acquire()
	fa_body = os.path.basename(fa).split('.')[0]
	seq_data = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))
	seq_count = len(seq_data)
	if seq_count<max_seq_count:
		shutil.copy(fa, f'{finaldir}/{fa_body}.fa')
	else:
		if os.path.exists(f'{outdir}/{fa_body}.fa'):
			os.remove(f'{outdir}/{fa_body}.fa')
		for t in range(int(10*seq_count/50)):
			random_seq = random.sample(list(seq_data.keys()),round(max_seq_count))
			tmp_file = open(f'{outdir}/tmp/{fa_body}.{t}.fa', 'w')
			for seqid in random_seq:
				tmp_file.write(f'>{seqid}\n{seq_data[seqid].seq}\n')
			seqAlign(f'{outdir}/tmp/{fa_body}.{t}.fa', f'{outdir}/{fa_body}.fa', outdir, fa_body, t)
	semaphore.release()

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

def combineSeq(outdir, indir):
	combine_file = open(f'{outdir}/final.combine.fasta', 'w')
	for file in glob.glob(f'{indir}/*fa*'):
		if file == f'{outdir}/final.combine.fasta':
			continue
		with open(file) as infile:
			for line in infile:
				combine_file.write(line)
	for file in glob.glob(f'{outdir}/*fa*'):
		if file == f'{outdir}/final.combine.fasta':
			continue
		with open(file) as infile:
			for line in infile:
				combine_file.write(line)
		shutil.move(file, f'{outdir}/tmp')
	combine_file.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-i","--indir", help="input fasta data path")
	parser.add_argument("-o","--outdir", help="output data path, do not same as input data path")
	parser.add_argument("-p","--parallel", type=int, help="max tasks analysis")
	parser.add_argument("-n","--maxnumber", help="max random sequence choosed each time [default=50,250]")
	args = parser.parse_args()
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)
		os.makedirs(f'{args.outdir}/step1/tmp')
		os.makedirs(f'{args.outdir}/step2/tmp')
		os.makedirs(f'{args.outdir}/final/tmp')
	maxtask = args.parallel
	maxnum1, maxnum2 = args.maxnumber.split(',')
	semaphore = multiprocessing.Semaphore(maxtask)
	step = 1
	#getIterFile(args.indir, args.iternumber, args.outdir, step)
	# first cycle
	for fa in glob.glob(f'{args.indir}/*fa*'):
		processes = []
		p = Process(target=getSeqAln, args=(fa, f'{args.outdir}/step1', f'{args.outdir}/final', int(maxnum1), semaphore))
		p.start()
		processes.append(p)
	for i in processes:
		p.join()
	# second cycle
	for fa in glob.glob(f'{args.outdir}/step1/*fa*'):
		processes = []
		p = Process(target=getSeqAln, args=(fa, f'{args.outdir}/step2', f'{args.outdir}/final', int(maxnum2), semaphore))
		p.start()
		processes.append(p)
	for i in processes:
		p.join()
	# final aln
	combineSeq(f'{args.outdir}/final', f'{args.outdir}/step2')
	seqAlign(f'{args.outdir}/final/final.combine.fasta', f'{args.outdir}/final/final.combine.consense.fasta', f'{args.outdir}/final/', 'final', 0)
	'''
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
	'''
