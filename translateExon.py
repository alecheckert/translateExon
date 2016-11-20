'''
translateExon.py -- translate exons with the appropriate phase information
'''
import math
import os
import argparse
import numpy as np
import pandas as pd
from codonTable import codonTable

def translate(cds, startPhase, endPhase, find_orfs=True):
	'''
	Translate a nucleotide sequence into a protein sequence.

	This program deals with three situations:
		1. Nonnegative start phase. Translate the sequence until the first
			stop codon.
		2. Negative start phase and nonnegative end phase. Translate the
			sequence in reverse until the first start codon.
		3. Negative start phase and negative end phase. If *find_orfs*, look
			for and translate the longest open reading frame. Otherwise return
			an empty string.

	INPUT
		cds : string, nucleotide sequence to be translated
		startPhase : int, starting phase (number of nucleotides on first codon
			that lie on the previous exon)
		endPhase : int, ending phase (number of nucleotides of last exon that
			lie on the current exon).

	RETURNS
		string, translated sequence

	'''
	if startPhase>=0:
		if startPhase==0:
			pass
		else:
			cds = cds[(3-startPhase):]
		result = ''
		for i in range(len(cds)/3):
			codon=cds[(i*3):((i*3)+3)]
			if codonTable[codon]=='X':
				break
			result += codonTable[codon]
		return result
	elif (startPhase<0) and (endPhase>=0):
		result = ''
		if endPhase==0:
			pass
		else:
			cds = cds[:-endPhase]
		while len(cds)>=3:
			aa = codonTable[cds[-3:]]
			result = aa + result
			if aa=='M':
				#break #trying something out
				pass
			cds = cds[:-3]
		result = result[result.find('M'):]
		return result
	else:
		if find_orfs:
			if findORF(cds):
				start_pos, stop_pos, orf_length = findORF(cds)
				translate(cds[start_pos:stop_pos], startPhase=0, endPhase=0)
			else:
				return ''
		else:
			return ''

def findORF(sequence):
	'''
	Find the start and stop of the longest ORF.

	INPUT
		sequence : string, the nucleotide sequence to be analyzed

	RETURNS
		(start, stop) : int, start and stop positions of the ORF.

	'''
	orfs = []
	for i in range(len(sequence)-3):
		if codonTable[sequence[i:(i+3)]]=='M':
			start_pos = i
			stop_pos = i
			new_seq = sequence[i:]
			for j in range(1, len(new_seq)/3-1):
				codon = new_seq[(j*3):((j+1)*3)]
				if codonTable[codon]=='X':
					stop_pos=i+(j*3)
					break
			if start_pos==stop_pos:
				pass
			else:
				orfs.append((start_pos, stop_pos, stop_pos-start_pos))
	if len(orfs)>0:
		orfs = sorted(orfs, key=lambda i: i[2], reverse=False)
		return orfs[0]
	else:
		return False

def newFindStartExon(transcript_df, start_phase='startPhase', end_phase='endPhase', rank='rank'):
	'''
	Get the rank of the exon containing the start codon. Revised version: the program looks for
	the first exon that satisfies the first criterion, then moves onto the next criterion:
		1. Exon has (start_phase, end_phase) = (-1, >=0)
		2. Exon has (start_phase, end_phase) = (>=0, x) and all previous exons are (-1, -1)
		3. All exons are (-1, -1) and exon contains the longest ORF in the transcript

	INPUT
		transcript_df : pandas DataFrame, with exons in rank order
		start_phase, end_phase, rank : names of the columns in *transcript_df*

	RETURNS
		int, the index of the start exon

	'''
	transcript_df = transcript_df.sort_values(by=rank)
	transcript_df = transcript_df.set_index(rank, drop=False)
	highest_rank = max(transcript_df.index)

	if transcript_df.ix[1,'startPhase']>=0:
		return 1
	elif all([transcript_df['startPhase']<0][0]) and all([transcript_df['endPhase']<0][0]):
		exon_orfs = [findORF(transcript_df.ix[i,'sequence']) for i in transcript_df.index]
		for i in range(len(exon_orfs)):
			if exon_orfs[i]==False:
				exon_orfs[i]=(0,0,0)
		if (not exon_orfs) or (not exon_orfs[0]) or (len(exon_orfs)==0):
			return False
		else:
			return exon_orfs.index(max(exon_orfs, key=lambda i: i[2]))+1
	else:
		for i in transcript_df.index:
			if (transcript_df.ix[i,'startPhase']<0) and (transcript_df.ix[i,'endPhase']>=0):
				return i
		return False

def newTranslateDF(transcript_df, start_phase='startPhase', end_phase='endPhase', rank='rank'):
	'''
	Adds a new column, ``protein'', to the DataFrame *transcript_df*.

	INPUT
		transcript_df : pandas DataFrame, indices are exons corresponding to that transcript. Must
			have phase and rank information.

		start_phase, end_phase, rank : strings, names of the corresponding columns in *transcript_df*

	RETURNS
		pandas DataFrame, a copy of *transcript_df* with the new ``protein'' column

	'''
	start_exon = newFindStartExon(transcript_df, start_phase=start_phase, end_phase=end_phase, rank=rank)
	if not start_exon:
		peptides = pd.Series(['' for i in transcript_df.index])
		transcript_df['protein']=peptides
		return transcript_df
	peptides = pd.Series([])
	for i in transcript_df.index:
		if i < start_exon:
			peptides.ix[i]=''
		elif i==start_exon:
			peptides.ix[i]=translate(transcript_df.ix[i,'sequence'], transcript_df.ix[i,start_phase], transcript_df.ix[i,end_phase])
		else:
			peptide = translate(transcript_df.ix[i,'sequence'], transcript_df.ix[i,start_phase], transcript_df.ix[i,end_phase])
			if i>1: #get the first codon, part of which lies on the previous exon
				if transcript_df.ix[i,start_phase]>0:
					_start_phase = transcript_df.ix[i,start_phase]
					prev_cds = transcript_df.ix[i-1,'sequence']
					first_codon=prev_cds[-_start_phase:]
					first_codon = first_codon + transcript_df.ix[i,'sequence'][:(3-_start_phase)]
					peptide = codonTable[first_codon] + peptide
				else:
					pass
			peptides.ix[i]=peptide
	transcript_df['protein']=peptides
	print transcript_df ##debugging
	print "START EXON: %d\n\n" % start_exon ##debugging
	return transcript_df

def translateTranscriptFile(transcriptFile, rank='rank', start_phase='startPhase', end_phase='endPhase'):
	'''
	Reads a file containing exon sequences corresponding to one transcript and
	translates each of the exons, adding a ``protein'' column to the file.

	INPUT
		transcriptFile : Pandas-type CSV with the exon sequences to be translated

	RETURNS
		<None> (writes to file)

	'''
	f = pd.read_csv(transcriptFile)
	f = f.sort_values(by=rank)
	f = f.set_index(rank, drop=False)
	f = newTranslateDF(f, rank=rank, start_phase=start_phase, end_phase=end_phase)
	f.to_csv(transcriptFile, index=False)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='given a CSV with exon sequences, translate them with the correct phase')
	parser.add_argument('infile', type=str, help='file or directory of files with exons corresponding to a single transcript')
	parser.add_argument('--startphase', type=str, help="name of the start phase column. default: ``startPhase''.", default='startPhase')
	parser.add_argument('--endphase', type=str, help="name of the end phase column. default: ``endPhase''.", default='endPhase')
	parser.add_argument('--rank', type=str, help="name of the rank column. default: ``rank''.", default='rank')

	args = parser.parse_args()
	if os.path.isdir(args.infile):
		files = [i for i in os.listdir(args.infile) if '.csv' in i]
		for u in files:
			try:
				translateTranscriptFile('%s/%s' % (args.infile,u), rank=args.rank, start_phase=args.startphase, end_phase=args.endphase)
			except pandas.io.common.CParserError:
				continue
	else:
		translateTranscriptFile(args.infile, rank=args.rank, start_phase=args.startphase, end_phase=args.endphase)
