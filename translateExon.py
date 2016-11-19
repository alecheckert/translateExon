'''
translateExon.py -- translate exons with the appropriate phase information
'''
import math
import os
import argparse
import numpy as np
import pandas as pd
from codonTable import codonTable

def translate(cds, startPhase, endPhase, find_orfs=False):
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
		cds = cds[:-endPhase]
		while len(cds)>=3:
			aa = codonTable[cds[-3:]]
			result = aa + result
			if aa=='M':
				break
			cds = cds[:-3]
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

def findStartExon(transcript_df, start_phase='startPhase', end_phase='endPhase', rank='rank'):
	'''
	Get the rank of the exon containing the start codon.

	INPUT
		transcript_df : pandas DataFrame, containing the exons corresponding
			to a single transcript. Must have the rank, startPhase, and endPhase
			columns.

		start_phase, end_phase, rank : strings, names of the columns containing
			the relevant information

	RETURNS
		int, the rank of the exon with the start codon, or 0 if start codon is not
			found.

	'''
	transcript_df = transcript_df.sort_values(by=rank)
	transcript_df = transcript_df.set_index(rank, drop=False)
	highest_exon_rank = max(transcript_df.index)

	if len(transcript_df)==1:
		return 1
	for i in transcript_df.index:
		if transcript_df.ix[i,start_phase]<0 and transcript_df.ix[i,end_phase]>=0:
			return i
		elif transcript_df.ix[i,start_phase]<0 and i<highest_exon_rank and transcript_df.ix[i+1,end_phase]>=0:
			return i+1
	if all([i>=0 for i in transcript_df[start_phase].values.tolist()]) and ((transcript_df.ix[1,'sequence'][0:3]=='ATG') or (transcript_df.ix[1,'sequence'][0:3]=='AUG')):
		return 1
	if all([i>=0 for i in transcript_df[start_phase].values.tolist()]):
		return 1
	elif all([i<0 for i in transcript_df[start_phase].values.tolist()]):
		orf_lengths = []
		for i in transcript_df.index:
			if findOpenReadingFrame(transcript_df.ix[i,'sequence']):
				orf_lengths.append(findLargestOpenReadingFrame(transcript_df.ix[i,'sequence'])[2])
			else:
				orf_lengths.append(0)
		return orf_lengths.index(max(orf_lenths))+1
	for i in transcript_df.index:
		if transcript_df.ix[i,start_phase]==0 and transcript_df.ix[i,end_phase]<0:
			return i
	return 0


def findOpenReadingFrame(cds):
    '''
    If an open reading frame - an in-frame start codon and stop codon- exist
    in the sequence, returns it and its length.
    '''
    for i in range(len(cds)-6):
        if (cds[i:(i+3)]=='AUG') or (cds[i:(i+3)]=='ATG'):
            putative_orf = cds[i:]
            for j in range(len(putative_orf)/3-1):
                if codonTable[putative_orf[(j*3):((j+1)*3)]]=='X':
                    return j*3
    return False
 
def findLargestOpenReadingFrame(cds):
    '''
    Looks for open reading frames and returns the start position, stop position,
    and length of the longest one it can find. If it cannot find any open reading
    frames, returns *False*.
    '''
    orfs = []
    for i in range(len(cds)-6):
        if (cds[i:(i+3)]=='AUG') or (cds[i:(i+3)]=='ATG'):
            putative_orf=cds[i:]
            for j in range(len(putative_orf)/3-1):
                if codonTable[putative_orf[(j*3):((j+1)*3)]]=='X':
                    orfs.append((i, i+(j*3), j*3))
    if len(orfs)==0:
        return False
    orfs=sorted(orfs, key=lambda u: u[2])
    return orfs[-1]
 
def translate_df(transcript_df, start_phase='startPhase', end_phase='endPhase', rank='rank'):
    '''
    Adds a new column, ``protein'', to the DataFrame *transcript_df* that
    consists of the translated sequence for that exon.
    '''
    start_exon = findStartExon(transcript_df, rank=rank, start_phase=start_phase, end_phase=end_phase)
    peptides = pd.Series([])
    for i in transcript_df.index:
        if i < start_exon:
            peptides.ix[i]=''
        elif i==start_exon:
            if transcript_df.ix[i,end_phase]==-1:
                cds=transcript_df.ix[i,'sequence']
                start_pos, stop_pos, length = findLargestOpenReadingFrame(cds)
                cds=cds[start_pos:]
                peptides.ix[i] = translate(cds, startPhase=0, endPhase=0)
            else:
                cds=transcript_df.ix[i,'sequence']
                peptides.ix[i]=translate(cds,startPhase=-1,endPhase=transcript_df.ix[i,end_phase])
        else:
            if transcript_df.ix[i-1,end_phase]<0:
                peptides.ix[i] = ''
            else:
                cds=transcript_df.ix[i,'sequence']
                _start_phase=transcript_df.ix[i,start_phase]
                prev_cds=transcript_df.ix[i-1,'sequence']
                if _start_phase == 0:
                    prev_cds=''
                else:
                    prev_cds = prev_cds[(-_start_phase):]
                cds = '%s%s' % (prev_cds, cds)
                peptides.ix[i] = translate(cds, startPhase=0, endPhase=0)
    transcript_df['protein']=peptides
    print transcript_df
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
	f = translate_df(f, rank=rank, start_phase=start_phase, end_phase=end_phase)
	f.to_csv(transcriptFile, index=False)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='given a CSV with exon sequences, translate them with the correct phase')
	parser.add_argument('-i', '--infile', type=str, help='file with exons corresponding to a single transcript', required=True)
	parser.add_argument('--startphase', type=str, help="name of the start phase column. default: ``startPhase''.", default='startPhase')
	parser.add_argument('--endphase', type=str, help="name of the end phase column. default: ``endPhase''.", default='endPhase')
	parser.add_argument('--rank', type=str, help="name of the rank column. default: ``rank''.", default='rank')

	args = parser.parse_args()
	translateTranscriptFile(args.infile, rank=args.rank, start_phase=args.startphase, end_phase=args.endphase)
