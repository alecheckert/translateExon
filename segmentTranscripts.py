'''
segmentTranscripts.py -- split a CSV containing sequences from
many genes into separate files, one for each transcript
'''
import argparse
import os
import numpy as np
import pandas as pd

def segmentTranscripts(filename, transcript_column='transcript_id', sort_by='rank', out_prefix=None):
	'''
	Split the lines in a sequence CSV into separate files, one for each transcript.
	Assumes that each index in the CSV has a single transcript in its ID field.

	INPUT
		filename : string, name of the file containing the sequence information
	
		transcript_column : string, name of the column containing the transcript
			IDs
		
		out_prefix : string, prefix to be appended to the transcript names in 
			generating the output filenames

	RETURNS
		<None>

	'''
	f = pd.read_csv(filename)
	if transcript_column not in f.columns:
		print "Cannot find the transcript ID column ``%s'' in the file. File has columns %r" % (transcript_column, f.columns)
		exit(1)
	else:
		f = splitIndices(f, transcript_column)
		unique_transcripts = np.unique(f[transcript_column])
		for transcript in unique_transcripts:
			print transcript, out_prefix
			transcript_df = f[f[transcript_column]==transcript]
			transcript_df = transcript_df.sort_values(by=sort_by)
			transcript_df = transcript_df.set_index(sort_by, drop=False)
			if (not out_prefix) or (len(out_prefix)==0):
				transcript_df.to_csv('%s.csv' % transcript, index=False)
			elif out_prefix[-1]=='/':
				transcript_df.to_csv('%s/%s.csv' % (out_prefix, transcript), index=False)
			else:
				transcript_df.to_csv('%s_%s.csv' % (out_prefix, transcript), index=False)

def splitIndices(df, column, delimiter=';'):
	'''
	For indices that have more than one value at the attribute *column*,
	splits these into separate indices.

	INPUT
		df : pandas DataFrame, containing the data to be split
		column : string, name of the column to look for splits

	RETURNS
		pandas DataFrame

	'''
	for i in df.index:
		if delimiter in df.ix[i, column]:
			new_row_transcripts = df.ix[i, column].split(delimiter)
			for t in new_row_transcripts[1:]:
				rowSeries = df.ix[i]
				rowSeries[column]=t
				df = df.append(rowSeries, ignore_index=True)
			df.ix[i, column]=new_row_transcripts[0]
	return df

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='split a sequence CSV into separate CSVs, one for each transcript')
	parser.add_argument('infile', type=str, help='Pandas-format CSV containing the sequence information')
	parser.add_argument('-c', '--column', type=str, help="name of the column to split the file by. Default: ``transcript_id''", default='transcript_id')
	parser.add_argument('-o', '--outdir', type=str, help='directory to output to', default=None)
	parser.add_argument('-s', '--sort', type=str, help="name of the column to sort the values by. Default is ``rank''", default='rank')
	args = parser.parse_args()

	if args.outdir:
		if not os.path.isdir(args.outdir):
			os.mkdir(args.outdir)
		out_prefix = '%s/' % args.outdir
	else:
		out_prefix = None

	f = segmentTranscripts(args.infile, transcript_column=args.column, sort_by=args.sort, out_prefix=out_prefix)

	print "Finished"
	
