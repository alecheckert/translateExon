'''
getLongestTranscript.py -- given a set of ENSEMBL gene IDs, transcript IDs, and transcript
lengths, takes only the longest transcript for each gene. Useful as a starting point in
exon analyses.
'''
import pandas as pd
import numpy as np
import argparse

def getLongestTranscript(filename, outname, columns=['gene_id', 'transcript_length', 'transcript_id']):
	'''
	INPUT
		filename : a CSV containing a list of transcripts, their corresponding genes, and their
			lengths. Must contain the column names in *columns*.

		columns : a list of the column names.

	RETURNS
		pandas DataFrame object
		writes to CSV

	'''
	print "Reading file..."
	f = pd.read_csv(filename)
	gene_id, transcript_length, transcript_id = columns[0], columns[1], columns[2]
	unique_genes = np.unique(f[gene_id])
	df = pd.DataFrame(columns=f.columns)
	print "Searching for unique transcripts..."
	for unique_gene in unique_genes:
		gene_df = f[f[gene_id]==unique_gene]
		max_transcript = gene_df[gene_df[transcript_length]==max(gene_df[transcript_length])]
		if len(max_transcript)>1:
			max_transcript=max_transcript.iloc[0]
		df = df.append(max_transcript, ignore_index=True)
	print "Writing output..."
	df.to_csv(outname,index=False)
	return df
	
if __name__=='__main__':
	parser = argparse.ArgumentParser(description='take the longest transcript for each gene')
	parser.add_argument('infile', type=str, help="CSV input file. must contain `transcript_id', `gene_id', and `transcript_length' columns")
	parser.add_argument('outfile', type=str, help="file to write results to")
	args = parser.parse_args()
	getLongestTranscript(args.infile, args.outfile)
	print "Finished"
