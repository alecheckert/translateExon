'''
printFullProtein.py -- given a translated transcript file, prints the
full amino acid sequence of the protein. The transcript file must be
translated.
'''
import pandas as pd
import argparse

def printFullProtein(filename, column='protein'):
	f = pd.read_csv(filename)
	full_protein = ''
	for i in f.index:
		if type(f.ix[i,column])==type(''):
			full_protein = full_protein + f.ix[i,column]
	print full_protein

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='given a transcript file, prints the sequence of the full protein')
	parser.add_argument('infile', type=str, help='name of the transcript file, which must contain a translated protein column and be sorted by rank.')
	parser.add_argument('-c', '--column', type=str, help="name of the column containing the protein sequences. default ``protein''.", default='protein')

	args = parser.parse_args()
	printFullProtein(args.infile, column=args.column)
