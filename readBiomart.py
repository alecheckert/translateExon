'''
readBiomart.py -- read files of the type output by ENSEMBL BioMart,
an online biological data mining application.

BioMart is a way to download a set of desired sequences with various header
information. These sequences are in modified FASTA format, with the header
information on the FASTA identity line. The headers are `|'-delimited. The problem
is that BioMart export files may not contain all of the header fields if the
last header fields are empty. For this reason, you should always specify
header fields with a universally populated field (such as transcript ID, or
exon ID) last. This guarantees that each FASTA header line contains all of the
header columns.

*readBiomart* takes a BioMart export file and a header file specifying the
header columns. The header file is a single-line `,'-delimited file with
the names of the headers. This program returns a dataframe containing the
information in the file.
'''
import pandas as pd
import argparse
import math

def readBiomart(biomart_file, header_file):
	'''
	Reads a BioMart export file.

	INPUT:
		biomart_file : (string) a file in ENSEMBL BioMart export format, assumed
			to be FASTA with `|'-delimited header columns on the identity line.

		header_file : (string) a file containing the names of the header columns,
			which should be a single-line, `,'-delimited list of the column names.
			Must match with BioMart file. Does not contain the ``sequence'' field,
			which is added separately.

	RETURNS:
		a pandas DataFrame object

	'''
	print "Reading BioMart file..."
	f=open(biomart_file,'r')
	flines=f.read().split('\n')
	f.close()

	print "Reading columns..."
	h=open(header_file,'r')
	columns=h.read().split('\n')[0].split(',')
	columns=[i for i in columns if len(i)>0]
	columns=columns+['sequence']
	h.close()

	print "Making dataframe..."
	out=pd.DataFrame(columns=columns)
	for i in range(len(flines)):
		if len(flines[i])>0 and flines[i][0]=='>':
			fline_0 = flines[i].replace('>','').split('|')
			cds = ''
			j = i+1
			while (len(flines[j])!=0) and (flines[j][0]!='>'):
				cds += flines[j]
				j+=1
			row=fline_0+[cds]
			rowSeries = pd.Series(row, index=columns)
			out=out.append(rowSeries, ignore_index=True)
	return out

def writeBiomart(df, outname):
	'''
	Writes data into a CSV.

	INPUT:
		df : a pandas DataFrame object
		outname : a string containing the name of the CSV to write to
	RETURNS:
		<None>

	'''
	print "Writing to CSV..."
	df.to_csv(outname, index=False)

def zeroColumn(df, columnNames):
	'''
	Replaces blank values in a DataFrame with zeros, and converts to int if
	possible.

	INPUT:
		df : a pandas DataFrame object
		columnNames : a list of the column names to correct
	RETURNS:
		a pandas DataFrame object

	'''
	columnNames = [i for i in columnNames if i in df.columns]
	for columnName in columnNames:
		for i in df.index:
			if df.ix[i,columnName]=='':
				df.ix[i,columnName]=0
		else:
			df.ix[i,columnName]=try_int(df.ix[i,columnName])
	return df
	
def try_int(arg):
	'''
	Attempts to convert an argument to integer.
	'''
	if (type(arg)==type('')) or (math.isnan(arg)):
		try:
			return int(arg)
		except (KeyError, ValueError) as e2:
			return arg
	else:
		return arg

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='read FASTA files in ENSEMBL BioMart output format and write to CSV')
	parser.add_argument('-i', '--infile', type=str, help='FASTA file containing the BioMart sequence information')
	parser.add_argument('-j', '--headerfile', type=str, help='CSV containing the names of the FASTA columns')
	parser.add_argument('-o', '--outfile', type=str, help='file to write results to')
	parser.add_argument('-z', '--zero', action='store_true', help='convert all values to int if possible and replace empty values with 0')
	args = parser.parse_args()
	f=readBiomart(args.infile, args.headerfile)
	if args.zero:
		f = zeroColumn(f, f.columns)
	writeBiomart(f,args.outfile)
	print "Finished"
