'''
concatenateCSV.py -- join a list of CSVs into a large CSV.
'''
import pandas as pd
import os
import argparse

def concatenateCSV(filenames, outfile):
	if len(filenames)==0:
		print "No files specified."
		exit(1)
	else:
		df = pd.read_csv(filenames[0])
		for f in filenames[1:]:
			try:
				df2 = pd.read_csv(f)
				df = df.append(df2,ignore_index=True)
			except (KeyError, ValueError, TypeError) as e3:
				continue
		df.to_csv(outfile,index=False)
		return df

def readFileList(infile, delimiter='\n'):
	'''
	Reads a file containing a list of files to be concatenated.

	INPUT
		infile: string, the file containing a list of files
		delimiter: string, character used to separate files in said master file

	RETURNS
		list of string, the filenames

	'''
	g = open(infile, 'r')
	glines = g.read().split(delimiter)
	glines = [i for i in glines if len(i)>0]
	g.close()
	return glines

def generateCSVListFromDirectory(dir_name):
	'''
	Given a directory, finds all of the CSVs in that directory and returns them as a 
	list. The directory must be in the path.

	INPUT
		dir_name : string, the directory containing the CSVs

	RETURNS
		list of string, the files in that directory

	'''
	x = os.listdir(dir_name)
	x = [i for i in x if '.csv' in i]
	x = ['%s/%s' % (dir_name, i) for i in x]
	return x

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='concatenate similar CSVs into a single CSV')
	parser.add_argument('-i', '--infile', type=str, help='file containing a list of CSVs to be concatenated')
	parser.add_argument('-o', '--outfile', type=str, help='CSV to write the concatenated information to', required=True)
	parser.add_argument('-s', '--delimiter', type=str, help='delimiter used in the infile. default newline.', default='\n')
	parser.add_argument('-d', '--directory', type=str, help='directory containing the files to be concatenated')

	args = parser.parse_args()
	if args.directory:
		try:
			target_files = generateCSVListFromDirectory(args.directory)
			concatenateCSV(target_files, args.outfile)
			print "Finished"
		except ValueError:
			print "Encountered directory error."
			exit(1)
	elif args.infile:
		try:
			target_files = readFileList(args.infile, delimiter=args.delimiter)
			concatenateCSV(target_files, args.outfile)
			print "Finished"
		except ValueError:
			print "Encountered target file list error."
			exit(1)
	else:
		print "Incorrect input; see usage."
		exit(1)
