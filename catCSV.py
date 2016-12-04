'''
catCSV.py -- concatenate multiple CSVs into a single CSV
'''
import pandas as pd
import argparse
import os

def catCSV(dir_name, outname, no_header=False, exclude=[]):
	'''
	First, checks to see if the directory *dir_name*. If so, finds all of the CSVs
	inside of it and opens them with pandas. The program then attempts to
	concatenate them. Files must all be of the same type.

	INPUT
		dir_name : string, name of the directory containing the CSVs
		outname : string, name of the CSV to write the result to
		no_header : bool, *True* if the target CSVs do not have a header line
		exclude : list of string, files to be excluded from the concatenation

	RETURNS
		Pandas DataFrame : the concatenated data

	'''
	if not os.path.isdir(dir_name):
		print "Could not find the directory %s in the current working directory." % dir_name
		exit(1)
	else:
		print "Found directory %s" % dir_name
		fs = ['%s/%s' % (dir_name, i) for i in os.listdir(dir_name) if ('.csv' in i) and (i not in exclude)]
		if len(fs) == 0:
			print "No CSV files found in directory %s." % dir_name
			exit(1)
		try:
			if no_header:
				f_data = [pd.read_csv(i, header=False) for i in fs]
			else:
				f_data = [pd.read_csv(i) for i in fs]
			result = pd.DataFrame(columns=f_data[0].columns)

			for f in range(len(f_data)):
				result = pd.concat([result, f_data[f]])
				print "Concatenated file %s" % fs[f]
			result.to_csv(outname,index=False)
			print "Writing to %s" % outname
			return result
		except (TypeError, ValueError, KeyError) as e3:
				print 'Could not concatenate the CSV files in directory %s' % dir_name
				exit(1)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='concatenate all CSVs in a given directory into a single CSV.')
	parser.add_argument('directory', type=str, help='name of the directory containing the CSVs to be concatenated. All CSVs must have the ``.csv'' suffix.')
	parser.add_argument('outfile', type=str, help='name of file to write concatenated CSVs to')
	parser.add_argument('-n', '--noheader', action='store_true', help='target CSVs do not have header lines.')
	parser.add_argument('-e', '--exclude', action='append', dest='files_to_exclude', help='files to exclude from the concatenation')
	args = parser.parse_args()
	catCSV(args.directory, args.outfile, no_header=args.noheader, exclude=args.files_to_exclude)
	print "Finished"
