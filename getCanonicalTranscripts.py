'''
getCanonicalTranscripts.py -- retrieve the ENSEMBL stable IDs for the ENSEMBL
canonical transcripts, given a list of ENSEMBL gene IDs. Note: the latest
release of ENSEMBL as of 11/19/2016 is 86.
'''
from cogent.db.ensembl import Genome
human=Genome('human', 86, None)
import pandas as pd
import argparse

def readGenes(filename):
	'''
	Reads a list of ENSEMBL genes. The file should be a `\n'-delimited
	list of ENSEMBL gene IDs.
	'''
	f = open(filename,'r')
	flines = f.read().split('\n')
	f.close()
	return [i for i in flines if len(i)>0]

def getCanonicalTranscripts(geneList, startIndex=0, stopIndex=None):
    '''
    Given a list of ENSEMBL gene IDs, returns a DataFrame that relates
    these gene IDs to the transcript IDs of their canonical transcripts.
    Assumes that the ENSEMBL stable IDs in *geneList* are human.
    '''
    if not startIndex:
        startIndex=0
    if (not stopIndex) or (stopIndex == 0):
        stopIndex=len(geneList)
    result=pd.DataFrame(columns=['geneID', 'transcriptID'])
    for geneID in geneList[startIndex:stopIndex]:
        try:
            geneObj=human.getGeneByStableId(StableId=geneID)
            transcriptID=geneObj.CanonicalTranscript.StableId
            print geneID, transcriptID
            newRow = pd.Series([geneID, transcriptID], index=['geneID', 'transcriptID'])
            result=result.append(newRow,ignore_index=True)
        except AttributeError:
            continue
    return result
 
def writeCanonicalTranscript(geneList, outName, startIndex=None, stopIndex=None):
    '''
    Performs *getCanonicalTranscripts* on *geneList*, then writes to a file.
    '''
    geneTranscriptDF=getCanonicalTranscripts(geneList, startIndex=startIndex, stopIndex=stopIndex)
    geneTranscriptDF.to_csv(outName,index=False)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='find the canonical transcripts for a list of ENSEMBL gene IDs')
	parser.add_argument('-i', '--infile', type=str, help='file containing newline-delimited list of genes', required=True)
	parser.add_argument('-o', '--outfile', type=str, help='file to write the canonical transcript list to', required=True)
	parser.add_argument('-s', '--startindex', type=int, help='position in the gene list to start at. Default 0.', default=0)
	parser.add_argument('-e', '--stopindex', type=int, help='position in the gene list to stop at. Default None.', default=0)

	args = parser.parse_args()
	f = readGenes(args.infile)
	writeCanonicalTranscript(f, args.outfile, startIndex=args.startindex, stopIndex=args.stopindex)
