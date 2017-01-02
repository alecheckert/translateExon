===============================================================================
translateExon: a Python utility for translating protein-encoding exons in phase
===============================================================================

Motivation.
~~~~~~~~~~~

Eukaryotic proteins are frequently encoded in small pieces of information - exons - scattered over large regions of genomic space. The phase of an exon determines
its reading frame. This tool is designed to translate independent 
exon sequences with the appropriate phase information.

Workflow.
~~~~~~~~~~~~~~~~~~

The following workflow is recommended in order to analyze exons from individual
proteins:

Download a set of gene IDs from ENSEMBL's BioMart service. Use filters for protein-
coding genes and protein-coding transcripts, returning a dataset that contains
gene IDs, transcript IDs, and transcript lengths. Use the function
``getLongestTranscript.py`` to find the longest transcript corresponding to each
gene, generating a nonredundant list of canonical transcripts.

Submit the list of canonical transcripts to ENSEMBL BioMart, downloading exon
sequences in FASTA format. The following information is required for the
``translateExon.py`` program: rank, start phase, end phase, and transcript ID; 
any other information can be included but is not used by the program.

Use the function ``readBiomart.py`` to convert the BioMart-type FASTA files to
more tractable CSV files. Segment this CSV into individual CSVs with 
``segmentTranscripts.py``. Finally, apply the translation function ``translateExon.py``
to the directory containing the transcript CSVs. The result is a directory
containing the canonical transcript files of interest, each containing a list
of exon sequences, their phase information, their encoded peptide sequences,
and any other information.

Component Programs.
~~~~~~~~~~~~~~~~~~~

<To be added.>
