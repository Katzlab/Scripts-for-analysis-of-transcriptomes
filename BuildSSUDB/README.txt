BuildSSUDB.py README

This script will do a text search for SSU sequences based on user input entrez query (i.e. taxon name)
Downloaded files will then be used to blast against environmental samples.  The original
file plus environmental sequences will then be concatenated and 'Ridded'.  Ridding means 
sequences that do not meet the identity min and max and length cut-offs (provided by the user) will 
be removed.

User is asked if they want to do a text search or they can provide a fasta file of SSUs to
be used to pull environmental seqeunces

User is asked for min identity match percent and minimum sequence length.

User is asked for the Rid percent.  Sequences with a greater percent identity will be purged.




#############################################################################
Useage is "python BuildSSUDB.py"
Rid4.py  and the input fasta file (if any) must be in the same directory
#############################################################################
There are several output files.  The useful ones are:
InputFileName_envseqs.fas: The fasta file of environmental sequences
InputFileName_sequence_names: A file with the gb number and complete definition (name) field
InputFileName_gb: The genbank formatted files for each blast hit from the env. seq. search
WhoRiddedWho: sequences and the list of sequences they ridded
Output99 (Output + rid%identity) Output file of ridded sequences ready to be aligned and treed out
#############################################################################
#############################################################################


