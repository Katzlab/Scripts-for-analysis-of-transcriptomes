renameGBfasta.py README

Files downloaded from genbank come with a messy name that includes the gi number, genbank 
number, and description which can include spaces, odd characters, etc.
It is not useful to use these fastas directly in any phylogenetic software.
This script takes the names and renames the sequences with the gi number, preceeded
by a taxon code (unique for each taxon.)  This naming structure is much more tractable
in downstream applications.

You are able to preceed this code with more information, either two letters from another 
part of the ncbi taxonomy or a constant string of your choice.


#############################################################################
Useage is "python renameGBfasta.py <yourFastaFile>
#############################################################################
output is renamed file and a file containing the code key
#############################################################################
This script takes a fasta file downloaded from genbank and renames the sequences with a standard naming format based on the NCBI taxonomy for that taxon.
The standard is to rename each sequences with a 4 letter species code plus the gi number unique to that sequence.
If this code is not unique, it adds a letter to the end of the name. If there are too many non-unique codes, this script will fail .
#############################################################################
#############################################################################


