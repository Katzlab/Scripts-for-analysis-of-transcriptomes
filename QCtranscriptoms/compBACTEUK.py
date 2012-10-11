
#########################################################################################
# useage python compBACTEUK.py <file_name>
# This script takes a fasta file and blasts against nr limited to bacteria and archaea and
# limited to eukaryotes.  It takes the ratio of evalues to determine whether the sequence
# is strongly euk, strongly bact/arch or undetermined.
# input - a fasta file
# output - blast result files bacterial.txt and eukaryotic.txt and reult file eval_comparison.csv
# This program blasts on line so it is vulnerable to time outs/other online issues.  Watch
# the terminal for sequences whose blasts fail!
#########################################################################################

import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import os



def blastBACTEUK(arg):
	out=open('bacterial.txt','a')
	out2=open('eukaryotic.txt','a')
	records = SeqIO.parse(open(arg), format="fasta")
	
	for record in records:
		try:
			name = record.id
			result_handleB = NCBIWWW.qblast("blastx", "nr", record.format("fasta"), ncbi_gi=False, descriptions= "1", alignments="1", format_type="XML", hitlist_size="1", entrez_query='(Bacteria[ORGN] OR Archaea[ORGN])')
			result_handleE = NCBIWWW.qblast("blastx", "nr", record.format("fasta"), ncbi_gi=False, descriptions= "1", alignments="1", format_type="XML", hitlist_size="1", entrez_query='(Eukaryota[ORGN])')

			blast_recordsB = NCBIXML.read(result_handleB)
			blast_recordsE = NCBIXML.read(result_handleE)

			if blast_recordsB.descriptions:
				print record.id

				name = record.id


				out.write(name + ',' + str(blast_recordsB.alignments[0].hsps[0].expect) + '\n')
			else:
				out.write(name + ', no hit'  + '\n')

			if blast_recordsE.descriptions:
				out2.write(name + ',' +  str(blast_recordsE.alignments[0].hsps[0].expect) + '\n')
			else:
				out2.write(name + ', no hit'  + '\n')
		except:
			errorout = open('errorlog.txt','a')
			error out.write('problem blasting ' + record.id + '\n')
			errorout.close()

	out.close()
	out2.close()


def compareEval():
	outfile = open('eval_comparison.csv','a')
	outfile.write('NAME,eval Bact,eval Euk,evalB/evalE,WHAT IS IT?\n')
	infileB=open('bacterial.txt','r').readlines()
	infileE=open('eukaryotic.txt','r').readlines()

	for lineB in infileB:

		name = lineB.split(',')[0]
		evalB = lineB.split(',')[1]
		try:
			evalB = float(evalB)
		except:
			evalB = evalB.strip()
		for lineE in infileE:
			if name == lineE.split(',')[0]:
				evalE = lineE.split(',')[1]
				
				try:
					evalE = float(evalE)
				except:
					evalE = evalE.strip()
				
				if evalE == 'no hit' and evalB == 'no hit':
					outfile.write(name + ',' + str(evalB) + ',' + str(evalE)  + ',N/A,UNDETERMINED\n')
				elif evalE == 'no hit' and evalB != 'no hit':
					outfile.write(name + ',' + str(evalB) + ',' + str(evalE)  + ',N/A,EUKARYOTIC\n')
				elif evalE != 'no hit' and evalB == 'no hit':
					outfile.write(name + ',' + str(evalB) + ',' + str(evalE)  + ',N/A,BACTERIAL\n')
				else:
					if evalE == 0.0:
						evalE = 1e-310 #smallest I can get without python reverting to 0.0 and not dividing
					if evalB == 0.0:
						evalB = 1e-310
					outfile.write(name + ',' + str(evalB) + ',' + str(evalE)  + ',' + str(evalB/evalE) + ',')
					if (evalB/evalE) > 10.00:
						outfile.write('EUKARYOTIC\n')
					elif (evalE/evalB) > 10.00:
						outfile.write('BACTERIAL\n')
					else:
						outfile.write('UNDETERMINED\n')
					break
def main():
	for arg in sys.argv[1:]:
		blastBACTEUK(arg)
		compareEval()


main()
