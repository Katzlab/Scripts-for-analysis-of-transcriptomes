'''
This script will search for SSU sequences from the taxa you input.
It will then use the results of the search to pull out closely related environmental seqeunces
and remove closely related sequences based on the RID percent you choose.

searchgb(stermList, mindate) does a text search for the terms given by the user
envSeqs() then blasts the sequences against environmental sequences and writes them to 

envseqs.fas also writes names out to the file 'sequence_names' along with their descriptions...so
you can see what they are.

It then concatenates and runs it through Rid to remove closely reated sequencs (not sure if
this is really working well!)
'''


from Bio import Entrez
Entrez.email = 'jgrant@smith.edu'
import string
import re
import sys
import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from datetime import date


def searchgb(stermList, mindate):
	usedID = []
	#try:
	out = open(stermList + '_searchdb.fas','a')
	print 'searching genbank'
	for sterm in stermList:
		searchTerm = '(' + sterm + '[Organism] OR '  + sterm +  '[All Fields]) AND ("' + str(mindate) + ' "[MDAT] : "' + str(date.today())+ ' "[MDAT])'
		handle = Entrez.esearch(db="nucleotide", retmax=10000, term= sterm)
	
		record = Entrez.read(handle)
		print record["Count"]
		for ID in  record["IdList"]:
			if ID not in usedID:
				usedID.append(ID)
			
				handle = Entrez.efetch(db="nucleotide", id=ID, rettype="fasta" )
				record = SeqIO.read(handle,"fasta")
				if re.search('18S', record.description):
					if re.search('5.8S', record.description):
						print "5.8S: %s" % record.description
					else:
					
						out.write('>' + record.id.split('|')[3] + '\n')
						out.write(str(record.seq) + '\n')
				elif re.search('small subunit', record.description):
					if re.search('5.8S', record.description):
						print "5.8S: %s" % record.description
					else:
					
						out.write('>' + record.id.split('|')[3] + '\n')
						out.write(str(record.seq) + '\n')

	
	out.close()
	
def envSeqs(inputfile):
	IDeDict = {}
	out2 = open(inputfile + '_envseqs.fas','a')
	out3 = open(inputfile + '_sequence_names','a')
	out4 = open(inputfile + '_gb','a')
	records = SeqIO.parse(open(inputfile),format="fasta")
	entrezQuery = '(environmental samples[filter] OR metagenomes[orgn])' # AND ("' + str(mindate) + ' "[MDAT] : "' + str(date.today())+ ' "[MDAT])'
	for record in records:
		IDlist = []
		#try:
		#result_handle = NCBIWWW.qblast("blastn", "nr", record.format("fasta"),  minidentity =  ".95")
		
		result_handle = NCBIWWW.qblast("blastn", "nr", record.format("fasta"), ncbi_gi=False, format_type="XML", entrez_query= entrezQuery ,expect = "2e-20", hitlist_size = "10000") #change evalue cut off to 2e-10

		#print result_handle.read()
		blast_records = NCBIXML.parse(result_handle)


		for blast_record in blast_records:
			#print blast_record.id
			if blast_record.descriptions:
				for alignment in blast_record.alignments:
					ID  = alignment.accession
					if ID not in IDlist:
						IDlist.append(ID)
						
						IDeDict[ID] = blast_record.alignments[0].hsps[0].expect ##
		#except:
		#	writelog("blast error: " +  record.id)
		out3.write( "**blast for " +  record.id + ':\n')
		#try:
		for ID in IDlist:	
			#handle = Entrez.efetch(db="nucleotide", id=ID, rettype="fasta")
			handle2 = Entrez.efetch(db="nucleotide", id=ID, rettype="gb")
			
			#record = SeqIO.read(handle,"fasta")
			record2 = SeqIO.read(handle2,"gb")
			
			name = record2.description.split()[1][0:2] + '_' + record2.description.split()[2][0:3] + "_" + record2.id
			
			out3.write( name + ':' + str(IDeDict[ID]) + ':' + record2.description + '\n') ##
			out2.write('>' + name + '\n')
			out2.write(str(record2.seq) + '\n')
			out4.write(str(record2))

		#except:
		#	writelog("efetch error: " +  record.id)
	
	
	out2.close()
	out3.close()
	out4.close()			
def writelog(input):
	outfile = open('log.txt','a')
	outfile.write(input + '\n')
	outfile.close

def main():
	print "This script will search for SSU sequences from the taxa you input." 
	print "It will then use the results of the search to pull out closely related environmental seqeunces."
	print "and remove closely related sequences based on the RID percent you choose. \n"
	search = raw_input('Do you need to do a text search? ')
	if search[0] == 'y':
	
		update = raw_input('Are you updating a file? If so, enter the date since your last update (YYYY/MM/DD) or else, hit enter. ')
		if update == "":
			mindate = 1986/01/01 

		else:
			mindate = update

		
		sterm = raw_input('Enter the search terms you would like to use, separated by commas (i.e. "oligotrichia,choreotrichia") ')	
		stermList = sterm.split(',')
		inputfile = stermList + '_searchdb.fas'
		searchgb(stermList, mindate)
	elif search[0] == 'n':
		inputfile = raw_input('What is the file of SSU sequences? ')
	else:
		print 'Please answer yes or no. '
		main()
	i = raw_input('What percent would you like to RID at? (hit return for default of 99%) ')
	try:
		num = float(i) + 1
	except TypeError:
		print 'Your input must be a number.  Try again. '
		main()
	except ValueError:
		i = ""	
	if i == "":
		RIDin = 99.00
	else:
		RIDin = float(i)
	
	x = raw_input('What is the minimum identity cut off? (hit return for default of 60%) ')
	try:
		num = float(x) + 1
	except TypeError:
		print 'Your input must be a number.  Try again. '
		main()
	except ValueError:
		x = ""
	if x == "":
		idmin = 60.00
	else:
		idmin = float(x)
	
	y = raw_input('What is the minimum length cut off? (hit return for default of 350 nt) ')
	try:
		num = float(y) + 1
	except TypeError:
		print 'Your input must be a number.  Try again. '
		main()
	except ValueError:
		y = ""	
	if y == "":
		lenmin = 350
	else:
		lenmin = int(y)
	z = raw_input('Hit return when you are ready to continue. ')
	try:
		writelog("Sterm = " + sterm)
		writelog("Update date = " + str(update))
	except:
		a='a'
	writelog("Rid identity = " + str(RIDin))
	writelog("min identity = " + str(idmin))
	writelog("min length = " + str(lenmin))
	
	
	
	#
	envSeqs(inputfile)
	os.system('cat ' + inputfile + ' ' + inputfile + '_envseqs.fas  > ready2RidFast.fas')
	

	import Rid4
	Rid4.rid(RIDin,idmin, lenmin)


main()