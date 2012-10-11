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
entrezQuery = '(environmental samples[filter] OR metagenomes[orgn])' # AND ("' + str(mindate) + ' "[MDAT] : "' + str(date.today())+ ' "[MDAT])'


def searchgb(stermList, mindate):
	usedID = []
	#try:
	out = open('seq_searchdb.fas','a')
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
	
def blast(record, IDlist, IDeDict,toBlast,listsize):	
	try:		
		result_handle = NCBIWWW.qblast("blastn", "nr", record.format("fasta"), ncbi_gi=False, format_type="XML", entrez_query= entrezQuery , expect = "2e-50", hitlist_size = listsize) #change evalue cut off to 2e-10
		blast_records = NCBIXML.parse(result_handle)

						
								
		for blast_record in blast_records:
			if blast_record.descriptions:
				for alignment in blast_record.alignments:
					ID  = alignment.accession
					if ID not in IDlist:
						IDlist.append(ID)
						#print '2: ' + ID
						IDeDict[ID] = float(alignment.hsps[0].identities)/float(alignment.hsps[0].align_length)*100.0 ##
						#print '3: ' + str(IDeDict[ID])
						
	except:
		writelog("blast error: " +  record.id)
		toBlast.append(record)
	if toBlast != []: #if there were problems blasting, try again
		x = toBlast.pop(0)
		blast(x, IDlist, IDeDict,toBlast,listsize)
	return IDlist, IDeDict,toBlast



def envSeqs(inputfile,listsize,idmin):
	toBlast = []
	IDeDict = {}
	out2 = open(inputfile + '_envseqs.fas','a')
	out3 = open(inputfile + '_sequence_names','a')
	out4 = open(inputfile + '_gb','a')
	records = SeqIO.parse(open(inputfile),format="fasta")
	IDlist = []
	for record in records:
		IDlist, IDeDict,toBlast = blast(record, IDlist, IDeDict,toBlast, listsize)
	
	
	print IDlist
	
	
	try:
		for ID in IDlist:
			print IDeDict[ID], idmin
			if 	IDeDict[ID] > idmin:
				handle2 = Entrez.efetch(db="nucleotide", id=ID, rettype="gb")
				record2 = SeqIO.read(handle2,"gb")
				name = record2.description.split()[1][0:2] + '_' + record2.description.split()[2][0:3] + "_" + record2.id
				out3.write( name + ':' + str(IDeDict[ID]) + ':' + record2.description + '\n') ##
				out2.write('>' + name + '\n')
				out2.write(str(record2.seq) + '\n')
				out4.write(str(record2))

	except:
		writelog("efetch error: " +  ID)
	
	
	out2.close()
	out3.close()
	out4.close()			
def writelog(input):
	outfile = open('log.txt','a')
	outfile.write(input + '\n')
	outfile.close

def reblast():
	
	for record in SeqIO.parse(open('RiddedSeqs.fasta','r'),'fasta'):
		try:
			result_handle2 = NCBIWWW.qblast("blastn", "nr", record.format("fasta"), ncbi_gi=False, format_type="XML", entrez_query= 'NOT ' + entrezQuery , expect = "2e-50", hitlist_size = 1) #change evalue cut off to 2e-10
			blast_records2 = NCBIXML.parse(result_handle2)
		
			for blast_record2 in blast_records2:
		
				if blast_record2.descriptions:
					for alignment in blast_record2.alignments:
						ID2  = alignment.accession
						print record.id, float(alignment.hsps[0].identities)/float(alignment.hsps[0].align_length)*100.0
						handle2 = Entrez.efetch(db="nucleotide", id=ID2, rettype="gb")
						record2 = SeqIO.read(handle2,"gb")
						print record2.annotations["taxonomy"] #and float(alignment.hsps[0].identities)/float(alignment.hsps[0].align_length) > 90: ??
						newname = record.id + "_" + record2.annotations["taxonomy"][2][0:3] + "_" + record2.annotations["taxonomy"][-1] + '_' + str(float(alignment.hsps[0].identities)/float(alignment.hsps[0].align_length)*100)[0:4]
				else:
					newname = record.id
		except:
			newname = record.id
		outfile = open('RiddedSeqs_rename.fasta','a')
		outfile.write('>' + newname + '\n' + str(record.seq) + '\n')
		outfile.close()


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
		if not re.search(',',sterm):
			stermList = []
			stermList.append(sterm)
		else:
			stermList = sterm.split(',')
		inputfile = 'seq_searchdb.fas'
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

	try:
		writelog("Sterm = " + sterm)
		writelog("Update date = " + str(update))
	except:
		a='a'
	writelog("Rid identity = " + str(RIDin))
	writelog("min identity = " + str(idmin))
	writelog("min length = " + str(lenmin))
	l = raw_input('What is the maximum return list length cut off? (hit return for default of 1000 seqs) ')
	try:
		num = float(l) + 1
	except TypeError:
		print 'Your input must be a number.  Try again. '
		main()
	except ValueError:
		l = ""	
	if l == "":
		listsize = 1000
	else:
		listsize  = int(l)
	#m = raw_input('Do you want to limit your returned environmental sequence by taxonomy of best blast hit? If yes, enter one term to limit (i.e. "Alveolata" or  "Ciliophora" or "Spirotrichea") '
	
	
	z = raw_input('Hit return when you are ready to continue. ')

	if search[0] == 'y':
		searchgb(stermList, mindate)
	
		
	
	#envSeqs(inputfile, listsize, idmin)
	#os.system('cat ' + inputfile + ' ' + inputfile + '_envseqs.fas  > ready2RidFast.fas')
	

	import rid_by_treev2
	rid_by_treev2.rid(RIDin)
	
	reblast()
	
	os.system('mafft --auto  --anysymbol RiddedSeqs_rename.fasta  > RiddedSeqs_renamemafftout')
	os.system('FastTree -nt RiddedSeqs_renamemafftout >  RiddedSeqs_renameOut.tree')

main()