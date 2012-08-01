'''
This script takes a fasta downloaded from genbank.  Sequences will be named like this:
'>gi|389607217|dbj|AB725339.1| Oxytrichidae sp. MP29 gene for 18S ribosomal RNA, partial sequence_389607217'

'''


try:
	from Bio import SeqIO
	from Bio import Entrez
except:
	print 'Biopython is needed for these scripts. Please make sure it is installed. '

import sys
import re
email = raw_input('An email address is necessary to use the ncbi entrez server.  Please input your email address. ')
while  not re.search('[.+@.+]',email):
	email = raw_input('Please input a valid email address. ')

Entrez.email = email


GspeList = {}
print '#############################################################################'
print 'Useage is "python renameGBfasta.py <yourFastaFile>'
print '#############################################################################'
print 'output is renamed file and a file containing the code key'
print '#############################################################################'
print 'This script takes a fasta file downloaded from genbank and renames the sequences with a standard naming format based on the NCBI taxonomy for that taxon.\n'
print 'The standard is to rename each sequences with a 4 letter species code plus the gi number unique to that sequence.\n'
print 'If this code is not unique, it adds a letter to the end of the name. If there are too many non-unique codes, this script will fail .\n'
print 'You may want to preceed this code with more information, i.e. the major clade of the sequence or another part of the ncbi taxonomy.\n'
print '#############################################################################'

	


def rename(arg, code, prec):
	if prec == '_':
		prec = ''
	additional = ''
	alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
	out = open(arg + '_renamed.fas','a')
	GspeList = []
	infile = open(arg,'r')
	for line in infile:
		if line[0] == '>':
			ginum = line.split('|')[1]
			taxon, taxonomy = getTax(ginum)
			#print taxon, taxonomy
			try:	
				Genus = taxon.split()[0]
				species = taxon.split()[1]
				if code == '':
					additional = ''
				else:
					additional = taxonomy[int(code)][0:2] + '_'

				Gspetemp = Genus[0] + species[0:3]
			except:
				errorOut = open('errorlog','a')
				errorOut.write('no taxon for ' + ginum + '\n')
				Gspetemp = 'none'			
			if checkGspe(Gspetemp, taxon) == True:
			 	Gspe = Gspetemp
			else:
				for i in alphabet:
                                        Gspetemp = Genus[0] + species[0:2] + i
                                        if checkGspe(Gspetemp, taxon) == True:
                                                Gspe = Gspetemp
                                                break
                                        elif checkTaxon(taxon) != False:
                                                Gspe = checkTaxon(taxon)
                                else:
                                        running = True
                                        while running:
                                                Gspetemp = raw_input('This genus and species code ' + Gspetemp + ' has been used.  Enter a code for ' + taxon)
                                                if checkGspe(Gspetemp, taxon) == True:
                                                        Gspe = Gspetemp
                                                        running = False
                                                else:
                                                        running = True
			#print Gspe

			out.write('>' + prec + additional + Gspe + '_' + ginum + '\n')
		else:
			out.write(line)




def getTax(ID):
	
	try:
		handle = Entrez.efetch(db="nucleotide", id=ID, rettype="gb")
		record = SeqIO.read(handle,"gb")

		return record.annotations["organism"], record.annotations["taxonomy"]
		#record_iterator = SeqIO.parse("ls_orchid.gbk", "genbank")
		#first_record = record_iterator.next()
		#print first_record
	except:
		out2=open('errorlog','a')
		out2.write(ID + '\n')
		out2.close()


def checkGspe(Gspetemp, taxon):
	
	
	if Gspetemp in GspeList.keys():
		if GspeList[Gspetemp] == taxon:
			return True
		else:
			return False
			
	else:
		outGspe = open('Genus_sp_codeList','w')
		GspeList[Gspetemp] = taxon
		outGspe.write(str(GspeList) + '\n')
		outGspe.close()
		return True

def get_key(taxon):
	return [key for key, value in GspeList.iteritems() if value == taxon][0]
	
def checkTaxon(taxon):
	
	if taxon in GspeList.values():
		Gspetemp = get_key(taxon)
		return Gspetemp
	else:
		return False

def getCode():
	print "the ncbi taxonomy looks something like this: ['Eukaryota', 'Alveolata', 'Ciliophora', 'Intramacronucleata', 'Spirotrichea', 'Choreotrichia', 'Tintinnida', 'Metacylididae', 'Climacocylis']"
	print '#############################################################################'

	b = raw_input('The standard is to rename each sequence with a 4 letter code made up from the genus and species names of the taxon. Would you like to add the first two letters of another taxonomic level? ')
	try:
		assert b[0] == 'y' or a[0] == 'n'
	except:
		print 'please enter y or n'		
		getCode()
	if b[0] == 'y':
		print '#############################################################################'
		rank = raw_input('Enter F for family, O for order, C for class,	or M for the second word in the taxonomy (i.e. Alveolata in the example above). If there is no entry for that level for a taxon, it will add "No" ')
		try:
			assert rank == 'F' or rank == 'O'  or rank == 'C' or rank == 'M'
		except:
			print 'please enter F, O, C or M'		
			getCode()		
		if rank == 'F':
			code = '-2'
		if rank == 'O':
			code = '-3'
		if rank == 'C':
			code = '-4'
		if rank == 'M':
			code = '2'
	else:
		code = ''		
	print '#############################################################################'

	c = raw_input('Would you like to add something standard to the beginning of every sequence? ')
	try:
		assert c[0] == 'y' or c[0] == 'n'
	except:
		print 'please enter y or n'		
		getCode()
	if c[0] == 'y':
		print '#############################################################################'

		print ('Remember, if your name is too long you might run into trouble later, as some software truncates the sequence names. ')
		prec = raw_input('Enter the code you would like added to the beginning of every sequence name. ')
		prec = prec + '_'
	else:
		prec = ''	
	return code , prec


def main():
	

	for arg in sys.argv[1:]:

		a = raw_input('Do you need to add to the code? ')
		try:
			assert a[0] == 'y' or a[0] == 'n'
		except:
			print 'please enter y or n'		
			main()
				
		if a[0] == 'y':
			code, prec = getCode()
		else:
			code = ''
			prec = ''
		
		rename(arg, code, prec)

main()
