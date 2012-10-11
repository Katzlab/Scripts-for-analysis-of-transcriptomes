"""
removeSSU.py
	usage: python removeSSU.py FILE-TO_ANALYZE
	FILE-TO_ANALYZE is a fasta file of ESTs
nec. files:  Database,  to blast against to remove SSUs
path /Users/katzlab/Files/ssulsu_website.fa (or enter your own)

output: OutSSUSeqs: SSU and LSU sequences found in the input file
		OutnoSSUSeqs.fas: sequences that are not SSU or LSU

"""



import os
import sys




def ridSSU(arg,input,path): #uses blastn to remove SSU sequences.  Output 2 fasta files, one all SSU, one no SSU
	placeholder = open(arg + '_out','a')
	OutSSU = arg + "OutSSUSeqs.fas"
	outSSU = open(OutSSU,'a')
	OutnoSSU = arg + "OutnoSSUSeqs.fas"
	outnoSSU = open(OutnoSSU,'a')
	toRID = []
	command = 'blastn  -evalue .00001 -num_descriptions 1 -outfmt 10 -db ' + path + ' -num_threads 2 -num_alignments 0  -query ' + arg + ' -out ' + arg + '_out'
	os.system(command)

	out = arg + '_out'
	file_in = open(out, 'r')
	f=file_in.readlines()
	file_in.close
	
	for line in f:
		contigName = line.split(',')[0]	
		toRID.append(contigName)
	
	
	estIn = open(arg, 'r')
	estFile = estIn.readlines()
	estIn.close

	for line in estFile:
		if (line[0] == '>'):
			seqName = line.split()[0][1:]
		if (seqName in toRID):
			outSSU.write(line)
		else:
			outnoSSU.write(line)
	if input == 'y':
		addSSU(arg)



def addSSU(arg):	
	#try:
	os.system('cp SSU/* SSU/SSU.fas')
	os.system('formatdb -i SSU/SSU.fas -p F -o T') 
	command = ('blastn  -evalue .00001  -outfmt 10 -db SSU/SSU.fas -num_threads 2 -num_alignments 0  -query '  + arg + 'OutnoSSUSeqs.fas -out ' + arg + '_outlist')
	os.system(command)
	name = arg + '_outlist'
	infile2=open(name,'r')
	inFile2=[]
	for line in infile2:
		inFile2.append(line.split(',')[0])
		
	
	
	outfile=('outlist','a')
	
	inname = arg + 'OutnoSSUSeqs.fas'
	infile = open(inname,'r')
	inFile = infile.readlines()
	infile.close()
	outFile = open(inname,'w')
	outFile2 = open('moreSSUs','a')
	for line in inFile:
		if (line[0] == '>'):
			flag = 0
			
			if (line.split()[0][1:] in inFile2):
				print line.split()[0][1:]
				flag = 1
				outFile2.write(line)
		if (flag == 0):
			outFile.write(line)
		if (flag == 1):
			outFile2.write(line) 	
	outFile.close()	
	os.system('rm SSU/SSU*')


def main():
	for arg in sys.argv[1:]:
		
		path = raw_input('Give the full path to the ssu/lsu data file you are blasting against. ')
		if path == '':
			path = '/Users/katzlab/Files/ssulsu_website.fa' #default for JG
		ans = raw_input('does this database need to be formatted for blast? (y/n) ')
		if ans[0] == 'y':
			os.system('formatdb -i ' + path + ' -p F -o T')
		
		input = raw_input("Would you like to add a taxon-specific SSU to the SSU removal pipeline? (y/n) ")
		if str.lower(input[0]) == 'y':
			x = raw_input('Put your ssu as a fasta in a folder called "SSU" in this directory.  Then hit return (there should only be one file in the SSU folder) ')
		ridSSU(arg,str.lower(input[0]),path)
		
		
main()
