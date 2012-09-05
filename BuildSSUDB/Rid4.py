"""
#
#RID4.py
#
#This sums up all scripts -- from large targeted protein list to a file that is ready to run RID on.  
#Input:  large file to blast against i.e. "parsedAllProt"
#		small file of genes to blast, i.e. "representative_actins.fas"
#Output:	file of all fastas that blast to that gene, i.e. "allActinsfas" that is ready to use in RID
#
#**************************************
#
#blast(): 
#	runs formatdb on parsedAllProt and blasts file your_file against parsedAllProt, output: blast_out in tab form
#	Calls  rmDups() to remove duplicate blast hits (since blasting from a file of representitives of 
#	all 6 supergroups, many duplicate hits in blast_out)
#parser():
#	reads in blast_out, pulls out the GBnum and the %identity.  Writes GBnum to a file (gbNumList) and passes it to 
#	blast2fasta(GBnum).
#blast2fasta(GBnum): 
#	gets the fasta file from the genbank number and writes each to a file.
#rmDups():
#	opens blast__out and rewrites file without duplicated lines.
#getFasta(): does the work of parser() but gets the fasta from the big file instead of from genbank.  May be faster.
#"""
#
#
#from Bio.Blast import NCBIStandalone
#from Bio import GenBank
import os
#import re
#import string
#import sys
#
def rid(RIDin,idmin, lenmin):	

	#x=[]
	outName = 'Output' + str(RIDin)
	ridformatdb()
	ridblast()
	
	ridorthology('ready2RidFast.fas',RIDin, outName)
	removeDups(idmin, lenmin)
#
#	
def ridformatdb(): #formats the fasta file for blasting  CHANGE YOUR-FILE-NAME
	os.system('formatdb -i ready2RidFast.fas -p F -o T')


def ridblast():  #blasts fasta file against itself  CHANGE YOUR-FILE-NAME
	os.system('blastn -db ready2RidFast.fas -query ready2RidFast.fas -num_descriptions 1000 -num_alignments 1 -evalue 1e-20  -num_threads 2 -outfmt 7 -out ridblast_out')
	


def ridorthology(arg,RIDin, outName): #extracts gb numbers and %identity from blast report - calls listManage and blast2fasta
	#Passes the numbers to listManage(GBnum, identity) gets back a list of numbers that are not orthologs
	#Passes these numbers to blast2fasta(GBnum) which writes a fasta file
	file_in = open('ridblast_out', 'r')
	#f=file_in.readlines()
	#file_in.close
	toRID = []
	toKeep = []
	ridDict = {}
	qGBnum = ""
	gbnumlist = []
	fileName = ('tmp_' + qGBnum)
	flag = 0
	out1 = open('WhoRiddedWho','a')
	for line in file_in:
	
		if line[0] != '#':
			qGBnum = line.split()[0]
				
	
			if toRID.count(qGBnum) == 0:
				if toKeep.count(qGBnum) == 0:
					toKeep.append(qGBnum)
				lineSplit = line.split()
				subject = lineSplit[1]
				sGBnum = lineSplit[1]
				identity = lineSplit[2]
				print qGBnum, sGBnum, identity	
				if qGBnum.strip() != sGBnum.strip():
					#print float(identity), float(RIDin)
					if float(identity) >= float(RIDin):
				#x = listManage2(qGBnum, sGBnum, identity, RIDin) #x == None or sGBnum
		
				#if x != None: #x is an subject to be ridded

						if toRID.count(sGBnum) == 0:
							toRID.append(sGBnum) #a list of things to be ridded.
							try:
								ridDict[qGBnum].append(sGBnum)
							except:
								ridDict[qGBnum] = []
								ridDict[qGBnum].append(sGBnum)
	for query in ridDict.keys():
		out1.write(query + ':' + str(ridDict[query]) + '\n')			
	print 'keeping ' + str(len(toKeep))
	print 'ridding ' + str(len(toRID))
	print 'getFastaFast(' + arg + outName + ')'
	getFastaFast(arg,toKeep, outName)

	
#	
#def listManage2(qGBnum, sGBnum, identity, RIDin): #REMOVED__returns sGBnum if it should be ridded, else returns None
#	fas_in = open('ready2RidFast.fas','r')
#	a = fas_in.readlines()
#	
#	if float(identity) >= float(RIDin):
#		
#		#get staxname and qtaxname
#		for line in a:
#			if line[0] == '>':
#				print line
#				if qGBnum in line.split()[0].split('|')[3]:
#					qtax = line.split('[')[1][:-1]
#				if sGBnum in line.split()[0].split('|')[3]:
#					stax = line.split('[')[1][:-1]	
#
#		if stax == qtax:
#			
#			return sGBnum
#		else:
#			return	""
#			
#	else:
#		return ""
#
#
#
#
def getFastaFast(arg,m, outName): # Gets fastas from file instead of Genbank
	flag = 0

	fasfile_in=open(arg, 'r')
	r=fasfile_in.readlines()
	fasfile_in.close()
	#print m
	
	for line in r:
		if line[0] == '>':
			#print line.split()[0].strip()[1:]	
			#print line.split('>')[1].strip()	
			if (line.split()[0].strip()[1:]) in m:
				flag = 1
		
					
			else:
				flag = 0
					
									
		if flag == 1:
			file_out= open(outName,'a')
			file_out.write(line)
			file_out.close
					




#def blast(arg, idmin, lenmin):
#	# blasts new sequences against the database and outputs a report in tabular form. Calls
#	#rmDups to remove duplicate blast hits (since blasting from a file of representitives of 
#	#all 6 supergroups, many duplicate hits in blast_out)
#	os.system('formatdb -i '+ arg + ' -p F -o T')
#	os.system('perl /usr/local/ncbi/blast/bin/legacy_blast.pl blastall -p blastn -m 8 -b 100000 -d '+ arg +  ' -v 100000 < '+ arg + ' > blast_out --path /usr/local/ncbi/blast/bin')
#	removeDups(idmin, lenmin) #Get rid of duplicate blast hits
#
#
def removeDups(idmin, lenmin): #cleans up blast_out so it only has one of each hit, 
	#also removes identiry < 60, length < 750
	x = []
	flag = 0
	blast_out = open('ridblast_out', 'r')
	b = blast_out.readlines()
	blast_out.close
	os.system('cp ridblast_out ridblast_out1')
	file_out = open('ridblast_out', 'w')		
	file_out.close
	
	for line in b:
		#print line
		if line[0] != '#' and len(line.split()) > 4:
			z=line.split()
			identity = z[2]
			length = z[3]
			GBnum=z[1]
			if x.count(GBnum) == 0:	

				if float(identity) >= idmin:
					
					if float(length) >= lenmin:					
						flag = 1
						x.append(GBnum)
					else:
						flag = 0
				else:
					flag = 0	
			else:
				flag = 0

			
			if flag == 1:	
				file_out = open('ridblast_out', 'a')		
				file_out.write(line)
				file_out.close

#def parser(arg,RIDin): # Gets fasta files from genbank
#	#reads in the file from the blast and pulls out the genbank numbers
#	m = []
#	
#	blast_out = open('blast_out', 'r')
#	
#	for line in blast_out.readlines():
#		x=line.split()
#		identity=x[2]
#		y=x[1].split('|')
#		GBnum=y[3]
#		print GBnum
#		if float(identity) <= float(RIDin):
#			print identity 
#			print GBnum
#			m.append(GBnum)
#	
#	outName=arg+'_'+ str(RIDin) + '_Ridded.fas'	
#	
#	getFastaFast(arg,m,outName) # pass the list of GBnumbers to the next function 
#					
#