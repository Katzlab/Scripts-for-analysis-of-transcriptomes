import re,os
from Bio import SeqIO
from Bio import Phylo
toKeep = []
toRemove = []
seqDict = {}
arg = 'ready2RidFast.fas' #fastafile name
infile=SeqIO.parse(arg,'fasta')

def rid(RIDin):
	for seq in infile:
		try:
			name = seq.id
			
		except:
			name = seq.id	
		seqDict[name] = seq.seq
	


#build a mafft guide tree.
		
	os.system('mafft --auto  --anysymbol ' +   arg + '  > mafftout')
	os.system('FastTree -nt mafftout >  tree')
	treefile ='tree'
	tree = Phylo.read(treefile,'newick')
	checkClade(tree, RIDin)
	#print toKeep, toRemove
	#for seqname in toKeep:
#		outfile = open('riddedSeqs.fas','a')
#		outfile.write('>' + seqname.strip('_') + '\n' + str(seqDict[seqname.strip('_')]) + '\n')
#		outfile.close()
#take a leaf, compare it to its sister
def get_parent(tree, child_clade):
		#############################################################
		#############################################################
	node_path = tree.get_path(child_clade)
	try:
		return node_path[-2]				
	except:
		return 'None' 	
def getName(cladename):
	newname = ""
	cladelist = cladename.split('_')
	for x in cladelist:
		if x != cladelist[0] and x != cladelist[-1]:
			newname = newname + str(x) + '_'
	name = newname.strip('_')
	return name
def getName2(cladename):
	newname = ""
	cladelist = cladename.split('_')
	for x in cladelist:
		if x != cladelist[0]:# and x != cladelist[-1]:
			newname = newname + str(x) + '_'
	name = newname.strip('_')
	return name
	
	
def checkNeedle(tree, parent, child, keepone, seen, RIDin):
	answer = []
	f = tree.name
	for leaf in parent.get_terminals():
		#newname = 
		out = open('seq1.fasta','w')
		try:
			out.write('>' + leaf.name  + '\n' + str(seqDict[leaf.name]) + '\n')
		except:
			out.write('>' + leaf.name  + '\n' + str(seqDict[leaf.name]) + '\n')
                
		out.close()

	for j in range(len(parent.get_terminals())):
		
		if parent.get_terminals()[j] != leaf:
			out2 = open('seq2.fasta','w')			
			try:
				out2.write('>' + parent.get_terminals()[j].name  + '\n' + str(seqDict[parent.get_terminals()[j].name]) + '\n')
			except:
				out2.write('>' + parent.get_terminals()[j].name  + '\n' + str(seqDict[parent.get_terminals()[j].name]) + '\n')
				
			out2.close()
#	
			cline = '/usr/local/bin/needle -outfile=' + f + 'needle.txt -asequence=seq1.fasta -bsequence=seq2.fasta -gapopen=10 -gapextend=0.5' #NeedleCommandline(asequence="../Temp/seq1.fasta", bsequence="../Temp/seq2.fasta", gapopen=10, gapextend=0.5, outfile=('../Temp/' + f + "needle.txt"))		
			os.system(cline)
		
			os.system('cat ' + f + 'needle.txt >> allneedleout.txt' )
		
			import NeedleRids6
			x = NeedleRids6.test(f, RIDin)			
			print 'x = ' + str(x)
			answer.append(x)
			
			if x == True:
				break
					
	#if all(x in cladeList for x in TestList)	
	if all(ans == False for ans in answer): #all are within cut off, o up a node and try again
		newparent = get_parent(tree, parent)
		if newparent != 'None':
			checkNeedle(tree, newparent, parent, keepone, seen, RIDin)
		else:
			keepone.append(child) #the parent has a node outside the cut off, but the child does not.  Keep one from child, 'rid' the rest.
			for leaf in child.get_terminals():
				seen.append(leaf)
	else:
		keepone.append(child) #the parent has a node outside the cut off, but the child does not.  Keep one from child, 'rid' the rest.
		for leaf in child.get_terminals():
			seen.append(leaf)
	return seen
		
			
def checkClade(tree, RIDin):
	keepone = []
	seen = []
	written = []
	for clade in tree.get_terminals():
		print clade.name
		if clade not in seen:
			#seen.append(clade)
			parent = get_parent(tree,clade)

			if parent == 'None':
				seqlist =  clade.name.split('_')
				seq1name = ""	
				for i in range(len(seqlist)-1):
					if i != 0:
						seq1name = seq1name + '_' + seqlist[i]
				toKeep.append(seq1name)	
				
			else:			

				seen = checkNeedle(tree, parent, clade, keepone, seen,RIDin)
				
	for clade in keepone:
		#print clade
		out = open('RiddedSeqs.fasta','a')
		out2 = open('WhoRiddedWho','a')
		out2.write('\n' + str(clade.get_terminals()[0].name) + ':\n')
		if clade.get_terminals()[0].name not in written:
			written.append(clade.get_terminals()[0].name)
			out.write('>' + str(clade.get_terminals()[0].name)  + '\n' + str(seqDict[clade.get_terminals()[0].name]) + '\n')
		out.close()
		
		
		for seq in clade.get_terminals():
			if seq != clade.get_terminals()[0]:
				out2.write('>' + str(seq) + ',')
		out2.close()
		
		
