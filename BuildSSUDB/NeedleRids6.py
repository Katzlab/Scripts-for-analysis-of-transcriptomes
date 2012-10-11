def test(f,ridid):
	outfinal = open(f + 'needlescores.txt','a')
	outscore = open(f + 'pairwise_out_scores.csv','a')
	from Bio import AlignIO
	import re
	flag = 1
	infile = open(f + 'needle.txt','r')
	midline = ""
	gap = 0.0
	ident = 3.0
	sim = 3.0
	diff = 0.0
	
	for line in infile:
		outfinal.write(line)
		if re.search('# 1: ', line):
			seq1 = str(line.split(':')[1].strip())
		if re.search('# 2: ', line):
			seq2 = str(line.split(':')[1].strip())		
		if re.search('# Length:',line):
			tlengthlist = re.findall('\d+', line)
			tlength = float(tlengthlist[0])
			
		#print line[0]
		if line.strip() == '#=======================================':
			flag = 0
		
		if line.strip() == '#---------------------------------------':
			flag = 1
			
		if flag == 0:
			if line[0] == ' ':
				midline=midline + line.strip()
	outfinal.write('\n')	
	midline = midline + '*'
	#print midline
	midline2 = ""
			
	id = 0 
	gp = 0
	
	
	
	
	for char in midline:
		if id < 4:	
			if char == '|':
				id = id + 1
			else: 
				id = 0
		else:
			midline2 = midline2 + char
	
	for char in midline2:
	
		if gp < 4:	
			if char == ' ':
				gp = gp + 1
				gap = gap + 1
			else:
				gp = 0
			if char == '|':
				ident = ident + 1
				sim = sim + 1
			if char == ':':
				sim = sim + 1
			if char == '.':
				diff = diff + 1
			if char == '*':
				gp = 10
	len = gap + sim + diff
	
	infile.close()
	outscore.write(seq1 + ',' + seq2 + ', len: ' + str(len) + ', tlength: ' + str(tlength) + ', identity: ' + str(ident) + ', similarity: ' + str(sim) + ', len/tlength: ' + str(len/tlength) + ', sim/len: ' + str(sim/len) + '\n')
	print float(sim)/float(len), ridid
	if float(sim)/float(len) * 100 < float(ridid):# and sim/ident > 1.02: #len/tlength > 0.67 and sim/len > 0.75 and 
		return True
	else:
		return False