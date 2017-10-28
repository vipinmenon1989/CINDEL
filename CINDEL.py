#### A Vipin Menon,BIG LAB,May 15th 2016
####////////////////////////////////////////////////////////////////////

# single code for calcualting score for CRISPR-Cpf1 sgRNA,CINDEL score 
# The code  comprises of Free Energy ,Mononucleotide ,Dinucleotide and Independent nucleotide composition
# This code generates two results one for whole batch and other one for single sequence,option is given 
# The total length of the sequence should be  27 bp , whereas the target sequence is 23 bp and 4bp would be PAM , Please comply with PAM : TTTV (V = A,G or C) 
import sys,os 
import math,RNA				 
def calculatescore(seq):
        					
        parameters = [(9,'A',0.037608778),(13,'AA',0.167340348),(25,'AA',0.364227683),(6,'AA',-0.37092459),(6,'AC',0.083159983),
(9,'AC',0.459109261),(16,'AG',-0.257994467),(22,'AT',0.469569195),(26,'C',0.166195809),(14,'CA',0.11154741),(5,'CA',-0.317283912),(23 ,'CC',0.130179283),(4,'CC',1.065828073),(24,'CG',-0.186173348),(6,'CG',-0.231597251),(7,'CG',0.275263279),(13,'CT',0.025086518),(19,'CT',0.164219823),(24,'CT',0.083624311),(6,'CT',0.120799685),(7,'CT',-0.336424167),(9,'CT',-0.094751841),(19,'G',-0.076228205),(22,'G',-0.066751042),(23,'G',-0.102884761),(26,'G',-0.043253074),(4,'G',0.489354734),(10,'GA',0.059114363),(18,'GA',0.077528126),(21,'GA',0.216211044),(4,'GA',0.077064537),(20,'GC',-0.123286737),(21,'GG',-0.124324225),(22,'GG',-0.358412408),(5,'GG',-0.058789262),(8,'GG',0.126113307),(23,'GT',-0.15567756),(21,'T',0.048550229),(4,'T',-0.781196013),(9,'T',-0.028800337),(23,'TA',0.444556704),(24,'TA',0.274164613),(7,'TA',0.391209275),(15,'TC',0.123925716),(25,'TC',0.194045073),(10,'TG',-0.190097449),(16,'TG',0.036768788),(6,'TG',0.12530035),(10,'TT',-0.173588498),(12,'TT',-0.035141319),(16,'TT',-0.127599558)] # Coefficients for dinucleotides and single nucleotides based on position

        intercept = -1.830804463
        Free_Energy = 0.11018728 
	A = 0.101900819
	AC = 0.038389043
	CG = 0.006835173
	CC = -0.24886262
	TA = 0.037868724
        score  = intercept
        #### feature extraction 
	##Free_energy
        Energy = seq[4:27]
        Energycal = RNA.fold(Energy)[-1]
        Energycal = round(Energycal,0)
        score +=(Energycal*Free_Energy)
	# independent nucleotide composition
	score  = score + A*(seq.count('A'))
	score = score + AC*(seq.count('AC'))
	score = score + CG*(seq.count('CG'))
	score = score + CC*(seq.count('CC'))
	score = score + TA*(seq.count('TA'))
	
        for position, bp, wt in parameters:
                subSeq = seq[position:position+len(bp)]
                if subSeq==bp:
                        score += wt
        return (1.0/(1.0 + math.exp(-score)))




### Batch processing for Cpf1 based ,Please provide gRNA ID and the sequence with  header  in csv format #########
def bashmethod (input):

	file = open('Score.csv','wb')
	sequence = open(input,'r')
	file.write('Score' + '\n')
	sequence = sequence.readlines()
	del sequence [0]
	set = {}
	#fal = []
	p = {}
	for line in sequence:
        	seq = line.strip().split(',')
        	idseq = seq[0]
		ids  = seq[1]
        	set[idseq] = ids
	#	fal.append(clas)
	for j in set:
		if set[j][0:3] == 'TTTA' or 'TTTG' or 'TTTC':
	       		n = calculatescore(set[j])
			p[j] = set[j],n
		else:
			print "please provide only 'TTTV' where V = 'A' ,'G' or 'C'"
	for x ,y in p.items():
		file.write(str(x) +  ',' + str(y[0]) + ',' + str(y[1]) + '\n')
	file.close()
	print file

def guidefinder(input):
	sequence = input
	ds = {}
	if (len(sequence) > 100):
		for n in xrange(len(sequence)): 
			st = sequence.find('TTT',n)
			if st == n:
				if int(0) <= int(st) < len(sequence)-int(27):
					nstart = int(st) - int(0)
					nend = int(st) + int(27)
					sequences = sequence[nstart:nend]
					score = calculatescore(sequences)
				#	pstream = sequences[0:17]
					PAM = sequences[0:4]
					target = sequences[4:27]
				#	target = sequence[st:]
					ds[sequences] = PAM,target,score
		print ds



### single sequence score ######
def sequencescore(input):
	sequence = input
	result = calculatescore(sequence)
	score = 'sgRNA score:' + str(result)
	print  score


func_arg = {"-a": bashmethod, "-b": sequencescore, "-c": guidefinder}


if __name__ == "__main__":
    func_arg[sys.argv[1]](sys.argv[2])

