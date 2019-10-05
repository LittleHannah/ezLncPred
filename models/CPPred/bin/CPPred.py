import os,sys
from feamodule import CTD
from feamodule import ProtParam as PP
from feamodule import ORF_length as len
import Bio.SeqIO as Seq
from feamodule import fickett
from feamodule  import FrameKmer
import argparse as agp

def coding_nocoding_potential(input_file):
	coding={}
	noncoding={}
	for line in open(input_file).readlines():
		fields = line.split()
		if fields[0] == 'hexamer':continue
		coding[fields[0]] = float(fields[1])
		noncoding[fields[0]] =  float(fields[2])
	return coding,noncoding

def output_feature(seq_file,hex_file,species):
	tmp = open('test.f_svm','w')
	feature = open('test.feature','w')
	out_label = 1
	coding,noncoding = coding_nocoding_potential(hex_file)
	if species == "Human":
		feature.write("\t".join(map(str,["#ID","ORF-integrity","ORF-coverage","Instability","T2","C0","PI","ORF-length","AC","T0","G0","C2","A4","G2","TG","A0","TC","G1","C3","T3","A1","GC","T1","G4","C1","G3","A3","Gravy","Hexamer","C4","AG","Fickett","A2","T4","C","G","A","T"]))+"\n")
	if species == "Integrated":
		feature.write("\t".join(map(str,["#ID","ORF-coverage","ORF-integrity","GC","Instability","ORF-length","T0","Fickett","G2","C3","PI","A3","C1","G3","Hexamer","TG","G1","TC","A0","A1","AC","C2","G0","T4","C0","A4","G","A2","T","T3","G4","C4","Grary","T2","AG","AT","T1","A","C"]))+"\n")
	for seq in Seq.parse(seq_file,'fasta'):
		seqid = seq.id
		A,T,G,C,AT,AG,AC,TG,TC,GC,A0,A1,A2,A3,A4,T0,T1,T2,T3,T4,G0,G1,G2,G3,G4,C0,C1,C2,C3,C4 = CTD.CTD(seq.seq)
		insta_fe,PI_fe,gra_fe = PP.param(seq.seq)
		fickett_fe = fickett.fickett_value(seq.seq)
		hexamer = FrameKmer.kmer_ratio(seq.seq,6,3,coding,noncoding)
		Len,Cov,inte_fe = len.len_cov(seq.seq)
		if species == "Human":
			tem = [inte_fe,Cov,insta_fe,T2,C0,PI_fe,Len,AC,T0,G0,C2,A4,G2,TG,A0,TC,G1,C3,T3,A1,GC,T1,G4,C1,G3,A3,gra_fe,hexamer,C4,AG,fickett_fe,A2,T4,C,G,A,T]
			feature.write("\t".join(map(str,[seqid,inte_fe,Cov,insta_fe,T2,C0,PI_fe,Len,AC,T0,G0,C2,A4,G2,TG,A0,TC,G1,C3,T3,A1,GC,T1,G4,C1,G3,A3,gra_fe,hexamer,C4,AG,fickett_fe,A2,T4,C,G,A,T]))+"\n")
		if species == "Integrated":
			tem = [Cov,inte_fe,GC,insta_fe,Len,T0,fickett_fe,G2,C3,PI_fe,A3,C1,G3,hexamer,TG,G1,TC,A0,A1,AC,C2,G0,T4,C0,A4,G,A2,T,T3,G4,C4,gra_fe,T2,AG,AT,T1,A,C]
			feature.write("\t".join(map(str,[seqid,Cov,inte_fe,GC,insta_fe,Len,T0,fickett_fe,G2,C3,PI_fe,A3,C1,G3,hexamer,TG,G1,TC,A0,A1,AC,C2,G0,T4,C0,A4,G,A2,T,T3,G4,C4,gra_fe,T2,AG,AT,T1,A,C]))+"\n")
		print >> tmp, out_label,
		for label,item in enumerate(tem):
			print >> tmp, str(label+1)+':'+str(item),
		print >> tmp
	tmp.close()

def predict(range_file,model_file):
	os.system('models/CPPred/libsvm-3.22/svm-scale -r '+ range_file + ' test.f_svm  > test.scaled ')
	os.system('models/CPPred/libsvm-3.22/svm-predict -b 1 test.scaled ' + model_file +' tmp.txt >  tmp2.txt')

	coding_poten = open('coding_potential','w')
	coding_poten.write("\t".join(map(str,["table","coding_potential"]))+"\n")

	for line in open('tmp.txt').readlines():
		if line[0]=="l":
			continue
		coding_potential=line.split(" ")[1]
		if line.split(" ")[0] == "1":
			coding_poten.write("\t".join(map(str,["coding",coding_potential]))+"\n")
		else:
			coding_poten.write("\t".join(map(str,["noncoding",coding_potential]))+"\n")
	

def merge(output_file):
	os.system("paste test.feature coding_potential >"  + output_file)


def deleted():
	os.system("rm test.*")
	os.system("rm coding_potential")
	os.system("rm tmp*")


def main(args):
	#parser = agp.ArgumentParser()
	#parser.add_argument('-i','--RNA_file',help="the input FASTA file of RNA sequence")
	#parser.add_argument('-hex','--hexmar',help="the input of hexmar")
	#parser.add_argument('-r','--range',help="the input file of range")
	#parser.add_argument('-mol','--model',help="the input file of model")
	#parser.add_argument('-spe','--species',help="the input species")
	#parser.add_argument('-o','--outfile',help="output file")
	#args = parser.parse_args()

	args.hexmar = 'models/CPPred/Hexamer/'+args.species+'_Hexamer.tsv'
	args.range = 'models/CPPred/'+args.species+'_Model/'+args.species+'.range'
	args.whichmodel = 'models/CPPred/'+args.species+'_Model/'+args.species+'.model'
	output_feature(args.fasta,args.hexmar,args.species)
	predict(args.range,args.whichmodel)
	merge(args.outfile)
	deleted()


if __name__ == '__main__':
	main()
