import sys
from string import maketrans
import ORF

def extract_feature_from_seq(seq,stt,stp):
	'''extract features of sequence from fasta entry'''
	
	stt_coden = stt.strip().split(',')
	stp_coden = stp.strip().split(',')
	transtab = maketrans("ACGTNX","TGCANX")
	mRNA_seq = seq.upper()
	mRNA_size = len(seq)
	tmp = ORF.ExtractORF(mRNA_seq)
	(CDS_size1, CDS_integrity, CDS_seq1) = tmp.longest_ORF(start=stt_coden, stop=stp_coden)
	return (mRNA_size, CDS_size1,CDS_integrity)

start_codons = 'ATG'
stop_codons = 'TAG,TAA,TGA'
Coverage = 0

def len_cov(seq):
	(mRNA_size, CDS_size,CDS_integrity) = extract_feature_from_seq(seq = seq, stt = start_codons,stp = stop_codons)
	mRNA_len = mRNA_size
	CDS_len = CDS_size
	Coverage = float(CDS_len)/mRNA_len
	Integrity = CDS_integrity
	return(CDS_len,Coverage,Integrity)
