#!/usr/bin/env python

'''
Extract the most probable ORF in a given sequence 
The most probable ORF is the longest open reading frame found in the sequence
When having same length, the upstream ORF is selected
modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/v1.2.1/
'''

class ExtractORF:
	def __init__(self,seq):
		self.seq = seq
		self.result=(0,0,0,0)
		self.longest = 0
	
	def codons(self,frame):
		start_coord = frame
		while start_coord + 3 <= len(self.seq):
			yield (self.seq[start_coord:start_coord+3],start_coord)
			start_coord += 3
	def longest_orf_in_seq(self,frame_number,start_codon,stop_codon):
		codon_posi = self.codons(frame_number)
		start_codons = start_codon
		stop_codons = stop_codon
		while True:
			try:
				codon,index = codon_posi.next()
			except StopIteration:
				break
			if codon in start_codons and codon not in stop_codons:
				ORF_start = index
				end = False
				while True:
					try:
						codon,index = codon_posi.next()
					except StopIteration:
						end = True
						integrity = -1
					if codon in stop_codons:
						integrity = 1
						end = True
					if end:
						ORF_end = index+3
						ORF_Length=(ORF_end-ORF_start)
						if ORF_Length > self.longest:
							self.longest = ORF_Length
							self.result = (integrity,ORF_start,ORF_end,ORF_Length)
						if ORF_Length == self.longest and ORF_start < self.result[1]:
							self.result = (integrity,ORF_start,ORF_end,ORF_Length)
						break
	def longest_ORF(self,start=['ATG'],stop=['TAA','TAG','TGA']):
		orf_seq = ""
		for frame in range(3):
			self.longest_orf_in_seq(frame,start,stop)
		orf_seq = self.seq[self.result[1]:self.result[2]]
		ORF_integrity = self.result[0]
		ORF_length = self.result[3]
		return ORF_length,ORF_integrity,orf_seq

