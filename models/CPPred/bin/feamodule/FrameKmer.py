#!/usr/bin/env python
'''the python script is downloaded from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/'''
'''deal with Kmer. DNA sequence should only A, C, G, T. python2.7 or newer'''

import os,sys
import math

def word_generator(seq,word_size,step_size,frame=0):
	'''generate DNA word from sequence using word_size and step_size. Frame is 0, 1 or2'''
	for i in xrange(frame,len(seq),step_size):
		word =  seq[i:i+word_size]
		if len(word) == word_size:
			yield word

def kmer_ratio(seq,word_size,step_size,coding,noncoding):
	if len(seq) < word_size:
		return 0
		
	sum_of_log_ratio_0 = 0.0
	sum_of_log_ratio_1 = 0.0
	sum_of_log_ratio_2 = 0.0	
	frame0_count=0.0
	frame1_count=0.0
	frame2_count=0.0
	for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=0):	
		if (not coding.has_key(k)) or (not noncoding.has_key(k)):
			continue
		if coding[k]>0 and noncoding[k] >0:
			sum_of_log_ratio_0  +=  math.log( coding[k] / noncoding[k])
		elif coding[k]>0 and noncoding[k] == 0:
			sum_of_log_ratio_0 += 1
		elif coding[k] == 0 and noncoding[k] == 0:
			continue
		elif coding[k] == 0 and noncoding[k] >0 :
			sum_of_log_ratio_0 -= 1
		else:
			continue
		frame0_count += 1
	'''	
	for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=1):
		if (not coding.has_key(k)) or (not noncoding.has_key(k)):
			continue
		if coding[k]>0 and noncoding[k] >0:
			sum_of_log_ratio_1  +=  math.log( coding[k] / noncoding[k])
		elif coding[k]>0 and noncoding[k] == 0:
			sum_of_log_ratio_1 += 1
		elif coding[k] == 0 and noncoding[k] == 0:
			continue
		elif coding[k] == 0 and noncoding[k] >0 :
			sum_of_log_ratio_1 -= 1
		else:
			continue
		frame1_count += 1
	
	for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=2):
		if (not coding.has_key(k)) or (not noncoding.has_key(k)):
			continue
		if coding[k]>0 and noncoding[k] >0:
			sum_of_log_ratio_2  +=  math.log( coding[k] / noncoding[k])
		elif coding[k]>0 and noncoding[k] == 0:
			sum_of_log_ratio_2 += 1
		elif coding[k] == 0 and noncoding[k] == 0:
			continue
		elif coding[k] == 0 and noncoding[k] >0 :
			sum_of_log_ratio_2 -= 1
		else:
			continue
		frame2_count += 1
	return max(sum_of_log_ratio_0/frame0_count, sum_of_log_ratio_1/frame1_count,sum_of_log_ratio_2/frame2_count)	
	'''
	try:
		return sum_of_log_ratio_0/frame0_count
	except:
		return -1
	

