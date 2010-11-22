#!/usr/bin/env python
from __future__ import division
import numpy as N
import sys, csv, optparse

import itertools

#global dtype
dtype_allele = N.dtype([('A','float64'),('T','float64'),('C','float64'),('G','float64'),('?','float64')])

#valid alphabetically arranged combinations
valid = ['AA','AC','AG','AT','CC','CG','CT','GG','GT','TT']

#Count allele
def allele_count(allele_array,nucleotide,normalize,missing):
	if missing == 'false':
		allele_array = filter('?'.__ne__,allele_array)
	
	if normalize:
		return allele_array.count(nucleotide)/len(allele_array)
	else:
		return allele_array.count(nucleotide)

#Create allele_count rec array
def allele_freq(allele_array, normalize, missing):
	cast = N.cast
	temp = []
	
	#find unique nucleotides in allele
	#unique = list(frozenset(allele_array))
	#dtype_allele = N.dtype([(unique[0],'float64'),(unique[1],'float64')])
	
	
	for i in xrange(5):
		temp.append(cast[dtype_allele[i]](allele_count(allele_array,dtype_allele.names[i],normalize,missing)))

	count = N.rec.array(temp,dtype = dtype_allele)
	
	return count

#parse map file
def parseMap(map_file):
	dtype = N.dtype([('Chromosome', N.str_,3), ('RS', N.str_,12), \
			 ('Genetic_Distance', 'float32'), ('Base-Pair_Position','int32')])
	cast = N.cast
	data = [[] for dummy in xrange(len(dtype))]	
	temp = []

	reader = csv.reader(map_file, delimiter='\t')
	for row in reader:
		temp.append(row)
	
	temp = map(list, zip(*temp))
	
	for i in xrange(len(dtype)):
		data[i] = cast[dtype[i]](temp[i])
	
	mapArray = N.rec.array(data, dtype=dtype)
	
	return mapArray

#class section
class parsePlink:
	def __init__(self,datafile):
		ped_file = open("%s.ped" % datafile,  'r')
		map_file = open("%s.map" % datafile,  'r')
		
		cast = N.cast
		self.data = []
		temp = []
		
		self.mapArray = parseMap(map_file)
		
		#SNP list
		self.snp_list = self.mapArray['RS']
			
		dtype = N.dtype([('Family_ID', N.str_,10), ('Individual_ID', N.str_,10), \
			 ('Paternal_ID', N.str_,10), ('Maternal_ID',N.str_,10), \
			('Sex','int32'),('Phenotype','int32')] + list((i,N.str_,2) for i in self.snp_list))

		reader = csv.reader(ped_file, delimiter=' ')
		
		for row in reader:
			#replace '' with '?', adds 5 seconds
			row = row[:6] + list(x if x in ['A','T','C','G'] else '?' for x in row[6:])
			temp.append(row[:6] + list((row[i]+row[i+1]) \
					for i in xrange(6,len(row),2)))

		temp = map(list, zip(*temp))

		number_of_columns = len(dtype)	
		for i in xrange(number_of_columns):
			self.data.append(cast[dtype[i]](temp[i]))
		
		self.pedArray = N.rec.array(self.data, dtype=dtype)
		
		#Total Males and Females
		self.total_male = self.pedArray['Sex'].tolist().count(1)
		self.total_female = self.pedArray['Sex'].tolist().count(2)
		self.total_unknown_sex = self.pedArray['Sex'].tolist().count(0)		
		
	#return genotype for given SNP
	def get_genotype(self, rs_number):
		return self.pedArray[rs_number]
	
	#splits allele in no specific order	
	def split(self, rs_number):
		#get genotype
		genotype = self.pedArray[rs_number]
		
		#arrange genotype alphabetically
		genotype_arranged = [i if i in valid else i[::-1] for i in genotype]
		
		alleles = zip(*genotype_arranged)
		return alleles[0],alleles[1]
		
	#returns common,minor allele in order
	def split_alleles(self, rs_number):
		#get genotype
		genotype = self.pedArray[rs_number]
		
		#arrange genotype alphabetically
		genotype_arranged = [i if i in valid else i[::-1] for i in genotype]
		
		alleles = zip(*genotype_arranged)
		
		#don't check for common allele if both alleles are the same	
		if alleles[0] == alleles[1]:
			return alleles[0],alleles[1]	
		
		#get allele frequencies
		allele1_frequencies = allele_freq(alleles[0],'false','true')
		allele2_frequencies = allele_freq(alleles[1],'false','true')
		
		#if allele has more '?' than nucleotides, return the other allele
		#set creates list of unique elements
		for i in set(alleles[0])-set(['?']):
			if allele1_frequencies['?']>= allele1_frequencies[i]:
				return alleles[1], alleles[0]

		diff_a1 = 0
		diff_a2 = 0		
		#get non_zero indices and
		#get difference accross frequencies
		a1_nonzero = N.nonzero(allele1_frequencies.item())[0]
		a2_nonzero = N.nonzero(allele2_frequencies.item())[0]

		for i in a1_nonzero:
			diff_a1 = abs(allele1_frequencies.item()[i] - diff_a1)
		
		for i in a2_nonzero:
			diff_a2 = abs(allele2_frequencies.item()[i] - diff_a2)
		
		#print allele1_frequencies.item(), allele2_frequencies.item()
		#print diff_a1, diff_a2
		
		if(diff_a1>diff_a2 or diff_a1==diff_a2):
			return alleles[0],alleles[1]
		else:
			return alleles[1],alleles[0]

	def tup_count(self,tup,item):
                count = 0
                for ind in tup:
                        if ind==item: count+=1
                return count
                
	#return minor allele, common allele and minor allele frequency
	def maf(self, rs_number):
		temp = []
		
		#Split allele for given RS
		A1,A2 = self.split(rs_number)
		
		#remove unknowns from data	
		if '?' in A1:
			A1 = filter('?'.__ne__,A1)
		if '?' in A2:
			A2 = filter('?'.__ne__,A2)		
		
		#get a list of unique nucleotides in bot alleles
		unique = list(frozenset(A1).union(frozenset(A2)))
		
		#In case only one nucleotide is present
		if len(unique) == 1:
			return 0.0, unique[0], 0.0
		
		for nucleotide in unique:
                        #print self.tup_count(A1,nucleotide)
                        #print self.tup_count(A2,nucleotide)
			temp.append(self.tup_count(A1,nucleotide)+self.tup_count(A2,nucleotide))
		
		index_minor = N.argmin(temp)	
		minor_allele = unique[index_minor]
		
		if index_minor == 0:
			common_allele = unique[1]
		else:
			common_allele = unique[0]
		
		minor_allele_freq = (self.tup_count(A1,minor_allele)+self.tup_count(A2,minor_allele))/(2*(len(A1)))
		
		return minor_allele, common_allele, minor_allele_freq
	
	def ct_allelic(self, rs_number, phenotype):
		
		#only works for Sex and Phenotype
		if phenotype not in ['Sex','Phenotype']:
			print 'invalid phenotype'
			return 0
		
		#Get major and minor alleles
		a1a2 = self.maf(rs_number)
		A1 = a1a2[0]
		A2 = a1a2[1]
		#get phenotype column
		phenotype = self.pedArray[phenotype]
		#split genotype into alleles
		allele_col1, allele_col2 = self.split(rs_number)
		
		#count total nucleotides, not including unknowns
		total_nucleotides = len(filter('?'.__ne__,allele_col1)) + len(filter('?'.__ne__,allele_col2))
		
		array = N.vstack((phenotype,allele_col1,allele_col2)).transpose()
		
		#init contingency matrix
		matrix = N.zeros((2,2))
		
		for row in array:
			#remove subjects with unknown phenotypes
			if row[0][0] not in ['1','2']:
				array = N.delete(array,0,axis=0)
				total_nucleotides-=2
			#Case
			if row[0][0] == '2':
				matrix[0][0] += row[1:].tolist().count(A1)
				matrix[0][1] += row[1:].tolist().count(A2)
			#Control
			if row[0][0] == '1':
				matrix[1][0] += row[1:].tolist().count(A1)
				matrix[1][1] += row[1:].tolist().count(A2)
		
		
		return matrix/total_nucleotides
			
if __name__ == '__main__':
	sys.exit(main(sys.argv))
