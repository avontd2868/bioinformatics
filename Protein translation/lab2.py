# Author: Noah M. Jorgenson
# Course: Bioinformatics
# Lab 2

import sys
import numpy as np
import matplotlib.pyplot as plt

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def product(*args, **kwds):
    """
    Returns an iterable with tuples of the pairs of each character
    in the *args string
    """
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

def generate_nuc_dicts():
    nuc_counts = {}
    nuc_probs = {}

    for nuc in 'ATGC':
        nuc_counts[nuc]=0
        nuc_probs[nuc]=0
    return nuc_counts, nuc_probs

def generate_dinuc_dicts():
    """
    Generates the dictionary for dinucs with the dinuc pair
    string as the key to the hash

    Initializes the values as 0
    """
    dinuc_counts = {}
    dinuc_probs = {}
    for pair in product('ATGC',repeat=2):
        dinuc_counts[''.join([pair[0],pair[1]])]=0
        dinuc_probs[''.join([pair[0],pair[1]])]=0
    return dinuc_counts, dinuc_probs
    
def read_fasta(file):
    seq=''
    header=''
    try:
        fasta_file=open(file,'r')
    except IOError:
        print 'Cannot find file ', file, '.'
        sys.exit(1)
        
    for line in fasta_file:
        if line.startswith(">"):
            header += line.strip()[1:]
        else:
            seq += line.strip()
    return header, seq

def dinuc_probs(dna_seq):
    dna_len = len(dna_seq)
    dinuc_counts, dinuc_probs = generate_dinuc_dicts()
    
    total_pairs = float(0)
    for key in dinuc_counts:
        count=sum([dna_seq.count(key,i,i+2) for i in range(dna_len-1)])
        dinuc_counts[key]=count
        total_pairs += count
    for key in dinuc_probs:
        prob=(dinuc_counts[key]/total_pairs)
        dinuc_probs[key]=prob

    return dinuc_counts, dinuc_probs

def dinuc_probs_by_nucprod(dna_seq):
    dna_len = float(len(dna_seq))
    nuc_counts, nuc_probs = generate_nuc_dicts()
    
    for nuc in nuc_counts:
        count=dna_seq.count(nuc)
        nuc_counts[nuc]=count
        nuc_probs[nuc]=(count/dna_len)

    trash, dinuc_probs_by_nucprod = generate_dinuc_dicts()

    for pair in dinuc_probs_by_nucprod:
        dinuc_probs_by_nucprod[pair]=(nuc_probs[pair[0]]*nuc_probs[pair[1]])

    return dinuc_probs_by_nucprod
    

def probs_to_scatter(dinuc_probs, dinuc_probs_by_nucprod):
    ordering = []
    for pair in product('ATGC',repeat=2):
        ordering.append(''.join([pair[0],pair[1]]))

    dinuc_probs_list = []
    dinuc_probs_by_nucprod_list = []
    for pair in ordering:
        dinuc_probs_list.append(dinuc_probs[pair])
        dinuc_probs_by_nucprod_list.append(dinuc_probs_by_nucprod[pair])

    x_axis = np.array(dinuc_probs_list)
    y_axis = np.array(dinuc_probs_by_nucprod_list)
    plt.scatter(x_axis, y_axis)
    plt.xlabel('Empirical')
    plt.ylabel('Pairwise')
    plt.title('Dinucleotide Probabilities')
    plt.show()

def translate_dna(dna_seq):
    translated_dna = ''
    for n in range(0,len(dna_seq),3):
        translated_dna=''.join([translated_dna,gencode[dna_seq[n:n+3]]])
    return translated_dna

def calc_translation_accuracy(trans_dna, BLAST_trans):
    if len(trans_dna)!=len(BLAST_trans):
        return 'ERROR'
    
    index = 0
    wrong = float(0)
    wrong_indeces = []
    for amino_acid in BLAST_trans:
        if amino_acid != trans_dna[index]:
            wrong+=1
            wrong_indeces.append(index)
        index+=1
    return (wrong/len(trans_dna)), wrong_indeces
            
        
        


########
# MAIN #
########

try:
    filename=sys.argv[1]
except IndexError:
    print 'Input filename not provided.'
    print 'usage: python ' + sys.argv[0] + ' WhatSequenceAmI.fasta'
    sys.exit(1)

#1
name, dna_sequence=read_fasta(filename)

#2
dinuc_counts, dinuc_probs = dinuc_probs(dna_sequence)
dinuc_probs_by_nucprod = dinuc_probs_by_nucprod(dna_sequence)
probs_to_scatter = probs_to_scatter(dinuc_probs, dinuc_probs_by_nucprod)

#3
translated_dna = translate_dna(dna_sequence)
BLAST_translation = 'MKAKLLVLLCAFTATYADTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCRLKGIAPLQLGNCSVAGWILGNPECESLFSKESWSYIAETPNPENGTCYPGYFADYEELREQLSSVSSFERFEIFPKESSWPNHTVTKGVTASCSHNGKSSFYRNLLWLTEKNGLYPNLIKSYVNNKEKEVLVLWGVHHPSNIGDQRAIYHTENAYVSVVSSHYSRRFTPEIAKRPKVRDQEGRINYYWTLLEPGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNASMNECDAKCQTPQGAINSSLPFQNVHPVTIGECPKYVRSTKLRMVTGLRNIPSIQSR'
accuracy,wrong_indeces = calc_translation_accuracy(translated_dna, BLAST_translation)
print ''.join(['Translation had ',str(accuracy*100), '% inaccuracies according to BLAST database.'])
if len(wrong_indeces)>0:
    wrong_indeces_str = ''
    for index in wrong_indeces:
        wrong_indeces_str = ''.join([wrong_indeces_str,str(index)])
    print ''.join(['The incorrect indeces were: ',wrong_indeces_str])







