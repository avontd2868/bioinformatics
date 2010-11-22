# Author: Noah M. Jorgenson
# Course: Bioinformatics
# Lab 2

import sys
from numpy import *
from decimal import *

INDEX = {'A':0, 'T':1, 'C':2, 'G':3}
NUCS = 'ATGC'
FILE_test1 = 'test1.fasta'
FILE_testN = 'testN.fasta' 

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

def generate_dinuc_dicts():
    """
    Generates the dictionary for dinucs with the dinuc pair
    string as the key to the hash

    Initializes the values as 0
    """
    dinuc_counts = {}
    dinuc_probs = {}
    for pair in product(NUCS,repeat=2):
        dinuc_counts[''.join([pair[0],pair[1]])]=0
        dinuc_probs[''.join([pair[0],pair[1]])]=0
    return dinuc_counts, dinuc_probs

def calc_nuc_probs(dna_seq):
    """
    Calculates th probability of a nucleotide occuring in
    a given dna sequence
    """
    nuc_probs = {}
    dna_len = float(len(dna_seq))
    for key in NUCS:
        count = dna_seq.count(key)
        prob = count/dna_len
        nuc_probs[key]=prob
    return nuc_probs

def calc_dinuc_probs(dna_seq):
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
        
        
def markov_chain(dna_seq):
    """
    Returns the Markov chain transition matrix trained on the DNA
    sequence provided as argument.

    Remember the conditional probability of P(A,C) for example,
    is actually the join probability (empirically observed) of
    the pair CA (reversed), divided by the total probability
    of the nucleotide C.
    """
    t = mat(zeros((4,4)))
    nuc_probs = calc_nuc_probs(dna_seq)
    dinuc_counts, dinuc_probs = calc_dinuc_probs(dna_seq) 

    for nuc1 in INDEX:
        for nuc2 in INDEX:
            #calculate the conditional probability
            joint_prob_pair = ''.join([nuc2,nuc1])
            joint_prob = dinuc_probs[joint_prob_pair]
            cond_p = (joint_prob/nuc_probs[nuc2])
            t[INDEX[nuc1],INDEX[nuc2]] = cond_p
    return t, nuc_probs

def seq_prob_calc(t, nuc_probs, dna_seq_filename):
    """
    Having trained the transitional matrix of conditional prbabilities
    calculated by the Markov chain, calculate the probability for a
    sequence in the argument dna sequence to occur.
    """
    #get the sequence we are calculating probability for
    header, dna_seq = read_fasta(dna_seq_filename)
    seq_len = len(dna_seq)
    print 'Calculating probability of sequence %s[Length: %s] occuring...'%(header,seq_len)

    zero_hit = False
    p_total = Decimal(1)
    for i in range(len(dna_seq)-1):
        t_prob = Decimal(str(t[INDEX[dna_seq[i]],INDEX[dna_seq[i+1]]]))
        p_total = p_total * t_prob
        if p_total == 0 and not zero_hit:
            zero_hit = i
    if zero_hit:
        print 'Hit zero at %sth nucleotide in sequence.\n'%str(zero_hit)

    #now account for that last marginal probability of p(x)
    p_total = p_total*Decimal(str(nuc_probs[dna_seq[-1]]))

    return p_total
        

########
# MAIN #
########

try:
    filename=sys.argv[1]
except IndexError:
    print 'Input filename not provided.'
    print 'usage: python ' + sys.argv[0] + ' WhatSequenceAmI.fasta'
    sys.exit(1)

#Part A
name, dna_seq=read_fasta(filename)

trans_matrix, nuc_probs = markov_chain(dna_seq)
print 'Conditional probability matrix via markov chain trained by input DNA sequence:'
print trans_matrix
print

#train markov 
test1_seq_prob = seq_prob_calc(trans_matrix, nuc_probs, FILE_test1)
testN_seq_prob = seq_prob_calc(trans_matrix, nuc_probs, FILE_testN)

print 'Probability %s: %s'%(FILE_test1,str(test1_seq_prob))
print 'Probability %s: %s'%(FILE_testN,str(testN_seq_prob))

#the p_total is most likely not sufficient since it is such a
#small probability, so we will take the log-odds ratio of the two
#probabilities

test1_seq_prob_ln = Decimal.ln(test1_seq_prob)
testN_seq_prob_ln = Decimal.ln(testN_seq_prob)
print 'Natural (base e) logarithm of probability for %s: %s'%(FILE_test1,str(test1_seq_prob_ln))
print 'Natural (base e) logarithm of probability for %s: %s'%(FILE_testN,str(testN_seq_prob_ln))

ratio_prob_ln = Decimal.ln(test1_seq_prob/testN_seq_prob)
print 'Natural (base e) logarithm of ratio of probability for %s and %s: %s'%(FILE_test1,FILE_testN,str(ratio_prob_ln))

