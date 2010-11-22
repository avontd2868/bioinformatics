#Author: Noah M. Jorgenson

import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


########
# Main #
########

####10-6-2010 Pre-Lab 7

##pdb1_filename = "1TGH.pdb"
##pdb1_code="1TGH"
##
##parser=PDBParser() #Instance of parser class
##structure=parser.get_structure(pdb1_code, pdb1_filename)
##model=structure[0]
##chainA=model["A"]
##TBPresidue219=chainA[219]
###print TBPresidue219.resname #GLY
###print TBPresidue219['CA'].coord #coordinates of chain CA
##
###C1 prime coords
##chainC = model['C']
##DNAresidue118 = chainC[118]
##C1prime = DNAresidue118["C1'"]
##print C1prime.coord
##
##
###euclidean distance function = dist(x,y)
##CA=TBPresidue219['CA']
##print CA.coord
##
##C1prime=DNAresidue118["C1'"]
##print C1prime.coord

def find_nearest_residue():
    pdb1_filename = "1TGH.pdb"
    pdb1_code="1TGH"

    parser=PDBParser()
    structure=parser.get_structure(pdb1_code, pdb1_filename)
    model=structure[0]
    chainA=model["A"]
    chainC=model["B"]  # DNA chain

    nearest_residues_dict = {}
    for dna_res in chainC.get_list():
        dna_resid = res_to_resid(dna_res)
        # find the distances of this nucleotide to all protein side chains
        if dna_resid:
            dna_resnum=dna_resid
            dna_resname=chainC[dna_resnum] 
            dna_c1p=dna_resname["C1'"]
            dna_coord=dna_c1p.coord

            # rows are protein residues, and columns are residue id and distance to dna nucl
            distance_array=np.ones((len(chainA),2))*999 
            i=0 # distance array row index
            for prot_res in chainA.get_list():
                prot_resid = res_to_resid(prot_res)
                if prot_resid:
                    prot_resnum=prot_resid
                    # don't want' W's                                                     
                    prot_resObj=chainA[prot_resid]
                    prot_CA=prot_resObj["CA"]
                    prot_coord=prot_CA.coord
                    distance_array[i,0]=prot_resnum
                    distance_array[i,1]=dist(dna_coord,prot_coord)
                i=i+1  # this is part of for loop
            min_dist = find_min_dist(distance_array)
            nearest_residues_dict[dna_resnum]=min_dist
    return nearest_residues_dict

def find_min_dist(distance_array):
    min_residue_mask= np.argmin(distance_array[:,1])
    min_dist_residue=distance_array[min_residue_mask]
    return min_dist_residue

def dist(dna, prot):
    """
    Find the distance between the two coordinates (numpy arrays)
    """
    diff = dna-prot
    dist = np.sqrt(np.sum(np.power(diff,2)))
    return dist

def res_to_resid(res):
    """
    Returns the resnum of a res where the hetfield is not 'W'
    """
    resid=res.get_id()
    resnum=resid[1]
    hetfield=resid[0]
    if hetfield[0]==' ':
        return resnum
    else:
        return False


## 10-8-2010 - Lab 7
def mckinney_complement(seq):
    complDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    complList= [complDict[base] for base in seq]
    complseq=''.join(complList)
    return complseq


def find_palindromes(seq,size):
    palindromeDict={}
    seq_len=len(seq)
    for i in range(seq_len-size+1):
        sub_seq=seq[i:i+size]
        biopy_seq=Seq(sub_seq, generic_dna)
        #if sub_seq[::-1] == mckinney_complement(sub_seq):
        if sub_seq[::-1] == biopy_seq.complement().tostring():
            palindromeDict[i]=sub_seq
    return palindromeDict

def find_palindromes_variable(seq):
    """
    Go through for each length (starting with 2) by steps of 2 (even numbers are
    the only lenghts capable of being palindromes) and determine if there are
    palindromes in seq.

    Return once you go through the entire sequence or once you reach a length with no palindromes
    (once you find no palindromes there won't be any in the future)
    """
    palindromeDict={}
    seq_len = len(seq)
    max_size = seq_len
    for size in range(2,max_size+1,2):
        palindromeDict[size]={}
        found = False
        for i in range(seq_len-size+1):
            sub_seq=seq[i:i+size]
            biopy_seq=Seq(sub_seq, generic_dna)
            #if sub_seq[::-1] == mckinney_complement(sub_seq):
            if sub_seq[::-1] == biopy_seq.complement().tostring():
                found=True
                palindromeDict[size][i]=sub_seq
        if not found:
            #return if no palindromes found
            return palindromeDict
    return palindromeDict
    

# now apply the function
dna='AACAATGCCATGATGATGATTATTACGACACAACAACACCGCGCTTGACGGCGGCGGATGGATGCCGCGATCAGACGTTCAACGCCCACGTAACGTAACGCAACGTAACCTAACGACACTGTTAACGGTACGAT'
#palDict= find_palindromes(dna,4)
palDict= find_palindromes_variable(dna)

print '\nPalindromes'
for location in sorted(palDict.keys()):
    print location, '\t', palDict[location]

#find the nearest residues for each palindrome nucleotide
nearest_residues_dict = find_nearest_residue()
print '\nNearest residues for DNA sequence chain'
for res in sorted(nearest_residues_dict.keys()):
    print res, '\t', nearest_residues_dict[res]

