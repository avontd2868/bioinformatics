####################
# Noah M. Jorgenson
# Bioinformatics
# Lab 4
#

import sys
import scipy
from numpy import *
import matplotlib.pyplot as plt



def gen_jukes_mat(a):
    M = mat([[1-a,a/3,a/3,a/3],[a/3,1-a,a/3,a/3],[a/3,a/3,1-a,a/3],[a/3,a/3,a/3,1-a]])
    return M

def generations(juke_m, num_gen, p0):
    p = p0
    for i in range(num_gen):
        p_next=juke_m*p
        p=p_next
    return p

def generations_convergence(juke_m, p0):
    max_gen=3e6
    tol=1e-8

    diff=sum(1)
    curr_gen=0
    p=p0
    #while we are not outside tolerance for accepting convergence
    while (diff>tol):
        p_next=juke_m*p
        diff = sum(abs(p_next-p)) #calculate difference between new prob vector and last prob vector
        p=p_next
        curr_gen+=1 #increment generation
        if curr_gen > max_gen:
            print 'No convergence, exiting via maximum generation limit.\n'
            break #break while loop if we are about to exceed max generations
    return p, curr_gen #return the prob vector and the number of generations it took

def genetic_drift(A,N,G):
    """
    Model a population of chromosomes with a single gene (boolean 1/0; A/a)

    A = initial number of A genes
    
    """
    #BLOCK 1: Input simulation parameters: A, N and G
    A_prop=zeros(G)     #Initialize the 1D array of allele A proportions at each generation
    pop=zeros(N)        #Initial population of fixed size N
    pop[0:A]=ones(A)    #A is the number of genes with allele A (1's)
    A_prop[0]=sum(pop)/N#Fraction of allele A in the population at gen 0
    pop_next=zeros(N)   #Initialize the population for the next generation

    #BLOCK 2: Carry out the simulation from generation 1 to G
    #Pick one of the alleles at random from the previous population
    for gen in range(1,G,1): #for each generation
        for i in range(N): #for each person in population
            pop_next[i]=pop[random.randint(0,N-1)] #assign the person the value of some previous person from last population
        A_prop[gen]=sum(pop_next)/N # count the number of 1s
        pop=pop_next #set the next population to the generated population

    return A_prop #return the number of genes that have the A gene (1's) after the population

def main():
    numA = 750   #Number of A's to initalize the population with
    popN = 1000   #Size of the population within a generation
    gens = 100  #Number of generations to run simulation over

    t=linspace(1,gens,gens)
    props=genetic_drift(numA,popN,gens)
    plt.plot(t,props)
    plt.show()
    

## PART A
##########

p0 = mat([.1,.4,.25,.25]).T

m_1=gen_jukes_mat(1e-1) #a=1e-1
m_2=gen_jukes_mat(1e-2) #a=1e-2
m_5=gen_jukes_mat(1e-5) #a=1e-5
m_6=gen_jukes_mat(1e-6) #a=1e-6

num_gen = 15
p_m1 = generations(m_1, num_gen, p0)
print ''.join(['[a=1e-1] Probability vector (transposed) at 15th generation:\n\t',str(p_m1.T)])

p_m1,gen_m1 = generations_convergence(m_1,p0)
print ''.join(['[a=1e-1] [Convergence at generation %s]\nResultant probability vector(transposed):\n\t'%gen_m1,str(p_m1)])

p_m2,gen_m2 = generations_convergence(m_2,p0)
print ''.join(['[a=1e-2] [Convergence at generation %s]\nResultant probability vector(transposed):\n\t'%gen_m2,str(p_m2)])

p_m5,gen_m5 = generations_convergence(m_5,p0)
print ''.join(['[a=1e-5] [Convergence at generation %s]\nResultant probability vector(transposed):\n\t'%gen_m5,str(p_m5)])

p_m6,gen_m6 = generations_convergence(m_6,p0)
print ''.join(['[a=1e-6] [Convergence at generation %s]\nResultant probability vector(transposed):\n\t'%gen_m6,str(p_m6)])

M = gen_jukes_mat(1e-1)
vals, vecs = linalg.eig(M)
print vals
print vecs

## PART B
########

main()










