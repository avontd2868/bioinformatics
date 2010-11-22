import sys
from numpy import *

def usage():
    """Help screen"""
    print "Eventually: Dynamic programming sequence alignment"
    print "Usage: dp.py top_seq side_seq gap-penalty"
    print "Top and side refer to the top and side of the alignment matrix."

def subst_score(a,b):
    # substitution score
    if a==b:
        #print '%s = %s'%(a,b)
        return 1 # match
    else:
        #print '%s != %s'%(a,b)
        return -1 # mismatch

def judge_alignment(top,side):
    if not len(top)==len(side):
        print 'ERROR in scoring alignment\n\tAlignment sequences are not\
        of the same length'
    else:
        score = 0.0
        for i in range(len(top)):
            t_char = top[i]
            s_char = side[i]
            if t_char == s_char:
                score+=1
        ratio = '%s/%s'%(str(int(score)),str(len(top)))
        percent = round(score/len(top)*100,4)
        return ratio, percent
        

def create_dp_matrices(top_seq,side_seq,gap_penalty):
    m=len(side_seq)
    n=len(top_seq)
    ###
    ### initialize matrices
    trace_matrix=zeros([m+1,n+1])
    score_matrix=zeros([m+1,n+1])
    for i in range(1,m+1): # intialize first column
        score_matrix[i,0]=score_matrix[i-1,0]+gap_penalty
        trace_matrix[i,0]=0
    for j in range(1,n+1): #initialize first row
        score_matrix[0,j]=score_matrix[0,j-1]+gap_penalty
        trace_matrix[0,j]=1
    # TRACEBACK MATRIX INITIALIZATION HERE
    #
    ### Now recursion
    ###
    for i in range (1, m + 1):
        for j in range (1, n + 1):
            local_scores=[score_matrix[i-1,j]+gap_penalty, \
                          score_matrix[i,j-1]+gap_penalty, \
                          score_matrix[i-1,j-1]+subst_score(side_seq[i-1],top_seq[j-1])]
            # ARGMAX USAGE HERE
            trace_matrix[i,j]=argmax(local_scores)
            score_matrix[i,j]=max(local_scores)
    return score_matrix, trace_matrix

def trace_back(top_seq, side_seq, trace_matrix):
    (m,n)=trace_matrix.shape
    top_align=[]
    side_align=[]
    
    #print 'm,n:', m,n
    while m > 1 and n > 1:
        #print top_seq[n-2]
        #print side_seq[m-2]
        
        #print trace_matrix[m-1,n-1]
        if trace_matrix[m-1,n-1] == 0:
            top_align.append('-')
            side_align.append(side_seq[m-2])
            #print 'm,n:', m,n
            #print top_align
            #print side_align
            m = m - 1
        elif trace_matrix[m-1,n-1] == 1:
            top_align.append(top_seq[n-2])
            side_align.append('-')
            #print 'm,n:', m,n
            #print top_align
            #print side_align
            n = n - 1
        elif trace_matrix[m-1,n-1] == 2:
            top_align.append(top_seq[n-2])
            side_align.append(side_seq[m-2])
            #print 'm,n:', m,n
            #print top_align
            #print side_align
            m = m - 1
            n = n - 1
        else:
            print 'ERROR in traceback matrix\n\t \
            Unrecognized value %s received at [%s,%s]\
            '%(str(trace_matrix[m-1,n-1]),str(m-1),str(n-1))
    #empty the rest of the values out (we reached the end of
    #one sequence but not necessarily both of them, yet
    while m > 1:
        top_align += top_seq[m-2]
        side_align += '-'
        #print 'm,n:', m,n
        #print top_align
        #print side_align
        m -= 1
    while n > 1:
        top_align += '-'
        side_align += side_seq[n-2]
        #print 'm,n:', m,n
        #print top_align
        #print side_align
        n -= 1
            
    #print 'm,n:', m,n
    top_align.reverse()
    side_align.reverse()
    return top_align, side_align

def lab5_pair_tracebacks(pair, gap_penalty):
    """
    Just a function to make my code prettier and facilitate different pairs without replicating
    all my print statements respectively.
    """
    #print pair[0]
    #print pair[1]
    sm,tm = create_dp_matrices(pair[0],pair[1],gap_penalty)
    print ''.join(['Score matrix for %s and %s with gap penalty of %s:\n'
                   %(pair[0],pair[1],str(gap_penalty)),str(sm),
                   '\nTraceback matrix:\n',str(tm)])
    tb_top, tb_side = trace_back(pair[0],pair[1],tm)
    print tb_top
    print tb_side
    ratio, percent = judge_alignment(tb_top,tb_side)
    print 'Alignment matches: %s [%s%%]'%(ratio,percent)

### MAIN ###
if __name__ == "__main__":


###
# Part A
###
##    if len(sys.argv) < 3:            
##        usage()
##        sys.exit()
##    else:
##        top_seq = sys.argv[1]
##        side_seq = sys.argv[2]
##        gap_penalty = float(sys.argv[3])
##
##        m=len(top_seq)-1
##        while m>=0:
##            if len(top_seq) != len(side_seq):
##                print 'need equal lengths in this example.'
##                sys.exit()
##            print top_seq[m], ' ', side_seq[m]
##            print subst_score(top_seq[m],side_seq[m])
##            m-=1

###
# Part B
###
    small_gap_penalty = -.5
    gap_penalty = -1
    large_gap_penalty = -2
    larger_gap_penalty = -3.5
    
    #pair1 = ('WHAT','WHY')
    #lab5_pair_tracebacks(pair1, gap_penalty)
    
    #pair2 = ('WHEY','WHY')
    #lab5_pair_tracebacks(pair2, gap_penalty)
    
    #pair3 = ('AAAC','AGC')
    #lab5_pair_tracebacks(pair3, gap_penalty)
    
    pair4 = ('IQIFSFIFRQEWNDA','QIFFIFRMSVEWND')
    lab5_pair_tracebacks(pair4, larger_gap_penalty)

    
