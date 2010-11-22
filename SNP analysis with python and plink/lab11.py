#Author: Noah M. Jorgenson
#Bioinformatics Lab 11

import plinkParser as pp
#print dir(pp)
"""
['N', '__builtins__', '__doc__', '__file__', '__name__',
'allele_count', 'allele_freq', 'csv', 'division', 'dtype_allele',
'itertools', 'optparse', 'parseMap', 'parsePlink', 'sys', 'valid']
"""

def calc_all_snp_quants(extra_data):
    """
    Calculates the following quantities for each SNP in *.ped data
    """
    for item in extra_data.mapArray:
        chromosome = item[0]
        snp = item[1]
        bp_pos = item[3]        
        A1, A2 = extra_data.split(snp)
        c_m_p = extra_data.ct_allelic(snp, 'Phenotype')
        freq_A1_cases = (c_m_p[0,0] / (c_m_p[0,0]+c_m_p[0,1]))
        freq_A1_contr = (c_m_p[1,0] / (c_m_p[1,0]+c_m_p[1,1]))
        odds_ratio_pheno = calc_odds_ratio(c_m_p)

        print '\t'.join([chromosome, snp, str(bp_pos), str(A1), str(freq_A1_cases),
                         str(freq_A1_contr), str(A2), str(odds_ratio_pheno)])
    

def calc_odds_ratio(cont_tab):
    """
    Calculate the odds ratio of a SNP form its phenotype/A1 contingency table array
    """
    c_m_p = cont_tab
    #frequence of cases with minor allele A1
    # freq_A = a / (a+b)
    # freq_U = c / (c+d)
    freq_A1_cases = (c_m_p[0,0] / (c_m_p[0,0]+c_m_p[0,1]))
    freq_A1_contr = (c_m_p[1,0] / (c_m_p[1,0]+c_m_p[1,1]))
    #print freq_A, freq_U

    #probability of event; p
    # odds = p/(1-p)
    A1_cases_odds = (freq_A1_cases/(1-freq_A1_cases))
    A1_contr_odds = (freq_A1_contr/(1-freq_A1_contr))
    #print A1_cases_odds, A1_contr_odds

    #odds ratio
    #odds of cases having A1 to odds of controls having A1
    A1_odds_ratio = A1_cases_odds/A1_contr_odds
    #print A1_odds_ratio
    return A1_odds_ratio




extra_data = pp.parsePlink('extra')
#print dir(extra_data)
"""
['__doc__', '__init__', '__module__', 'ct_allelic', 'data',
'get_genotype', 'maf', 'mapArray', 'pedArray', 'snp_list',
'split', 'split_alleles', 'total_female', 'total_male', 'total_unknown_sex']
"""

##tot_fem = extra_data.total_female
##tot_mal = extra_data.total_male
##tot_unk = extra_data.total_unknown_sex
###print tot_fem, tot_mal, tot_unk
##
##col_ids = extra_data.pedArray.dtype #column ids in pedArray ped data
###print col_ids
##
##phenotype_states = extra_data.pedArray['Phenotype']
###print phenotype_states
##
##snp_rs_ids = extra_data.mapArray['RS'] #SNP rs ids
###print snp_rs_ids
##
##snp = 'rs17121574'
##print extra_data.get_genotype(snp) #== extra_data.pedArray[snp]
##
##A1, A2 = extra_data.split(snp) #split genotype into alleles
###print A1, A2
##
##freqs = extra_data.maf(snp) #count all occurrences of alles in all genotypes
###print freqs
##print('minor allele = ', freqs[0], ' with frequency ', freqs[2])
##print('common allele = ', freqs[1], ' with frequency ', 1-freqs[2])
##
#########################################
###calculate A1 odds ratio for Phenotype
##
###numpy array of contingency table for frequence of minor allele (A1)
###and common allele (A2) in the cases and controls
##conting_mat_pheno = extra_data.ct_allelic(snp, 'Phenotype')
##odds_ratio_pheno = calc_odds_ratio(conting_mat_pheno)
##print odds_ratio_pheno
##
#########################################
###calculate A1 odds ratio for Sex
##
###numpy array of contingency table for frequence of minor allele (A1)
###and common allele (A2) in the cases and controls
##conting_mat_sex = extra_data.ct_allelic(snp, 'Sex')
##odds_ratio_sex = calc_odds_ratio(conting_mat_sex)
##print odds_ratio_sex

calc_all_snp_quants(extra_data)









