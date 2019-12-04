
import argparse
import pandas as pd
import numpy as np
import random
from tanimoto import Tanimoto
import itertools

class Pvalue:
    def __init__(self,n, r, arg3,arg4,arg5,arg6):
        self.N = n
        self._random_seed = r
        random.seed(self._random_seed)
        self.drug = arg3
        self.targets = arg4
        self.A = arg5
        self.B = arg6
        self.tanimoto = Tanimoto(self.drug, self.targets, "outputfile_notused.txt")
        self.init()
        self.check={}
    
    def init(self):
        if len(self.A)>0 and len(self.B)>0:
            self.p()
        

    def p(self):
        pval = self.pb(self.A,self.B)
        print("pval",pval)
        return pval

    def tsummary_AB(self,A,B):
        '''
        input: A,B uniprot id, prot_to_drug dict of prot to all drugs which have this protein as target, 
        dbtp: db to protein dictionary, tc tanimoto coeffieicnt of all db pairs
        outptu: tsummary(A,B)
        self.tanimoto._prot_to_drug, self.tanimoto._db_to_prot(never used)
        '''
        db_A = self.tanimoto._prot_to_drug[A]
        db_B = self.tanimoto._prot_to_drug[B]
        a = list(db_A)
        b = list(db_B)
        sum_a_b=0.
        num=0
        for idx in range(0,len(a)):
            for jdx in range(0,len(b)):
                tc_ab = self.tanimoto.singleTc_drug(a[idx],b[jdx])
                if tc_ab > 0.5:
                    sum_a_b +=tc_ab
                num +=1
        return sum_a_b

    def tsummary_B(self,A,B):
        ligA = self.get_ligandSet(A)
        ligB = self.get_ligandSet(B)
        tb_sum = 0
        
        for idx in range(0,len(ligA)):
            for jdx in range(0,len(ligB)):
                tc_b = self.tanimoto.singleTc_drug(ligA[idx],ligB[jdx])
                if tc_b > 0.5:
                    tb_sum +=tc_b
        return tb_sum

    #random.choice()
    def get_ligandSet(self,A):
        list_drugs = list(self.tanimoto._drug_dict.keys())
        lenA = len(self.tanimoto._prot_to_drug[A])
        sample_rep  = [random.choice(list_drugs) for _ in range(lenA)]
        return sample_rep


    def pb(self,A,B):
        '''
        sum Tb>Tsummary
        '''
        tot = 0
        tsum=self.tsummary_AB(A,B)
        for _ in range(0,self.N):
            tb_val = self.tsummary_B(A,B)
            if tb_val >=tsum:
                 tot +=1        
        return tot/self.N


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, default=500, help='n, number of iterations')
    parser.add_argument('-r', type=int, default=214, help='random seed')
    parser.add_argument('arg3', type=str, help='drugs.csv location')
    parser.add_argument('arg4', type=str, help='targets.csv location')
    parser.add_argument('arg5', type=str, default='P21918', help='protein A')
    parser.add_argument('arg6', type=str, default='P18089', help='protein B')
    
    args = parser.parse_args()
    pvalue = Pvalue(args.n, args.r, args.arg3,args.arg4,args.arg5,args.arg6)

    