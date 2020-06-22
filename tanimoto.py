
import argparse
import pandas as pd
import numpy as np
import os
import itertools

class Tanimoto:
    def __init__(self,arg1, arg2, arg3):
        #nothing here so far
        self.drugs = arg1
        self.targets = arg2
        self.outputfile = arg3
        self._drug_dict = {}
        self._pairTc = None
        self._pairs = None
        self._db_to_prot = {}
        self._prot_to_drug = {}
        self.init()

    def init(self):
        #print("init tanimoto.py")
        self.make_drugdict(self.drugs)
        self.make_dbtoprot_and_prottodrug(self.targets)
        #print("DB01019:",self._db_to_prot['DB01019'])
        self.write_tc(self.outputfile)
        
    def make_drugdict(self,filename):
        '''
        input: drugs.csv
        output: dict key=drugs, value=maccs score
        '''
        df_drugs = pd.read_csv(filename)
        df_numpy = df_drugs.to_numpy()
        for idx in range(0,df_numpy.shape[0]):
            name = df_numpy[idx,:][0]
            #convert to int. Is maccs a set?
            score = set([int(x) for x in df_numpy[idx,:][2].split()])
            self._drug_dict[name] = score
    

    def make_dbtoprot_and_prottodrug(self,filename):
        '''
        input: target.csv
        output: dict key=drug, value=[list of proteins]
        We need both lookups drugs to list of proteins and proteins to list of drugs. Disk IO is expensive
        so once the file is read into memory create all datastructures for Tanimoto and pvalue. 
        '''
        df_targets = pd.read_csv(filename)
        for idx in range(0,df_targets.to_numpy().shape[0]):
            drug = df_targets.to_numpy()[idx][0].strip()
            protein = df_targets.to_numpy()[idx][1].strip()
            if drug not in self._db_to_prot:
               self._db_to_prot[drug] = set([protein])
            else:
                self._db_to_prot[drug].add(protein)
            if protein not in self._prot_to_drug:
                self._prot_to_drug[protein] = set([drug])
            else:
                self._prot_to_drug[protein].add(drug)
        #initialized both datastructures _db_to_prot and _prot_to_drug

    def tanimoto(self,a,b):
        '''
        input: 2 drugs, a and b
        output: tanimoto score of drugs a and b from MACCSKeys
        '''
        return round(float(len(a.intersection(b)))/float(len(a.union(b))),6)

    def singleTc_drug(self, db_a,db_b):
        '''
        input: db_a, db_b 2 drugs 
        output: Tanimoto for both drugs
        '''
        return self.tanimoto(self._drug_dict[db_a], self._drug_dict[db_b])
    
    def singleTc_protein(self,A,B):
        '''
        tanimoto for 2 proteins, this is Tsummary not Tc. 
        '''
        #we should put this in pvalue
        return None

    def write_tc(self,filename):
        self._pairs = list(itertools.combinations(self._drug_dict.keys(),2))
        self._pairTc = {(x[0],x[1]): self.tanimoto(self._drug_dict[x[0]],self._drug_dict[x[1]])  for x in self._pairs}
        fh = open(filename,"w")
        k =list(self._pairTc.keys())
        for idx in range(0,len(k)):
            first = k[idx][0]
            second = k[idx][1]
            tc_score = self._pairTc[(first,second)]
            fh.write(first+","+second+","+str(tc_score)+","+str(self.shared_targets(first,second,self._db_to_prot))+"\n")
        fh.close()


    def shared_targets(self,db1,db2,dbtp):
        if db1 in dbtp.keys():
            target_db1 = dbtp[db1]
        else:
            target_db1=set([])

        if db2 in dbtp.keys():
            target_db2 = dbtp[db2]
        else:
            target_db2 = set([])

        if len(target_db1.intersection(target_db2)) > 0:
            return 1
        return 0
    

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('arg1', type=str, help='drugs.csv location')
    parser.add_argument('arg2', type=str,  help='targets.csv location')
    parser.add_argument('arg3', type=str,  help='outputfile.csv location')
    args = parser.parse_args()
    t = Tanimoto(args.arg1,args.arg2,args.arg3)
    