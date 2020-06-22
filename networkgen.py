from pvalue import Pvalue
import argparse
import pandas as pd
import itertools
class NetworkGen:
    def __init__(self,arg1, arg2, arg3):
        self.drugs = arg1
        self.targets = arg2
        self.protein_nodes = arg3
        self.pvalue = Pvalue(500 ,214,self.drugs,self.targets,"","")
        self.init()
        
    def init(self):
        df_nodes = pd.read_csv(self.protein_nodes)
        proteins = []
        for row in df_nodes.iterrows():
            #print(type(row),len(row))
            #print(row[0],row[1].to_dict())
            u_acc = row[1].to_dict()['uniprot_accession']
            u_id = row[1].to_dict()['uniprot_id']
            ind = row[1].to_dict()['indications']
            #print(u_acc, u_id, ind)
            proteins.append(u_acc)

        #print(len(proteins))    
        pairs = list(itertools.combinations(proteins,2))
        #print(len(pairs))
        #x = pairs[0]
        #pb = self.pvalue.pb(x[0],x[1])
        #print(pb)
        check={}
        fh = open("network_edgelist.txt","w")
        for x in pairs:
            pb = self.pvalue.pb(x[0],x[1])
            if (x[0],x[1]) in check:
                print('DUPLICATE')
            check[(x[0],x[1])] = 1
            check[(x[1],x[0])] = 1
            #print("pb:",pb)
            if pb<=0.05:
                fh.write(x[0]+" "+x[1]+"\n")
        fh.close()

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('arg1', type=str, help='drugs.csv location')
    parser.add_argument('arg2', type=str, help='targets.csv location')
    parser.add_argument('arg3', type=str, help='protein_nodes location')
    
    args = parser.parse_args()
    #print(args.arg1,args.arg2,args.arg3)
    ngen  = NetworkGen(args.arg1,args.arg2,args.arg3)
    

