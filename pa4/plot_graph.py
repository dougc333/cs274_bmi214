
import argparse
from networkx import 
from networkgen import NetworkGen

class PlotGraph:
    def __init__(self,arg1, arg2, arg3):
        self.edgelist = arg1
        self.pnodes = arg2
        self.outputdir = arg3
        
        self.init()
        
    def init(self):
        '''
        draw and do shit on adj list sing protein_nodes.csv 
        '''
        with open(self.edgelist) as fh:
            lines=h=fh.readlines()
            for l in lines:
                print(l.split()[0],l.split()[0])


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('arg1', type=str, help='network_edgelist.txt')
    parser.add_argument('arg2', type=str, help='protein_nodes.csv')
    parser.add_argument('arg3', type=str, help='file path to output figure location')
    
    args = parser.parse_args()
    #print(args.arg1,args.arg2,args.arg3)
    pg  = PlotGraph(args.arg1,args.arg2,args.arg3)
    