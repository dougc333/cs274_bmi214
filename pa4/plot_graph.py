
import argparse
import networkx as nx
from networkgen import NetworkGen
import matplotlib.pyplot as plt
import os
class PlotGraph:
    def __init__(self,arg1, arg2, arg3):
        self.edgelist_file = arg1
        self.pnodes_file = arg2
        self.outputdir = arg3
        self.G=nx.Graph()
        self.H=None
        self.init()
        
    def init(self):
        '''
        draw and do shit on adj list sing protein_nodes.csv 
        '''
        name,color = self.read_protein(self.pnodes_file)

        with open(self.edgelist_file) as fh:
            lines = fh.readlines()
            for l in lines:
                if(len(l)>0):
                    p = l.split()
                    first=p[0]
                    second = p[1]
                    self.G.add_node(first)
                    self.G.add_node(second)
                    self.G.add_edge(first,second)
            mapping={}

        #rename the nodes and color them. 
        for x in list(self.G.nodes):
            #make mapping
            mapping[x]=name[x]
        self.H=nx.relabel_nodes(self.G,mapping)
        colormap=[]
        for node in self.G:
            colormap.append(color[node])
        nx.draw(self.H,node_color = colormap)
        plt.savefig(os.path.join(self.outputdir,"network.png"))
        #plt.plot()
        plt.figure(figsize=(8, 8),dpi=150)
        
    def read_protein(self,filename):
        name={}
        color={}
        with open(filename) as fh:
            lines = fh.readlines()
            for idx in range(0,len(lines)):
                uni_acc=lines[idx].split(sep=",")[0]
                uni_id=lines[idx].split(sep=",")[1]
                ind=lines[idx].split(sep=",")[2].strip()
                name[uni_acc]=uni_id
                if ind=="bp":
                    color[uni_acc] = "red"
                elif ind=="bp;diabetes":
                    color[uni_acc] = "purple"
                elif ind=="bp;cholesterol":
                    color[uni_acc] = "green"
                elif ind=="bp;cholesterol;diabetes":
                    color[uni_acc] = "blue"
        #print(name['P21918'])
        #print(color['P21918'])
        return name,color

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('arg1', type=str, help='network_edgelist.txt')
    parser.add_argument('arg2', type=str, help='protein_nodes.csv')
    parser.add_argument('arg3', type=str, help='file path to output figure location')
    
    args = parser.parse_args()
    #print(args.arg1,args.arg2,args.arg3)
    pg  = PlotGraph(args.arg1,args.arg2,args.arg3)
    