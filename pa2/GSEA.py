import sys
import pandas as pd



class GSEA:
  def __init__(self):
    self.exp_data = None
    self.samp_data = None
    self.gene_data = None

  def load_data(self,expfile,sampfile,genesets):
    df_exp = pd.read_csv(expfile,sep='\t')
    df_samp = pd.read_csv(expfile,sep='\t')
    df_geneset = pd.read_csv(expfile,sep='\t')
    self.exp_data = df_exp.to_numpy()
    self.samp_data = df_samp.to_numpy()
    self.gene_data = df_geneset.to_numpy()
  def get_gene_rank_order(self):
    return
  def get_enrichment_score(self):
    return
  def get_sig_sets(self,p):
    '''
    return list of string at corrected thrholsd. 
    '''
    print("threshold:",p)
    return