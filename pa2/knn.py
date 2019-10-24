
import sys
import pandas as pd

class KNN(object):
  def __init__(self):
    self.k = 5
    self.fn = None
    self.exp_columns = None
    self.exp_data = None
    self.samp_columns = None
    self.samp_data = None

  def load_data(self,expfile, sampfile):
    #process exp file
    df_exp = pd.read_csv(expfile,sep="\t")
    self.exp_columns = df_exp.columns()
    self.exp_data = df_exp.to_numpy()
    #process sample file
    df_samp = pd.read_csv(sampfile,sep="\t")
    self.samp_columns = df_samp.columns()
    self.samp_data = df_samp.to_numpy()
  
  def get_assignments(self,k, fn):
    '''
    input:
    output:return list of integer 0 and 1. 
    '''
    self.k = k
    self.fn = fn
    
  def kNN(self,training_set, test):
    '''
    store distance in list, sort, and return closest self.k
    '''
    distance = []
    
  def distance(self,x,y):
    return np.linalg.norm(x,y)

  def calc_metrics(self,k, fn): 
    #
    return None


if __name__ == "__main__":
  exp_file = sys.argv[1]
  samp_file = sys.argv[2]
  knn=KNN()
  knn.load_data(exp_file,samp_file)