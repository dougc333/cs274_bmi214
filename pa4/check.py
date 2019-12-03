import pandas as pd
import numpy as np


df = pd.read_csv("outputfile2.csv",header=None).to_numpy()
#print(df.head())
print(df.shape)
for i in range(0,df.shape[0]):
    if df[i,:][0]==df[i,:][1]:
        print("******************")
'''
for row in df.iterrows():
    print("row:",row)
    #print(row[0],"**",row[1])
'''