#quizcode

import random
import numpy as np
from pvalue import Pvalue
import pickle

seeds = list(np.ndarray.flatten(np.random.randint(1000, size=(10,1))))

print("num seeds:",len(seeds))
p_100 = []
p_500=[]
p_1000=[]
for s in seeds:
    pv1 = Pvalue(100, s, "data/drugs.csv", "data/targets.csv", "P54577", "Q7RTX0").p()
    print("pv1:",pv1)
    p_100.append(pv1)

for s in seeds:
    pv2 = Pvalue(500, s, "data/drugs.csv", "data/targets.csv", "P54577", "Q7RTX0").p()
    print("pv2:",pv2)
    p_500.append(pv2)

for s in seeds:
    pv3 = Pvalue(1000, s, "data/drugs.csv", "data/targets.csv", "P54577", "Q7RTX0").p()
    print("pv3:",pv3)
    p_1000.append(pv3)

print(len(p_100))
print(len(p_500))
print(len(p_1000))

with open("p100",'wb') as fh100:
    pickle.dump(p_100,fh100)
with open("p500",'wb') as fh500:
    pickle.dump(p_500,fh500)
with open("p1000",'wb') as fh1000:
    pickle.dump(p_1000,fh1000)


print("mean 100:",np.mean(p_100))
print("mean 500:",np.mean(p_500))
print("mean 1000:",np.mean(p_1000))

print("std dev 100:",np.std(p_100))
print("std dev 500:",np.std(p_500))
print("std dev 1000:",np.std(p_1000))

