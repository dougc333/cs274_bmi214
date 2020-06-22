#!/usr/bin/env python
# coding: utf-8

# In[3]:


from sklearn.neighbors import KNeighborsClassifier
import pandas as pd
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import LeaveOneOut
import numpy as np

df = pd.read_csv("leukemia.csv")
#print(df.head())
#print(df.columns)
#X = df["Zyxin,"PRG1 Proteoglycan 1, secretory granule",CD33 CD33 antigen (differentiation antigen),DF D component of complement (adipsin),RNS2 Ribonuclease 2 (eosinophil-derived neurotoxin; EDN),CST3 Cystatin C (amyloid angiopathy and cerebral hemorrhage),APLP2 Amyloid beta (A4) precursor-like protein 2,"GLUTATHIONE S-TRANSFERASE, MICROSOMAL",CTSD Cathepsin D (lysosomal aspartyl protease),ATP6C Vacuolar H+ ATPase proton channel subunit,PROTEASOME IOTA CHAIN,LYN V-yes-1 Yamaguchi sarcoma viral related oncogene homolog,PPGB Protective protein for beta-galactosidase (galactosialidosis),HKR-T1,SPI1 Spleen focus forming virus (SFFV) proviral integration oncogene spi1,LYZ Lysozyme1,Lysozyme gene (EC 3.2.1.17),LYZ Lysozyme2,FAH Fumarylacetoacetate,RABAPTIN-5 protein,Leukotriene C4 synthase (LTC4S) gene,MYL1 Myosin light chain (alkali),ALDOA Aldolase A,"ARHG Ras homolog gene family, member G (rho G)",Liver mRNA for interferon-gamma inducing factor(IGIF),TCF3 Transcription factor 3 (E2A immunoglobulin enhancer binding factors E12/E47)1,MPO Myeloperoxidase,MB-1 gene,"FTL Ferritin, light polypeptide",Azurocidin gene,TIMP2 Tissue inhibitor of metalloproteinase 2,Nucleoside-diphosphate kinase,Macmarcks,ALDR1 Aldehyde reductase 1 (low Km aldose reductase),ME491  gene extracted from H.sapiens gene for Me491/CD63 antigen,GPX1 Glutathione peroxidase 1,LAMP2 Lysosome-associated membrane protein 2 {alternative products},"C-myb gene extracted from Human (c-myb) gene, complete primary cds, and five complete alternatively spliced cds","Quiescin (Q6) mRNA, partial cds",SELL Leukocyte adhesion protein beta subunit2,CCND3 Cyclin D3,TCF3 Transcription factor 3 (E2A immunoglobulin enhancer binding factors E12/E47)2,"KAI1 Kangai 1 (suppression of tumorigenicity 6, prostate; CD82 antigen (R2 leukocyte antigen, antigen detected by monoclonal and antibody IA4))","ITGAX Integrin, alpha X (antigen CD11C (p150), alpha polypeptide)","PFC Properdin P factor, complement",PROBABLE G PROTEIN-COUPLED RECEPTOR LCR1 HOMOLOG,Transcriptional activator hSNF2b1,INTERLEUKIN-8 PRECURSOR,Oncoprotein 18 (Op18) gene,"Clone 22 mRNA, alternative splice variant alpha-1",CYSTATIN A,14-3-3 PROTEIN TAU,Hunc18b2,PTH2 parathyroid hormone receptor mRNA,MSH2 DNA repair protein MSH2,MANB Mannosidase alpha-B (lysosomal),KIAA0022 gene,"SPTAN1 Spectrin, alpha, non-erythrocytic 1 (alpha-fodrin)",VIL2 Villin 2 (ezrin),ALDH7 Aldehyde dehydrogenase 7,"ZNF91 Zinc finger protein 91 (HPF7, HTF10)",Keratin 10 type I intermediate filament (KRT10) mRNA,"GB DEF = Topoisomerase type II (Topo II) mRNA, partial cds",INTERFERON GAMMA UP-REGULATED I-5111 PROTEIN PRECURSOR,CDC10 Cell division cycle 10 (homologous to CDC10 of S. cerevisiae,26-kDa cell surface protein TAPA-1 mRNA,TOP2B Topoisomerase (DNA) II beta (180kD),SMT3B protein,"Calcineurin A catalytic subunit [human, testis, mRNA, 2134 nt]","PRKCD Protein kinase C, delta",Inducible protein mRNA,"GTF2E2 General transcription factor TFIIE beta subunit, 34 kD",NADH:ubiquinone oxidoreductase subunit B13 (B13) mRNA,Thymopoietin beta mRNA,"TCF12 Transcription factor 12 (HTF4, helix-loop-helix transcription factors 4)",LEUKOCYTE ELASTASE INHIBITOR,GRN Granulin,CRYZ Crystallin zeta (quinone reductase),PERIPHERAL-TYPE BENZODIAZEPINE RECEPTOR,PCBD 6-pyruvoyl-tetrahydropterin synthase/dimerization cofactor of hepatocyte nuclear factor 1 alpha (TCF1),"Heat Shock Protein, 70 Kda (Gb:Y00371)",ERC-55 mRNA,Terminal transferase mRNA,PLECKSTRIN,"RHD Rhesus blood group, D antigen",Transcriptional activator hSNF2b2,"CHRNA7 Cholinergic receptor, nicotinic, alpha polypeptide 7","GB DEF = Retinoblastoma susceptibility protein (RB1) gene, with a 3 bp deletion in exon 22 (L11910 bases 161855-162161)",IL7R Interleukin 7 receptor,PHOSPHATIDYLINOSITOL,PGD Phosphogluconate dehydrogenase,GB DEF = GTP-binding protein (RAB3B) mRNA,CAB3b mRNA for calcium channel beta3 subunit,"CD36 CD36 antigen (collagen type I receptor, thrombospondin receptor)","KIAA0128 gene, partial cds",SELL Leukocyte adhesion protein beta subunit1,Nuclear Factor Nf-Il6,Phosphotyrosine independent ligand p62 for the Lck SH2 domain mRNA,AKT1 V-akt murine thymoma viral oncogene homolog 1,LRPAP1 Low density lipoprotein-related protein-associated protein 1 (alpha-2-macroglobulin receptor-associated protein 1,SNRPN Small nuclear ribonucleoprotein polypeptide N,Amyloid precursor protein-binding protein 1 mRNA,DCK Deoxycytidine kinase,ADPRT ADP-ribosyltransferase (NAD+; poly (ADP-ribose) polymerase),GB DEF = P85 beta subunit of phosphatidyl-inositol-3-kinase,Clone 23721 mRNA sequence,Thymopoietin (TMPO) gene,GB DEF = DNA for cellular retinol binding protein (CRBP) exons 3 and 4,"ACADM Acyl-Coenzyme A dehydrogenase, C-4 to C-12 straight chain",LMNA Lamin A,"LMP2 gene extracted from H.sapiens genes TAP1, TAP2, LMP2, LMP7 and DOB","SERINE/THREONINE PROTEIN PHOSPHATASE 2B CATALYTIC SUBUNIT, BETA ISOFORM",GB DEF = Homeodomain protein HoxA9 mRNA,CSF3R Colony stimulating factor 3 receptor (granulocyte),Hepatocyte growth factor-like protein gene,Tetracycline transporter-like protein mRNA,"Putative enterocyte differentiation promoting factor mRNA, partial cds","PI Protease inhibitor 1 (anti-elastase), alpha-1-antitrypsin",Epb72 gene exon 1,Interleukin 8 (IL8) gene,GB DEF = Myosin-IE,Activin type II receptor,Skeletal muscle LIM-protein SLIM1 mRNA,TALDO Transaldolase,Transmembrane protein,HMG1 High-mobility group (nonhistone chromosomal) proteinleukemia_type"]
y = df["leukemia_type"] 
X=df.loc[:, df.columns != 'leukemia_type']
y_np = y.to_numpy()
print(y_np.shape)
print(y_np)
y_all = np.where(y_np=="ALL", 1, y_np) 
y_all= np.where(y_all=="AML",0,y_all)
y_all=y_all.astype('int')
#print("y_all replace ALL with 1 and AML 0")
#print(y_all)
y_aml=np.where(y_np=="ALL", 0, y_np)
y_aml=np.where(y_aml=="AML",0,y_aml)
print("y_aml replace AML with 1 and ALL 0")
#print(y_aml)
knn = KNeighborsClassifier(n_neighbors=5)

loo = LeaveOneOut()
for train_index, test_index in loo.split(X):
    X_train, X_test = X.to_numpy()[train_index], X.to_numpy()[test_index]
    #y_train, y_test = y.to_numpy()[train_index], y.to_numpy()[test_index]
    y_train, y_test = y_all[train_index], y_all[test_index]
    print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
    probas_ = knn.fit(X_train, y_train).predict(X_test)
    # Compute ROC curve and area the curve
    #print(type(probas_),probas_)
    #print(y_test)
    fpr, tpr, thresholds = roc_curve(y_test, probas_,pos_label=1)
    print("fpr:",fpr,' tpr:',tpr)
    


# In[5]:


from sklearn.model_selection import LeaveOneOut
loo = LeaveOneOut()
#loo.get_n_splits(X,y)
print('====')
print(type(loo))
print('====')
print(X.shape,y.shape)
tp_ALL=0
tn_ALL=0
fp_ALL=0
fn_ALL=0
tp_AML=0
tn_AML=0
fp_AML=0
fn_AML=0

for train_index, test_index in loo.split(X):
    #print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X.to_numpy()[train_index], X.to_numpy()[test_index]
    y_train, y_test = y.to_numpy()[train_index], y.to_numpy()[test_index]
    #print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
    knn = KNeighborsClassifier(n_neighbors=5)
    knn.fit(X_train, y_train) 
    res = knn.predict(X_test)
    print(y.to_numpy()[test_index],res)
    if y.to_numpy()[test_index]=="ALL" and res=="ALL":
        tp_ALL +=1
        tn_AML+=1
    elif(y.to_numpy()[test_index]=="ALL" and res=="AML"):
        fn_ALL+=1
        fp_AML+=1
    elif(y.to_numpy()[test_index]=="AML" and res=="AML"):
        tn_ALL+=1
        tp_AML+=1
    elif(y.to_numpy()[test_index]=="AML" and res=="ALL"):
        fn_ALL+=1
        fn_AML+=1
    else:
        print('not see')

print("tp_ALL:",tp_ALL," tn_ALL:",tn_ALL," fp_ALL:",fp_ALL," fn_ALL:",fn_ALL, "tpr_all:",(tp_ALL)/(tp_ALL+fn_ALL),"fpr_all:",fp_ALL/(fp_ALL+tn_ALL))
print("tp_AML:",tp_AML," tn_AML:",tn_AML," fp_AML:",fp_AML," fn_AML:",fn_AML, "tpr_aml:",(tp_AML)/(tp_AML+fn_AML),"fpr_aml:",fp_AML/(fp_AML+tn_AML))


# In[6]:


from sklearn.naive_bayes import GaussianNB

from sklearn.model_selection import LeaveOneOut
loo = LeaveOneOut()
loo.get_n_splits(X,y)
print('====')
print(type(loo))
print('====')
print(X.shape,y.shape)
tp_ALL=0
tn_ALL=0
fp_ALL=0
fn_ALL=0
tp_AML=0
tn_AML=0
fp_AML=0
fn_AML=0

for train_index, test_index in loo.split(X):
    #print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X.to_numpy()[train_index], X.to_numpy()[test_index]
    y_train, y_test = y.to_numpy()[train_index], y.to_numpy()[test_index]
    #print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
    nb = GaussianNB()
    nb.fit(X_train, y_train) 
    res = nb.predict(X_test)
    print(y.to_numpy()[test_index],res)
    if y.to_numpy()[test_index]=="ALL" and res=="ALL":
        tp_ALL +=1
        tn_AML+=1
    elif(y.to_numpy()[test_index]=="ALL" and res=="AML"):
        fn_ALL+=1
        fp_AML+=1
    elif(y.to_numpy()[test_index]=="AML" and res=="AML"):
        tn_ALL+=1
        tp_AML+=1
    elif(y.to_numpy()[test_index]=="AML" and res=="ALL"):
        fn_ALL+=1
        fn_AML+=1
    else:
        print('not see')
print("tp_ALL:",tp_ALL," tn_ALL:",tn_ALL," fp_ALL:",fp_ALL," fn_ALL:",fn_ALL, "tpr_all:",(tp_ALL)/(tp_ALL+fn_ALL),"fpr_all:",fp_ALL/(fp_ALL+tn_ALL))
print("tp_AML:",tp_AML," tn_AML:",tn_AML," fp_AML:",fp_AML," fn_AML:",fn_AML, "tpr_aml:",(tp_AML)/(tp_AML+fn_AML),"fpr_aml:",fp_AML/(fp_AML+tn_AML))


# In[ ]:




