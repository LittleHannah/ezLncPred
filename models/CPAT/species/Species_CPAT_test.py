import numpy as np
import pandas as pd
import os
import random
from Bio import SeqIO
from sklearn.metrics import roc_auc_score,accuracy_score,recall_score,f1_score,precision_score
from sklearn.ensemble import RandomForestClassifier

def model_performance(model_proba,model_predict,test_y):
    print(len(test_y))
    print(len(model_proba))
    AUC = roc_auc_score(test_y, model_proba)
    Accuracy = accuracy_score(test_y, model_predict)
    Sensitivity = recall_score(test_y, model_predict)
    Specificity=(Accuracy*len(test_y)-Sensitivity*sum(test_y))/(len(test_y)-sum(test_y))
    F1 = f1_score(test_y, model_predict)
    Precision = precision_score(test_y, model_predict)
    return AUC,Accuracy,Sensitivity,Specificity,F1,Precision

lncRNA = list(SeqIO.parse("/home/ls1/data/CPPred_data/S.cerevisiae/S.cerevisiae_ncrna.fa","fasta"))
pcts = list(SeqIO.parse("/home/ls1/data/CPPred_data/S.cerevisiae/S.cerevisiae_coding_RNA.fa","fasta"))

print(len(lncRNA))
print(len(pcts))
global_y = []
for num in range(len(lncRNA)):
        global_y.append(1)
for num in range(len(pcts)):
        global_y.append(0)
global_y = np.asarray(global_y)
all_X = np.concatenate((lncRNA, pcts), axis=0)
index = [i for i in range(len(all_X))] 
random.shuffle(index)
all_X = all_X[index]
SeqIO.write(all_X,"all_X_test","fasta")
test_file = "all_X_test"
global_y = global_y[index]
output_file = "CPAT_result"
species = 'Human'

print("CPAT is Running ... ")
os.system("python2 /home/ls1/CPAT-1.2.4/bin/cpat.py -g "+test_file+" -d "+species+".logit.RData -x "+species+"_Hexamer.tsv -o "+output_file)
result = pd.read_csv(output_file,sep="\t")
proba = list(1 - result['coding_prob'])
label = list(result['coding_prob'].apply(lambda x:1 if x<=0.5 else 0))


	
AUC,Accuracy,Sensitivity,Specificity,F1,Precision=model_performance(proba,label,global_y)
print("AUC:",AUC)
print("Accuracy:",Accuracy)
print("Sensitivity:",Sensitivity)
print("Specificity:",Specificity)
print("F1:",F1)
print("Precision:",Precision)


