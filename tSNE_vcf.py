import msprime
import argparse
import numpy as np
import math
import os
import time
import re
import random
import numba
import umap

from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
from bokeh.io import output_notebook, show
from bokeh.plotting import figure
from bokeh.models import TapTool, CustomJS, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.io import output_file
from bokeh import colors




parser = argparse.ArgumentParser()  
parser.add_argument('--vcf',nargs='+',type=str)
args = parser.parse_args()

vcf=open(args.vcf[0])





def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
    
    return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]

def get_colors(n):
    colorz=[]
    for k in range(0,n):
        c1=round(random.random(), 3)
        c2=round(random.random(), 3)
        c3=round(random.random(), 3)
        colorz.append((c1,c2,c3))
    return colorz
#####IIIINCOOOMPLEEEETE




for line in vcf:
    if line[0]=='#' and line[1]!='#':
        population_labels=line.strip().split()[9:]
        true_labels=[]
        for k in population_labels:
            familyname=re.search(r'([A-Z]+)_[A-Z]+',k)
            true_labels.append(familyname.group(1))
        genotypes=[[] for x in range(0,len(population_labels))]
    if line[0]!='#':
        line=line.strip().split()
        counter=0
        for j in line[9:]:
            if j=='0/0':
                genotypes[counter].append(0)
            if j=='1/0' or j=='0/1':
                genotypes[counter].append(1)
            if j=='1/1':
                genotypes[counter].append(2)
            if j=='./.':
                genotypes[counter].append(9)
            counter+=1

print(len(genotypes),len(genotypes[0]))     




####################################################################################################################################################
####################################################################################################################################################
       
#TRANSFORM TO tSNE
X = np.asarray(genotypes)  
pca_for_tSNE = PCA(n_components=20).fit_transform(genotypes)

WEIGHTS=(PCA(n_components=20).fit(genotypes).explained_variance_)
print(WEIGHTS)


def weighted_dist(a,b):

    distance = math.sqrt(sum([((a[x] - b[x])*WEIGHTS[x]) ** 2 for x in range(0,len(a))]))
    return distance
    
print(weighted_dist(pca_for_tSNE[0],pca_for_tSNE[50]))    
print(weighted_dist(pca_for_tSNE[0],pca_for_tSNE[1]))    



X_embedded = TSNE(verbose=0,n_components=2,learning_rate=200.0,n_iter=1000,perplexity=10.0).fit_transform(pca_for_tSNE) 
#print(X_embedded.shape)   





#PLOTING
plt.figure(figsize=(100, 60))

COLORPALLETE=get_colors(len(set(true_labels)))
COLORZ_TO_LABELS={}

uniquelabels=[x for x in set(true_labels)]
for j in range(0,len(uniquelabels)):
    COLORZ_TO_LABELS[uniquelabels[j]]=COLORPALLETE[j]

colors=[ COLORZ_TO_LABELS[x] for x in true_labels]
   
plt.scatter([x[0] for x in X_embedded],[x[1] for x in X_embedded],label=true_labels,c=colors)
plt.show()



X_embedded = TSNE(verbose=0,n_components=2,learning_rate=200.0,n_iter=1000,perplexity=10.0).fit_transform(pca_for_tSNE) 
#print(X_embedded.shape)   




####################################################################################################################################################
####################################################################################################################################################
# Second version of tSNE, weighted distance

X = np.asarray(genotypes)  
pca_for_tSNE = PCA(n_components=20).fit_transform(genotypes)

X_embedded = TSNE(verbose=0,n_components=2,learning_rate=200.0,n_iter=1000,perplexity=10.0,metric=weighted_dist()).fit_transform(pca_for_tSNE) 
#print(X_embedded.shape)   







#PLOTING
plt.figure(figsize=(100, 60))

COLORPALLETE=get_colors(len(set(true_labels)))
COLORZ_TO_LABELS={}

uniquelabels=[x for x in set(true_labels)]
for j in range(0,len(uniquelabels)):
    COLORZ_TO_LABELS[uniquelabels[j]]=COLORPALLETE[j]

colors=[ COLORZ_TO_LABELS[x] for x in true_labels]
   
plt.scatter([x[0] for x in X_embedded],[x[1] for x in X_embedded],label=true_labels,c=colors)
plt.show()




####################################################################################################################################################
####################################################################################################################################################
# Normal PCA





pca = PCA(n_components=2).fit_transform(genotypes)


plt.figure(figsize=(100, 60))

   
plt.scatter([x[0] for x in pca],[x[1] for x in pca],label=population_labels,c=colors)
plt.show()
        
        