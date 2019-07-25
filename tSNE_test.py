import msprime
import numpy as np
import math
import os
import time
import re
import random
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA


start_time = time.time()


reps=1
for REPS in range(0,reps):

    
    N_locals=10000
    N_metropolis=10000
    
    generation_time = 20
    T_COLONIZATION=700/generation_time
    
    
    COLONIZER=random.randint(0,1)
    COLONIZER=1
    if COLONIZER==0:
        N_initial_colony=int(round(random.uniform(200.0,float(N_locals))))
        while N_initial_colony>N_metropolis:
            N_initial_colony=int(round(random.uniform(200.0,float(N_metropolis))))
    if COLONIZER==1:
        N_initial_colony=1000



    r_locals=0.000
    r_metropolis=0.000
    r_colony=0.000
    


    N_finale_colony=N_initial_colony / (math.exp(-r_colony * T_COLONIZATION))

    ###############################################################################################################################
    


    population_configurations = [
        msprime.PopulationConfiguration(initial_size=N_locals,growth_rate=r_locals),
        msprime.PopulationConfiguration(initial_size=N_metropolis, growth_rate=r_metropolis),
        msprime.PopulationConfiguration(initial_size=N_finale_colony, growth_rate=r_colony)
    ]



    migration_matrix = [
        [0,0.00001,0.0001],
        [0.000001,0,0.0001],
        [0.0001,0.0001,0]]

    N1=30
    N2=30
    N3=30
    POPS=[N1,N2,N3]
    samples=[msprime.Sample(0,0)]*N1 + [msprime.Sample(1,0)]*N2 + [msprime.Sample(2,0)] *N3

    demographic_events = [
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(0, 2)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(2, 0)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(1, 2)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(2, 1)),
    msprime.PopulationParametersChange(time=T_COLONIZATION, initial_size=N_initial_colony, growth_rate=0, population_id=2),
    msprime.MassMigration(time=T_COLONIZATION, source=2, destination=COLONIZER, proportion=1.0),
    msprime.MassMigration(time=1000/20, source=1, destination=0, proportion=1.0),
    
    ]

    






######################################################################################################################################################
#RUN the simulation and output genotypes in vcfs and ms format files, one for each chrom 
    variantinfo=open('variants_info.txt','w')
    variantinfo.write('CHROM\tVARIANT\tPOSITION\n')
    variants=[]
    for j in range(0,1):
    
        dd = msprime.simulate(samples=samples,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,mutation_rate=1e-8,
            demographic_events=demographic_events,length=10000000,recombination_rate=2e-8)

        #for the ms file
        
        outfile=open('ms_prime_{}'.format(j),'w')   
        for var in dd.variants():
            variants.append([var.position,var.index])
            variantinfo.write('{}\t{}\t{}\n'.format(j,var.index,var.position))
            for genotype in var.genotypes:
                outfile.write(str(genotype))
            outfile.write('\n')
        outfile.close()    
        #for the vcf file    
        wow=open('mynewvcf{}.vcf'.format(j),'w')
        dd.write_vcf(wow,2,str(j))
        wow.close()

        population_labels= ["locals"]*int(N1/2) + ["metropolis"]*int(N2/2) + ["colony"]*int(N3/2)
        d=0
        newlabels=[]
        for i in range(0,len(population_labels)):
            newlabels.append(population_labels[i]+str(d))
            d+=1
            if i==len(population_labels)-2:
                newlabels.append(population_labels[i]+str(d))
                break
            if population_labels[i]!=population_labels[i+1]:
                d=0
        population_labels=newlabels
        wow=open('mynewvcf{}.vcf'.format(j))
        wowzers=open('myvcf{}.vcf'.format(j),'w')
        for line in wow:
            line=line.strip().split()
            if line[0]=='#CHROM':
                line[9:]=population_labels
            wowzers.write("\t".join(line))
            wowzers.write("\n")
        wow.close()
        wowzers.close()

        vcf=open('mynewvcf{}.vcf'.format(j))
        
        population_labels=population_labels
        genotypes=[[] for x in range(0,len(population_labels))]
        
        for line in vcf:
            if line[0]!='#':
                line=line.strip().split()
                counter=0
                for j in line[9:]:
                    if j=='0|0':
                        genotypes[counter].append(0)
                    if j=='1|0' or j=='0|1':
                        genotypes[counter].append(1)
                    if j=='1|1':
                        genotypes[counter].append(2)
                    counter+=1
        
        #TRANSFORM TO tSNE
        from sklearn.manifold import TSNE
        X = np.asarray(genotypes)    
        X_embedded = TSNE(n_components=2,learning_rate=250.0,n_iter=2000,perplexity=15.0).fit_transform(X)    
        #print(X_embedded.shape)   
        #perplexity=> keep it near sample size
        
        
        
        
   
        #PLOTING
        plt.figure(figsize=(100, 60))
        colors=['red']*int(N1/2) + ['blue']*int(N2/2) + ['yellow']*int(N3/2)
           
        plt.scatter([x[0] for x in X_embedded],[x[1] for x in X_embedded],label=population_labels,c=colors)
        plt.show()
        
        
        pca = PCA(n_components=10).fit_transform(genotypes)
        print(len(pca))
        
        plt.figure(figsize=(100, 60))
        colors=['red']*int(N1/2) + ['blue']*int(N2/2) + ['yellow']*int(N3/2)
           
        plt.scatter([x[0] for x in pca],[x[1] for x in pca],label=population_labels,c=colors)
        plt.show()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        
            