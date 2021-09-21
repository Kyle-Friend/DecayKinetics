# -*- coding: utf-8 -*-
"""
Created by Yavuz Durmaz and Kyle Friend at Washington and Lee University in May 2019.
"""
import math
import pickle
import statistics
import numpy as np
from scipy.optimize.minpack import curve_fit

"""This starting file correlates the samples with the number of hours they
were chased. The samples are also group according to experimental procedure
in the block after that."""
SAM_Times = {}
temp = []
with open(r'E:/DecayKinetics/SAMfiles.txt','r') as infile:
    for line in infile:
        line = line.split()
        line[1] = int(line[1])
        SAM_Times[line[0]] = line[1]
        temp.append(line[0])

SAMfiles_clustered = [] # list of samples grouped by 6
for i in range(len(temp)):
    if i % 6 == 0:
        SAMfiles_clustered.append(temp[i:i + 6])

"""Loads the dictionary and extracts a non-redundant list of genes from the
mouse genome."""
Time_Point_data_Dict = pickle.load(open(r'E:/DecayKinetics/Kinet_Dict.p','rb'))
joined_sample_gene = list(Time_Point_data_Dict.keys())
allgenes = []
for a in joined_sample_gene:
    allgenes.append(a[0])
allgenes = list(set(allgenes))

outfile = open('mrnahalflife.txt', 'w')

"""The first piece of this code block extracts the relevant values for the
downstream half-life calculation. These are percent converted, total number
of reads and number of hours (chase). The second section then uses the data
to calculate half-lives."""
for gene in allgenes:
    loop_breaker = 1
    outfile.write(str(gene) + '\t')
    results = []
    for group in SAMfiles_clustered:
        clust_temp = []

        for j in group:
            temp = []
            val = Time_Point_data_Dict[(gene, j)]
            temp.append(val[2])
            temp.append(val[1])
            time = SAM_Times[j]
            temp.append(time)
            clust_temp.append(temp)

        times = []
        amounts = []
        read_count = []
        for comb in clust_temp:
            times.append(comb[2]*60)
            amounts.append(comb[0])
            read_count.append(comb[1])
        x = np.asarray(times)
        y = np.asarray(amounts)
        guess = [1, -0.05, 0.5]
        exp_decay = lambda x, A, t, y0: A*2**(x*t) + y0
        try:
            params, cov = curve_fit(exp_decay, x, y, p0=guess)
            A, t, y0 = params
            perr = np.sum(np.sqrt(np.diag(cov)))
        except RuntimeError:
            A, t, y0 = 0, 0, 0
            perr = 'NA'
        if A <= 0:
            A, t, y0 = 0, 0, 0
        else:
            t = -1/t

        read_count = statistics.mean(read_count)
        results.append(read_count)
        results.append(t)
        results.append(A)
        results.append(y0)
        results.append(perr)
    for j in results:
        outfile.write(str(j) + '\t')
    outfile.write('\n')
outfile.close()
