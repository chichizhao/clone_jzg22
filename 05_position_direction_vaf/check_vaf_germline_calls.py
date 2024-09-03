#!/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: based on the vcf in the gatk pipeline, plot the afs in each sample
# usage: python check_vaf_germline_calls.py

# vcf file path/jzg_22_old/2022_snp.vcf

import pandas as pd

change_head = pd.read_csv("path/script/change_head.csv", sep="\t")
print(change_head)

afs ={}
old_sample = []
with open("path/jzg_22_old/2022_snp.vcf") as f:
    for line in f:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            line = line.strip().split("\t")
            for i in range(9, len(line)):
                if line[i] not in afs:
                    old_sample.append(line[i])
                    afs[line[i]] = []
        else:
            line = line.strip().split("\t")
            for i in range(9, len(line)):
                AD = line[i].split(":")[1].split(",")
                DP = line[i].split(":")[2]
                if DP != "." and AD[0] != "." and AD[1] != "." and int(DP) > 0:
                    vaf = float(AD[1])/float(DP)
                    afs[old_sample[i-9]].append(vaf)
#print(afs)

import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
fig.set_size_inches(6.5,2)
ax.axis('off')
for i in range(0,3):
    x =i
    y = 0
    box = plt.Rectangle((x*1.1+0.12, y*1.1+0.1), 1, .92, edgecolor="none", facecolor="white", linewidth=1)
    ax.add_artist(box)
    box = plt.Rectangle((x*1.1+0.12, y*1.1+0.1), 1, 0.001, edgecolor="black", facecolor="black", linewidth=1)
    ax.add_artist(box)
    box = plt.Rectangle((x*1.1+0.12, y*1.1+0.1), 0.001, .92, edgecolor="black", facecolor="black", linewidth=1)
    ax.add_artist(box)
    ax.set_xlim(0.1, 3.35)
    ax.set_ylim(0, 1.05)
    # plot the distribution of the snps and indels
    sample_name = change_head[change_head["change"] == i+1]["old"].values[0]
    data = np.array(afs[sample_name], dtype=float)
    # remove the data equal to 0
    data = data[data != 0]   
    # plot the data into the box
    counts = [0]*101
    for j in data:
        if j * 100 % 1 == 0:
            counts[int(j*100)] += 1
        else:
            counts[int(j*100)+1] += 1
    for j in range(101):
        if j < 100:
            box = plt.Rectangle((x*1.1+0.1+j*0.01, y*1.1+0.1), 1/100, counts[j]/400*0.9, edgecolor=None, facecolor="#447244", alpha=0.5)
            ax.add_artist(box)
        if j == 17:
            ax.text(x*1.1+0.1+j*0.01, y*1.1+0.1+counts[j]/400*0.9+0.09, str("1/6"), fontsize=7, verticalalignment='center', horizontalalignment='center')
            arrow = plt.Arrow(x*1.1+0.1+j*0.01+1/400, y*1.1+0.1+counts[j]/400*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            ax.add_artist(arrow)
        if j == 34:
            ax.text(x*1.1+0.1+j*0.01, y*1.1+0.1+counts[j]/400*0.9+0.09, str("2/6"), fontsize=7, verticalalignment='center', horizontalalignment='center')
            arrow = plt.Arrow(x*1.1+0.1+j*0.01+1/400, y*1.1+0.1+counts[j]/400*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            ax.add_artist(arrow)
        if j == 50:
            ax.text(x*1.1+0.1+j*0.01, y*1.1+0.1+counts[j]/400*0.9+0.09, str("3/6"), fontsize=7, verticalalignment='center', horizontalalignment='center')
            arrow = plt.Arrow(x*1.1+0.1+j*0.01+1/400, y*1.1+0.1+counts[j]/400*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            ax.add_artist(arrow)
        if j == 67:
            ax.text(x*1.1+0.1+j*0.01, y*1.1+0.1+counts[j]/400*0.9+0.09, str("4/6"), fontsize=7, verticalalignment='center', horizontalalignment='center')
            arrow = plt.Arrow(x*1.1+0.1+j*0.01+1/400, y*1.1+0.1+counts[j]/400*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            ax.add_artist(arrow)
        if j == 84:
            ax.text(x*1.1+0.1+j*0.01, y*1.1+0.1+counts[j]/400*0.9+0.09, str("5/6"), fontsize=7, verticalalignment='center', horizontalalignment='center')
            arrow = plt.Arrow(x*1.1+0.1+j*0.01+1/400, y*1.1+0.1+counts[j]/400*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            ax.add_artist(arrow)
        if j == 100:
            box = plt.Rectangle((x*1.1+0.1+j*0.01, y*1.1+0.1), 1/100, counts[j]/3000*0.9, edgecolor=None, facecolor="#447244", alpha=0.5)
            ax.add_artist(box)
            ax.text(x*1.1+0.1+j*0.01, y*1.1+0.1+counts[j]/3000*0.9+0.09, str(counts[j]), fontsize=7, verticalalignment='center', horizontalalignment='right')
            arrow = plt.Arrow(x*1.1+0.1+j*0.01+1/400, y*1.1+0.1+counts[j]/3000*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            ax.add_artist(arrow)
            
            
            
    # add the sample name
    ax.text(x*1.1+0.55, y*1.1+0.85, 'sample '+str(i+1), fontsize=7)
    
    # add the x and y label
    for k in range(1, 6):
        box = plt.Rectangle((x*1.1+0.12, y*1.1+0.1+k*0.1*1.8), 0.01, 0.001, edgecolor="black", facecolor="black", linewidth=0.4)
        ax.add_artist(box)
        ax.text(x*1.1+0.132, y*1.1+0.08+ k*0.1*1.8, str(k*80), fontsize=5)
 
    ax.text(x*1.1+0.08, y*1.1+0.55, "Alleles Number", fontsize=7, rotation=90, verticalalignment='center', horizontalalignment='center')
    for k in range(0, 6):
        box = plt.Rectangle((x*1.1+0.12+k*0.2, y*1.1+0.09), 0.001, 0.01, edgecolor="black", facecolor="black", linewidth=0.4)
        ax.add_artist(box)
        ax.text(x*1.1+0.1+k*0.2, y*1.1+0.05, str(k/5), fontsize=5)
 
    ax.text(x*1.1+0.65, y*1.1, "Alleles Frequency", fontsize=7, rotation=0, verticalalignment='center', horizontalalignment='center')   


plt.savefig("afs_gatk3.pdf", dpi=600, bbox_inches='tight')