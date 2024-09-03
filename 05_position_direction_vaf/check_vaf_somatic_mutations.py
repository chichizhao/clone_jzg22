#!/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: plot the somatic mutations vaf distribution of each sample
# usage: python check_vaf_somatic_mutations.py

import pandas as pd

change_head = pd.read_csv("path/script/change_head.csv", sep="\t")
print(change_head)

# read the af distribution of each sample
snps = {}
snps_dp = {}
indels = {}
indels_dp = {}

with open("path/mutation_site_count3.csv") as f:
    for line in f:
        if line.startswith("snp_vaf"):
            orgin_head = line.strip().split(",")[0].split("_")[2]
            change_header = change_head[change_head["old"] == orgin_head]["change"].values[0]
            # add the value to the snps
            snps[change_header] = line.strip().split(",")[1:]
        elif line.startswith("indel_vaf"):
            orgin_head = line.strip().split(",")[0].split("_")[2]
            #print(orgin_head)
            change_header = change_head[change_head["old"] == orgin_head]["change"].values[0]
            # add the value to the indels
            indels[change_header] = line.strip().split(",")[1:]
        elif line.startswith("snp_dp"):
            orgin_head = line.strip().split(",")[0].split("_")[2]
            change_header = change_head[change_head["old"] == orgin_head]["change"].values[0]
            # add the value to the snps_dp
            snps_dp[change_header] = line.strip().split(",")[1:]
        elif line.startswith("indel_dp"):
            orgin_head = line.strip().split(",")[0].split("_")[2]
            change_header = change_head[change_head["old"] == orgin_head]["change"].values[0]
            # add the value to the indels_dp
            indels_dp[change_header] = line.strip().split(",")[1:]
    
# plot the af distribution of each sample
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
fig.set_size_inches(6.5, 2)
ax.axis('off')
for i in range(0,3):
    x =i
    y =0
    box = plt.Rectangle((x*1.1+0.12, y*1.1+0.1), 1, .92, edgecolor="none", facecolor="white", linewidth=0.4)
    ax.add_artist(box)
    box = plt.Rectangle((x*1.1+0.12, y*1.1+0.1), 1, 0.001, edgecolor="black", facecolor="black", linewidth=0.4)
    ax.add_artist(box)
    box = plt.Rectangle((x*1.1+0.12, y*1.1+0.1), 0.001, .92, edgecolor="black", facecolor="black", linewidth=0.4)
    ax.add_artist(box)
    ax.set_xlim(0.1, 3.35)
    ax.set_ylim(0, 1.05)
    # plot the distribution of the snps and indels
    data_snps = np.array(snps[i+1], dtype=float)
    data_indels = np.array(indels[i+1], dtype=float)
    data = np.concatenate((data_snps, data_indels))
    data_snps_dp = np.array(snps_dp[i+1], dtype=float)
    data_indels_dp = np.array(indels_dp[i+1], dtype=float)
    data_dp = np.concatenate((data_snps_dp, data_indels_dp))
    
    # plot the data into the box
    counts = [0]*501
    dp_counts_10 = [0]*501
    dp_counts_15 = [0]*501
    dp_counts_20 = [0]*501
    dp_counts_25 = [0]*501
    number = 0
    for j in data:
        number += 1
        if j * 500 % 1 == 0:
            counts[int(j*500)] += 1
            if data_dp[number-1] < 10:
                dp_counts_10[int(j*500)] += 1
            elif data_dp[number-1] < 15:
                dp_counts_15[int(j*500)] += 1
            elif data_dp[number-1] < 20:
                dp_counts_20[int(j*500)] += 1 
            elif data_dp[number-1] < 25:
                dp_counts_25[int(j*500)] += 1   
        else:
            counts[int(j*500)+1] += 1
            if data_dp[number-1] < 10:
                dp_counts_10[int(j*500)+1] += 1
            elif data_dp[number-1] < 15:
                dp_counts_15[int(j*500)+1] += 1
            elif data_dp[number-1] < 20:
                dp_counts_20[int(j*500)+1] += 1
            elif data_dp[number-1] < 25:
                dp_counts_25[int(j*500)+1] += 1
                
    for j in range(501):
        if j < 500:
            box = plt.Rectangle((x*1.1+0.115+j*0.002, y*1.1+0.1), 1/500, counts[j]/100*0.9, edgecolor=None, facecolor="#447244", alpha=0.5)
            ax.add_artist(box)
            # add the dp less than 10
            box = plt.Rectangle((x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9), 1/500, -dp_counts_15[j]/100*0.9, edgecolor=None, facecolor="red", alpha=0.5)
            ax.add_artist(box)
        # and mark the sample that cross all the 72 samples 209
        if j == 500:
            box = plt.Rectangle((x*1.1+0.115+j*0.002, y*1.1+0.1), 1/500, (counts[j]-209)/300*0.9, edgecolor=None, facecolor="#447244", alpha=0.5)
            ax.add_artist(box)
            box = plt.Rectangle((x*1.1+0.115+j*0.002, y*1.1+0.1+(counts[j]-209)/300*0.9), 1/500, (-dp_counts_15[j])/300*0.9, edgecolor=None, facecolor="red", alpha=0.5)
            ax.add_artist(box)
            ax.text(x*1.1+0.115+j*0.002, y*1.1+0.1+(counts[j]-209)/300*0.9+0.09, str(counts[j]-209), fontsize=7, ha='right', va='center')
            arrow = plt.Arrow(x*1.1+0.115+j*0.002, y*1.1+0.1+(counts[j]-209)/300*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.4)
            plt.gca().add_patch(arrow)
            #box = plt.Rectangle((x*1.1+0.11+j*0.002, y*1.1+0.1+counts[j]/500*0.9-dp_counts_15[j]/500*0.9), 1/500, -209/500*0.9, edgecolor=None, facecolor="blue", alpha=0.5)
            #ax.add_artist(box)
        if j == 83:
            ax.text(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.09, str("1/6"), fontsize=7, ha='center', va='center')
            arrow = plt.Arrow(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            plt.gca().add_patch(arrow)

        if j == 167:
            ax.text(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.09, str("2/6"), fontsize=7, ha='center', va='center')
            arrow = plt.Arrow(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            plt.gca().add_patch(arrow)

        if j == 250:
            ax.text(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.09, str("3/6"), fontsize=7, ha='center', va='center')
            arrow = plt.Arrow(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            plt.gca().add_patch(arrow)

        if j == 334:
            ax.text(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.09, str("4/6"), fontsize=7, ha='center', va='center')
            arrow = plt.Arrow(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            plt.gca().add_patch(arrow)

        if j == 417:
            ax.text(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.09, str("5/6"), fontsize=7, ha='center', va='center')
            arrow = plt.Arrow(x*1.1+0.115+j*0.002, y*1.1+0.1+counts[j]/100*0.9+0.06, 0, -0.05, width=0.01, edgecolor="black", facecolor="black", linewidth=0.1)
            plt.gca().add_patch(arrow)

    ax.text(x*1.1+0.52, y*1.1+0.8, 'sample '+str(i+1), fontsize=7)
    
    # add the x and y label
    for k in range(1, 6):
        box = plt.Rectangle((x*1.1+0.12, y*1.1+0.1+k*0.1*1.8), 0.01, 0.001, edgecolor="black", facecolor="black", linewidth=0.4)
        ax.add_artist(box)
        ax.text(x*1.1+0.132, y*1.1+0.08+ k*0.1*1.8, str(k*20), fontsize=5)

    ax.text(x*1.1+0.08, y*1.1+0.55, "Alleles Number", fontsize=7, rotation=90, verticalalignment='center', horizontalalignment='center')
    for k in range(0, 6):
        box = plt.Rectangle((x*1.1+0.12+k*0.2, y*1.1+0.09), 0.001, 0.01, edgecolor="black", facecolor="black", linewidth=0.4)
        ax.add_artist(box)
        ax.text(x*1.1+0.1+k*0.2, y*1.1+0.05, str(k/5), fontsize=5)
   
    ax.text(x*1.1+0.65, y*1.1, "Alleles Frequency", fontsize=7, rotation=0, verticalalignment='center', horizontalalignment='center')   
    
plt.savefig("afs4.pdf", dpi=600, bbox_inches='tight')