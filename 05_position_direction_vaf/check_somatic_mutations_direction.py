#/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: calculat the mutation change direction when compare to the reference genome
# usage: python check_somatic_mutations_direction.py
# goal, for the snps in each sample we need calculat the mutation change type of each snps

# read the snps file of sample
# and for the snp we got the direction of the mutation
# A->T, A->C, A->G, T->A, T->C, T->G, C->A, C->T, C->G, G->A, G->T, G->C
# total 12 types of the mutation
import pandas as pd
import gzip
from matplotlib import rcParams

# Update the font properties
rcParams['font.family'] = 'Arial'


# read the vcfs of samples

def read_vcf(vcf):
    mutation = {}
    with gzip.open(vcf, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            elif line.startswith('Chr'):
                #print(line)
                line = line.strip().split('\t')
                #pos = line[0]+'_'+line[1]
                ref = line[3]
                alt = line[4]
                quality = line[6]
                if quality != 'REJECT' and len(ref) == 1 and len(alt) == 1:
                    if ref + '->' + alt not in mutation:
                        mutation[ref + '->' + alt] = 1
                    else:
                        mutation[ref+'->'+alt] = mutation.get(ref+'->'+alt, 0) + 1
    return mutation

# read the change header file
change_head = pd.read_csv("path/script/change_head.csv", sep="\t")

#vcf1 = 'path/somaticseq_data/somatic_1/indels.snvs_filter.vcf.gz'
#mutation1 = read_vcf(vcf1)
#print(mutation1)
total_mutation = {}

# plot the muatation direction of each sample 
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
fig.set_size_inches(6.5, 8)
ax.axis('off')
# the mutation direction order of types
# order = ['C_to_T', 'G_to_A','A_to_G', 'T_to_C', 'G_to_T', 'C_to_A', 'T_to_G', 'A_to_C', 'G_to_C', 'C_to_G', 'T_to_A', 'A_to_T']
order = ['C->T', 'G->A','A->G', 'T->C', 'G->T', 'C->A', 'T->G', 'A->C', 'G->C', 'C->G', 'T->A', 'A->T']
# in here C_to_T', 'G_to_A','A_to_G', 'T_to_C', 'G_to_T' are the transition mutation
# while others are the transversion mutation
for i in range(0,72):
    x =i//12
    y = i%12
    #box = plt.Rectangle((x*1.1, y*1.1), 1.1, 1.1, edgecolor="black", facecolor="white", linewidth=0.4)
    #ax.add_artist(box)
    box = plt.Rectangle((x*1.1+0.08, y*1.1+0.1), 0.004, 0.9, edgecolor="black", facecolor="black", linewidth=0.2)
    ax.add_artist(box)
    for k in range(0, 4):
        box = plt.Rectangle((x*1.1+0.07, y*1.1+0.1+k*100 *0.0025), 0.01, 0.001, edgecolor="black", facecolor="black", linewidth=0.2)
        ax.add_artist(box)
        plt.text(x*1.1+0.06-0.02, y*1.1+0.1+k*100 *0.0025+0.001, str(k*100), fontsize=2.5, color='black', ha='center', va='center',rotation=90)
    ax.add_artist(box)
    ax.set_xlim(0, 6.45)
    ax.set_ylim(0, 13.2)
    # plot the distribution of the snps and indel
    # read vcf file
    change_header = change_head[change_head["change"] == i+1]["old"].values[0]
    vcf = 'path/somaticseq_data/somatic_'+str(change_header)+'/indels.snvs_filter.vcf.gz'
    mutation = read_vcf(vcf)
    total_mutation[i+1] = mutation
    print(mutation)
    total_mutation_num = 0
    transition_num = 0
    transversion_num = 0
    
    for j in range(12):
        type = order[j]
        num = mutation.get(type, 0)
        if j <4 :
            box = plt.Rectangle((x*1.1+0.12+j*0.08, y*1.1+0.1), 0.04, num*0.0025, edgecolor=None, facecolor="#79BA83", linewidth=0.4, alpha=0.8)
            transition_num += num
        else:
            box = plt.Rectangle((x*1.1+0.12+j*0.08, y*1.1+0.1), 0.04, num*0.0025, edgecolor=None, facecolor="#79BA83", linewidth=0.4, alpha=0.4)
            transversion_num += num
        ax.add_artist(box)
        if num > 65:
            text = plt.text(x*1.1+0.12+j*0.08+0.02, y*1.1+0.1+num*0.0025 +0.05, str(num), fontsize=2.4, color='black', ha='center', va='center')
        else:
            text = plt.text(x*1.1+0.12+j*0.08+0.02, y*1.1+0.1+65*0.0025 +0.05, str(num), fontsize=2.4, color='black', ha='center', va='center')
        text = plt.text(x*1.1+0.12+j*0.08+0.03, y*1.1+0.1+0.1, type, fontsize=2.5, color='black', ha='center', va='center', rotation=90)
        total_mutation_num += num
    # add the mark of all mutation type
    box = plt.Rectangle((x*1.1+0.12, y*1.1+0.075), 0.28, 0.015, edgecolor=None, facecolor="#79BA83", linewidth=0.4, alpha=0.8)
    ax.add_artist(box)
    text = plt.text(x*1.1+0.26, y*1.1+0.02, 'Transition: '+str(format((transition_num/total_mutation_num), '.0%')), fontsize=2.5, color='black', ha='center', va='center')
    box = plt.Rectangle((x*1.1+0.44, y*1.1+0.075), 0.6, 0.015, edgecolor=None, facecolor="#79BA83", linewidth=0.4, alpha=0.4)
    ax.add_artist(box)
    text = plt.text(x*1.1+0.72, y*1.1+0.02, 'Transversion: '+str(format((transversion_num/total_mutation_num), '.0%')), fontsize=2.5, color='black', ha='center', va='center')
    # and the snp numbers 
    text = plt.text(x*1.1+0.5, y*1.1+0.8, 'sample '+str(i+1), fontsize=4)
    text = plt.text(x*1.1-0.02, y*1.1+0.6, 'SNP Numbers', fontsize=2.5, color='black', ha='center', va='center', rotation=90)
# plt.savefig('mutation_direction.png', dpi=600, bbox_inches='tight')
plt.savefig('mutation_direction.pdf', dpi=600, bbox_inches='tight')
# and we need plot a overall number of the snps

plt.close()

fig, ax = plt.subplots()
fig.set_size_inches(2.4, 2.4)
plt.axis('off')
ax.set_xlim(0.05, 1.05)
ax.set_ylim(0, 1.0)
#box = plt.Rectangle((0, 0), 1.1, 1.1, edgecolor="black", facecolor="white", linewidth=2)
#ax.add_artist(box)
order = ['C->T', 'G->A','A->G', 'T->C', 'G->T', 'C->A', 'T->G', 'A->C', 'G->C', 'C->G', 'T->A', 'A->T']
# in here C_to_T', 'G_to_A','A_to_G', 'T_to_C', 'G_to_T' are the transition mutation
# while others are the transversion mutation
total_mutation_num = 0
transition_num = 0
transversion_num = 0

for j in range(12):
    type = order[j]
    number_data = []
    for i in range(0, 72):
        number_data.append(total_mutation[i+1].get(type, 0))
    average = sum(number_data)/72
    sd = (sum([(x-average)**2 for x in number_data])/72)**0.5
    if j <4 :
        box = plt.Rectangle((0.12+j*0.08, 0.1), 0.04, average*0.0025, edgecolor=None, facecolor="#79BA83", linewidth=0.4, alpha=0.8)
        transition_num += average
    else:
        box = plt.Rectangle((0.12+j*0.08, 0.1), 0.04, average*0.0025, edgecolor=None, facecolor="#79BA83", linewidth=0.4, alpha=0.4)
        transversion_num += average
    ax.add_artist(box)
    # add the sd
    box = plt.Rectangle((0.12+j*0.08+0.019, 0.1+average*0.0025-sd*0.0025), 0.002, sd*0.0025*2, edgecolor=None, facecolor="black", linewidth=0.4, alpha=0.5)
    ax.add_artist(box)  
    box = plt.Rectangle((0.12+j*0.08+0.016, 0.1+average*0.0025-sd*0.0025), 0.008, 0.00002, edgecolor='black', facecolor="black", linewidth=1 , alpha=0.5)
    ax.add_artist(box)
    box = plt.Rectangle((0.12+j*0.08+0.016, 0.1+average*0.0025+sd*0.0025), 0.008, 0.00002, edgecolor='black', facecolor="black", linewidth=1 , alpha=0.5)
    ax.add_artist(box)
    #plt.text(0.12+j*0.08+0.02, 0.1+average*0.0025+sd*0.0025+0.03, str(round(average, 2)), fontsize=3, color='black', ha='center', va='center')
    
    if average > 35:
        text = plt.text(0.12+j*0.08+0.02, 0.1+(average+sd)*0.0025+0.03, str(int(round(average, 0))), fontsize=5, color='black', ha='center', va='center')
    else:
        text = plt.text(0.12+j*0.08+0.02, 0.1+(35+sd)*0.0025+0.03, str(int(round(average, 0))), fontsize=5, color='black', ha='center', va='center')
    text = plt.text(0.12+j*0.08+0.02, 0.1+0.05, type, fontsize=5, color='black', ha='center', va='center', rotation=90, )
    total_mutation_num += average
    
box = plt.Rectangle((0.10, 0.1), 0.003, 0.9, edgecolor="black", facecolor="black", linewidth=0.4, alpha=0.8)
ax.add_artist(box)
for k in range(4):
    box =plt.Rectangle((0.09, 0.1+k*100*0.0025), 0.01, 0.002, edgecolor="black", facecolor="black", linewidth=0.4, alpha=0.8)
    ax.add_artist(box)
    text =text = plt.text(0.08, 0.1+k*100*0.0025, str(k*100), fontsize=5, color='black', ha='right', va='center')
    
# add the mark of all mutation type
box = plt.Rectangle((0.12, 0.075), 0.28, 0.015, edgecolor=None, facecolor="#79BA83", linewidth=0.4, alpha=0.8)
plt.gca().add_artist(box)
text = plt.text(0.26, 0.04, 'Transition: '+str(format((transition_num/total_mutation_num), '.0%')), fontsize=6, color='black', ha='center', va='center')
box = plt.Rectangle((0.44, 0.075), 0.6, 0.015, edgecolor=None, facecolor="#79BA83", linewidth=0.4, alpha=0.4)
plt.gca().add_artist(box)
text = plt.text(0.76, 0.04, 'Transversion: '+str(format((transversion_num/total_mutation_num), '.0%')), fontsize=6, color='black', ha='center', va='center')
# and the snp numbers
#text = plt.text(0.3, 0.8, 'All samples', fontsize=12)
text = plt.text(0.15, 1.03, 'SNP number', fontsize=5, color='black', ha='center', va='center')
# ax.text(1.02, 1.02, "b", fontsize=40, rotation=0, verticalalignment='center', horizontalalignment='center', color="black", weight="bold", family="Arial")
#plt.savefig('mutation_direction_all.png', dpi=600, bbox_inches='tight')
plt.savefig('mutation_direction_all.pdf', dpi=800, bbox_inches='tight')




