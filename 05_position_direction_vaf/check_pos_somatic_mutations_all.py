#!/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: check the mutation position in gene feature in four categories: Promoter, gene, Transposable elements, Intergenic region of all the 72 samples
# usage: python3 path/check_pos_somatic_mutations_all.py
import matplotlib.pyplot as plt
import numpy as np

path = 'path/somatic_information.txt'


Promoter = []
Gene = []
Transposon = []
mutation = []
TEs = []
IGR = []
Promoter_Gene_TE = []
Gene_TE = []
Promoter_Gene = []
Promoter_TE = []
TEprotein = []
with open(path, 'r') as f:
    k=0
    lines = f.readlines()
    #lines = lines.split('\n')
    sample_mutation = []
    for line in lines:
        k = k + 1
        if k%31 == 7:
            Gene.append(line.strip())
        if k%31 == 11:
            Promoter.append(line.strip())
        if k%31 == 14:
            Transposon.append(line.strip())
        if k%31 == 15:
            mutation.append(line.strip())
        if k%31 == 18:
            Gene_TE.append(line.strip())
        if k%31 == 19:
            Promoter_TE.append(line.strip())
        if k%31 == 20:
            Promoter_Gene.append(line.strip())
        if k%31 == 21:
            Promoter_Gene_TE.append(line.strip())
        if k%31 == 22:
            TEprotein.append(line.strip())
#IGR = list(set(mutation) - set(Promoter) - set(Gene) - set(TEs) + set(Gene_TE) + set(Promoter_TE) + set(Promoter_Gene) + set(Promoter_Gene_TE) +set(Promoter_Gene_TE))
#IGR = mutation-Promoter-Gene-TEs+Gene_TE+Promoter_TE+Promoter_Gene+Promoter_Gene_TE*2
for i in range(len(mutation)):
    TEs.append(int(Transposon[i])+int(TEprotein[i]))
    #igr = mutation[i] - Promoter[i] - Gene[i] - TEs[i] + Gene_TE[i] + Promoter_TE[i] + Promoter_Gene[i] + Promoter_Gene_TE[i]+Promoter_Gene_TE[i]
    igr = int(mutation[i]) - int(Promoter[i]) - int(Gene[i]) - int(TEs[i]) + int(Gene_TE[i]) + int(Promoter_TE[i]) + int(Promoter_Gene[i]) - int(Promoter_Gene_TE[i])
    IGR.append(igr)


# we get the 72 samples' mutation information
# next we need to plot the mutation information

fig, ax = plt.subplots()
# the size of the figure
fig.set_size_inches(2.4, 2.4)

mean_promoter = 0
std_promoter = 0
for i in range(len(Promoter)):
    mean_promoter = mean_promoter + int(Promoter[i])
mean_promoter = mean_promoter/len(Promoter)
for i in range(len(Promoter)):
    std_promoter = std_promoter + (int(Promoter[i])-mean_promoter)**2
std_promoter = np.sqrt(std_promoter/len(Promoter))
ax.bar(1, mean_promoter, yerr=std_promoter, align='center', alpha=0.5, ecolor='black', capsize=2, label='Promoter',fc='#cccc00',width=0.4)
box = plt.Rectangle((0.8, -20), 0.9, 10, color='#cccc00',alpha=0.5)
ax.add_patch(box)
ax.text(1, mean_promoter+std_promoter+80, str(format(mean_promoter, '.1f')), ha='center', va='center', fontsize=6)
mean_gene = 0
std_gene = 0
for i in range(len(Gene)):
    mean_gene = mean_gene + int(Gene[i])
mean_gene = mean_gene/len(Gene)
for i in range(len(Gene)):
    std_gene = std_gene + (int(Gene[i])-mean_gene)**2
std_gene = np.sqrt(std_gene/len(Gene))
ax.bar(2, mean_gene, yerr=std_gene, align='center', alpha=0.5, ecolor='black', capsize=2, label='Gene',fc='green',width=0.4)
box = plt.Rectangle((1.8, -20), 0.9, 10, color='green',alpha=0.5)
ax.add_patch(box)
ax.text(2, mean_gene+std_gene+80, str(format(mean_gene, '.1f')), ha='center', va='center', fontsize=6)
mean_tes = 0
std_tes = 0
for i in range(len(TEs)):
    mean_tes = mean_tes + int(TEs[i])
mean_tes = mean_tes/len(TEs)
for i in range(len(TEs)):
    std_tes = std_tes + (int(TEs[i])-mean_tes)**2
std_tes = np.sqrt(std_tes/len(TEs))
ax.bar(3, mean_tes, yerr=std_tes, align='center', alpha=0.5, ecolor='black', capsize=2, label='TEs',fc='red',width=0.4)
box = plt.Rectangle((2.8, -20), 0.9, 10, color='red',alpha=0.5)
ax.add_patch(box)
ax.text(3, mean_tes+std_tes+80, str(format(mean_tes, '.1f')), ha='center', va='center', fontsize=6)
mean_igr = 0
std_igr = 0
for i in range(len(IGR)):
    mean_igr = mean_igr + int(IGR[i])
mean_igr = mean_igr/len(IGR)
for i in range(len(IGR)):
    std_igr = std_igr + (int(IGR[i])-mean_igr)**2
std_igr = np.sqrt(std_igr/len(IGR))
ax.bar(4, mean_igr, yerr=std_igr, align='center', alpha=0.5, ecolor='black', capsize=2, label='IGR',fc='gray',width=0.4)
box = plt.Rectangle((3.8, -20), 0.9, 10, color='gray',alpha=0.5)
ax.add_patch(box)
ax.text(4, mean_igr+std_igr+80, str(format(mean_igr, '.1f')), ha='center', va='center', fontsize=6) 
#
ax.axis('off')

#total mutation

mean_mutation = 0
std_mutation = 0
for i in range(len(mutation)):
    mean_mutation = mean_mutation + int(mutation[i])
mean_mutation = mean_mutation/len(mutation)
for i in range(len(mutation)):
    std_mutation = std_mutation + (int(mutation[i])-mean_mutation)**2
std_mutation = np.sqrt(std_mutation/len(mutation))

ax.text(1, 100, str(format(mean_promoter/mean_mutation*100, '.0f'))+'%', ha='center', va='center', fontsize=5, rotation=90)
ax.text(2, 100, str(format(mean_gene/mean_mutation*100, '.0f'))+'%', ha='center', va='center', fontsize=5,rotation=90)
ax.text(3, 100, str(format(mean_tes/mean_mutation*100, '.0f'))+'%', ha='center', va='center', fontsize=5,rotation=90)
ax.text(4, 100, str(format(mean_igr/mean_mutation*100, '.0f'))+'%', ha='center', va='center', fontsize=5,rotation=90)
# add the information of the genome
genome_size = 216955151
gene = 99495975
promoter = 38963587
TEs = 73155645
gene_TE = 12776936
promoter_TE = 5524050
promoter_gene = 10782633
promoter_gene_TE = 350126
IGRS = genome_size - gene - promoter - TEs + gene_TE + promoter_TE + promoter_gene - promoter_gene_TE
ax.bar(1.5, promoter/genome_size*mean_mutation, align='center', alpha=0.3, ecolor='black', capsize=12, label='Promoter',fc='#cccc00',width=0.4)
ax.text(1.5,100, str(format(promoter/genome_size*100, '.0f'))+'%', ha='center', va='center', fontsize=5, rotation=90)
ax.bar(2.5, gene/genome_size*mean_mutation,align='center', alpha=0.3, ecolor='black', capsize=12, label='Gene',fc='green',width=0.4)
ax.text(2.5,100, str(format(gene/genome_size*100, '.0f'))+'%', ha='center', va='center', fontsize=5, rotation=90)
ax.bar(3.5, TEs/genome_size*mean_mutation,align='center', alpha=0.3, ecolor='black', capsize=12, label='TEs',fc='red',width=0.4)
ax.text(3.5,100, str(format(TEs/genome_size*100, '.0f'))+'%', ha='center', va='center', fontsize=5, rotation=90)
ax.bar(4.5, IGRS/genome_size*mean_mutation, align='center', alpha=0.3, ecolor='black', capsize=12, label='IGR',fc='gray',width=0.4)
ax.text(4.5,100, str(format(IGRS/genome_size*100, '.0f'))+'%', ha='center', va='center', fontsize=5, rotation=90)


# add the text

ax.text(1.25, -80, 'Promoter', ha='center', va='center', fontsize=6)
ax.text(2.25, -80, 'Gene', ha='center', va='center', fontsize=6)
ax.text(3.25, -80, 'TEs', ha='center', va='center', fontsize=6)
ax.text(4.25, -80, 'IGR', ha='center', va='center', fontsize=6)

ax.add_patch(plt.Rectangle((0.44, 0), 0.02, 2050,fc='black'))
for i in range(5):
    ax.add_patch(plt.Rectangle((0.34, i*500), 0.1, 5,fc='black'))
    ax.text(0.32, i*500, str(i*500), ha='right', va='center', fontsize=5,rotation=90)



ax.set_ylim(-150, 2100)
ax.set_xlim(-0.2, 4.9)
ax.text(-0.1, 1000, 'Number of mutations', ha='center', va='center', fontsize=6, rotation=90)
#box = plt.Rectangle((-0.39, -150), 5.30, 2250, fill=None, edgecolor='black', linewidth=1)
#ax.add_patch(box)




#ax.text(-0.2, 2000, "a", fontsize=34, rotation=0, verticalalignment='center', horizontalalignment='center', color="black", weight="bold", family="Arial")

plt.savefig('mutation_information.pdf', dpi=300, bbox_inches='tight')    