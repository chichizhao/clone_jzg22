# !/py3
#-*- coding:utf-8 -*-
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: based on the genome and the Sanger sequencing results, plot the 1/1 mutation sites and 1/0 mutation sites
# files 1. align fasta file 2. Sanger sequencing results file(ab1)
# usage: python plot_alignment.py -align chr3.aln -sanger jzg_Chr03_*.ab1_path -output chr3
import sys
import os
import matplotlib.pyplot as plt
import argparse
from Bio import SeqIO
from collections import defaultdict 
from statistics import mean
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="plot the 1/1 mutation sites and 1/0 mutation sites")
    parser.add_argument('-align', '--align', type=str, help='align fasta file')
    parser.add_argument('-sanger', '--sanger', type=str, help='sanger sequencing results file(ab1)')
    parser.add_argument('-o', '--output', type=str, help='output file name')
    args = parser.parse_args()
    return args
    
def read_align(align):
    align_dict = {}
    with open(align, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            align_dict[record.id] = str(record.seq)
    seq_id = list(align_dict.keys())[0] 
    align_start = 0
    align_end = 0
    for i in range(len(align_dict[seq_id])):
        if align_dict[seq_id][i] != '-':
            test = 0
            for key in align_dict.keys():
                if align_dict[key][i] == '-':
                    test = 1
                    break
            if test == 0:
                for j in range(i,i+35):
                    for key in align_dict.keys():
                        if align_dict[key][j] == '-':
                            test = 1
                            break
                    if test == 1:
                        break
                if test == 0:
                    align_start = i
                    break

    for i in range(len(align_dict[seq_id])):
        i = len(align_dict[seq_id]) - i - 1
        if align_dict[seq_id][i] != '-':
            test = 0
            for key in align_dict.keys():
                if align_dict[key][i] == '-':
                    test = 1
                    break
            if test == 0:
                for j in range(i-40,i):
                    for key in align_dict.keys():
                        if align_dict[key][j] == '-':
                            test = 1
                            break
                    if test == 1:
                        break
                if test == 0:
                    align_end = i
                    break
    align_dict['region'] = [align_start, align_end]                
    return align_dict

def read_ab1(path):
    ab1_list = []
    record_list = {}
    for record in os.listdir(path):
        if record.endswith('.ab1'):
            ab1_list.append(record)
    for file in ab1_list:
        record = SeqIO.read(path + '/' + file, "abi")
        channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
        trace = defaultdict(list)
        fasta = record.annotations["abif_raw"]["PBAS2"]
        ploc = record.annotations["abif_raw"]["PLOC2"]
        trace['fasta'] = fasta
        trace['ploc'] = ploc
        for channel in channels:
            trace[channel] = record.annotations["abif_raw"][channel]
        record_list[file] = trace
    return record_list
def get_mutation_or_snp(align_dict, record_list):
    record_list_new = {}
    ref_key = list(align_dict.keys())[0]
    region_key = 'region'
    files_key = list(record_list.keys())
    remove_before_start_list = {}
    for file in files_key:
        file = str(file.split('.')[0])
        remove_before_start_list[file] = 0
        for i in range(align_dict[region_key][0]):
            if align_dict[file][i] != '-':
                remove_before_start_list[file] += 1
    for i in range(align_dict[region_key][0]+10, align_dict[region_key][1]-10):
        ref_base = align_dict[ref_key][i]
        if ref_base != '-':
            test = 0
            # here we check if the base is the same as the referene base
            # here we based the fasta result.
            for file in files_key:
                file = str(file.split('.')[0])
                if align_dict[file][i] != ref_base:
                    record_list_new[str(i) + '.'+ str(ref_base)] = {}
                    for file in files_key:  
                        file = str(file.split('.')[0])
                        pos_ploc = remove_before_start_list[file] + i - align_dict[region_key][0]
                        pos_base = record_list[file+str('.ab1')]['ploc'][pos_ploc]
                        record_list_new[str(i) + '.'+ str(ref_base)][file] = pos_base
                    test = 1
                    break
            # continue if the base is the same as the reference base
            # & check if the base is snp site
            if test == 0:
                for file in files_key:
                    file = str(file.split('.')[0])
                    # here we need compare the base data to check if it is snp site
                    pos_ploc = remove_before_start_list[file] + i - align_dict[region_key][0]
                    pos_base = record_list[file+str('.ab1')]['ploc'][pos_ploc]
                    T_base_max = max(record_list[file+str('.ab1')]['DATA9'][pos_base-3:pos_base+4])
                    C_base_max = max(record_list[file+str('.ab1')]['DATA10'][pos_base-3:pos_base+4])
                    G_base_max = max(record_list[file+str('.ab1')]['DATA11'][pos_base-3:pos_base+4])
                    A_base_max = max(record_list[file+str('.ab1')]['DATA12'][pos_base-3:pos_base+4])
                    mx = max(A_base_max, T_base_max, C_base_max, G_base_max)
                    top2 = sorted([A_base_max, T_base_max, C_base_max, G_base_max])[-2]
                    peaks = int(A_base_max/mx +0.3) + int(T_base_max/mx +0.3) + int(C_base_max/mx +0.3) + int(G_base_max/mx +0.3)
                    if peaks > 1:
                        record_list_new[str(i) + '.'+ str(ref_base)] = {}
                        for file in files_key:  
                            file = str(file.split('.')[0])
                            pos_ploc = remove_before_start_list[file] + i - align_dict[region_key][0]
                            pos_base = record_list[file+str('.ab1')]['ploc'][pos_ploc]
                            record_list_new[str(i) + '.'+ str(ref_base)][file] = pos_base
                        break
    return record_list_new

def plot_snp_mut(record_list_new, record_list, output):
    # here we base on the record_list_new, which record the snp data9-12 ploc and fasta
    with open(output+'.txt', 'w') as f:
        f.write(output)
        color_group = ['#FF4136','#0074D9','#2ECC40','#FFDC00']
        fig, ax = plt.subplots()
        fig.set_size_inches(18.5, 10.5)
        ax.set_aspect('equal')
        ax.set_xlim(-1, 77)
        ax.set_ylim(-3, 25)
        snp_mut_list = list(record_list_new.keys())
        file_list = list(record_list_new[snp_mut_list[0]].keys())
        for i in range(len(snp_mut_list)):
            f.write('\n'+str(snp_mut_list[i])+'.'+str(6329308 - int(snp_mut_list[i].split('.')[0]) +1))
            for j in range(len(file_list)):
                ploc_base = record_list_new[snp_mut_list[i]][file_list[j]]
                ref_base = snp_mut_list[i].split('.')[1]
                G_line = list(record_list[file_list[j]+str('.ab1')]['DATA9'][ploc_base-5:ploc_base+6])
                A_line = list(record_list[file_list[j]+str('.ab1')]['DATA10'][ploc_base-5:ploc_base+6])
                T_line = list(record_list[file_list[j]+str('.ab1')]['DATA11'][ploc_base-5:ploc_base+6])
                C_line = list(record_list[file_list[j]+str('.ab1')]['DATA12'][ploc_base-5:ploc_base+6])
                f.write('\t')
                f.write(str([G_line,A_line,T_line,C_line]))
                for k in range(len(G_line)):
                    G_line[k] = G_line[k]/1500 + i 
                    A_line[k] = A_line[k]/1500 + i
                    T_line[k] = T_line[k]/1500 + i
                    C_line[k] = C_line[k]/1500 + i
                x = []
                for k in range(len(G_line)):
                    x.append(j+k/12+0.1)
                
                g_line = ax.plot(x, G_line, color=color_group[0], linewidth=0.5, alpha=0.8)   
                a_line = ax.plot(x, A_line, color=color_group[1], linewidth=0.5, alpha=0.8)
                t_line = ax.plot(x, T_line, color=color_group[2], linewidth=0.5, alpha=0.8)
                c_line = ax.plot(x, C_line, color=color_group[3], linewidth=0.5, alpha=0.8)
                box = plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=0.2)
                ax.add_artist(box)
        # add position label
        for i in range(len(snp_mut_list)):
            if snp_mut_list[i].split('.')[1] == 'a':
                box = plt.Rectangle((72, i), 5, 1, fill=False, edgecolor='grey', lw=0.1, alpha=0.8)
                ax.add_artist(box)
                ax.text(72.05, i+0.31, 'A', fontsize=8, color=color_group[1])
                ax.text(72.31, i+0.01, 'T', fontsize=4, color=color_group[2])
            elif snp_mut_list[i].split('.')[1] == 'c':
                box = plt.Rectangle((72, i), 5, 1, fill=False, edgecolor='grey', lw=0.1, alpha=0.8)
                ax.add_artist(box)
                ax.text(72.05, i+0.31, 'C', fontsize=8, color=color_group[3])
                ax.text(72.31, i+0.01, 'G', fontsize=4, color=color_group[0])
            elif snp_mut_list[i].split('.')[1] == 'g':
                box = plt.Rectangle((72, i), 5, 1, fill=False, edgecolor='grey', lw=0.1, alpha=0.8)
                ax.add_artist(box)
                ax.text(72.05, i+0.31, 'G', fontsize=8, color=color_group[0])
                ax.text(72.31, i+0.01, 'C', fontsize=4, color=color_group[3])
            elif snp_mut_list[i].split('.')[1] == 't':
                box = plt.Rectangle((72, i), 5, 1, fill=False, edgecolor='grey', lw=0.1, alpha=0.8)
                ax.add_artist(box)
                ax.text(72.05, i+0.31, 'T', fontsize=8, color=color_group[2])
                ax.text(72.31, i+0.01, 'A', fontsize=4, color=color_group[1])
            pos = 6329308 - int(snp_mut_list[i].split('.')[0]) +1
            ax.text(72.5, i+0.1, ':'+str(pos), fontsize=8)
        box = plt.Rectangle((0, 0), 77, len(snp_mut_list), fill=False, edgecolor='black', lw=0.4)
        ax.add_artist(box)
        box = plt.Rectangle((76.33, 0), 0.65,len(snp_mut_list), facecolor='white', edgecolor='none', alpha=0.8)
        ax.add_artist(box)
        ax.text(76.4, 0.3, output, fontsize=8, rotation=90)
        # clock wise 90 degree
        x_label = [str(i) for i in range(1, 73)]
        y_label = [str(i) for i in range(1, 11)]
        ax.set_xticks(np.arange(0.5, 72.5, 1))
        ax.set_yticks(np.arange(0.5, 10.5, 1))
        ax.tick_params(axis='both', which='both', length=0)
        ax.set_xticklabels(x_label, fontsize=8)
        ax.set_yticklabels(y_label, fontsize=8)
        ax.axis('off')
        plt.savefig(output+'.png', dpi=600)
        plt.close()
def main():
    args = get_args()
    align_dict = read_align(args.align)
    record_list = read_ab1(args.sanger)
    snp_mut_dict = get_mutation_or_snp(align_dict, record_list)
    plot_snp_mut(snp_mut_dict, record_list, args.output)

if __name__ == '__main__':
    main()