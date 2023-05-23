#!/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# the mutation have its specific character,it mutates in on specific individual, not in the whole population.
# And the mutation has the mutation order, which the mutation happen in the early time in the time scale should be checked first.
# here we based on the mutation site time and the population vcf file to check the mutation site whether it is a non-synonymous mutation
# FUNCTION: check the mutation site individual by individual, and check whether it is a non-synonymous mutation
# INPUT: the population vcf file, the mutation time tree file, the gene annotation file, reference genome
# OUTPUT: the non-synonymous mutation site file individual by individual in column
# Usage: python3 non_synonymous_check.py -vcf solid_somatic_mutation.vcf -tree typha_somatic_mutation_tree_time.txt -gff annotation.gff3 -mut_type mutation_type.txt -ref typha.fa -out non_synonymous_mutation

import argparse
import re
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

codon_table_dna = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
                'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

def get_args():
    parser = argparse.ArgumentParser(description="get the non-synonymous mutation site")
    parser.add_argument('-vcf', '--vcf_file', type=str, help='the population vcf file')
    parser.add_argument('-gff', '--gff_file', type=str, help='the gene annotation file')
    parser.add_argument('-ref', '--ref_seq', type=str, help='the reference genome')
    parser.add_argument('-tree', '--mutation_time_tree', type=str, help='the mutation time tree file')
    parser.add_argument('-mut_type', '--mutation_type', type=str, help='the mutation type file')
    parser.add_argument('-out', '--out_file', type=str, help='the output file')
    args = parser.parse_args()
    return args
def get_tree(tree_file):
    tree = ""
    with open(tree_file) as f:
        for line in open(tree_file):
            if line[0] == '(':
                tree = line.strip()
                break
    return tree

def read_vcf(vcf_file):
    vcf_dict = {}
    sample_list = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                line = line.strip().split('\t')
                sample_list = line[9:]
            else:
                line = line.strip().split('\t')
                chr = line[0]
                start = line[1]
                ref = line[3]
                alt = line[4]
                vcf_dict[chr + '_' + start] = [ref, alt]
    return vcf_dict
def get_vcf_record(vcf_file):
    vcf = open(vcf_file, 'r')
    vcf_record = {}
    sample_name = []
    for line in vcf:
        if line[0] == '#':
            if line[1] == '#':
                continue
            else:
                sample_name = line.strip().split()
                sample_name = sample_name[9:]
                continue
        else:
            line = line.strip().split()
            chr = line[0]
            pos = line[1]
            gene_type = []
            for i in range(9,len(line)):
                if line[i].split(':')[0] == '0/0' or line[i].split(':')[0] == './.':
                    gene_type.append('0')
                else:
                    gene_type.append('1')
            vcf_record[chr + '_' + pos] = gene_type
            
    return vcf_record, sample_name

def get_the_gene_time(tree):
    gene_time = {}
    new_tree = ""
    if tree[-1] == ';':
        tree = tree[:-1]
    before_node_time = 0
    after_node_time = 0
    node_time_length = 0
    time_order = []
    time_order = get_estimate_time_order(tree,time_order,1)
    k = 2 
    for i in range(len(tree)):
        if tree[i] == '(':
            k = k - 1
            if k == 0:
                lft = i
                break
    k = 2
    for i in range(len(tree)):
        i = len(tree) - 1 - i
        if tree[i] == ')':
            k = k - 1
            if k == 0:
                rgt = i
                break
    root = tree[lft:rgt+1]
    if root not in gene_time:
        gene_time[root] = [1194.63,1194.63+410.64]
    for i in range(len(time_order)):
        after_node_time = 0
        node_time_length = 0
        for j in range(len(time_order[i])):
            j = len(time_order[i]) - 1 - j
            if time_order[i][j] == ':':
                node_time_length = float(time_order[i][j+1:])
                newtree = time_order[i][:j]
                break
        after_node_time = get_after_node_time(newtree)
                
        if time_order[i] not in gene_time:
            gene_time[time_order[i]] = [after_node_time,after_node_time+node_time_length]
    return gene_time

def get_after_node_time(tree):
    after_node_time = 0
    new_tree = ""
    comma_or_left_cal = 0
    node = 0
    for i in range(len(tree)):
        i = len(tree) - 1 - i
        if tree[i] == ')':
            node += 1
        elif node == 0 and (tree[i] == ',' or tree[i] == '('):
            new_tree = tree[i+1:]
            break
        elif tree[i] == ',':
            comma_or_left_cal += 1
        elif tree[i] == '(':
            node -= 1
            comma_or_left_cal += 1
        elif node == 0 and comma_or_left_cal != 0:
            new_tree = tree[i+1:]
            break
    new_tree = tree
    num = len(get_samples(new_tree))
    if num > 2:
        fake_left = 0
        node = 0
        k = 0
        for i in range(len(new_tree)):
            if new_tree[i] == '(':
                if k == 1:
                    fake_left += 1
            elif new_tree[i] == ')':
                node -= 1
                if k == 1 and fake_left > 0:
                    fake_left -= 1
            elif new_tree[i] == ',' :
                node += 1
            if new_tree[i] == ':' and k == 0:
                k = 1
                for j in range(i,len(new_tree)):
                    if new_tree[j] == ',' or new_tree[j] == ')':
                        after_node_time += float(new_tree[i+1:j])
                        break
            elif new_tree[i] == ':' and k == 1 and fake_left == 0 and node == 0:
                for j in range(i,len(new_tree)):
                    if new_tree[j] == ',' or new_tree[j] == ')':
                        after_node_time += float(new_tree[i+1:j])
                        break
    elif num == 2:
        for i in range(len(new_tree)):
            if new_tree[i] == ':':
                for j in range(i,len(new_tree)):
                    if new_tree[j] == ',' or new_tree[j] == ')':
                        after_node_time += float(new_tree[i+1:j])
                        break
                break
    
    return after_node_time
                   
def get_samples(node):
    samples = []
    node = node.strip().split(',')
    for i in node:
        if i[0] == '(':
            samples.append(i.strip().split(':')[0].split('(')[-1])
        else:
            samples.append(i.strip().split(':')[0])
    return samples 

def get_mutation_type(mutation_type_file):
    mutation_type = {}
    mutation_type_file = open(mutation_type_file, 'r')
    for line in mutation_type_file:
        line = line.strip().split(',')
        chr = line[0]
        pos = line[1]
        type = line[2]
        mutation_type[chr + '_' + pos] = type
    return mutation_type
def get_estimate_time_order(tree,time_order,k):
    if k == 1 :
        time_order = []
        tree = tree[1:-1]
    node = 0
    l = 0
    for i in range(len(tree)):
        i = len(tree) - 1 - i
        if tree[i] == ')':
            node += 1
        elif tree[i] == '(':
            node -= 1
            l += 1
        if node == 1 and l != 0 and tree[i] == ',':
            for j in range(i):
                if tree[j] == '(':
                    tree_left = tree[j+1:i]
                    break
            for j in range(len(tree)):
                j = len(tree) - 1 - j
                if tree[j] == ')':
                    tree_right = tree[i+1:j]
                    break
            samples_left = []
            samples_right = []
            samples_left = get_samples(tree_left)
            samples_right = get_samples(tree_right)
            if len(samples_left) == 1 and len(samples_right) == 1:
                time_order.append(tree_right)
                time_order.append(tree_left)
            elif len(samples_left) > 1 and len(samples_right) == 1:
                time_order.append(tree_right)
                time_order.append(tree_left)
                time_order_next = get_estimate_time_order(tree_left,time_order,k+1)
                time_order = time_order_next
            elif len(samples_left) == 1 and len(samples_right) > 1:
                time_order.append(tree_left)
                time_order.append(tree_right)
                time_order_next = get_estimate_time_order(tree_right,time_order,k+1)
                time_order = time_order_next
            elif len(samples_left) > 1 and len(samples_right) > 1:
                time_order.append(tree_left)
                time_order.append(tree_right)
                time_order_next = get_estimate_time_order(tree_left,time_order,k+1)
                time_order = time_order_next
                time_order_next = get_estimate_time_order(tree_right,time_order,k+1)
                time_order = time_order_next
            break
        elif node == 1 and l == 0 and tree[i] == ',':
            tree_left = ""
            tree_right = ""
            for j in range(i):
                if tree[j] == '(':
                    tree_left = tree[j+1:i]
                    break
            for j in range(len(tree)):
                j = len(tree) - 1 - j
                if tree[j] == ')':
                    tree_right = tree[i+1:j]
                    break
            samples_left = []
            samples_right = []
            samples_left = get_samples(tree_left)
            samples_right = get_samples(tree_right)
            if len(samples_left) == 1 and len(samples_right) == 1:
                time_order.append(tree_left)
                time_order.append(tree_right)
            elif len(samples_left) > 1 and len(samples_right) == 1:
                time_order.append(tree_right)
                time_order.append(tree_left)
                time_order_next = get_estimate_time_order(tree_left,time_order,k+1)
                time_order = time_order_next
            else:
                print("error")
            break
        
    return time_order


def get_after_node_time(tree):
    after_node_time = 0
    new_tree = ""
    comma_or_left_cal = 0
    node = 0
    for i in range(len(tree)):
        i = len(tree) - 1 - i
        if tree[i] == ')':
            node += 1
        elif node == 0 and (tree[i] == ',' or tree[i] == '('):
            new_tree = tree[i+1:]
            break
        elif tree[i] == ',':
            comma_or_left_cal += 1
        elif tree[i] == '(':
            node -= 1
            comma_or_left_cal += 1
        elif node == 0 and comma_or_left_cal != 0:
            new_tree = tree[i+1:]
            break
    new_tree = tree
    num = len(get_samples(new_tree))
    if num > 2:
        fake_left = 0
        node = 0
        k = 0
        for i in range(len(new_tree)):
            if new_tree[i] == '(':
                if k == 1:
                    fake_left += 1
            elif new_tree[i] == ')':
                node -= 1
                if k == 1 and fake_left > 0:
                    fake_left -= 1
            elif new_tree[i] == ',' :
                node += 1
            if new_tree[i] == ':' and k == 0:
                k = 1
                for j in range(i,len(new_tree)):
                    if new_tree[j] == ',' or new_tree[j] == ')':
                        after_node_time += float(new_tree[i+1:j])
                        break
            elif new_tree[i] == ':' and k == 1 and fake_left == 0 and node == 0:
                for j in range(i,len(new_tree)):
                    if new_tree[j] == ',' or new_tree[j] == ')':
                        after_node_time += float(new_tree[i+1:j])
                        break
    elif num == 2:
        for i in range(len(new_tree)):
            if new_tree[i] == ':':
                for j in range(i,len(new_tree)):
                    if new_tree[j] == ',' or new_tree[j] == ')':
                        after_node_time += float(new_tree[i+1:j])
                        break
                break
    
    return after_node_time
                   
def get_samples(node):
    samples = []
    node = node.strip().split(',')
    for i in node:
        if i[0] == '(':
            samples.append(i.strip().split(':')[0].split('(')[-1])
        else:
            samples.append(i.strip().split(':')[0])
    return samples 

def annotate_mutationsite_time(vcf_dict, vcf_record, gene_time, mutation_type, sample_name):
    mutation_time = {}
    for i in sample_name:
        mutation_time[i] = []
    for i in vcf_record:
        if i in mutation_type:
            for j in gene_time:
                samples = get_samples(j)
                test = 0
                pos_s = []
                sample_name_list = []
                for k in range(len(samples)):
                    pos = sample_name.index(samples[k])
                    if vcf_record[i][pos] == '1':
                        test = 2
                        sample_name_list.append(samples[k])
                        pos_s.append(pos+1)
                        continue
                    else:
                        test = 1
                        break
                if test == 2: 
                    chr_pos = i
                    ref = vcf_dict[i][0]
                    alt = vcf_dict[i][1]
                    start = int(gene_time[j][0])
                    end = int(gene_time[j][1])
                    type = mutation_type[i]
                    for k in pos_s:
                        mutation_time[sample_name[k-1]].append([chr_pos,ref,alt,start,end,type])
                    for k in range(len(pos_s)):
                        vcf_record[i][pos_s[k]-1] = '0'
    return mutation_time

def get_ref_seq(genome_fasta_file):
    ref_seq = {}
    with open(genome_fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chromosome = line.strip().split('>')[1]
                ref_seq[chromosome] = []
            else:
                ref_seq[chromosome].append(line.strip())
    for chromosome in ref_seq:
        ref_seq[chromosome] = ''.join(ref_seq[chromosome])
    return ref_seq

def read_gff(gff_file):
    gff_dict = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line = line.strip().split('\t')
                if line[2] == 'gene':
                    ID = line[8].split(';')[0].split('=')[1]
                    start = line[3]
                    end = line[4]
                    gff_dict[ID] = [start,end]   
                elif line[2] == 'CDS':
                    Parent = line[8].split(';')[1].split('=')[1]
                    start = line[3]
                    end = line[4]
                    gff_dict[ID].append([start,end])
    return gff_dict
def check_non_synonymous(ref_seq, mutation_time, gff_dict):
    non_synonymous = {}
    samples = mutation_time.keys()
    for sample in samples:
        non_synonymous_sample = {}
        non_synonymous_sample_seq = {}
        sample_mutation = {}
        for i in mutation_time[sample]:
            pos = i[0]
            sample_mutation[pos] = i[1:]
        # sort the sample_mutation by the change start time, which is the 3rd element in the list in descending order
        sample_mutation = sorted(sample_mutation.items(), key=lambda x:x[1][2], reverse=True)
        #based on the mution time, here we can orderly chech the mutation whether it occurs in non-synonymous type
        for i in sample_mutation:
            if i[1][4] != 'noncoding':
                mut_pos = int(i[0].split('_')[1])
                ref = i[1][0]
                mut = i[1][1]
                start_time = int(i[1][2])
                end_time = int(i[1][3])
                cds_pos = gff_dict[i[1][4]][2:]
                for j in range(int(len(cds_pos))):
                    if int(mut_pos) >= int(cds_pos[j][0]) and int(mut_pos) <= int(cds_pos[j][1]):
                        if len(ref) != len(mut):
                            if non_synonymous_sample == {} or str(i[1][4])not in non_synonymous_sample:
                                start_pos = int(gff_dict[i[1][4]][0])
                                end_pos = int(gff_dict[i[1][4]][1])
                                non_synonymous_sample_seq[str(i[1][4])]=ref_seq[i[0].split('_')[0]][start_pos:mut_pos]+mut+ref_seq[i[0].split('_')[0]][mut_pos+len(ref):end_pos+1]
                                non_synonymous_sample[str(i[0])] = [gff_dict[i[1][4]],mut_pos,ref,mut,start_time,end_time]
                                if str(i[1][4])not in non_synonymous_sample:
                                    non_synonymous_sample[str(i[1][4])] = [i[0]]
                                else:
                                    non_synonymous_sample[str(i[1][4])].append(i[0])
                            else:
                                move_pos = 0
                                for k in non_synonymous_sample[str(i[1][4])]:
                                    if int(mut_pos) > int(non_synonymous_sample[k][1]):
                                        move_pos += len(mut)-len(ref)
                                start_pos = int(gff_dict[i[1][4]][0])
                                end_pos = int(gff_dict[i[1][4]][1])
                                non_synonymous_sample_seq[str(i[1][4])]=non_synonymous_sample_seq[str(i[1][4])][:mut_pos-start_pos+move_pos]+mut+non_synonymous_sample_seq[str(i[1][4])][mut_pos-start_pos+len(ref)+move_pos:end_pos+move_pos+1]
                                non_synonymous_sample[str(i[0])] = [gff_dict[i[1][4]],mut_pos,ref,mut,start_time,end_time]
                                if str(i[1][4])not in non_synonymous_sample:
                                    non_synonymous_sample[str(i[1][4])] = [i[0]]
                                else:
                                    non_synonymous_sample[str(i[1][4])].append(i[0])
                            break
                        else:
                            if non_synonymous_sample == {} or str(i[1][4])not in non_synonymous_sample:
                                sample_ref_seq = ref_seq[i[0].split('_')[0]][int(cds_pos[j][0]):int(cds_pos[j][1])+1]
                                sample_mut_seq = ref_seq[i[0].split('_')[0]][int(cds_pos[j][0]):mut_pos]+mut+ref_seq[i[0].split('_')[0]][mut_pos+len(ref):int(cds_pos[j][1])+1]
                                start_pos = int(gff_dict[i[1][4]][0])
                                end_pos = int(gff_dict[i[1][4]][1])
                                for k in range(int(len(sample_ref_seq)/3)):
                                    if codon_table_dna[sample_ref_seq[3*k:3*k+3]] != codon_table_dna[sample_mut_seq[3*k:3*k+3]]:
                                        # here is non-synonymous mutation
                                        non_synonymous_sample_seq[str(i[1][4])]=ref_seq[i[0].split('_')[0]][start_pos:mut_pos]+mut+ref_seq[i[0].split('_')[0]][mut_pos+len(ref):end_pos+1]
                                        non_synonymous_sample[str(i[0])] = [gff_dict[i[1][4]],mut_pos,ref,mut,start_time,end_time]
                                        if str(i[1][4])not in non_synonymous_sample:
                                            non_synonymous_sample[str(i[1][4])] = [i[0]]
                                        else:
                                            non_synonymous_sample[str(i[1][4])].append(i[0])
                                        break
                                    else:
                                        continue
                            else:
                                move_pos = 0
                                for k in non_synonymous_sample[str(i[1][4])]:
                                    if int(mut_pos) > int(non_synonymous_sample[k][1]):
                                        move_pos += len(mut)-len(ref)
                                start_pos = int(gff_dict[i[1][4]][0])
                                end_pos = int(gff_dict[i[1][4]][1])
                                sample_ref_seq = non_synonymous_sample_seq[str(i[1][4])][int(cds_pos[j][0])-start_pos+move_pos:int(cds_pos[j][1])-start_pos+move_pos+1]
                                sample_mut_seq = sample_ref_seq[:mut_pos-int(cds_pos[j][0])+move_pos]+mut+sample_ref_seq[mut_pos-int(cds_pos[j][0])+move_pos+len(ref):]
                                for k in range(int(len(sample_ref_seq)/3)):
                                    if codon_table_dna[sample_ref_seq[3*k:3*k+3]] != codon_table_dna[sample_mut_seq[3*k:3*k+3]]:
                                        # here is non-synonymous mutation
                                        non_synonymous_sample_seq[str(i[1][4])]=non_synonymous_sample_seq[str(i[1][4])][:mut_pos-start_pos+move_pos]+mut+non_synonymous_sample_seq[str(i[1][4])][mut_pos-start_pos+len(ref)+move_pos:end_pos+move_pos+1]
                                        non_synonymous_sample[str(i[0])] = [gff_dict[i[1][4]],mut_pos,ref,mut,start_time,end_time]
                                        if str(i[1][4])not in non_synonymous_sample:
                                            non_synonymous_sample[str(i[1][4])] = [i[0]]
                                        else:
                                            non_synonymous_sample[str(i[1][4])].append(i[0])
                                        break
                                    else:
                                        continue
                            break
                    else:
                        continue
        non_synonymous[str(sample)] = []
        for key in non_synonymous_sample.keys():
            if key.startswith('e'):
                mutation_in_sample = []
                mutation_in_sample.append(key)
                for i in non_synonymous_sample[key]:
                    mutation_in_sample.append(i)
                    mutation_in_sample.append(non_synonymous_sample[i][-2])
                    mutation_in_sample.append(non_synonymous_sample[i][-1])
                non_synonymous[str(sample)].append(mutation_in_sample)
            else:
                continue
    return non_synonymous

def non_synonymous_mutation_site_matrix(non_synonymous,gene_time,outfile):
    non_synonymous_matrix = {}
    # here we want plot a matrix that with the time as the y axis samples as the x axis
    # so here we need give y value  to each sample
    # here we conut the mutation sites number which locate in the same time
    # here we set the time is 2000 years, and 1 year as the unit
    # and we have 72 samples
    # so the matrix is 2000*72 
    non_synonymous_num = {}
    non_synonymous_num2 = {}
    non_synonymous_num_sample = {}
    for sample in non_synonymous.keys():
        non_synonymous_num_sample[str(sample)] = 0
        for i in non_synonymous[sample]:
            non_synonymous_num_sample[str(sample)] += 1
            start_time = int(i[2])
            end_time = int(i[3])
            if start_time == 0:
                if str(start_time)+'_'+str(end_time) + '_'+str(sample) not in non_synonymous_num:
                    non_synonymous_num[str(start_time)+'_'+str(end_time) + '_'+str(sample)] = 1
                else:
                    non_synonymous_num[str(start_time)+'_'+str(end_time) + '_'+str(sample)] += 1
            else:
                if str(start_time)+'_'+str(end_time) not in non_synonymous_num:
                    non_synonymous_num[str(start_time)+'_'+str(end_time)] = 1
                else:
                    non_synonymous_num[str(start_time)+'_'+str(end_time)] += 1
    print(non_synonymous_num)
    print(non_synonymous_num_sample)
    for sample in non_synonymous_num.keys():
        if len(sample) >2:
                # as the calculate unit is the sample by sample 
                # here we need use the number of samples to get the real num
                for key in gene_time.keys():
                    if int (gene_time[key][0]) > 0 and int(gene_time[key][0]) == int(sample.split('_')[0]) and int(gene_time[key][1]) == int(sample.split('_')[1]):
                        non_synonymous_num2[str(sample)] = non_synonymous_num[str(sample)]/len(get_samples(key))
                        break
                    else:
                        non_synonymous_num2[str(sample)] = non_synonymous_num[str(sample)]
    print(non_synonymous_num2)
    write_non_synonymous_num_for_plot(non_synonymous_num2,outfile)   
    for sample in non_synonymous.keys():
        non_synonymous_matrix[str(sample)] = [0]*2000
        for i in non_synonymous[sample]:
            start_time = int(i[2])
            end_time = int(i[3])
            for j in range(start_time,end_time):
                non_synonymous_matrix[str(sample)][j] += 1
    # we sort the matrix by sample name
    new_matrix = {}
    for i in range(1,73):
        new_matrix[str(i)] = non_synonymous_matrix[str(i)]
    #print(new_matrix)
    return new_matrix
def write_non_synonymous_num_for_plot(non_synonymous_num,outfile):
    with open(outfile+'_non_synonymous_num_for_plot.txt','w') as f:
        for i in non_synonymous_num.keys():
            if len(i) > 2:
                f.write(i.split('_')[0]+'\t'+i.split('_')[1]+'\t'+str(non_synonymous_num[i])+'\n')
    

def plot_heatmap(non_synonymous_matrix):
    df = pd.DataFrame(non_synonymous_matrix)
    path = os.getcwd()
    data = df.values
    # Define the figure and axis objects
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 8)
    ax.set_xlim([-70, 1700])
    ax.set_ylim([-5, 75])
    max_value = 20
    for j in range(0,72):
        for i in range(0,1600):
            if data[i][j] != 0:
                color = plt.cm.Reds(data[i][j]/max_value)
                box = plt.Rectangle((i,j), 1, 1, facecolor=color, edgecolor=color, linewidth=0.1)
                ax.add_artist(box)
        box = plt.Rectangle((0, j), 1600, 1, fill=False, edgecolor='grey', lw=0.1)
        ax.add_artist(box)
    box = plt.Rectangle((-30, 0), 1630, 72, fill=False, edgecolor='black', lw=0.4)
    ax.add_artist(box)
    for i in range(72):
        ax.text(-15, i+0.4, str(i+1), ha='center',va='center', fontsize=5.5)
        box = plt.Rectangle((-30, i), 1630, 1, fill=False, edgecolor='grey', lw=0.1)
        ax.add_artist(box)
    box = plt.Rectangle((-30,0), 30, 72, fill=False, edgecolor='black', lw=0.2)
    ax.add_artist(box)
    ax.text(-65, 36, 'Samples', ha='center',va='center', fontsize=12,rotation=90)
    for i in range(0,1700,400):
        ax.text(i, -1.5, str(i), ha='center', fontsize=9)
        #box = plt.Rectangle((i, -0.5), 1, 72.5, fill='grey', edgecolor='grey', lw=0.1)
        #ax.add_artist(box)
        box = plt.Rectangle((i, -0.5), 1, 0.5, facecolor='black', edgecolor='black', linewidth=0.4)
        ax.add_artist(box)
    box = plt.Rectangle((0, -0.1), 1600, 0.1, facecolor='black', edgecolor='black', linewidth=0.4)
    ax.add_artist(box)
    ax.text(800, -3.2, 'Time', ha='center', fontsize=12)
    ax.text(1500, -3, 'years', ha='center', fontsize=12)
    for i in range(0,13):
        box = plt.Rectangle((1630, 5+i*4), 20, 4, facecolor=plt.cm.Reds(i/20), edgecolor='black', linewidth=0.1)
        ax.add_artist(box)
    box = plt.Rectangle((1630, 1), 20, 4, fill=False, edgecolor='black', lw=0.2)
    ax.add_artist(box)
    box = plt.Rectangle((1630, 1), 20, 56, fill=False, edgecolor='black', lw=0.2)
    ax.add_artist(box)
    ax.text(1655, 2.5, '0', ha='left', fontsize=9)
    ax.text(1655, 54.5, '12', ha='left', fontsize=9)
    ax.text(1685, 10, 'Non-synonymous mutation number', ha='center', fontsize=12,rotation=90)
    plt.tight_layout()  
    ax.axis('off')
    plt.savefig(path+'/mutation_num_heatmap.png',dpi=600)
    
def write_non_synonymous_mutation_site_matrix(non_synonymous_matrix):
    with open('non_synonymous_mutation_site_matrix.txt','w') as f:
        for i in non_synonymous_matrix.keys():
            f.write(i+'\t'+'\t'.join([str(j) for j in non_synonymous_matrix[i]])+'\n')
def write_non_synonymous(non_synonymous):
    with open('non_synonymous_samples.txt','w') as f:
        for i in non_synonymous.keys():
            f.write(i+'\t'+'\t'.join([str(j) for j in non_synonymous[i]])+'\n')

                                                 
def main():
    args = get_args()
    vcf_file = args.vcf_file
    gff_file = args.gff_file
    ref_seq = args.ref_seq
    tree_file = args.mutation_time_tree
    mut_type = args.mutation_type
    out_file = args.out_file
    tree = get_tree(tree_file)
    vcf_record , sample_name = get_vcf_record(vcf_file)
    mutation_type = get_mutation_type(mut_type)
    gene_time = get_the_gene_time(tree)
    vcf_dict = read_vcf(vcf_file)
    mutation_time=annotate_mutationsite_time(vcf_dict, vcf_record, gene_time, mutation_type,sample_name)
    ref_seq = get_ref_seq(ref_seq)
    gff_dict = read_gff(gff_file)
    # here we need notice that the site are two strings in the dna sequence
    # so we need check it by strings while the vcf file do not have the information
    # so here we think all the mutation happend in the same single string
    non_synonymous=check_non_synonymous(ref_seq, mutation_time, gff_dict)
    write_non_synonymous(non_synonymous)
    non_synonymous_matrix = non_synonymous_mutation_site_matrix(non_synonymous,gene_time,out_file)
    write_non_synonymous_mutation_site_matrix(non_synonymous_matrix)
    plot_heatmap(non_synonymous_matrix)
    

if __name__ == '__main__':
    main()