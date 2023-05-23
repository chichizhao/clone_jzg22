#!py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# estaminat the time scale of the tree, with the mutation rate of 5.5e-9
# we assume the tree is rooted, and different branches have different mutation rates
# aslo, we assume the leaf has different mutation rates
# usage: python3 evolution_tree_time.py -tree tree.nwk -mut 5.52e-9 -ml mutational_length.txt -ref typha.fa -out out.txt
# gemome_length is the length of the reference genome: 216955151
import sys
import argparse
import scipy.stats as stats
import numpy as np
from scipy import linalg
from scipy.linalg import solve

def argparse_init():
    parser = argparse.ArgumentParser(description='estaminat the time scale of the tree, with the mutation rate of 5.5e-9')
    parser.add_argument('-tree', help='the tree file in newick format')
    parser.add_argument('-mut', help='the mutation rate')
    parser.add_argument('-ml', help='the mutational length of each branch')
    parser.add_argument('-out', help='the output file')
    args = parser.parse_args()
    return args

def read_mutational_length(mutational_length_file):
    mutational_length = {}
    with open(mutational_length_file) as f:
        for line in f:
            line = line.strip()
            if line:
                line = line.split('\t')
                mutational_length[line[0]] = float(line[1])
    return mutational_length

def read_tree(tree_file):
    tree = ""
    with open(tree_file) as f:
        for line in open(tree_file):
            if line[0] == '(':
                tree = line.strip()
                break
    return tree

def get_divison_time(mutation_rate, mutational_length):
    # as the average mutation rate is 5.5e-9, we can estimate the time scale of the tree
    # we assume all the mutations are neutral, and obey the mutation rate as average mutation rate normal distribution
    # so we first fit mutational length to normal distribution, and then estimate average mutation rate
    genome_length = 216955151
    mutational_length = list(mutational_length.values())
    u = sum(mutational_length)/len(mutational_length)
    std = (sum([(x-u)**2 for x in mutational_length])/len(mutational_length))**0.5
    D,p = stats.kstest(mutational_length, 'norm', args=(u, std))
    if p > 0.05:
        print("the mutational length is normal distribution")
    else:
        print("the mutational length is not normal distribution")
    # D =  0.09882778339598808 p =  0.4539586354351576
    # here we got the mutational length is  normal distribution
    # so here the average mutation rate is 5.5e-9,and average mutation length is u 
    print(u)
    average_mutation_time = int(u/(float(mutation_rate)*genome_length))
    print("the average mutation time is %s year" % average_mutation_time)
    # the average mutation time is 1194 year
    # so we think the begining of the tree is 1194 year ago and all the samples are collected in 2022 august 
    return average_mutation_time

def calcuate_time(tree, average_mutation_time):
    gemome_length = 216955151
    mut = 5.5e-9
    
def get_node_and_branch(tree):
    node_length = {}
    leaf_length = {}
    node_name = []
    leaf_name = []
    samples = []
    samples_retriv_node = {}
    node = 0
    for i in range(len(tree)):
        if tree[i] == '(':
            node += 1
        elif tree[i] == ')':
            node -= 1
        elif tree[i] == ':':
            node_group = ""
            leaf_group = ""
            if tree[i-1] == ')':
                r_node = 1
                for j in range(i-2):
                    j = i-2-j
                    if tree[j] == ')':
                        r_node += 1
                    elif tree[j] == '(':
                        r_node -= 1
                    if r_node == 0:
                        node_group = tree[j:i]
                        break
            else:
                for j in range(i-1):
                    j = i-1-j
                    if tree[j] == ',' or tree[j] == '(':
                        leaf_group = tree[j+1:i]
                        break
            if node_group != "":
                node_name.append(node_group)
            elif leaf_group != "":
                leaf_name.append(leaf_group)
            
            for j in range(i+1, len(tree)):
                if tree[j] == ',' or tree[j] == ')':
                    if node_group != "":
                        node_length[node_group] = float(tree[i+1:j])
                    if leaf_group != "":
                        leaf_length[leaf_group] = float(tree[i+1:j])
                    break
    total_length = 0
    if len(node_name) != 0:
        for i in leaf_name:
            total_length += leaf_length[i]
            #print(total_length)
            #print(i)
    if len(leaf_name) != 0:
        for i in node_name:
            samples = get_samples(i)
            total_length += node_length[i] * len(samples)
    average_branch_length = total_length / len(leaf_name)
    return node_name, leaf_name, node_length, leaf_length, average_branch_length
def get_root_length(tree,node_tree):
    node_length = {}
    leaf_length = {}
    node_name = []
    leaf_name = []
    samples = []
    root_length = []
    root_nodes = []
    node = 0
    for i in range(len(tree)):
        if tree[i] == '(':
            node += 1
        elif tree[i] == ')':
            node -= 1
        elif tree[i] == ':':
            node_group = ""
            leaf_group = ""
            if tree[i-1] == ')':
                r_node = 1
                for j in range(i-2):
                    j = i-2-j
                    if tree[j] == ')':
                        r_node += 1
                    elif tree[j] == '(':
                        r_node -= 1
                    if r_node == 0:
                        node_group = tree[j:i]
                        break
            else:
                for j in range(i-1):
                    j = i-1-j
                    if tree[j] == ',' or tree[j] == '(':
                        leaf_group = tree[j+1:i]
                        break
            if node_group != "":
                node_name.append(node_group)
            elif leaf_group != "":
                leaf_name.append(leaf_group)
            
            for j in range(i+1, len(tree)):
                if tree[j] == ',' or tree[j] == ')':
                    if node_group != "":
                        node_length[node_group] = float(tree[i+1:j])
                    if leaf_group != "":
                        leaf_length[leaf_group] = float(tree[i+1:j])
                    break
    for i in node_name:
        node_tree_samples = get_samples(node_tree)
        samples = get_samples(i)
        for j in samples:
            if j in node_tree_samples:
                node_tree_samples.remove(j)
        if len(node_tree_samples) == 0:
            root_length.append(node_length[i])
            root_nodes.append(i)
    return root_length, root_nodes
    
def get_time( tree, mut, average_mutation_time):
    new_tree = ""
    if tree[-1] == ';':
        tree = tree[:-1]
    root_length = 0
    average_branch_length = 0
    node = 0
    estaminat_time_order = []
    get_the_time_order = get_estimate_time_order(tree,estaminat_time_order,1)
    mutation_time_each_branch_n_leaf = {}
    mutation_rate_each_branch_n_leaf = {}
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
    mutation_time_each_branch_n_leaf[root] = float(490/(5.5e-9 * 216955151))
    mutation_rate_each_branch_n_leaf[root] = float(5.5e-9)
    mutation_time_nodes_name = []
    mutation_time_nodes_name.append(root)
    gemome_length = 216955151
    u0 = 5.5e-9
    L = 1425.5 + 490
    T =  (1425.5+490)/(5.5e-9 * 216955151)
    for i in get_the_time_order:
        node_length = 0
        for j in range(len(i)):
            j = len(i) - 1 - j
            if i[j] == ':':
                new_tree = i[:j]
                node_length = float(i[j+1:])
                new_tree,avrage_branch_length = get_average_branch_length(new_tree)
                break
        root_length, root_nodes = get_root_length(tree,new_tree)
        T1 = 0
        for j in range(len(root_nodes)):
            for k in range(len(mutation_time_nodes_name)):
                if root_nodes[j] == mutation_time_nodes_name[k]:
                    T1 += mutation_time_each_branch_n_leaf[mutation_time_nodes_name[k]]
        z= (avrage_branch_length+node_length)/(T -T1)
        z = z/216955151
        t = node_length/z
        t = t/216955151
        mutation_time_nodes_name.append(new_tree)
        mutation_time_each_branch_n_leaf[new_tree] = t
        mutation_rate_each_branch_n_leaf[new_tree] = z
    estimate_time, mutation_rate = retrive_tree(tree,mutation_time_each_branch_n_leaf,mutation_rate_each_branch_n_leaf,mutation_time_nodes_name)
    return estimate_time, mutation_rate
def retrive_tree(tree,mutation_time_each_branch_n_leaf,mutation_rate_each_branch_n_leaf,mutation_time_nodes_name):
    k = 0
    new_tree = ""
    new_tree2 = ""
    for i in range(len(tree)):
        if tree[i] == ':':
            if tree[i-1] != ')':
                for j in range(i-1):
                    j = i-1-j
                    if tree[j] == ',' :
                        sample1 = tree[j+1:i]
                        for t in range(len(mutation_time_nodes_name)):
                            samples = get_samples(mutation_time_nodes_name[t])
                            if len(samples) == 1 and sample1 == samples[0]:
                                # set the time with 2 number after the point
                                if new_tree == "":
                                    new_tree = sample1 + ':' + str(round(mutation_time_each_branch_n_leaf[mutation_time_nodes_name[t]],2))
                                    new_tree2 = sample1 + ':' + str(round(mutation_rate_each_branch_n_leaf[mutation_time_nodes_name[t]]*1000000000,2))
                                else:
                                    new_tree = new_tree + ',' + sample1 + ':' + str(round(mutation_time_each_branch_n_leaf[mutation_time_nodes_name[t]],2)) + ')'
                                    new_tree2 = new_tree2 + ',' + sample1 + ':' + str(round(mutation_rate_each_branch_n_leaf[mutation_time_nodes_name[t]]*1000000000,2)) + ')'
                    if tree[j] == '(':
                        p1 = 0
                        for t in range(i):
                            t = i-1-t
                            if tree[t] == ',':
                                p1 = t 
                                break
                        sample1x = tree[p1+1:i]
                        sample1 = tree[j+1:i]
                        for t in range(len(mutation_time_nodes_name)):
                            samples = get_samples(mutation_time_nodes_name[t])
                            if len(samples) == 1 and sample1 == samples[0]:
                                # set the time with 2 number after the point
                                new_tree = new_tree+','+sample1x + ':' + str(round(mutation_time_each_branch_n_leaf[mutation_time_nodes_name[t]],2))
                                new_tree2 = new_tree2+','+sample1x + ':' + str(round(mutation_rate_each_branch_n_leaf[mutation_time_nodes_name[t]]*1000000000,2))
            else:
                node = ""
                k = 0
                for j in range(i-1):
                    j = i-1-j
                    if tree[j] == ')':
                        k += 1
                    elif tree[j] == '(':
                        k -= 1
                    if k == 0:
                        node = tree[j:i]
                        break
                sample1 = ""
                sample2 = ""
                sample1 = get_samples(node)
                for j in range(len(mutation_time_nodes_name)):
                    sample2 = get_samples(mutation_time_nodes_name[j])
                    s1 = 0
                    s2 = 0
                    for t in range(len(sample1)):
                        for k in range(len(sample2)):
                            if sample1[t] == sample2[k]:
                                s1 += 1
                                s2 += 1
                    if s1 == len(sample1) and s2 == len(sample2):
                        # set the time with 2 number after the point
                        new_tree = new_tree + ':' + str(round(mutation_time_each_branch_n_leaf[mutation_time_nodes_name[j]],2) )
                        new_tree2 = new_tree2 + ':' + str(round(mutation_rate_each_branch_n_leaf[mutation_time_nodes_name[j]]*1000000000,2))
                        if i != len(tree)-1 :
                            for t in range(i+1,len(tree)):
                                if tree[t] == ',':
                                    break
                                if tree[t] == ')':
                                    new_tree = new_tree + ')'
                                    new_tree2 = new_tree2 + ')'
                                    break
    return new_tree , new_tree2
            
def get_estimate_time_order(tree,estaminat_time_order,k):
    if k == 1 :
        estaminat_time_order = []
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
                estaminat_time_order.append(tree_right)
                estaminat_time_order.append(tree_left)
            elif len(samples_left) > 1 and len(samples_right) == 1:
                estaminat_time_order.append(tree_right)
                estaminat_time_order.append(tree_left)
                estaminat_time_order_next = get_estimate_time_order(tree_left,estaminat_time_order,k+1)
                estaminat_time_order = estaminat_time_order_next
            elif len(samples_left) == 1 and len(samples_right) > 1:
                estaminat_time_order.append(tree_left)
                estaminat_time_order.append(tree_right)
                estaminat_time_order_next = get_estimate_time_order(tree_right,estaminat_time_order,k+1)
                estaminat_time_order = estaminat_time_order_next
            elif len(samples_left) > 1 and len(samples_right) > 1:
                estaminat_time_order.append(tree_left)
                estaminat_time_order.append(tree_right)
                estaminat_time_order_next = get_estimate_time_order(tree_left,estaminat_time_order,k+1)
                estaminat_time_order = estaminat_time_order_next
                estaminat_time_order_next = get_estimate_time_order(tree_right,estaminat_time_order,k+1)
                estaminat_time_order = estaminat_time_order_next
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
            #print(tree_left)
            #print(tree_right)
            samples_left = []
            samples_right = []
            samples_left = get_samples(tree_left)
            samples_right = get_samples(tree_right)
            if len(samples_left) == 1 and len(samples_right) == 1:
                estaminat_time_order.append(tree_left)
                estaminat_time_order.append(tree_right)
            elif len(samples_left) > 1 and len(samples_right) == 1:
                estaminat_time_order.append(tree_right)
                estaminat_time_order.append(tree_left)
                estaminat_time_order_next = get_estimate_time_order(tree_left,estaminat_time_order,k+1)
                estaminat_time_order = estaminat_time_order_next
            else:
                print("error")
            break
        
    return estaminat_time_order
                      
def get_average_branch_length(tree):
    branch_length = []
    new_tree = ""
    node = 0
    average_branch_length = 0
    comma_or_left_cal = 0
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
        node_name, leaf_name, node_length, leaf_length, average_branch_length = get_node_and_branch(new_tree)
        return new_tree, average_branch_length
    elif num == 2:
        l1 = 0
        l2 = 0
        l1 = new_tree.split(',')[0].split(':')[1]
        l2 = new_tree.split(',')[1].split(':')[1].split(')')[0]
        average_branch_length = (float(l1) + float(l2)) / 2
        return new_tree, average_branch_length
    else:
        return new_tree, average_branch_length

    
def get_samples(node):
    samples = []
    node = node.strip().split(',')
    for i in node:
        if i[0] == '(':
            samples.append(i.strip().split(':')[0].split('(')[-1])
        else:
            samples.append(i.strip().split(':')[0])
    return samples 
def write_result(out, estimate_time, mutation_rate):
    with open(out+"mutation_rate.txt", "w") as f:
        f.write(str(mutation_rate))
    with open(out+"estimate_time.txt", "w") as f:
        f.write(str(estimate_time))     
    
    
def main():
    args = argparse_init()
    mutational_length = read_mutational_length(args.ml)
    tree = read_tree(args.tree)
    average_mutation_time = get_divison_time(args.mut, mutational_length)
    estimate_time,mutation_rate = get_time(tree, args.mut, average_mutation_time)
    write_result(args.out, estimate_time, mutation_rate)
    
if __name__ == "__main__":
    main()

    
    