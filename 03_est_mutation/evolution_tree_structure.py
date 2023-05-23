#/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: built the evolution tree of the jzg_22 with snps and indels
# usage: python3 evolution_tree.py -i input.vcf -o output

import argparse
def get_args():
    parser = argparse.ArgumentParser(description="built the evolution tree of the jzg_22 with snps and indels")
    parser.add_argument("-i", "--input", help="input vcf file", required=True)
    parser.add_argument("-o", "--output", help="output file", required=True)
    args = parser.parse_args()
    return args
def count_the_same_n_diff_site_num(vcf_file):
    same_num = {}
    diff_num = {}
    sample_name = []
    sample_vcf = {}
    same_vcf = {}
    diff_vcf = {}
    same_vcf_pop = {}
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                # get the sample name
                if line.startswith("#CHROM"):
                    line = line.strip().split("\t")
                    sample_name = line[9:]
                    sample_num = len(sample_name)
                    for i in range(sample_num):
                        for j in range(i, sample_num):
                            same_num[sample_name[i] + "_" + sample_name[j]] = 0
                            diff_num[sample_name[i] + "_" + sample_name[j]] = 0
                            same_vcf[sample_name[i] + "_" + sample_name[j]] = []
                            diff_vcf[sample_name[i] + "_" + sample_name[j]] = []        
                else:
                    continue
            else:
                line = line.strip().split("\t")
                # get the number of different sites
                gts = []
                for i in range(sample_num):
                    if line[i+9].startswith("./.") or line[i+9].startswith("0/0"):
                        gts.append(0)
                    else:
                        gts.append(1)
                chr = line[0]
                pos = line[1]
                if sum(gts) == len(gts):
                    same_vcf_pop[chr + "_" + pos] = gts
                else:
                    sample_vcf[chr + "_" + pos] = gts
                for i in range(sample_num):
                    for j in range(i, sample_num):
                        if gts[i] == gts[j] and gts[i] == 1:
                            same_num[sample_name[i] + "_" + sample_name[j]] += 1
                            same_vcf[sample_name[i] + "_" + sample_name[j]].append(chr + "_" + pos)
                        elif gts[i] == gts[j] and gts[i] == 0:
                            continue
                        else:
                            diff_num[sample_name[i] + "_" + sample_name[j]] += 1
                            diff_vcf[sample_name[i] + "_" + sample_name[j]].append(chr + "_" + pos)
    return same_num, diff_num, sample_name, same_vcf, diff_vcf, same_vcf_pop,sample_vcf

def write_to_file(same_num, diff_num, output_file, sample_name):
    with open((output_file + "_same.csv"), "w") as same:
        str_sample_name = ",".join(sample_name)
        same.write('sample_name' + ',' + str_sample_name + '\n')
        for i in range(len(sample_name)):
            for j in range(len(sample_name)):
                if j == 0:
                    same.write(sample_name[i] +',')
                    if i == j:
                        same.write(str(same_num[sample_name[i] + "_" + sample_name[j]]) + ',')
                    else:
                        same.write(str(same_num[sample_name[j] + "_" + sample_name[i]]) + ',')
                elif j == len(sample_name) - 1:
                    same.write(str(same_num[sample_name[i] + "_" + sample_name[j]]) + '\n')
                elif j > i :
                    same.write(str(same_num[sample_name[i] + "_" + sample_name[j]]) + ',')
                else:
                    same.write(str(same_num[sample_name[j] + "_" + sample_name[i]]) + ',')                   
    with open((output_file + "_diff.csv"), "w") as diff:
        diff.write('sample_name' + ',' + str_sample_name + '\n')
        for i in range(len(sample_name)):
            for j in range(len(sample_name)):
                if j == 0:
                    diff.write(sample_name[i] +',')
                    if i == j:
                        diff.write(str(diff_num[sample_name[i] + "_" + sample_name[j]]) + ',')
                    else:
                        diff.write(str(diff_num[sample_name[j] + "_" + sample_name[i]]) + ',')
                elif j == len(sample_name) - 1:
                    diff.write(str(diff_num[sample_name[i] + "_" + sample_name[j]]) + '\n')
                else:
                    if int(i) < int(j):
                        diff.write(str(diff_num[sample_name[i] + "_" + sample_name[j]]) + ',')
                    else:
                        diff.write(str(diff_num[sample_name[j] + "_" + sample_name[i]]) + ',')

def build_tree(same_num, diff_num, sample_name, same_vcf, diff_vcf, same_vcf_pop, sample_vcf):
    # the 2 samples with the fewest number of same sites are the farthest
    # we need find the fewest number of same sites and the sample name
    sample_a, sample_b = get_first_and_second_sample_name(same_num)
    root_length = get_root_length(same_vcf_pop)
    #print(root_length)
    # build the tree
    tree_strcture = " "
    tree_strcture = "(" + sample_a  + "," + sample_b + ")"
    leaf_length ,pos_sample_a, pos_sample_b = get_first_and_second_sample_leaf_length(sample_a, sample_b, sample_name, sample_vcf)        
    node_list_dict, node_list_len = get_node_list(tree_strcture, root_length)
    # we generally divide the tree into 2 parts, left and right
    group_a, group_b, a_sample_vcf, b_sample_vcf = get_groupa_n_groupb(sample_a, sample_b, sample_name, sample_vcf,root_length,same_num, pos_sample_a, pos_sample_b)
    #print(b_sample_vcf)
    tree =[]
    # k equals to 2^n, n is the number of samples
    k = 1
    #print(len(group_a))
    sample_num = len(group_a)
    tree_left =get_structrue_b(a_sample_vcf, group_a, group_a,k,tree,sample_num)
    #print(tree_left)
    k = 1
    sample_num = len(group_b)
    tree_right = get_structrue_b(b_sample_vcf, group_b, group_b,k,tree,sample_num)
    #print(tree_right)
    tree= "(" + tree_left + "," + tree_right + ")"+ ":" + str(int(root_length)) 
    #print(tree)
    tree, mutational_length, leaf_length = get_leaf_length(tree, sample_name, sample_vcf,root_length)
    tree = "(" + tree + ");"
    print(tree)
    
    return tree, mutational_length, leaf_length
    
def get_leaf_length(tree, sample_name, sample_vcf,root_length_x):
    sample_root_length = {}
    nodes = 0 
    for i in range(len(tree)):
        if tree[i] == "(" :
            nodes = nodes + 1
        elif tree[i] == ")":
            real_node = nodes
            fake_node = 0
            root_length = 0
            sample_name1 = ""
            sample_name2 = ""
            t = 1
            for j in range(i-1):
                j = i - j-1
                if tree[j] == ":": 
                    t = 2
                elif tree[j] == "," and t == 1:
                    sample_name1 = tree[j+1:i]
                    p = j
                elif tree[j] == "(" and t == 1:
                    sample_name2 = tree[j+1:p]
                    t = 2
            for j in range(i):
                if tree[j] == ")":
                    real_node = real_node - 1
                    
            t =1
            for j in range(i+1, len(tree)):
                if tree[j] == "(": 
                    fake_node = fake_node + 1
                if fake_node > 0 :
                    if tree[j] == ")":
                        fake_node = fake_node - 1
                        t = 0
                if fake_node == 0 and real_node > 0:
                    if tree[j] == ":" and t == 0:
                        t = 1
                    elif tree[j] == ":" and t == 1:
                        real_node = real_node - 1
                        for k in range(j+1, len(tree)):
                            if tree[k] == "," or tree[k] == ")":
                                #print(tree[j+1:k])
                                root_length += int(tree[j+1:k])
                                break
                if real_node == 0 and fake_node == 0 and t == 1:
                    if sample_name1 != "":
                        sample_root_length[sample_name1] = root_length
                        #print("sample_name1")
                        #print(sample_name1)
                        #print(root_length)
                    if sample_name2 != "":
                        sample_root_length[sample_name2] = root_length
                    t = 2
    mutational_length = {}
    for i in range(len(sample_name)):
        mutational_length[sample_name[i]] = 0
    for vcf in sample_vcf:
        for i in range(len(sample_name)):
            if sample_vcf[vcf][i] == 1:
                mutational_length[sample_name[i]] += 1
    #print(sample_root_length)
    #print(mutational_length)
    leaf_length = {}
    #print(root_length_x)
    for i in range(len(sample_name)):
        leaf_length[sample_name[i]] = mutational_length[sample_name[i]] - sample_root_length[sample_name[i]] 
    #print(leaf_length)
    new_tree= ""
    for i in range(len(tree)):
        if tree[i] == ",":
            p = i
            break
    sample = tree[:p].replace("(","").replace(")","").replace(":","").replace(" ","")
    new_tree= tree[:p] + ":" + str(leaf_length[sample])+ ',' 
    #print(new_tree)
    K =0
    for i in range(len(tree)):
        comma1 = 0
        comma2 = 0 
        if tree[i] == ",":
            comma1 = i
            for j in range(i+1, len(tree)):
                if tree[j] == ":":
                    k = 1
                if tree[j] == ",":
                    comma2 = j
                    break
        if comma1 > 0 and comma2 == 0:
            # the last one
            if k == 1:
                sample_name1 = tree[comma1+1:].split(":")[0].replace("(","").replace(")","").replace(":","").replace(" ","")
            else:
                sample_name1 = tree[comma1+1:].replace("(","").replace(")","").replace(":","").replace(" ","")
        
            for t in range(comma1+1, len(tree)):
                if sample_name1 == tree[t:t+len(sample_name1)]:
                    new_tree2 = tree[comma1+1:t]+ sample_name1 + ":" + str(leaf_length[sample_name1]) + str(tree[t+len(sample_name1):len(tree)])
                    break
            new_tree = new_tree + new_tree2
        elif comma1 > 0 and comma2 > 0:
            if k == 1:
                sample_name1 = tree[comma1+1:comma2].split(":")[0].replace("(","").replace(")","").replace(":","").replace(" ","")           
            else:
                sample_name1 = tree[comma1+1:comma2].replace("(","").replace(")","").replace(":","").replace(" ","")
            for t in range(comma1+1, comma2):
                if sample_name1 == tree[t:t+len(sample_name1)]:
                    new_tree2 = tree[comma1+1:t]+ sample_name1 + ":" + str(leaf_length[sample_name1]) + str(tree[t+len(sample_name1):comma2])
                    break
            new_tree = new_tree + new_tree2+ ","
    new_tree = new_tree.replace(" ","").replace("[","").replace("]","").replace("'","")
    return new_tree,mutational_length,leaf_length

def get_structrue_b(b_sample_vcf,group_b,sample_b,k,tree,sample_num):
    if k == 1:
        tree=[]
    b_sample_vcf, b_same_num, b_same_pop, b_sample_name = get_gropu_same_diff_num(b_sample_vcf,group_b,sample_b)
    group_a, group_b, a_sample_vcf, b_sample_vcf,sample_a ,sample_b = split_to_2_group(b_sample_vcf, b_same_num,b_same_pop, b_sample_name)
    if len(group_a) == 1 and len(group_b) == 1:
        k = k+1
        tree_strcture = "("+str(group_a)+","+str(group_b)+")" + ":" + str(b_same_pop)
        tree.append(tree_strcture)
        return tree
    elif len(group_a) == 1 and len(group_b) > 1:
        k = k+1
        tree_strcture = "("+ str(group_a)+","+str(group_b)+")"+ ":" + str(b_same_pop)
        tree.append(tree_strcture)
        tree_next = get_structrue_b(b_sample_vcf,group_b,sample_b,k,tree,sample_num)
 
    elif len(group_a) > 1 and len(group_b) == 1:
        k = k+1
        tree_strcture = "("+str(group_a)+","+str(group_b)+")"+ ":" + str(b_same_pop)
        tree.append(tree_strcture)
        tree_next = get_structrue_b(a_sample_vcf,group_a,sample_a,k,tree,sample_num)
    else:
        k = k+1
        tree_strcture = "("+str(group_a)+","+str(group_b)+")"+ ":" + str(b_same_pop)
        tree.append(tree_strcture)
        tree_next = get_structrue_b(a_sample_vcf,group_a,sample_a,k,tree,sample_num)
        tree_next = get_structrue_b(b_sample_vcf,group_b,sample_b,k+1,tree,sample_num)
    #print(tree)
    if k == 2:
        #print(tree)
        tree = retrive_tree(tree)
        #print(tree)
        return tree
def retrive_tree(tree):
    tree_strcture = tree[-1]
    length = len(tree_strcture)
    time1= tree_strcture.split(":")[1]
    tree_strcture = tree_strcture.split(":")[0]
    #print(time1)
    tree_strcture = tree_strcture.replace("'","").replace("[","").replace("]","").replace(" ","")
    samples=[]
    numbers = tree_strcture.split(",")
    for i in range(len(numbers)):
        sample = numbers[i].replace("(","").replace(")","").replace("'","").replace("[","").replace("]","").replace(" ","")
        samples.append(sample)
    tree_strcture = tree_strcture+":"+str(int(time1))
    #print(samples)
    #print(tree_strcture)
    #print(len(tree))
    for i in range(len(tree)-1):
        i = len(tree)-i-2
        #print(tree[i])
        time2 = tree[i].split(":")[1]
        tree_strcture2 = tree[i].split(":")[0]
        new_tree_strcture2 = ""
        samples2 = []
        numbers = tree_strcture2.split(",")
        for j in range(len(numbers)):
            sample = numbers[j].replace("(","").replace(")","").replace("'","").replace("[","").replace("]","").replace(" ","")
            samples2.append(sample)
        #print(samples2)
        k = 0
        t = 1
        for j in range(len(samples)):
            if samples[j] in samples2:
                samples2.remove(samples[j])
                t = t+1
                k = 1
        #print(samples2)
        #print(k)
        if k == 1:
            #print(samples2)
            if len(samples2) == 0:
                t = t - 1
                tree_strcture =tree_strcture + ")"
                tree_strcture = tree_strcture.replace("'","")
                length = len(tree_strcture)
                #print(tree_strcture)
                for j in range(length):
                    j = len(tree_strcture)-j-1
                    if tree_strcture[j] == ",":
                        t= t -1
                    #print("t")
                    #print(t)
                    #print(tree_strcture[j])
                    if t ==1 and tree_strcture[j] == "(":
                        #print (j)
                        tree_strcture = tree_strcture[:j]+"("+tree_strcture[j:]
                        tree_strcture = tree_strcture.replace("'","")
                        t = t - 1   
            elif len(samples2) == 1:
                tree_strcture =tree_strcture + "," + samples2[0] + ")"
                tree_strcture = tree_strcture.replace("'","")
                for j in range(len(tree_strcture)):
                    j = len(tree_strcture)-j-1
                    if tree_strcture[j] == ",":
                        t= t -1
                    if t ==1 and tree_strcture[j] == "(":
                        tree_strcture = tree_strcture[:j]+"("+tree_strcture[j:]
                        tree_strcture = tree_strcture.replace("'","")
                        t = t - 1  
        if k == 0 :
            tree_strcture = tree_strcture+"," + tree_strcture2.replace("[","").replace("]","").replace("'","")
            tree_strcture = tree_strcture.replace("'","")
        #print(tree_strcture)
        #print(k)
        #print(samples2)
        tree_strcture = tree_strcture + ":" + str(int(time2))
        samples = []
        numbers = tree_strcture.split(",")
        for j in range(len(numbers)):
            p = 0
            for k in range(len(numbers[j])):
                if numbers[j][k] == ":":
                    p = k
                    break
            if p == 0:
                sample = numbers[j].replace("(","").replace(")","").replace("'","").replace("[","").replace("]","").replace(" ","")
            else:
                numbers[j] = numbers[j][:p]
                sample = numbers[j].replace("(","").replace(")","").replace("'","").replace("[","").replace("]","").replace(" ","")
            samples.append(sample)
        #print(samples)
        #print(tree_strcture)
    return tree_strcture
        
def split_to_2_group(a_sample_vcf, a_same_num,a_same_pop, a_sample_name):
    sample_a, sample_b = get_a_b_name(a_same_num)
    #print(sample_a, sample_b)
    root_length = a_same_pop
    leaf_length ,pos_sample_a, pos_sample_b = get_first_and_second_sample_leaf_length(sample_a, sample_b, a_sample_name, a_sample_vcf)
    group_a, group_b, a_sample_vcf, b_sample_vcf = get_groupa_n_groupb(sample_a, sample_b, a_sample_name, a_sample_vcf,root_length,a_same_num, pos_sample_a, pos_sample_b)
    return group_a, group_b, a_sample_vcf, b_sample_vcf,sample_a, sample_b
def get_a_b_name(a_same_num):
    sample_a = ""
    sample_b = ""
    #print(a_same_num)
    for i in a_same_num:
        if a_same_num[i] == min(a_same_num.values()):
            sample_a = i.split("_")[0]
            sample_b = i.split("_")[1]
    return sample_a, sample_b

def get_gropu_same_diff_num(a_sample_vcf,group_a,sample_a):
    same_num = {}
    sample_name = group_a
    sample_vcf = {}
    same_pop = 0
    for i in range(len(group_a)):
        #print(group_a[i])
        for j in range(i+1,len(group_a)):
            same_num[group_a[i] + "_" + group_a[j]] = 0
    for vcfs in a_sample_vcf:
        if sum(a_sample_vcf[vcfs]) == len(group_a):
            same_pop += 1
        else:
            sample_vcf[vcfs] = a_sample_vcf[vcfs]
        for i in range(len(group_a)):
            if a_sample_vcf[vcfs][i] == 1:
                for j in range(i+1,len(group_a)):
                    if a_sample_vcf[vcfs][j] == 1:
                        same_num[group_a[i] + "_" + group_a[j]] += 1
    return sample_vcf, same_num, same_pop ,sample_name
        
def get_groupa_n_groupb(sample_a, sample_b, sample_name, sample_vcf,root_length,same_num, pos_sample_a, pos_sample_b):
    groups={}
    namepos={}
    groups[sample_a] = [sample_a]
    namepos[sample_a] = [pos_sample_a]
    groups[sample_b] = [sample_b]
    namepos[sample_b] = [pos_sample_b]
    for i in range(len(sample_name)):
        if sample_name[i] != sample_a and sample_name[i] != sample_b:
            # here we need find the sample closest branch hierarchically, and add the sample to the branch
            if pos_sample_a < i:
                same_num_left = same_num[sample_a + "_" + sample_name[i]] - root_length
            else:
                same_num_left = same_num[sample_name[i] + "_" + sample_a] - root_length
            if pos_sample_b < i:
                same_num_right = same_num[sample_b + "_" + sample_name[i]] - root_length
            else:
                same_num_right = same_num[sample_name[i] + "_" + sample_b] - root_length
            if same_num_left < same_num_right:
                groups[sample_b].append(sample_name[i])
                namepos[sample_b].append(i)
            else:
                groups[sample_a].append(sample_name[i])
                namepos[sample_a].append(i)
    group_a=groups[sample_a]
    group_b=groups[sample_b]
    namepos_a=namepos[sample_a]
    namepos_b=namepos[sample_b]
    a_sample_vcf={}
    for i in sample_vcf:
        g_vcf=[]
        for j in namepos_a:
            g_vcf.append(sample_vcf[i][j])
        if sum (g_vcf) > 0:
            a_sample_vcf[i]=g_vcf
    b_sample_vcf={}
    for i in sample_vcf:
        g_vcf=[]
        for j in namepos_b:
            g_vcf.append(sample_vcf[i][j])
        if sum (g_vcf) > 0:
            b_sample_vcf[i]=g_vcf
    return group_a, group_b, a_sample_vcf, b_sample_vcf
            
def get_first_and_second_sample_leaf_length(sample_a, sample_b, sample_name, sample_vcf):  
    leaf_length = {}
    leaf_length[sample_a] = 0
    leaf_length[sample_b] = 0
    pos_sample_a = 0
    #print(sample_a)
    for i in range(len(sample_name)):
        if sample_name[i] == sample_a:
            pos_sample_a = i
            break
    for i in sample_vcf:
        if sample_vcf[i][pos_sample_a] == 1:
            leaf_length[sample_a] += 1
    pos_sample_b = 0
    for i in range(len(sample_name)):
        if sample_name[i] == sample_b:
            pos_sample_b = i
            break
    for i in sample_vcf:
        if sample_vcf[i][pos_sample_b] == 1:
            leaf_length[sample_b] += 1
    return leaf_length, pos_sample_a, pos_sample_b
            
def get_root_length(same_vcf_pop):
    root_length = 0
    for i in same_vcf_pop:
        root_length += 1
    return root_length

def get_node_list(tree_strcture, root_length):
    node_list = []
    node_list_dict = {}
    # add the root node to the node_list_dict
    node_list_dict["root"] = "root"
    for i in range(len(tree_strcture)):
        if tree_strcture[i] == "(":
            node_list.append(i)
    node_list_len = {}
    node_list_len["root"] = root_length
    for i in range(len(node_list), 0, -1):
        i -= 1
        left_pos = node_list[i]
        m = len(node_list)- i 
        for j in range(left_pos,len(tree_strcture)):
            if tree_strcture[j] == "(":
                m += 1
            elif tree_strcture[j] == ")":
                m -= 1
                if m == 0:
                    node = tree_strcture[left_pos:j+1]
                    node_list_dict[node] = node
    return node_list_dict, node_list_len
 
def get_first_and_second_sample_name(same_num):
    sample_a = ""
    sample_b = ""
    farthest_same_num = min(same_num.values())
    for key, value in same_num.items():
        if value == farthest_same_num:
            sample_a = key.split("_")[0]
            sample_b = key.split("_")[1]
            break
    return sample_a, sample_b
 
def write_tree_to_file(tree, output):
    with open(output+".tree", "w") as f:
        f.write(tree)          
def write_mutation_length_to_file(mutational_length, output):
    with open(output+".mutational_length.csv", "w") as f:
        for i in mutational_length:
            f.write(i + "\t" + str(mutational_length[i]) + "\n")
def write_leaf_length_to_file(leaf_length, output):
    with open(output+".leaf_length.csv", "w") as f:
        for i in leaf_length:
            f.write(i + "\t" + str(leaf_length[i]) + "\n")
    
def main():
    args = get_args()
    same_num, diff_num, sample_name, same_vcf, diff_vcf, same_vcf_pop, sample_vcf = count_the_same_n_diff_site_num(args.input)
    tree, mutational_length, leaf_length, = build_tree(same_num, diff_num, sample_name, same_vcf, diff_vcf, same_vcf_pop, sample_vcf)
    write_tree_to_file(tree, args.output)
    write_mutation_length_to_file(mutational_length, args.output)
    write_leaf_length_to_file(leaf_length, args.output)
    write_to_file(same_num, diff_num, args.output, sample_name)
    
    
if __name__ == "__main__":
    main()
    
