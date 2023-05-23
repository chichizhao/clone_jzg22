#!/python
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: plot the same site number between two samples 
# usage: python plot_same_site.py -i1 same_site.csv -i2 different_site.csv -o same_site.png

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def arg_parse():
    import argparse
    parser = argparse.ArgumentParser(description='plot the same site number between two samples')
    parser.add_argument('-i1', '--input', help='input file', required=True)
    input2 = parser.add_argument('-i2', '--input2', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    args = parser.parse_args()
    return args

def get_data(data1):
    data = {}
    with open(data1, 'r') as f:
        for line in f:
            if line.startswith('sample_name'):
                data['sample_name'] = line.strip().split(',')[1:]
            else:
                data[line.strip().split(',')[0]] = line.strip().split(',')[1:]
                for i in range(len(data[line.strip().split(',')[0]])):
                    data[line.strip().split(',')[0]][i] = int(data[line.strip().split(',')[0]][i])
    # we need sort the data by the sample name in ascending order
    order = []
    for i in range(len(data['sample_name'])):
        for j in range(len(data['sample_name'])):
            if i == int(data['sample_name'][j])-1:
                order.append(j)
    sorted_data = {}
    sorted_data['sample_name'] = [str(i+1) for i in range(len(order))]
    #print(sorted_data)
    for i in range(1,len(order)+1):
        sorted = []
        for j in range(len(order)):
            #print(i, order[j])
            sorted=sorted+[data[str(i)][order[j]]]
        sorted_data[str(i)] = sorted
    #print(sorted_data)        
    return sorted_data

def plot_heatmap(data1, data2, output):
    name1 = data1['sample_name']
    df1 = pd.DataFrame(data1)
    df1 = df1.set_index('sample_name')
    name2 = data2['sample_name']
    df2 = pd.DataFrame(data2)
    df2 = df2.set_index('sample_name')
    #print(name1)
    data1_value = df1.values
    data2_value = df2.values
    # Define the figure and axis objects
    fig, ax = plt.subplots()
    # Set the aspect ratio to 1 to make the circles circular
    ax.set_aspect('equal')
    # Set the axis limits
    ax.set_xlim([-1, 75])
    ax.set_ylim([-1, 75])
    max_same = 0
    min_same = 10000
    max_diff = 0
    min_diff = 10000
    max_sample = 0
    min_sample = 10000
    
    for i in range(72):
        for j in range(72):
            if i > j :
                if data2_value[i][j] > max_diff:
                    max_diff = data2_value[i][j]
                if data2_value[i][j] < min_diff:
                    min_diff = data2_value[i][j]
                color = plt.cm.Reds((data2_value[i][j]-800)/1800)
                #size = data2_value[i][j]/6000
                box = plt.Rectangle((i, j), 1, 1, facecolor=color, edgecolor='black', linewidth=0.1)
                #circle = plt.Circle((i+0.5, j+0.5), size, facecolor=color, edgecolor=color)
                ax.add_artist(box)
                #ax.add_artist(circle)
            elif i < j :
                if data1_value[i][j] > max_same:
                    max_same = data1_value[i][j]
                if data1_value[i][j] < min_same:
                    min_same = data1_value[i][j]
                color = plt.cm.Blues((data1_value[i][j]-800)/800)
                #size = data1_value[i][j]/6000
                box = plt.Rectangle((i, j), 1, 1, facecolor=color, edgecolor='black', linewidth=0.1)
                #circle = plt.Circle((i+0.5, j+0.5), size, facecolor=color, edgecolor=color)
                ax.add_artist(box)
                #ax.add_artist(circle)
            else:
                if data1_value[i][j] > max_sample:
                    max_sample = data1_value[i][j]
                if data1_value[i][j] < min_sample:
                    min_sample = data1_value[i][j]
                color = plt.cm.Greens((data1_value[i][j])/3000)
                #size = data1_value[i][j]/8000
                box = plt.Rectangle((i, j), 1, 1, facecolor=color, edgecolor='black', linewidth=0.1)
                #circle = plt.Circle((i+0.5, j+0.5), size, facecolor=color, edgecolor=color)
                ax.add_artist(box)
                #ax.add_artist(circle)
    for i in range(72):
        ax.text(i+0.5, 72.5, name1[i], ha='center', va='center', fontsize=2.5)
        ax.text(72.5, i+0.5, name1[i], ha='center', va='center', fontsize=2.5)
        box = plt.Rectangle((i, 72), 1, 1, facecolor='none', edgecolor='black', linewidth=0.1)
        ax.add_artist(box)
        box = plt.Rectangle((72, i), 1, 1, facecolor='none', edgecolor='black', linewidth=0.1)
        ax.add_artist(box)
    box = plt.Rectangle((0, 0), 73, 73, facecolor='none', edgecolor='black', linewidth=0.4)
    ax.add_artist(box)
    box = plt.Rectangle((0, 0), 72, 72, facecolor='none', edgecolor='black', linewidth=0.4)
    ax.add_artist(box)
    box = plt.Rectangle((72, 72), 1, 1, facecolor='none', edgecolor='black', linewidth=0.4)
    ax.add_artist(box)
    print(max_same, min_same, max_diff, min_diff, max_sample, min_sample)
    blues = plt.cm.Blues
    step_same = (max_same - 800)/100
    for i in range(100):
        color = plt.cm.Blues((step_same*i)/1800)
        box = plt.Rectangle((1+ i/5,73.5), 0.2, 0.5, facecolor=color, edgecolor=color)
        ax.add_artist(box)
    box = plt.Rectangle((0.9,73.4), 20.2, 0.7, facecolor='none', edgecolor='grey', linewidth=0.4)
    ax.add_artist(box)
    box = plt.Rectangle((0.9,73.4), 0.1, 0.9, facecolor='black', edgecolor='black', linewidth=0.2)
    ax.add_artist(box)
    ax.text(0.9, 74.6, '800', ha='left', fontsize=5)
    box = plt.Rectangle((21,73.4), 0.1, 0.9, facecolor='black', edgecolor='black', linewidth=0.2)
    ax.add_artist(box)
    ax.text(21, 74.6, max_same, ha='right', fontsize=5)
    ax.text(11, 74.6, 'Same sites', ha='center', fontsize=5)
    step_diff = (max_diff - 800)/100
    for i in range(100):
        color = plt.cm.Reds((step_diff*i)/1800)
        box = plt.Rectangle((22+ i/5,73.5), 0.2, 0.5, facecolor=color, edgecolor=color)
        ax.add_artist(box)
    box = plt.Rectangle((21.9,73.4), 20.2, 0.7, facecolor='none', edgecolor='grey', linewidth=0.4)
    ax.add_artist(box)
    box = plt.Rectangle((21.9,73.4), 0.1, 0.9, facecolor='black', edgecolor='black', linewidth=0.2)
    ax.add_artist(box)
    ax.text(21.9, 74.6, '800', ha='left', fontsize=5)
    box = plt.Rectangle((42,73.4), 0.1, 0.9, facecolor='black', edgecolor='black', linewidth=0.2)
    ax.add_artist(box)
    ax.text(42, 74.6, max_diff, ha='right', fontsize=5)
    ax.text(32, 74.6, 'Different sites', ha='center', fontsize=5)
    step_sample = (max_sample)/100
    for i in range(100):
        color = plt.cm.Greens((step_sample*i)/3000)
        box = plt.Rectangle((43+ i/5,73.5), 0.2, 0.5, facecolor=color, edgecolor=color)
        ax.add_artist(box)
    box = plt.Rectangle((42.9,73.4), 20.2, 0.7, facecolor='none', edgecolor='grey', linewidth=0.4)
    ax.add_artist(box)
    box = plt.Rectangle((42.9,73.4), 0.1, 0.9, facecolor='black', edgecolor='black', linewidth=0.2)
    ax.add_artist(box)
    ax.text(42.9, 74.6, '0', ha='left', fontsize=5)
    box = plt.Rectangle((63,73.4), 0.1, 0.9, facecolor='black', edgecolor='black', linewidth=0.2)
    ax.add_artist(box)
    ax.text(63, 74.6, max_sample, ha='right', fontsize=5)
    ax.text(53, 74.6, 'Samples variants', ha='center', fontsize=5)
    
    ax.text(68, 73.5, 'Samples', ha='center', fontsize=8)
    
           
    
    plt.tight_layout()  
    ax.axis('off')
    plt.savefig('same_site2.png', dpi=600)

def main():
    args = arg_parse()
    same_site = get_data(args.input)
    same_site2 = get_data(args.input2)
    plot_heatmap(data1=same_site, data2=same_site2, output=args.output)

if __name__ == '__main__':
    main()