#!/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: calculate the distance metrix between the samples
# the distance metrix is the distance between the samples
x=[103.763654,103.764301,103.764889,103.764608,103.764502,103.764288,103.763775,103.763702,103.763756,103.763755,103.763725,103.760635,103.760357,103.757515,103.755514,103.758088,103.757893,103.758538,103.758825,103.758796,103.758735,103.758984,103.759034,103.7591,103.759195,103.759447,103.75934,103.759392,103.759477,103.759592,103.759652,103.759689,103.759797,103.759816,103.759765,103.759335,103.760224,103.761106,103.761708,103.75755,103.758766,103.758753,103.758743,103.758949,103.759238,103.759181,103.759483,103.759584,103.759651,103.761979,103.760714,103.759915,103.759687,103.759427,103.759348,103.759335,103.759199,103.759423,103.759511,103.763726,103.763759,103.754918,103.758533,103.758205,103.758118,103.757918,103.75687,103.756137,103.756045,103.755898,103.756904,103.756943]
y=[33.261663,33.262041,33.264221,33.26095,33.260906,33.260767,33.260143,33.25943,33.258929,33.257167,33.256892,33.249185,33.249178,33.249165,33.248849,33.249935,33.250169,33.250454,33.25082,33.251882,33.252507,33.25319,33.25333,33.253459,33.253657,33.253664,33.253872,33.254014,33.25416,33.254221,33.254389,33.254553,33.254757,33.254989,33.255163,33.256076,33.259091,33.259659,33.260041,33.250075,33.250185,33.252133,33.252726,33.252911,33.253307,33.25318,33.253822,33.253931,33.254038,33.260085,33.259293,33.258951,33.258801,33.258327,33.257764,33.257608,33.256847,33.25587,33.255692,33.259084,33.258738,33.248196,33.249152,33.250251,33.24941,33.249793,33.249847,33.249668,33.248747,33.24923,33.248951,33.249454]
import matplotlib.pyplot as plt
import numpy as np
distance_metrix = np.zeros((len(x),len(x)))
for i in range(len(x)):
    for j in range(len(x)):
        # here the x is the longitude and the y is the latitude
        # so we need to calculate the distance between the two points but the uint is meter
        distance_metrix[i][j] =np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)*111000
tree ='((((((((((72:521.12,34:521.12):172.29,55:693.41):70.79,71:764.2):111.17,(((3:458.76,36:458.76):174.55,32:633.31):106.08,49:739.38):135.99):37.4,29:912.77):87.38,(((((47:515.21,30:515.21):136.87,51:652.08):95.29,57:747.37):91.42,(45:525.69,59:525.69):313.1):48.32,28:887.11):113.03):13.17,21:1013.31):26.93,22:1040.24):154.39,((((((69:574.84,31:574.84):459.1,((((53:594.47,18:594.47):331.69,(((27:552.52,6:552.52):168.05,66:720.57):144.28,(41:628.39,39:628.39):236.46):61.31):32.83,(68:584.6,63:584.6):374.39):33.93,(54:598.59,7:598.59):394.33):41.01):14.27,61:1048.2):102.59,(((((((8:524.89,43:524.89):150.61,9:675.5):144.21,(1:539.65,10:539.65):280.07):44.24,2:863.95):46.38,64:910.33):62.57,((67:559.99,50:559.99):174.67,40:734.66):238.24):140.86,(((65:592.68,13:592.68):168.93,16:761.61):317.78,((((((20:512.53,35:512.53):306.92,((58:520.8,37:520.8):165.9,46:686.7):132.74):123.62,((((25:523.09,70:523.09):155.09,17:678.18):75.83,52:754.01):55.11,19:809.12):133.94):50.56,(((((15:476.42,14:476.42):153.85,24:630.26):93.43,60:723.69):62.71,48:786.4):48.26,44:834.66):158.96):9.99,((12:567.32,26:567.32):147.42,56:714.74):288.87):12.37,23:1015.98):63.41):34.37):37.04):5.72,(((33:559.48,11:559.48):175.55,4:735.03):106.84,42:841.87):314.65):6.57,((5:582.23,38:582.23):165.17,62:747.41):415.68):31.55):410.64)'
import re
divsion_time_metrix = np.zeros((len(x),len(x)))
tree_sample = tree.split(',')
def get_sample(tree):
    sample = 0
    sample= int(tree.split(':')[0].replace('(','').replace(')',''))
    return sample
samples =[]
for i in range(len(tree_sample)):
    samples.append(get_sample(tree_sample[i]))
def get_division_time(tree):
    division_time = 0
    for i in range(len(tree)):
        i = len(tree)-i-1
        if tree[i] == ':':
            division_time += float(tree[i:len(tree)].split(':')[1].replace(')',''))
        elif tree[i] == ',':
            break
    return division_time



for i in range(len(x)):
    for j in range(len(x)):
        # here the x is the longitude and the y is the latitude
        # so we need to calculate the distance between the two points but the uint is meter
        sample1 = i
        sample2 = j
        p1 = 0
        p2 = 0
        for k in range(len(samples)):
            if samples[k] == sample1:
                p1 = k
            if samples[k] == sample2:
                p2 = k
        if p1 == p2:
            divsion_time_metrix[i][j] = 0
        elif p1 < p2:
            new_tree =''
            start1 = 0
            start2 = 0
            for k in range(p1+1):
                start1 = len(str(tree_sample[k]))+1+start1
            for k in range(p2):
                start2 = len(str(tree_sample[k]))+1+start2
            for k in range(start2,len(tree)):
                if tree[k] == ')':
                    # here we chech this struture if complete
                    node = 0
                    test = 0 
                    check_point = k
                    for l in range(k):
                        l = k-l
                        if tree[l] == '(':
                            test = 1
                            node -= 1
                        elif tree[l] == ')':
                            test = 1
                            node += 1
                        if node == 0 and test == 1:
                            check_point = l
                            break
                    if check_point<start1:
                        new_tree = tree[check_point+1:k+1]
                        break
            divsion_time_metrix[i][j] = get_division_time(new_tree)
            divsion_time_metrix[j][i] = get_division_time(new_tree)

from math import log 
expand = np.zeros((len(x),len(x)))
for i in range(len(x)):
    for j in range(len(x)):
        if divsion_time_metrix[i][j] == 0:
            expand[i][j] = 0
        else:
            expand[i][j] = distance_metrix[i][j]/divsion_time_metrix[i][j]/2

fig, ax = plt.subplots(figsize=(20, 20))
ax.set_aspect('equal')
ax.set_xlim(-2, 75)
ax.set_ylim(-2, 85)
# here we plot the distance metrix and the expnading matrix into one plot
box = plt.Rectangle((0, 76), 72, 3, facecolor='none', edgecolor='black', linewidth=0.8, alpha=0.1)
ax.add_artist(box)
box = plt.Rectangle((0, 79), 72, 3, facecolor='none', edgecolor='black', linewidth=0.8, alpha=0.1)
ax.add_artist(box)
max_distance = 0
max_expand = 0
for i in range(len(x)):
    for j in range(len(x)):
        if i> j :
            if distance_metrix[i][j] > max_distance:
                max_distance = distance_metrix[i][j]
            color = plt.cm.Blues(((log(distance_metrix[i][j])/log(10))/6))
            size = float((log(distance_metrix[i][j])/log(10))/6)
            box = plt.Rectangle((i, j), 1, 1, facecolor=color, edgecolor='black', linewidth=0.2)
            #circle = plt.Circle((i+0.5, j+0.5), size, facecolor=color, edgecolor='black', linewidth=0.1)
            ax.add_artist(box)
            #ax.add_artist(circle)
        elif i< j:
            if expand[i][j] > max_expand:
                max_expand = expand[i][j]
            color = plt.cm.Reds((expand[i][j]/2))
            size = expand[i][j]/2
            box = plt.Rectangle((i, j), 1, 1, facecolor=color, edgecolor='black', linewidth=0.2)
            #circle = plt.Circle((i+0.5, j+0.5), size, facecolor=color, edgecolor='black', linewidth=0.1)
            ax.add_artist(box)
            #ax.add_artist(circle)

    #box = plt.Rectangle((i, i), 1, 11, facecolor='none', edgecolor='black', linewidth=0.1)
    #ax.add_artist(box)
    # here we want add a boxplot of each column into the box
    box_data = []
    for j in range(len(x)):
        if i != j:
            box_data.append(expand[i][j]*8+73.1)
    ax.boxplot(box_data, positions=[i+0.5], widths=0.7, showfliers=False)
    box = plt.Rectangle((i, 72), 1, 1, facecolor='none', edgecolor='black', linewidth=0.2)
    ax.add_artist(box)
    ax.text(i+0.5, 72.5, str(i+1), fontsize=8, horizontalalignment='center', verticalalignment='center')
    box = plt.Rectangle((0,72), 72, 1, facecolor='none', edgecolor='black', linewidth=0.4)
    ax.add_artist(box)
    box = plt.Rectangle((-1, i), 1, 1, facecolor='none', edgecolor='black', linewidth=0.2)
    ax.add_artist(box)
    ax.text(-0.5, i+0.5, str(i+1), fontsize=8, horizontalalignment='center', verticalalignment='center')
box = plt.Rectangle((0, 0), 72, 72, facecolor='none', edgecolor='black', linewidth=0.8)
ax.add_artist(box)
box = plt.Rectangle((0, 72), 72, 12, facecolor='none', edgecolor='black', linewidth=0.8)
ax.add_artist(box)
box = plt.Rectangle((-1, 0), 1, 84, facecolor='none', edgecolor='black', linewidth=0.8)
ax.add_artist(box)
max_value = np.max(expand)
box = plt.Rectangle((0,79), 0.2, 0.1, facecolor='black', edgecolor='black', linewidth=0.8)
ax.add_artist(box)
box = plt.Rectangle((0,76), 0.2, 0.1, facecolor='black', edgecolor='black', linewidth=0.8)
ax.add_artist(box)
box = plt.Rectangle((0,73), 0.2, 0.1, facecolor='black', edgecolor='black', linewidth=0.8)
ax.add_artist(box)
box = plt.Rectangle((0,82), 0.2, 0.1, facecolor='black', edgecolor='black', linewidth=0.8)
ax.add_artist(box)
ax.text(-0.55,82.1, str('1.5'), fontsize=8, horizontalalignment='center', verticalalignment='center')
ax.text(-0.55,79.1, str('1'), fontsize=8, horizontalalignment='center', verticalalignment='center')
ax.text(-0.55,76.1, str('0.5'), fontsize=8, horizontalalignment='center', verticalalignment='center')
ax.text(-0.55,73.1, str('0'), fontsize=8, horizontalalignment='center', verticalalignment='center')
ax.text( 36, 83, str('Expanding rate (m/year)'), fontsize=16, horizontalalignment='center', verticalalignment='center')
ax.text(2,-1.5, str('Sample'), fontsize=16, horizontalalignment='center', verticalalignment='center')
for i in range(99):
    box = plt.Rectangle((45+i/5,-1.0 ), 0.5, 0.5, facecolor=plt.cm.Blues((log(max_distance)/log(10))/6*i/100), edgecolor=plt.cm.Blues((log(max_distance)/log(10))/6*i/100), linewidth=0.2)
    ax.add_artist(box)
    box = plt.Rectangle((10+i/5,-1.0 ), 0.5, 0.5, facecolor=plt.cm.Reds((max_expand)/2*i/100), edgecolor=plt.cm.Reds((max_expand)/2*i/100), linewidth=0.2)
    ax.add_artist(box)
box = plt.Rectangle((45,-1.0 ),20, 0.5, facecolor='none', edgecolor='black', linewidth=0.2)
ax.add_artist(box)
box = plt.Rectangle((10,-1.0 ),20, 0.5, facecolor='none', edgecolor='black', linewidth=0.2)
ax.add_artist(box)
box = plt.Rectangle((45,-1.2 ),0.05, 0.7, facecolor='black', edgecolor='black', linewidth=0.2)
ax.add_artist(box)
box = plt.Rectangle((10,-1.2 ),0.05, 0.7, facecolor='black', edgecolor='black', linewidth=0.2)
ax.add_artist(box)
box = plt.Rectangle((65,-1.2 ),0.05, 0.7, facecolor='black', edgecolor='black', linewidth=0.2)
ax.add_artist(box)
box = plt.Rectangle((30,-1.2 ),0.05, 0.7, facecolor='black', edgecolor='black', linewidth=0.2)
ax.add_artist(box)
ax.text( 55.5, -2.0, str('Distance (m)'), fontsize=16, horizontalalignment='center', verticalalignment='center')
ax.text( 20.5, -2.0, str('Expanding rate (m/year)'), fontsize=16, horizontalalignment='center', verticalalignment='center')
ax.text( 10, -2, str('0'), fontsize=10, horizontalalignment='center', verticalalignment='center')
ax.text( 30, -2., str(format(max_expand, '.2f')), fontsize=14, horizontalalignment='center', verticalalignment='center')
ax.text( 45, -2, str('0'), fontsize=10, horizontalalignment='center', verticalalignment='center')
ax.text( 65, -2, str(format(max_distance, '.0f')), fontsize=14, horizontalalignment='center', verticalalignment='center')
ax.axis('off')

plt.savefig('expand.png', dpi=600, bbox_inches='tight')
