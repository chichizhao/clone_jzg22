#!/python3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: plot the phylogenetic tree of jzg population
# target: plot the phylogenetic tree of jzg population
# 1. plot the time tree and the branch length is the mutation sites
# 2. give the branch a color according to the mutation tree. the colors are based on the mutation rate of each branch
# 3. add a circle to the branch and the size of the circle is the non-synonymous mutation number of each branch

from matplotlib import pyplot as plt
import numpy as np

def get_samples(node):
    samples = []
    #print(node)
    node = node.strip().split(',')
    for i in node:
        if i[0] == '(':
            samples.append(i.strip().split(':')[0].split('(')[-1])
        else:
            samples.append(i.strip().split(':')[0])
    #print(samples)
    return samples 

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
    #print(get_samples(new_tree))
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

def get_node_time(tree,node):
    # total time is 2434
    #print(node)
    node_time = 0
    # find the node position in the tree
    node_pos = 0
    length_node = len(node)
    for i in range(len(tree)):
        if str(tree[i:i+length_node]) == str(node):
            node_pos = i
            break
    print(node)
    print(tree)
    #print(node_pos)
    branch_length = 0

    # calculate the length to the very outside
    distance2out = 0
    k = 1
    m = 0
    n = 0
    
    for i in range(node_pos+len(node)-1,len(tree)):
        if tree[i] == ')':
            if n > 0:
                n -= 1
            elif n == 0 and k == 0:
                k = 1
        elif tree[i] == '(':
            n += 1
            k = 0
        elif tree[i] == ',':
            m += 1
        elif tree[i] == ':' and n == 0 and k == 1 and tree[i-1] == ')':
            for j in range(i,len(tree)):
                if tree[j] == ',' or tree[j] == ')':
                    distance2out = distance2out + float(tree[i+1:j])
                    break
        #print(distance2out)
        
    for i in range(len(node)):
        i = len(node) - 1 - i
        if node[i] == ':':
            print("ture")
            print(i)
            print(node[i+1:len(node)])
            branch_length = float(node[i+1:len(node)])

            break
    #
    print(distance2out)
    print(branch_length)
    print(2434-distance2out-branch_length,2434-distance2out)
    return 2434-distance2out-branch_length,2434-distance2out        

time_tree='(((((((((49:614,55:582):311,(71:555,52:718):297):92,(32:567,3:764):413):54,29:755):154,(((((37:571,47:741):257,(51:544,30:796):278):149,(((36:765,34:861):295,72:781):144,45:806):178):43,57:1012):45,28:753):118):16,21:785):33,(22:306,59:870):426):176,((((((69:570,1:757):197,31:743):348,(((41:628,27:648):173,68:753):267,((((18:617,33:855):218,6:1009):123,53:906):136,((63:705,39:1022):273,4:695):280):62):68):26,(61:613,38:601):486):115,((((((10:691,50:748):222,64:810):91,5:863):279,((23:464,17:917):521,(((((70:680,14:945):249,46:861):198,(19:639,44:692):319):173,((((60:737,24:744):222,20:829):104,35:832):118,(58:569,15:804):379):140):23,(12:632,56:613):448):29):91):57,((((11:720,8:890):253,66:787):241,((67:662,13:614):185,54:803):182):89,(65:616,16:727):390):152):20,(((((7:783,43:841):360,(25:903,26:836):385):92,48:993):154,((2:612,9:782):202,40:678):254):35,42:856):207):46):2,62:980):46):525)'
target_mark_tree = '((((49:614,55:582):311,(71:555,52:718):297):92,(32:567,3:764):413):54,29:755):154'
G1 = [1,2,3,4,5,6,7,8,9,60,61,10,11]
G2 = [50,39,38,51,37,52,53,54,55,56,57]
G3 = [36,58,59,35,34,33,32,31,30,29,28,27,25,26,24,23,22,49,48,47,26,45,46,44,43,21,42,20]
G4 = [12,13,19,18,41,64,17,16,40,66,65,63,67,72,71,14,68,70,69,15,62]
group = {}
group[1]= G1
group[2]= G2
group[3]= G3
group[4]= G4
# here we need four color to color the group 
group_color = ['#dcd116', '#6495ED', '#32CD32', '#b627b1']
#
# group_color = ['None', '#6495ED', 'None', 'None']
# group_color = ['None', 'None', '#32CD32', 'None']
# group_color = ['None', 'None', 'None', '#b627b1']
#group_color = ['#dcd116', 'None', 'None', 'None']
non_synonymous ={'1944_2434': 2.0, '1863_1909': 1.0, '1459_1738': 2.0, '1368_1459': 1.0, '1146_1368': 2.0, '455_1146': 5.0, '1643_1795': 1.0, '1313_1554': 3.0, '1060_1313': 2.0, '340_1060': 4.0, '1647_1738': 4.0, '1170_1618': 2.0, '538_1170': 6.0, '1372_1554': 2.0, '1187_1372': 2.0, '573_1187': 5.0, '1224_1422': 1.0, '975_1224': 1.0, '30_975': 4.0, '1076_1455': 5.0, '272_1076': 8.0, '1253_1643': 4.0, '526_1253': 8.0, '1126_1647': 2.0, '209_1126': 5.0, '1746_1861': 1.0, '1454_1590': 1.0, '1331_1454': 1.0, '1113_1331': 1.0, '496_1113': 4.0, '1103_1422': 7.0, '464_1103': 2.0, '1372_1720': 2.0, '418_1175': 8.0, '404_1233': 7.0, '1733_1909': 4.0, '1700_1733': 1.0, '915_1700': 7.0, '1307_1733': 2.0, '662_1126': 3.0, '1011_1233': 1.0, '267_1011': 4.0, '1608_1815': 1.0, '1419_1573': 1.0, '942_1327': 7.0, '39_942': 1.0, '106_942': 5.0, '1385_1652': 1.0, '564_1212': 4.0, '813_1566': 8.0, '1530_1684': 1.0, '775_1530': 6.0, '1329_1478': 1.0, '1051_1329': 1.0, '255_1051': 6.0, '1319_1573': 1.0, '1117_1319': 1.0, '505_1117': 3.0, '629_1372': 4.0, '1063_1476': 5.0, '496_1063': 4.0, '258_1113': 10.0, '1300_1478': 1.0, '1156_1300': 1.0, '861_1156': 2.0, '0_861_34': 7, '505_1337': 6.0, '96_861': 5.0, '1072_1329': 4.0, '501_1072': 2.0, '1260_1746': 1.0, '659_1260': 5.0, '1310_1590': 1.0, '1037_1310': 4.0, '15_1037': 7.0, '641_1319': 7.0, '299_1063': 5.0, '584_1212': 11.0, '752_1608': 7.0, '967_1327': 1.0, '126_967': 1.0, '411_1103': 7.0, '494_1300': 9.0, '363_1224': 6.0, '331_1072': 6.0, '426_1419': 13.0, '1073_1384': 2.0, '459_1073': 2.0, '615_1310': 6.0, '398_1146': 2.0, '507_1051': 1.0, '1087_1384': 3.0, '369_1087': 4.0, '548_1454': 3.0, '569_1372': 6.0, '491_1073': 5.0, '557_1170': 3.0, '509_1521': 5.0, '507_1076': 5.0, '437_1307': 8.0, '596_1459': 11.0, '274_1011': 4.0, '647_1260': 4.0, '883_1863': 9.0, '322_1331': 9.0, '332_1037': 3.0, '558_1368': 10.0, '637_1253': 3.0, '526_1313': 6.0, '525_1187': 2.0, '632_1385': 4.0, '605_1175': 6.0, '295_975': 10.0, '184_967': 9.0, '532_1087': 2.0, '375_1156': 11.0, '170_1060': 5.0, '335_1117': 7.0}

fig, ax = plt.subplots(figsize=(6.5, 6))
ax.set_xlim(-1, 24.5)
ax.set_ylim(0, 220)
#ax.set_xticks([])
#ax.set_yticks([])
#ax.set_title('Tree', fontsize=20)
num = 0
x_value = [3]
x2_value = []
num_count = 0
for i in range(len(time_tree)):
    if time_tree[i] == ':':
        num_count += 1
        sample = ''
        value = ''
        start = 0
        for j in range(i):
            j = i - j - 1
            if time_tree[j] == ',':
                sample = time_tree[j+1:i]
                for k in range(i, len(time_tree)):
                    if time_tree[k] == ')' or time_tree[k] == ',':
                        value = float(time_tree[i+1:k])
                        break
                break
            elif time_tree[j] == ')':
                sample = ''
                for k in range(i, len(time_tree)):
                    if time_tree[k] == ',' or time_tree[k] == ')':
                        value = float(time_tree[i+1:k])
                        break
                break
            elif time_tree[j] == '(':
                sample = time_tree[j+1:i]
                #print(sample)
                for k in range(i, len(time_tree)):
                    if time_tree[k] == ')' or time_tree[k] ==',':
                        value = float(time_tree[i+1:k])
                        break
                break
        for j in range(i):
            j = i - j - 1
            if time_tree[j] == ',' or time_tree[j] == '(':
                # start = 0
                type = 0
                # add a x_value[-1]+3 to x_value
                x_value.append(x_value[-1]+3)
                #print(x_value)
                for k in range(i, len(time_tree)):
                    if time_tree[k] == ')' or time_tree[k] == ',':
                        break
                node_tree = time_tree[j+1:k]
                start, end = get_node_time(time_tree,node_tree)
                #print(value)
                #print(node_tree)
                #print(time_tree)
                #print(start,end)
                break
            
            
            
            elif time_tree[j] == ')':
                type = 1
                tree = ''
                node_right = 1
                node_left=0
                for k in range(j):
                    k = j-k-1
                    if time_tree[k] == ')':
                        node_right +=1
                    elif time_tree[k] == '(':
                        node_left +=1
                    if node_left == node_right:
                        break
                tree = time_tree[k:i]
                for m in range(i,len(time_tree)):
                    if time_tree[m] == ',' or time_tree[m] == ')':
                        break
                node_tree = time_tree[k:m]
                sample_num = 0
                sample_num = len(get_samples(tree))
                if sample_num <= 2:
                    x3_start = x_value[-2]
                    x3_width = x_value[-1] - x_value[-2]+0.1
                    x2_value.append((x_value[-1]+x_value[-2])/2)
                elif sample_num > 2:
                    x2_tree_right = 0
                    for k in range(len(tree)):
                        k = len(tree) - k - 1
                        if tree[k] == ',':
                            break
                        elif tree[k] == ')': 
                            x2_tree_right += 1
                    if x2_tree_right == 1:
                        x3_start = x2_value[-1]
                        x3_width = x_value[-1] - x2_value[-1] +0.1
                        x2_value[-1] = (x2_value[-1]+x_value[-1])/2
                    else:
                        x3_start = x2_value[-2]
                        x3_width = x2_value[-1] - x2_value[-2] +0.2
                        x2_value[-2] = (x2_value[-2]+x2_value[-1])/2
                        x2_value = x2_value[:-1]

                start , end = get_node_time(time_tree,node_tree)
                print(value)
                print(node_tree)
                print(time_tree)
                print(start,end) 
                break
        if str(node_tree) in str(target_mark_tree):
            target = 1
        else:
            target = 0 

        if type == 0:
            color = 'black'
            for k in range(1,5):
                if int(sample) in group[k]:
                    color = group_color[k-1]
                    break
            if color == 'black':
                continue
                #box = plt.Rectangle((start,x_value[-1]), value/100, 0.1, facecolor='black', edgecolor='black', lw=0.1)
                #ax.add_artist(box)
            else:

                box = plt.Rectangle((-0.8,x_value[-1]-1.5), 1, 3, facecolor=color, edgecolor='none', lw=0.1)
                ax.add_artist(box)
                box = plt.Rectangle((start/100,x_value[-1]), value/100, 0.1, facecolor='black', edgecolor='black', lw=0.5)
                ax.add_artist(box)
                box = plt.Rectangle((start/100,x_value[-1]-0.1), 0.05, 0.4, facecolor='black', edgecolor='black', lw=0.1)
                ax.add_artist(box)
                box = plt.Rectangle((0,x_value[-1]), start/100, 0.1, facecolor='black', edgecolor='black', lw=0.1, alpha=0.1)
                ax.add_artist(box)
                if target == 1:
                    box = plt.Rectangle((start/100,x_value[-1]), value/100, 0.1, facecolor='#0ED3FA', edgecolor='#0ED3FA', lw=0.5)
                    ax.add_artist(box)


            for key in non_synonymous:
                #print(key)
                if int(key.split('_')[0]) == start and int(key.split('_')[1]) == end:
                    for i in range(int(non_synonymous[key])):
                        if i%2 == 0:
                            box = plt.Rectangle(((start+value/2)/100-0.1+i//2* 0.25,x_value[-1]+0.2), 0.2, 1, facecolor='#B5E4D0', edgecolor='grey', lw=0.1)
                            ax.add_artist(box)
                        else:
                            box = plt.Rectangle(((start+value/2)/100-0.1+i//2* 0.25,x_value[-1]+1.4), 0.2, 1, facecolor='#F3AFA7', edgecolor='grey', lw=0.1)
                            ax.add_artist(box)
            # add sample name background box
            box = plt.Rectangle((-0.5,x_value[-1]-1.25), 0.5, 2.5, facecolor='white', alpha =0.5)
            ax.add_artist(box)
            # add sample name
            text = plt.text(-0.5-0.076, x_value[-1]-0.5-0.95, str(sample), fontsize=5,fontname= 'Arial')
            ax.add_artist(text)
            # add sample node
            box = plt.Rectangle((0,x_value[-1]-0.3),0.05,0.6,facecolor='black',)
            ax.add_artist(box)
        elif start != 0:
            box2 = plt.Rectangle(((start)/100,x2_value[-1]), value/100, 0.2, facecolor='black', edgecolor='black', lw=0.5)
            ax.add_artist(box2)
            box = plt.Rectangle(((start)/100,x3_start), 0.01, x3_width, facecolor='black', edgecolor='black', lw=0.1)
            ax.add_artist(box)
            # add a node symbol
            box= plt.Rectangle(((start-2)/100,x2_value[-1]-0.2), 0.05, 0.6, facecolor='black', edgecolor='black', lw=0.1)
            ax.add_artist(box)
            if target == 1:
                box = plt.Rectangle(((start)/100,x2_value[-1]), value/100, 0.2, facecolor='#0ED3FA', edgecolor='#0ED3FA', lw=0.5)
                ax.add_artist(box)
                box = plt.Rectangle(((start)/100,x3_start), 0.01, x3_width, facecolor='#0ED3FA', edgecolor='#0ED3FA', lw=0.1)
                ax.add_artist(box)
                box= plt.Rectangle(((start-2)/100,x2_value[-1]-0.2), 0.05, 0.6, facecolor='#0ED3FA', edgecolor='#0ED3FA', lw=0.1)
                ax.add_artist(box)
            for key in non_synonymous:
                if int(key.split('_')[0]) == int(start) and int(key.split('_')[1]) == int(value+start):
                    for i in range(int(non_synonymous[key])):
                        if i%2 == 0:
                            box = plt.Rectangle(((start+value/2)/100-0.1+i//2* 0.25,x2_value[-1]+0.2), 0.2, 1, facecolor='#B5E4D0', edgecolor='grey', lw=0.1)
                            ax.add_artist(box)
                        else:
                            box = plt.Rectangle(((start+value/2)/100-0.1+i//2* 0.25,x2_value[-1]+1.4), 0.2, 1, facecolor='#F3AFA7', edgecolor='grey', lw=0.1)
                            ax.add_artist(box)
# add the root box
for key in non_synonymous:
    if int(key.split('_')[0]) == int(1944) and int(key.split('_')[1]) == int(2434):
        for i in range(int(non_synonymous[key])):
            if i%2 == 0:
                box = plt.Rectangle(((1944+2434)/200-0.1+i//2* 0.25,38 * 3+0.2), 0.2, 1, facecolor='#B5E4D0', edgecolor='grey', lw=0.1)
                ax.add_artist(box)
            else:
                box = plt.Rectangle(((1944+2434)/200-0.1+i//2* 0.25,38 * 3+1.4), 0.2, 1, facecolor='#F3AFA7', edgecolor='grey', lw=0.1)
                ax.add_artist(box)
            
num = 0
x_value = [3]
x2_value = []
num_count = 0 

# add variants scale
t = (2500-2434)/100
box = plt.Rectangle((0-t,4.4),25,0.4,facecolor='black')
ax.add_artist(box)
box = plt.Rectangle((0-t,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(0-t, -0.5, str(2500), fontsize=5,fontname= 'Arial',ha ='center')
ax.add_artist(text)

box = plt.Rectangle((5-t,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(5-t, -0.5, str(2000), fontsize=5,fontname= 'Arial',ha ='center')
ax.add_artist(text)

box = plt.Rectangle((10-t,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(10-t, -0.5, str(1500), fontsize=5,fontname= 'Arial',ha ='center')
ax.add_artist(text)
box = plt.Rectangle((15-t,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(15-t, -0.5, str(1000), fontsize=5,fontname= 'Arial',ha ='center')
ax.add_artist(text)
box = plt.Rectangle((20-t,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(20-t, -0.5, str(500), fontsize=5,fontname= 'Arial',ha ='center')
ax.add_artist(text)
box = plt.Rectangle((24.95-t,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(25-t, -0.5, str(0), fontsize=5,fontname= 'Arial',ha ='center')
ax.add_artist(text)


ax.add_artist(box)

# add the x axies mark
text = plt.text((25-t)/2, -3, "Variants", fontsize=6,fontname= 'Arial',ha ='center')
# off the axis
ax.axis('off')
# add the Nonsynonymous mutation
box = plt.Rectangle((21.5,15),0.5,4,facecolor='#B5E4D0')
ax.add_artist(box)
box = plt.Rectangle((21.5,20),0.5,4,facecolor='#F3AFA7')
ax.add_artist(box)
box = plt.Rectangle((22.1,15),0.5,4,facecolor='#B5E4D0')
ax.add_artist(box)

text = plt.text(18, 17, "Non-synonymous \n Mutations", fontsize=6,fontname= 'Arial',ha ='left')
#plt.show()
plt.savefig('tree.pdf', dpi=300,bbox_inches='tight')
