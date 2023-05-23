#!/python3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: plot the phylogenetic tree of jzg population
# time_tree =((((((((((72:521.12,34:521.12):172.29,55:693.41):70.79,71:764.2):111.17,(((3:458.76,36:458.76):174.55,32:633.31):106.08,49:739.38):135.99):37.4,29:912.77):87.38,(((((47:515.21,30:515.21):136.87,51:652.08):95.29,57:747.37):91.42,(45:525.69,59:525.69):313.1):48.32,28:887.11):113.03):13.17,21:1013.31):26.93,22:1040.24):154.39,((((((69:574.84,31:574.84):459.1,((((53:594.47,18:594.47):331.69,(((27:552.52,6:552.52):168.05,66:720.57):144.28,(41:628.39,39:628.39):236.46):61.31):32.83,(68:584.6,63:584.6):374.39):33.93,(54:598.59,7:598.59):394.33):41.01):14.27,61:1048.2):102.59,(((((((8:524.89,43:524.89):150.61,9:675.5):144.21,(1:539.65,10:539.65):280.07):44.24,2:863.95):46.38,64:910.33):62.57,((67:559.99,50:559.99):174.67,40:734.66):238.24):140.86,(((65:592.68,13:592.68):168.93,16:761.61):317.78,((((((20:512.53,35:512.53):306.92,((58:520.8,37:520.8):165.9,46:686.7):132.74):123.62,((((25:523.09,70:523.09):155.09,17:678.18):75.83,52:754.01):55.11,19:809.12):133.94):50.56,(((((15:476.42,14:476.42):153.85,24:630.26):93.43,60:723.69):62.71,48:786.4):48.26,44:834.66):158.96):9.99,((12:567.32,26:567.32):147.42,56:714.74):288.87):12.37,23:1015.98):63.41):34.37):37.04):5.72,(((33:559.48,11:559.48):175.55,4:735.03):106.84,42:841.87):314.65):6.57,((5:582.23,38:582.23):165.17,62:747.41):415.68):31.55):410.64)
# mutation_tree=((((((((((72:5.18,34:8.94):7.06,55:5.02):6.38,71:4.9):6.01,(((3:6.2,36:7.8):7.0,32:5.04):6.34,49:5.24):6.07):6.04,29:3.91):5.8,(((((47:6.01,30:6.85):6.43,51:4.98):5.95,57:5.14):5.75,(45:5.09,59:5.69):5.39):5.63,28:3.63):5.34):5.6,21:3.48):5.48,22:3.1):5.34,((((((69:4.92,31:4.68):4.8,((((53:5.02,18:5.34):5.18,(((27:4.48,6:6.77):5.62,66:4.95):5.4,(41:4.09,39:8.66):6.37):5.79):5.62,(68:4.2,63:6.61):5.4):5.57,(54:4.3,7:7.38):5.84):5.62):5.49,61:4.7):5.44,(((((((8:7.14,43:7.24):7.19,9:5.85):6.74,(1:5.71,10:5.31):5.51):6.25,2:5.12):6.06,64:4.88):5.89,((67:5.05,50:5.71):5.38,40:4.27):5.01):5.63,(((65:4.73,13:4.99):4.86,16:5.17):4.96,((((((20:5.77,35:5.01):5.39,((58:5.12,37:5.04):5.08,46:6.08):5.42):5.41,((((25:7.82,70:5.97):6.9,17:6.63):6.81,52:5.67):6.52,19:5.22):6.26):5.83,(((((15:5.93,14:8.33):7.13,24:6.31):6.86,60:6.18):6.69,48:5.24):6.4,44:5.58):6.26):6.0,((12:4.35,26:8.21):6.28,56:4.76):5.78):5.96,23:4.5):5.89):5.77):5.72):5.64,(((33:6.64,11:6.28):6.46,4:4.68):5.87,42:3.67):5.32):5.61,((5:5.36,38:4.69):5.02,62:3.46):4.5):5.55):5.5)
# target: plot the phylogenetic tree of jzg population
# 1. plot the time tree and the branch length is the time
# 2. give the branch a color according to the mutation tree. the colors are based on the mutation rate of each branch
# 3. add a circle to the branch and the size of the circle is the non-synonymous mutation number of each branch

from matplotlib import pyplot as plt
import numpy as np

def get_samples(node):
    samples = []
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
        

time_tree ='((((((((((72:521.12,34:521.12):172.29,55:693.41):70.79,71:764.2):111.17,(((3:458.76,36:458.76):174.55,32:633.31):106.08,49:739.38):135.99):37.4,29:912.77):87.38,(((((47:515.21,30:515.21):136.87,51:652.08):95.29,57:747.37):91.42,(45:525.69,59:525.69):313.1):48.32,28:887.11):113.03):13.17,21:1013.31):26.93,22:1040.24):154.39,((((((69:574.84,31:574.84):459.1,((((53:594.47,18:594.47):331.69,(((27:552.52,6:552.52):168.05,66:720.57):144.28,(41:628.39,39:628.39):236.46):61.31):32.83,(68:584.6,63:584.6):374.39):33.93,(54:598.59,7:598.59):394.33):41.01):14.27,61:1048.2):102.59,(((((((8:524.89,43:524.89):150.61,9:675.5):144.21,(1:539.65,10:539.65):280.07):44.24,2:863.95):46.38,64:910.33):62.57,((67:559.99,50:559.99):174.67,40:734.66):238.24):140.86,(((65:592.68,13:592.68):168.93,16:761.61):317.78,((((((20:512.53,35:512.53):306.92,((58:520.8,37:520.8):165.9,46:686.7):132.74):123.62,((((25:523.09,70:523.09):155.09,17:678.18):75.83,52:754.01):55.11,19:809.12):133.94):50.56,(((((15:476.42,14:476.42):153.85,24:630.26):93.43,60:723.69):62.71,48:786.4):48.26,44:834.66):158.96):9.99,((12:567.32,26:567.32):147.42,56:714.74):288.87):12.37,23:1015.98):63.41):34.37):37.04):5.72,(((33:559.48,11:559.48):175.55,4:735.03):106.84,42:841.87):314.65):6.57,((5:582.23,38:582.23):165.17,62:747.41):415.68):31.55):410.64)'
mutation_tree ='((((((((((72:5.18,34:8.94):7.06,55:5.02):6.38,71:4.9):6.01,(((3:6.2,36:7.8):7.0,32:5.04):6.34,49:5.24):6.07):6.04,29:3.91):5.8,(((((47:6.01,30:6.85):6.43,51:4.98):5.95,57:5.14):5.75,(45:5.09,59:5.69):5.39):5.63,28:3.63):5.34):5.6,21:3.48):5.48,22:3.1):5.34,((((((69:4.92,31:4.68):4.8,((((53:5.02,18:5.34):5.18,(((27:4.48,6:6.77):5.62,66:4.95):5.4,(41:4.09,39:8.66):6.37):5.79):5.62,(68:4.2,63:6.61):5.4):5.57,(54:4.3,7:7.38):5.84):5.62):5.49,61:4.7):5.44,(((((((8:7.14,43:7.24):7.19,9:5.85):6.74,(1:5.71,10:5.31):5.51):6.25,2:5.12):6.06,64:4.88):5.89,((67:5.05,50:5.71):5.38,40:4.27):5.01):5.63,(((65:4.73,13:4.99):4.86,16:5.17):4.96,((((((20:5.77,35:5.01):5.39,((58:5.12,37:5.04):5.08,46:6.08):5.42):5.41,((((25:7.82,70:5.97):6.9,17:6.63):6.81,52:5.67):6.52,19:5.22):6.26):5.83,(((((15:5.93,14:8.33):7.13,24:6.31):6.86,60:6.18):6.69,48:5.24):6.4,44:5.58):6.26):6.0,((12:4.35,26:8.21):6.28,56:4.76):5.78):5.96,23:4.5):5.89):5.77):5.72):5.64,(((33:6.64,11:6.28):6.46,4:4.68):5.87,42:3.67):5.32):5.61,((5:5.36,38:4.69):5.02,62:3.46):4.5):5.55):5.5)'
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
group_color = ['#FFA07A', '#6495ED', '#32CD32', '#FFC0CB']
non_synonymous ={'1194_1605': 2.0, '1163_1194': 1.0, '863_910': 1.0, '819_863': 1.0, '539_819': 4.0, '0_539_10': 4, '841_1156': 5.0, '559_735': 4.0, '0_559_11': 1, '1079_1113': 2.0, '1015_1079': 2.0, '714_1003': 1.0, '567_714': 2.0, '0_567_12': 5, '761_1079': 2.0, '0_592_13': 6, '834_993': 1.0, '630_723': 1.0, '476_630': 1.0, '0_476_14': 3, '0_476_15': 10, '0_761_16': 9, '0_678_17': 7, '594_926': 2.0, '0_594_18': 6, '0_809_19': 9, '0_539_1': 5, '819_943': 1.0, '0_512_20': 6, '1040_1194': 4.0, '1013_1040': 1.0, '0_1013_21': 7, '0_1040_22': 2, '0_1015_23': 5, '0_630_24': 3, '523_678': 3.0, '0_523_25': 3, '0_567_26': 7, '720_864': 1.0, '552_720': 1.0, '0_552_27': 4, '0_887_28': 8, '912_1000': 1.0, '0_912_29': 6, '652_747': 1.0, '515_652': 3.0, '0_515_30': 4, '0_863_2': 5, '574_1033': 3.0, '0_574_31': 4, '633_739': 3.0, '0_633_32': 6, '0_559_33': 5, '764_875': 1.0, '693_764': 1.0, '521_693': 1.0, '0_521_34': 7, '0_512_35': 5, '458_633': 2.0, '0_458_36': 3, '686_819': 2.0, '520_686': 2.0, '0_520_37': 2, '747_1163': 2.0, '582_747': 1.0, '0_582_38': 4, '628_864': 7.0, '0_628_39': 6, '734_972': 2.0, '0_734_40': 7, '0_458_3': 5, '0_628_41': 6, '0_841_42': 3, '524_675': 1.0, '0_524_43': 1, '0_834_44': 13, '525_838': 5.0, '0_525_45': 5, '0_686_46': 4, '0_515_47': 7, '0_786_48': 10, '0_739_49': 4, '0_735_4': 3, '559_734': 1.0, '0_559_50': 4, '0_652_51': 2, '0_754_52': 8, '0_594_53': 3, '598_992': 3.0, '0_598_54': 6, '0_693_55': 5, '0_714_56': 4, '0_747_57': 5, '0_520_58': 5, '0_525_59': 4, '0_582_5': 10, '0_723_60': 4, '0_1048_61': 6, '0_747_62': 7, '0_552_6': 10, '584_958': 2.0, '0_584_63': 7, '0_910_64': 13, '0_592_65': 4, '0_720_66': 9, '0_559_67': 4, '0_584_68': 4, '0_574_69': 6, '0_523_70': 9, '0_598_7': 9, '0_764_71': 4, '0_521_72': 9, '0_524_8': 8, '0_675_9': 8}

fig, ax = plt.subplots(figsize=(6, 10))
ax.set_xlim(-2, 17)
ax.set_ylim(0, 240)
#ax.set_xticks([])
#ax.set_yticks([])
ax.set_title('Tree', fontsize=20)
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
                start = 0
                # add a x_value[-1]+3 to x_value
                x_value.append(x_value[-1]+3)
                #print(x_value)
                break
            elif time_tree[j] == ')':
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
                #print(x2_value)
                start = get_after_node_time(tree)
                break
        color_base = mutation_tree.split(':')[num_count]
        for k in range(len(color_base)):
            if color_base[k] == ',' or color_base[k] == ')':
                color_base = color_base[:k]
                break
        # here color we give a double color
        if ((float(color_base)-5.5)/10) >= 0:
            color = plt.cm.Reds((float(color_base)-5.5)/5)
        else:
            color = plt.cm.Blues(abs((float(color_base)-5.5)/5))
        if start == 0:
            color = 'black'
            for k in range(1,5):
                if int(sample) in group[k]:
                    color = group_color[k-1]
                    break
            if color == 'black':
                box = plt.Rectangle((0,x_value[-1]), value/100, 0.1, facecolor='black', edgecolor='black', lw=0.1)
                ax.add_artist(box)
            else:
                box = plt.Rectangle((1,x_value[-1]), value/100-1, 0.1, facecolor='black', edgecolor='black', lw=0.1)
                ax.add_artist(box)
                box = plt.Rectangle((-0.5,x_value[-1]-1.5), 1.5, 3, facecolor=color, edgecolor='none', lw=0.1)
                ax.add_artist(box)
                box = plt.Rectangle((0,x_value[-1]), value/100, 0.1, facecolor='black', edgecolor='black', lw=0.1)
                ax.add_artist(box)
                # here we need add a mutation rate box above the line
                box = plt.Rectangle((1.2,x_value[-1]+0.2), 1, 1.5, facecolor='none', edgecolor='black', lw=0.2 )
                ax.add_artist(box)
                color_base_list =[]
                for i in range(int(float(color_base)*10)):
                    color_base_list.append(plt.cm.Reds(float(i)/100)) 
                for i in range(len(color_base_list)):
                    box = plt.Rectangle((1.2+i/100,x_value[-1]+0.2), 0.010, 1.45, facecolor=color_base_list[i], edgecolor='none', lw=0.1)
                    ax.add_artist(box)
                box = plt.Rectangle((1.75,x_value[-1]+0.2), 0.01, 1.5, facecolor='white', edgecolor='black', lw=0.1 )
                ax.add_artist(box)
                box = plt.Rectangle((1.2,x_value[-1]+0.2), 1, 1.5, facecolor='none', edgecolor='black', lw=0.2 )
                ax.add_artist(box)
            #print(sample)
            #print(int(value))
            for key in non_synonymous:
                #print(key)
                if int(key.split('_')[0]) == 0 and int(key.split('_')[1]) == int(value) and str(key.split('_')[2]) == str(sample):
                    for i in range(int(non_synonymous[key])):
                        if i%2 == 0:
                            box = plt.Rectangle((2.5+i//2* 0.25,x_value[-1]+0.2), 0.2, 1, facecolor='#B5E4D0', edgecolor='grey', lw=0.1)
                            ax.add_artist(box)
                        else:
                            box = plt.Rectangle((2.5+i//2* 0.25,x_value[-1]+1.4), 0.2, 1, facecolor='#F3AFA7', edgecolor='grey', lw=0.1)
                            ax.add_artist(box)
            # add sample name background box
            box = plt.Rectangle((-0.5,x_value[-1]-1.25), 0.5, 2.5, facecolor='white', alpha =0.5)
            ax.add_artist(box)
            # add sample name
            text = plt.text(-0.3, x_value[-1]-0.5, str(sample), fontsize=9,fontname= 'Arial')
            ax.add_artist(text)
            # add sample node
            box = plt.Rectangle((0,x_value[-1]-0.3),0.05,0.6,facecolor='black',)
            ax.add_artist(box)
        elif start != 0:
            box2 = plt.Rectangle(((start)/100,x2_value[-1]), value/100, 0.2, facecolor='black', edgecolor='black', lw=0.1)
            ax.add_artist(box2)
            box = plt.Rectangle(((start)/100,x3_start), 0.01, x3_width, facecolor='black', edgecolor='black', lw=0.1)
            ax.add_artist(box)
            # add a node symbol
            box= plt.Rectangle(((start-2)/100,x2_value[-1]-0.2), 0.05, 0.6, facecolor='black', edgecolor='black', lw=0.1)
            ax.add_artist(box)
            for key in non_synonymous:
                if int(key.split('_')[0]) == int(start) and int(key.split('_')[1]) == int(value+start):
                    for i in range(int(non_synonymous[key])):
                        if i%2 == 0:
                            box = plt.Rectangle(((start+value/2)/100-0.1+i//2* 0.25,x2_value[-1]-1), 0.2, 1, facecolor='#B5E4D0', edgecolor='grey', lw=0.1)
                            ax.add_artist(box)
                        else:
                            box = plt.Rectangle(((start+value/2)/100-0.1+i//2* 0.25,x2_value[-1]-2.2), 0.2, 1, facecolor='#F3AFA7', edgecolor='grey', lw=0.1)
                            ax.add_artist(box)
            
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
                start = 0
                # add a x_value[-1]+3 to x_value
                x_value.append(x_value[-1]+3)
                #print(x_value)
                break
            elif time_tree[j] == ')':
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
                #print(x2_value)
                start = get_after_node_time(tree)
                break
        color_base = mutation_tree.split(':')[num_count]
        for k in range(len(color_base)):
            if color_base[k] == ',' or color_base[k] == ')':
                color_base = color_base[:k]
                break
        # here color we give a double color
        if ((float(color_base)-5.5)/10) >= 0:
            color = plt.cm.Reds((float(color_base)-5.5)/5)
        else:
            color = plt.cm.Blues(abs((float(color_base)-5.5)/5))
        if start != 0:
            #print(int(start))
            color_base_list =[]
            box = plt.Rectangle(((start+value/2)/100-0.1,x2_value[-1]+0.3), 0.2, 4, facecolor='white', edgecolor='black', lw=0.2 )
            ax.add_artist(box)
            for i in range(int(float(color_base)*10)):
                color_base_list.append(plt.cm.Reds(i/100))
            for i in range(len(color_base_list)):
                #box = plt.Rectangle(((start+value/2)/100-0.1,x2_value[-1]+0.3+i/100*4), 0.2, 0.04, facecolor=color_base_list[i], edgecolor='none', lw=0.1)
                box = plt.Rectangle(((start+value/2)/100-0.1,x2_value[-1]+0.3+i/100*4), 0.2, 0.04, facecolor=color_base_list[i], edgecolor='none', lw=0.1)
                ax.add_artist(box)
            box = plt.Rectangle(((start+value/2)/100-0.1,x2_value[-1]+0.3+2), 0.2, 0.1, facecolor='white', edgecolor='black', lw=0.2 )
            ax.add_artist(box)
            box = plt.Rectangle(((start+value/2)/100-0.1,x2_value[-1]+0.3), 0.2, 4, facecolor='none', edgecolor='black', lw=0.3 )
            ax.add_artist(box) 
# add time scale
box = plt.Rectangle((0,4.4),16.5,0.4,facecolor='black')
ax.add_artist(box)
box = plt.Rectangle((0,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(0, 0.5, str(0), fontsize=14,fontname= 'Arial',ha ='center')
ax.add_artist(text)

box = plt.Rectangle((5,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(5, 0.5, str(500), fontsize=14,fontname= 'Arial',ha ='center')
ax.add_artist(text)

box = plt.Rectangle((10,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(10, 0.5, str(1000), fontsize=14,fontname= 'Arial',ha ='center')
ax.add_artist(text)
box = plt.Rectangle((15,3.4),0.05,1,facecolor='black')
ax.add_artist(box)
text=plt.text(15, 0.5, str(1500), fontsize=14,fontname= 'Arial',ha ='center')
ax.add_artist(text)

#box = plt.Rectangle((13.5,20),3.3,80,facecolor='none',edgecolor='black',lw =0.2)
#ax.add_artist(box)
text = plt.text(15.7, 95, str('Legend'), fontsize=14,fontname= 'Arial',ha ='right')
ax.add_artist(text)

ax.add_artist(box)

for i in range(100):
    box =plt.Rectangle((15.8,72+i/10*2),0.2,0.2,facecolor=plt.cm.Reds(i/100))
    ax.add_artist(box)
box = plt.Rectangle((15.8,83),0.2,0.2,facecolor='white',edgecolor='grey',lw =0.2)
ax.add_artist(box)
text = plt.text(16.05, 82.5, str('5.52'), fontsize=8,fontname= 'Arial',ha ='left')
ax.add_artist(text)
box = plt.Rectangle((15.8,72),0.25,0.1,facecolor='black',edgecolor='black',lw =0.2)
ax.add_artist(box)
text = plt.text(16.06,71.6, str('0'), fontsize=8,fontname= 'Arial',ha ='left')
ax.add_artist(text)
box = plt.Rectangle((15.8,91.9),0.25,0.1,facecolor='black',edgecolor='black',lw =0.2)
ax.add_artist(box)
text = plt.text(16.05, 91.6, str('10'), fontsize=8,fontname= 'Arial',ha ='left')
ax.add_artist(text)
box = plt.Rectangle((16,72),0.01,20,facecolor='black',edgecolor='black',lw =0.1)
ax.add_artist(box)

# off the axis
ax.axis('off')
plt.show()
