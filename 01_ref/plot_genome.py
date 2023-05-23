#!/py
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: plot the genomic information by circles

import sys
import os
import pycircos
import collections
import matplotlib
import math
# as the python package is too many, use the following code to import the 2 classes
from pycircos.pycircos import Garc
from pycircos.pycircos import Gcircle 

circle = Gcircle(figsize=(10,10)) 
with open("path/chr_genral.csv") as f:
    f.readline()
    for line in f:
        line   = line.rstrip().split("\t") 
        name   = line[0]
        length = int(line[-1]) 
        arc    = Garc(arc_id=name, size=length, interspace=2, raxis_range=(800,850), labelposition=60, label_visible=True, labelsize=14)
        circle.add_garc(arc) 
circle.set_garcs(-90,90) 
for arc_id in circle.garc_dict:
    circle.tickplot(arc_id, raxis_range=(850,853), tickinterval=2500000)

#bar plot1
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("path/trf_num.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["values"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(line[-1]))
        values_all.append(float(line[-1]))
vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:  
    circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], base_value=0.0, rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)],
                   raxis_range=[650,700], facecolor="#8ECFC9", spine=True)

values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("path/TE_Copia_num.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["values"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(line[-1]))
        values_all.append(float(line[-1]))
vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:  
    circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], base_value=0.0, rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)],
                   raxis_range=[600,650], facecolor="#55A5E9", spine=True)

values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("path/TE_Gypsy_num.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["values"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(line[-1]))
        values_all.append(float(line[-1]))
vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:  
    circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], base_value=0.0, rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)],
                   raxis_range=[550,600], facecolor="#FA7F6F", spine=True)

#bar plot1
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("path /gene_density_per100k.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split("\t")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["values"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(0-float(line[-1]))
        values_all.append(0-float(line[-1]))
vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:  
    circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], base_value=0.0, rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)],
                   raxis_range=[700,800], facecolor="#2878b5", spine=True)



#
#bar plot3
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("path/repeat.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split("\t")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["values"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(1))
        values_all.append(float(1))
vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:  
    circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], base_value=0.0, rlim=[vmin, vmax],
                   raxis_range=[520,550], facecolor="#12853A", spine=True)

# link plot
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("path/links.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(" ")
        name1  = line[0]     
        start1 = int(line[1])-1
        end1   = int(line[2])
        name2  = line[3]     
        start2 = int(line[4])-1
        end2   = int(line[5])
        source = (name1, start1, end1, 520)
        destination = (name2, start2, end2, 520)
        circle.chord_plot(source, destination, facecolor=circle.garc_dict[name1].facecolor)

circle.figure.savefig("path/genome.png", dpi=600)


