#!/usr/bin/env python

# Python 3.7 or later
import ete3
import matplotlib.pyplot as plt
from matplotlib import colors


f_nwk='test.nwk'
f_clade='test.txt'

colors=['green', 'red', 'orange', 'blue', 'purple']

def layout(node):
    if node.is_leaf():
        ete3.faces.add_face_to_node(ete3.TextFace(node.name, fsize=10), node, column=0)

# read clades and determine colors
name_to_clade_color={}
current_clade=''
color_index=0
with open(f_clade) as infile:
    next(infile)
    for line in infile:
        ls=line.split()
        if ls[0] == 'Not_clustered':
            color = 'silver'
        elif not ls[0] == current_clade:
            color_index= (color_index + 1) % len(colors)
            color=colors[color_index]
        name_to_clade_color[ls[1]]=(ls[0], color)
        current_clade=ls[0]


# read nwk
t=ete3.Tree(f_nwk)
ts=ete3.TreeStyle()
#ts.mode='c'  # draw in circular
ts.show_branch_support=False
ts.layout_fn = layout
ts.show_leaf_name = False


# plot
for node in t.traverse():
    if node.is_leaf() is True:
        name=node.name
        clade,color=name_to_clade_color[name]
        node.name= '%s:%s' % (clade, node.name)
        nstyle=ete3.NodeStyle()
        nstyle['fgcolor']=color
        nstyle['size']=10
        node.set_style(nstyle)

t.render('test.pdf', w=1000, units='mm', tree_style=ts)

