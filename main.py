#!/usr/bin/env python

"""
Author: Shohei Kojima @ RIKEN
Copyright (c) 2023 RIKEN
All Rights Reserved
See file LICENSE for details.
"""

import os,sys,argparse
sys.setrecursionlimit(10000)

# version
version='2023/1/25'

# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-i', metavar='str', type=str, help='Input phylogenetic tree (must be newick format).', required=True)
parser.add_argument('-o', metavar='str', type=str, help='Output file name.', required=True)
parser.add_argument('-t', metavar='float', type=float, help='Threshold for distance to make clusters.', required=True)
parser.add_argument('-m', metavar='int', type=int, help='Minimum count of members to make clusters. Default = 1.', default = 1)
parser.add_argument('-v', '--version', action='version', version=version)
args=parser.parse_args()
args.version=version


class Clade:
    def __init__(self):
        self.all_set=set()
        self.components=[]
        self.terminal=False
    
    def add_clade(self, clade, length):
        self.all_set |= clade.all_set
        self.components.append((clade, length))
    
    def add_terminal(self, leaf_name, length):
        self.terminal = True
        self.all_set = {leaf_name}
        self.components.append((leaf_name, length))
    
    def find_path_to_element(self, target):
        if not target in self.all_set:
            return [],[]
        if self.terminal:
            return [],[]
        for clade,length in self.components:
            if target in clade.all_set:
                path,lengths=clade.find_path_to_element(target)
                path.append(self)
                lengths.append(length)
                return path,lengths


def nwk_to_clades(ls, s, e):
    # find nested clade
    nested_found=False
    start_index=[]
    end_index=[]
    depth = 0
    for i in range(s + 1, e):
        if ls[i] == '(':
            depth += 1
            if depth == 1:
                start_index.append(i)
                nested_found = True
        elif ls[i] == ')':
            depth -= 1
            if depth == 0:
                end_index.append(i)
    if nested_found:
        clades=Clade()
        # add nested
        for si,ei in zip(start_index, end_index):
            clade,length=nwk_to_clades(ls, si, ei)
            clades.add_clade(clade, length)
        # add non nested
        depth = 0
        for i in range(s + 1, e):
            if ls[i] == '(':
                depth += 1
            elif ls[i] == ')':
                depth -= 1
            else:
                if depth > 0:
                    continue
                if ls[i-1] == ')':
                    continue
                s=ls[i]
                name,length=s.split(':')
                length = float(length)
                term_clade=Clade()
                term_clade.add_terminal(name, 0)
                clade=Clade()
                clade.add_clade(term_clade, 0)
                clades.add_clade(clade, length)
        length = float(ls[e+1].split(':')[1])
        return clades, length
    else:
        clades=Clade()
        for i in range(s + 1, e):
            if ls[i-1] == ')':
                continue
            s=ls[i]
            name,length=s.split(':')
            length = float(length)
            clade=Clade()
            clade.add_terminal(name, 0)
            clades.add_clade(clade, length)
        length = float(ls[e+1].split(':')[1])
        return clades, length


def nwk_split(line):
    ls=[]
    tmp=[]
    inside_name=0
    for c in line:
        if c == '"':
            inside_name = (inside_name + 1) % 2
            tmp.append(c)
        else:
            if inside_name == 0:
                if c == '(' or c == ')':
                    if len(tmp) >= 1:
                        ls.append(''.join(tmp))
                        tmp=[]
                    ls.append(c)
                elif c == ',':
                    ls.append(''.join(tmp))
                    tmp=[]
                elif c == ';':
                    if len(tmp) >= 1:
                        ls.append(''.join(tmp))
                else:
                    tmp.append(c)
            else:
                tmp.append(c)
    if ls[-1] == ')':
        # unrooted tree
        ls.append('_:0')
    # list up names
    names=[]
    for i in range(1, len(ls)):
        prev = ls[i - 1]
        current = ls[i]
        if not current == '(' and not current == ')':
            if prev == ')':
                continue
            name,_ = current.split(':')
            names.append(name)
    return ls,names


def nwk_parser(f):
    clusters=[]
    with open(f) as infile:
        line=next(infile).strip()
    ls,leaf_names=nwk_split(line)
    # find isolated clusters
    clades,length = nwk_to_clades(ls, 0, len(ls)-2)
    return clades,length,leaf_names


def calc_distance(clades, target_leaf, seed_path, seed_lengths):
    target_path,target_lengths = clades.find_path_to_element(target_leaf)
    min_length = min(len(seed_lengths), len(target_lengths))
    seed_path=seed_path[::-1]
    target_path=target_path[::-1]
    seed_lengths=seed_lengths[::-1]
    target_lengths=target_lengths[::-1]
    shared = 0
    for i in range(min_length):
        sp = seed_path[i]
        tp = target_path[i]
        if sp == tp:
            shared += 1
        else:
            break
    shared -= 1
    distance = sum(seed_lengths[shared:]) + sum(target_lengths[shared:])
    return distance


def find_close_leaves(clades, target_name, classified, threshold):
    res=[]
    classified.add(target_name)
    path,lengths = clades.find_path_to_element(target_name)
    total_length=0
    for clade,length in zip(path, lengths):
        total_length += length
        for name in clade.all_set:
            if name in classified:
                continue
            distance = calc_distance(clades, name, path, lengths)
            if distance <= threshold:
                res.append(name)
                classified.add(name)
                tmp=find_close_leaves(clades, name, classified, threshold)
                res.extend(tmp)
        if total_length > threshold:
            break
    return res


def sort_leaf_names(leaf_names, old):
    # sort
    old=set(old)
    new=[]
    for name in leaf_names:
        if name in old:
            new.append(name)
    return new


def main(args):
    threshold = args.t
    member_threshold = args.m
    clustered = 0
    not_clustered = 0
    
    # parse nwk
    clades,length,leaf_names=nwk_parser(args.i)
    classified=set()
    # find clusters
    clade_num=0
    out=['clade\tleaf_name\n']
    for name in leaf_names:
        if name in classified:
            continue
        classified.add(name)
        close_leaves=[name]
        tmp=find_close_leaves(clades, name, classified, threshold)
        close_leaves.extend(tmp)
        close_leaves=sort_leaf_names(leaf_names, close_leaves)
        n_member = len(close_leaves)
        if n_member >= member_threshold:
            clade_name = 'cluster%d' % clade_num
            clade_num += 1
            clustered += n_member
        else:
            clade_name = 'Not_clustered'
            not_clustered += n_member
        for name in close_leaves:
            out.append('%s\t%s\n' % (clade_name, name))
        total_processed = clustered + not_clustered
        print('\rLeaves processed: %d' % total_processed, end = '')
    
    with open(args.o, 'w') as outfile:
        outfile.write(''.join(out))
    
    print('\n%d leaves were clustered into %s clusters.' % (clustered, clade_num))
    print('%d leaves were not clustered.' % not_clustered)

if __name__ == '__main__':
    main(args)
