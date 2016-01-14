# -*- coding: utf-8 -*-
"""Disjoint paths using path-removal and Dijkstra's algorithm.
"""

#    Copyright (C) 2015-2016 by
#    Niels L. M. van Adrichem <n.l.m.vanadrichem@tudelft.nl>
#    All rights reserved.
#    BSD license.

__author__ = """\n""".join(['Niels L. M. van Adrichem <n.l.m.vanadrichem@tudelft.nl>'])
__all__ = ['simple_disjoint_paths']

import networkx as nx

def simple_disjoint_paths(G, source, target, weight=None, k=2, node_disjoint=False):
    """ Simplest implementation ever by removing the shortest path and computing the next shortest path, no guarantees of finding anymore than the shortest path, even if other possiblities exist.


    """

    if k < 2:
        raise nx.NetworkXUnfeasible("You need at least k>=2 paths to be disjoint, k=%d"%(k))
        
    if source == target:
        raise nx.NetworkXUnfeasible("There is no such thing as a disjoint path to oneself, as oneself has to excluded from the disjoint path to be able to exist")
        
    paths = []
    G_copy = G.copy()

    for i in range(0, k):

        #plt.figure(i)
        #nx.draw(G_copy)            
        
        try:
            path = nx.shortest_path(G_copy, source, target, weight)
        except nx.NetworkXNoPath:
            raise nx.NetworkXNoPath(
                "Cannot find more than %d disjoint path(s)"%(i))
            
        
        paths.append(path)
        
        if i == k:
            break
        
        if node_disjoint:        
            for j in range(1, len(path)-1):
                G_copy.remove_node(path[j])
                
        else:
            for j in range(0, len(path)-1):
                G_copy.remove_edge(path[j], path[j+1])
                
                
    return paths

        
