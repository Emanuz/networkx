# -*- coding: utf-8 -*-
"""Disjoint paths using path-removal and Dijkstra's algorithm.
"""

#    Copyright (C) 2015-2016 by
#    Niels L. M. van Adrichem <n.l.m.vanadrichem@tudelft.nl>
#    All rights reserved.
#    BSD license.

__author__ = """\n""".join(['Niels L. M. van Adrichem <n.l.m.vanadrichem@tudelft.nl>'])
__all__ = ['extended_disjoint']

import networkx as nx
from collections import defaultdict

def extended_disjoint(G, weight=None, node_disjoint=False, edge_then_node_disjoint=False):
    if node_disjoint == True and edge_then_node_disjoint == True:
        raise nx.NetworkXUnfeasible("Edge-then-node disjointness overrides node disjointness")
        
        
    if G.is_multigraph():
        raise nx.NetworkXUnfeasible(
            "Apply link-splitting before calling this algorithm to allow multigraphs. Note that multigraphs are not useful when searching for guaranteed node disjointness, a graph containing the minimum-weight edges suffices.")
        
    if weight == None:
        weight = "_weight"
        while any( weight in d for u, v, d in G.edges(data = True) ):
            weight = "_"+weight
    
    succ,dist = nx.floyd_warshall_successor_and_distance(G, weight=weight)
    
    #print succ
    
    succ = defaultdict(dict, succ)
    
    #Create forwarding-matrix
    

    #For memory overhead, we should try to make a more shallow copy that only stores the one edge that gets removed from the original
    G_copy = G.copy(with_data=False)
    
    if node_disjoint == True or edge_then_node_disjoint == True:        
        for u in G:
            for v in G.neighbors(u):
            
                #print "Removing node %s"%(u,)
                G_copy.remove_node(v)
                
                _pred,_dist = nx.bellman_ford_predecessor_and_distance(G_copy, u, weight=weight)
                dst_affected = [n for n in succ[u] if succ[u][n] == v]
                #print dst_affected
                
                for n in dst_affected:
                    if n in _dist: #Check if node is reachable at all through another path
                        next = n
                        while next != u:
                            prev = _pred[next][0]
                            if succ[prev].get(n) != next:
                                succ[(prev,v)][n] = next
                            next = prev
                    else:
                        succ[(u,v)][n] = None
                
                #Restore the copy
                G_copy.add_node(v, G.node[v])
                G_copy.add_edges_from(G.edges(nbunch=G[v], data=True))
            
    if node_disjoint == False or edge_then_node_disjoint == True:
        for u in G:
            for v in G.neighbors(u):
                #print "Removing edge %s-%s"%(u,v)
                G_copy.remove_edge(u, v)
    
                _pred,_dist = nx.bellman_ford_predecessor_and_distance(G_copy, u, weight=weight)
                
                dst_affected = [n for n in G if succ[u][n] == v]
                #print dst_affected
                
                for n in dst_affected:
                    if n in _dist:
                        next = n                
                        while next != u:
                            prev = _pred[next][0]
                            if (prev not in succ or succ[prev].get(n) != next) and (edge_then_node_disjoint == False or (prev,v) not in succ or succ[(prev,v)].get(n) != next):
                                succ[(prev,(u,v))][n] = next
                            next = prev
                    elif edge_then_node_disjoint == False or (u,v) not in succ: # succ[(u,v)].get(n) != None), if there is no edge-failure-disjoint detour, neither is there a node-failure-disjoint one, no need to check
                        succ[(u,(u,v))][n] = None
    
                #Restore the copy
                G_copy.add_edge(u, v, G[u][v])
            
    return dict(succ)