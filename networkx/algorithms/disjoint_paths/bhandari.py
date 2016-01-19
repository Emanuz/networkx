# -*- coding: utf-8 -*-
"""Disjoint paths using  Bhandari's algorithm based on the Bellman-Ford algorithm.
"""

#    Copyright (C) 2015-2016 by
#    Niels L. M. van Adrichem <n.l.m.vanadrichem@tudelft.nl>
#    All rights reserved.
#    BSD license.

__author__ = """\n""".join(['Niels L. M. van Adrichem <n.l.m.vanadrichem@tudelft.nl>'])
__all__ = ['bhandari_paths']

import networkx as nx

def bhandari_paths(G, source, target, weight=None, k=2, node_disjoint=False, force_maximally_disjoint=False, node_over_edge_disjointness=True, cutoff_disjointness=False):

    if k < 2:
        raise nx.NetworkXUnfeasible("You need at least k>=2 paths to be disjoint, k=%d"%(k))
    
    if source == target:
        raise nx.NetworkXUnfeasible("There is no such thing as a disjoint path to oneself, as oneself has to excluded from the disjoint path to be able to exist")

    if G.has_negative_edges():
        raise nx.NetworkXUnfeasible("Bhandari's algorithm may not give correct results for negative edges")
    
    #if node_disjoint == True:
    #    raise NotImplementedError(
    #        "To Be Done: Node-disjoint paths are not yet implemented")

    if force_maximally_disjoint == True:
        raise NotImplementedError(
            "To Be Done: Maximally-disjoint paths not yet implemented")
            
    if cutoff_disjointness == True:
        raise NotImplementedError(
            "To Be Done: Cutoff for maximal disjointness not yet implemented")    
    
    if G.is_multigraph():
        raise NotImplementedError(
            "Implement link-splitting to allow multigraphs, impossible in the case of node-disjointness")
    
    #Find an unused parameter so we don't mixup our internal weights with
    #Possible existing parameters named "weight" that should not be considered.
    if weight == None:
        weight = "_weight"
        while any( weight in d for u, v, d in G.edges(data = True) ):
            weight = "_"+weight
    
    #We will use node-splitting to assure node_disjointness, hence the
    #"original" graph may change.
    if node_disjoint == True:
        G_orig = nx.DiGraph(G)
        for node in G_orig.nodes():
    else:
        G_orig = G

    G_copy = nx.DiGraph(G_orig)
    path_init = []
    edges = {}
   
    for i in range(0, k):
        (pred, dist) = nx.bellman_ford_predecessor_and_distance(G_copy, source, target=target, weight=weight)
        if target not in dist:
            raise nx.NetworkXNoPath(
                "Cannot find more than %d disjoint path(s)"%(i))
                
        next = None
        v = target
        u = pred[v][0]
        
        while v != source:            
            prev = pred[u][0] if u in pred else None
            
            if (v,u) not in edges:                
                #Remove edge from graph and add negatively-weighted arc in opposite direction
                weight_val = G_copy[u][v].get(weight, 1)
                G_copy.remove_edge(u, v)

                #If reverse arc exists, replace it
                if G_copy.has_edge(v, u):
                    G_copy.remove_edge(v, u)
                
                G_copy.add_edge(v, u, {weight : -weight_val})

                #Add to internal path
                edges[(u, v)] = (prev, next)
                
                #Set value for next iteration
                next = v
                
            else: #(v,u) in edges, thus also already splitted where applicable
                #Restore original edge in graph, if it exists:
                weight_val = G_copy[u][v].get(weight, 1)                
                G_copy.remove_edge(u,v)
                
                #If original arc existed, restore it
                if G_orig.has_edge(u, v):
                    G_copy.add_edge(u, v, {weight : G_orig[u][v].get(weight, 1)})
                    
                #Add original arc weight
                G_copy.add_edge(v, u, {weight : -weight_val})

                #Switch paths
                #Connect encountered path to our tail
                (sPrev, sNext) = edges[(v, u)]
                del edges[(v, u)]
                
                (sPrevPrev, _) = edges[(sPrev, v)]
                edges[(sPrev, v)] = (sPrevPrev, next)
                
                (_, nextNext) = edges[(v, next)]
                edges[(v, next)] = (sPrev, nextNext)
                
                #Connect found tail to our path in next round
                next = sNext
                
            v = u
            u = prev
                
        path_init.append( (source, next) )
        
    paths = []
    for (u, v) in path_init:
        path = [u]
        
        while v != None:
            path.append(v)
            (_, next) = edges[(u, v)]
            u = v
            v = next
            
        paths.append(path)
        
    return paths