# -*- coding: utf-8 -*-
"""Disjoint paths using  Bhandari's algorithm based on the Bellman-Ford algorithm.
"""

#    Copyright (C) 2015-2016 by
#    Niels L. M. van Adrichem <n.l.m.vanadrichem@tudelft.nl>
#    All rights reserved.
#    BSD license.

__author__ = """\n""".join(['Niels L. M. van Adrichem <n.l.m.vanadrichem@tudelft.nl>'])

import networkx as nx

def bhandari(G, source, target, weight=None, k=2, node_disjoint=False, force_maximally_disjoint=False, node_over_edge_disjointness=False, edge_disjointness_penalty=None, node_disjointness_penalty=None):

    if k < 2:
        raise nx.NetworkXUnfeasible("You need at least k>=2 paths to be disjoint, k=%d"%(k))
    
    if source == target:
        raise nx.NetworkXUnfeasible("There is no such thing as a disjoint path to oneself, as oneself has to excluded from the disjoint path to be able to exist")

    if G.has_negative_edges():
        raise nx.NetworkXUnfeasible("Bhandari's algorithm is not correct for negative edges as it may insert negative cycles.")

    if force_maximally_disjoint == True or node_over_edge_disjointness == True or edge_disjointness_penalty != None or node_disjointness_penalty != None:
        raise NotImplementedError(
            "To Be Done: Maximally-disjoint paths not yet implemented")
    
    if G.is_multigraph():
        raise nx.NetworkXUnfeasible(
            "Apply link-splitting before calling Bhandari's to allow multigraphs. Note that multigraphs are not useful when searching for non-maximal node disjointness, a graph containing the minimum-weight edges suffices.")
    
    #Find an unused parameter 'weight' so we don't mixup our internal weights with
    #Possible existing parameters named "weight" that should not be considered.
    if weight == None:
        weight = "_weight"
        while any( weight in d for u, v, d in G.edges(data = True) ):
            weight = "_"+weight
    
    
    if G.is_multigraph():
        raise NotImplementedError(
            "Implement link-splitting to allow multigraphs, impossible in the case of node-disjointness")
    
    def split_node_criteria(n):
        #Don't split source and destination nodes
        if n == source or n == target:
            return False
        #Already splitted by us
        if n not in G:
            return False
        
        #An absolute degree <= 3 does not need to be splitted, neither do nodes with just 1 in- or out-going edge.
        if G.is_directed():
            if G.out_degree(n) == 1:
                return False
            if G.in_degree(n) == 1:
                return False
            #Distinct degree of incoming and outgoing node
            if ( len(G.succ[n]) + sum(1 for _node in G.pred[n] if _node not in G.succ[n]) ) <= 3:
                return False
        else:
            if G.degree(u) <= 3:
                return False
                
        return True
           
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
                
                if node_disjoint == True and split_node_criteria(u):
                    (u_in, u_out) = G_copy.split_node(u, weight=weight)
                    
                    #Swap direction and weight
                    G_copy.remove_edge(u_in, u_out)
                    G_copy.add_edge(u_out, u_in, {weight : -0} )
                    
                    G_copy.add_edge(v, u_out, {weight : -weight_val})

                    #Add to internal path
                    edges[(u_out, v)] = (u_in , next)
                    edges[(u_in, u_out)] = (prev, v)
                    
                    #Set values for next iteration
                    u = u_in
                    v = u_out                    
                    
                else: #Being in edges implies that u and v are already node-splitted if necessary
                    #Swap direction and weight
                    G_copy.add_edge(v, u, {weight : -weight_val})
                    
                    #Add to internal path
                    edges[(u, v)] = (prev, next)
                    
                #Set value for next iteration
                next = v
                
            else: #(v,u) in edges, thus also already splitted where applicable
                #Restore original edges in graph, if they exists:
                weight_val = -G_copy[u][v].get(weight, 1)                
                G_copy.remove_edge(u,v)

                #Add original arc weight
                G_copy.add_edge(v, u, {weight : weight_val})
                
                if node_disjoint == True:
                    
                    if u not in G_orig:
                        _u = G_copy.node[u]["split_node"]
                    else:
                        _u = u
                    
                    if v not in G_orig:
                        _v = G_copy.node[v]["split_node"]
                    else:
                        _v = v
                    #If original arc existed, restore it (possibly a self-loop, but that seems to be correct)
                    if G_orig.has_edge(_u, _v): 
                        G_copy.add_edge(u, v, {weight : G_orig[_u][_v].get(weight, 1)})
                    
                else:
                    #If original arc existed, restore it
                    if G_orig.has_edge(u, v):
                        G_copy.add_edge(u, v, {weight : G_orig[u][v].get(weight, 1)})
                    
                

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
    #Compute paths and true lengths
    paths = []
    dists = []
    
    for (u, v) in path_init:
        
        path = [u]
        dist = 0
        
        if node_disjoint == True:
            u_split = u
            
        while v != None:

            if node_disjoint == True:
                if v not in G:
                    v_split = G_copy.node[v].get("split_node")
                    if u_split != v_split:
                        path.append(v_split)
                        dist += G[u_split][v_split].get(weight, 1)
                        u_split = v_split
                else:
                    path.append(v)
                    dist += G[u_split][v].get(weight, 1)
                    u_split = v
                    
            else:
                path.append(v)
                dist += G[u][v].get(weight, 1)
 
            (_, next) = edges[(u, v)]
            u = v
            v = next

        dists.append(dist)
        paths.append(path)
        
    dists,paths = zip(*sorted(zip(dists,paths)))        
        
    return dists,paths