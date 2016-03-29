# -*- coding: utf-8 -*-
"""Floyd-Warshall algorithm for shortest paths.
"""
#    Copyright (C) 2004-2016 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
import networkx as nx
from collections import defaultdict
__author__ = """Aric Hagberg <aric.hagberg@gmail.com>"""
__all__ = ['floyd_warshall',
           'floyd_warshall_predecessor_and_distance',
           'floyd_warshall_successor_and_distance',
           'floyd_warshall_predecessor_successor_and_distance',
           'floyd_warshall_numpy']

def floyd_warshall_numpy(G, nodelist=None, weight='weight'):
    """Find all-pairs shortest path lengths using Floyd's algorithm.

    Parameters
    ----------
    G : NetworkX graph

    nodelist : list, optional
       The rows and columns are ordered by the nodes in nodelist.
       If nodelist is None then the ordering is produced by G.nodes().

    weight: string, optional (default= 'weight')
       Edge data key corresponding to the edge weight.

    Returns
    -------
    distance : NumPy matrix
        A matrix of shortest path distances between nodes.
        If there is no path between to nodes the corresponding matrix entry
        will be Inf.

    Notes
    ------
    Floyd's algorithm is appropriate for finding shortest paths in
    dense graphs or graphs with negative weights when Dijkstra's
    algorithm fails.  This algorithm can still fail if there are
    negative cycles.  It has running time O(n^3) with running space of O(n^2).
    """
    try:
        import numpy as np
    except ImportError:
        raise ImportError(\
          "to_numpy_matrix() requires numpy: http://scipy.org/ ")

    # To handle cases when an edge has weight=0, we must make sure that
    # nonedges are not given the value 0 as well.
    A = nx.to_numpy_matrix(G, nodelist=nodelist, multigraph_weight=min,
                              weight=weight, nonedge=np.inf)
    n,m = A.shape
    I = np.identity(n)
    A[I==1] = 0 # diagonal elements should be zero
    for i in range(n):
        A = np.minimum(A, A[i,:] + A[:,i])
    return A

def floyd_warshall_predecessor_and_distance(G, weight='weight'):
    """Find all-pairs shortest path lengths using Floyd's algorithm.

    Parameters
    ----------
    G : NetworkX graph

    weight: string, optional (default= 'weight')
       Edge data key corresponding to the edge weight.

    Returns
    -------
    predecessor,distance : dictionaries
       Dictionaries, keyed by source and target, of predecessors and distances
       in the shortest path.

    Notes
    ------
    Floyd's algorithm is appropriate for finding shortest paths
    in dense graphs or graphs with negative weights when Dijkstra's algorithm
    fails.  This algorithm can still fail if there are negative cycles.
    It has running time O(n^3) with running space of O(n^2).

    See Also
    --------
    floyd_warshall
    floyd_warshall_numpy
    all_pairs_shortest_path
    all_pairs_shortest_path_length
    """
    pred = defaultdict(dict)
    dist = _floyd_warshall(G, weight=weight, pred=pred)
    return dict(pred), dict(dist)

def floyd_warshall_successor_and_distance(G, weight='weight'):
    """Find all-pairs shortest path lengths using Floyd's algorithm.

    Parameters
    ----------
    G : NetworkX graph

    weight: string, optional (default= 'weight')
       Edge data key corresponding to the edge weight.

    Returns
    -------
    successor,distance : dictionaries
       Dictionaries, keyed by source and target, of successors and distances
       in the shortest path.

    Notes
    ------
    Floyd's algorithm is appropriate for finding shortest paths
    in dense graphs or graphs with negative weights when Dijkstra's algorithm
    fails.  This algorithm can still fail if there are negative cycles.
    It has running time O(n^3) with running space of O(n^2).

    See Also
    --------
    floyd_warshall
    floyd_warshall_numpy
    all_pairs_shortest_path
    all_pairs_shortest_path_length
    """
    succ = defaultdict(dict)
    dist = _floyd_warshall(G, weight=weight, succ=succ)
    return dict(succ), dict(dist)

def floyd_warshall_predecessor_successor_and_distance(G, weight='weight'):
    """Find all-pairs shortest path lengths using Floyd's algorithm.

    Parameters
    ----------
    G : NetworkX graph

    weight: string, optional (default= 'weight')
       Edge data key corresponding to the edge weight.

    Returns
    -------
    predecessor,successor,distance : dictionaries
       Dictionaries, keyed by source and target, of predecessors, successors
       and distances in the shortest path.

    Notes
    ------
    Floyd's algorithm is appropriate for finding shortest paths
    in dense graphs or graphs with negative weights when Dijkstra's algorithm
    fails.  This algorithm can still fail if there are negative cycles.
    It has running time O(n^3) with running space of O(n^2).

    See Also
    --------
    floyd_warshall
    floyd_warshall_numpy
    all_pairs_shortest_path
    all_pairs_shortest_path_length
    """
    succ = defaultdict(dict)
    pred = defaultdict(dict)
    dist = _floyd_warshall(G, weight=weight, pred=pred, succ=succ)
    return dict(pred), dict(succ), dict(dist)    
    
def _floyd_warshall(G, weight, pred=None, dist=None, succ=None):
    
    # dictionary-of-dictionaries representation for dist and pred
    # use some defaultdict magick here
    # for dist the default is the floating point inf value
    
    if dist == None:
        dist = defaultdict(dict)

    #if pred == None:
    #    pred = defaultdict(dict)
 
    #if succ == None:
    #    succ = defaultdict(dict)

    inf = float('inf')
    
    for u in G:
        dist[u][u] = 0
        if pred != None:
            pred[u][u] = None
        if succ != None:
            succ[u][u] = None
    
    # initialize path distance dictionary to be the adjacency matrix
    # also set the distance to self to 0 (zero diagonal)
    undirected = not G.is_directed()
    for u,v,d in G.edges(data=True):
        e_weight = d.get(weight, 1.0)
        dist[u][v] = min(e_weight, dist[u].get(v, inf))
        if pred != None:
            pred[u][v] = u
        if succ != None:
            succ[u][v] = v
        if undirected:
            dist[v][u] = min(e_weight, dist[v].get(u, inf))
            if pred != None:
                pred[v][u] = v
            if succ != None:
                succ[v][u] = u
    for w in G:
        for u in G:
            for v in G:
                if dist[u].get(v, inf) > dist[u].get(w, inf) + dist[w].get(v, inf):
                    dist[u][v] = dist[u][w] + dist[w][v]
                    if pred != None:
                        pred[u][v] = pred[w][v]
                    if succ != None:
                        succ[u][v] = succ[u][w]

    return dist


def floyd_warshall(G, weight='weight'):
    """Find all-pairs shortest path lengths using Floyd's algorithm.

    Parameters
    ----------
    G : NetworkX graph

    weight: string, optional (default= 'weight')
       Edge data key corresponding to the edge weight.


    Returns
    -------
    distance : dict
       A dictionary,  keyed by source and target, of shortest paths distances
       between nodes.

    Notes
    ------
    Floyd's algorithm is appropriate for finding shortest paths
    in dense graphs or graphs with negative weights when Dijkstra's algorithm
    fails.  This algorithm can still fail if there are negative cycles.
    It has running time O(n^3) with running space of O(n^2).

    See Also
    --------
    floyd_warshall_predecessor_and_distance
    floyd_warshall_numpy
    all_pairs_shortest_path
    all_pairs_shortest_path_length
    """
    # could make this its own function to reduce memory costs
    return _floyd_warshall(G, weight=weight)

# fixture for nose tests
def setup_module(module):
    from nose import SkipTest
    try:
        import numpy
    except:
        raise SkipTest("NumPy not available")
