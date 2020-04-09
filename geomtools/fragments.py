#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains anything related to analysis of distances, ordering, and combining/dividing fragments
"""
import numpy as np
#import sys
from geomtools.geom import geom

def mix_geoms(g1, g2, ratio=0.5):
    """
    Note
    ----
    Yields a midpoint between two geometries with a specific ratio. Very useful if a geometry optimisation "oscillates" between two structures
    Parameters
    ----------
    g1 : geom
        geometry 1
    g2 : geom
        geometry 2
    ratio : float
        weight of g1, default is 0.5. 
        
    Returns
    -------
    geom
        mixed geometry
    """
    from geomtools.geom import check_same_atoms,check_same_unit
    check_same_atoms(g1,g2)
    check_same_unit(g1,g2)
    return geom(g1.atoms,np.add(np.multiply(ratio,g1.coords),np.multiply(1-ratio,g2.coords)))
                
def dist_order(g):
    """
    Note
    ----
    d[order] yields the distances in increasing order
    
    Parameters
    ----------
    g : geom
        geometry to reorder
    
    Returns
    -------
    tuple
        (array of distances in the order of g.atoms, order of atom indexes from the closest to the furthest)
    """
    d=np.sqrt(np.add(np.add((g.coords[:,0])**2,(g.coords[:,1])**2),(g.coords[:,2])**2))
    order=d.argsort()
    return (d, order)

def check_geoms_order(g1, g2, thresh=0.00005):
    """
    Note
    ----
    Checks if the atoms in two geometries are in the same order.
    This is done via ordering according to the distance from the com.
    In case of symmetry a new random point near the com is used to distinguish atoms.
    
    Parameters
    ----------
    g1 : geom
        geometry 1
    g2 : geom
        geometry 2
    thresh : float
        threshold distance below which two points are considered equivalent, default is 0.00005
        
    Returns
    -------
    bool
        False if the geometry are ordered differently
        True if they are ordered in the same way
    """
    g1.check_same_atoms(g2)
    g1.check_same_units(g2)
    d1, o1 = dist_order(g1)
    d2, o2 = dist_order(g2)
    if (g1.atoms[o1]!=g2.atoms[o1]).any() or (d1[o1] - d2[o1] > thresh).any():
        return False
    else:
        ncm=np.add(np.multiply(g1.get_com_coords()[o1][0],np.multiply(3,thresh)),np.multiply(np.multiply(1.5,thresh),np.random.rand(3)))
        from transformations import center
        g1a,g2a = center(g1), center(g2)
        g1a.transl(ncm), g2a.transl(ncm) 
        da1, oa1 = dist_order(g1a)
        da2, oa2 = dist_order(g2a)
        return (da1[oa1] - da2[oa1] < thresh).all()
    
def geoms_reorder_equally(g1, g2, thresh=0.00005):
    """
    Note
    ----
    Reorders two molecules to the same ordering.
    This is done via ordering according to the distance from the com.
    In case of symmetry a new random point near the com is used to distinguish atoms.
    Parameters
    ----------
    g1 : geom
        geometry 1
    g2 : geom
        geometry 2
    thresh : float
        threshold distance below which two points are considered equivalent, default is 0.00005
        
    Returns
    -------
    tuple
        (geometry 1 reordered, geometry 2 reordered)
    """
    from geomtools.geom import check_same_atoms,check_same_unit
    check_same_atoms(g1,g2)
    check_same_unit(g1,g2)
    d1, o1 = dist_order(g1)
    d2, o2 = dist_order(g2)
    if (d1[o1] - d2[o1] > thresh).any():
        raise geom.otherError("These are not the same geometry!!!")
    else:
        ncm=np.add(np.multiply(g1.get_com_coords()[o1][0],np.multiply(3,thresh)),np.multiply(np.multiply(1.5,thresh),np.random.rand(3)))
        g1a, g2a = geom(g1.atoms,g1.get_com_coords() - ncm), geom(g2.atoms,g2.get_com_coords() - ncm)  #ctodo: use transl
        da1, oa1 = dist_order(g1a)
        da2, oa2 = dist_order(g2a)   
        return (geom(g1.atoms[oa1],g1.coords[oa1]),geom(g2.atoms[oa2],g2.coords[oa2]))  #TODO other attributes

def get_ordd_dist_mat(g):
    """
    Parameters
    ----------
    g : geom
        geometry to obtain the distance matrix of
    
    Returns
    -------
    tuple
        (distance matrix, ordering of the matrix rows and columns)
        the distance matrix is ordered so that the first line's values are in increasing order
    """
    from scipy.spatial import distance_matrix
    d, o = dist_order(g)
    DM1=distance_matrix(g.coords,g.coords)
    row_ord=DM1[o[-1]].argsort()
    P = np.identity(len(row_ord))[:,row_ord]
    return np.dot(np.dot(P.T,DM1),P), row_ord

def get_bond_matrix(g): #todo: sum of vdw radii
    """
    Parameters
    ----------
    g : geom
        geometry to obtain the bond matrix of
    
    Returns
    -------
    tuple
        (bond matrix, ordering of the matrix rows and columns)
        the bond matrix has a 0 if the atoms are not bound and a 1 if they are
    """
    import mendeleev as md
    DM, row_ord =get_ordd_dist_mat(g)
    BM=np.zeros([len(row_ord),len(row_ord)])
#    for i in range(len(row_ord)):
#        for j in range(len(row_ord)):
#            R_i, R_j = np.multiply(0.01,md.element(g.atoms[row_ord][i]).vdw_radius), np.multiply(0.01,md.element(g.atoms[row_ord][j]).vdw_radius)
#            BM[i][j]=(DM[i][j] < R_i + R_j)
    for i in range(len(row_ord)):
        R_i = np.multiply(0.01,md.element(g.atoms[row_ord][i]).vdw_radius)
        BM[i]=(DM[i] < R_i)
    BM=np.maximum(BM,BM.T)
    return (BM, row_ord)

def comb_geoms(*args, List=False, keep_charges=False, keep_identifier=True):
    """
    Parameters
    ----------
    the geometry to combine can be given as:
        separately=> args=geom1,geom2...
        list=> args=[geom1,geom2...]
    
    Return
    ------
    geom
        combined geometry
    """
    if len(args)==1:
        geoms=args[0]
    else:
        geoms=args

    g_comb  = geoms[0].copy()
    for i in geoms[1:]:
        g_comb.check_same_unit(i)
        g_comb.add_atoms(i, keep_charges=keep_charges, keep_identifier=keep_identifier)
    return g_comb

def get_fragments(g):  
    """
    Parameters
    ----------
    g : geom
        geometry to divide into fragments
    
    Returns
    -------
    list[geoms]
        list of fragment geometries
    """
    import itertools as it
    BM, o = get_bond_matrix(g)
    Natms=len(o)
    bonds = np.argwhere(BM==1)
    bonds = bonds[bonds[:,0] != bonds[:,1]]
    j=0  # counter of fragments found
    f=[]  # will be list of lists, NB, just indexes
    while len(list(it.chain.from_iterable(f))) < Natms:  # obtains the list f[0]+f[1]+...+f[n]
        f.append([])  # New empty fragment
        f[j]=[min([i for i in range(Natms) if i not in list(it.chain.from_iterable(f))])]  # puts first unassigned atom in newly formed frag 
        cnt=0  #counter of atoms checked
        while cnt<len(f[j]):
            for i in f[j][cnt:]:  # only checking new ones
                cnt=len(f[j])  #updating "checked" counter
                f[j].extend([z for z in bonds[bonds[:,0]==i][:,1] if z not in f[j]])  #adding any new bound atom
        j+=1  # counter of fragments found
    frags=[]
    for i in f:
        frags.append(geom(g.atoms[o][i],g.coords[o][i]))  # going from indexes to real frags(geom objects)
    return frags