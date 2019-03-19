#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains anything related to analysis of distances, ordering, and combining/dividing fragments
"""
import numpy as np
import sys
from geom import geom

def dist_order(g, which="com"):
    """
    Note
    ----
    d[order] yields the distances in increasing order
    
    Parameters
    ----------
    g : geom
        geometry to reorder
    which : {"com", "inp"}
        coordinates to order, default is com
    
    Returns
    -------
    tuple
        (array of distances in the order of g.atoms, order of atom indexes from the closest to the furthest)
    """
    d=np.sqrt(np.add(np.add((g.coords(which)[:,0])**2,(g.coords(which)[:,1])**2),(g.coords(which)[:,2])**2))
    order=d.argsort()
    return (d, order)

def check_geoms_order(g1, g2, thresh = 0.00005):
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
    d1, o1 = dist_order(g1)
    d2, o2 = dist_order(g2)
    if (g1.atoms[o1]!=g2.atoms[o1]).any() or (d1[o1] - d2[o1] > thresh).any():
        return False
    else:
        ncm=np.add(np.multiply(g1.get_com_coords()[o1][0],np.multiply(3,thresh)),np.multiply(np.multiply(1.5,thresh),np.random.rand(3)))
        g1a, g2a = g1.transl(ncm, which="com"), g2.transl(ncm, which="com") 
        da1, oa1 = dist_order(g1a,which="inp")
        da2, oa2 = dist_order(g2a, which="inp")
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
    d1, o1 = dist_order(g1)
    d2, o2 = dist_order(g2)
    if (g1.atoms[o1]!=g2.atoms[o1]).any() or (d1[o1] - d2[o1] > thresh).any():
        print("These are not the same geometry!!!")
        sys.exit(0)
    else:
        ncm=np.add(np.multiply(g1.get_com_coords()[o1][0],np.multiply(3,thresh)),np.multiply(np.multiply(1.5,thresh),np.random.rand(3)))
        g1a, g2a = geom(g1.atoms,g1.get_com_coords() - ncm), geom(g2.atoms,g2.get_com_coords() - ncm) #todo: use transl
        da1, oa1 = dist_order(g1a,which="inp")
        da2, oa2 = dist_order(g2a, which="inp")   
        return (geom(g1.atoms[oa1],g1.coords("inp")[oa1]),geom(g2.atoms[oa2],g2.coords("inp")[oa2]))

def get_cov_radius(El):
    """
    Parameters
    ----------
    El : str
        element whose radius is desired
        
    Returns
    -------
    float
        the first radius available from [Bragg, Slater, Cordero, Pyykko]
        
    """
    import mendeleev as md
    if md.element(El).covalent_radius_bragg!=None:
        return np.multiply(0.01,md.element(El).covalent_radius_bragg)
    elif md.element(El).covalent_radius_slater!=None:
        return np.multiply(0.01,md.element(El).covalent_radius_slater)
    elif md.element(El).covalent_radius_cordero!=None:
        return np.multiply(0.01,md.element(El).covalent_radius_cordero)
    elif md.element(El).covalent_radius_pyykko!=None:
        return np.multiply(0.01,md.element(El).covalent_radius_pyykko)
    else:
        print("Cannot find radius!")
        return None

def get_ordd_dist_mat(g):
    """
    Parameters
    ----------
    g : geom
        geometry to obtain the distance matrix of
    
    Returns
    -------
    tuple
        (density matrix, ordering of the matrix rows and columns)
        the density matrix is ordered so that the first line's values are in increasing order
    """
    from scipy.spatial import distance_matrix
    d, o = dist_order(g)
    DM1=distance_matrix(g.inp_coords,g.inp_coords)
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
    for i in range(len(row_ord)):
        for j in range(len(row_ord)):
            R_i, R_j = np.multiply(0.01,md.element(g.atoms[row_ord][i]).vdw_radius), np.multiply(0.01,md.element(g.atoms[row_ord][j]).vdw_radius)
            BM[i][j]=(DM[i][j] < R_i + R_j)
    return (BM, row_ord)

def comb_geoms(*args,List=False):
    """
    Parameters
    ----------
    List : bool
        whether the geometries will be given individually(False) or as a list(True), default is False
        False=> args=geom1,geom2...
        True=> args=[geom1,geom2...]
    
    Return
    ------
    geom
        combined geometry
    """
    if List:
       g_comb  = args[0][0]
       for i in args[0][1:]:
            g_comb = geom(np.append(g_comb.atoms,i.atoms),np.append(g_comb.inp_coords,i.inp_coords,axis=0))
    else:
        g_comb  = args[0]
        for i in args[1:]:
            g_comb = geom(np.append(g_comb.atoms,i.atoms),np.append(g_comb.inp_coords,i.inp_coords,axis=0))
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
    j=0
    f=[]
    while len(list(it.chain.from_iterable(f))) < Natms:
        f.append([])
        f[j]=[min([i for i in range(Natms) if i not in list(it.chain.from_iterable(f))])]
        cnt=0
        while cnt<len(f[j]):
            for i in f[j][cnt:]:
                cnt=len(f[j])
                f[j].extend([z for z in bonds[bonds[:,0]==i][:,1] if z not in f[j]])
        j+=1
    frags=[]
    for i in f:
        frags.append(geom(g.atoms[o][i],g.inp_coords[o][i]))
    return frags