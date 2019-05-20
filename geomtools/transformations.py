#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains all geometric transformations: translations, rotations,etc...
"""
import numpy as np
from geomtools.geom import geom

def find_center(geom):
    """
    Parameters
    ----------
    g : geom
        geometry to find the center of
        
    Returns
    -------
    array(3)
        center of mass of the geometry
    """
    import mendeleev as md
    tot=0.0
    num=np.zeros(3)
    for i in range(len(geom.atoms)):
        atom_contrib=md.element(geom.atoms[i]).atomic_weight
        tot=np.add(tot,atom_contrib)
        num=np.add(num,np.multiply(atom_contrib,geom.inp_coords[i]))
    return np.divide(num,tot)

def transl(g, vect, which="inp"):
    """
    Parameters
    ----------
    g : geom
        geometry to apply the transition to
    vect : array(3)
        translation vector to apply
    which : {"inp", "com"}
        coordinates to translate, default is inp
    Returns
    -------
    geom
        geometry object where inp_coords are the results of the translation
    """
    return geom(g.atoms, np.add(g.coords(which),vect))

def uvec(vec):
    """
    Parameters
    ----------
    vec : array(3)
        vector
    
    Returns
    -------
    array(3)
        unitary vector with the same direction as "vec"
    """
    return np.divide(vec,np.linalg.norm(vec))

def get_angle(v1, v2, unit="deg"):
    """
    Parameters
    ----------
    v1: array(3)
        vector1
    v2: array(3)
        vector2
    unit: {"deg", "rad"}
        desired unit for the angle, default is "deg"
    
    Returns
    -------
    float
        angle between v1 and v2        
    """
    dot = np.dot(v1, v2)
    den = np.multiply(np.linalg.norm(v1), np.linalg.norm(v2))
    frac=np.divide(dot, den)
    if frac<-1.0:
        frac=-1.0
    angle = np.arccos(frac)
    if unit=="deg":
        angle=np.divide(np.multiply(angle,180),np.pi)
    elif unit=="rad":
        pass
    else:
        print("This unit for angles is not implemented yet. Why don't you do it, champ?")
    return angle

def transl_back(g2A, g2B, g1):     #todo: Maybe redo with Kabsch?
    """
    Parameters
    ----------
    g2A : geom
        geometry of fragment A in the roto-translated complex
    g2B : geom
        geometry of fragment B in the roto-translated complex
    g1 : geom
        original/reference geometry of fragment B
    Returns
    -------
    tuple
        (geometry of fragment A translated so that geom(B)==g1, geometry of B (translated back, i.e. geom(B)==g1))
    """
    if (g1.atoms!=g2B.atoms).all():
        print("The two geometries do not have the same atoms!!")    #Check if geometries have the same atom list
    else:
        g2Bt=transl(g2B,np.subtract(g1.get_com(),g2B.get_com()))
        g2At=transl(g2A,np.subtract(g1.get_com(),g2B.get_com()))
    return (g2At,g2Bt)

def rot_ax(g, ax, x, which = "com", angle_unit = "deg"): #todo: check if np.asmatrix can be avoided, implement for vector
    """
    Parameters
    ----------
    g : geom
        geometry to rotate
    ax : array(3)
        axis(vector) to rotate around
    x : float
        angle to rotate of, in deg or grad
    which : {"com", "inp"}
        coordinates to rotate, default is com
    angle_unit : {"deg", "grad"}
        unit for the angle, default is degree
    Returns
    geom
        rotated geometry
    -------
    """
    if np.linalg.norm(ax)!=1:
        ax=uvec(ax)
    if angle_unit=="deg":
        x=np.divide(np.multiply(x,np.pi),180)
    elif angle_unit=="rad":
        pass
    else:
        print("This unit for angles is not implemented yet. Why don't you do it, champ?")
    out=np.add(np.multiply(g.coords(which),np.cos(x)),
    np.add((np.multiply(np.asmatrix(ax).T,np.multiply(np.asmatrix(np.dot(g.coords(which),ax)),np.subtract(1,np.cos(x))))).T,
    np.multiply(np.cross(ax,g.coords(which)),np.sin(x))))
    return geom(g.atoms,np.asarray(out))

def align_vec(v1,v2,unit="deg"):
    """
    Parameters
    ----------
    v1: array(3)
        vector1
    v2: array(3)
        vector2
    unit: {"deg", "rad"}
        desired unit for the angle, default is "deg"
    
    Returns
    -------
    tuple
        ax: axis for the rotation, x: angle for the rotation
    """
    x = get_angle(v1, v2,unit=unit)
    ax = np.cross(v1,v2)
    return (ax,x)
    
def find_rot(g2, g1, thresh=0.00000001,angle_unit="deg"):
    """
    Parameters
    ----------
    g2 : geom
        final geometry
    g1 : geom
        initial geometry
    thresh : float
        threshold distance below which two points are considered equivalent, default is 0.00000001
    angle_unit :{"deg", "grad"}
        unit for the angles, default is degree
        
    Returns
    -------
    tuple
        (angle for the first rotation, angle for the second rotation, axis for the first rotation, axis for the second rotation)
    """
    if (g1.atoms!=g2.atoms).all():    #Check if geometries have the same atom list
        print("The two geometries do not have the same atoms!!")
    else:
        first=0
        add=1
        while np.linalg.norm(g1.get_com_coords()[first]) < thresh : #Skip atoms if they are "on" the center of mass
            first+=1
        #Skip atoms if they are "on" the center of mass
        while np.linalg.norm(g1.get_com_coords()[first+add]) < thresh \
        or np.linalg.norm(np.cross(g1.get_com_coords()[first],g1.get_com_coords()[first+add])) < thresh: #Skip aligned atoms
            add+=1 
        A1, A2 = uvec(g1.get_com_coords()[first]), uvec(g2.get_com_coords()[first])
        cross1=np.cross(A1,A2)
        #TODO check if they are very close
        if np.linalg.norm(cross1) < thresh: #if symmetric with respect to c.o.m. 
            Th_a = np.pi
            Va = uvec(np.cross(A1,g1.get_com_coords()[first+1]))
        else:
            Va = uvec(cross1)
            Th_a=np.arccos(np.dot(A1,A2))
        
        g1or1 =rot_ax(g1, Va, Th_a) #geometry after first rotation
        if (g1or1.get_com_coords() != g2.get_com_coords()).all():        
            A1r=g1or1.get_com_coords()[first] #A1 after first rotation
            s1, s2 = uvec(g1or1.get_com_coords()[first+1]), uvec(g2.get_com_coords()[first+1])  #auxiliary vectors to obtain B
            B1, B2 = uvec(np.cross(A1r,s1)), uvec(np.cross(A2,s2))  #vector B, made so to be perpendicular to A
            cross2=np.cross(B1,B2)
            if np.linalg.norm(cross2) < thresh:
                Vb=uvec(A1r)
                Th_b=np.pi
            else:
                Vb = uvec(cross2)
                Th_b=np.arccos(np.dot(B1,B2))
        else: #if one rotation is enough, then second rotation is 0°
            Vb=Va
            Th_b=0
        if angle_unit=="deg":   #adjust units
            Th_a=np.divide(np.multiply(Th_a,np.pi),180)
            Th_b=np.divide(np.multiply(Th_b,np.pi),180)
    return (Th_a, Th_b, Va, Vb)

def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two vectors
    """
    D = len(V[0])
    N = len(V)
    result = 0.0
    for v, w in zip(V, W):
        result += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(result/N)


def kabsch_rmsd(g1, g2,which="com"):
    """
    Rotate matrix P unto Q using Kabsch algorithm and calculate the RMSD.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    which : {"inp", "com"}
        coordinates to use, default is com

    Returns
    -------
    rmsd : float
        root-mean squared deviation
    """
    g1 = kabsch_rotate(g1, g2,which)
    return rmsd(g1.coords(which), g2.coords(which))


def kabsch_rotate(g1, g2,which="com"):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    which : {"inp", "com"}
        coordinates to use, default is com
        
    Returns
    -------
    P : array
        (N,D) matrix, where N is points and D is dimension,
        rotated

    """
    U = kabsch(g1, g2,which)

    # Rotate P
    P = np.dot(g1.coords(which), U)
    return geom(g1.atoms,P)


def kabsch(g1, g2, which="com"):
    """
    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.

    The algorithm works in three steps:
    - a centroid translation of P and Q (assumed done before this function
      call)
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    which : {"inp", "com"}
        coordinates to use, default is com
        
    Returns
    -------
    U : matrix
        Rotation matrix (D,D)
    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(g1.coords(which)), g2.coords(which))

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U

def apply_planar_symm(g,plane):
    s=geom(g.atoms,np.zeros([len(g.atoms),3]))
    for n,atom in enumerate(g.inp_coords):
        s.inp_coords[n]=plane.symmetric_point(atom)
    return s

def apply_point_symm(g,point):
    s=geom(g.atoms,np.zeros([len(g.atoms),3]))
    for n,atom in enumerate(g.inp_coords):
        s.inp_coords[n]=np.subtract(np.multiply(2,point),atom)
    return s    