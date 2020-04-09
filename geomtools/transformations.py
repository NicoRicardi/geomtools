#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains all geometric transformations: translations, rotations,etc...
"""
import numpy as np
from geomtools.geom import geom
import sys

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
    tot = 0.0
    num = np.zeros(3)
    for n,i in enumerate(geom.atoms):
        if "x" in i or "X" in i:  # ignoring ghost atoms
            continue
        atom_contrib=md.element(i).atomic_weight
        tot=np.add(tot,atom_contrib)
        num=np.add(num,np.multiply(atom_contrib,geom.coords[n]))
    tot = 1 if tot==0 else tot  # if all ghosts avoid division by 0 and return [0,0,0]
    return np.divide(num,tot)

def transl(g, vect):
    """
    Parameters
    ----------
    g : geom
        geometry to apply the transition to
    vect : array(3)
        translation vector to apply
        
    Returns
    -------
    geom
        geometry object where inp_coords are the results of the translation
    """
    gt=g.copy()
    gt.transl(vect)
    return gt

def center(g):
    """
    Parameters
    ----------
    g : geom
        geometry to center
        
    Returns
    -------
    geom
        centered geometry
    """
    copy=g.copy()
    copy.center()
    return copy

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
    from geomtools.geom import check_same_atoms,check_same_unit
    check_same_atoms(g1,g2B)
    check_same_unit(g1,g2B)
    g2Bt=transl(g2B,np.subtract(g1.get_com(),g2B.get_com()))
    g2At=transl(g2A,np.subtract(g1.get_com(),g2B.get_com()))
    return (g2At,g2Bt)

def rot_ax(g, ax, x, angle_unit = "deg"): #todo: check if np.asmatrix can be avoided, implement for vector
    """
    Parameters
    ----------
    g : geom
        geometry to rotate
    ax : array(3)
        axis(vector) to rotate around
    x : float
        angle to rotate of, in deg or rad
    angle_unit : {"deg", "rad"}
        unit for the angle, default is degree
        
    Returns
    -------
    geom
        rotated geometry
    """
    c=g.copy()
    c.rot_ax(ax, x, angle_unit = angle_unit)
    return c

def rot_matrix(g,M):
    """
    Rotate geometry g with matrix M

    Parameters
    ----------
    g : geometry object
        geometry to rotate
    M : np.array(3,3)
        rotation matrix
        
    Returns
    -------
    geom
        rotated geometry
    """
    c=g.copy()
    c.rot_matrix(M)
    return c

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
    Note
    ----
    Works if geometries are exactly the same. Consider using kabsch_find_rot
    
    Parameters
    ----------
    g2 : geom
        final geometry
    g1 : geom
        initial geometry
    thresh : float
        threshold distance below which two points are considered equivalent, default is 0.00000001
    angle_unit :{"deg", "rad"}
        unit for the angles, default is degree
        
    Returns
    -------
    tuple
        (angle for the first rotation, angle for the second rotation, axis for the first rotation, axis for the second rotation)
    """
    from geomtools.geom import check_same_atoms,check_same_unit
    check_same_atoms(g1,g2)
    check_same_unit(g1,g2)
    first=0
    add=1
    while np.linalg.norm(g1.coords[first]) < thresh : #Skip atoms if they are "on" the center of mass
        first+=1
    #Skip atoms if they are "on" the center of mass
    while np.linalg.norm(g1.coords[first+add]) < thresh \
    or np.linalg.norm(np.cross(g1.coords[first],g1.coords[first+add])) < thresh: #Skip aligned atoms
        add+=1 
    A1, A2 = uvec(g1.coords[first]), uvec(g2.coords[first])
    cross1=np.cross(A1,A2)
    #TODO check if they are very close
    if np.linalg.norm(cross1) < thresh: #if symmetric with respect to c.o.m. 
        Th_a = np.pi
        Va = uvec(np.cross(A1,g1.coords[first+1]))
    else:
        Va = uvec(cross1)
        Th_a=np.arccos(np.dot(A1,A2))
    
    g1or1 =rot_ax(g1, Va, Th_a,angle_unit="rad") #geometry after first rotation
    if (g1or1.coords != g2.coords).all():        
        A1r=g1or1.coords[first] #A1 after first rotation
        s1, s2 = uvec(g1or1.coords[first+1]), uvec(g2.coords[first+1])  #auxiliary vectors to obtain B
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
    Note
    ----
    Calculate Root-mean-square deviation between two geometries.
    NB reads both np.array(coords) or geom.

    Parameters
    ----------
    V : array or geom
        (N,D) matrix, where N is points and D is dimension.
    W : array or geom
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    float
        Root-mean-square deviation between the two vectors
    """
    for i in [V,W]: #TODO change!!!
        try:
            i=i.coords
        except:
            pass
    D = len(V[0])
    N = len(V)
    result = 0.0
    for v, w in zip(V, W):
        result += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(result/N)


def kabsch_rmsd(g1, g2):
    """
    Note
    ----
    Rotate geometry P unto Q using Kabsch algorithm and calculate the RMSD.
    NB reads both np.array(coords) or geom
    
    Parameters
    ----------
    P : array or geom
        (Natoms,3) coord array
    Q : array or geom
        (Natoms,3) coord array 
        
    Returns
    -------
    float
        root-mean squared deviation
    """
    for i in [g1,g2]: # TODO change!!
        try:
            i=i.coords
        except:
            pass
    g1 = kabsch_rotate(g1, g2)
    return rmsd(g1, g2)


def kabsch_rotate(g1, g2):
    """
    Note
    ----
    Rotate geometry P unto matrix Q using Kabsch algorithm.
    NB reads both np.array(coords) or geom
    
    Parameters
    ----------
    P : array or geom
        (Natoms,3) coord array
    Q : array or geom
        (Natoms,3) coord array
        
    Returns
    -------
    array(N,D) or geom
        (Natoms,3) coord array,
        rotated

    """
    c=[]
    for i in [g1,g2]:
        try:
            c.append(i.coords)
        except:
            c.append(i)
    U = kabsch_find_rot(c[0], c[1])

    # Rotate P
    try:
        out=rot_matrix(g1,U)
    except:
        out = np.dot(c[0], U)
    return out


def kabsch_find_rot(g1, g2):
    """
    Note
    ----
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
        (Natoms,3) coord array
    Q : array
        (Natoms,3) coord array
        
    Returns
    -------
    array(3,3)
        Rotation matrix (3,3)
    """
    try:
        g1=g1.coords
    except:
        pass
    try:
        g2=g2.coords
    except:
        pass
    # Computation of the covariance matrix
    C = np.dot(np.transpose(g1), g2)

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
    """
    Parameters
    ----------
    g : geom object
        the geometry to apply the transformation to
    plane : plane object
        the symmetry plane 
    
    Returns
    -------
    geom object
        the transformed geometry
    """
    s=geom(g.atoms,np.zeros([len(g.atoms),3]))
    for n,atom in enumerate(g.coords):
        s.coords[n]=plane.symmetric_point(atom)
    return s

def apply_point_symm(g,point):
    """
    Parameters
    ----------
    g : geom object
        the geometry to apply the transformation to
    point : array(3) or list of lenght 3
        coordinates of the inversion center 
    
    Returns
    -------
    geom object
        the transformed geometry
    """
    if type(point)==list:
        point=np.asarray(point)
    elif type(point)==np.ndarray:
        if len(point.shape)==1:
            pass
        elif len(point.shape)==2:
            if 1 in point.shape:
                a=point.shape.index(1)
                point.reshape(point.shape[0 if a==1 else 1])
            else:
               print("2D array!! only give the charges as a list, (n,) or (n,1) or (1,n) array!")
               sys.exit() #TODO change with exceptions
        else:
            raise geom.otherError("Array with more than 2 axes!! Give point as listor as array of shape: (3,),(3,1),(1,3)")
    s=geom(g.atoms,np.zeros([len(g.atoms),3]))
    for n,atom in enumerate(g.coords):
        s.coords[n]=np.subtract(np.multiply(2,point),atom)
    return s    