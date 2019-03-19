#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains all input-output functions
"""
import numpy as np
from geom import geom
    
def read_xyz(fnm):
    """
    Parameters
    ----------
    fnm : string
        name or path of the .xyz file
        
    Returns
    -------
    geom 
        geometry object from the .xyz file
    """
    with open(fnm,"r") as f:
        rl=f.readlines()
        Natoms=int(rl[0])
        atoms=[]
        coords=[]
        for i in range(2,Natoms+2):
            atoms.append(rl[i].split()[0])
            coords.append(np.asarray([np.float64(i) for i in rl[i].split()[1:]]))
        Geom=geom(np.array(atoms), np.array(coords))
    return Geom

def read_zr(fnm):
    """
    Parameters
    ----------
    fnm : string
        name or path of the .zr file
        
    Returns
    -------
    geomA 
        geometry object of fragment A from the .zr file
     geomB 
        geometry object of fragment B from the .zr file   
    """
    with open(fnm,"r") as f:
        rl=f.readlines()
        Natoms=len(rl)-1
        atoms=[[],[]]
        coords=[[],[]]
        AorB=0
        for i in range(0,Natoms+1):
            if rl[i]=="----\n":
                AorB=1
            else:
                atoms[AorB].append(rl[i].split()[0])
                coords[AorB].append(np.asarray([np.float64(i) for i in rl[i].split()[1:]]))
        GeomA=geom(np.array(atoms[0]), np.array(coords[0]))
        GeomB=geom(np.array(atoms[1]), np.array(coords[1]))
    return GeomA,GeomB

def read_charge_txt(fnm):
    """
    Note
    ----
    The file should have the charges as 1 value per line, as many values as the atoms.
    Handles empty lines at the end.
    
    Parameters
    ----------
    fnm : str
        name or path of the text file to read.
    
    Returns
    -------
    list[floats]
        list of the point atomic charges in the file
    """
    with open(fnm,"r") as f:
        rl=f.readlines()
        charge_list=[]
        i=0
        while len(rl[i].split())==1:
            charge_list.append(float(rl[i]))
            i+=1
            if i==len(rl):
                break
    return charge_list

def write_xyz(g, fnm, which="inp"):
    """
    Note
    ----
    writes a geometry as a .xyz file
    
    Parameters
    ----------
    g : geom
        geometry to write
    fnm : str
        name or path of the output file
    which : {"inp", "com"}
        coordinates to write
    """
    with open(fnm,"w") as out:
        out.write(" "+str(len(g.atoms))+"\n\n")
        for i in range(len(g.atoms)):
            out.write(g.atoms[i]+"    "+"    ".join(map(str,g.coords(which)[i]))+"\n")   

def write_zr(gA, gB, fnm):
    """
    Note
    ----
    writes 2 geomtries as a .zr file. inp_coords of both gA and gB are taken!
    
    Parameters
    ----------
    gA : geom
        geometry of fragment A
    BA : geom
        geometry of fragment B
    fnm : str
        name or path of the output file
    """
    with open(fnm,"w") as out:
            for i in range(len(gA.atoms)):
                out.write(gA.atoms[i]+"    "+"    ".join(map(str,gA.inp_coords[i]))+"\n")
            out.write("----\n")
            for i in range(len(gB.atoms)):
                out.write(gB.atoms[i]+"    "+"    ".join(map(str,gB.inp_coords[i]))+"\n")
 
def write_frag_file(fnm,*args,Type="calculate"):
    """
    Note
    ----
    Writes fragments separated by "--\n"
    Parameters
    ----------
    fnm : str
        name or path of the fragment file to write
    Type : {"calculate", "list", "individual"}
        default is calculate
        calculate=> args=geom(all fragments)
        list=> args=[frag1,frag2...]
        individual=> args=frag1,frag2...
    """
    if Type=="calculate":
        frag_list=get_fragments(args[0])
    elif Type=="list":
        frag_list=args[0]
    elif Type=="individual":
        frag_list=args
    else:
        print("Go home bro, you're drunk")
    sl=[]
    for i in frag_list:
        sl.append("\n".join([i.atoms[j]+"    "+"    ".join(map(str,i.inp_coords[j])) for j in range(len(i.atoms))]))
    with open(fnm,"w") as out:
        out.write(("\n--\n").join(sl))
        