#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains all input-output functions
"""
import numpy as np
from geomtools.geom import geom
    
def read_xyz(fnm, identifier=""):
    """
    Parameters
    ----------
    fnm : string
        name or path of the .xyz file
    identifier: str or array(Natoms)
        identifier for the geom object
            
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
        Geom=geom(np.array(atoms), np.array(coords), identifier=identifier, coord_unit="Angstrom", charges_dict={})
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
        GeomA=geom(np.array(atoms[0]), np.array(coords[0]), identifier="-A", coord_unit="Angstrom", charges_dict={})
        GeomB=geom(np.array(atoms[1]), np.array(coords[1]), identifier="-B",coord_unit="Angstrom", charges_dict={})
    return GeomA,GeomB

def read_coords(fnm, identifier="", inp="Angstrom", out="Angstrom"):
    """
    Parameters
    ----------
    fnm : string
        name or path of the coordinate file (any coord-only file, as xyz without header, in any unit)
    identifier: str or array(Natoms)
        identifier for the geom object
    inp : {"Angstrom","au","a.u.","bohr"}
        unit of the input, default is Angstrom. NB case insensitive
    out : {"Angstrom","au","a.u.","bohr"}
        unit of the output, default is Angstrom. NB case insensitive
        
    Returns
    -------
    geom 
        geometry object from the coord file
    """
    dict_ = {"angstrom":"angstrom","au":"au","a.u.":"au","bohr":"au"}
    with open(fnm,"r") as f:
        rl=f.readlines()
        Natoms=int(rl[0])
        atoms=[]
        coords=[]
        for i in range(Natoms):
            atoms.append(rl[i].split()[0])
            coords.append(np.asarray([np.float64(i) for i in rl[i].split()[1:]]))
        if dict_[inp.lower()]==dict_[out.lower()]:
            pass
        elif dict_[inp.lower()]=="au" and dict_[out.lower()]=="angstrom":
            coords=np.multiply(0.529177,coords)
        elif dict_[inp.lower()]=="angstrom" and dict_[out.lower()]=="au":
            coords=np.divide(coords,0.529177)
        else:
            print("combination of units of measure not implemented yet. Why don't you do it, champ?")
        Geom=geom(np.array(atoms), np.array(coords), identifier=identifier, coord_unit=out, charges_dict={})
    return Geom

def read_string(coord_string, identifier="", inp="Angstrom", out="Angstrom"):
    """
    Parameters
    ----------
    coord_string : string
        string with the coordinates
    identifier: str or array(Natoms)
        identifier for the geom object
    inp : {"Angstrom","au","a.u.","bohr"}
        unit of the input, default is Angstrom. NB case insensitive
    out : {"Angstrom","au","a.u.","bohr"}
        unit of the output, default is Angstrom. NB case insensitive
        
    Returns
    -------
    geom 
        geometry object from the coord file
    """
    import sys
    if sys.version_info[0] < 3:
        from StringIO import StringIO
    else:
        from io import StringIO
    tofeed = StringIO(coord_string)
    return read_coords(tofeed, identifier=identifier, inp=inp, out=out)
    
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

def write_xyz(g, fnm, decs=6, spacing=4):
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
    decs : int
        number of desired decimal digits
    spacing : int
        number of empty spaces between coordinates.NB: the "-" sign will take one of these spaces, so values <2 are unadvisable
    """
    with open(fnm,"w") as out:
        out.write(" "+str(len(g.atoms))+"\n\n")
        for i in range(len(g.atoms)):
            out.write(g.atoms[i]+' {:{w}.{p}f} {:{w}.{p}f} {:{w}.{p}f}\n'.format(g.coords[i][0],g.coords[i][1],g.coords[i][2],w=spacing+decs+1, p=decs))
            
def write_zr(gA, gB, fnm, decs=6, spacing=4):
    """
    Note
    ----
    writes 2 geomtries as a .zr file.
    
    Parameters
    ----------
    gA : geom
        geometry of fragment A
    BA : geom
        geometry of fragment B
    fnm : str
        name or path of the output file
    decs : int
        number of desired decimal digits
    spacing : int
        number of empty spaces between coordinates.NB: the "-" sign will take one of these spaces, so values <2 are unadvisable
    """
    with open(fnm,"w") as out:
            for i in range(len(gA.atoms)):
                out.write(gA.atoms[i]+' {:{w}.{p}f} {:{w}.{p}f} {:{w}.{p}f}\n'.format(gA.coords[i][0],gA.coords[i][1],gA.coords[i][2],w=spacing+decs+1, p=decs))
            out.write("----\n")
            for i in range(len(gB.atoms)):
                out.write(gB.atoms[i]+' {:{w}.{p}f} {:{w}.{p}f} {:{w}.{p}f}\n'.format(gB.coords[i][0],gB.coords[i][1],gB.coords[i][2],w=spacing+decs+1, p=decs))
            

def write_coords(g, fnm, inp="Angstrom", out="Angstrom", decs=6, spacing=4):
    """
    Note
    ----
    writes a geometry as coord-only file, in any unit
    
    Parameters
    ----------
    g : geom
        geometry to write
    fnm : str
        name or path of the output file
    inp : {"Angstrom","au","a.u.","bohr"}
        unit of the input, default is Angstrom. NB case insensitive
    out : {"Angstrom","au","a.u.","bohr"}
        unit of the output, default is Angstrom. NB case insensitive
    decs : int
        number of desired decimal digits
    spacing : int
        number of empty spaces between coordinates.NB: the "-" sign will take one of these spaces, so values <2 are unadvisable
    """
    dict_ = {"angstrom":"angstrom","au":"au","a.u.":"au","bohr":"au"}
    if dict_[inp.lower()]==dict_[out.lower()]:
        factor=1
    elif dict_[inp.lower()]=="au" and dict_[out.lower()]=="angstrom":
        factor=0.529177
    elif dict_[inp.lower()]=="angstrom" and dict_[out.lower()]=="au":
        factor=1.88973
    else:
        print("combination of units of measure not implemented yet. Why don't you do it, champ?")
    coords_to_print=np.multiply(factor,g.coords)
    with open(fnm,"w") as out:
        for i in range(len(g.atoms)):
             out.write(g.atoms[i]+' {:{w}.{p}f} {:{w}.{p}f} {:{w}.{p}f}\n'.format(coords_to_print[i][0],coords_to_print[i][1],coords_to_print[i][2],w=spacing+decs+1, p=decs))
           
           
def write_frag_file(fnm, *args, Type="calculate", decs=6, spacing=4):
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
    decs : int
        number of desired decimal digits
    spacing : int
        number of empty spaces between coordinates.NB: the "-" sign will take one of these spaces, so values <2 are unadvisable
    """
    from fragments import get_fragments
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
        sl.append("\n".join([i.atoms[j]+' {:{w}.{p}f} {:{w}.{p}f} {:{w}.{p}f}'.format(i.coords[j][0],i.coords[j][1],i.coords[j][2],w=spacing+decs+1, p=decs) for j in range(len(i.atoms))]))
#        sl.append("\n".join([i.atoms[j]+"    "+"    ".join(map(str,i.coords[j])) for j in range(len(i.atoms))]))
    with open(fnm,"w") as out:
        out.write(("\n--\n").join(sl))
        