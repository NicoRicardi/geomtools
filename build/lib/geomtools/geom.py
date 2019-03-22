#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains the class "geom"
"""
import numpy as np
import sys

class geom:
    """
    Note
    ----
    Geometry object that contains atom labels and one or more set of coordinates for the molecule(s)
    """
    def __init__(self, atoms, inp_coords, coord_unit="Angstrom"):
        """
        Parameters
        ----------
        atoms : array(Natoms)
            Element symbols for the atoms
        inp_coords : array(3,Natoms) *check*
            x,y,z coordinates for each atom
        coord_unit : str
            unit for coordinates, defaults is Angstrom 
        """
        self.atoms = atoms
        self.inp_coords = inp_coords
        self.dipole_from_qm = {}
        self.coord_unit = coord_unit
        
    def get_com(self):
        """
        Note
        ----
            The center of mass is calculated at most once per geom object
        
        Returns
        -------
        array(3)
            the center of mass of the geometry
        """
        from transformations import find_center
        if not hasattr(self, "com"):
            self.com = find_center(self)
        return self.com
    
    def get_com_coords(self):
        """      
        Returns
        -------
        array(3,Natoms)
            translates the geometry so that the center of mass is in the origin
        """       
        if not hasattr(self, "com_coords"):
            self.com_coords = self.inp_coords - self.get_com()
        return self.com_coords    
    
    def add_charges(self, method, charge_list):
        """
        Note
        ----
            stores the atomic point charges obtained with a specific method
        
        Parameters
        ----------
        method : str
            method used to obtain the atomic point charges
        charge_list : list[floats]
            list of the atomic point charges        
        """
        if not hasattr(self,"charges"):
                self.charges = {}
        self.charges[method] = np.array(charge_list)
        
    def transl(self, vect, which = "inp"):
        """
        Parameters
        ----------
        vect : array(3)
            translation vector to apply
        which : {"inp", "com"}
            coordinates to translate, default is inp
            
        Returns
        -------
        geom
            geometry object where inp_coords are the results of the translation
        """
        return geom(self.atoms, np.add(self.coords(which), vect))
    
    def mix(self, g, ratio=0.5):
        """
        Note
        ----
        useful tool for "oscillating" geometry optimisations. NB uses input coordinates
        
        Parameters
        ----------
        g : geom
            geometry to mix with
        ratio : float
            ratio of self
            
        Returns
        -------
        geom
            geometry obtained from mixing self with g with the desired ratio 
        """
        from fragments import mix_geoms
        return mix_geoms(self,g,ratio)
    
    def coords(self, which):
        """
        Parameters
        ----------
        which : {"inp", "com"}
            "inp" for the input coordinates, "com" for coordinates translated to have the com at the origin
        
        Returns
        -------
        coords : array(3,Natoms)
            Either inp_coords or com_coords
        """
        if which=="com":
            return self.get_com_coords()
        if which=="inp":
            return self.inp_coords
        else:
            print("ERROR. Use either \"inp\" for input coordinates or \"com\" for coordinates centered at the center of mass")
            
    def add(self, g):
        """
        Parameters
        ----------
        g : geom
            geometry to add
        
        Returns
        -------
        geom
            geometry object containing self+g
        """
        from fragments import comb_geoms
        return comb_geoms(self, g)
    
    def get_charge_dipole(self, method, charge=0, dip_unit="au"):
        """
        Parameters
        ----------
        method : str
            method for the charges
        charge : int
            charge of the molecule, default is 0
        dip_unit : str
            unit for the dipole moment, default is au
        
        Returns
        -------
        array(3)
            dipole obtained from the atomic point charges
        """
        if not hasattr(self,"charge_dipoles"):
                    self.charge_dipoles = {}
                    self.charge_dipoles_units={}
        self.charge_dipoles[method]=calculate_charge_dipole(self,method,charge=0,self.coord_unit,out=dip_unit)
        self.charge_dipoles_units[method]=dip_unit
        return self.charge_dipoles[method]

        
def calculate_charge_dipole(g, method, charge=0, coords="Angstrom", out="au"):
    """
    Note
    ----
    Handles different names au/a.u./bohr
    
    Parameters
    ----------
    g : geom
        geometry 
    method : str
        method identifying the atomic point charges
    charge : int
        the charge of the molecule(s), default is 0
    coords : {"Angstrom", "bohr", "au", "a.u."}
        unit of measure for the coordinates, default is Angstrom. NB case insensitive
    out : {"Debye", "au", "a.u."}
        desired unit for the output, default is au
    
    Returns
    -------
    array(3)
        dipole obtained from the atomic point charges
    """
    if charge!=0:
        coc=np.divide(np.sum(g.inp_coords,axis=0),np.multiply(charge,len(g.atoms)))
    else:
        coc=np.zeros(3)
    d=np.sum(np.multiply(g.charges[method].reshape(len(g.atoms),1),g.inp_coords-coc),axis=0)
    dict_={"au":"au", "a.u.":"au", "bohr":"au", "angstrom":"angstrom", "debye":"debye"}
    if coords.lower() not in dict_.keys() or out not in dict_.keys():
        print("combination of units of measure not implemented yet. Why don't you do it, champ?")
        sys.exit(0)
    if dict_[coords.lower()]=="angstrom":
        d=np.divide(d,0.529177)
    if dict_[out.lower()]=="debye":
        d=np.divide(d,0.393456)
    return d 