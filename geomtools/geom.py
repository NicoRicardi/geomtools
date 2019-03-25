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
        self.coord_unit = coord_unit
        
    def from_xyz(fnm):
        """
        Parameters
        ----------
        fnm : str
            name or path to the .xyz file
        
        Returns
        -------
        geom
            geometry object from the .xyz file
        """
        from geomtools.io import read_xyz
        return read_xyz(fnm)
    
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
        from geomtools.transformations import find_center
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
        from geomtools.fragments import mix_geoms
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
        from geomtools.fragments import comb_geoms
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
        self.charge_dipoles[method]=calculate_charge_dipole(self,method,coords=self.coord_unit,out=dip_unit,charge=0)
        self.charge_dipoles_units[method]=dip_unit
        return self.charge_dipoles[method]
    
    def change_coord_unit(self,out):
        """
        Note
        ----
        Changes the unit of inp_coords
        
        Parameters
        ----------
        out : str
            desired unit for the output. Handles cases properly and au/a.u./bohr
        """
        dict_={"au":"au", "a.u.":"au", "bohr":"au", "angstrom":"angstrom"}
        if dict_[self.coord_unit.lower()]==dict_[out.lower()]:
            print("The coordinates are already in "+out)
        elif dict_[self.coord_unit.lower()]=="angstrom" and dict_[out.lower()]=="au":
            self.inp_coords=np.multiply(1,88973,self.inp_coords)
            self.coord_unit="au"
            print("coordinates changes from "+self.coord_unit+"to "+out)
        elif dict_[self.coord_unit.lower()]=="au" and dict_[out.lower()]=="angstrom":
            self.inp_coords=np.multiply(0.529177,self.inp_coords)
            self.coord_unit="angstrom"
            print("coordinates changes from "+self.coord_unit+"to "+out)
        else:
            print("Unit combination not implemented yet. Why don't you do it?")
            
    def reorder(self,o):
        """
        Parameters
        ----------
        o : list
            new order of the geometry
        
        Returns
        -------
        geom
            geometry in the new order
        """
        return geom(self.atoms[o],self.inp_coords[o])
    
    def xyz(self,fnm, which="inp", decimals=6, spaces=4):
        """
        Does not work yet. Do not know why
        """
        from geomtools.io import write_xyz
        write_xyz(self,which,decs=decimals,spacing=spaces)
           
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

class plane:
    """
    Note
    ----
    Plane object that stores the parameters for the plane equation
    """
    def __init__(self,A,B,C,D):
        """
        Parameters
        ----------
        A : float 
            or anything that can be turned into a float
        B : float 
            or anything that can be turned into a float
        C : float 
            or anything that can be turned into a float
        D : float 
            or anything that can be turned into a float
        Note
        ----
        Uses equation form Ax + By + Cz + D = 0
        """
        self.A=float(A)
        self.B=float(B)
        self.C=float(C)
        self.D=float(D)
        
    def eq(self):
        """
        Note
        ----
        prints the equation of the plane object
        """
        print("{:.2f}*x + {:.2f}*y + {:.2f}*z + {:.2f} = 0".format(self.A,self.B,self.C,self.D))
    
    @classmethod
    def normal_to(cls,v,P=[0,0,0]):
        """
        Parameters
        ----------
        v : array(3)
            vector the plane needs to be perpendicular to
        P : array(3)
            point that the plane needs to pass by. default is the origin
        
        Returns
        ------
        plane object perpendicular to v and passing by P
        """
        A=v[0]
        B=v[1]
        C=v[2]
        D=-A*P[0]-B*P[1]-C*P[2]
        return cls(A,B,C,D)
    
    @classmethod
    def by_three_points(cls,pointarray):
        """
        Parameters
        ----------
        pointarray : array(3,3)
        
        Returns
        ------
        plane object passing by the three points
        """
        n=np.cross(pointarray[1]-pointarray[0],pointarray[1]-pointarray[2])
        A=n[0]
        B=n[1]
        C=n[2]  
        P=pointarray[1]
        D=-A*P[0]-B*P[1]-C*P[2]
        return cls(A,B,C,D)
    
    def get_norm(self):
        """
        Returns
        -------
        array(3)
            the norm of the plane
        """
        if not hasattr(self,"norm"):
            self.norm=[self.A,self.B,self.C]
        return self.norm
    
    def get_unorm(self):
        """
        Returns
        -------
        array(3)
            the norm of the plane, as a unitary vector
        """
        from geomtools.transformations import uvec
        if not hasattr(self,"unorm"):
            self.unorm=uvec([self.A,self.B,self.C])
        return self.unorm
        
    def dist_from_point(self,p,signed=False):
        """
        Parameters
        ----------
        p : array(3)
            point to calculate the distance from
        signed : bool
            wether the distance should be signed or absolute
        
        Returns
        -------
        float
            the distance from p to the plane
        """
        num=np.add(np.add(np.add(np.multiply(self.A,p[0]),np.multiply(self.B,p[1])),np.multiply(self.C,p[2])),self.D)
        if not signed:
            num=abs(num)
        den=np.linalg.norm([self.A,self.B,self.C])
#        d=np.divide(np.add(np.add(np.add(np.multiply(self.A,p[0]),np.multiply(self.B,p[1])),np.multiply(self.C,p[2])),self.D),np.linalg.norm([self.A,self.B,self.C]))
        return np.divide(num,den)
   
    def symmetric_point(self,p):
        """
        p : array(3)
            point we need the symmetric of
        
        Returns
        -------
        array(3)
            the symmetric point
        """
        return np.add(p,np.multiply(self.get_unorm(),np.multiply(-2,self.dist_from_point(p,signed=True))))
