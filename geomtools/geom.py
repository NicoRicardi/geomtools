#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains the class "geom"
"""
import numpy as np
#import sys

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class NatomsError(Error):
    """Mismatch in the number of atoms.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self):
        self.message = "Mismatch in atom number!!"
        print(self.message)
        
class atomsError(Error):
    """Mismatch in the atom labels.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self):
        self.message = "Mismatch in atom labels!!"
        print(self.message)
        
class unitError(Error):
    """Mismatch in the units.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self):
        self.message = "Mismatch in unit!!"
        print(self.message)
        
class sizemismatchError(Error):
    """Mismatch in the units.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self,message):
        self.message = message
        print(self.message)      
        
class otherError(Error):
    """Other error.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self,message):
        self.message = message
        print(self.message)

class geom:
    """
    Note
    ----
    Geometry object that contains atom labels and one or more set of coordinates for the molecule(s)
    """
    def __init__(self, atoms, coords, identifier= "", coord_unit="Angstrom", charges_dict={}):
        """
        Parameters
        ----------
        atoms : array(Natoms)
            Element symbols for the atoms
        coords : array(Natoms,3)
            x,y,z coordinates for each atom
        identifier: str or array(Natoms)
            identifier for this molecule to add after every atom. Useful for pyscf, but use wisely!
        coord_unit : str
            unit for coordinates, defaults is Angstrom 
        charges_dict: dict
            keys=method(s) values=array of charges
        """
        if type(atoms) == str:
            atoms = np.array([atoms])
        if type(atoms) in [list, tuple]:
            atoms = np.array(atoms)
        self.atoms = atoms
        if identifier:
            print("Beware of risky identifiers ('a','b','c','d') and third party modules! 'C'+'a/A' could be read as calcium!")
        if type(identifier) not in [str, np.ndarray]:  # if non-standard type
            if type(identifier) in [tuple, list]:  # equivalent to array
                identifier = np.asarray(identifier)
            else:  # int or other stuff
                try:
                    identifier = str(identifier)  # try turn to string 
                except:
                    raise TypeError("identifier is an array(Natoms), lists/tuples can be converted, int/float/str are converted to array of strings")
        ##Now identifier is either np.ndarray or str
        if type(identifier) == np.ndarray and type(identifier[0]) != np.str_:  # array of int/float => array of str
            identifier = np.array([str(i) for i in identifier])
        if type(identifier) == str:  # str => array of str
            identifier = np.array([identifier for i in range(len(self.atoms))])
        self.identifier = identifier
        if type(coords) != np.ndarray:
            raise TypeError("coords must be a numpy array")
        if len(coords.shape) == 1:  # 1D array 
            coords = coords.reshape(-1,3)
        self.coords = coords
        if not (np.array([len(self.atoms),len(self.identifier)])== len(self.coords)).all():
            raise sizemismatchError("atoms, identifier, and coords must have the same lenght!")
        self.coord_unit = coord_unit.lower()
        self.charges = charges_dict
    
    def __str__(self, identifier = True, spacing=4, decimals=6):
        """
        Note
        ----
        Allows easy printing
        
        Parameters
        ----------
        spacing: int
            desired spacing between columns
        decimals: int
            number of decimal digits desired
        """
        if identifier:
            return "\n".join([self.atoms[i]+self.identifier[i]+' {:{w}.{p}f} {:{w}.{p}f} {:{w}.{p}f}'.format(self.coords[i][0],self.coords[i][1],self.coords[i][2],w=spacing+decimals+1, p=decimals) for i in range(len(self.atoms))])
        else:
            return "\n".join([self.atoms[i]+' {:{w}.{p}f} {:{w}.{p}f} {:{w}.{p}f}'.format(self.coords[i][0],self.coords[i][1],self.coords[i][2],w=spacing+decimals+1, p=decimals) for i in range(len(self.atoms))])
        
    def __repr__(self, spacing=4, decimals=6):
        """
        Note
        ----
        Allows easy calling
        
        Parameters
        ----------
        spacing: int
            desired spacing between columns
        decimals: int
            number of decimal digits desired
        """
        return "\n".join([self.atoms[i]+' {:{w}.{p}f} {:{w}.{p}f} {:{w}.{p}f}'.format(self.coords[i][0],self.coords[i][1],self.coords[i][2],w=spacing+decimals+1, p=decimals) for i in range(len(self.atoms))])
    
    def __add__(self, add):
        """
        Note
        ----
        Allows addition of atoms (e.g. geomC=geomA+geomB). Nb: it returns geomC!
        
        Parameters
        ----------
        add: geom or vector(array/list(3))
            geometry to add or vector to translate
        
        Returns
        -------
        geom
            the total geometry
        """
        if type(add) in [list, np.array]:
            gt=self.copy()
            gt.transl(add)
            return gt
        elif type(add) == geom:
            from geomtools.fragments import comb_geoms #TODO check identifier
            return comb_geoms(self,add)
        else:
            print(type(add))
            print(isinstance(add,geom))
            raise TypeError("Use + to combine geometry (e.g. g = gA + gB) or to translate (e.g. gt = g + [1,1,1])")
    def __iadd__(self, add):
        """
        Note
        ----
        Allows addition of atoms (e.g. geomA+=geomB). NB: it changes geomA, returns nothing!
        
        Parameters
        ----------
        g: geom
            geometry to add
        """
        if type(add) in [list, np.array]:
            self.transl(add)
            return self
        elif type(add) == geom:
            self.add_atoms(add)
            return self
        else:
            print(type(add))
            print(isinstance(add,geom))
            raise TypeError("Use + to combine geometry (e.g. g = gA + gB) or to translate (e.g. gt = g + [1,1,1])")
        
    def __sub__(self,vect):
        """
        Note
        ----
        Allows quick translation (e.g. geomAt=geomA-vect). Nb: it returns geomAt!
        
        Parameters
        ----------
        vect: array (or list)
            vector to subtract
        
        Returns
        -------
        geom
            the translated geometry
        """
        gt=self.copy()
        gt.transl(-vect)
        return gt
    
    def __isub__(self,vect):
        """
        Note
        ----
        Allows quick translation (e.g. geomA-=vect). Nb: it changes geomA and returns nothing!
        
        Parameters
        ----------
        vect: array (or list)
            vector to subtract
        """
        self.transl(-vect)
        return self
        
    def __eq__(self,g, thresh=1e-6):
        """
        Note
        ----
        Allows quick comparison (==) of geometries within threshold
        
        Parameters
        ----------
        g: geom
            geometry to compare to
        
        Returns
        -------
        bool
            whether they are equal or not
        """
        if type(self) != type(g):
            return False
        elif not self.have_same_unit(g):
            return False
        elif len(self.atoms) != len(g.atoms):#Check if geometries have the same atom list
            return False
        elif (self.atoms != g.atoms).all(): #Check if geometries have the same atom list
            return False
        else:
            return (abs(self.coords - g.coords) < thresh).all()
        
    def __neq__(self,g, thresh=1e-6):
        """
        Note
        ----
        Allows quick comparison (!=) of geometries within threshold
        
        Parameters
        ----------
        g: geom
            geometry to compare to
        
        Returns
        -------
        bool
            whether they are different or not
        """
        if type(self) != type(g):
            return True
        elif not self.have_same_unit(g):
            return True
        elif len(self.atoms) != len(g.atoms):#Check if geometries have the same atom list
            return True
        elif (self.atoms != g.atoms).all(): #Check if geometries have the same atom list
            return True
        else:
            return (abs(self.coords - g.coords) > thresh).any()
        
    def strip_identifier(self):
        """
        Note
        ----
        sets the identifier to array of empty strings
        """
        self.identifier = np.array(["" for i in range(len(self.atoms))])
            
    def add_identifier(self, identifier):
        """
        Note
        ----
        Sets the identifier.
        
        Parameters
        ----------
        identifier: str 
            desired identifier. if not string tries to convert
        """
        if (self.identifier!="").any():
            raise otherError("Geometry already has an identifier! use change_identifier!")
        print("Beware of risky identifiers ('a','b','c','d') and third party modules! 'C'+'a/A' could be read as calcium!")
        if type(identifier) not in [str, np.ndarray]:  # if non-standard type
            if type(identifier) in [tuple, list]:  # equivalent to array
                identifier = np.asarray(identifier)
            else:  # int or other stuff
                try:
                    identifier = str(identifier)  # try turn to string 
                except:
                    raise TypeError("identifier is an array(Natoms), lists/tuples can be converted, int/float/str are converted to array of strings")
        ##Now identifier is either np.ndarray or str
        if type(identifier) == np.ndarray and type(identifier[0]) != np.str_:  # array of int/float => array of str
            identifier = np.array([str(i) for i in identifier])
        if type(identifier) == str:  # str => array of str
            identifier = np.array([identifier for i in range(len(self.atoms))])
        if len(identifier) != len(self.atoms):
            raise sizemismatchError("your identifier does not match natoms!")
        self.identifier = identifier
        
    def change_identifier(self, identifier):
        """
        Note
        ----
        Changes the identifier.
        
        Parameters
        ----------
        identifier: str 
            desired identifier. if not string tries to convert
        """
        self.strip_identifier()
        self.add_identifier(identifier)
        
    def number_in_group(self):
        """
        Sets
        ----
        self.identifier
            sets it to the numbering for each element
        """
        idn = self.atoms.copy()
        elems, count = np.unique(idn, return_counts=True)
        for n in range(len(elems)):
            idn[idn==elems[n]] =  np.arange(1, count[n]+1)
        self.identifier = idn
        
    def ghostify(self):
        """
        Note
        ----
        Center of mass ignores ghost atoms (they have no mass)
        Sets
        ----
        self.atoms
            all atoms become their ghost version
        """
        self.atoms = np.array(["X"+i for i in self.atoms])
        
    def have_same_atoms(self,g):
        """
        Note
        ----
        Useful for scripts where different operations are to be done if geomA.atoms==geomB.atoms or geomA.atoms!=geomB.atoms
        
        Parameters
        ----------
        g: geom
            the geometry to compare with
            
        Returns
        -------
        bool
            whether they have the same atoms or not
        """
        if len(self.atoms) != len(g.atoms):#Check if geometries have the same atom list
            return False
        elif (self.atoms != g.atoms).all(): #Check if geometries have the same atom list
            return False
        else:
            return True
        
    def check_same_atoms(self,g):
        """
        Note
        ----
        Returns nothing! Just raises specific errors in case
        
        Parameters
        ----------
        g: geom
            the geometry to compare with
        """
        if len(self.atoms) != len(g.atoms):#Check if geometries have the same atom list
            raise NatomsError()
        elif (self.atoms != g.atoms).all(): #Check if geometries have the same atom list
            raise atomsError()
            
    def have_same_unit(self,g):
        """
        Note
        ----
        Useful for scripts where different operations are to be done if geomA.coord_unit==geomB.coord_unit or geomA.coord_unit!=geomB.coord_unit (e.g. convert)
        
        Parameters
        ----------
        g: geom
            the geometry to compare with
            
        Returns
        -------
        bool
            whether they have the same unit or not
        """
        dict_ = {"angstrom":"angstrom","au":"au","a.u.":"au","bohr":"au"}
        if dict_[self.coord_unit.lower()] != dict_[g.coord_unit.lower()]:
            return False
        else:
            return True
        
    def check_same_unit(self,g):
        """
        Note
        ----
        Returns nothing! Just raises specific errors in case
        
        Parameters
        ----------
        g: geom
            the geometry to compare with
        """
        dict_ = {"angstrom":"angstrom","au":"au","a.u.":"au","bohr":"au"}
        if dict_[self.coord_unit.lower()] != dict_[g.coord_unit.lower()]:
            raise unitError()
            
    def from_xyz(fnm, identifier=""):
        """
        Parameters
        ----------
        fnm : str
            name or path to the .xyz file
        identifier: str or array(Natoms)
            identifier for the geom object
        
        Returns
        -------
        geom
            geometry object from the .xyz file
        """
        from geomtools.geom_io import read_xyz
        return read_xyz(fnm, identifier=identifier)
    
    def from_coordfile(fnm, identifier="", inp="Angstrom", out="Angstrom"):
        """
        Parameters
        ----------
        fnm : str
            name or path to the .xyz file
        identifier: str or array(Natoms)
            identifier for the geom object
        inp : {"Angstrom","au","a.u.","bohr"}
            unit of the input, default is Angstrom. NB case insensitive
        out : {"Angstrom","au","a.u.","bohr"}
            unit of the output, default is Angstrom. NB case insensitive            
        
        Returns
        -------
        geom
            geometry object from the .xyz file
        """
        from geomtools.geom_io import read_coords
        return read_coords(fnm, identifier=identifier, inp=inp, out=out)  
    
    def from_string(coord_string, identifier="", inp="Angstrom", out="Angstrom"):
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
        from geomtools.geom_io import read_string
        return read_string(coord_string, identifier=identifier, inp=inp, out=out)    
    
    def get_com(self):
        """
        Note
        ----
            The center of mass is calculated at most once per geom object.
            Ghost atoms have no mass, hence they are ignored
        
        Returns
        -------
        array(3)
            the center of mass of the geometry
        """
        from geomtools.transformations import find_center #TODO check identifier
        if not hasattr(self, "com"):
            self.com = find_center(self)
        return self.com
    
    def get_com_coords(self):
        """      
        Note
        ----
            com_coords are calculated at most once per geom object.
            Ghost atoms have no mass, hence they are ignored
            
        Returns
        -------
        array(3,Natoms)
            the coordinates translated so that the center of mass is in the origin
        """       
        if not hasattr(self, "com_coords"):
            self.com_coords = self.coords - self.get_com()
        return self.com_coords    
    
    def center(self):
        """
        Note
        ----
            com_coords are calculated at most once per geom object.
            Ghost atoms have no mass, hence they are ignored
            
        Translates the geometry so that the center of mass is in the origin
        """
        self.coords = self.get_com_coords()
        self.com = np.array([0.0,0.0,0.0])
    
    def add_charges(self, method, charges):
        """
        Note
        ----
            stores the atomic point charges obtained with a specific method
        
        Parameters
        ----------
        method : str
            method used to obtain the atomic point charges
        charges : array or list/tuple[floats]
            array/list/tuple of the atomic point charges        
        """
        if type(charges) in [list,tuple]:
            charges=np.asarray(charges)
        elif type(charges) == np.ndarray:
            if len(charges.shape) == 1:
                pass
            elif len(charges.shape) == 2:
                if 1 in charges.shape:
                    a=charges.shape.index(1)
                    charges.reshape(charges.shape[0 if a==1 else 1])
                else:
                    raise otherError("2D array!! only give the charges as a list, (n,) or (n,1) or (1,n) array!")
            elif len(charges.shape)>2:
                raise otherError("Your array has 3 or more axis!! only give the charges as a list, (n,) or (n,1) or (1,n) array!")
        else:
            raise otherError("Weird type: only give the charges as a list, (n,) or (n,1) or (1,n) array!")
        if len(charges) != len(self.atoms):
            raise NatomsError()
        self.charges[method] = charges
        
    def copy(self):
        """
        Returns
        -------
        geom
            a copy of self
        """
        import copy as c
        return c.deepcopy(self)
    
    def transl(self, vect):
        """
        Parameters
        ----------
        vect : array(3)
            translation vector to apply
        """
        self.coords = np.add(self.coords, vect)
        
    def rot_matrix(self,M):
        """
        Rotate with matrix M. If the geometry is not centered it deletes com and com_coords
    
        Parameters
        ----------
        M : np.array(3,3)
            rotation matrix
            
        """
        self.coords = np.dot(self.coords, M)
        if hasattr(self,"com"):
            if not hasattr(self,"com_coords"):
                if not (self.com != np.zeros(3)).all():
                    del self.com
            elif not (self.coords==self.com_coords).all():
                del self.com
                del self.com_coords
        
    def rot_ax(self, ax, x, angle_unit = "deg"): #todo: check if np.asmatrix can be avoided, implement for vector
        """
        Parameters
        ----------
        ax : array(3)
            axis(vector) to rotate around
        x : float
            angle to rotate of, in deg or rad
        angle_unit : {"deg", "rad"}
            unit for the angle, default is degree
            
        Returns
        geom
            rotated geometry
        -------
        """
        if np.linalg.norm(ax)!=1:
            from geomtools.transformations import uvec
            ax=uvec(ax)
        if angle_unit=="deg":
            x=np.divide(np.multiply(x,np.pi),180)
        elif angle_unit=="rad":
            pass
        else:
            print("This unit for angles is not implemented yet. Why don't you do it, champ?")
        out=np.add(np.multiply(self.coords,np.cos(x)),
        np.add((np.multiply(np.asmatrix(ax).T,np.multiply(np.asmatrix(np.dot(self.coords,ax)),np.subtract(1,np.cos(x))))).T,
        np.multiply(np.cross(ax,self.coords),np.sin(x))))
        self.coords=np.asarray(out)
    
    def mix(self, g, ratio=0.5):
        """
        Note
        ----
        useful tool for "oscillating" geometry optimisations. 
        NB. it wipes charges!
        
        Parameters
        ----------
        g : geom
            geometry to mix with
        ratio : float
            ratio of self
        """
        self.check_same_atoms(g)
        self.check_same_unit(g)
        self.coords=np.add(np.multiply(ratio,self.coords),np.multiply(1-ratio,g.coords))
        self.charges = {}
    
    def inherit_optionals(self,g):
        """
        Note
        ----
        copies all attributes but atoms, coords, and coord_units from g. Useful for FF charges or custom attributes
        
        Parameters
        ----------
        g: geom
            geometry to inherit from
        """
        optionals=[i for i in self.__dict__.keys() if i not in ["atoms","coords", "coord_unit"]]
        for i in optionals:
            setattr(self,i,getattr(g,i))
        
    def add_atoms(self, g, keep_charges=False, keep_identifier=False):
        """
        Note
        ----
        Any custom attribute of the geometry you add is lost. 
        If keep_charges == True it will combine the charges of the two geometries for any method both have available.
        If keep_charges == False it will also delete the charges of self.
        If keep_identifier == False all atoms will have the identifier of self.
        If keep_identifier == True every atom will have its specific identifier,
        and self.identifier will be as many spacing as the length of the different identifiers.
        Identifiers of different length will raise an error.
        
        Parameters
        ----------
        g : geom
            geometry to add
        keep_charges : bool
            whether to keep charges from the 2 geometries and combine them
        keep_identifier: bool
            whether to keep the identifier of the new geometry for its atoms or not.
        """
        self.check_same_unit(g)
        if not keep_identifier:
            if (self.identifier==self.identifier[0]).all():
                idf = self.identifier[0]
            else:
                idf = ""
            g = change_identifier(g, idf)
        self.atoms = np.append(self.atoms,g.atoms)
        self.identifier = np.append(self.identifier,g.identifier)
        self.coords = np.append(self.coords,g.coords,axis=0)
        todel=[]
        if keep_charges:
            for m in self.charges.keys():
                if m in g.charges.keys():
                    self.charges[m]=np.append(self.charges[m],g.charges[m])
                else:
                    todel.append(m)
            for m in todel:
                del self.charges[m]
        else:
            self.charges={}
    
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
        self.charge_dipoles[method]=calculate_charge_dipole(self,method,coords=self.coord_unit,out=dip_unit,charge=charge)
        self.charge_dipoles_units[method]=dip_unit
        return self.charge_dipoles[method]
    
    def change_coord_unit(self,out):
        """
        Note
        ----
        Changes the unit of coords
        
        Parameters
        ----------
        out : str
            desired unit for the output. Handles cases properly and au/a.u./bohr
        """
        dict_={"au":"au", "a.u.":"au", "bohr":"au", "angstrom":"angstrom"}
        if dict_[self.coord_unit.lower()]==dict_[out.lower()]:
            print("The coordinates are already in "+out)
        elif dict_[self.coord_unit.lower()]=="angstrom" and dict_[out.lower()]=="au":
            print("changing coordinates from "+self.coord_unit+" to "+out)
            self.coords=np.multiply(1.88973,self.coords)
            self.coord_unit="au"
        elif dict_[self.coord_unit.lower()]=="au" and dict_[out.lower()]=="angstrom":
            print("changing coordinates from "+self.coord_unit+" to "+out)
            self.coords=np.multiply(0.529177,self.coords)
            self.coord_unit="angstrom"
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
        return geom(self.atoms[o],self.coords[o])
    
    def to_xyz(self, fnm, decimals=6, spacing=4):
        """
        Parameters
        ----------
        fnm: str
            filename to write
        decimals: int
            number of decimal digits desired
        spacing: int
            number of spaces desired
        """
        from geomtools.geom_io import write_xyz
        write_xyz(self, fnm, decimals=decimals, spacing=spacing)
        
    def to_coordfile(self, fnm, unit="Angstrom", decimals=6, spacing=4):
        """
        Parameters
        ----------
        fnm: str
            filename to write
        decimals: int
            number of decimal digits desired
        spacing: int
            number of spaces desired
        """
        from geomtools.geom_io import write_coords
        write_coords(self, fnm, unit=unit, decimals=decimals, spacing=spacing)
        
    def to_string(self, identifier = True, decimals=6, spacing=4):
        """
        Parameters
        ----------
        identifier: bool
            whether to print the identifier
        decimals: int
            number of decimal digits desired
        spacing: int
            number of spaces desired
        """
        return self.__str__(identifier=identifier, spacing=spacing, decimals=decimals)

def ghostify(g):
    """
    Note
    ----
    Center of mass ignores ghost atoms (they have no mass)
    
    Parameters
    ----------
    g
        geometry to ghostify
    
    Returns
    -------
    geom 
        ghostified geometry
    """
    c = g.copy()
    c.ghostify()
    return c

def strip_identifier(g):
    """
    Note
    ----
    Strips the identifier.
    
    Parameters
    ----------
    g: geom 
        geometry to strip the identifier of
    """
    c = g.copy()
    c.strip_identifier()
    return c

def add_identifier(g, identifier):
    """
    Note
    ----
    Adds an identifier.
    
    Parameters
    ----------
    g: geom 
        geometry to add the identifier to
    identifier: str 
        desired identifier. if not string tries to convert
    """
    c= g.copy()
    c.add_identifier(identifier)
    return c

def change_identifier(g, identifier):
    """
    Note
    ----
    Changes the identifier.
    
    Parameters
    ----------
    g: geom 
        geometry to change the identifier of
    identifier: str 
        desired identifier. if not string tries to convert
    """
    c= g.copy()
    c.change_identifier(identifier)
    return c        

def have_same_atoms(g1,g2):
    """
    Note
    ----
    Useful for scripts where different operations are to be done if geomA.atoms==geomB.atoms or geomA.atoms!=geomB.atoms
    
    Parameters
    ----------
    g1: geom
        the geometry to compare
    g2: geom
        the geometry to compare to
        
    Returns
    -------
    bool
        whether they have the same atoms or not
    """
    if len(g1.atoms) != len(g2.atoms):#Check if geometries have the same atom list
        return False
    elif (g1.atoms!=g2.atoms).all(): #Check if geometries have the same atom list
        return False
    else:
        return True
    
def check_same_atoms(g1,g2):
    """
    Note
    ----
    Returns nothing! Just raises specific errors in case
    
    Parameters
    ----------
    g1: geom
        the geometry to compare
    g2: geom
        the geometry to compare to
        
    """
    if len(g1.atoms) != len(g2.atoms):#Check if geometries have the same atom list
        raise NatomsError()
    elif (g1.atoms!=g2.atoms).all(): #Check if geometries have the same atom list
        raise atomsError()
        
def have_same_unit(g1,g2):
    """
    Note
    ----
    Useful for scripts where different operations are to be done if geomA.coord_unit==geomB.coord_unit or geomA.coord_unit!=geomB.coord_unit (e.g. convert)
    
    Parameters
    ----------
    g1: geom
        the geometry to compare
    g2: geom
        the geometry to compare to
        
    Returns
    -------
    bool
        whether they have the same unit or not
    """
    if g1.coord_unit != g2.coord_unit:#Check if geometries have the same atom list
        return False
    else:
        return True
    
def check_same_unit(g1,g2):
    """
    Note
    ----
    Returns nothing! Just raises specific errors in case
    
    Parameters
    ----------
    g1: geom
        the geometry to compare
    g2: geom
        the geometry to compare to

    """
    if g1.coord_unit != g2.coord_unit:#Check if geometries have the same atom list
        raise unitError() 
        
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
        coc=np.divide(np.sum(g.coords,axis=0),np.multiply(charge,len(g.atoms)))
    else:
        coc=np.zeros(3)
    d=np.sum(np.multiply(g.charges[method].reshape(len(g.atoms),1),g.coords-coc),axis=0)
    dict_={"au":"au", "a.u.":"au", "bohr":"au", "angstrom":"angstrom", "debye":"debye"}
    if coords.lower() not in dict_.keys() or out not in dict_.keys():
        print("combination of units of measure not implemented yet. Why don't you do it, champ?")
        raise NotImplementedError
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
