'''
Module: ediff.gcryst
--------------------

Collection of algorithms from general/geometric crystallography.
'''

import numpy as np
from pymatgen.core import Lattice
from math import gcd
from functools import reduce

class CrystalCell(Lattice):
    '''
    Class defining CrystalCell object, which contains:
    
    * unit cell parameters
    * unit cell-related functions from geometrical crystallography 
    
    Simple usage of the class:
        
    * The following example shows key functions of CrystalCell objects.
    * Only if you need more, continue reading the *Technical details* section.
    
    >>> # CrystalCell class :: simple usage
    >>> import ediff as ed
    >>>
    >>> # Define cubic unit cell and calculate interplanar distance
    >>> cell = ed.gcryst.CrystalCell.cubic(a=4.08)
    >>> print(cell.dhkl([0,2,0]))  # 2.04
    >>>
    >>> # Define hexagonal unit cell and calculate:
    >>> # 1) normal vector [uvw] of lattice plane (hkl)
    >>> # 2) lattice plane (hkl) of normal vector [uvw]
    >>> cell = ed.gcryst.CrystalCell.hexagonal(a=5.91, c=3.49)
    >>> print(cell.normal_of_plane([1,0,0]))  # (2,1,0)
    >>> print(cell.plane_of_normal([2,1,0]))  # (1,0,0)
    
    Technical details:
    
    * CrystalCell is a subclass of pymatgen.core.lattice.Lattice
        - It contains all props and methods of pymatgen.core.lattice.Lattice
        - https://pymatgen.org/pymatgen.core.html#module-pymatgen.core.lattice
    * CrystalCell class offers the following features:
        - Initialization with unit cell definition = using Lattice methods
        - All properties and methods from Lattice class are fully inherited
        - Several user-friendly function for Xtallographic calculations added
    '''

    def dhkl(self, hkl):
        '''
        Compute interplanar distance for (hkl) plane.

        Parameters
        ----------
        hkl : array-like
            Miller indices of the plane

        Returns
        -------
        dhkl : float
            Interplanar distance for the (hkl) plane.
        '''
        return( self.d_hkl(hkl) )

    
    @staticmethod
    def _integerize(vec, tol=1e-6):
        """
        Round a fractional vector to small integers for (hkl) or [uvw].

        Parameters
        ----------
        vec : array-like
            Fractional vector
        tol : float
            Small values below tol are treated as zero

        Returns
        -------
        tuple of ints
        """
        v = np.array(vec, dtype=float)
        v[np.abs(v) < tol] = 0.0
        if np.all(np.abs(v) < tol):
            return tuple(int(0) for _ in v)
        # Scale to smallest nonzero component = 1
        scale = 1.0 / np.min([abs(x) for x in v if abs(x) > tol])
        int_vec = tuple(int(round(x * scale)) for x in v)
        # Reduce by GCD of nonzero elements
        nz = np.array([abs(x) for x in int_vec if x != 0], dtype=int)
        if len(nz) > 0:
            g = np.gcd.reduce(nz)
            int_vec = tuple(int(x // g) for x in int_vec)
        return int_vec


    def normal_of_plane(self, hkl, tol=1e-6):
        """
        Compute the [uvw] direction, which is normal to a given (hkl) plane.

        Parameters
        ----------
        hkl : array-like
            Miller indices of the plane
        tol : float
            Tolerance for rounding small values

        Returns
        -------
        tuple of ints
            Miller indices [uvw] of the direction,
            which normal to (hkl) plane.
        """
        # Cartesian vector of plane normal (reciprocal lattice)
        cart_normal = self. \
            reciprocal_lattice_crystallographic.get_cartesian_coords(hkl)
        # Fractional coordinates in direct lattice
        uvw_frac = self.get_fractional_coords(cart_normal)
        # Round to small integers
        return self._integerize(uvw_frac, tol=tol)


    def plane_of_normal(self, uvw, tol=1e-6):
        """
        Compute the (hkl) plane, which is normal to a given [uvw] direction.

        Parameters
        ----------
        uvw : tuple of ints
            Direction indices [uvw]
        tol : float
            Tolerance for rounding small values

        Returns
        -------
        tuple of ints
            Miller indices (hkl) of the plane,
            which is normal to [uvw] direction.
        """
        # Cartesian vector of the direction
        cart_dir = self.get_cartesian_coords(uvw)
        # Fractional coordinates in reciprocal lattice
        hkl_frac = self. \
            reciprocal_lattice_crystallographic.get_fractional_coords(cart_dir)
        # Round to small integers
        return self._integerize(hkl_frac, tol=tol)


class HexLattice:
    '''
    Class with functions converting *Miller* - *Bravais* - *Weber* indices.
    
    >>> # HexLattice class :: simple usage
    >>> import ediff as ed
    >>>
    >>> # Conversions between Miller (hkl) and Bravais (hkil) indices
    >>> hkil = ed.gcryst.HexLattice.hkl_to_hkil([0,1,0])
    >>> hkl  = ed.gcryst.HexLattice.hkil_to_hkl([0,1,-1,0])
    >>> print(f'{hkl} <-> {hkil}')  # (0, 1, 0) <-> (0, 1, -1, 0)
    >>>
    >>> # Conversions between Miller [uvw] and Weber [uvtw] indices
    >>> uvtw = ed.gcryst.HexLattice.uvw_to_uvtw([2,1,0])
    >>> uvw  = ed.gcryst.HexLattice.uvtw_to_uvw([1,0,-1,0])
    >>> print(f'{uvw} <-> {uvtw}')  # (2, 1, 0) <-> (1, 0, -1, 0)
    '''

    def _convert_to_smallest_integers(uvtw):
        # Functional programming style
        # (why functional programming? => finding gcd of more numbers 
        # (source/inspiration: https://stackoverflow.com/q/16628088/
        # (1) map: for each element - multiply by 3 and round
        uvtw = list(map(lambda x: round(x*3), uvtw))
        # (2) reduce: for all elements - calculate greatest common divisor
        my_gcd = reduce(gcd, uvtw)
        # (3) map: for each element - divide by gcd and round
        uvtw = list(map(lambda x: round(x/my_gcd) , uvtw))
        # Return the modified list
        return(uvtw)

    
    def hkl_to_hkil(hkl):
        '''
        Convert hexagonal indices: Miller (hkl) to Bravais (hkil).
    
        Parameters
        ----------
        hkl : tuple or list or array of three integers
            Miller indices (hkl)
            defining a crystallographic plane.
    
        Returns
        -------
        hkil : four integers
            Bravais indices (hkil) of the crystallographic plane.
        '''
        h,k,l = hkl
        i = -(h+k)
        return(h,k,i,l)
    
    
    def hkil_to_hkl(hkil):
        '''
        Convert hexagonal indices: Bravais (hkil) to Miller (hkl).

        Parameters
        ----------
        hkil : tuple or list or array of four integers
            Bravais indices (hkil)
            defining a crystallographic plane.

        Returns
        -------
        hkl : tuple of three integers
            Miller indices (hkl) of the crystallographic plane.
        '''
        h,k,i,l = hkil
        return(h,k,l)
    
    
    def uvw_to_uvtw(uvw):
        '''
        Convert hexagonal indices: Miller [uvw] to Weber [uvtw]

        Parameters
        ----------
        uvw : tuple or list or array of three integers
            Miller indices [uvw]
            defining a crystallographic direction.

        Returns
        -------
        uvtw : tuple of four integers
            Weber indices [uvtw] of the crystallographic direction.
        '''
        u,v,w = uvw
        U = (2*u - v) / 3
        V = (2*v - u) / 3
        T = - (u + v) / 3
        W = w
        U,V,T,W = HexLattice._convert_to_smallest_integers([U,V,T,W])
        return(U,V,T,W)

    
    def uvtw_to_uvw(uvtw):
        '''
        Convert hexagonal indices: Weber [uvtw] to Miller [uvw].

        Parameters
        ----------
        uvtw : tuple or list or array of four integers
            Weber indices [uvtw]
            defining a crystallographic direction.

        Returns
        -------
        uvw : tuple of three integers
            Miller indices [uvw] of the crystallographic direction.
        '''
        U,V,T,W = uvtw
        u = U - T
        v = V - T
        w = W
        u,v,w = HexLattice._convert_to_smallest_integers([u,v,w])
        return(u,v,w)
        

class ZoneAxis:
    '''
    Class with functions to determine *zone axis* of a 2D diffractogram.
    '''

    def from_two_diffs(diff1, diff2):
        '''
        Calculate zone axis [uvw] from (hkl) of 2 diffractions in diff.pattern.

        Parameters
        ----------
        diff1 : tuple or list or array of three integers
            Miller indices (hkl) of the 1st diffraction within given zone 
        diff2 : TYPE
            Miller indices (hkl) of the 2nd diffraction within given zone.

        Returns
        -------
        zone_axis : tuple of three integers
            Miller indices [uvw] of the zone axis
            calculated from *diff1* and *diff2*.
        
        Technical notes
        ---------------
        * Zone axis should be perpendicular
          to any two diffractions of given zone.
        * Therefore, it can be calculated as a vector cross-product
          of any two diffractoins of given zone.
        * We use np.cross to calculate the cross product
          and convert the result to tuple for the sake of consistency
          as all other functions in this module return vectors as tuples.
        '''
        zone_axis = np.cross(diff1,diff2)
        return( tuple(zone_axis) )


    def from_multiple_diffs(diffs):
        '''
        Calculate zone axis [uvw] from (hkl) of N diffractions in diff.pattern.

        Parameters
        ----------
        diffs : list of 3-item lists = list of (hkl) of multiple diffractions
            The list containing (hkl) indices of multiple diffractions
            that are supposed to correspond to one zone axis [uvw].

        Returns
        -------
        None
            The list of vector cross-products is printed on the screen.
            If all cross-products are the same (up to scale/sign),
            then the zone axis [uvw] explains the presense of all *diffs*
            in the observed diffraction pattern.
        '''
        for d1 in (diffs):
            for d2  in (diffs):
                index_d1 = diffs.index(d1) + 1
                index_d2 = diffs.index(d2) + 1
                if index_d1 < index_d2:
                    print(f'd{index_d1} x d{index_d2}', end='')
                    print(f' = {np.cross(d1,d2)}')

    def weiss_zone_law(zone_axis):
        '''
        Using Wess Zone Law, calculate diffractions (hkl) for zone axis [uvw]. 

        Parameters
        ----------
        zone_axis : tuple or list or array of three integers
            Crystallographic direction [uvw] that defines zone axis.

        Returns
        -------
        WZL : string
            A condition for (hkl) indices of diffractions
            that satisfy Weis Zone Law for given zone axis/direction [uvw].
        
        Technical note
        --------------
        * The condition is a simple string.
        * It would be more straightforward to use a SymPy expression.
        * Nevertheless, in this simple case the "manual symbolic calculation"
          is quite easy and so we do not need the heavy-weight SymPy library.
        '''
        # WZL = Weiss Zone Law = h*u + k*v + l*w = 0
        WZL = ''
        u,v,w = zone_axis
        # Process h*u part of WZL
        if   u == 1     : WZL += 'h'
        elif u == -1    : WZL += '-h'
        elif abs(u) > 1 : WZL += str(u) + 'h'
        # Process k*v part of WZL
        if   v == 1     : WZL += ' + k'
        elif v == -1    : WZL += ' - k'
        elif v > 1      : WZL += ' + ' + str(v) + 'k'
        elif v < -1     : WZL += ' - ' + str(v) + 'k'
        # Process l*w part of WZL
        if   w == 1     : WZL += ' + l'
        elif w == -1    : WZL += ' - l'
        elif w > 1      : WZL += ' + ' + str(w) + 'l'
        elif w < -1     : WZL += ' - ' + str(w) + 'l'
        # Add right-hand side of WZL
        WZL += ' = 0'
        # Return the result
        return(WZL)
        
    