'''
Module: ediff.gcryst
--------------------

Collection of formulas from general/geometric crystallography.
'''

from pymatgen.core import Lattice
import numpy as np

class CrystalCell(Lattice):
    """
    Class defining CrystalCell object, which contains:
    
    * unit cell parameters
    * unit-cell related functions from geometrical crystallography. 
    
    Recommendation:
        
    * Start with examples within this documentation.
    * Only if you need more, continue reading the *Technical details* section.
    
    Technical details:
    
    * CrystalCell is a subclass of pymatgen.core.lattice.Lattice
    * it contains all properties and methods of pymatgen.core.lattice.Lattice
    * doc: https://pymatgen.org/pymatgen.core.html#module-pymatgen.core.lattice
    * CrystalCell unit cell can be defined/initialized with Lattice methods
    * CrystalCell properties and methods from Lattice class are fully inherited
    * CrystalCell adds user-friendly extentsions for Xtallographic calculations
    """
    
    @staticmethod
    def _integerize(vec, tol=1e-6):
        """
        Round a fractional vector to small integers for [uvw] or (hkl).

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
        Compute the [uvw] direction normal to a given (hkl) plane.

        Parameters
        ----------
        hkl : tuple of ints
            Miller indices of the plane
        tol : float
            Tolerance for rounding small values

        Returns
        -------
        tuple of ints
            [uvw] direction normal to the plane
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
        Compute the (hkl) plane normal to a given [uvw] direction.

        Parameters
        ----------
        uvw : tuple of ints
            Direction indices [uvw]
        tol : float
            Tolerance for rounding small values

        Returns
        -------
        tuple of ints
            Miller indices (hkl) of the plane normal to the direction
        """
        # Cartesian vector of the direction
        cart_dir = self.get_cartesian_coords(uvw)
        # Fractional coordinates in reciprocal lattice
        hkl_frac = self. \
            reciprocal_lattice_crystallographic.get_fractional_coords(cart_dir)
        # Round to small integers
        return self._integerize(hkl_frac, tol=tol)
