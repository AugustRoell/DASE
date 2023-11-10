"""
Helper functions to define standard dome profiles
"""

# import packages
from numpy import (
    cos,
    exp,
    pi,
    sqrt,
)



def sphere(
    r : float,
    R : float = 11.521e-3,
    h : float = 0.434e-6,
): 
    """Returns the profile of a plano-convex spherical cap HBAR.
    
    Parameters:
    -----------
    r: float
        r-coordinate (in meters)
    R: float
        Radius of curvature (in meters)
    h: float
        Height at apex of dome (in meters)
    """

    if r**2 <= 2*R*h - h**2:
        return(h - R + sqrt(R**2 - r**2))
    return 0.


def SBQ(
        r : float,
        R : float = 100.0,
        h : float = 0.434e-6,
        a : float = 0.061375350123207426,
):
    """Returns the profile of a plano-convex SBQ HBAR.
    
    Parameters:
    -----------
    r: float
        r-coordinate (in meters)
    R: float
        Radius of curvature at top of dome (in meters)
    h: float
        Height at apex of dome (in meters)
    a : float
        Positive anharmonicity parameter (in meters^-2)
    """

    if r**2 < h*(sqrt(4*a*R**2 + 1) - 1)/(2*a*R):
        return(sqrt(h**2 - (h/R)*(r**2) - a*(r**2)**2))
    return 0.


def pnorm(
        r : float,
        R : float = 0.002695958723694581,
        h : float = 0.434e-6,
        p : int = 6,
):
    """Returns the profile of a plano-convex p-norm HBAR.
    
    Parameters:
    -----------
    r: float
        r-coordinate (in meters)
    R: float
        Radius of curvature parameter. For p=2 it is the sphere radius, for 
        p->inf it approaches the cylinder radius (in meters)
    h: float
        Height at apex of dome (in meters)
    p : int
        P-value of the p-norm. p=2 gives the spherical cap, p-> converges to 
        the cyliner (unitless)
    """

    if r < (R**p - (h - R)**p)**(1/p):
        return(h - R + (R**p - r**p)**(1/p))
    return 0.


def cylinder(
        r : float,
        h : float = 0.434e-6,
        r_cyl : float = 337.22e-6,
):
    """Returns the profile of a plano-convex cylindrical HBAR.
    
    Parameters:
    -----------
    r: float
        r-coordinate (in meters)
    h: float
        Height at apex of dome (in meters)
    r_cyl : float
        Cylinder radius (in meters)
    """

    if r <= r_cyl:
        return(h)
    return(0.)


def gaussian(
    r : float,
    s : float = 7.07e-5,
    h : float = 0.434e-6,
):
    """Returns the profile of a plano-convex gaussian HBAR.
    
    Parameters:
    -----------
    r: float
        r-coordinate (in meters)
    s: float
        Standard deviation  (in meters)
    h: float
        Height at apex of dome (in meters)
    """

    return(h*exp(-(r**2)/(2*s**2)))


def cosine(
    r : float, 
    p : float = 10e-6, 
    h : float = 0.434e-6,
    ):
    """Returns the profile of a plano-convex cosine HBAR.
    
    Parameters:
    -----------
    r: float
        r-coordinate (in meters)
    p: float
        Period  (in meters)
    h: float
        Height at apex of dome (in meters)
    """

    if r < p/2.:
        return((h/2)*(1+cos(2*pi*r/p)))
    return 0.0


def sphere_clipped(
        r : float,
        R : float = 20.0e-3,
        h : float = 0.434e-6,
        r_clip : float = 50e-6,
):
    """Returns the profile of a plano-convex a discontinuously clipped 
    spherical cap HBAR.
    
    Parameters:
    -----------
    r: float
        r-coordinate (in meters)
    R: float
        Radius of curvature (in meters)
    h: float
        Height at apex of dome (in meters)
    r_clip : float
        r-coordinate at which the profile discontinuously 
        goes to zero (in meters)
    """

    if r_clip > sqrt(2*R*h-h**2):
        raise ValueError("The clip radius r_clip should be smaller than or "
        "equal to the cutoff radius of the dome.\nThe values of the clip and "
        "cutoff radii respectively were: "
        f"{r_clip*1e6:9.3f} um, {sqrt(2*R*h-h**2)*1e6:9.3f} um")
    if r <= clip_radius:
        return(h - R + sqrt(R**2 - r**2))
    return(0.)
