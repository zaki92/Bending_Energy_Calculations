# -*- coding: utf-8 -*-
"""
These functions can be used to calculate the surface areas and bending energies 
of oblate spheroids, elliptic hyperboloids and combinations of oblate spheroids
and elliptic hyperboloids, as described in the publication 'A Biophysical 
Model for Plant Cell Plate Development' by Jawaid et al. Further explanations  
about how to use this codeare given in the functions.'

Authored by Muhammad Zaki Jawaid
Contact via email (zjawaid@ucdavis.edu) for questions
Uploaded to GitHub 8/23/2020
"""

import numpy as np
from scipy import integrate

def surface_area_hyperboloid_generic(a,b,c,h):
    
    """
    Input:
        a (float): major axis radius at the center of the 
                   hyperboloid
        b (float): minor axis radius at the center of the 
                   hyperboloid
        c (float): measure of elongation of the hyperboloid along the
                   polar axis
        h (float): height of the hyperboloid

    Output: 
        Surface area of a elliptical one sheeted hyperboloid with
        given a,b,c,h (float)
    """
    
    #simplifyng variables
    e = (b**2)*(c**2)
    f = (a**2)*(c**2)
    g = (a**2)*(b**2)
    
    #u parameterizes the surface along the z-axis
    #v goes from 0 to 2pi
    
    integrand = lambda v,u : np.sqrt(e*(u**2+1)*(np.cos(v)**2)+f*(u**2+1)*
                                     (np.sin(v)**2)+g*u*u)
    
    i = lambda v : 0.0
    j = lambda v : np.pi/2
    
    return (8*integrate.dblquad(integrand,0,h/(2*c),i,j)[0])

def get_pd(a1,c1,a2,c2):
    
    """
    This function is used to calculate pd. 
    
    pd is described as the penetration distance, or the shared distance of the 
    hyperboloid and the oblate spheroid, and is used to calculate the height 
    of the connected hyperboloid.

    Input:
        a1 (float) : major axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        c1 (float) : minor axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        a2 (float) : major axis radius of the hyperboloid at the center of the
                     hyperboloid.
        c2 (float) : minor axis radius of the hyperboloid at the center of the
                     hyperboloid.   

    Output: 
        pd (float) : penetration distance, or the shared distance of the 
                     hyperboloid and the oblate spheroid, and is used to 
                     calculate the height of the connected hyperboloid.
    """
    
    d = get_d(a1,c1,a2,c2)
    pd = a1-d*(a2**2/(a2**2+c2**2))
    return pd

def get_h(a1,c1,a2,c2):

    """
    This function is used to calculate the height of a hyperboloid the is 
    connected to an oblate spheroid. 
    
    Input:
        a1 (float) : major axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        c1 (float) : minor axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        a2 (float) : major axis radius of the hyperboloid at the center of the
                     hyperboloid.
        c2 (float) : minor axis radius of the hyperboloid at the center of the
                     hyperboloid.   

    Output: 
        h (float) : height of a hyperboloid connected to an oblate spheroid.
    """
    
    d = get_d(a1,c1,a2,c2)
    h = 2*d*(1-a2**2/(a2**2+c2**2))
    return h


def surface_area_hyperboloid(a1,c1,a2,c2):
    
    """
    This functions returns the surface area of a one-sheeted elliptical 
    hyperboloid that is connected to an oblate spheroid such that curvature 
    continuity is maintained. Because of this constraint, Eqs. S2-S4 from the 
    supplemental information can be invoked. 
    
    Input: 
        a1 (float) : major axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        c1 (float) : minor axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        a2 (float) : major axis radius of the hyperboloid at the center of the
                     hyperboloid.
        c2 (float) : minor axis radius of the hyperboloid at the center of the
                     hyperboloid.   
                     
    Output: 
        Surface area of a one-sheeted elliptical hyperboloid 
        that is connected to an oblate spheroid (float).
                     
    """
    
    # the height, h, and the minor axis radius, b, of the hyperboloid is 
    # constrained by the oblate spheroid it is connected to because of 
    # continuity considerations, as explained in Eqs. S2 - S4.
    # the get_h function calculates this h, and the expression for b is given below.
    
    h = get_h(a1,c1,a2,c2)
    b = a2*c1/a1  
    
    #simplifying variables
    a = a2
    c = c2
    
    e = (b**2)*(c**2)
    f = (a**2)*(c**2)
    g = (a**2)*(b**2)
    
    # u,v, parameterize the surface
    #u is along z axis
    #v goes from 0 to 2pi
    
    integrand = lambda v,u : np.sqrt(e*(u**2+1)*(np.cos(v)**2)+f*(u**2+1)*
                                    (np.sin(v)**2)+g*u*u)
    
    i = lambda v : 0.0
    j = lambda v : np.pi/2
    
    return (8*integrate.dblquad(integrand,0,h/(2*c),i,j)[0])

def get_d(a1,c1,a2,c2):

    """
    This function is used to calculate d. 
    
    d is described in Eq. S1 as the distance from the the center point of the 
    hyperboloid on the polar axis from the center of the connected oblate 
    spheroid, and is used to calculate the penetration depth 
    as well as the height of the connected hyperboloid.

    Input:
        a1 (float) : major axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        c1 (float) : minor axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        a2 (float) : major axis radius of the hyperboloid at the center of the
                     hyperboloid.
        c2 (float) : minor axis radius of the hyperboloid at the center of the
                     hyperboloid.   

    Output: 
        d (float) : distance from the the center point of the 
                    hyperboloid on the polar axis from the center of the 
                    connected oblate spheroid.
    """
    
    d = (np.sqrt(a1**2-a2**2))*(np.sqrt(a2**2+c2**2))/a2
    return d

def surface_area_os(a,c):

    """
    This function returns the surface area of an oblate spheroid, 
    given that a >= c.
    
    Input: 
        a (float) : major axis of the oblate spheroid
        c (float) : minor axis of the oblate spheroid
        
    Output:
        S (float) : surface area of the oblate spheroid
    """
    
    if a == c:
        return 4*np.pi*a*c
    else:
        ecc = np.sqrt(1.0-((c**2)/(a**2)))
        S = 2*np.pi*(a**2)+(np.log((1+ecc)/(1-ecc)))*np.pi*(c**2)/(ecc)
        return S


def surface_area_correction_os(a1,c1,a2,c2):
    
    """
    This function calculates the area of within the intersection created
    when an oblate spheroid and a hyperboloid are connected with the 
    constraint of curvature continuity.
    
    Input:
        a1 (float) : major axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        c1 (float) : minor axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        a2 (float) : major axis radius of the hyperboloid at the center of the
                     hyperboloid.
        c2 (float) : minor axis radius of the hyperboloid at the center of the
                     hyperboloid.   

    Output: 
        area within intersection (float)
    """
    
    #center offset distance as calculated in get_d
    
    d = get_d(a1,c1,a2,c2) 
    
    #simplifying variables
    
    cons = 1-(((d*(a2**2))/(a2**2+c2**2))**2)/(a1**2)
    
    integrand = lambda u,v : (a1/np.sqrt(2))*np.sqrt(a1**2+c1**2+(a1**2-c1**2)*
                                                     np.cos(2*v))*np.sin(v)
    i = lambda u : 0
    j = lambda v : np.arcsin(np.sqrt((cons-np.cos(v)**2)/(np.sin(v)**2)))
    v_i = np.arccos(np.sqrt(cons))
    
    return (4*integrate.dblquad(integrand,v_i,np.pi/2,i,j)[0])

def bending_energy_os(a,c,c0,K):
    
    """
    This functions calculates the mean curvature energy of an oblate spheriod
    with major axis radius 'a', and minor axis radius 'c', given that a >= c. 
    
    Input: 
        a (float) : major axis of the oblate spheroid
        c (float) : minor axis of the oblate spheroid
        c0 (float): spontaneous curvature as described in Eq. 2
        K (float) : bending rigidity, as described in Eq. 2
        
    Output:
        mean curvature energy of oblate spheroid (float)
        
    """
    
    E = lambda x: ((K/2)*((np.sqrt(2)*c*(3*a**2+c**2+(a**2-c**2)*np.cos(2*x))/
                  (a*(a**2+c**2+(a**2-c**2)*np.cos(2*x))**1.5) -2*c0)**2)*
                  np.sqrt(2)*np.pi*a*(np.sqrt(a**2+c**2+(a**2-c**2)*np.cos(2*x)))
                  *np.sin(x))
    
    return integrate.quad(E,0,np.pi)[0]

def bending_energy_hyperboloid(a1,c1,a2,c2,c0,K):
    
    """
    This function calculates the mean curvature energy of an elliptical hyperboloid
    that is connected to an oblate spheroid with curvature continuity. 
    
    Because of this constraint, Eqs. S2-S4 from the 
    supplemental information can be invoked. 
    
    Input: 
        a1 (float) : major axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        c1 (float) : minor axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        a2 (float) : major axis radius of the hyperboloid at the center of the
                     hyperboloid.
        c2 (float) : minor axis radius of the hyperboloid at the center of the
                     hyperboloid. 
        
        c0 (float): spontaneous curvature as described in Eq. 2
        
        K (float) : bending rigidity, as described in Eq. 2
                     
    Output: 
        bending energy of the hyperboloid (float)

    """
    
    #calculating the height and the minor axis radius
    h = get_h(a1,c1,a2,c2)
    b = a2*c1/a1
    
    
    #simplifying variables
    
    a = a2
    c = c2
    e = (b**2)*(c**2)
    f = (a**2)*(c**2)
    g = (a**2)*(b**2)
    
    integrand = lambda v,u : (((((a*b*c)*((b*b*u*u-a*a)*(np.sin(v)**2)+
                             (a*a*u*u-b*b)*(np.cos(v)**2)+
                             (u*u+1)*(c*c))/((e*(u**2+1)*
                             (np.cos(v)**2)+f*(u**2+1)*(np.sin(v)**2)+g*u*u)**1.5))-2*c0)**2)*
                             (np.sqrt(e*(u**2+1)*(np.cos(v)**2)+
                             f*(u**2+1)*(np.sin(v)**2)+g*u*u)))
                                              
    i = lambda v : 0.0
    j = lambda v : np.pi/2
    
    return (0.5*K*(8*integrate.dblquad(integrand,0,h/(2*c),i,j)[0]))

def bending_energy_os_correction(a1,c1,a2,c2,c0,K):
    
    """
    This function calculates the contribution to the mean curvature energy by
    the surface area within the intersection that is created
    when an oblate spheroid and a hyperboloid are connected with the 
    constraint of curvature continuity.
    
    Input:
        a1 (float) : major axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        c1 (float) : minor axis radius of the oblate spheroid that the 
                     hyperboloid is connected to.
        a2 (float) : major axis radius of the hyperboloid at the center of the
                     hyperboloid.
        c2 (float) : minor axis radius of the hyperboloid at the center of the
                     hyperboloid.   
        c0 (float): spontaneous curvature as described in Eq. 2

        K (float) : bending rigidity, as described in Eq. 2

    Output: 
        mean curvature energy cortribution due to
        the area within intersection (float)
    
    """
    #calculating offset distance
    d=get_d(a1,c1,a2,c2)
    
    #x is the theta angle (in the xy plane)
    
    #we do the phi integral as a function of theta, then we do the theta integral
    integrand = lambda x,u: (((K/2)*((np.sqrt(2)*c1*
                            (3*a1**2+c1**2+(a1**2-c1**2)*np.cos(2*x))/
                            (a1*(a1**2+c1**2+(a1**2-c1**2)*np.cos(2*x))**1.5) -2*c0)**2)*
                             np.sqrt(2)*np.pi*a1*(np.sqrt(a1**2+c1**2+(a1**2-c1**2)*np.cos(2*x)))*
                             np.sin(x))/(2*np.pi))
    
    cons=1-(((d*(a2**2))/(a2**2+c2**2))**2)/(a1**2)
    
    i = lambda u : 0
    j = lambda x : np.arcsin(np.sqrt((cons-np.cos(x)**2)/(np.sin(x)**2)))
    
    v_i = np.arccos(np.sqrt(cons))
    
    return (4*integrate.dblquad(integrand,v_i,np.pi/2,i,j)[0])







   



