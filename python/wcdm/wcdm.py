# Python implementation of formulae in Eisenstein (2015).
# 
# Written by Daniel Eisenstein, save for adaptation of GSL code in 
# beta_cont_frac_gsl(a, b, x).  

'''
wcdm(z,om,w):  Compute the proper motion distance in (c/H_0) units.
	Allows vector z, om, and/or w
	Can accept any om>0.
	Requires w<-1/3 if one is using the SciPy beta function, w<0 if using GSL

wcdm_time(z, om, w): Compute t(z) in (1/H_0) units.
	Allows vector z and w
    	Requires om<=1, but can accept a vector.
	Use om=0.99999 if you want om=1.
	Requires w<-1 if one is using the SciPy beta function, w<0 if using GSL

owcdm(z,om,w,ok): Compute the comoving coordinate distance and proper motion distance,
	both in (c/H_0) units.
	Allows z, om, or w as vectors; will produce a vector output.
	ok must be a scalar.
	Requires om<1 and ox=1-om-ok>0.  
	Use om=0.99999 if you want om=1.
	Requires w<-1/3 if one is using the SciPy beta function, w<0 if using GSL
 	The curvature is treated only perturbatively; it will be fine for |Ok|<0.1,
 	but will start to produce noticeable offsets for Ok ~ Om.

coorddist(z, om, w, ok): owcdm() driver for comoving coordinate distance
propmotdis(z, om, w, ok): owcdm() driver for proper motion distance
angdist(z, om, w, ok): owcdm() driver for angular diameter distance
lumdist(z, om, w, ok): owcdm() driver for luminosity distance

wcdm_rad(z, om, w, rad=0.0): Compute the proper motion distance in (c/H_0) units,
	including an approximation for radiation terms.
	Allows vector z, om, and/or w.
	Can accept 0<om<1.
        Requires w<-1/3 if one is using the SciPy beta function, w<0 if using GSL
	This uses the beta functions with arguments that may not be allowed in 
	all implementations, but appears to work for the GSL code.

test(): Produces a set of accuracy tests relative to numerical evaluation
	of the integrals being used.

-----------------

Eisenstein is distributing the code in this file under the BSL license.

However, note that the GSL code in beta_cont_frac_gsl(a, b, x) is
under the Gnu Public License.  So if you want to use that particular
function, supplied in the file incomplete_beta.py, then you need
to treat this as bound by GPL, which is stricter.

Copyright (c) 2015, Daniel Eisenstein
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
'''

import numpy as np
import scipy.special
from incomplete_beta import beta_cont_frac_gsl

############### Implementation of the Incomplete Beta Function ###########

def incomplete_beta(a,b,x1,x2):
    # Compute B_x2(a,b) - B_x1(a,b).
    # This is the integral from x1 to x2 of x^(a-1) (1-x)^(b-1) dx
    # Normally this requires a>0, b>0, and x1,x2 in the range [0..1]
    # However, the b>0 requirement is relaxed if one uses the modified GSL code.

    # To use the modified GSL code.
    # Here we avoid the pre-factor normalization with the complete beta function.
    # This allows one to use b<0 without problem.
    return beta_cont_frac_gsl(a,b,x2)-beta_cont_frac_gsl(a,b,x1)

    # If you don't want the modified GSL code, here's the standard SciPy call.
    # If you use SciPy, then you must have b>0, which limits these routines 
    # to w < -1/3.  This is because larger w have a divergent horizon in the 
    # future, which breaks the normalization used internally to the betainc() code.
    # return scipy.special.beta(a,b) * \
    #	    (scipy.special.betainc(a,b,x2)- scipy.special.betainc(a,b,x1))


############# Dark Energy models #########################################

def wcdm(z, om, w):
    # For redshifts z and Omega_m om and constant w, 
    # compute D_M(z)=r(z) in Hubble distance (c/H0) units for a flat Universe.
    # Can give z or w as vectors; will produce a vector output.
    # om must be a scalar; see owcdm() to allow a vector
    # This requires w<-1/3 if one is using the SciPy beta function, w<0 if using GSL
    ox = 1-om
    m = 1.0/(-6.0*w)
    # Return the EdS simple case or compute the incomplete beta function
    return np.where(ox==0, 2.0*(1-(1+z)**-0.5),
	2.0*m/np.sqrt(om) * np.abs(om/np.where(ox!=0,ox,1.0))**m * \
	incomplete_beta(m, 
	    np.where(ox>0, 0.5-m, 0.5),
	    np.where(ox>0, ox/(om*(1+z)**(-3.0*w)+ox), -ox/om*(1+z)**(3.0*w)),
	    np.where(ox>0, ox/(om+ox), -ox/om)))
    # This is handling all three ox branches at once.
    # Not the prettiest code, but it avoids floating point errors on the 
    # conditional branch executions in the vector case.


def wcdm_time(z, om, w):
    # For redshifts z and Omega_m om and constant w, for a flat Universe,
    # compute the age of the Universe at the given z in Hubble (1/H0) units. 
    # Can give z or w as vectors; will produce a vector output.
    # This requires om<=1 and ox=1-om>=0.  
    # This requires w<-1 if one is using the SciPy beta function, w<0 if using GSL
    ox = 1-om
    if (np.min(ox)<0):
	print "Can't evaluate negative dark energies"
	exit
    xz = ox/(om*(1+z)**(-3.0*w)+ox)
    m = 1.0/(-2.0*w)
    return np.where(ox==0, 2.0/3.0*np.power(1+z,-1.5), 
	2.0/3.0*m/np.sqrt(om) * np.power(om/np.where(ox==0,1.0,ox),m) * 
	incomplete_beta(m,0.5-m,0.0,xz))
    # Separately compute the flat and open cases, avoiding 
    # the divide-by-zero possibility on the second branch.


def owcdm(z, om, w, ok=0.0):
    # For redshifts z and Omega_m om and constant w and a small ok 
    # compute r(z) and D_M(z) in Hubble distance (c/H0) units for LCDM.
    # Can give z, om, or w as vectors; will produce a vector output.
    # ok must be a scalar, because of our if statements.
    # We must have om<1 and ox=1-om-ok>0.  Use om=0.99999 if you want om=1
    # This requires w<-1/3 if one is using the SciPy beta function, w<0 if using GSL
    ox = 1-om-ok
    xz = ox/(om*(1+z)**(-3.0*w)+ox)
    x0 = ox/(om+ox)
    # Make the flat limit
    m = 1.0/(-6.0*w)
    c = 2.0*m/np.sqrt(om)
    rz = c * (om/ox)**m * incomplete_beta(m,0.5-m,xz,x0)
    if (ok==0):
	propmotdis = rz
    else:
	# Now make the curvature corrections, up to 10th order in Ok/Om
	# Stop when the contributions are small
	for order in range(1,11):
	    m = (1.0+2*order)/(-6.0*w)
	    c *= -(2.0*order-1.0)/(2.0*order)*ok/om
	    delrz = c*(om/ox)**m * incomplete_beta(m,0.5+order-m,xz,x0)
	    rz += delrz
	    if (np.max(np.abs(delrz))<1e-7): break
	# Compute the proper motion distance
	if (ok>0): propmotdis = np.sinh(np.sqrt(ok)*rz)/np.sqrt(ok)
	else: propmotdis = np.sin(np.sqrt(-ok)*rz)/np.sqrt(-ok)
    # Now return the answer
    return rz,propmotdis

############# Drivers ####################################################

def coorddist(z, om, w, ok): return owcdm(z,om,w,ok)[0]

def propmotdis(z, om, w, ok): return owcdm(z,om,w,ok)[1]

def angdist(z, om, w, ok): return owcdm(z,om,w,ok)[1]/(1+z)

def lumdist(z, om, w, ok): return owcdm(z,om,w,ok)[1]*(1+z)

############# Radiation perturbations #########################################

def incomplete_beta_approx(a,b,x1,x2):
    # Compute B_x2(a,b) - B_x1(a,b).
    # This is the integral from x1 to x2 of x^(a-1) (1-x)^(b-1) dx
    # However, this is a crude approximation; only suitable for 
    # very small perturbations, such as for the radiation term.
    i1 = x1**a/a*(1.0-a*(b-1.0)/(a+1.0)*x1)
    i2 = x2**a/a*(1.0-a*(b-1.0)/(a+1.0)*x2)
    return i2-i1

def wcdm_rad(z, om, w, rad=0.0):
    # For redshifts z and Omega_m om and constant w and a small ok 
    # compute r(z) and D_M(z) in Hubble distance (c/H0) units for LCDM.
    # Can give z, om, or w as vectors; will produce a vector output.
    # We must have om<1 and ox=1-om>0.
    # This requires w<-1/3 if one is using the SciPy beta function, w<0 if using GSL
    ox = 1-om-rad
    xz = ox/(om*(1+z)**(-3.0*w)+ox)
    x0 = ox/(om+ox)
    # Make the radiation-free limit
    m = 1.0/(-6.0*w)
    c = 2.0*m/np.sqrt(om)
    rz = c * (om/ox)**m * incomplete_beta(m,0.5-m,xz,x0)
    # Now add in the radiation part
    m = 1.0/(6.0*w)
    c = m*rad/om/np.sqrt(om)
    delrz = c*(om/ox)**m * incomplete_beta(m,1.5-m,xz,x0)
    # This uses a B(a,b,x) implementation with a<0, which is not formally allowed
    # and may be fragile.  But the GSL version seems to behave well enough.
    # Could use the following instead.
    #delrz = c*(om/ox)**m * incomplete_beta_approx(m,1.5-m,xz,x0)
    print incomplete_beta(m,1.5-m,xz,x0), incomplete_beta_approx(m,1.5-m,xz,x0)
    rz += delrz
    return rz

################### Direct integration codes for testing ###################

import scipy.integrate

def wcdm_romberg(z, om, w):
    integrand = lambda z: 1.0/np.sqrt(om*(1.0+z)**3.0+(1.0-om)*(1.0+z)**(3.0*(1.0+w)))
    return scipy.integrate.romberg(integrand, 0.0, z)

def owcdm_romberg(z, om, w, ok=0.0):
    integrand = lambda z: 1.0/np.sqrt(om*(1.0+z)**3.0+ok*(1.0+z)**2.0+(1.0-om-ok)*(1.0+z)**(3.0*(1.0+w)))
    return scipy.integrate.romberg(integrand, 0.0, z)

def time_romberg(z, om, w):
    integrand = lambda a: np.sqrt(a)/np.sqrt(om+(1.0-om)*a**(-3.0*w))
    return scipy.integrate.romberg(integrand, 0.0, 1.0/(1.0+z), tol=1e-5)

def wcdmrad_romberg(z, om, w, rad):
    integrand = lambda z: 1.0/np.sqrt(rad*(1.0+z)**4.0+om*(1.0+z)**3.0+(1.0-om-rad)*(1.0+z)**(3.0*(1.0+w)))
    return scipy.integrate.romberg(integrand, 0.0, z)


def rad_test(xmin):
    #integrand = lambda x: x**(-1.0/6.0-1) * (1-x)**(5.0/3.0-1.0)
    #print scipy.integrate.romberg(integrand, xmin, 0.7, tol=1e-5, divmax=13)
    integrand = lambda y: -6.0*(1-y**-6.0)**(5.0/3.0-1.0)
    print scipy.integrate.romberg(integrand, xmin**(-1.0/6.0), 0.7**(-1.0/6.0), tol=1e-5, divmax=13)
    print incomplete_beta_nrcf(-1.0/6.0,5.0/3.0,xmin,0.7)

########################  The test driver ################################
# To run the tests, execute test()

def test():
    print "Testing the base wcdm code"
    x = wcdm(1.0,0.3,-1) 
    print "LCDM, om=0.3, z=1: r(z) = ", x, " Err = ",x-wcdm_romberg(1.0,0.3,-1) 

    x = wcdm(1.0,0.3,-0.2)
    print "WCDM, om=0.3, w=-0.2, z=1: r(z) = ", x, " Err = ",x-wcdm_romberg(1.0,0.3,-0.2) 

    x = wcdm(1.0,0.3,-0.4)
    print "WCDM, om=0.3, w=-0.4, z=1: r(z) = ", x, " Err = ",x-wcdm_romberg(1.0,0.3,-0.4) 

    x = wcdm(1.0,0.3,-1.4) 
    print "WCDM, om=0.3, w=-1.4, z=1: r(z) = ", x, " Err = ",x-wcdm_romberg(1.0,0.3,-1.4) 

    x = wcdm(1.0,1.3,-0.4) 
    print "WCDM, om=1.3, w=-0.4, z=1: r(z) = ", x, " Err = ",x-wcdm_romberg(1.0,1.3,-0.4) 

    print "\nTesting the owcdm code in the flat limit"
    x = owcdm(1.0,0.3,-1)[0]
    print "Non-flat LCDM, om=0.3, ok=0, z=1: r(z) = ", x, " Err = ",x-wcdm_romberg(1.0,0.3,-1) 

    x = owcdm(1.0,0.3,-0.4)[0]
    print "Non-flat WCDM, om=0.3, w=-0.4, ok=0, z=1: r(z) = ", x, " Err = ",x-wcdm_romberg(1.0,0.3,-0.4) 

    x = owcdm(1.0,0.3,-1.4)[0]
    print "Non-flat WCDM, om=0.3, w=-1.4, ok=0, z=1: r(z) = ", x, " Err = ",x-wcdm_romberg(1.0,0.3,-1.4) 

    print "\nSlightly open universes"
    x = owcdm(1.0,0.3,-1,0.05)[0]
    print "Non-flat LCDM, om=0.3, ok=0.05, z=1: r(z) = ", x, " Err = ",x-owcdm_romberg(1.0,0.3,-1,0.05) 

    x = owcdm(1.0,0.3,-0.2,0.05)[0]
    print "Non-flat WCDM, om=0.3, w=-0.2, ok=0.05, z=1: r(z) = ", x, " Err = ",x-owcdm_romberg(1.0,0.3,-0.2,0.05) 

    x = owcdm(1.0,0.3,-0.4,0.05)[0]
    print "Non-flat WCDM, om=0.3, w=-0.4, ok=0.05, z=1: r(z) = ", x, " Err = ",x-owcdm_romberg(1.0,0.3,-0.4,0.05) 

    x = owcdm(1.0,0.3,-1.4,0.05)[0]
    print "Non-flat WCDM, om=0.3, w=-1.4, ok=0.05, z=1: r(z) = ", x, " Err = ",x-owcdm_romberg(1.0,0.3,-1.4,0.05) 

    print "\nSlightly closed univeses"
    x = owcdm(1.0,0.3,-1,-0.05)[0]
    print "Non-flat LCDM, om=0.3, ok=-0.05, z=1: r(z) = ", x, " Err = ",x-owcdm_romberg(1.0,0.3,-1,-0.05) 

    x = owcdm(1.0,0.3,-0.2,-0.05)[0]
    print "Non-flat WCDM, om=0.3, w=-0.2, ok=-0.05, z=1: r(z) = ", x, " Err = ",x-owcdm_romberg(1.0,0.3,-0.2,-0.05) 

    x = owcdm(1.0,0.3,-0.4,-0.05)[0]
    print "Non-flat WCDM, om=0.3, w=-0.4, ok=-0.05, z=1: r(z) = ", x, " Err = ",x-owcdm_romberg(1.0,0.3,-0.4,-0.05) 

    x = owcdm(1.0,0.3,-1.4,-0.05)[0]
    print "Non-flat WCDM, om=0.3, w=-1.4, ok=-0.05, z=1: r(z) = ", x, " Err = ",x-owcdm_romberg(1.0,0.3,-1.4,-0.05) 


    print "\nTesting the t(z) code (expect O(1e-5) because control version is not so accurate)"
    x = wcdm_time(0.0,1.0,-1.4)
    print "CDM, om=1.0, w=-1.4, z=0: r(z) = ", x, " Err = ",x-time_romberg(0.0,1.0,-1.4) 

    x = wcdm_time(0.0,0.3,-0.2)
    print "WCDM, om=0.3, w=-0.2, z=0: r(z) = ", x, " Err = ",x-time_romberg(0.0,0.3,-0.2) 

    x = wcdm_time(0.0,0.3,-1.4)
    print "WCDM, om=0.3, w=-1.4, z=0: r(z) = ", x, " Err = ",x-time_romberg(0.0,0.3,-1.4) 

    x = wcdm_time(1.0,1.0,-1.4)
    print "CDM, om=1.0, w=-1.4, z=1: r(z) = ", x, " Err = ",x-time_romberg(1.0,1.0,-1.4) 

    x = wcdm_time(1.0,0.3,-1.4)
    print "WCDM, om=0.3, w=-1.4, z=1: r(z) = ", x, " Err = ",x-time_romberg(1.0,0.3,-1.4) 


    print "\nTesting the radiation code"
    x = wcdm_rad(1.0,0.3,-1.0,8e-5)
    print "WCDM, om=0.3, w=-1.0, or=8e-5, z=1: r(z) = ", x, " Err = ",x-wcdmrad_romberg(1.0,0.3,-1.0,8e-5) 
    x = wcdm_rad(1.0,0.3,-0.4,8e-5)
    print "WCDM, om=0.3, w=-0.4, or=8e-5, z=1: r(z) = ", x, " Err = ",x-wcdmrad_romberg(1.0,0.3,-0.4,8e-5) 
    x = wcdm_rad(10.0,0.3,-1.0,8e-5)
    print "WCDM, om=0.3, w=-1.0, or=8e-5, z=10: r(z) = ", x, " Err = ",x-wcdmrad_romberg(10.0,0.3,-1.0,8e-5) 


    print "\nTesting the vectorization"
    z = np.arange(0.1,1.1,0.1)
    w = np.arange(-1.0,-0.4,0.1)
    om = np.arange(1.0,0.2,-0.1)
    om[0] -= 1e-7
    om2 = np.arange(1.4,0.2,-0.2)
    om2[2] = 1.0

    print "z: ", z
    print "owcdm(): ", owcdm(z, 0.3, -1, 0)[0]
    print "wcdm(): ", wcdm(z, 0.3, -1)

    print 
    print "w: ", w
    print "owcdm(): ", owcdm(1.0, 0.3, w, 0)[0]
    print "wcdm(): ", wcdm(1.0, 0.3, w)

    print 
    print "om: ", om
    print "wcdm(): ", wcdm(1.0,om,-1)
    print "owcdm(): ", owcdm(1.0,om,-1,0)[0]
    print "owcdm(): ", propmotdis(1.0,om,-1,0)
    print "om: ", om2
    print "wcdm(): ", wcdm(1.0,om2,-1)

    print
    print "LCDM D_A: ", angdist(z, 0.3, -1, 0)*3000.0
    print "OCDM D_A: ", angdist(z, 0.3, -1, 0.1)*3000.0
    print "OCDM D_M: ", propmotdis(z, 0.3, -1, 0.1)

    print
    print "LCDM H*t(z): ", wcdm_time(z, 0.3, -1)
    print "wCDM H*t(z=1): ", wcdm_time(1, 0.3, w)
    print "wCDM H*t(z=0): ", wcdm_time(0, 0.3, w)
    print "LCDM H*t(z=0) om: ", wcdm_time(0, om, -1)
    print "SCDM H*t(z=0): ", wcdm_time(0, 1, -1)




if __name__ == "__main__":
    test()
