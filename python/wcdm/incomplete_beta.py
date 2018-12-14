# This is an adaptation by Daniel Eisenstein of GSL code in 
# beta_cont_frac_gsl(a, b, x).  

# This file is distributed under the GNU Public License

import numpy as np

######### An adaptation of the GSL beta function code

def beta_cont_frac_gsl(a, b, x):
    # This computes B_x(a,b) using a continued fraction approximation.
    # We do require a>0, but b can be negative.
    # Having b<0 and x very near 1 can cause of precision (not enough iterations)
    # However, x near 1 only occurs in the far future for our cosmology application
    # Require 0<=x<1.
    #
    # This python subroutine is adapted from the 
    # Gnu Science Library (GSL) specfunc/beta_inc.c code
    # by Daniel Eisenstein (July 2015).
    # Changes were generally to strip down to the case of interest, removing
    # the pre-factor from the complete beta function.  Also vectorized.
    #
    # Original GSL header:
    # Copyright (C) 2007 Brian Gough
    # Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
    #
    # This program is free software; you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation; either version 3 of the License, or (at
    # your option) any later version.
    # 
    # This program is distributed in the hope that it will be useful, but
    # WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    # General Public License for more details.
    # 
    # You should have received a copy of the GNU General Public License
    # along with this program; if not, write to the Free Software
    # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
    #
    # Author:  G. Jungman 

    x = np.array(x, copy=False, ndmin=1)
    if (np.min(x)<0 or np.max(x)>=1):
    	print "Illegal entry in beta_cont_frac_gsl()\n"
	exit
    cutoff = 1e-30	#  control the zero cutoff

    # standard initialization for continued fraction 
    num_term = 1.0
    den_term = 1.0 - (a+b)*x/(a+1.0)
    den_term[np.where(np.abs(den_term)<cutoff)] = cutoff
    den_term = 1.0/den_term
    cf = den_term

    for k in range(1,200):
	# first step 
	coeff = k*(b-k)*x/(((a-1.0)+2*k)*(a+2*k))
	den_term = 1.0 + coeff*den_term
	num_term = 1.0 + coeff/num_term
	den_term[np.where(np.abs(den_term)<cutoff)] = cutoff
	num_term[np.where(np.abs(num_term)<cutoff)] = cutoff
	den_term  = 1.0/den_term
	cf *= den_term*num_term

	# second step 
	coeff = -(a+k)*(a+b+k)*x/((a+2*k)*(a+2*k+1.0))
	den_term = 1.0 + coeff*den_term
	num_term = 1.0 + coeff/num_term
	den_term[np.where(np.abs(den_term)<cutoff)] = cutoff
	num_term[np.where(np.abs(num_term)<cutoff)] = cutoff
	den_term = 1.0/den_term
	cf *= den_term*num_term

	# Are we done?
	if (np.max(np.abs(den_term*num_term-1))<1e-12): break
    # End k loop
    # If this ends, we're just accepting the answer even if we haven't converged

    # Include the prefactor 
    # We need a>0 so that x=0 doesn't crash.
    cf *= np.power(x,a)*np.power(1-x,b)/a
    if (len(cf)==1): return cf[0]     # Get back to a scalar
    else: return cf



