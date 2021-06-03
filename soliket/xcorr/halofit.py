import time
import numpy as np

class HalofitPowerSpectrum():
    """A class which returns P(k) from CAMB linear P(k)"""
    # The main purpose of this class is to smoothly handle the
    # interpolation in z.  If we have a small-ish number of points
    # Lagrange interpolating polynomials are ideal for vectorized ops.
    # Currently P(k) is computed from linear P(k) from CAMB.
    def lagrange_spam(self,z):
        """Returns the weights to apply to each z-slice to interpolate to z."""
        dz = self.zlist[:,None] - self.zlist[None,:]
        singular = (dz == 0)
        dz[singular] = 1.0
        fac = (z - self.zlist) / dz
        fac[singular] = 1.0
        return(fac.prod(axis=-1))
    def __init__(self,k,Pk_interpolator,halofit,zlist):
        """Makes the P(k) tables for some k points,
        a Pk_interpolator object that calls the linear power spectrum,
        a list of redshifts. one_loop and shear tell you whether
        to include 1-loop and shear terms in the LPT expansion.
        Also adds a halofit object for the magnification term;
        halofit is a P(k) interpolator object with nonlinear=True."""
        self.halofit = halofit


    def __call__(self,kval,zz,b1_HF,b2_HF,one_loop=True,shear=True):
        """Returns power spectra at k=kval and z=zz.
        pars is a dict telling you which method to use (halozeldovich, zeft, halofit)
        and the bias/1-halo parameters for 2 tracers.  
        ex: pars = {'method': 'halofit', 'tracer1': {'b': lambda z: 1.5},
        'tracer2': {'b': lambda z: 2.0}}.
        Note that biases must be functions of redshift.
        This function then returns all 4 spectra needed for CMB lensing:
        tracer1 x tracer2 , tracer1 x matter, tracer2 x matter, and matter x matter.
        If you want auto-spectra, just use the same bias coefficients for tracer1 and tracer2."""
        # Get the interpolating coefficients then just weigh each P(k,z).
        t0 = time.time()
        
        p_mm = np.zeros((1,len(kval)))
        p_gm_HF = np.zeros((1,len(kval)))
        p_gg_HF = np.zeros((1,len(kval)))
        #print('coef',coef)

        #print('pmm',p_mm)
        #print('pgm1_HF',p_gm1_HF)
        p_mm[0,:] = self.halofit(zz, kval)
        p_gm_HF[0,:] = self.halofit(zz, kval)
        p_gg_HF[0,:] = self.halofit(zz, kval)
                
        #print('kval',np.shape(kval))
        #print('kk',np.shape(kk))

        return(p_gg_HF, p_gm_HF, p_mm)