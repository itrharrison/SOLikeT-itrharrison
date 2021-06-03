"""
.. module:: xcorr

:Synopsis: Definition of cross-correlations likelihood for Simons Observatory
:Authors: Ian Harrison. From orginal xcorr by Andrina Nicola, Alex Krolewski, Emmanuel Schaan

"""
from ..gaussian import GaussianData, GaussianLikelihood
from ..ps import PSLikelihood
from .halofit import HalofitPowerSpectrum
from .utils import bin_cl, apply_shift, apply_width

class XcorrLikelihood(GaussianLikelihood):

    def initilize(self):
        name: str = "Xcorr"
        self.log.info('Initialising.')

        x, y, dy = self._get_data()
        if self.covpath is None:
            self.log.info('No Xcorr covariance specified. Using diag(dy^2).')
            cov = np.diag(dy**2)
        else:
            cov = self._get_cov()
        self.data = GaussianData(self.name, x, y, cov)

    def get_requirements(self):


    def _get_data(self, **params_values):
        data_auto = np.loadtxt(self.auto_file)
        data_cross = np.loadtxt(self.cross_file)

        # Get data
        self.ell_auto = data_auto[0]
        cl_auto = data_auto[1]
        cl_auto_err = data_auto[2]

        self.ell_cross = data_cross[0]
        cl_cross = data_cross[1]
        cl_cross_err = data_cross[2]  

        x = np.concatenate([self.ell_auto, self.ell_cross])
        y = np.concatenate([cl_auto, cl_cross])
        dy = np.concatenate([cl_auto_err, cl_cross_err])

        return x, y, dy

    def _get_theory(s_wise=0.4,
                    b1=1.0,
                    b2=0.0,
                    bs=0.0,
                    alpha_auto=0.0,
                    alpha_cross=0.0,
                    alpha_matter=0.0,
                    SN=1e-7,
                    shift=0.0,
                    width=1.0,
                    _theory={'Pk_interpolator': {'z': np.linspace(0,zmax,Nz), 
                             'k_max': kmax, 'nonlinear': False,
                             'hubble_units': False,'k_hunit': False, 
                             'vars_pairs': [[var,var]]}}):
        '''Function for the likelihood.'''
        t0 = time.time()
        Omegam = (_theory.get_param('ombh2') + _theory.get_param('omch2'))/((_theory.get_param('H0')/100.)**2.)
        
        Pk_interpolator = _theory.get_Pk_interpolator()['delta_nonu_delta_nonu'].P

        # Set up halofit object
        
        halofit = _setup_halofit(_theory, _theory.get_param('H0'), _theory.get_param('ombh2'),
                                 _theory.get_param('omch2'), _theory.get_param('nnu'),
                                 _theory.get_param('mnu'),_theory.get_param('num_massive_neutrinos'),
                                 _theory.get_param('As'),_theory.get_param('ns'))
        
        SN_auto = SN
        SN_cross = 0.   

        lim_cross, lim_auto = _get_angular_power_spectra(halofit, _theory, s_wise, Omegam, b1, b2, bs, alpha_auto, alpha_cross, alpha_matter, SN, shift, width)

        #print('lim_auto immediately before bin',lim_auto)

        # Bin
        lim_bin = bin_cl(all_ell,lim_cross,lmin[low_ind_kg:high_ind],lmax[low_ind_kg:high_ind]) + SN_cross
        # Correct for namaster binning
        lim_bin_cross = lim_bin #* clkg_corr[low_ind:high_ind]
        
        # Bin and add SN
        lim_bin = bin_cl(all_ell,lim_auto,lmin[low_ind_gg:high_ind],lmax[low_ind_gg:high_ind]) + SN_auto
        # Correct for namaster binning
        lim_bin_auto = lim_bin #* clgg_corr[low_ind:high_ind]

        return np.concatenate([lim_bin_auto, lim_bin_cross])

    def _get_data(self):

def _get_angular_power_spectra(halofit,_theory,
                                  s_wise=0.4,
                                  Omegam=0.3,
                                  b1=1.0,
                                  b2=0.0,
                                  bs=0.0,
                                  alpha_auto=0.0,
                                  alpha_cross=0.0,
                                  alpha_matter=0.0,
                                  SN=1e-7,
                                  shift=0.0,
                                  width=1.0,autoCMB=False):
    '''Returns angular power spectra'''
    t0 = time.time()
    Pk_interpolator = _theory.get_Pk_interpolator()['delta_nonu_delta_nonu'].P
    #print('pk interpolator',time.time()-t0)

    auto, cross = _wrap_limber(Pk_interpolator, halofit, s_wise, _theory.get_param('As'), _theory.get_param('H0'),
                                _theory.get_param('omch2'), _theory.get_param('ombh2'), _theory.get_param('mnu'),
                                _theory.get_param('nnu'), _theory.get_param('ns'), _theory.get_param('tau'),
                                _theory.get_param('num_massive_neutrinos'), shift, width, autoCMB)
        

    # Multiply by bias  
    lim_clkg, alpha_cross_term, lim_clkmu = cross
    lim_cross = b1* lim_clkg + lim_clkmu

    # Multiply by bias
    lim_clgg, alpha_auto_term, lim_clgmu, lim_clmumu = auto
    lim_auto = b1**2 * lim_clgg + 2*b1*lim_clgmu + lim_clmumu
    
    return np.squeeze(lim_cross), np.squeeze(lim_auto)


@cached(cache={}, key=lambda theory, H0, ombh2, omch2, nnu, mnu, num_massive_neutrinos, As, ns: 
        hashkey(H0, ombh2, omch2, nnu, mnu, num_massive_neutrinos, As, ns))
def _setup_halofit(theory, H0, ombh2, omch2, nnu, mnu, num_massive_neutrinos, As, ns):
    pars = theory.camb.CAMBparams()

    pars.set_cosmology(H0=H0,
        ombh2=ombh2,
        omch2=omch2,
        nnu=nnu,
        mnu=mnu,
        num_massive_neutrinos=num_massive_neutrinos)
    pars.InitPower.set_params(As=As,
        ns=ns)
    from camb import model as CAMBModel
    pars.NonLinear = CAMBModel.NonLinear_both
    
    halofit = theory.camb.get_matter_power_interpolator(pars,
        nonlinear=True,hubble_units=False,k_hunit=False,kmax=kmax,
        var1=var,var2=var).P
    return halofit

@cached(cache={}, key=lambda Pk_interpolator, halofit, s_wise, As, H0, omch2, ombh2, mnu, 
    nnu, ns, tau, num_massive_neutrinos, shift, width, autoCMB: hashkey(s_wise, As, H0, omch2, ombh2, mnu, nnu, ns, tau, num_massive_neutrinos, shift, width, autoCMB))
def _wrap_limber(Pk_interpolator, halofit, s_wise, As, H0, omch2, ombh2, mnu, nnu, ns, tau, num_massive_neutrinos, shift, width, autoCMB):
    t0 = time.time()
    halofit_pk = _wrap_halofit(Pk_interpolator, halofit, As, H0, omch2, ombh2, mnu, nnu, ns, tau, num_massive_neutrinos)
    #print('cleft pk',time.time()-t0)
    h = H0/100.
    #omch2 = h**2 * (Om0 - (mnu/93.14)/h**2.) - ombh2
    Om0 = (omch2 + ombh2)/(h**2.) + (mnu/93.14)/h**2.
    Ob0 = (ombh2)/h**2.
    Neff = nnu

    cosmo = FlatLambdaCDM(H0=H0,Om0=Om0,Ob0=Ob0,Tcmb0=2.7255,Neff=Neff) 

    # Modify dn/dz by a shift and a width
    out = apply_shift(np.array([dndz_xcorr[:,0],0.0035*np.ones_like(dndz_xcorr[:,0]),dndz_xcorr[:,1]]).T, color, shift)
    out = apply_width(np.array([dndz_xcorr[:,0], 0.0035*np.ones_like(dndz_xcorr[:,0]), out]).T, color, width)   
    dndz_xcorr_modified = np.array([dndz_xcorr[:,0], out]).T        

    out = apply_shift(np.array([dndz_xmatch[:,0],0.0035*np.ones_like(dndz_xmatch[:,0]),dndz_xmatch[:,1]]).T, color, shift)
    out = apply_width(np.array([dndz_xmatch[:,0], 0.0035*np.ones_like(dndz_xmatch[:,0]), out]).T, color, width) 
    dndz_xmatch_modified = np.array([dndz_xmatch[:,0], out]).T      

    # Pre-define chi- grid to save evaluation time
    setup_chi_out = limber.setup_chi(cosmo, dndz_xcorr_modified, dndz_xmatch_modified, Nchi, Nchi_mag)
    
    alpha_cross = lambda z: 1.0
    alpha_auto = lambda z: 1.0
    
    out = limber.do_limber(all_ell, cosmo, dndz_xcorr_modified, dndz_xcorr_modified, s_wise, s_wise, halofit_pk, b1_HF, b1_HF, alpha_auto, alpha_cross, Nchi=Nchi, autoCMB=autoCMB, use_zeff=False, dndz1_mag=dndz_xmatch, dndz2_mag=dndz_xmatch,setup_chi_flag=True,setup_chi_out=setup_chi_out)   
        
    return out

@cached(cache={}, key=lambda Pk_interpolator, halofit, As, H0, omch2, ombh2, mnu, 
    nnu, ns, tau, num_massive_neutrinos: hashkey(As, H0, omch2, ombh2, mnu, nnu, ns, tau, num_massive_neutrinos))
def _wrap_halofit(Pk_interpolator, halofit, As, H0, omch2, ombh2, mnu, nnu, ns, tau, num_massive_neutrinos):
    # This is a wrapper for the limber call of the Halofit portion of the power spectrum
    # The idea here is to ensure that at fixed cosmology, the limber integral is cached
    # so that varying the bias and shot noise is as simple as multiplying by a constant
    # and adding another constant.
    # This requires the use of the "@cached" decorator from cachetools
    # see https://stackoverflow.com/questions/30730983/make-lru-cache-ignore-some-of-the-function-arguments
    # The idea is that you pass all of the cosmological parameters to this function and the power spectrum interpolator
    # and it caches the output with a key given by all of the parameters
    # If you add more parameters into your cosmology, don't forgot to add them to this block of code!
    t0 = time.time()
    h = H0/100.
    #omch2 = h**2 * (Om0 - (mnu/93.14)/h**2.) - ombh2
    Om0 = (omch2 + ombh2)/(h**2.) + (mnu/93.14)/h**2.
    Ob0 = (ombh2)/h**2.
    Neff = nnu

    # This cosmology is needed for distance-redshift conversions
    # since we are at low redshift, I include neutrinos in the matter density (assume they are non-relativistic)
    cosmo = FlatLambdaCDM(H0=H0,Om0=Om0,Ob0=Ob0,Tcmb0=2.7255,Neff=Neff) 


    setup_chi_out = limber.setup_chi(cosmo, dndz_xcorr, dndz_xcorr, Nchi, Nchi_mag)

    
    # Only the *internal* number of k points affects the runtime
    # These points are just interpolated between, so the runtime doesn't care about them
    # these should just be something reasonable, where I get reasonable answers out
    k = np.logspace(-4,np.log10(10),1000)
    
    halofit_pk = HalofitPowerSpectrum(k, Pk_interpolator, halofit, zlist) #, zelda)
    
    return halofit_pk


