# Binning function
def bin_cl(ell,theory_cl,lmin,lmax):
	binned_theory_cl = np.zeros_like(lmin)
	for i in range(len(lmin)):
		binned_theory_cl[i] = np.mean(theory_cl[(ell >= lmin[i]) & (ell < lmax[i])])
	return binned_theory_cl

def apply_shift(bdndz, color, deltaz):
	z = bdndz[:,0]
	dz = bdndz[:,1]
	data = bdndz[:,2]
	
	dndz = normalize_dndz(bdndz, color)
		
	shifted_dndz = np.interp(z+deltaz,z,dndz,left=0,right=0)
	
	new_bdndz = np.array([z, dz, shifted_dndz]).T
	
	out = normalize_dndz(new_bdndz, color)
	#out[0] = 0
	#out[-1] = 0
	
	return out
	
def apply_width(bdndz, color, width):
	z = bdndz[:,0]
	dz = bdndz[:,1]
	data = bdndz[:,2]
	
	dndz = normalize_dndz(bdndz, color)
	
	zbar = np.sum(z * dndz * dz)/np.sum(dndz * dz)
		
	shifted_dndz = np.interp(width * (z-zbar) + zbar,z,dndz,left=0,right=0)
	
	new_bdndz = np.array([z, dz, shifted_dndz]).T
	
	out = normalize_dndz(new_bdndz, color)
	#out[0] = 0
	#out[-1] = 0
	
	return out