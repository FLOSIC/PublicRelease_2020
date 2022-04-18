from pyscf import gto,dft
from scipy.optimize import least_squares
import numpy as np
import pandas

tmpl = ['X only', 'XC']
spind = {'N': 3, 'Ne': 0}
bas_l = ['ccpvdz','ccpvtz','ccpvqz','ccpv5z']
z = np.array([2,3,4,5])

def fit_obj(c,zmax,en):
    ffn = c[0] + c[1]*np.exp(-c[2]*zmax)
    #ffn = c[0] + c[1]*(zmax+1)*np.exp(-c[2]*(zmax)**(0.5))
    return en - ffn

tmpv = np.zeros(z.shape[0]+1)

xclw = pandas.ExcelWriter('./atom_calcs/r2scan_cbs_extrap.xlsx',\
    engine='xlsxwriter')

for at in ['N','Ne']:

    tdict = {'Basis': bas_l+['CBS'], 'Energy (X only, Ha)': tmpv,
        'Energy (XC, Ha)': tmpv}

    for ic in range(2):

        if ic == 0:
            wentry = 'Energy (X only, Ha)'
        elif ic == 1:
            wentry = 'Energy (XC, Ha)'

        for ibas,bas in enumerate(bas_l):

            srep = gto.M(atom='{:} 0 0 0'.format(at), basis=bas, symmetry=False,\
                spin = spind[at],cart=True)

            srep.cart2sph_coeff(normalized='sp')

            wrep = dft.UKS(srep)
            wrep.conv_tol = 1.e-8
            wrep.grids.level = 9
            dfa_str = '1.0*MGGA_X_R2SCAN,{:}*MGGA_C_R2SCAN'.format(1.0*ic)
            dft.libxc.define_xc_(wrep._numint,dfa_str)
            #dft.libxc.define_xc_(wrep._numint,'1.0*GGA_X_PBE,{:}*GGA_C_PBE'.format(1.0*ic))
            wrep.verbose = 1
            tdict[wentry][ibas] = wrep.kernel()

        lsd = least_squares(fit_obj,np.ones(3),args=(z,tdict[wentry][:-1]))
        tdict[wentry][-1] = lsd.x[0]

    tmpdt = pandas.DataFrame(tdict)
    tmpdt.to_excel(xclw,sheet_name=at,index=False)

xclw.save()
