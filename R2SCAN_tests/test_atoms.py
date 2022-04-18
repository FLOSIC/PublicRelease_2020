from pyscf import gto,dft
from os import system,path
import pandas

if not path.isdir('./atom_calcs/'):
    system('mkdir ./atom_calcs')

# Reference values from r2 SCAN paper (DOI: 10.1021/acs.jpclett.0c02405) or the code repo:
# https://gitlab.com/dhamil/r2scan-subroutines

ref_d = {
    'R2SCAN': {'N': [-54.39774875432, -54.57900797069], # X only, then XC
        'Ne': [-128.5698768302, -128.9168416529]},
    'SCAN' : {'N': [-54.40541174994,-54.58565736367],
        'Ne': [-128.5891914977,-128.9340794821]}
}

tmpl = ['X only', 'XC']
spind = {'N': 3, 'Ne': 0}
bas_l = ['3-21G','6-311Gss','cc-pvDZ','cc-pvTZ']
key_order = ['Atom','Basis','E_tot (X only, Ha)','Reference (X only, cc-pvTZ, Ha)',\
'Diff (X only, Ha)','E_tot (XC, Ha)','Reference (XC, cc-pvTZ, Ha)','Diff (XC, Ha)']


for lcart in [False]:#,False]:

    if lcart:
        auxstr = '_cart'
    else:
        auxstr = '_sph'

    for cnorm in ['sp']:#['all','sp']:

        auxstr += '_{:}_norm'.format(cnorm)

        xclw = pandas.ExcelWriter('./atom_calcs/atom_calcs{:}.xlsx'.format(auxstr),\
            engine='xlsxwriter')

        for dfa in ['PBE','SCAN','R2SCAN']:

            tdict = {'Atom': [], 'Basis': [], 'E_tot (X only, Ha)': [], 'E_tot (XC, Ha)': []}
            if dfa in ref_d:
                tstr = 'Atom, Basis, X or XC, Energy (Ha), Reference (cc-pvTZ), Diff\n'
                tdict['Reference (X only, cc-pvTZ, Ha)'] = []
                tdict['Diff (X only, Ha)'] = []
                tdict['Reference (XC, cc-pvTZ, Ha)'] = []
                tdict['Diff (XC, Ha)'] = []
                tdict = {akey: tdict[akey] for akey in key_order}
            else:
                tstr = 'Atom, Basis, X or XC, Energy (Ha)\n'

            for at in ['N','Ne']:

                for bas in bas_l:

                    tdict['Atom'].append(at)
                    tdict['Basis'].append(bas)

                    for ic in range(2):

                        srep = gto.M(atom='{:} 0 0 0'.format(at), basis=bas, symmetry=False,\
                            spin = spind[at],cart=lcart)

                        srep.cart2sph_coeff(normalized=cnorm)

                        wrep = dft.UKS(srep)
                        wrep.conv_tol = 1.e-8
                        wrep.grids.level = 9
                        if dfa == 'PBE':
                            dfa_str = '1.0*GGA_X_{:},{:}*GGA_C_{:}'.format(dfa,1.0*ic,dfa)
                        else:
                            dfa_str = '1.0*MGGA_X_{:},{:}*MGGA_C_{:}'.format(dfa,1.0*ic,dfa)
                        dft.libxc.define_xc_(wrep._numint,dfa_str)
                        #dft.libxc.define_xc_(wrep._numint,'1.0*GGA_X_PBE,{:}*GGA_C_PBE'.format(1.0*ic))
                        wrep.verbose = 1
                        en = wrep.kernel()

                        if ic == 0:
                            tdict['E_tot (X only, Ha)'].append(en)
                        elif ic == 1:
                            tdict['E_tot (XC, Ha)'].append(en)

                        if dfa in ref_d:
                            rdiff = en - ref_d[dfa][at][ic]
                            rval = ref_d[dfa][at][ic]
                            tstr += ('{:},'*5 + '{:}\n').format(at, bas, tmpl[ic], en, \
                                rval, rdiff)
                            if ic == 0:
                                tdict['Reference (X only, cc-pvTZ, Ha)'].append(rval)
                                tdict['Diff (X only, Ha)'].append(rdiff)
                            elif ic == 1:
                                tdict['Reference (XC, cc-pvTZ, Ha)'].append(rval)
                                tdict['Diff (XC, Ha)'].append(rdiff)
                        else:
                            tstr += ('{:},'*3 + '{:}\n').format(at, bas, tmpl[ic], en)

            tmpdt = pandas.DataFrame(tdict)
            tmpdt.to_excel(xclw,sheet_name=dfa,index=False)

            with open('./atom_calcs/{:}_N_Ne_ref_calc{:}.csv'.format(dfa,auxstr),'w+') as tfl:
                tfl.write(tstr)

        xclw.save()
