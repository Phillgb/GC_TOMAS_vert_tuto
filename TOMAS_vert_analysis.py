#!/usr/bin/python3
# Author: Phillipe Gauvin-Bourdon
# Created: April 03, 2023
# ------------------------------------------------------------------------------
# Import modules
# ------------------------------------------------------------------------------
import numpy as np
import xarray as xr 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from GC_TOMAS_functions import GC_sig2press
# ------------------------------------------------------------------------------
# Setup of inputs variables
# ------------------------------------------------------------------------------
# OUTPUT files parameters
dir = './'
prefix = 'GEOSChem.SpeciesConc.'
mns = ['20200301', '20200401', '20200501', '20200601', '20200701']
sdhead = ['NK','SF','SS','ECIL','ECOB','OCIL','OCOB','DUST','AW']
nbins = 15  # Number of TOMAS bins
nlevs = 47  # Number of pressure levels

# Molecular weights (SF, SS, ECIL, ECOB, OCIL, OCOB, DUST)
molwgt = [96., 58.5, 12., 12., 12., 12., 100., 18.]      # [g/mol]

# Dry densities of the aerosol species
dens = [1770.,2000.,2000.,1400.,1500.]      #kg m-3 [so4, ss, bc, oa, dust]

# Retrieve GC outputs
print('Importing data...')
dist = []
for i, m in enumerate(mns):
    fname = dir+prefix+m+'_0000z.nc4'
    nc_file = xr.open_dataset(fname)
    data = []
    for ic, c in enumerate(sdhead):
        sizedist = []
        for k in range(0, nbins):
            sizedist.append(nc_file['SpeciesConcVV_'+c+str(k+1)])
        data.append(sizedist)
    dist.append(data)

# Store all distributions concentrations in a numpy array with 
# Dimensions [month, composition, size bin, level, lat, lon]
# Composition are in the order defined by 'sdhead'
dist = np.array(dist)[:, :, :, 0, :, :, :]

# Retreive GC sigma pressure levels
# for the midpoints
etam = nc_file['lev']
# for the edges of the levels
etae = nc_file['ilev']
psurf = nc_file['P0'] # mb
ptop = 0.01 # mb

# Extract coordinates of the netcdf grid
lat = np.array(nc_file['lat'])
long = np.array(nc_file['lon'])

# ------------------------------------------------------------------------------
# GEOS-Chem outputs conversions
# ------------------------------------------------------------------------------
print('GEOS-Chem output convertion...')
# Convert pressure levels
pressm = GC_sig2press(etam, psurf, ptop)
presse = GC_sig2press(etae, psurf, ptop)

# Convert GC particle mixing ratio to concentrations [cm^-3] STP
dist[:, 0, :, :, :, :] = dist[:, 0, :, :, :, :] / 22.4E6 * 273 / 293   #comp, bin, level, lat, lon

# Convert mass species to µg m-3
for ic, c in enumerate(molwgt):
    dist[:, ic+1, :, :, :, :] = dist[:, ic+1, :, :, :, :] * c / 22.4E6 * 1E6 * 1E9

# Convert TOMAS bins limits to diameter instead of mass
# Recreation of the mass bin edges
# This is defined for all TOMAS simulations
xk = np.zeros(nbins+1)  # particle mass cutoffs [kg] ... 15 bins makes 16 bins edges
if nbins == 15:
    xk[0] = 1.6E-23         # avg mass per particle in lowest bin edge
    for i in range(1, nbins+1):
        if i < nbins-1:
            xk[i] = xk[i-1]*4
        else:
            xk[i] = xk[i-1]*32
elif nbins == 40:
    for i in range(1, nbins+1):
        xk[i] = xk[i-1]*2

rho_assumed = 1400.     # [kg m^-3] assumed density

vol = xk/rho_assumed
Dpk1 = (6.*vol/np.pi)**(1./3.)*1E6      # particle diameter cutoffs [µm]
Dp_lower = Dpk1[:-1]                    # Extract lower limits
Dp_upper = Dpk1[1:]                     # Extract upper limits

# Calculation of diameter for the center of the bin
x = np.sqrt(xk[:-1]*xk[1:])
vol = x/rho_assumed
Dp1 = (6.*vol/np.pi)**(1./3.)*1E6       # particle diameter center [µm]

# aero_num = dist[:, 0, :, :, :, :]*1E6                   # [m^-3] -- first composition is number

# SO4 = dist[:, 1,:,:,:,:]*1E-9                           # [kg m-3] -- SO4
# SS = dist[:, 2,:,:,:,:]*1E-9                            # [kg m-3] -- sea salt
# BC = (dist[:, 3,:,:,:,:] + dist[:, 4,:,:,:,:])*1E-9     # [kg m-3] -- internally mixed and externally mixed BC
# OA = (dist[:, 5,:,:,:,:] + dist[:, 6,:,:,:,:])*1E-9     # [kg m-3] -- hydrophilic and hydrophobic OA
# DU = dist[:, 7,:,:,:,:]*1E-9                            # [kg m-3] -- dust

# dry_mass = SO4[:, :, :, :, :] + SS[:, :, :, :, :] + BC[:, :, :, :, :] + \
#             OA[:, :, :, :, :] + DU[:, :, :, :, :]                # total dry mass [µg m^-3]
# tot_pvol = (SO4[:, :, :, :, :]/aero_num[:, :, :, :, :]/dens[0]) +\
#             (SS[:, :, :, :, :]/aero_num[:, :, :, :, :]/dens[1]) +\
#             (BC[:, :, :, :, :]/aero_num[:, :, :, :, :]/dens[2]) +\
#             (OA[:, :, :, :, :]/aero_num[:, :, :, :, :]/dens[3]) +\
#             (DU[:, :, :, :, :]/aero_num[:, :, :, :, :]/dens[4])     # total volume per particle

# dry_mass_per = dry_mass/aero_num            # dry mass per particle

# Dp = (6.*tot_pvol/np.pi)**(1/3)             # [m] average particle diamter of bin
# Dp_lower = Dp[:, :, :, :, :] * (xk[:-1, np.newaxis, np.newaxis, np.newaxis, np.newaxis]/dry_mass_per[:, :, :, :, :])**(1/3) # [m]
# Dp_upper = Dp[:, :, :, :, :] * (xk[1:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]/dry_mass_per[:, :, :, :, :])**(1/3)  # [m]

# Convert to micron
# Dp = Dp*1E9             # [nm]
# Dp_lower = Dp_lower*1E9 # [nm]
# Dp_upper = Dp_upper*1E9 # [nm]

# ------------------------------------------------------------------------------
# GEOS-Chem outputs analysis
# ------------------------------------------------------------------------------
print('Preparing arrays for PEARL plotting...')
# Find index of cell corresponding to PEARL coord
ilat = np.where(lat == 82)[0]
ilong = np.where(long == -85)[0]

# Normalize the particle distributions concentrations at PEARL
PEARL_dist = dist[:, :, :, :, ilat, ilong]      # Extracting PEARL cell data
PEARL_spec_tot = np.array(PEARL_dist).sum(2)    # Calcul the total mass concentration for all species    
dry_mass_tot = PEARL_spec_tot[:, 1:-1, :, :].sum(1)   # Sum all dry species mass concentrations

# Normalize the mass concentration by the diameter range
dry_dMdlogDp_tot = dry_mass_tot / np.log10(Dp_upper[-1]/Dp_lower[0])    

# Create a continuous colorscheme for the TOMAS bins
lcolors = plt.cm.jet(np.linspace(0, 1, nbins))

# Create a continuous colorscheme for the TOMAS dry species
scolor = plt.cm.jet(np.linspace(0, 1, len(sdhead)-2))

# Plotting vertical distribution at PEARL from each TOMAS species
for i, m in enumerate(mns):
    for j, c in enumerate(sdhead):
        print(f'Plotting vertical distribution of {c} during {m[:-2]}...')

        plt.figure(figsize=(10, 6))
        for n in range(0, nbins):
            vert_dist = np.array(PEARL_dist[i, j, n, :]).flatten(order='C')
            plt.plot(vert_dist, pressm, color=lcolors[n],
                    label=f'{Dp_lower[n]:.2e} - {Dp_upper[n]:.2e} µm')
        plt.plot(PEARL_spec_tot[i, j, :], pressm, linestyle='--', color='k', label='Total')
        plt.gca().invert_yaxis()
        plt.title(f'Vertical distribution of {c} at PEARL on {m}')
        if c == 'NK':
            plt.xlabel(r'Particle Number Concentration [$cm^{-3}$]')
        else:
            plt.xlabel(r'Particle Mass Concentration [$µg\ m^{-3}$]')
        plt.ylabel('Pressure (mb)')
        plt.legend(loc='lower right')

        plt.savefig(f'spec_vertprof/{m}_{c}_vert_profile.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Verticale distribution color map of the particle distribution
        print(f'Plotting vertical colormap of {c} concentration during {m[:-2]}...')
        
        vert_dist = PEARL_dist[i, j, :, :].T

        fig = plt.figure(figsize=(10, 6))
        cmesh = plt.pcolormesh(Dp1, pressm, vert_dist[0], norm=LogNorm(), cmap='jet')
        plt.xscale('log')
        plt.gca().invert_yaxis()
        plt.xlabel('Particle size (µm)')
        plt.ylabel('Pressure (mb)')
        plt.minorticks_on()

        # Create a new axes for the colorbar
        cbar_ax = fig.add_axes([0.905, 0.1, 0.025, 0.78])
        # Create the colorbar and configure it
        cb = plt.colorbar(cmesh, norm=LogNorm(), cmap='jet', cax=cbar_ax)
        if c == 'NK':
            cb.set_label(r'dNdlogDp ($cm^{-3}$)', rotation=270, labelpad=20, fontsize=14)
        else:
            cb.set_label(r'dmdlogDp ($µg\ m^{-3}$)', rotation=270, labelpad=20, fontsize=14)
        plt.savefig(f'colormap/{m}_{c}_vertprof_cmap.png', dpi=300, bbox_inches='tight')
        plt.close()


    # Plotting the total mass concentrations observed at PEARL by TOMAS
    print(f'Plotting the vertical total mass concentrations during {m[:-2]}...')
    plt.figure(figsize=(10, 6))
    plt.plot(dry_mass_tot[i], pressm, color='k', linestyle='--', label='Total')
    for j, c in enumerate(sdhead[1:-1]):
        dMdlog10Dp = PEARL_spec_tot[i, j+1, :]
        plt.plot(dMdlog10Dp, pressm, color=scolor[j], label=c)
    plt.gca().invert_yaxis()
    plt.title(f'Vertical distribution of the total mass concentrations observed at PEARL')
    plt.xlabel(r'Particle mass concentration [$µg\ m^{-3}$]')
    plt.ylabel('Pressure (mb)')
    plt.legend(loc='lower right')

    plt.savefig(f'total_mass/{m}_masstot_vertprof.png', dpi=300, bbox_inches='tight')
    plt.close()