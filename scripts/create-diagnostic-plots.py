from astropy.table import Table
import random
import matplotlib.pyplot as plt
import numpy as np
import os
from settings import *

from helpers_plots import compare_spectra_template_vs_calibrated, compare_spectra_groupedcoadded_vs_not

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--outdir',
                  dest='outdir', type='str',
                  help='path to directory where things should be saved.')
parser.add_option('--exposures-path',
                  dest='exposures_path', type='str',
                  help='path to directory with all exposures folder with night sub-folders.')
parser.add_option('--zcat-path',
                  dest='path_zcat', type='str',
                  help='path to zca file.')
parser.add_option('--spectra-path',
                  dest='path_spec16', type='str',
                  help='path to spectra-<nside> folder.')
parser.add_option('--nside',
                  dest='nside', type=int,
                  help='HEALPix resolution param')
parser.add_option('--rand',
                  action='store_true', dest='rand_obj', default=False,
                  help='pick a random object to plot some details for' +
                        'otherwise going for a high-z Lya one.')

(options, args) = parser.parse_args()
outdir = options.outdir
path_exposures = options.exposures_path
path_zcat = options.path_zcat
path_spec16 = options.path_spec16
nside = options.nside
randobj = options.rand_obj

# ------------------------------------------------------------------------------
# lets deal with the zcatalog stuff first
zcat = Table.read(path_zcat)
zcat = zcat[zcat['HPIXELNUM'] != 0]  # dealing with duplicates -- needs to be fixed
grps_spectype = zcat.group_by('SPECTYPE')
grps_subtype = zcat.group_by('TEMPLATESUBTYPE')

for i, grp in enumerate(grps_spectype.groups.keys['SPECTYPE']):
    print(f'group = {grp}\n  {len(grps_spectype.groups[i])} objects\n')

print('--')
for i, grp in enumerate(grps_subtype.groups.keys['TEMPLATESUBTYPE']):
    print(f'group = {grp}\n  {len(grps_subtype.groups[i])} objects\n')

# qso catalog
zcat_qso = grps_spectype.groups[np.where(grps_spectype.groups.keys['SPECTYPE'] == 'QSO')]
# lya catalog
zcat_lya = grps_subtype.groups[np.where(grps_subtype.groups.keys['TEMPLATESUBTYPE'] == 'LYA')]
# ----------------------------------
# n(z) plot
bins = np.arange(0, 3.8, 0.05)
plt.hist(zcat_qso['Z'], histtype='step', bins=bins, lw=4, label='QSO (%s objects)' % len(zcat_qso['Z']))
plt.hist(zcat_lya['Z'], histtype='step', bins=bins, lw=2, label='Lya (%s objects)' % len(zcat_lya['Z']))
plt.title('QSO catalog')
plt.legend(loc='best')
plt.xlabel('z')
plt.ylabel('counts')
# save fig
fname = f'{outdir}/plot_qso-cat_z-hist.png'
plt.savefig(fname, format='png', bbox_inches='tight')
plt.close()
print(f'## saved {fname}')
# ----------------------------------
# now the ra-dec plot
fig, axes = plt.subplots(1, 2)
plt.subplots_adjust(wspace=0.1)
axes[0].plot(zcat_qso['RA'], zcat_qso['DEC'], '.')
axes[1].plot(zcat_lya['RA'], zcat_lya['DEC'], '.')
# plot details
axes[0].set_title('QSO (%s objects)' % len(zcat_qso['Z']))
axes[1].set_title('Lya (%s objects)' % len(zcat_lya['Z']))
# labels
for nrow in range(2):
    axes[nrow].set_xlabel('RA (deg)')
axes[0].set_ylabel('DEC (deg)')
fig.set_size_inches(20, 5)
# save fig
fname = f'{outdir}/plot_qso-cat_ra-dec.png'
plt.savefig(fname, format='png', bbox_inches='tight')
plt.close()
print(f'## saved {fname}')
# ----------------------------------
# now the supposed flux plot
bins = np.arange(0, 40, 0.5)
upp = max(bins)

fig, axes = plt.subplots(1, 2)
plt.subplots_adjust(wspace=0.1)
for band in 'GRZ':
    flux = zcat_qso[f'FLUX_{band}_1']
    axes[0].hist(flux, histtype='step', bins=bins, lw=2,
                 label=f'{band} ({len(np.where(flux > upp)[0])} objects with flux > {upp}')
    flux = zcat_lya[f'FLUX_{band}_1']
    axes[1].hist(flux, histtype='step', bins=bins, lw=2,
                 label=f'{band} ({len(np.where(flux > upp)[0])} objects with flux > {upp}')
# plot details
axes[0].set_title('QSO (%s objects)' % len(zcat_qso['Z']))
axes[1].set_title('Lya (%s objects)' % len(zcat_lya['Z']))
for nrow in range(2):
    axes[nrow].legend(loc='best')
    axes[nrow].set_xlabel('flux')
axes[0].set_ylabel('counts')
fig.set_size_inches(20, 5)
# save fig
fname = f'{outdir}/plot_qso-cat_fluxes.png'
plt.savefig(fname, format='png', bbox_inches='tight')
plt.close()
print(f'## saved {fname}')
# ------------------------------------------------------------------------------
# now plot things for the spectra
print('\n## setting up for the spectra plots ..')
if randobj:
    hpix = None
    targetid = None
else:
    # read  in zcat - and randomly choose a target id
    zcat_lya = zcat_lya[zcat_lya['Z'] > 2.5]
    targetid = random.choice(np.unique(zcat_lya['TARGETID'].value))
    hpix = zcat_lya['HPIXELNUM'][np.where(zcat_lya['TARGETID'].data == targetid)[0]].value[0]

print(f'targetid = {targetid}; hpix = {hpix}')
# ----------------------------------
# plot the simspectra vs the grouped/coadded ones
out = compare_spectra_groupedcoadded_vs_not(exposures_path=path_exposures,
                                            coadd_spectra_path=path_spec16,
                                            nside=nside, hpixnum=hpix, targetid=targetid,
                                            outdir=outdir, savefig=True,
                                            return_objdata=True
                                            )
simspec_ind, night, expid = out
print('etwet', simspec_ind)
# ----------------------------------
# now plot the template vs. calibrated plot for this object
compare_spectra_template_vs_calibrated(exposures_path=path_exposures,
                                       night=night, expid=expid,
                                       outdir=outdir, savefig=True,
                                       simspec_ind=simspec_ind
                                       )
