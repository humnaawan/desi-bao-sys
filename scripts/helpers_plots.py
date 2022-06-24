from astropy.table import Table
import random
from desispec.io import read_spectra, read_frame
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os
from settings import *

__all__ = ['compare_spectra_template_vs_calibrated', 'compare_spectra_groupedcoadded_vs_not']

# ------------------------------------------------------------------------
# function to plot template spectra and calibration spectra (From the Spectra files and the cframe files)
def compare_spectra_template_vs_calibrated(exposures_path, night, expid,
                                           showfig=True, savefig=False, outdir=None, simspec_ind=None):
    path = f'{exposures_path}/exposures/{night}/{expid}'
    data_simspec = Table.read(f'{path}/simspec-{expid}.fits', format='fits')
    simspec_hdul = fits.open(f'{path}/simspec-{expid}.fits')
    data_fibermap = Table.read(f'{path}/fibermap-{expid}.fits', format='fits')
    
    # print out some stats
    print(f'\nkeys in the sim-spec file: {data_simspec.keys()}\n')
    for key in ['TRUESPECTYPE', 'TEMPLATETYPE', 'TEMPLATESUBTYPE']:
        print(f'{key}: {np.unique(data_simspec[key].value.data, return_counts=True)}')

    print(f'\nkeys in the fibermap file: {data_fibermap.keys()}\n')
    
    # lets now plot the spectrum of a random object
    if simspec_ind is None:
        simspec_ind = random.choice(range(len(data_simspec)))
        
    # sanity check
    if data_simspec['TARGETID'][simspec_ind] != data_fibermap['TARGETID'][simspec_ind]:
        raise ValueError(f'simspec and fibermaps indices do NOT match: simspec-index = {simspec_ind}')
    else:
        targetid = data_simspec['TARGETID'][simspec_ind]
        
    # first plot the simspec spectrum
    plt.clf()
    fig, axes = plt.subplots(2, 1)
    axes[0].plot(simspec_hdul[0].data, simspec_hdul[1].data[simspec_ind, :], 'kx-', label='simspec spectrum')
    axes[1].plot(simspec_hdul[0].data, simspec_hdul[1].data[simspec_ind, :], 'kx-', label='simspec spectrum')
    
    # now the calibrated one - from the Spectra file
    # lets now extract the petal info
    petal_ind = data_fibermap['PETAL_LOC'][simspec_ind]
    data_spectra = read_spectra(f'{path}/spectra-{petal_ind}-{expid}.fits')
    ind_spectra = np.where(data_spectra.fibermap['TARGETID'] == targetid)[0]
    for band in 'brz':
        nwave = len(data_spectra.wave[band])
        axes[0].plot(data_spectra.wave[band], data_spectra.flux[band][ind_spectra, :].reshape(nwave),
                     label=f'calibrated: Spectra: {band}')

    # and now the calibrated one - from the cframe files
    for band in 'brz': 
        data_frame = read_frame(f'{path}/cframe-{band}{petal_ind}-{expid}.fits')
        ind_frame = np.where(data_frame.fibermap['TARGETID'] == targetid)[0]
        axes[1].plot(data_frame.wave, data_frame.flux[ind_frame][0, :], label=f'calibrated: cframe: {band}')
    
    # set up the title
    title = ''
    for key in ['TRUESPECTYPE', 'TEMPLATETYPE', 'TEMPLATESUBTYPE']:
        title += f'{key} {data_simspec[key][simspec_ind]} ; '
    title = title[0:-2]
    plt.suptitle(f'night {night} ; expid {expid}\n' +
                 f'target ID {targetid} ; petal {petal_ind} ; simspec-ind {simspec_ind}\n{title}', y=0.99)
    # pot details
    for ax in axes:
        ax.legend(bbox_to_anchor=(1,1))
        ax.set_ylabel('flux')
    axes[1].set_xlabel('wavelength')
    fig.set_size_inches(10, 6)
    # save plot and/or show
    if showfig:
        plt.show()
    if savefig:
        fname = f'plot_spectra_template-vs-calibrated_targetid{targetid}_simspec-ind{simspec_ind}.png'
        fig.savefig('%s/%s' % (outdir, fname), format='png', bbox_inches='tight')
        plt.close()
        print('## saved %s' % fname )

# ------------------------------------------------------------------------
# function to compare ungrouped/uncoadded spectra with grouped/coadded ones
def compare_spectra_groupedcoadded_vs_not(exposures_path, nside,
                                          showfig=True, savefig=False, outdir=None, hpixnum=None, targetid=None):
    path = f'{exposures_path}/exposures/spectra-{nside}'
    # pick one of the healpix numbers randomly
    if hpixnum is None:
        hpixnum = int(random.choice(os.listdir(path)))
    # now read the grouped data
    grouped = read_spectra(f'{path}/{hpixnum}/grouped-{nside}-{hpixnum}.fits')
    
    # now pick a random targetid
    if targetid is None:
        targetid = random.choice(np.unique(grouped.fibermap['TARGETID'].value))
    
    ind_grouped = np.where(grouped.fibermap['TARGETID'].value == targetid)[0]
    # now plot things out
    nrows, ncols = len(ind_grouped), 4
    print(f'{nrows} spectra found to be grouped')
    fig, axes = plt.subplots(nrows, ncols)
    plt.subplots_adjust(hspace=0.6, wspace=0.25)
    
    for i, ind in enumerate(ind_grouped):
        # lets plot the grouped spectra
        for band in 'brz':
            nwave = len(grouped.wave[band])
            if nrows == 1:
                ax = axes[2]
            else:
                ax = axes[i, 2]
            ax.plot(grouped.wave[band], grouped.flux[band][ind, :].reshape(nwave))
            ax.set_title(f'grouped observation; ind {ind}\nflux from grouped-spectra file', fontsize=14)
            
        # now plot the spectra file from disc that were presumbly used to create the grouped spectra
        petal = grouped.fibermap[ind]['PETAL_LOC']
        night, expid = grouped.fibermap['NIGHT'][ind], grouped.fibermap['EXPID'][ind]
        expid = f'00000{expid}'
        subdir = f'{exposures_path}/exposures/{night}/{expid}'
        
        # first the simspec file
        fname = f'simspec-{expid}.fits'
        print(f'reading in {fname}')
        simspec_hdul = fits.open(f'{subdir}/{fname}')
        simspec_ind = np.where(simspec_hdul[13].data['TARGETID'] == targetid)[0]
        if nrows == 1:
            ax = axes[0]
        else:
            ax = axes[i, 0]
        ax.plot(simspec_hdul[0].data, simspec_hdul[1].data[simspec_ind, :][0, :])
        ax.set_title(f'simspec template; ind {simspec_ind}\nnight {night}; expid {expid}', fontsize=14)

        # now the Spectra files   
        fnames = [f for f in os.listdir(subdir) if f.startswith(f'spectra-{petal}-')]
        if len(fnames) != 1:
            raise ValueError(f'somethings weird: have {len(fnames)} Spectra fits files in {subdir}')
        else:
            print(f'reading in {fnames[0]}')
            # now read in the spectra
            spectra = read_spectra(f'{subdir}/{fnames[0]}')
            ind_petal_spec = np.where(spectra.fibermap['TARGETID'] == targetid)[0]
            if nrows == 1:
                ax = axes[1]
            else:
                ax = axes[i, 1]
            for band in 'brz':
                nwave = len(spectra.wave[band])
                ax.plot(spectra.wave[band], spectra.flux[band][ind_petal_spec, :].reshape(nwave))
            ax.set_title(f'calibrated Spectra; ind {ind_petal_spec}\npetal {petal}; night {night}; expid {expid}', fontsize=14)
            
    grouped = []
    # now lets read the coadded spectrum
    coadded = read_spectra(f'{path}/{hpixnum}/spectra-{nside}-{hpixnum}.fits')
    ind_coadded = np.where(coadded.fibermap['TARGETID'].value == targetid)[0]
    if len(ind_coadded) != 1:
        raise ValueError(f'somethings wrong - expect 1 coadd entry but got {len(ind_coadded)}; targetid = {targetid}; hpix = {hpixnum}')
    else:
        ind = ind_coadded[0]
        if nrows == 1:
            ax = axes[3]
        else:
            ax = axes[0, 3]
            for i in range(1, nrows):
                axes[i, 3].axis('off')
        for band in 'brz':
            nwave = len(coadded.wave[band])
            ax.plot(coadded.wave[band], coadded.flux[band][ind, :].reshape(nwave))
            ax.set_title(f'coadded observation; ind {ind}\nflux from coadded-spectra file', fontsize=14)
    # plot details
    fig.set_size_inches(20, 3*nrows)
    if nrows <= 2:
        y = 1.25
    else:
        y = 1.01
    plt.suptitle(f'targetID {targetid}\n healpix pixel {hpixnum} (nside {nside})', fontsize=18, y=y, fontweight="bold")  
    if nrows == 1:
        for col in range(ncols):
            axes[col].set_xlabel('wavelength')
        axes[0].set_ylabel('flux')
    else:
        for col in range(ncols):
            axes[-1, col].set_xlabel('wavelength')
        axes[0, -1].set_xlabel('wavelength')
        for row in range(nrows):
            axes[row, 0].set_ylabel('flux')
    # save fig and/or show
    if showfig:
        plt.show()
    if savefig:
        fname = f'plot_spectra_template-to-coadd_targetid{targetid}_hpix{hpixnum}_nside{nside}.png'
        fig.savefig('%s/%s' % (outdir, fname), format='png', bbox_inches='tight')
        plt.close()
        print('## saved %s' % fname )