from astropy.table import Table
import random
from desispec.io import read_spectra, read_frame
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os
from settings import *
import fitsio
import scipy as sp
import h5py
from picca.wedgize import wedge

__all__ = ['compare_spectra_template_vs_calibrated',
           'compare_spectra_groupedcoadded_vs_not',
           'plot_delta_attributes', 'plots_cfs', 'plot_baofit'
           ]

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

def plot_delta_attributes(data_path, outdir):
    """
    data_path: str: path to directory that contains the delta_attributes.fits.gz file.
    """
    # adapted version of https://desi.lbl.gov/trac/attachment/wiki/LymanAlphaWG/how_to_run_picca/plot_delta_attributes.py
    data = fitsio.FITS(f'{data_path}/delta_attributes.fits.gz')

    f, (axs) = plt.subplots(nrows=2, ncols=3, figsize=(12,5))
    axs[-1, -1].axis('off')

    ### Stack
    loglam = data[1]['LOGLAM'][:]
    stack  = data[1]['STACK'][:]
    cut = (stack!=0.) & (data[1]['WEIGHT'][:]>0.)
    loglam = loglam[cut]
    stack  = stack[cut]
    axs[0][0].plot(10.**loglam, stack, linewidth=1)
    axs[0][0].grid()
    axs[0][0].set_xlabel(r'$\lambda_{\mathrm{Obs.}} \, [\AA]$',fontsize=20)
    axs[0][0].set_ylabel(r'$\mathrm{\overline{Deltas}}$',fontsize=20)

    ### ETA
    loglam    = data[2]['LOGLAM'][:]
    eta       = data[2]['ETA'][:]
    nb_pixels = data[2]['NB_PIXELS'][:]
    cut = (nb_pixels>0.)&(eta!=1.)
    loglam = loglam[cut]
    eta    = eta[cut]
    axs[0][1].errorbar(10.**loglam, eta, linewidth=1)
    axs[0][1].grid()
    axs[0][1].set_xlabel(r'$\lambda_{\mathrm{Obs.}} \, [\AA]$',fontsize=20)
    axs[0][1].set_ylabel(r'$\eta$',fontsize=20)

    ### VAR_LSS
    loglam    = data[2]['LOGLAM'][:]
    var_lss   = data[2]['VAR_LSS'][:]
    nb_pixels = data[2]['NB_PIXELS'][:]
    cut       = (nb_pixels>0.)&(var_lss!=0.1)
    loglam    = loglam[cut]
    var_lss   = var_lss[cut]
    axs[0][2].errorbar(10.**loglam, var_lss, linewidth=1)
    axs[0][2].grid()
    axs[0][2].set_xlabel(r'$\lambda_{\mathrm{Obs.}} \, [\AA]$',fontsize=20)
    axs[0][2].set_ylabel(r'$\sigma^{2}_{\mathrm{LSS}}$',fontsize=20)

    ### FUDGE
    loglam    = data[2]['LOGLAM'][:]
    fudge     = data[2]['FUDGE'][:]
    nb_pixels = data[2]['NB_PIXELS'][:]
    cut       = (nb_pixels>0.)&(fudge!=1.e-7)
    loglam    = loglam[cut]
    fudge     = fudge[cut]
    axs[1][0].errorbar(10.**loglam, fudge, linewidth=1)
    axs[1][0].grid()
    axs[1][0].set_xlabel(r'$\lambda_{\mathrm{Obs.}} \, [\AA]$',fontsize=20)
    axs[1][0].set_ylabel(r'$\mathrm{Fudge}$',fontsize=20)

    ### Mean cont
    loglam_rest = data[3]['LOGLAM_REST'][:]
    mean_cont   = data[3]['MEAN_CONT'][:]
    cut = (mean_cont!=0.) & (data[3]['WEIGHT'][:]>0.)
    loglam_rest = loglam_rest[cut]
    mean_cont   = mean_cont[cut]
    axs[1][1].plot(10.**loglam_rest, mean_cont, linewidth=1)
    axs[1][1].grid()
    axs[1][1].set_xlabel(r'$\lambda_{\mathrm{R.F.}} \, [\AA]$', fontsize=20)
    axs[1][1].set_ylabel(r'$\mathrm{\overline{Flux}}$', fontsize=20)

    plt.tight_layout()
    fname = f'{outdir}/plot_delta_attributes.png'
    plt.savefig(fname, format='png', bbox_inches='tight')
    plt.close()

    return fname

#
def plots_cfs(data_path, outdir, data_tag):
    # adapted version of https://desi.lbl.gov/trac/attachment/wiki/LymanAlphaWG/how_to_run_picca/plot_xcf.py
    """
    data_path: str: path to directory that contains e.g. the cf_lyalya_lyalya.fits.gz file.
    data_tag: str: e.g. lyalya_lyalya
    """
    data = fits.open(f'{data_path}/cf_lyalya_lyalya.fits.gz')

    #-- these are 100x50=5000 bins vectors containing the separations of each bin
    rp = data[1].data.RP   #-- this is radial separation r_parallel
    rt = data[1].data.RT   #-- this is transverse separation r_perp
    z = data[1].data.Z     #-- this is the mean redshift (almost constant)
    we = data[2].data.WE   #-- this is the weight associated to each bin, and each pixel
    cf = data[2].data.DA   #-- this is the correlation function in each bin and each pixel

    n_p = data[1].header['NP']
    n_t = data[1].header['NT']
    rpmin = data[1].header['RPMIN']
    rpmax = data[1].header['RPMAX']
    rtmax = data[1].header['RTMAX']

    rper = (np.arange(n_t)+0.5)*rtmax/n_t
    rpar = np.arange(rpmin+2, rpmax+2, (rpmax - rpmin)/n_p)

    #--  number of pixels and number of bins
    npix, nbins = cf.shape

    #-- doing a weighted average to get the final correlation function
    mcf = np.sum(cf*we, axis=0)/np.sum(we, axis=0)
    mwe = np.sum(we, axis=0)

    r = np.sqrt(rp**2+rt**2)

    #-- making a 2D plot of the correlation function
    plt.figure(figsize=(8, 7))
    plt.ion()
    step = int(len(rper))
    plt.pcolormesh( rper, rpar, (r*mcf).reshape(int(len(r)/step), step), vmin=-0.4, vmax=0.4, cmap=plt.cm.seismic)
    plt.xticks(np.arange(0, int(rpmax)+1, step=step))
    plt.yticks(np.arange(int(rpmin), int(rpmax)+1, step=step))
    plt.xlabel(r'$r_\perp [h^{-1} Mpc]$', fontsize=14)
    plt.ylabel(r'$r_\parallel [h^{-1} Mpc]$', fontsize=14)
    plt.colorbar()
    plt.title('cf_lyalya_lyalya')
    fname = f'{outdir}/plot_cf_{data_tag}.png'

    plt.savefig(fname, format='png', bbox_inches='tight')
    plt.close()

    return fname

def plot_baofit(path_cf, outdir, data_tag, path_baofit=None, power=2,
                 mus=[1., 0.95, 0.8, 0.5, 0], fit_rmax=200.,
                 absMus=True):
    # need to add doc strings for the input params
    # adapted from https://desi.lbl.gov/trac/attachment/wiki/LymanAlphaWG/how_to_run_picca/plot_wedges.py
    #- Read correlation function and covariance matrix
    h = fitsio.FITS(path_cf)
    da = h[1]['DA'][:]
    co = h[1]['CO'][:]
    hh = h[1].read_header()
    rpmin = hh['RPMIN']
    rpmax = hh['RPMAX']
    rtmin = 0
    rtmax = hh['RTMAX']
    nrp = hh['NP']
    nrt = hh['NT']
    h.close()

    #-- Read fit h5 file if there is one
    if path_baofit:
        ff = h5py.File(path_baofit, 'r')
        keys = list(ff.keys())
        for k in keys:
            if k != 'best fit' and k != 'minos':
                base = k
                break
        fit = ff[base+'/fit'][...]
        attr = dict(ff['best fit'].attrs)
        chi2 = attr['fval']
        ndata = attr['ndata']
        npars = attr['npar']
        rchi2 = chi2/(ndata-npars)
        print(f'chi2/(ndata-npars) = {chi2}/({ndata}-{npars}) = {rchi2}')
        ff.close()

    f, (axs) = plt.subplots(nrows=2, ncols=2, figsize=(12,6))

    j=0

    for i, (mumax,mumin) in enumerate(zip(mus[:-1],mus[1:])):
        b = wedge(mumin=mumin, mumax=mumax,
                                rpmin=rpmin, rpmax=rpmax,
                                rtmin=rtmin, rtmax=rtmax,
                                nrt=nrt, nrp=nrp, absoluteMu=absMus,
                                rmin=0., rmax=min(rpmax, rtmax),
                                nr=min(nrt, nrp))
        r,d,c = b.wedge(da,co)

        nrows = 2

        #-- Wedges and best model
        y = d*r**power
        dy = np.sqrt(c.diagonal())*r**power
        if absMus:
            axs[j//2][j%2].errorbar(
                r, y, dy, fmt="o",
                label=r"${}<|\mu|<{}$".format(mumin,mumax))
        else:
            axs[j//2][j%2].errorbar(
                r, y, dy, fmt="o",
                label=r"${}<\mu<{}$".format(mumin,mumax))
        if path_baofit:
            r, model, _ = b.wedge(fit, co)
            ym = model*r**power
            w = r < fit_rmax
            r = r[w]
            ym = ym[w]
            axs[j//2][j%2].plot(r, ym, linewidth=2, color="red", label='Fit')

        axs[j//2][j%2].set_ylabel(r"$r^{power}\xi(r)$".format(power=power))
        if j//2==1:
            axs[j//2][j%2].set_xlabel(r"$r \, [h^{-1}\, \mathrm{Mpc}]$")
        axs[j//2][j%2].legend(loc="upper right", fontsize=12)
        axs[j//2][j%2].grid(True)
        j+=1

        plt.tight_layout()

    fname = f'{outdir}/plot_wedges_{data_tag}.png'
    plt.savefig(fname, format='png', bbox_inches='tight')
    plt.close()

    return fname
