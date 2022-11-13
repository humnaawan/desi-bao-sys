# ------------------------------------------------------------------------
# this script takes the simspec files that have contributed to the coadded
# spectra for each healpix pixel, and saves an analog of zbest.fits
# file; here we are not running redrock so this just contains the true zs.
# the goal is to then run desi_zcatalog to create the catalog that can then
# go to picca.
# ------------------------------------------------------------------------
import os
import sys
import time
import datetime
import logging
from astropy.io import fits
from astropy.table import Table, vstack, join
import fitsio
import numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--simspec-path',
                  dest='sims_path', type='str',
                  help='path to directory with all the simulated data \
                        (i.e., where folder name = night)')
parser.add_option('--coadds-path',
                  dest='coadds_path', type='str',
                  help='path to directory that contains the spectra-<nside> folder\
                        (i.e., where subfolders {hpix//100}{hpix})')
parser.add_option('--nside',
                  dest='nside', type=int,
                  help='HEALPix resolution param')
(options, args) = parser.parse_args()
sims_path = options.sims_path
coadds_path = options.coadds_path
nside = options.nside
# ------------------------------------------------------------------------
start0 = time.time()

# set up the logger
temp = f'save-zfile_{datetime.datetime.now()}'.replace(' ', '_').split('.')[0]
logging.basicConfig(filename=f'{coadds_path}/log_{temp}.log',
                    level=logging.DEBUG, filemode='w',
                   )
logging.info(f'running {sys.argv[0]}')
logging.info(f'\n ## inputs: {options}\n')

# some params for the output file
# minimal column headers according to https://github.com/desihub/desispec/blob/720153babcf85dd93530252b0c1f631d48edfc0d/bin/desi_zcatalog#L226
colnames = ['TARGETID', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z']
# add some more
colnames += ['Z', 'MOCKID', 'SPECTYPE', 'TEMPLATETYPE', 'TEMPLATESUBTYPE', 'TEMPLATEID']

for hpix_mod in os.listdir(f'{coadds_path}/spectra-{nside}/'):
    for hpix in os.listdir(f'{coadds_path}/spectra-{nside}/{hpix_mod}'):
        logging.info(f' ## working with healpix pixel= {hpix}')
        # read in the coadded spectra
        fname = f'{coadds_path}/spectra-{nside}/{hpix_mod}/{hpix}/spectra-{nside}-{hpix}.fits'
        hdul_coadded = fits.open(fname)
        # get the exposure fibermap - has NIGHT, EXPID
        fibermap_exp = hdul_coadded['EXP_FIBERMAP'].data
        # get the fibermap - has TARGET_RA, TARGET_DEC, FLUX_G, FLUX_R, FLUX_Z
        fibermap = hdul_coadded['FIBERMAP'].data

        # set things up for the output file
        data4later = {}
        for col in colnames:
            data4later[col] = []
        fibermap_all = None
        # now lets loop over all the different nights and exposures,
        # find the simpec file for all targets in this coadd, and extract the redshift.
        # another way to this will be to loop over all the targetid, find the
        # night + exposure and then extract the reshift. grouping things by night +
        # exposure should be faster.
        for night in np.unique(fibermap_exp['NIGHT']):
            logging.info(f' -- working with night = {night}')
            ind_night = np.where( fibermap_exp['NIGHT'] == night)[0]
            for expid in np.unique(fibermap_exp['EXPID']):
                logging.info(f' -- working with expid = {expid}')
                # find all the exposures
                inds_to_pull = np.where(fibermap_exp['EXPID'] == expid)[0]
                targetids = fibermap_exp['TARGETID'][inds_to_pull]
                if len(inds_to_pull) != len(np.unique(targetids)):
                    message = f' ## somethings wrong: n-exp ({len(inds_to_pull)}) '
                    message += f' != unique targets ({len(np.unique(targetids))})'
                    raise ValueError(message)

                # note: we dont have to worry about the petal info since we can
                #       just go directly to simspec file
                expid = f'00000{expid}'
                subdir = f'{sims_path}/{night}/{expid}'
                # read the file
                logging.info(f' -- reading in {fname}')
                fname = f'simspec-{expid}.fits'
                simspec_hdul = fits.open(f'{subdir}/{fname}')

                joined = join({'TARGETID': targetids}, simspec_hdul['TRUTH'].data, keys='TARGETID')
                # join reorders things so lets pull the correctly ordered targetids
                targetids = joined['TARGETID'].data
                # see if any of the targetids here are in the data4later
                intersection = set(targetids).intersection(set(data4later['TARGETID']))
                if len(intersection) != 0:
                    # this means that data4later has targetid(s) we see here
                    logging.info(f' -- {len(intersection)} targetids seen before')
                    # this is a problem if the redshift between now and before differ
                    for tid in intersection:
                        ind_here = np.where( joined['TARGETID'].data == tid)[0].astype(int)
                        ind_before = np.where(np.array(data4later['TARGETID']) == tid)[0].astype(int)

                        if (joined['TRUEZ'][ind_here].data != np.array(data4later['Z'])[ind_before]).any():
                            message = f' -- somethings wrong:\nz_here = {joined["TRUEZ"][ind_here].data } while '+ \
                                             f'z_before = {np.array(data4later["Z"])[ind_before]} ' + \
                                             f'for targetid = {tid}; inds = {ind_here} vs {ind_before}'
                            cols_ = ['TARGETID', 'Z', 'TRUESPECTYPE', 'TEMPLATETYPE', 'TEMPLATESUBTYPE', 'TEMPLATEID']
                            message += f'\nmore-data from before:\n{Table(data4later)[cols_][ind_before]}\n'
                            cols_ = ['TARGETID', 'TRUEZ', 'TRUESPECTYPE', 'TEMPLATETYPE', 'TEMPLATESUBTYPE', 'TEMPLATEID']
                            message += f'more-data from now:\n{joined[cols_][ind_here]}'
                            logging.critical(message)
                            raise ValueError(message)
                # okay things are okay - lets add to the data4later object
                data4later['TARGETID'] += list(targetids)
                data4later['Z'] += list(joined['TRUEZ'])
                data4later['SPECTYPE'] += list(joined['TRUESPECTYPE'])
                # lets also add so other stuff from the simspec table
                for col in ['MOCKID', 'TEMPLATETYPE', 'TEMPLATESUBTYPE', 'TEMPLATEID']:
                    data4later[col] += list(joined[col])

                # now need to pull in data from the fibermap
                # have to loop over the targetids
                inds_fibermap = []
                for targetid in targetids:
                    ind = np.where(fibermap['TARGETID'] == targetid)[0]
                    if len(ind) != 1:
                        message = f' -- somethings wrong - found {len(ind)} ' + \
                        f'entries for targetid = {targetid} in the fibermap'
                        logging.critical(message)
                        raise ValueError(message)
                    else:
                        inds_fibermap.append(ind)

                # lets do one more check - confirm that the fluxes are consistent
                # between the fibermap and the truth table
                for col in ['FLUX_G', 'FLUX_R', 'FLUX_Z']:
                    if (fibermap[col][inds_fibermap].flatten() != joined[col].data).any():
                        message = f' ## {col} not consistent between fibermap and truth: ' + \
                                         f'{fibermap[col][inds_fibermap].data} \nvs\n{joined[col].data}'
                        logging.critical(message)
                        raise ValueError(message)

                for col in ['RA', 'DEC']:
                    data4later[col] += list(fibermap[f'TARGET_{col}'][inds_fibermap])
                for col in ['FLUX_G', 'FLUX_R', 'FLUX_Z']:
                    data4later[col] += list(fibermap[col][inds_fibermap])

                # now store the fibermap; append additional ones as they come along.
                if fibermap_all is None:
                    fibermap_all = Table(fibermap)
                else:
                    fibermap_all = vstack([fibermap_all, Table(fibermap)])
        # since our redshifts are the true redshifts, zerr = 0
        data4later['ZERR'] = np.zeros_like(data4later['Z'])
        # picca seems to want to this - i think 0 is okay?
        data4later['ZWARN'] = np.zeros_like(data4later['Z'])
        # now lets add the hpixel for later
        data4later['HPIXELNUM'] = np.array([hpix] * len(data4later['Z'])).astype(int)
        # now lets drop duplicate targetids and sort by the targetids so things match between the ztrue
        # data and the fibermap (checked in desi_zcatalogs)
        # ztrue data
        _, inds_to_keep = np.unique(data4later['TARGETID'], return_index=True)
        logging.info(f' -- dropping duplicate targetids in ztrue data: going from {len(data4later["TARGETID"])} rows to {len(inds_to_keep)}')
        data4later = Table(data4later)[inds_to_keep].group_by('TARGETID')
        # stacked fibermaps
        _, inds_to_keep = np.unique(fibermap_all['TARGETID'].data, return_index=True)
        logging.info(f' -- dropping duplicate targetids in fibermap: going from {len(fibermap_all["TARGETID"].data)} rows to {len(inds_to_keep)}')
        fibermap_all = fibermap_all[inds_to_keep].group_by('TARGETID')

        # write this data out
        fname = f'{coadds_path}/spectra-{nside}/{hpix_mod}/{hpix}/ztrue.fits'
        fitsio.write(fname, data4later.as_array(),  extname='ZBEST', )  # table name is named so that desi_zcatalogs can recognize it
        #fitsio.write(fname, fibermap_exp,  extname='FIBERMAP_EXP', )
        fitsio.write(fname, fibermap_all.as_array(),  extname='FIBERMAP', )
        logging.info(f' -- wrote out {fname}')
        logging.info(f' ---------\n')

logging.info(' ## all done')
logging.info(f' ## time taken: {(time.time() - start0)/60: .2f} min')
