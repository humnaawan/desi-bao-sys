# ------------------------------------------------------------------------
# this script takes in all exposures and creates a list of all the QSOs;
# needed for per-healpix-pixel grouping and coaddition.
# ------------------------------------------------------------------------
import os
import time
import datetime
import sys
import logging
from desimodel import footprint as desifootprint
from astropy.io import fits
from astropy.table import Table, join
import pandas as pd
import numpy as np

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--data-path',
                  dest='data_path', type='str',
                  help='path to directory with all the data \
                        (i.e., where folder name = night)')
parser.add_option('--outdir',
                  dest='outdir', type='str',
                  help='path to the directory where the output csv file is to saved.')
parser.add_option('--nside',
                  dest='nside', type=int,
                  help='nside resolution param')
(options, args) = parser.parse_args()
print('\n## inputs: %s' % options)
data_path = options.data_path
outdir = options.outdir
nside = options.nside
# ------------------------------------------------------------------------
start0 = time.time()

# set up the logger
temp = f'get-qso-list_{datetime.datetime.now()}'.replace(' ', '_').split('.')[0]
logging.basicConfig(filename=f'{outdir}/log_{temp}.log',
                    level=logging.DEBUG, filemode='w',
                    #format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p'
                   )
logging.info(f'running {sys.argv[0]}')
logging.info('\n## inputs: %s\n' % options)

explist = []

for night in [f for f in os.listdir(data_path) if os.path.isdir(f'{data_path}/{f}')]:
    logging.info(f'## looking at exposures for night {night}')
    for expid in [f for f in os.listdir(f'{data_path}/{night}') if os.path.isdir(f'{data_path}/{night}/{f}')]:
        logging.info(f'## looking at data for expid {expid}')
        subdir = f'{data_path}/{night}/{expid}'
        fnames = [f for f in os.listdir(subdir) if f.startswith('simspec') and f.endswith('fits')]
        if len(fnames) > 0: # doing this since not all exposure folders have data and we dont need to break because of that
            if len(fnames) > 1:
                raise ValueError(f'somethings weird: have {len(fnames)} simspec fits files in {subdir}')
            else:
                # i.e. have exactly only simspec file
                data_simspec = Table.read(f'{subdir}/{fnames[0]}', hdu='TRUTH', format='fits')
                # keep only the QSOs
                data_simspec = data_simspec[data_simspec['TRUESPECTYPE'] == 'QSO']

            for petal in range(10):
                logging.info(f'## dealing with petal {petal}')
                fname = f'{subdir}/spectra-{petal}-{expid}.fits'
                if os.path.exists(fname):
                    data_spectra = Table.read(fname, hdu='FIBERMAP', format='fits')

                    joined = join(data_spectra, data_simspec, keys='TARGETID')

                    # lets pull the tileid
                    tileid = np.unique(joined['TILEID'].value)
                    if len(tileid) != 1:
                        raise ValueError('dont have the functionality to handle multiple TILEIDs')
                    tileid = tileid[0]

                    # lets check the petal number
                    petal_ = np.unique(joined['PETAL_LOC'].value)
                    if len(petal_) != 1:
                        raise ValueError('should not have multiple petals here')
                    if petal != petal_:
                        raise ValueError(f'somethings wrong with the petal number here: {petal}, {petal_}')

                    # lets get the HEALPix number and store it
                    for hppixel in np.unique(desifootprint.radec2pix(nside=nside,
                                                                           ra=joined['TARGET_RA'].value,
                                                                           dec=joined['TARGET_DEC'].value
                                                                    )
                                                #hp.ang2pix(nside=nside,
                                                #        theta=np.pi / 2 - np.radians(joined['TARGET_DEC'].value),
                                                #        phi=np.radians(joined['TARGET_RA'].value))
                                            ):
                        explist.append([night, expid, tileid, petal, hppixel])

# covert the list to dataframe
explist = pd.DataFrame(explist, columns=['NIGHT', 'EXPID', 'TILEID', 'SPECTRO', 'HEALPIX'])
# save the data
fname = f'qso-exposures-list_nside{nside}.csv'
explist.to_csv(f'{outdir}/{fname}', index=False)

logging.info(f'## saved {fname} in {data_path}')
logging.info('## all done')
logging.info(f'## time taken: {(time.time() - start0)/60: .2f} min')