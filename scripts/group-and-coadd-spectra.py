# ------------------------------------------------------------------------
# this script takes in the qso exposure list and uses it to group and coadd
# spectra in each healpix pixel.
# ------------------------------------------------------------------------
import os
import time
import datetime
import sys
import logging
import subprocess
import pandas as pd
import numpy as np

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--exposures-path',
                  dest='exposures_path', type='str',
                  help='path to directory with all exposures folder with night sub-folders.')
parser.add_option('--exposures-list-path',
                  dest='exposures_list_path', type='str',
                  help='path to csv file with the exposures and their corresponding healpix numbers.')
parser.add_option('--outdir',
                  dest='outdir', type='str',
                  help='path to the directory with the spectra-{nside} folder would sit.')
parser.add_option('--nside',
                  dest='nside', type=int,
                  help='HEALPix resolution param')
parser.add_option('--debug',
                  action='store_true', dest='debug', default=False,
                  help='run things in debug mode => some consider only one HEALPix pixel.')
(options, args) = parser.parse_args()
print('\n## inputs: %s' % options)
exposures_path = options.exposures_path
exposures_list_path = options.exposures_list_path
main_outdir = options.outdir
nside = options.nside
debug = options.debug
# ------------------------------------------------------------------------
start0 = time.time()

# set up the logger
temp = f'group+coadd_{datetime.datetime.now()}'.replace(' ', '_').split('.')[0]
logging.basicConfig(filename=f'{main_outdir}/log_{temp}.log',
                    level=logging.DEBUG, filemode='w', 
                    #format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p'
                   )
logging.info(f'running {sys.argv[0]}')
logging.info(f'\n## inputs: {options}\n')

subdir = f'{main_outdir}/spectra-{nside}'
if not os.path.exists(subdir):
    os.makedirs(subdir, exist_ok=True)
    print(f'## created the directory {subdir}')

df = pd.read_csv(exposures_list_path)
pix_list = np.unique(df['HEALPIX'])
if debug: pix_list = pix_list[0:1]

for hp_pixel in pix_list:
    logging.info(f'## running things for pixelnum = {hp_pixel}')
    subdir = f'{main_outdir}/spectra-{nside}/{hp_pixel//100}'
    if not os.path.exists(subdir):
        if not os.path.exists(subdir):
            os.makedirs(subdir, exist_ok=True)
            logging.info(f'## created the directory {subdir}')
    subdir = f'{main_outdir}/spectra-{nside}/{hp_pixel//100}/{hp_pixel}'
    if not os.path.exists(subdir):
        if not os.path.exists(subdir):
            os.makedirs(subdir, exist_ok=True)
            logging.info(f'## created the directory {subdir}')

    bash_command = 'desi_group_spectra ' + \
                    f'--expfile {exposures_list_path} ' + \
                    f'--outfile {subdir}/grouped-{nside}-{hp_pixel}.fits ' + \
                    f'--reduxdir {exposures_path} ' + \
                    f'--nside {nside} --healpix {hp_pixel}'
    logging.info(f'running bash_command = {bash_command}')
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    if error is not None:
        print(error)
    logging.info(f'command output: {output}')
    logging.info(f'command error: {error}\n')


    bash_command = 'desi_coadd_spectra ' + \
                    f'--infile {subdir}/grouped-{nside}-{hp_pixel}.fits ' + \
                    f'--outfile {subdir}/spectra-{nside}-{hp_pixel}.fits '
    logging.info(f'running bash_command = {bash_command}')
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    if error is not None:
        print(error)
    logging.info(f'command output: {output}')
    logging.info(f'command error: {error}\n')

logging.info('\n## all done')
logging.info(f'## time taken: {(time.time() - start0)/60: .2f} min')