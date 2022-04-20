import os
import time
import datetime
import sys
import logging

from desispec.scripts import group_spectra, coadd_spectra

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--data-path',
                  dest='data_path', type='str',
                  help='path to directory with all the data \
                        (i.e., where folder name = night)')
parser.add_option('--night',
                  dest='night', type='str',
                  help='night to consider')
parser.add_option('--expid',
                  dest='expid', type='str',
                  help='exposure ID to consider')
(options, args) = parser.parse_args()
print('\n## inputs: %s' % options)
data_path = options.data_path
night = options.night
expid = options.expid
# ------------------------------------------------------------------------
start0 = time.time()

path = f'{data_path}/{night}/{expid}/'

# set up the logger
temp = f'{datetime.datetime.now()}'.replace(' ', '_').split('.')[0]
logging.basicConfig(filename=f'{path}/log_{temp}.log',
                    level=logging.DEBUG, filemode='w', 
                    #format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p'
                   )
logging.info(f'running {sys.argv[0]}')
logging.info('\n## inputs: %s' % options)


nside = 64
npetals = 1 #10

# loop over all spectrographs
for npetal in range(npetals):
    # group things - pass all cframes of a given spectrograph
    cframes = [f'{path}/{f}' for f in os.listdir(path) if f.__contains__(f'{npetal}-{expid}.fits') and f.startswith('cframe')]           
    outfile = f'{path}/grouped-{npetal}-{expid}.fits'
    options = [f'--nside={nside}',
               f'--outfile={outfile}',
               '--inframes'] + cframes
    print(f'\nOPTIONS:  {options}\n\n')
    args = group_spectra.parse(options)
    group_spectra.main(args)
    
    # now coadd things
    infile = outfile
    outfile = f'{path}/coadd-{npetal}-{expid}.fits'
    options = ['--infile', infile,
               '--outfile', outfile]

    args = coadd_spectra.parse(options)
    coadd_spectra.main(args)

logging.info('\n## all done')
logging.info(f'## time taken: {(time.time() - start0)/60: .2f} min')