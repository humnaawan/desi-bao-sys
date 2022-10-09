# ------------------------------------------------------------------------
# this script takes the zcatalog (zcat.fits) and creates the drq catalog
# needed for picca. following https://desi.lbl.gov/trac/wiki/LymanAlphaWG/how_to_run_picca
# ------------------------------------------------------------------------
from picca import converters

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--zcat-path',
                  dest='zcat_path', type='str',
                  help='path to zcat.fits file.')
parser.add_option('--outdir',
                  dest='outdir', type='str',
                  help='path to directory where zcat_drq should be saved.')
(options, args) = parser.parse_args()
zcat_path = options.zcat_path
out_path = f'{options.outdir}/zcat_drq.fits'

converters.desi_from_ztarget_to_drq(in_path=zcat_path, out_path=out_path,
                                    spec_type='QSO',
                                    downsampling_z_cut=2.1,
                                    downsampling_num=100000
                                    )
