from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
"""

Example table i/o for votable with timing comparison

"""

import os
import sys
import time

import numpy as np

t0 = time.time()
import astropy
print('Elapsed time(secs):', time.time() - t0)
print()
print(astropy.__version__)

from astropy.table import Table
from astropy.io.votable import from_table, writeto

nrows = int(1e6)

t0 = time.time()
print('Create data:',  nrows, 'rows')
col0 = np.linspace(1, nrows, num=nrows, dtype=np.int32)
col1 = np.linspace(1, nrows, num=nrows, dtype='float32')
col2 = np.linspace(1, nrows, num=nrows, dtype='float32')
print('Elapsed time(secs):', time.time() - t0)
print()
print(col0[0], col0[-1])
print(col1[0], col1[-1])

# table = Table([col1, col2, col3])

table = Table([col0])
table = Table([col0, col1, col2])

print('Elapsed time(secs):', time.time() - t0)
print()
table.info()
print('Elapsed time(secs):', time.time() - t0)
print()
# table['col1'].unit = 'deg'

t0 = time.time()
outfile = 'table.fits'
print('Write:', outfile)
table.write(outfile, overwrite=True)
print('Elapsed time(secs):', time.time() - t0)
print()


t0 = time.time()
infile = outfile
print('Read:', infile)
input = Table.read(infile)
print('Elapsed time(secs):', time.time() - t0)
input.info()
print()


t0 = time.time()
outfile = 'table_binary.vot'
print('Write:', outfile)
table.write(outfile, table_id='example_table',
            format='votable',  tabledata_format='binary',
            overwrite=True)
print('Elapsed time(secs):', time.time() - t0)
print()

t0 = time.time()
infile = outfile
print('Read:', infile)
input = Table.read(infile)
print('Elapsed time(secs):', time.time() - t0)
print('Number of rows:', len(input))
input.info()
print()


t0 = time.time()
outfile = 'table_binary2.vot'
print('Write:', outfile)
table.write(outfile, table_id='example_table',
            format='votable',  tabledata_format='binary2',
            overwrite=True)
print('Elapsed time(secs):', time.time() - t0)
print()

t0 = time.time()
infile = outfile
print('Read:', infile)
input = Table.read(infile)
print('Elapsed time(secs):', time.time() - t0)
print('Number of rows:', len(input))
input.info()
print()


#
t0 = time.time()
outfile = 'table_ascii.vot'
print('Write:', outfile)
table.write(outfile, table_id='example_table',
            format='votable',
            overwrite=True)
print('Elapsed time(secs):', time.time() - t0)
print()

t0 = time.time()
infile = outfile
print('Read:', infile)
input = Table.read(infile)
print('Elapsed time(secs):', time.time() - t0)
print('Number of rows:', len(input))
input.info()
print()

# hdf5
t0 = time.time()
outfile = 'table.hdf5'
print('Write:', outfile)
table.write(outfile, path='data',
            overwrite=True)
print('Elapsed time(secs):', time.time() - t0)
print()

t0 = time.time()
infile = outfile
print('Read:', infile)
input = Table.read(infile, path='data')
print('Elapsed time(secs):', time.time() - t0)
print('Number of rows:', len(input))
input.info()
print()

#hdf with compression
t0 = time.time()
outfile = 'table_compressed.hdf5'
print('Write:', outfile)
table.write(outfile, path='data',
            compression=True,
            overwrite=True)
print('Elapsed time(secs):', time.time() - t0)
print()

t0 = time.time()
infile = outfile
print('Read:', infile)
input = Table.read(infile, path='data')
print('Elapsed time(secs):', time.time() - t0)
print('Number of rows:', len(input))
print()


# csv
t0 = time.time()
outfile = 'table.csv'
print('Write:', outfile)
table.write(outfile,
            overwrite=True)
print('Elapsed time(secs):', time.time() - t0)
print()

t0 = time.time()
infile = outfile
print('Read:', infile)
input = Table.read(infile)
print('Elapsed time(secs):', time.time() - t0)
print('Number of rows:', len(input))
input.info()
print()
