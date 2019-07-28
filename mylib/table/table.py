"""
Convenient tools of crossmatching astronomical tables.
"""

import os, sys
import numpy as np
from astropy.table import Table, vstack, hstack
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.constants as const

class astrotable(object):
    def __init__(self):
        self.masks = {}

    def _join_best(self, ext_table, radius):
        coords1 = self.coords
        coords2 = ext_table.coords

        idx, d2d, d3d = coords1.match_to_catalog_sky(coords2)
        idx_masked = idx[d2d.arcsecond<radius]

        tbl1 = self.table
        tbl2 = ext_table.table

        tbl2_matched = tbl2[:][idx]

        newtbl = hstack([tbl1, tbl2_matched])
        return newtbl

    def _join_all(self, ext_table, radius):
        do_something = 1

    def aread(self, filename, ra_col='ra', dec_col='dec', coord_unit='deg',\
              *args, **kwargs):
        tbl = Table.read(filename, *args, **kwargs)
        ra = tbl[ra_col]
        dec = tbl[dec_col]

        coords = SkyCoord(ra=ra, dec=dec, unit=coord_unit)

        self.table = tbl.copy()
        self.coords = coords

    def awrite(self, filename, *args, **kwargs):
        self.table.write(filename, *args, **kwargs)

    def join(self, ext_table, radius=2, type='best'):
        do_something = 1

    def add_mask(self, mask_name, mask_value):
        self.masks[mask_name] = mask_value
