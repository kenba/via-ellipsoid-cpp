#!/usr/bin/env python

# Copyright (c) 2025 Ken Barker
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#  @file verify_GeodTest_file_initial_azimuths.py
#  @brief Reads and verifys the initial azimuths in the GeodTest.dat file at
#  https://geographiclib.sourceforge.io/C++/doc/geodesic.html#testgeod

import argparse
import polars as pl
from enum import Enum

from via_angle import Angle, Degrees
from via_sphere import calculate_gc_azimuth

# The columns of the GeodTest.dat file
class Column(Enum):
    latitude_1   = 0
    longitude_1  = 1
    azimuth_1    = 2
    latitude_2   = 3
    longitude_2  = 4
    azimuth_2    = 5
    distance_m   = 6
    distance_deg = 7
    m12          = 8
    area         = 9

# Read the first 6 columns of the GeodTest.dat file into a polar LazyFrame
def read_geotest_file():
    # Select all the test data entries
    pathname = 'https://sourceforge.net/projects/geographiclib/files/testdata/GeodTest.dat.gz/download'
    lf = pl.scan_csv(pathname, separator=' ', has_header=False).select(
        ['column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6']
    ).collect()

    return lf

def verify_geotest_initial_azimuths(output_values):
    lf = read_geotest_file()

    suspicious_southbound_geodesics = 0

    if output_values:
        print("line,lat_1,lat_2,lon_2,azi_1,azi_2,gc_azi")

    i = 0
    for row in lf.rows():
        # Get departure and arrival latitudes and departure azimuth from lf
        lat_1 = row[Column.latitude_1.value]
        azi_1 = row[Column.azimuth_1.value]
        lat_2 = row[Column.latitude_2.value]
        lon_2 = row[Column.longitude_2.value]
        azi_2 = row[Column.azimuth_2.value]

        # A valid non-antipodal southbound_equator_crossing should have an azimuth
        # in the same quadrant as a great circle arc between the positions
        delta_lat_abs = abs(lat_1 + lat_2)
        southbound_equator_crossing = (lat_1 > 0.0) and (lat_2 < 0.0)
        if southbound_equator_crossing and lon_2 < 179.0 and delta_lat_abs < 0.01 and azi_1 < 90.0:
            # Test whether in same quadrant as great circle
            gc_azi = calculate_gc_azimuth(Angle(Degrees(lat_1)), Angle(Degrees(lat_2)), Angle(Degrees(lon_2))).to_degrees().v()
            if gc_azi > 90.0:
                suspicious_southbound_geodesics += 1
                if output_values:
                    print(f"{i},{lat_1},{lat_2},{lon_2},{azi_1},{azi_2},{gc_azi}")

        i += 1

    if not output_values:
        print(f"total geodesics: {i}")
        print(f"suspicious_southbound_geodesics: {suspicious_southbound_geodesics}")

if __name__ == '__main__':
    # Initialize parser
    parser = argparse.ArgumentParser(description = 'Verfiy GeodTest.dat initial azimuths.')
    parser.add_argument("-l", "--Lines", help = "show suspicious lines, default False")
    args = parser.parse_args()

    output_values = args.Lines
    verify_geotest_initial_azimuths(output_values)
