"""

http://docs.astropy.org/en/v1.3.3//coordinates/skycoord.html#example-2-plotting-star-positions-in-bulge-and-disk

http://docs.astropy.org/en/stable/coordinates/skycoord.html#example-2-plotting-star-positions-in-bulge-and-disk



"""
import sys
import time
import numpy as np

import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import GeocentricTrueEcliptic
from astropy.coordinates import BarycentricTrueEcliptic
from astropy.coordinates import HeliocentricTrueEcliptic
from astropy.coordinates import Angle
from astropy.coordinates import Longitude, Latitude

import corner

def mk_radec(ndata=10):

    ndata =  int(ndata)
    ra = np.random.uniform(0.0, 360.0, ndata)
    cosdec = np.random.uniform(-1.0, 1.0, ndata)
    dec = np.rad2deg(np.arccos(cosdec)) - 90.0

    return ra, dec

ra, dec = mk_radec()

coords_radec = SkyCoord(ra, dec, unit='deg')

print(coords_radec)

print(coords_radec.galactic)

t0 = time.time()
coords_galactic = coords_radec.transform_to(Galactic)
print(coords_galactic)
print(time.time() - t0)
print(coords_galactic.l)
print(coords_galactic.b)


t0 = time.time()
coords_ecliptic = coords_radec.transform_to(BarycentricTrueEcliptic)
print(coords_ecliptic)
print(time.time() - t0)
print(coords_ecliptic.lon)
print(coords_ecliptic.lat)


def mk_random_galaxy(size=(1000, 1000)):

    disk = np.random.multivariate_normal(mean=[0, 0, 0],
                                         cov=np.diag([1.0, 1.0, 0.1]),
                                         size=size[0])

    print(type(disk), len(disk), disk.shape, disk.size)

    bulge = np.random.multivariate_normal(mean=[0,0,0],
                                          cov=np.diag([1.0, 1.0, 1.0]),
                                          size=size[1])
    return disk, bulge

disk, bulge = mk_random_galaxy()

make_cornerplot = False
if make_cornerplot:
    figure = corner.corner(disk,
                       labels=[r"$u$", r"$v$", r"$w$"])
    figure.suptitle('corner.py: mk_random_galaxy: disk')

    plt.savefig('disk_corner.png')
    plt.show()

    figure = corner.corner(bulge,
                       labels=[r"$u$", r"$v$", r"$w$"])
    figure.suptitle('corner.py: mk_random_galaxy: bulge')

    plt.savefig('bulge_corner.png')
    plt.show()



disk, bulge = mk_random_galaxy(size=(5000, 2000))

galaxy = np.concatenate([disk, bulge]) * u.rad
disk = disk * u.rad
bulge = bulge * u.rad
print(len(galaxy), type(galaxy), galaxy.shape, galaxy.size)
print(np.min(galaxy), np.max(galaxy), np.median(galaxy))

# coordinates are transformed into an astropy.coordinates SkyCoord object.
coords_disk = SkyCoord(disk,
                     representation='cartesian',
                     frame='galactic')
print(coords_disk)
print(coords_disk.u)
# print(coords_disk.lon)

coords_disk.representation = 'spherical'
print(coords_disk)
print(coords_disk.l)
print(coords_disk.b)

coords_disk_ecliptic = coords_disk.transform_to(BarycentricTrueEcliptic)
print(coords_disk_ecliptic)
print(time.time() - t0)
print(coords_disk_ecliptic.lon)
print(coords_disk_ecliptic.lat)




coords_disk_icrs = coords_disk.icrs
print(coords_disk.icrs)
coords_disk_icrs = coords_disk.galactic
#SkyCoord_ICRS = SkyCoord.transform_to(ICRS)
coords_disk_icrs_ra = coords_disk.icrs.ra
coords_disk_icrs_dec = coords_disk.icrs.dec

# print(coords_disk_icrs.representation_component_names)
# print(ICRS().representation_info)
# print(Galactic().representation_info)


# help(coords_disk_icrs)

coords_disk_galactic = coords_disk.galactic
print(coords_disk_galactic.representation_component_names)

# help(coords_disk_galactic)

coords_bulge = SkyCoord(bulge,
                     representation='cartesian',
                     frame='galactic')
coords_bulge_icrs = coords_bulge.icrs
coords_bulge_galactic = coords_bulge.galactic


# convert coordinates to radians and make sure they are between -pi and +pi
ra_disk_rad = coords_disk.icrs.ra.wrap_at(180 * u.deg).radian
dec_disk_rad = coords_disk.icrs.dec.radian

ra_bulge_rad = coords_bulge_icrs.ra.wrap_at(180 * u.deg).radian
dec_bulge_rad = coords_bulge_icrs.dec.radian


plt.figure(figsize=(8, 4.5))
plt.subplot(111, projection="aitoff")
plt.suptitle("Aitoff projection of random data")

plt.grid(True)
plt.plot(ra_disk_rad, dec_disk_rad, 'ob',
         markersize=2, alpha=0.5,
         label='disk')
plt.plot(ra_bulge_rad, dec_bulge_rad, 'or',
         markersize=2, alpha=0.5,
         label='bulge')

plt.xlabel('Right Ascension (ra)')
plt.ylabel('Declination (dec)')
plt.legend()

plt.subplots_adjust(top=0.90, bottom=0.0)

plt.savefig('plot_galactic_radec.png')
plt.show()


# plot in Galactic coordinates
coords_disk.representation = 'spherical'
coords_bulge.representation = 'spherical'

# convert coordinates to radians and make sure they are between -pi and +pi
l_disk_rad = coords_disk.l.wrap_at(180 * u.deg).radian
b_disk_rad = coords_disk.b.radian

l_bulge_rad = coords_bulge.galactic.l.wrap_at(180 * u.deg).radian
b_bulge_rad = coords_bulge.galactic.b.radian

plt.figure(figsize=(8, 4.5))
plt.subplot(111, projection="aitoff")
plt.suptitle("Aitoff projection of our random data")

plt.grid(True)
plt.plot(l_disk_rad, b_disk_rad, 'ob', markersize=1, alpha=1.0)
plt.plot(l_bulge_rad, b_bulge_rad, 'or', markersize=1, alpha=1.0)

plt.xlabel('Galactic Longitude (l)')
plt.ylabel('Galactic Latitute (b)')
plt.legend()


plt.subplots_adjust(top=0.90, bottom=0.0)

plt.savefig('plot_galactic_galactic.png')
plt.show()

# now do Ecliptc
coords_disk_ecliptic = coords_disk.transform_to(BarycentricTrueEcliptic)
print(coords_disk_ecliptic)
print(time.time() - t0)
print(coords_disk_ecliptic.lon)
print(coords_disk_ecliptic.lat)
coords_bulge_ecliptic = coords_bulge.transform_to(BarycentricTrueEcliptic)

# plot in Galactic coordinates
# coords_disk.representation = 'spherical'
# coords_bulge.representation = 'spherical'

# convert coordinates to radians and make sure they are between -pi and +pi
l_disk_rad = coords_disk_ecliptic.lon.wrap_at(180 * u.deg).radian
b_disk_rad = coords_disk_ecliptic.lat.radian

l_bulge_rad = coords_bulge_ecliptic.lon.wrap_at(180 * u.deg).radian
b_bulge_rad = coords_bulge_ecliptic.lat.radian

plt.figure(figsize=(8, 4.5))
plt.subplot(111, projection="aitoff")
plt.suptitle("Aitoff projection of our random data")

plt.grid(True)
plt.plot(l_disk_rad, b_disk_rad, 'ob', markersize=1, alpha=1.0)
plt.plot(l_bulge_rad, b_bulge_rad, 'or', markersize=1, alpha=1.0)

plt.xlabel('Ecliptic Longitude (lon)')
plt.ylabel('Ecliptic Latitute (lat)')
plt.legend()


plt.subplots_adjust(top=0.90, bottom=0.0)

plt.savefig('plot_galactic_ecliptic.png')
plt.show()
