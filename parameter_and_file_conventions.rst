v0
--

Telescope aperture
==================

Parameter menu
---------------

Symmetry (defined by coronagraph subclass): ``'full'``, ``'half'``, ``'quart'``

Primary mirror configuration (key ``'pm'``): ``'hex1'``, ``'hex2'``, ``'hex3'``, ``'hex4'``, ``'key24'``, ``'pie12'``, ``'pie8'``

Support strut configuration (key ``'ss'``): ``'y60'``, ``'y60off'``, ``'x'``, ``'cross'``, ``'t'``, ``'y90'``

Support strut thickness, either 2.5 or 10 cm (key ``'sst'``): ``'025'``, ``'100'``

Secondary obscuration present? (key ``'so'``): ``True``, ``False``

File name format
----------------

Example: ``'TelAp_half_hex3crosst025so1_N0125.dat'``

``'TelAp_{0:s}_{1:s}{2:s}t{3:s}so{4:d}_N{5:04d}.dat'``

0 - Symmetry, string
1 - Primary mirror, string
2 - Support strut config, string
3 - Support stut thickness, string
4 - Secondary obscuration, single digit integer (Boolean)
5 - Pupil quadrant width, zero-padded 4-digit integer

Focal plane mask
----------------



Lyot stop
----------
