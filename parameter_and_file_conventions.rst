=====================
1. Telescope aperture
=====================
Design parameter menu
---------------------
With the exception of symmetry, all parameter keys are nested under the ``'Pupil'`` key category. See https://github.com/neilzim/SCDA/blob/master/scda_demo.ipynb for usage when initializing a design survey.

- Symmetry (defined by coronagraph subclass): ``'full'``, ``'half'``, ``'quart'``

- Primary mirror configuration (key ``'pm'``): ``'hex1'``, ``'hex2'``, ``'hex3'``, ``'hex4'``, ``'key24'``, ``'pie12'``, ``'pie08'``

- Secondary support strut configuration (key ``'ss'``): ``'y60'``, ``'y60off'``, ``'x'``, ``'cross'``, ``'t'``, ``'y90'``

- Secondary support strut thickness, either 2.5 or 10 cm (key ``'sst'``): ``'025'``, ``'100'``

- Secondary mirror obscuration present? (key ``'sm'``): ``True``, ``False``

File name format
----------------
Format spec: ``'TelAp_{0:s}_{1:s}{2:s}t{3:s}sm{4:d}_N{5:04d}.dat'``

0. Symmetry, string
1. Primary mirror, string
2. Support strut config, string
3. Support stut thickness, string
4. Secondary mirror obscuration, single-digit Boolean integer
5. Pupil array quadrant width, zero-padded 4-digit integer

Example: ``'TelAp_half_hex3crosst025so1_N0125.dat'``

===================
2. Focal plane mask
===================
Design parameter menu
---------------------
For now, only the occulting spot is supported. All keys are nested under the 'FPM' category.

- Spot radius (key ``rad''): float

- Mask array quadrant width (key ``M``): integer

File name format
----------------
FPM files are agnostic to physical radius.

Format spec: ``'FPM_occspot_M{0:03d}.dat'``

0. Mask array quadrant width, zero-padded 3-digit integer

Example: ``'FPM_occspot_M050.dat'``

=============
3. Lyot stop
=============
Design parameter menu
---------------------
For now, only an annular stop is supported, with and without secondary support struts.

- Symmetry (defined by coronagraph subclass): ``'full'``, ``'half'``, ``'quart'``

- Inner diameter of annular aperture, percentage of telescope aperture diameter (key ``'id'``, integer)

- Outer diameter of annular aperture, percentage of telescope aperture diameter (key ``'od'``, integer)

- Secondary support strut configuration, if telescope features are mimicked by the stop (key ``'ss'``, string, defaults to the telescope aperture secondary support strut configuration). If the stop is a simple annulus with no struts, the string is '0'.

- Secondary support strut thickness, if telescope features are mimicked by the stop (key ``'st'``, string, defaults to the telescope aperture secondary support strut thickness). If the stop is a simple annulus with no struts, the string is '0'.

- Padding of obscuration features, if present (key ``'pad'``, integer). The padding parameter is specified in units of telescope aperture diameter percentage. Padding is applied in an omindirectial sense by a shift-and-combine-and-mask routine, so it increases thickness on all sides of a given obscuration feature, and the thickness of all features increases by the same absolute propportion of the telescope aperture. This parameter is zero if obscuration features are not mimicked by the Lyot stop, or if they are mimicked but not padded.

- Alignment tolerance, percentage of telescope aperture diameter (key ``'altol'``, integer). Does not affect the Lyot stop file name, but  modifies the AMPL optimization program. Not yet supported.

File name format
----------------
Format spec: ``'LS_{0:s}_ann{1:02d}D{2:02d}_{3:s}t{4:s}p{5:02d}_N{6:04d}.dat'``

0. Symmetry, string
1. Inner diameter, zero-padded 2-digit integer
2. Outer diameter, zero-padded 2-digit integer
3. Secondary support stut config, string
4. Secondary support stut thickness in telescope aperture, string
5. Obscuration padding, zero-padded 2-digit integer
6. Pupil array quadrant width, zero-padded 4-digit integer
