=====================
1. Telescope aperture
=====================
Design parameter input menu
---------------------
With the exception of symmetry, all parameter keys are nested under the ``'Pupil'`` key category. See https://github.com/neilzim/SCDA/blob/master/scda_demo.ipynb for usage when initializing a design survey.

- Symmetry (defined by coronagraph subclass): ``'full'``, ``'half'``, ``'quart'``

- Primary mirror configuration (key ``'pm'``): ``'hex1'``, ``'hex2'``, ``'hex3'``, ``'hex4'``, ``'key24'``, ``'pie12'``, ``'pie08'``

- Secondary support strut configuration (key ``'ss'``): ``'Y60d'``, ``'Yoff60d'``, ``'X'``, ``'Cross'``, ``'T'``, ``'Y90d'``

- Secondary support strut thickness, either 2.5 or 10 cm (key ``'sst'``): ``'025'``, ``'100'``

- Secondary mirror obscuration present? (key ``'sm'``): ``True``, ``False``

File name format
----------------
Format spec: ``'TelAp_{0:s}_{1:s}{2:s}{3:s}sm{4:d}_N{5:04d}.dat'``

0. Symmetry, string
1. Primary mirror, string
2. Support strut config, string
3. Support stut thickness, string
4. Secondary mirror obscuration, single-digit Boolean integer
5. Pupil array quadrant width, zero-padded 4-digit integer

Example: ``'TelAp_half_hex3Cross025sm1_N0125.dat'``

===================
2. Focal plane mask
===================
Design parameter input menu
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
Design parameter input menu
---------------------
For now, only an annular stop is supported, with and without secondary support struts.

- Symmetry (defined by coronagraph subclass): ``'full'``, ``'half'``, ``'quart'``

- Inner diameter of annular aperture, percentage of telescope aperture diameter (key ``'id'``, integer, default ``90``)

- Outer diameter of annular aperture, percentage of telescope aperture diameter (key ``'od'``, integer, default ``20``)

- Obscuration switch (key 'obscure', integer, default ``0``): ``0``, ``1``, ``2``. If ``0``, the stop is a clear annulus. If ``1``, the stop mimics the secondary obscuration configuration (secondary mirror plus support struts) of the telescope aperture. Otherwise, the stop is a clear annulus. If ``2``, the stop mimics the complete re-imaged telescope pupil, including secondary obscuration and mirror gaps.

- Padding of secondary obscuration features, if present (key ``'spad'``, integer in the range 0 to 100, default ``0``). The padding parameter is specified as a percentage of telescope pupil diameter [1]_.

- Padding of primary mirror gap features, if present (key ``'ppad'``, integer, default ``0``). The padding parameter is specified as a percentage of telescope pupil diameter. Padding is applied in an omindirectial sense by a shift-and-combine-and-mask routine, so it increases thickness on all sides of a given obscuration feature, and the thickness of all features increases by the same absolute propportion of the pupil diameter. This parameter remains zero if gaps are not mimicked by the Lyot stop, or if they are mimicked but not padded.

- Alignment tolerance, percentage of telescope aperture diameter (key ``'altol'``, integer). Does not affect the Lyot stop file name, but  modifies the AMPL optimization program. Not yet supported.

File name format
----------------
Format spec

A. When telescope obscurations are ommitted and the stop is a clear annulus:

Format spec: ``'LS_{0:s}_ann{1:02d}D{2:02d}_clear_N{3:04d}.dat'``

0. Symmetry, string
1. Inner diameter, zero-padded 2-digit integer
2. Outer diameter, zero-padded 2-digit integer
3. Pupil array quadrant width, zero-padded 4-digit integer

Example: ``'LS_quart_ann15D80_clear_N0300.dat'``

B. When telescope obscurations are mimicked by the stop, the relevant design parameters from the telescope aperture are included:

Format spec: ``'LS_{0:s}_ann{1:02d}D{2:02d}_{3:s}{4:s}sm{5:d}p{6:02d}_N{7:04d}.dat'``

0. Symmetry, string
1. Inner diameter, zero-padded 2-digit integer
2. Outer diameter, zero-padded 2-digit integer
3. Support stut config, string
4. Support stut thickness in telescope aperture, string
5. Secondary mirror obscuration, single-digit Boolean integer
6. Obscuration padding, zero-padded 2-digit integer
7. Pupil array quadrant width, zero-padded 4-digit integer

Examples: ``'LS_quart_ann20D85_X100sm1p08_N0300.dat'``


..  [1] Padding is applied in an omindirectial sense by a shift-and-combine-and-mask routine, so it increases thickness on all sides of a given obscuration feature, and the thickness of all features increases by the same absolute propportion of the pupil diameter. This parameter remains zero if obscuration features are not mimicked by the Lyot stop, or if they are mimicked but not padded.
