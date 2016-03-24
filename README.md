# RG_test

This is a part of the code that I used in my thesis, to determine the absolute magnitude of small 
amplitude red giants. Details regarding the procedure can be found here:

https://obswww.unige.ch/~palavers/Thesis/thesis_col.pdf

(Chapter 3, and in particular 3.6 and 3.8).

Briefly, this particular example uses only OGLE III LMC small amplitude RG's. Code does the following:

1) Picks a "test" star.
2) Searches all the remaining stars that have the same combo of top three periods (within +-2.5% 
in logP)
3) Peak of the distribution of the selected stars (i.e. those with the same periods) is considered as
the indicator of the "true" absolute magnitude.

Black points in the left panel represent the entire sample. Yellow crosses represent the position of 
the test star in PL diagram. Transparent (R,G,B) bins indicate +-2.5% around the three periods of the
test star (logP).

Black points in the middle panel represent only the stars with the "three-period combo" matching that
of the test star. Yellow dashed line is the true mag of the test star, cyan dashed line is the 
inferred mag (Ideally they should be overlapping since the code is "reconstructing" LMC).

Right panel shows the difference between true and inferred mag (blue patch between yellow and cyan 
lines), and the normalized smoothed histogram and KDE of the magnitude distribution of the stars from the middle plot.

