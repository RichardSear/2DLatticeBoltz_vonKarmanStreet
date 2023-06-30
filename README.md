# 2DLatticeBoltz_vonKarmanStreet
Lattice Boltzmann code to simulate flow past disc in 2D - shows von Karman streets of vorticities at Re around 100

This is just an edited version of a code from Palabos group at Uni Geneve (https://palabos.unige.ch/get-started/lattice-boltzmann/lattice-boltzmann-sample-codes-various-other-programming-languages/
). Note that this code does not quite implement Zou-He boundary conditions at left hand edge (to fix flow speed) and this code does (I hope!) but that does not seem to make any real difference.

It has a logical variable to either throw images to screen as simulation runs, or to save set of image (in specified subdirectory that you will need to create), which you can make into a movie. Images show both vorticity (as red/blue) and flow speed as arrow plot.

Code varies Reynolds number Re by varying LB's kinematic viscosity. von Karman street for Re around 100 or 200. For Re small get relatively viscous flow around disc. LB code has limitations and system is small, so will not work at Re of 1000s (except perhaps with larger system and larger discs) and finite size effects will influence results. 
