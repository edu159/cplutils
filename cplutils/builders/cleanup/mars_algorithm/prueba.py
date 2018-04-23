#!/usr/bin/env python
import numpy as np
from ase.lattice.hexagonal import *
from ase.lattice.compounds import *
import ase.io as io
from ase import Atoms, Atom

index1=1
index2=1
index3=1
mya = 5.038
myb = 5.038
myc = 13.772
myalpha = 90
mybeta = 90
mygamma = 120
gra = HEX_Fe2O3(symbol = ('Fe', 'O'),
# gra = HEX_Fe2O3(directions=[[0,0,0,1],[0,1,-1,0],[-2,1,1,0]],symbol = ('Fe', 'O'),
latticeconstant={'a':mya,'b':myb, 'c':myc,
'alpha':myalpha,
'beta':mybeta,
'gamma':mygamma},
size=(index1,index2,index3))

xmin, ymin, zmin = np.min(gra.arrays['positions'], axis=0)
gra.arrays['positions'][:,0] -= xmin
gra.arrays['positions'][:,1] -= ymin
gra.arrays['positions'][:,2] -= zmin
xmax, ymax, zmax = np.max(gra.arrays['positions'], axis=0)
print xmax, ymax, zmax
io.write('hexaFe2O3.xyz', gra, format='xyz')
