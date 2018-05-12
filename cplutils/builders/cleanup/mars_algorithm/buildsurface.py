from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import add_adsorbate
from ase.lattice.cubic import FaceCenteredCubic
from ase.io import read, write
from ase.visualize import view
import numpy as np
import time
# from joblib import Parallel, delayed
# from joblib.pool import has_shareable_memory

from ase.data import reference_states, atomic_numbers
from ase.lattice.hexagonal import *
from ase.lattice.compounds import *
import ase.io as io
from ase import Atoms, Atom

if False:
    ## HCP IRON OXIDE
    index1=4
    index2=4
    index3=4
    mya = 5.038
    myb = 5.038
    myc = 13.772
    myalpha = 90
    mybeta = 90
    mygamma = 120

    # slab = HEX_Fe2O3(directions=[[0,0,0,1],[0,-1,1,0],[2,-1,-1,0]]  ,symbol = ('Fe', 'O'),
    slab = HEX_Fe2O3(directions=[[0,0,0,1],[0,1,-1,0],[-2,1,1,0]],symbol = ('Fe', 'O'),
    latticeconstant={'a':mya,'b':myb, 'c':myc,
    'alpha':myalpha,
    'beta':mybeta,
    'gamma':mygamma},
    size=(index1,index2,index3), pbc=(1,1,1))
    # io.write('hexaFe2O3.xyz', gra, format='xyz')
    ### HCP IRON OXIDE
    # write("surface.xyz", slab)

Z = atomic_numbers['Cu']
print reference_states[Z]['a']
slab = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                         size=(100,20,100), symbol='Cu', pbc=(1,0,1))


a = []
xmin, ymin, zmin = np.min(slab.arrays['positions'], axis=0)
print xmin, ymin, zmin
slab.arrays['positions'][:,0] -= xmin
slab.arrays['positions'][:,1] -= ymin
slab.arrays['positions'][:,2] -= zmin
xmax, ymax, zmax = np.max(slab.arrays['positions'], axis=0)
print xmax, ymax, zmax



## READ ROUGHNESS
fname = "./heightmap.dat"
with open(fname, 'r') as f:
   header = [next(f) for x in xrange(4)]

Nx = int(header[0].replace("\n","").strip())
Ny = int(header[1].replace("\n","").strip())
Lx = float(header[2].replace("\n","").strip())
Ly = float(header[3].replace("\n","").strip())
print Nx, Ny

data = np.genfromtxt(fname, skip_header=4)

X = data[:,0].reshape(Nx, Ny)*xmax
Y = data[:,1].reshape(Nx, Ny)*zmax
heightmap = data[:,2].reshape(Nx, Ny) 
max_y, min_y= np.max(heightmap), np.min(heightmap)
heightmap = (heightmap - min_y) / (max_y - min_y) * ymax
# READ ROUGHNESS END

# molecule = read("squalane.pdb")

def keep_atom(roughness, atom_pos, bounds):
    xpoints, zpoints = roughness.shape
    dx, dz = bounds[0] / (xpoints-1), bounds[2] / (zpoints-1)
    i, k = int(atom_pos[0] / dx), int(atom_pos[2] / dz)
#    print i, k, roughness[i, k], atom_pos[1]
    if roughness[i, k] >= atom_pos[1]:
        return True
    return False

#TODO MAKE it 1D roughness 
def spectral_surface(nrange, bounds, res=5., nx=600, ny=600):

    x = np.linspace(0.0, bounds[0], nx)
    y = np.linspace(0.0, bounds[2], ny)
    z = np.zeros([nx,ny])

    for level in range(10):
        rand = np.random.random(3)
        m = np.ceil(rand[0]*res/(0.1*level+1))
        n = np.ceil(rand[1]*res/(0.1*level+1))
        A = (rand[2]-0.5)/float(n*m)
        z += A*np.outer(np.sin(np.pi*n*x/bounds[0])+np.cos(np.pi*n*x/bounds[0]), np.sin(m*np.pi*y/bounds[2])+np.cos(m*np.pi*y/bounds[2]))

    #scale
    z = nrange*np.abs(z/np.max(z))

    return z
    
print "paso 0"

#slab.set_calculator(EMT())
#e_slab = slab.get_potential_energy()

#molecule.set_calculator(EMT())
#e_N2 = molecule.get_potential_energy()

import matplotlib.pyplot as plt
plt.pcolormesh(X, Y, heightmap, cmap=plt.cm.RdYlBu_r)
plt.axis("tight")
plt.colorbar()

keep_list = [atom.index for atom in slab if keep_atom(heightmap, atom.position, (xmax, ymax, zmax))]
slab = slab[keep_list]
write('slab_test.xyz', slab)
# plt.show()



from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
print "ymax ", ymax 
coeffs = [(1.0, 1.0), (0.1, 50.0)]
# Add resx and resz to resolve correctly both dimensions
def periodic1Dsin_surface(coeffs, scale, Lx, Lz, resx=1000, resz=1000):
    height = 0.0
    X = np.linspace(0.0, Lx, resx)
    Z = np.linspace(0.0, Lz, resz)
    X, Z = np.meshgrid(X, Z)
    for amp, l in coeffs:
        height += amp*np.sin(l*np.pi*Z/zmax)
    height = (scale)*abs(height/np.max(height))
    return (height, X, Z)


#roughness0 = spectral_surface(ymax, (xmax, ymax,zmax), res=5.0, nx=1000, ny=1000)
# roughness,X,Z = periodic1Dsin_surface(coeffs, ymax, xmax, zmax)
# print "H: ", len(roughness[0]), len(roughness[1])
#roughness = (roughness + roughness0)/2
# Plot the surface.

# PLOT
# surf = ax.plot_surface(X, Z, heightmap, cmap=cm.coolwarm,
#                                linewidth=0, antialiased=False)
#
# # Customize the z axis.
# ax.set_zlim(0.0, ymax)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))


# Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

# plt.show()
def f(k):
    return k
t0 = time.clock()
#keep_list = Parallel(n_jobs=8, max_nbytes=1e9)(delayed(f)(atom.index) for atom in has_shareable_memory(slab)) #if keep_atom(roughness, atom.position, (xmax, ymax, zmax)))
print "Time parallel: ", time.clock() - t0
t0 = time.clock()
print "Time serial: ", time.clock() - t0
print "paso 2"

overage=0.5
nmols= 100
offset = 5
#for m in xrange(nmols):
#    add_adsorbate(slab, molecule, 5, (offset*m,offset*m))
#constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
#slab.set_constraint(constraint)
#dyn = QuasiNewton(slab, trajectory='N2Cu.traj')
#dyn.run(fmax=0.05)
#view(slab, viewer="VMD")
#view(slab)


#print('Adsorption energy:', e_slab + e_N2 - slab.get_potential_energy())
