from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import add_adsorbate
from ase.lattice.cubic import FaceCenteredCubic
from ase.io import read, write
from ase.visualize import view
import numpy as np
import time
from joblib import Parallel, delayed
from joblib.pool import has_shareable_memory
from ase.data import reference_states, atomic_numbers
Z = atomic_numbers['Cu']
print reference_states[Z]['a']
slab = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                         size=(10,2,10), symbol='Cu', pbc=(1,0,1))
write("surface.xyz", slab)
a = []
xmax, ymax, zmax = np.amax(slab.arrays['positions'], axis=0)
print xmax, ymax, zmax
#molecule = read("squalane.pdb")

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
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
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
roughness,X,Z = periodic1Dsin_surface(coeffs, ymax, xmax, zmax)
#roughness = (roughness + roughness0)/2
# Plot the surface.
surf = ax.plot_surface(X, Z, roughness, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(0.0, ymax)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
def f(k):
    return k
t0 = time.clock()
#keep_list = Parallel(n_jobs=8, max_nbytes=1e9)(delayed(f)(atom.index) for atom in has_shareable_memory(slab)) #if keep_atom(roughness, atom.position, (xmax, ymax, zmax)))
print "Time parallel: ", time.clock() - t0
t0 = time.clock()
keep_list = [atom.index for atom in slab if keep_atom(roughness, atom.position, (xmax, ymax, zmax))]
print "Time serial: ", time.clock() - t0
slab = slab[keep_list]
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

write('slab_test.xyz', slab)

#print('Adsorption energy:', e_slab + e_N2 - slab.get_potential_energy())
