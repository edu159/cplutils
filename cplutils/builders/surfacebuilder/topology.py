from ase.lattice.cubic import FaceCenteredCubic
from ase.io import write, read
import numpy as np
from ase.data import reference_states, atomic_numbers
import mars

def mask_atom(roughness, atom_pos, bounds):
    xpoints, zpoints = roughness.shape
    dx, dz = bounds[0] / (xpoints-1), bounds[2] / (zpoints-1)
    i, k = int(atom_pos[0] / dx), int(atom_pos[2] / dz)
    if roughness[i, k] >= atom_pos[1]:
        return True
    return False

def heightmap_generator_MARS(fname, slab):
    with open(fname, 'r') as f:
       header = [next(f) for x in xrange(4)]
    max_coords = np.max(slab.arrays['positions'], axis=0)
    Nx = int(header[0].replace("\n","").strip())
    Ny = int(header[1].replace("\n","").strip())
    Lx = float(header[2].replace("\n","").strip())
    Ly = float(header[3].replace("\n","").strip())
    data = np.genfromtxt(fname, skip_header=4)
    X = data[:,0].reshape(Nx, Ny)*max_coords[0]
    Z = data[:,1].reshape(Nx, Ny)*max_coords[2]
    heightmap = data[:,2].reshape(Nx, Ny) 
    # Scale heightmap
    max_y, min_y= np.max(heightmap), np.min(heightmap)
    heightmap = (heightmap - min_y) / (max_y - min_y) * max_coords[1]
    return X, Z, heightmap

# def heightmap_generator_sine1D(amp, periods, phase, slab):
def heightmap_generator_SINE2D(amp, period_x, shift_x, period_z, shift_z, slab):
    Ncx, Ncy, Ncz = slab.ncells
    max_coords = np.max(slab.arrays['positions'], axis=0)
    # assert Ncx % period_x == 0 and  Ncz % period_z  == 0
    assert amp <= 1.0 and amp >= 0.0
    ppc = 10 # points per cell
    X = np.linspace(0.0, max_coords[0], Ncx*ppc)
    Z = np.linspace(0.0, max_coords[2], Ncz*ppc) 
    X_grid, Z_grid = np.meshgrid(X,Z)
    Lx, Lz = max_coords[0], max_coords[2]
    dx, dz = Lx / Ncx, Lz / Ncz
    period_x_len = period_x * dx
    period_z_len = period_z * dz
    heightmap = (0.5*(np.sin(2*np.pi/period_x_len*X_grid)) + (0.5*np.sin(2*np.pi/period_z_len*Z_grid)) + 1)/2.0 # * max_coords[1]
    heightmap *= max_coords[1] * amp 
    heightmap += max_coords[1] * (1.0 - amp)
    heightmap = np.roll(heightmap, shift_x*ppc, axis=1)
    heightmap = np.roll(heightmap, shift_z*ppc, axis=0)
    return X_grid, Z_grid, heightmap


# def heightmap_generator_sine1D(amp, periods, phase, max_coords):
def heightmap_generator_BLOCKS(amp, ncx_block, ncx_gap, shift_x, ncz_block, ncz_gap, shift_z, invert, slab):
    Ncx, Ncy, Ncz = slab.ncells
    heightmap = np.ones((Ncx, Ncz)) * (1.0-amp)
    cell_size_x = ncx_block + ncx_gap
    cell_size_z = ncz_block + ncz_gap
    max_coords = np.max(slab.arrays['positions'], axis=0)
    assert amp <= 1.0 and amp >= 0.0
    assert Ncx % cell_size_x == 0 and Ncz % cell_size_z == 0
    blocks_x = Ncx / cell_size_x
    blocks_z = Ncz / cell_size_z
    for i in xrange(blocks_x):
        px = i*cell_size_x
        begin_x, end_x = px, px+ncx_block-1
        for k in xrange(blocks_z):
            pz = k*cell_size_z
            begin_z, end_z = pz, pz+ncz_block-1
            heightmap[begin_x:end_x,begin_z:end_z] = 1.0
    # print heightmap
    if invert:
        heightmap = 1 - heightmap + (1 - amp)
    heightmap *= max_coords[1]
    heightmap = np.roll(heightmap, shift_x, axis=0)
    heightmap = np.roll(heightmap, shift_z, axis=1)
    X = np.linspace(0.0, max_coords[0], Ncx)
    Z = np.linspace(0.0, max_coords[2], Ncz) 
    X_grid, Z_grid = np.meshgrid(X,Z)
    return X_grid, Z_grid, heightmap

    

HEIGHTMAP_GENERATORS = {"MARS": heightmap_generator_MARS,
                        "BLOCKS": heightmap_generator_BLOCKS,
                        "SINE2D": heightmap_generator_SINE2D,}

def build_surface_topology(slab, ncells, heightmap_function=None, heightmap_params=None, plot_heightmap=True):
    slab.ncells = ncells
    # slab = fcc111('Cu', size=(10,3,10))
    # Slab min max dimensions, shift to 0
    xmin, ymin, zmin = np.min(slab.arrays['positions'], axis=0)
    slab.arrays['positions'][:,0] -= xmin
    slab.arrays['positions'][:,1] -= ymin
    slab.arrays['positions'][:,2] -= zmin
    xmax, ymax, zmax = np.max(slab.arrays['positions'], axis=0)
    # print slab.size
    # Show info
    print "Slab dimensions: {} x {} x {}.".format(xmax, ymax, zmax)
    if  heightmap_function is not None:
        assert heightmap_params is not None
        hf = HEIGHTMAP_GENERATORS[heightmap_function]
        X, Z, heightmap = hf(*(heightmap_params + (slab,)))

    keep_list = [atom.index for atom in slab if mask_atom(heightmap, atom.position, (xmax, ymax, zmax))]
    slab = slab[keep_list]
    return X, Z, heightmap, slab


def create_mars_heightmap(Lx, Lz, nx, nz, ncorrx, ncorrz, cutoff, fout="heightmap.dat"):
        ppc = 1
        dx = Lx / (nx*ppc)
        dz = Lz / (nz*ppc)
        s= mars.surface(ncorrx*ppc,ncorrz*ppc,nx*ppc,nz*ppc,cutoff,dx,dz)
        # Step 1: Specify ACF
        acf= s.acf()
        # Step 2: Assemble & solve nonlinear system of equations
        guess= s.f0()
        alpha= s.ncgm(guess)
        # Step 3: Generate a random number matrix
        rand= s.eta()
        # Step 4: Generate the heightmap
        hmap= s.heightmap(alpha,rand)
        # Step 5: Save the surface
        s.save(fout)


