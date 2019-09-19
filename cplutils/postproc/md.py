import numpy as np
from cplutils.task import Task
import signal
import sys
from multiprocessing import Pool, TimeoutError

#TODO: Can be added support to do VACF between different molecules?
class Correlation(Task):
    def __init__(self, atoms, corr_length=None, orig_points=10, norm=True):
        Task.__init__(self)
        self.atoms = atoms
        self.corr_length = corr_length
        self.orig_points = orig_points
        self.norm = norm
        self.corr_func = None

    def _run(self):
        # Format of atoms ndarray should be 'C' and 'row' to be efficient.
        max_corr_length = self.atoms.shape[1] - self.orig_points + 1
        if self.corr_length is None:
            self.corr_length = max_corr_length
            if self.corr_length < 1:
                raise Exception("Number of origin points is larger or equal to number of time frames.")
        if self.corr_length > max_corr_length:
            raise Exception("There are '{}' samples available. Correlation length cannot be greater than '{}' for orig_points={}."\
                            .format(self.atoms.shape[1], max_corr_length, self.orig_points))
        self.corr_func = np.zeros(self.corr_length)
        corr_func_t = np.zeros(self.corr_length)
        N = self.atoms.shape[0]
        for orig in xrange(self.orig_points):
            self.progress = float(orig)/float(self.orig_points)*100.0
            if self.print_progress and orig % self.print_every == 0:
                print "Progress %s%%..." % self.progress
            corr_func_t[:] = 0.0
            end = orig+self.corr_length
            for atom_vels in np.rollaxis(self.atoms[:,orig:end,:], 0):
                vel0 = atom_vels[0,:]
                corr_func_t += np.matmul(atom_vels, vel0)
            corr_func_t /= N
            # Normalize by default
            if self.norm:
                self.corr_func += corr_func_t/corr_func_t[0]
            else:
                self.corr_func += corr_func_t
        self.corr_func = self.corr_func/self.orig_points


def correlation(atoms, corr_length=None, orig_points=10, norm=True):
    corr = Correlation(atoms, corr_length, orig_points, norm)
    corr.run()
    return corr.corr_func


def rdf_sub(args):
    atoms1, atoms2, nbins, cutoff, box, auto = args
    gauss_sum = lambda n: n*(n+1)/2
    N1 = atoms1.shape[1]
    N2 = atoms2.shape[1]
    time_steps = atoms1.shape[0]
    print "Time_steps", time_steps
    assert atoms1.shape[0] == atoms2.shape[0]
    if auto:
        no_pairs = gauss_sum(N1-1)
    else:
        no_pairs = N1*N2
    distances = np.zeros(no_pairs)
    histo, axis = np.histogram(distances, nbins, (1e-6, cutoff))
    histo *= 0
    def distance(x0, x1, dimensions):
        delta = np.abs(x0 - x1)
        delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
        return np.sqrt((delta ** 2).sum(axis=-1))
    for t in xrange(time_steps):
        for i in xrange(N1):
            if auto:
                orig = i*N1-gauss_sum(i)
                end = (i+1)*N1-gauss_sum(i+1)
                distances[orig:end] = distance(atoms1[t,i,:], atoms2[t,i+1:N1,:], np.array(box))
            else:
                orig = i*N2
                end = orig + N2
                distances[orig:end] = distance(atoms1[t,i,:], atoms2[t,:,:], np.array(box))
        histo += np.histogram(distances, nbins, (1e-6, cutoff))[0]
    return histo, axis 


def rdf(atoms1, atoms2, nbins, cutoff, box=(1.0,1.0,1.0), auto=False, no_procs=1):
    assert no_procs > 0
    assert atoms1.shape[0] == atoms2.shape[0]
    time_steps = atoms1.shape[0]
    N1 = atoms1.shape[1]
    N2 = atoms2.shape[1]
    chunks_perproc = int(time_steps / no_procs)
    remaining = time_steps % no_procs
    chunks_bounds = [[c*chunks_perproc, (c+1)*chunks_perproc] for c in xrange(no_procs)]
    # The last processor gets the remaining chunks
    chunks_bounds[-1][1] += remaining
    chunks = []
    for begin, end in chunks_bounds:
        chunks.append((atoms1[begin:end,:,:], atoms2[begin:end,:,:], nbins, cutoff, box, auto))
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = Pool(no_procs)
    signal.signal(signal.SIGINT, original_sigint_handler)
    result = []
    # result = pool.map(rdf_sub, chunks)
    try:
        partial_rdfs = pool.map_async(rdf_sub, chunks)
        # Keep running if timeout happens. Timeout is necessary to
        # be able to pick up the SIGINT. It is how it works get()...
        while True:
            try:
                result = partial_rdfs.get(600)
            except TimeoutError:
                continue
            break
    except KeyboardInterrupt:
        pool.terminate()
        pool.join()
        sys.exit(0)
    pool.close()
    pool.join()
    histo, axis = result.pop()
    for rdf in result:
        histo += rdf[0]
    shell_vol = 4.0/3.0 * np.pi * (np.power(axis[1:], 3) - np.power(axis[:-1], 3))
    box_vol = box[0] * box[1] * box[2]
    den = float(N1)/box_vol
    norm_factor = N1 * den * shell_vol
    # The 2 factor is for the double counting avoided
    if auto:
        norm_factor /= 2.0
    g = histo/(norm_factor*float(time_steps))
    axis = axis[1:]-cutoff/(2.0*nbins)
    return g, axis

def msd(atoms):
    #Needs 'frame' ordering
    return np.average(atoms, axis=1)

def atomic_rdf():
    import MDAnalysis
    u = MDAnalysis.Universe("system.data", "traj.dcd", format="LAMMPS", in_memory=True)
    u.trajectory.timeseries().shape
    g = u.select_atoms("type 1")
    return u.trajectory.timeseries()
