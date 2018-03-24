import md
import numpy as np
from cplutils.readers.lammps import chunk_reader
import time
import matplotlib.pyplot as plt

# a = np.array([[(1,2,3), (2,3,4), (2,2,3)], [(2,3,4), (1,1,1), (2,3,1)]])
# corr = md.correlation(a, orig_points=1, norm=True)
# print corr

# print "Loading..."
# a = time.time()
# output = chunk_reader("com.out", "molecule", order="frame")
# print "Time:", a - time.time()
#
# a = time.time()
# g1, axis1 = md.rdf(output["data"], output["data"], 100, 2.5, no_procs=1)
# print "Time1:", a - time.time()
# print g1
# a = time.time()
# g2, axis2 = md.rdf(output["data"], output["data"], 100, 2.5, no_procs=2)
# print "Time2:", a - time.time()
# a = time.time()
# g3, axis3 = md.rdf(output["data"], output["data"], 100, 2.5, no_procs=3)
# print "Time3:", a - time.time()
# a = time.time()
# g4, axis4 = md.rdf(output["data"], output["data"], 100, 2.5, no_procs=4)
# print "Time4:", a - time.time()
#
# if (g1==g2).all() and (g1==g3).all() and (g1==g4).all():
#     print "SON iguales!"


# corr = md.Correlation(output["data"], orig_points=1000)
# corr.run()
# print "time:",  corr.compute_time
# plt.plot(corr.corr_func)
# plt.show()

g_vmd = np.loadtxt("multiplot.dat")
system = md.atomic_rdf()
system = np.asfortranarray(np.transpose(system,(1,0,2)))
system2 = chunk_reader("traj.mpiio.out", "trajectory-xyz", order="frame")
# print system2["data"].shape, system.shape
# print system2["data"]
# print "----------------"
# print system
if np.allclose(system,system2["data"], rtol=1e-5):
    print "Son iguales."
g, axis = md.rdf(system2["data"], system2["data"], 100, 8.0, (17.0291, 17.0291, 17.0291), no_procs=4, auto=True)
plt.plot(g_vmd[:,0], g_vmd[:,1],"--")
plt.plot(axis, g, "--")
plt.show()
