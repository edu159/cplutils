import numpy as np

def timestep_generator(time_steps, interval, twrite):
    if isinstance(interval, (float, int)):
        interval_list = [interval] * len(time_steps)
    elif isinstance(interval, (list, tuple)):
        if len(time_steps) != len(interval):
            raise Exception("List of time_steps and interval has to be the same length.")
        else:
            interval_list = list(interval)
    for t in interval_list:
        if float(t)%float(twrite) > 1e-15:
            raise Exception("Interval is not an integer number of twrite.")
    tv_it = (t/2.0 for t in interval_list)
    nsteps_it = (t/twrite+1 for t in interval_list)
    tsteps_out = []
    for t in time_steps:
        tv = tv_it.next()
        nsteps = nsteps_it.next()
        tsteps_out += list(np.linspace(t-tv, t+tv, nsteps))
    tsteps_out = list(set(tsteps_out))
    tsteps_out.sort()
    return tsteps_out
