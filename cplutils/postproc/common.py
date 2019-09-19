import numpy as np
import os
import pickle

unit_labels = {"lj" : {
                  "velocity": r"(\sigma\tau^{-1})",
                  "density": r"(m\sigma^{-3})",
                  "pressure": r"(\epsilon\sigma^{-3})",
                  "stress":r"(\epsilon\sigma^{-3})",
                  "temperature": r"(\epsilon K_{b}^{-1})",
                  "time": r"(\tau)",
                  "length": r"(\sigma)", 
                  "nu": r"(\sigma^{2}\tau^{-1})", }, 
              "real" : { 
                  "velocity": r"(m/s)",
                  "density": r"(Kg/m^3)",
                  "pressure": r"(MPa)",
                  "stress": r"(MPa)",
                  "temperature": r"(K)",
                  "time": r"(ps)",
                  "length": r"(nm)",
                  "nu": r"(m^{2}s^{-1})", }, 
        }

unit_lammps2si_factors = {"lj" : {
                  "velocity": 1.0, 
                  "density": 1.0,
                  "pressure": 1.0,
                  "stress": 1.0,
                  "temperature": 1.0,
                  "time": 1.0,
                  "nu": 1.0,
                  "length": 1.0, }, 
              "real" : { 
                  "velocity": 1e5, # (Angs/fs -> m/s)
                  "density": 1e3, # (g/cm3 -> kg/m3)
                  "pressure": 0.101335, # (Atm -> MPa)
                  "stress": 0.101335, # (Atm -> MPa)
                  "temperature": 1.0,
                  "time": 1e-3, # (fs -> ps)
                  "nu": 1e-5, # (Angs^2/fs -> m2/s)
                  "length": 0.1, }, # (Angs -> nm)
        }



field_labels = { "velocity-x": r"V_{x}",
                  "density": r"\rho",
                  "pressure": r"P",
                  "stress-xy":r"\tau_{xy}",
                  "temperature": r"T",
                  "time": r"t",
                  "length": r"L", 
                  "nu": r"\nu",
               }

def get_conversion_factor(field, units):
    strip_fn = field.split("-")[0]
    return unit_lammps2si_factors[units][strip_fn]

def get_field_label(field_name, units, label=""):
    fn = ""
    if label:
        fn = label
    else:
        fn = field_labels[field_name]
    strip_fn = field_name.split("-")[0]
    return r"$%s\ %s$" % (fn, unit_labels[units][strip_fn])


def compute_subarray(x, data, dx, xmin, xmax, incxmax=False):
    dx2 = dx / 2.0
    if incxmax:
        xtop = xmax
    else:
        xtop = xmax - dx2
    xbot = xmin + dx2 
    ncells = int((xtop - xbot)/dx) 
    x_new = np.linspace(xbot, xbot + ncells * dx, ncells + 1)
    y_new = np.zeros((data.shape[0], len(x_new)))
    for t in xrange(data.shape[0]):
        y_new[t, :] = np.interp(x_new, x, data[t, :])
    return y_new, x_new

error_methods = ["rel-abs", "abs", "diff", "rel-diff"]
def compute_error(d1, d2, t, y, method="abs"):
    # Avoid nans
    eps = 1e-15
    d1 += eps
    d2 += eps
    fmax = np.maximum(np.absolute(d1), np.absolute(d2))
    if method == "rel-abs":
        error = np.divide(np.absolute(d1 - d2), fmax) * 100
    elif method == "abs":
        error = np.absolute(d1 - d2)
    elif method == "diff":
        error = d1 - d2
    elif method == "rel-diff":
        error = np.divide(d1 - d2, fmax) * 100
    return error

def compute_mean_field(data, tidx, tavg, dt=None, times=None):
    assert dt is not None or times is not None
    if times is None:
        tfin = max(tidx.keys()) - tavg/2.0
        #NOTE: time 0 is not used as there is a spike in the analytical solution
        tini = tavg/2.0 + dt
        # Correct to the closest value 
        tfin = tfin - (tfin-tini) % dt
        nsteps = int((tfin-tini) / dt) + 1
        times = np.linspace(tini, tfin, nsteps)
    else:
        nsteps = len(times)
    data_out = np.zeros((nsteps, ) + data.shape[1:])
    for i, t in enumerate(times):
        begin = t - tavg/2.0
        end = t + tavg/2.0
        if begin == end:
            data_out[i, ...] = data[tidx[begin]]
        else:
            data_out[i, ...] = np.mean(data[tidx[begin]:tidx[end],...], axis=0)
    return times, data_out

def load_fields(path, domain_fields, units, axis=["x", "y", "z"] , ignore_errors=False):
    domain_names = domain_fields.keys()
    fields = {}
    for dom in domain_names:
        fields[dom] = {"fields": {}}
        try:
            for xi in axis:
                axis_fname = "%s-%s.npy" % (dom, xi)
                fields[dom][xi] = np.load(os.path.join(path, axis_fname))
        except Exception:
            raise Exception("Problem loading axis '%s' for domain '%s'." % (axis_fname, dom))
        tidx_fname = "%s-tidx" % dom
        try:
            with open(os.path.join(path, tidx_fname)) as f:
                fields[dom]["tidx"] = pickle.load(f)
        except Exception:
            raise Exception("Problem loading time index '%s' for domain '%s'." % (tidx_fname, dom))
        field_names = domain_fields[dom]
        for f in field_names:
            if type(f) is tuple:
                f_name = f[0]
                indexes = f[1]
                # Eat-away the last dimension when it is of size 1
                if len(indexes) == 1:
                    indexes = indexes[0]
            else:
                f_name = f
                indexes = slice(None)
            try:
                field_fname = "%s-%s.npy" % (dom, f_name)
                field_path = os.path.join(path, field_fname)
                fields[dom]["fields"][f_name] = np.load(field_path)[..., indexes]
                # Remove last dimension if the size is 1. Makes comparison consistent in later operations.
                fshape = fields[dom]["fields"][f_name].shape
                if fshape[-1] == 1:
                    fields[dom]["fields"][f_name] = np.squeeze(fields[dom]["fields"][f_name],
                                                               axis=len(fshape)-1)
            except Exception as e:
                if not ignore_errors:
                    raise Exception("Problem loading field '%s' for domain '%s'." % (field_fname, dom))
                else:
                    print "Warning: Field %s not loaded." % field_fname
    return fields

def compute_fields_error(fields, error_fields, f1, f2, dx, xmin, xmax, save_path=None, err_name=None):
    error_methods = ["rel-abs", "abs", "diff", "rel-diff"]
    for field_name in error_fields:
        interp_field_1, interp_y_1 = compute_subarray(fields[f1]["y"],
                                                    fields[f1]["fields"][field_name],
                                                    dx, xmin, xmax)
        interp_field_2, interp_y_2 = compute_subarray(fields[f2]["y"],
                                                    fields[f2]["fields"][field_name],
                                                    dx, xmin, xmax)
        for em in error_methods:
            error = compute_error(interp_field_1, interp_field_2, fields[f1]["t"], interp_y_1, method=em)
            if save_path is not None:
                if err_name is None:
                    err_name = "%s-%s" % (f1, f2)
                np.save(os.path.join(save_path, "error_%s_%s_%s" % (err_name, field_name, em)), error)
        if save_path is not None:
            np.save(os.path.join(save_path, "error_%s_y" % err_name), interp_y_1)

def stats_subsamples(Ns, std_samples, mean_samples):
    if type(Ns) == int:
        samples = np.ones(len(mean_samples)) * Ns
    else:
        samples = Ns
    meanN = np.sum(mean_samples*samples) / np.sum(samples)
    m = np.mean(mean_samples)
    N = np.sum(samples)
    K = np.sum((samples - 1)*std_samples**2 + samples*(mean_samples - meanN)**2)
    stdN = np.sqrt(K/(N - 1))
    semN = stdN/np.sqrt(N-1)
    return  meanN, stdN, semN



