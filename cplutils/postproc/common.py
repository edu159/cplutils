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
              "length": r"(\sigma)", }, 
          "real" : { },
        }

field_labels = { "velocity-x": r"V_{x}",
                  "density": r"\rho",
                  "pressure": r"P",
                  "stress-xy":r"\tau_{xy}",
                  "temperature": r"T",
                  "time": r"t",
                  "length": r"L", 
               }

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
        nsteps = int((tfin-tini) / dt) + 1
        times = np.linspace(tini, tfin, nsteps)
    else:
        nsteps = len(times)
    data_out = np.zeros((nsteps, ) + data.shape[1:])
    for i, t in enumerate(times):
        begin = t - tavg/2.0
        end = t + tavg/2.0
        data_out[i, ...] = np.mean(data[tidx[begin]:tidx[end],...], axis=0)
    return times, data_out

def load_fields(path, domain_fields, axis=["x", "y", "z"], ignore_errors=False):
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
            except Exception:
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


