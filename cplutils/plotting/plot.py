import numpy as np
from cplutils.postproc.common import compute_mean_field
from cplutils.plotting.extras import DiscreteSlider
import matplotlib.pyplot as plt
from cplutils.postproc.common import field_labels, get_field_label

def plot_coupled(fields_in, labels=None, tidx=None, tavg=None, dt=None, times=None, slider=False, savefig=None, limits=None, units="lj", plot_opts=None):
    if slider:
        assert (dt is not None or times is not None) and tavg is not None and tidx is not None
        if times is None:
            tfin = max(tidx.keys()) - tavg/2.0
            tini = tavg/2.0 + dt
            nsteps = int((tfin-tini) / dt) + 1
            times = np.linspace(tini, tfin, nsteps)
        else:
            tfin = max(times)
            tini = min(times)


    field_names = list(set([f2 for f1 in fields_in.values() for f2 in f1["fields"].keys()]))
    no_fields = len(field_names)
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif', size=12)
    fig, axarr = plt.subplots(no_fields, sharex=True, squeeze=False)
    fig.subplots_adjust(hspace=0.1)
    
    handles = {f:v for f,v in zip(fields_in.keys(), [[None]*no_fields for f in fields_in.keys()])}
    if slider:
        plt.subplots_adjust(left=0.15, bottom=0.15)

    def update(val):
        plot([val], update=True)
        fig.canvas.draw_idle()

    def plot(times, update=False):
        axarr[-1, 0].set_xlabel(get_field_label("length", units)) 
        for domain, fields in fields_in.items():
            dom_opts = plot_opts[domain]
            for field_name, data in fields["fields"].items():
                n = field_names.index(field_name)
                times, fdf_data = compute_mean_field(data, fields["tidx"],
                                                     tavg, dt=dom_opts["twrite"],
                                                     times=times)
                # fdf_data = data # For tavg=0

                for it, t in enumerate(times):
                    # Search indexes withing the range
                    if limits is not None:
                        i = np.where((fields["y"]+dom_opts["shift"] > limits[0]) &\
                                     (fields["y"]+dom_opts["shift"] < limits[1]) )[0]
                    else:
                        i = slice(None)
                    if update:
                        # Only one handler expected (plot only one at a time)
                        handles[domain][n].set_ydata(fdf_data[:,i])
                        linfit = np.polyfit(fields["y"][i]+dom_opts["shift"], fdf_data[it,i], 1)
                        print "Intersection(%s): y = %f" % (domain, np.poly1d(linfit).roots)
                        fig.canvas.draw_idle()
                    else:
                        handle = axarr[n, 0].plot(fields["y"][i]+dom_opts["shift"],\
                                               fdf_data[it][i], dom_opts["style"], label=dom_opts["tag"])
                        handles[domain][n] = handle[0]
                try:
                    flabel = get_field_label(field_name, units)
                except KeyError:
                    try:
                        flabel = labels[field_name]
                    except Exception:
                        raise Exception("No label found for field '%s'." % field_name) 
                axarr[n, 0].set_ylabel(flabel)
        for ax in axarr[:-1, 0]:
            ax.get_xaxis().set_visible(False)

        for ax in axarr[:, 0]:
            ax.relim()
            ax.autoscale_view()

    if slider:
        axcolor = 'gray'
        stride = plt.axes([0.15, 0.05, 0.75, 0.03], facecolor=axcolor)
        slider_vals = np.array(times)
        stride_slider = DiscreteSlider(stride, 'Time:', tini, tfin, valinit=tini,\
                                       allowed_vals = slider_vals, valfmt="%0.0f")
        stride_slider.on_changed(update)
        plot([times[0]])
    else:
        plot(times)
    axarr[0, 0].set_title(r'Coupled solution')
    if savefig is not None:
        plt.savefig(savefig, dpi=600)
    plt.show()


def plot_error_contour(error, t, y):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    T, Y = np.meshgrid(y,t)
    # contour = ax.contourf(T,Y, error, 10, cmap='jet')
    contour = ax.pcolormesh(T,Y, error, 10, cmap='jet')
    plt.colorbar(contour, shrink=0.8, extend='both')
    plt.show()

def plot_error(plot_title, error_label, error, t, mode="mean", error_bars=False, print_stats=False, plot_type="2d", units="lj"):
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif', size=12)
    if plot_type == "2d":
        plot_rows = 2
        fig, axarr = plt.subplots(2, len(error.keys()), sharex=False, sharey=False, squeeze=False)
    elif plot_type == "contour":
        plot_rows = len(error_label)
        fig, axarr = plt.subplots(plot_rows, len(error.keys()), sharex=False, sharey=False, squeeze=False)
    else:
        raise Exception("Plot type '%s' not recognised." % plot_type)
    plt.suptitle("Error (%s) - %s" % (mode.title(), plot_title))
    for i, label in enumerate(error.keys()):
        assert len(error_label) == len(error[label])
        if print_stats:
            print "Errors - %s" % label
        if plot_type == "2d":
            axarr[0][i].set_title(label)
            axarr[0][i].set(xlabel=get_field_label("time", units), ylabel="$Error$")
            axarr[1][i].set(xlabel=get_field_label("length", units), ylabel="$Error$")
            no_rows = 2
        elif plot_type == "contour":
            axarr[0][i].set_title(label)
            no_rows = len(error_label)

        for l, (e, y) in enumerate(error[label]):
            no_cols = len(error.keys()) - 1
            if print_stats:
                print "%s:" % error_label[l]
                print "   error(mean): mean->", np.mean(e), "stderr->", np.std(e)
                print "   error(sum): ", np.sum(e)
                print "   error(l2): ", np.sqrt(np.sum(e**2))
            if plot_type == "2d":
                if mode == "mean":
                    error_y = np.mean(e, axis=0)
                    error_y_std = np.std(e, axis=0)
                    error_t = np.mean(e, axis=1)
                    error_t_std = np.std(e, axis=1)
                elif mode == "l2":
                    error_y = np.sum(e**2, axis=0)
                    error_t = np.sum(e**2, axis=1)
                elif mode == "mse":
                    error_y = np.sqrt(np.mean(e**2, axis=0))
                    error_t = np.sqrt(np.mean(e**2, axis=1))
                if mode == "mean" and error_bars:
                    axarr[1][i].errorbar(y, error_y, yerr=error_y_std, mec='black', fmt='--o',capsize=5, label=error_label[l])
                    axarr[0][i].errorbar(t, error_t, yerr=error_t_std, fmt='-o', ms=3, capsize=2, label=error_label[l])
                else:
                    if "results" in error_label[l]:
                        axarr[1][i].plot(y, error_y, "-o", label=error_label[l])
                        axarr[0][i].plot(t, error_t, "-o", mfc='r', ms=3, label=error_label[l])
                    else:
                        axarr[1][i].plot(y, error_y, "-x", label=error_label[l])
                        axarr[0][i].plot(t, error_t, "-x", mfc='r', ms=3, label=error_label[l])

                axarr[1][i].legend(loc="upper right")
                axarr[0][i].legend(loc="upper right")
            elif plot_type == "contour":
                T, Y = np.meshgrid(y,t)
                contour = axarr[l][i].contourf(T,Y, e, 10, cmap='jet')
                if i == no_cols:
                    fig.colorbar(contour, shrink=0.8, extend='both', ax=axarr[l][no_cols])
                # Place axis labels proprerly
                ylabel = xlabel = ""
                if i == 0:
                    ylabel = "%s\n" % error_label[l] + get_field_label("time", units)
                if l == 1:
                    xlabel = get_field_label("length", units)
                axarr[l][i].set(xlabel=xlabel, ylabel=ylabel)

        # Make ticks invisible
        for i in xrange(0, no_rows):
            for ax in axarr[i, 1:no_cols+1]:
                ax.get_yaxis().set_visible(False)
         
    fig.savefig(plot_title + ".jpeg")
    plt.show()


