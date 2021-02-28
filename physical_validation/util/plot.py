###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (C) 2012 University of Virginia                            #
#              (C) 2017 University of Colorado Boulder                    #
#                                                                         #
#    This library is free software; you can redistribute it and/or        #
#    modify it under the terms of the GNU Lesser General Public           #
#    License as published by the Free Software Foundation; either         #
#    version 2.1 of the License, or (at your option) any later version.   #
#                                                                         #
#    This library is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    Lesser General Public License for more details.                      #
#                                                                         #
#    You should have received a copy of the GNU Lesser General Public     #
#    License along with this library; if not, write to the                #
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     #
#    Boston, MA 02110-1301 USA                                            #
#                                                                         #
###########################################################################

import warnings

import numpy as np


def plot(
    res,
    legend=None,
    title=None,
    xlabel=None,
    ylabel=None,
    xlim=None,
    ylim=None,
    inv_x=False,
    inv_y=False,
    sci_x=False,
    sci_y=False,
    axtext=None,
    annotation_location=None,
    percent=False,
    filename=None,
    screen=True,
):

    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.ticker import AutoMinorLocator, FuncFormatter
    except ImportError:
        warnings.warn("Install matplotlib to enable plotting.")
        return

    def to_percent(y_ticks, _):
        # Adapted from https://matplotlib.org/examples/pylab_examples/histogram_percent_demo.html
        return "{:g}".format(100 * y_ticks)

    font = {"family": "serif", "weight": "normal", "size": 16}
    mpl.rc("font", **font)

    plt.ioff()
    fig, ax = plt.subplots()
    xmin = float("inf")
    xmax = float("-inf")
    for r in res:
        if "args" in r:
            args = r["args"]
        else:
            args = dict()
        if "name" in r:
            # backwards compatibility
            args["label"] = r["name"]

        if "hist" in r:
            y = r["y"]
            _, x, _ = ax.hist(y, r["hist"], **args)
        else:
            x = r["x"]
            y = r["y"]
            if xlim is not None:
                x = x[(r["x"] >= xlim[0]) & (r["x"] <= xlim[1])]
                y = y[(r["x"] >= xlim[0]) & (r["x"] <= xlim[1])]
            if "y_err" in r:
                dy = r["y_err"]
                if xlim is not None:
                    dy = dy[(r["x"] >= xlim[0]) & (r["x"] <= xlim[1])]
                ax.errorbar(x, y, yerr=dy, **args)
            else:
                ax.plot(x, y, **args)

        xmin = min(np.min(x), xmin)
        xmax = max(np.max(x), xmax)

    if legend is not None:
        ax.legend(loc=legend)
    box = ax.get_position()
    if title is not None:
        ax.set_title(title, y=1.05)
        box = box.from_bounds(box.x0, box.y0, box.width, box.height * 0.95)
    if xlabel is not None:
        ax.set_xlabel(xlabel, labelpad=5)
        box = box.from_bounds(
            box.x0, box.y0 + 0.05 * box.height, box.width, box.height * 0.95
        )
    if ylabel is not None:
        ax.set_ylabel(ylabel, labelpad=10)
        box = box.from_bounds(
            box.x0 + 0.05 * box.width, box.y0, box.width * 0.95, box.height
        )
    ax.set_position([box.x0, box.y0, box.width, box.height])
    ax.axis("auto")
    if xlim is not None:
        ax.set_xlim(xlim)
    elif np.isfinite(xmin) and np.isfinite(xmax):
        ax.set_xlim([xmin, xmax])
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    if inv_x:
        ax.invert_xaxis()
    if inv_y:
        ax.invert_yaxis()

    if axtext is not None:
        if isinstance(axtext, str):
            axtext = [axtext]
        if annotation_location is None:
            annotation_location = [None for _ in axtext]
        if isinstance(annotation_location, tuple):
            annotation_location = [annotation_location]
        for t, loc in zip(axtext, annotation_location):
            bbox = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
            if loc is None:
                ax.text(
                    0.95,
                    0.05,
                    t,
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    bbox=bbox,
                )
            else:
                ax.text(loc[0], loc[1], t, bbox=bbox)

    if percent:
        formatter = FuncFormatter(to_percent)
        ax.yaxis.set_major_formatter(formatter)

    if sci_x:
        ax.ticklabel_format(style="sci", axis="x", scilimits=(-3, 4))
    if sci_y:
        ax.ticklabel_format(style="sci", axis="y", scilimits=(-3, 4))
    ax.xaxis.major.formatter._useMathText = True

    if filename is not None:
        fig.savefig(filename, dpi=300)
    if screen:
        fig.show()

    plt.ion()
    if not screen:
        plt.close(fig)
