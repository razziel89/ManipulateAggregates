"""Class definition to control Gnuplot from Python.
"""

# This file is part of ManipulateAggregates.
#
# Copyright (C) 2016 by Torsten Sachse
#
# ManipulateAggregates is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ManipulateAggregates is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ManipulateAggregates.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import sys
import itertools
from subprocess import Popen, PIPE
from shutil import move
import tempfile as tf
import distutils.spawn as ds

from . import gpcsv
from . import postprocess

# default values for gnuplot
GNUPLOT_DEFAULT = {
    "dash": "AUTO",
    "color": "AUTO",
    "point": "AUTO",
    "lines": True,
    "points": False,
    "title": "AUTO",
    "xcol": 1,
    "ycol": 2,
    "zcol": 3,
    "type": "UNKNOWN",
    "head": False,
    "linewidth": 1.0,
    "pointsize": 1.0,
    "solid": True,
    "opacity": 1.0,
}


def _which(executable):
    """Return whether or not an executable was found.

    Args:
        executable: (string) the name of the executable

    Returns:
        whether or not the specified executable was found
    """
    return ds.find_executable(executable) is not None


def _mv(src, dest):
    """Move a file from src to dest.

    Args:
        src: (string) source path of the file
        dest: (string) destination path of the file

    """
    move(src, dest)


class gnuplot:
    """Controller class for gnuplot.

    Attributes:
        correct_mark: (bool) whether or not marks should be slightly displaced if they
            overlap
        correct_mark_dist: (float) by how much to displace overlaping marks
        dict: (dictionary) holds configuration parameters
        dimensions: (int, 2 or 3): the dimensionality of the plot
        f: (pipe) STDIO of GP to which the controlling information is written.
        font: (string) font to use
        fontsize: (int) fontsize to use
        GP: (process) an instance of Gnuplot opened via Popen
        rectanglecount: (int) how many rectangles were already drawn
        tempfiles: list of temporary files and whether they shall be auto-deleted
        xmarks: (dictionary) used to store marks in x direction - check for overlaps
        xrange: (tuple of 2 floats) range for x axis
        ymarks: (dictionary) used to store marks in y direction - check for overlaps
        yrange: (tuple of 2 floats) range for y axis
        zrange: (tuple of 2 floats) range for y axis if plot in 3 dimensions
    """

    class gp_dict:
        """A dummy class to hold the description of the config dictionaries.

        All members of this class are actually keys (of type string) that can
        be in each config dictionary and the given type is that of the
        associated value. The two letter abbreviations lc, lw and dt are
        Gnuplot commands. Please see Gnuplot's documentation.
        
        Attributes:
            lines: (none) if the key is set, plot a line
            points: (none) if the key is set, plot points
            type: (string) type of plot: 'file' (raw data) or 'function' (analytic function)
            function: (string) description of the analytic function
            filename: (string) the file that contains the raw data
            xcol: (int) coloumn in the file that contains x data
            ycol: (int) coloumn in the file that contains y data
            zcol: (int) coloumn in the file that contains z data
            border: (none) when plotting a rectange, it will have a border if this is set
            color: (int) colour as in 'lc INT'. Special value: 'AUTO'
            dash: (int) dash type as in 'dt INT'. Special value: 'AUTO'
            point: (int) point type as in 'dt INT'. Special value: 'AUTO'
            head: (none) whether arrows shall have heads
            linewidth: (float) line width as in 'lw FLOAT'
            opacity: (float) opacity of the rectangle (if solid is set)
            pointsize: (float) size of the points (if points is declared)
            solid: (none) if set, plot a solid rectangle (not just border)
            title: (string) the name of the plotted function/data
        """

        def __init__(self):
            """Dummy constructor, do not use."""
            raise Exception("This dummy class shall not be used directly")
            self.lines = ""
            self.points = ""
            self.type = ""
            self.function = ""
            self.filename = ""
            self.xcol = ""
            self.ycol = ""
            self.zcol = ""
            self.border = ""
            self.color = ""
            self.dash = ""
            self.point = ""
            self.head = ""
            self.linewidth = ""
            self.opacity = ""
            self.pointsize = ""
            self.solid = ""
            self.title = ""

    def __init__(
        self,
        filename,
        size=(20, 12),
        linewidth=1,
        xrange=None,
        yrange=None,
        correct_mark=True,
        correct_mark_dist=0.001,
        fontsize=10,
        xlog=False,
        ylog=False,
        classic_colors=True,
        dimensions=2,
        font="Helvetica",
        verbose=False,
    ):
        """Constructor.

        Args:
            filename: (string) name of the output file that will be created
            size: (tuple of 2 ints) the (x,y) size of the output in cm
            linewidth: (float) lines are scaled by this factor
            xrange: (tuple of 2 floats) range for x axis
            yrange: (tuple of 2 floats) range for y axis
            correct_mark: (bool) whether or not marks should be slightly displaced if
                they overlap
            correct_mark_dist: (float) by how much to displace overlaping marks
            fontsize: (int) fontsize to use
            xlog: (bool) whether or not the x axis is on a logarithmic scale
            ylog: (bool) whether or not the y axis is on a logarithmic scale
            classic_colors: (bool) whether or not to use Gnuplot's classic color scheme
            dimensions: (int, 2 or 3): the dimensionality of the plot
            font: (string) font to use
            verbose: (bool) if True, echo everything that is being passed to Gnuplot to
                stderr

        """

        self.tempfiles = []
        self.GP = Popen(["gnuplot"], stdin=PIPE, stdout=sys.stderr, stderr=sys.stderr)
        self.cat = Popen(["cat"], stdin=PIPE, stdout=sys.stderr, stderr=sys.stderr)
        self.f = self.GP.stdin
        self.rectanglecount = 1
        self.dimensions = dimensions
        self.fontsize = fontsize
        self.font = font
        self.filename = filename
        self.verbose = verbose
        if self.dimensions not in (2, 3):
            raise ValueError("Wrong number of dimensions provided.")
        self._write(
            'set term pdf enhanced colour font "%s,%d" size %dcm,%dcm linewidth %d\n'
            % (font, fontsize, size[0], size[1], linewidth)
        )
        if classic_colors:
            self._write("set colors classic\n")
        if xlog:
            self._write("set logscale x\n")
        if ylog:
            self._write("set logscale y\n")
        self.correct_mark = correct_mark
        if correct_mark:
            self.xmarks = {}
            self.ymarks = {}
            self.correct_mark_dist = correct_mark_dist
        self.xrange = xrange
        if self.xrange is not None:
            self._write("set xrange [%f:%f]\n" % (xrange[0], xrange[1]))
        self.yrange = yrange
        if self.yrange is not None:
            self._write("set yrange [%f:%f]\n" % (yrange[0], yrange[1]))
        self._write('set output "%s.pdf"\n' % (filename))

    def _write(self, s):
        """Write something to Gnuplot but also print it to
        stderr if verbose output is requested.

        Args:
            s (string): string to be passed to Gnuplot
        """
        self.f.write(s)
        if self.verbose:
            self.cat.stdin.write(s)

    def _set_dict(self, dict):
        """Set the dictionary holding config parameters.

        Each function to be plotted has a certain set of possible config
        options. See gnuplot.gp_dict for all possible options.

        Args:
            dict: (dictionary) holds configuration parameters in the form of
                key-value pairs.
        """
        self.dict = dict

    def _get(self, *args, **kwargs):
        """Retrieve a value from the dictionary.

        Args:
            `*args`: (strings) the config options whose associated value shall be
                retrieved.
            `**kwargs`: (dictionary) key strict: whether or not to raise an Error
                if the key cannot be found in the current dictionary or the
                default one. If False then None is returned in such cases.
        """
        for k in args:
            result = self.dict.get(k, GNUPLOT_DEFAULT.get(k, None))
            if result is not None:
                break
        if kwargs.get("strict", False) and result is None:
            raise KeyError(
                "Key %s not provided by the current dictionary and no default set."
                % (key)
            )
        return result

    def set_title(self, title):
        """Set the title of the plot.

        Args:
            title: (string) the title of the generated plot
        """
        self._write('set title "%s"\n' % (title))

    def emptyplot(self):
        """Create an empty plot.

        This is useful if only marks or arrows shall be plotted. Normally,
        Gnuplot would not create a plot under such conditions.

        Requires xrange and yrange to be set using set_xrange and
        set_yrange.
        """
        if self.xrange is None or self.yrange is None:
            raise RuntimeError(
                "Cannot perform emtpy plot if either xrange or yrange is not set."
            )
        self._write("plot NaN notitle\n")
        self.postplot()

    def postplot(self):
        """Unset arrows and labels.

        This is required when creating a plot that contains multiple pages since
        otherwise the labels and arrows/marks would be repeated on every page.
        """
        self._write("unset arrow\n")
        self._write("unset label\n")

    def unset(self, prop, oneprop=True):
        """Send an arbitrary 'unset' command to Gnuplot.

        Args:
            prop: (string or iterable of strings) if oneprop is True (the
                default), unset the one property given via prop. Otherwise
                unset all properties in the iterable prop.
            oneprop: (bool) whether or not prop is not an iterable
        """
        if oneprop:
            iterable = [prop]
        else:
            iterable = prop
        for p in iterable:
            self._write("unset %s\n" % (p))

    def set(self, prop, oneprop=True):
        """Send an arbitrary 'set' command to Gnuplot.

        Args:
            prop: (string or iterable of strings) if oneprop is True (the
                default), et the one property given via prop. Otherwise set
                all properties in the iterable prop.
            oneprop: (bool) whether or not prop is not an iterable
        """
        if oneprop:
            iterable = [prop]
        else:
            iterable = prop
        for p in iterable:
            self._write("set %s\n" % (p))

    def lineplot(self, data):
        """Plot one or several functions (can also be raw data).

        Each function has a certain set of possible config options. See
        gnuplot.gp_dict for all possible options.

        Args:
            data: (dictionary or list of dictionaries) each dictionary contains
                a set of key-value pairs that define the function to be plotted
                and how it shall be formated.
        """
        if isinstance(data, dict):
            dict_list = [data]
        else:
            try:
                if False in (isinstance(d, dict) for d in data):
                    raise TypeError("")
                else:
                    dict_list = data
            except TypeError:
                raise TypeError(
                    "Data for lineplot is neither a dictionary nor a list of dictionaries."
                )
        if len(dict_list) == 0:
            print("WARNING: cannot plot since no data was passed over", file=sys.stderr)
            return
        breakchar = ", "
        if self.dimensions == 2:
            self._write("plot ")
        elif self.dimensions == 3:
            self._write("splot ")
        count = 1
        for f in dict_list:
            self._set_dict(f)
            if f == dict_list[-1]:
                breakchar = "\n"

            if not self._get("lines") and not self._get("points"):
                raise ValueError(
                    "At least one of 'lines' or 'points' has to be declared otherwise nothing would be plotted."
                )

            if self._get("type") == "function":
                self._write("%s " % (self._get("function")))
            elif self._get("type") == "filename" or self._get("type") == "file":
                self._write('"%s" u ' % (self._get("filename", "file", strict=False)))
                # x coloumn
                if isinstance(self._get("xcol"), int):
                    self._write("($%d):" % (self._get("xcol")))
                else:
                    self._write("(%s):" % (self._get("xcol")))
                # y coloumn
                if isinstance(self._get("ycol"), int):
                    self._write("($%d)" % (self._get("ycol")))
                else:
                    self._write("(%s)" % (self._get("ycol")))
                # z coloumn, if present
                if self.dimensions == 3:
                    if isinstance(self._get("zcol"), int):
                        self._write(":($%d) " % (self._get("zcol")))
                    else:
                        self._write(":(%s) " % (self._get("zcol")))
                self._write(" ")
            elif self._get("type") == "UNKNOWN":
                raise ValueError(
                    "No plot type provided. Missing key 'type' from dictionary."
                )
            else:
                raise ValueError("Unknown plot type: %s" % (f["type"]))

            if self._set_style(count, "line"):
                count += 1
            self._write(breakchar)
        if self.correct_mark:
            self.xmarks = {}
            self.ymarks = {}
        self.postplot()

    def data_to_file(self, data, formatstring=None, delete=True):
        """Convert some data (given as x-y value pairs) to a format for Gnuplot.

        Args:
            data: (list of pairs of floats) the data to be plotted
            formatstring: (string) a printf-type string that will be used to
                convert each element of data to a string. Gnuplot will plot
                what's left after the conversion.
            delete: (bool) whether or not the temporary file that is created
                shall be deleted when finalize is called
        """
        tempfile = tf.NamedTemporaryFile(delete=False)
        if formatstring is None:
            for datum in data:
                tempfile.write("\t".join(map(str, datum)) + "\n")
        else:
            for datum in data:
                tempfile.write(formatstring % tuple(datum))
        tempfile.close()
        self.tempfiles.append((delete, tempfile.name))
        return tempfile.name

    def _set_style(self, count, type):
        """Create s string that Gnuplot understands and that describes a plot's style.

        Args:
            count: (int) how many times already the automatic style generation
                has been used
            type: (string) what kind of thing shall be plotted. Can be:
                "lines", "rectangle" or "arrow".
        """
        used_count = False
        style = ""
        if type == "line":
            style += "w "
            if self._get("lines"):
                style += "l"
            if self._get("points"):
                style += "p "
                style += "ps %f " % (self._get("pointsize"))
            style += " "
            dash = self._get("dash")
            if dash == "AUTO":
                style += "dt %d " % (count)
                used_count = True
            else:
                style += "dt %d " % (dash)
            point = self._get("point")
            if point == "AUTO":
                style += "pt %d " % (count)
                used_count = True
            else:
                style += "pt %d " % (point)
            color = self._get("color")
            if color == "AUTO":
                style += "lc %d " % (count)
                used_count = True
            else:
                style += "lc %d " % (color)
            width = self._get("linewidth")
            style += "lw %f " % (width)

            title = self._get("title")
            if title == "AUTO":
                pass
            elif title is None:
                style += "notitle "
            else:
                style += 'title "%s" ' % (title)
        elif type == "rectangle":
            color = self._get("color")
            if color == "AUTO" or color is None:
                style += "fc "
            else:
                style += "fc %s " % (color)
            style += "lt -1 lw 0 "
            if self._get("solid"):
                style += "fs solid %f " % (self._get("opacity"))
            if self._get("border") is None:
                style += "noborder "
        elif type == "arrow":
            dash = self._get("dash")
            if dash == "AUTO":
                style += "dt %d " % (count)
                used_count = True
            else:
                style += "dt %d " % (dash)
            color = self._get("color")
            if color == "AUTO":
                style += "lc %d " % (count)
                used_count = True
            else:
                style += "lc %d " % (color)
            width = self._get("linewidth")
            style += "lw %f " % (width)
            if not self._get("head"):
                style += "nohead "
        self._write(style)
        return used_count

    def set_xrange(self, start, stop):
        """Set the range of the plot in x-direction.

        Args:
            start: (float) start of the x range
            stop: (float) end of the x range
        """
        self._write("set xrange [%f:%f]\n" % (start, stop))
        self.xrange = (start, stop)

    def set_yrange(self, start, stop):
        """Set the range of the plot in y-direction.

        Args:
            start: (float) start of the y range
            stop: (float) end of the y range
        """
        self._write("set yrange [%f:%f]\n" % (start, stop))
        self.yrange = (start, stop)

    def set_zrange(self, start, stop):
        """Set the range of the plot in z-direction if the plot is 3D.

        Args:
            start: (float) start of the z range
            stop: (float) end of the z range
        """
        if self.dimensions == 3:
            self._write("set zrange [%f:%f]\n" % (start, stop))
            self.zrange = (start, stop)
        else:
            raise ValueError("Cannot set zrange for non-3d plot.")

    def set_stick(self, pos, height, color, base=0.0, width=0.5):
        """Create a vertical line of a certain height (i.e., stick).

        Args:
            pos: (float) the x position of the stick
            height: (float) the height in the y direction of the stick
            color: (int) the color if the stick as in 'lc INT'
            base: (float) where the stick shall start (defaults to x axis)
            width: (float) the width of the stick
        """
        try:
            if len(pos) != len(height):
                print(
                    "To print several sticks, the positions list and the height list must have the same number of elements."
                )
            else:
                gen = ((p, h) for p, h in zip(pos, height))
        except TypeError:
            gen = [(pos, height)]
        for p, h in gen:
            self.mark(
                p,
                "x",
                width=width,
                color=color,
                rectangle=False,
                opacity=1.0,
                center=True,
                extent=(base, base + h),
            )

    def set_background(
        self, opacities, colors=None, nr_fields=None, direction="x", extent=None
    ):
        """Create a non-white background for the plot.

        You can create backgrounds of areas with alternating colors or even
        checked backgrounds. This is realized as repeating a given pattern of
        opacities and colours until the entirety of the plot is filled in a
        certain direction. A checked pattern can be obtained by repeating a
        black-and-white pattern multiple times as black-and-white ->
        white-and-black -> black-and-white etc. These backgrounds are realized
        via rectangles, so they support all the properties of Gnuplot's
        rectangles.

        Args:
            opacities: (iterable of floats) pattern of opacities
            colors: (iterable of ints) pattern of colors. Defaults to black for
                all pattern elements.
            nr_fields: (int) the number of sections in which to partition the
                plot. E.g., if given a black-and-white pattern, a value of 5
                would result in black->white->black->white->black.
            direction: ("x" or "y") the direction of the pattern (defaults to "x")
            extent: (tuple of 2 floats) how far in the other direction (that
                not specified by direction) the pattern shall extent

        """
        if direction == "x":
            samerange = self.xrange
            otherrange = self.yrange
        elif direction == "y":
            samerange = self.yrange
            otherrange = self.xrange
        else:
            raise ValueError('Unknown direction "%s", must be x or y.' % (direction))
        if self.dimensions != 2:
            raise RuntimeError("Cannot create background for non-2d plot.")
        if extent is None:
            if otherrange is None:
                raise ValueError(
                    "Cannot create background in %s-direction without other range being set."
                    % (direction)
                )
            else:
                extent = otherrange
        if samerange is None:
            raise ValueError(
                "Cannot create background in %s-direction without same range being set."
                % (direction)
            )
        try:
            if colors is None:
                colors = [None] * len(opacities)
        except TypeError:
            opacities = [opacities]
            colors = [None]
        try:
            if len(colors) != len(opacities):
                raise ValueError(
                    "Cannot create background, colors and opacities do not have the same number of elements."
                )
            else:
                iterable = [(c, o) for c, o in zip(colors, opacities)]
        except TypeError:
            iterable = [(colours, opacities)]

        if nr_fields is None:
            nr_fields = len(colors)

        result = []
        width = 1.0 * (samerange[1] - samerange[0]) / (nr_fields)
        pos = samerange[0]
        count = 0
        for color, opacity in itertools.cycle(iterable):
            self.mark(
                pos,
                direction,
                width=width,
                color=color,
                rectangle=True,
                center=False,
                opacity=opacity,
                extent=extent,
            )
            result.append((pos, pos + width))
            pos += width
            count += 1
            if count == nr_fields:
                break
        return result

    def mark(
        self,
        pos,
        direction,
        width=0.5,
        color=None,
        rectangle=False,
        opacity=1.0,
        center=True,
        extent=None,
        label=None,
        zpos=None,
        dash=None,
    ):
        """Create vertival or horizontal line on the plot.

        If the plot is 3D, the position in the 3rd direction is also required.
        However, the plots are still in planes parallel to the x-y plane.
        
        Args:
            pos: (float) x or y position of the mark (depending on direction)
            direction: ("x" or "y") the direction of the line
            width: (float) the line width
            color: (int) colour as in 'lc INT'
            rectangle: (bool) whether the mark is not just a line but a rectangle
            opacity: (float) opacity of the mark (only used if rectangle)
            center: (bool) whether or not the given position is the mark's center.
                Otherwise, the pos is considered to be the left border (only
                used if rectangle)
            extent: (tuple of 2 floats) the startpoint and endpoint in the
                direction of the line (defaults to: entire plot)
            label: (dictionary) an optional description of an optional label. See
                description of set_label for required and optional keys.
            zpos: (float) position of the mark in a 3D plot. Required if the
                dimensionality of the plot is 3.
        """
        if direction == "x":
            hererange = self.yrange
            heremarks = self.xmarks
            startpos = lambda p, e: (p, e[0])
            endpos = lambda p, e: (p, e[1])
        elif direction == "y":
            hererange = self.xrange
            heremarks = self.ymarks
            startpos = lambda p, e: (e[0], p)
            endpos = lambda p, e: (e[1], p)
        else:
            raise ValueError('Unknown direction "%s", must be x or y.' % (direction))
        if self.dimensions != 2:
            if rectangle:
                raise RuntimeError(
                    "Cannot set %smark as rectangle for non-2d plot." % (direction)
                )
            elif self.dimensions == 3:
                if zpos is None:
                    raise RuntimeError(
                        "Cannot set %smark as arrow for non-2d plot without zpos defined."
                        % (direction)
                    )
            else:
                raise RuntimeError(
                    "Fatal internal error: wrong number of dimensions set. However that happened."
                )
        if extent is None:
            if hererange is None:
                raise ValueError(
                    "Cannot create %smark without other range being set." % (direction)
                )
            else:
                extent = hererange
        if not rectangle:
            if self.correct_mark:
                if pos in heremarks:
                    heremarks[pos] += 1
                    _pos = pos + self.correct_mark_dist * (heremarks[pos])
                else:
                    heremarks[pos] = 0
                    _pos = pos
            self._write("set arrow from %8.7E,%8.7E" % startpos(_pos, extent))
            if self.dimensions == 3:
                self._write(",%8.7E" % (zpos))
            self._write(" to %8.7E,%8.7E" % endpos(_pos, extent))
            if self.dimensions == 3:
                self._write(",%8.7E" % (zpos))
            self._write(" ")
            self._set_dict(
                {
                    "head": False,
                    "color": color if color is not None else 0,
                    "linewidth": width,
                    "dash": 1 if dash is None else dash,
                }
            )
            self._set_style(heremarks[pos], "arrow")
            if opacity != 1.0:
                print(
                    "WARNING: opacity unequal 100% set, but is ignored for xmark that is not a rectangle.",
                    file=sys.stderr,
                )
        else:
            if center:
                pos -= 0.5 * width
            self._write("set obj %d rect from " % (self.rectanglecount))
            self._write("%f,%f to " % startpos(pos, extent))
            self._write("%f,%f " % endpos(pos + width, extent))
            self._set_dict({"color": color, "opacity": opacity, "border": None})
            self._set_style(self.rectanglecount, "rectangle")
            self.rectanglecount += 1
        self._write("\n")
        if label is not None:
            if isinstance(label, dict):
                label = [label]
            for l in label:
                if "where" in l:
                    where = l["where"]
                else:
                    where = "tl"
                if where == "tl":
                    labelpos = (
                        startpos(pos, extent)[0] + l["offset"][0],
                        startpos(pos, extent)[1] + l["offset"][1],
                    )
                    l["pivot"] = "left"
                elif where == "bl":
                    labelpos = (
                        startpos(pos, extent)[0] + l["offset"][0],
                        endpos(pos, extent)[1] + l["offset"][1],
                    )
                    l["pivot"] = "left"
                elif where == "tr":
                    labelpos = (
                        endpos(pos, extent)[0] + l["offset"][0],
                        startpos(pos, extent)[1] + l["offset"][1],
                    )
                    l["pivot"] = "right"
                elif where == "br":
                    labelpos = (
                        endpos(pos, extent)[0] + l["offset"][0],
                        endpos(pos, extent)[1] + l["offset"][1],
                    )
                    l["pivot"] = "right"
                elif where == "c":
                    labelpos = (
                        0.5 * (startpos(pos, extent)[0] + endpos(pos, extent)[0])
                        + l["offset"][0],
                        0.5 * (startpos(pos, extent)[1] + endpos(pos, extent)[1])
                        + l["offset"][1],
                    )
                    l["pivot"] = "center"
                else:
                    raise ValueError(
                        'Wrong value for "where" in dictionary. Must be one of ["tl","bl","tr","br","c"] but is %s.'
                        % (where)
                    )
                l["position"] = labelpos
                self.set_label(l)

    def set_label(self, label):
        """Set a label on a plot.

        The argument label is a dictionary whose entries for "font",
        "fontsize", "pivot", "rotation" and "depth" can overwrite the defaults.
        Needs to have entries for "text" (the label's text) and "position"
        (tuple of floats describing the position).

        Args:
            label: (dictionary) a description of the label. See description of
                set_label for required and optional keys.
        """
        if "font" in label:
            font = label["font"]
        else:
            font = self.font
        if "fontsize" in label:
            fontsize = label["fontsize"]
        else:
            fontsize = self.fontsize
        if "pivot" in label:
            pivot = label["pivot"]
        else:
            pivot = "center"
        if "rotation" in label:
            rotation = label["rotation"]
        else:
            rotation = 0.0
        if "depth" in label:
            depth = label["depth"]
        else:
            depth = "front"
        self._write(
            'set label "%s" at %f,%f font "%s,%d" %s rotate by %.2f %s\n'
            % (
                label["text"],
                label["position"][0],
                label["position"][1],
                font,
                fontsize,
                depth,
                rotation,
                pivot,
            )
        )

    def finalize(self, delete=True, convert=False):
        """Finalize the plot.

        This calls a set of routines that finish the plotting procedure.
        Without calling this, the plot will never actually be created.

        Args:
            delete: (bool) whether or not temporary files that were declared as
                "to-be-deleted" shall actually be deleted.
            convert: (bool) whether or not to convert the eps file to a pdf
                file if the required software is installed
        """
        if not self.f.closed:
            self.f.close()
            rc = self.GP.wait()
            if delete:
                for d, filename in self.tempfiles:
                    if d:
                        os.remove(filename)
        return rc
