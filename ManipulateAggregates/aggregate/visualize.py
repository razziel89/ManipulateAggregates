# -*- coding: utf-8 -*-
"""A collection of visualization routines for aggregates.

This package requires GLUT (tested with freeGLUT). You only have to call the
function ManipulateAggregates.aggregate.visualize.visualize, the rest will be
taken care of.
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
import sys
import re
import os
import io
import copy

try:
    import cPickle as pickle

    pickleload = lambda h: pickle.load(h)
except ImportError:
    import pickle

    # The following code should ensure that nested dictionaries pickled
    # in Python 2 can also be openend in Python 3 (unless of course there
    # is ever a non-ascii character anywhere).
    def _convert(o):
        if isinstance(o, dict):
            return _convertdict(o)
        elif isinstance(o, bytes):
            return tobasestring(o)
        else:
            return o

    _convertdict = lambda d: {_convert(k): _convert(v) for k, v in d.items()}
    pickleload = lambda h: _convertdict(pickle.load(h, encoding="bytes"))
try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
    from ..collection.opengl import *

    gl_imported = True
except ImportError:
    gl_imported = False

import logging

logger = logging.getLogger(__name__)
try:
    import numpy
except ImportError:
    logger.warning("Could not import numpy")
try:
    from ..collection.p2p3IO import open, writeto, close, hashstring, tobasestring
except ImportError:
    logger.warning("Could not import p2p3IO")

org_translation = [0.0, 0.0, -260.0]

# A collection dictionary for all visualization options.
#  It contains all variables that need to be global for OpenGL to be able to
#  display stuff. Not to be changed by the user.
global gl_c
gl_c = {
    "fullscreen": False,  # whether fullscreen shall be activated right from the start
    "resolution": (1024, 768),  # default resolution in pixels
    "maxextent": -1.0,  # how to scale the view by default (default -1.0 means do not scale)
    "translation": [
        0.0,
        0.0,
        -260.0,
    ],  # position of camera with respect to drawn structures
    "look_at": [0.0, 0.0, 0.0],  # the camera will look at this point
    "angles": [0.0, -10.0, 0.0],  # rotational angles around x,y and z axes
    "globalscale": 10,  # global scale for whole plot, i.e., zooming in or out
    "faces": [],  # will contain all the triangles that make up the surface and shall be drawn
    "face_colourscale": [],  # will contain all the minimum and maximum z value that is associated with
    "spheres": [],  # will contain all the spheres
    "sphere_colours": [],  # will contain colours for all spheres
    "snap_count": 0,  # the counter for snapped images
    "snap_title": "snap",  # the title for the snapped images
    "keys": {},  # will contain all the keys pressed
    "firstrun": True,  # thether or not this is the first time the main function was executed
    "colours": [],  # color scheme - color values
    "borders": [],  # color scheme - transition values
    "high_contrast": False,  # select color scheme (True: red-black-blue, False: blue-turquiose-green-yellow-red)
    "povray": 0,  # scale for povray image and whether or not it's enabled
    "povray_data": None,  # the data povray requires to render this
    "povray_count": 0,  # how often povray images were rendered
    "povray_transform": "",  # transformation information for povray
    "graphical_output": True,  # whether graphical output is desired or not
    "window": 1,  # the number of the window used for the main display
    "savefile": None,  # the name of the file where visualization shall be saved to
    "savecount": 0,  # the number of visualization states already saved
    "additional": [],  # additional data. Only used to save the visualization
}

_keybindings = {
    hashstring("\033"): "quit",
    hashstring("+"): "zoom+",
    hashstring("="): "zoom+",
    hashstring("-"): "zoom-",
    hashstring("w"): "up",
    hashstring("s"): "down",
    hashstring("a"): "left",
    hashstring("d"): "right",
    hashstring("q"): "front",
    hashstring("e"): "back",
    hashstring("i"): "rot1+",
    hashstring("k"): "rot1-",
    hashstring("j"): "rot2+",
    hashstring("l"): "rot2-",
    hashstring("u"): "rot3+",
    hashstring("o"): "rot3-",
    hashstring("."): "snap",
    hashstring(","): "savevis",
    hashstring("p"): "povray",
}


class ParseError(Exception):
    """Raised when parsing a renderpath fails."""

    pass


def _ReSizeGLScene(Width, Height):
    global gl_c
    if Height == 0:  # Prevent A Divide By Zero If The Window Is Too Small
        Height = 1
    glViewport(
        0, 0, Width, Height
    )  # Reset The Current Viewport And Perspective Transformation
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60.0, float(Width) / float(Height), 0.01, 10000.0)
    glMatrixMode(GL_MODELVIEW)
    gl_c["resolution"] = (Width, Height)


def _expand_surface_data(data, face_indices):
    return [[data[i] for i in face] for face in face_indices]


# this function is repeatedly being called by OpenGL
def _main_control():
    global gl_c
    # draw everything
    GLAdjustCamera(gl_c["angles"], gl_c["translation"])
    if len(gl_c["faces"]) > 0:
        if gl_c["firstrun"]:
            DrawGLTrimesh(
                0,
                gl_c["faces"],
                gl_c["face_colourscale"],
                globalscale=gl_c["globalscale"],
                ccol=3,
                colours=gl_c["colours"],
                borders=gl_c["borders"],
            )
            gl_c["firstrun"] = False
        else:
            DrawGLTrimesh(0, globalscale=gl_c["globalscale"])
    if len(gl_c["spheres"]) > 0:
        DrawGLSpheres(
            gl_c["spheres"],
            (0, 1),
            globalscale=gl_c["globalscale"],
            sphere_elements=50,
            colour_list=gl_c["sphere_colours"],
        )
    glutSwapBuffers()


def _keyPressed(*args):
    global gl_c
    keys = gl_c["keys"]
    keyfunc = _keybindings.get(args[0], None)
    if keyfunc is not None:
        if keyfunc not in ("povray", "savevis"):
            keys[keyfunc] = True
        else:
            if keyfunc == "povray":
                if gl_c["povray"] > 0:
                    keys["povray"] = True
                else:
                    print(
                        "WARNING: PovRay support has either not been activated or is not supported by this type of visualization.",
                        file=sys.stderr,
                    )
            elif keyfunc == "savevis":
                if gl_c["savefile"] is not None:
                    gl_c["savecount"] += 1
                    SaveVisualizationState(
                        gl_c, gl_c["savefile"], prefix=str(gl_c["savecount"] - 1) + "_"
                    )
            else:
                raise Exception("Unhandled internal error")
        _evaluateKeyPressed()


def _evaluateKeyPressed():
    global gl_c
    keys = gl_c["keys"]
    # If escape is pressed, kill everything.
    if keys["quit"]:  # this is the escape sequence for the ESC key
        glutLeaveMainLoop()
        gl_c["running"] = False
    if keys["zoom+"]:
        gl_c["globalscale"] = gl_c["globalscale"] + 0.1
    if keys["zoom-"]:
        gl_c["globalscale"] = gl_c["globalscale"] - 0.1
    if keys["up"]:
        gl_c["translation"][1] += 1.0
    if keys["down"]:
        gl_c["translation"][1] -= 1.0
    if keys["left"]:
        gl_c["translation"][0] -= 1.0
    if keys["right"]:
        gl_c["translation"][0] += 1.0
    if keys["front"]:
        gl_c["translation"][2] -= 1.0
    if keys["back"]:
        gl_c["translation"][2] += 1.0
    if keys["rot1+"]:
        gl_c["angles"][0] -= 1.5
    if keys["rot1-"]:
        gl_c["angles"][0] += 1.5
    if keys["rot2+"]:
        gl_c["angles"][1] -= 1.5
    if keys["rot2-"]:
        gl_c["angles"][1] += 1.5
    if keys["rot3+"]:
        gl_c["angles"][2] += 1.5
    if keys["rot3-"]:
        gl_c["angles"][2] -= 1.5
    if keys["snap"]:
        snap(
            gl_c["resolution"],
            gl_c["snap_title"] + "_",
            "%3d",
            gl_c["snap_count"],
            "png",
        )
        gl_c["snap_count"] += 1
    if keys["povray"]:
        povray(
            [gl_c["povray"] * i for i in gl_c["resolution"]],
            gl_c["snap_title"] + "_",
            "%3d",
            gl_c["povray_count"],
            gl_c["angles"],
            [t - t0 for t, t0 in zip(gl_c["translation"], org_translation)],
            gl_c["povray_data"],
            gl_c["face_colourscale"],
            globalscale=gl_c["globalscale"],
            colours=gl_c["colours"],
            borders=gl_c["borders"],
            arrow_transform=gl_c["povray_transform"],
        )
        gl_c["povray_count"] += 1


def _keyReleased(*args):
    global gl_c
    keys = gl_c["keys"]
    keyfunc = _keybindings.get(args[0], None)
    if keyfunc is not None:
        if keyfunc not in ("savevis",):
            keys[keyfunc] = False
        else:
            if keyfunc == "savevis":
                pass
            else:
                raise Exception("Unhandled internal error")


def _initializeKeys(keys):
    for key in (
        "quit",
        "zoom+",
        "zoom-",
        "zoom-",
        "up",
        "down",
        "left",
        "right",
        "front",
        "back",
        "rot1+",
        "rot1-",
        "rot2+",
        "rot2-",
        "rot3+",
        "rot3-",
        "snap",
        "povray",
    ):
        keys[key] = False


def _TopLevelGlInitialization(
    gl_c, zoom, resolution, title="Molecule Visualization", use_light=False, hide=False
):
    _initializeKeys(gl_c["keys"])
    gl_c["running"] = True
    gl_c["globalscale"] *= zoom
    gl_c["resolution"] = resolution
    if not title == "Molecule Visualization":
        gl_c["snap_title"] = title
    if not gl_imported:
        return False
    try:
        os.environ["DISPLAY"]
    except KeyError:
        return False
    if len(os.environ["DISPLAY"]) == 0:
        return False
    if bool(glutInit):
        glutInit()
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
        # glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST) #this should enable antialiasing
        glutInitWindowSize(*gl_c["resolution"])
        glutInitWindowPosition(0, 0)
        gl_c["window"] = glutCreateWindow(title)
        # glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_CONTINUE_EXECUTION)
        glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS)
        InitGL(*gl_c["resolution"], use_light=use_light)
        glutDisplayFunc(_main_control)
        glutIdleFunc(_main_control)
        glutReshapeFunc(_ReSizeGLScene)
        glutKeyboardUpFunc(_keyReleased)
        glutKeyboardFunc(_keyPressed)
        if hide:
            glutHideWindow()
        return True
    else:
        return False


def _subcommandRangeNew(sc):
    s = sc[0].split("-")
    ranges = list(map(float, s))
    nsteps = int(sc[1])
    if nsteps > 0:
        return [[1.0 * r / nsteps for r in ranges]] * nsteps
    elif nsteps == 0:
        return []


def _parseTrajectory(trajectory):
    actions = trajectory.split(",")
    commands = [a.split("|") for a in actions]
    parsed = []
    for c in commands:
        p = []
        if re.match("(r[123][\+-]|t[123][\+-]|z[\+-])+", c[0]):
            nr_commands = sum(c[0].count(s) for s in ["r", "t", "z"])
            nr_ranges = c[1].count("-") + 1
            if not nr_ranges == nr_commands:
                raise ParseError(
                    "Could not parse "
                    + c[0]
                    + " and "
                    + c[1]
                    + " for an interlaced command."
                )
        for match in re.findall("r[123][\+-]|t[123][\+-]|z[\+-]", c[0]):
            subp = []
            if re.match("r", match):
                subp.append("angles")
            elif re.match("t", match):
                subp.append("translation")
            elif re.match("z", match):
                subp.append("globalscale")
            if re.match(".1", match):
                subp.append(0)
            elif re.match(".2", match):
                subp.append(1)
            elif re.match(".3", match):
                subp.append(2)
            else:
                subp.append(3)
            if re.match("(.|..)\+", match):
                subp.append("+")
            elif re.match("(.|..)-", match):
                subp.append("-")
            p.append(subp)
        for sc in _subcommandRangeNew(c[1].split("/")):
            parsed.append([p, sc])
    return parsed


def _set_high_contrast():
    if gl_c["povray"] > 0:
        middlecolour = [0.1, 0.1, 0.1]
    else:
        middlecolour = [0.0, 0.0, 0.0]
    sidecolours = ([0.3, 0.6, 1.0], [1.0, 0.4, 0.4])
    return [sidecolours[0], middlecolour, sidecolours[1]]


def _set_low_contrast():
    middlecolour = [0.0, 0.5, 0.0]
    sidecolours = ([0.75, 0.0, 0.0], [0.0, 0.0, 0.75])
    betweencolors = ([0.75, 0.75, 0.0], [0.0, 0.75, 0.75])
    return [
        sidecolours[0],
        betweencolors[0],
        middlecolour,
        betweencolors[1],
        sidecolours[1],
    ]
    # middlecolour=[0.2,0.2,0.2]
    # sidecolours=([0.0,0.2,1.0],[1.0,0.0,0.0])
    # return [sidecolours[0],middlecolour,sidecolours[1]]


def _svg_colorscale(filename, scale, high_contrast):
    svgfile = open(filename, "w")
    if high_contrast:
        s = r"""      <stop
         style="stop-color:#4d9aff;stop-opacity:1;"
         offset="0"
         id="stop3594" />
      <stop
         id="stop3629"
         offset="0.5"
         style="stop-color:#1a1a1a;stop-opacity:1;" />
      <stop
         style="stop-color:#ff6666;stop-opacity:1;"
         offset="1"
         id="stop3596" />
"""
    else:
        s = r"""      <stop
         style="stop-color:#C00000;stop-opacity:1;"
         offset="0"
         id="stop3594" />
      <stop
         id="stop3718"
         offset="0.25"
         style="stop-color:#C0C000;stop-opacity:1;" />
      <stop
         id="stop3629"
         offset="0.5"
         style="stop-color:#008000;stop-opacity:1;" />
      <stop
         id="stop3954"
         offset="0.75"
         style="stop-color:#00C0C0;stop-opacity:1;" />
      <stop
         style="stop-color:#0000C0;stop-opacity:1;"
         offset="1"
         id="stop3596" />
"""

    writeto(
        svgfile,
        r"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->

<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   width="259.62811"
   height="444.28571"
   id="svg2"
   version="1.1"
   inkscape:version="0.47 r22583"
   sodipodi:docname="New document 1">
  <defs
     id="defs4">
    <linearGradient
        id="linearGradient3592">
""",
    )
    writeto(svgfile, s)
    writeto(
        svgfile,
        """
</linearGradient>
    <inkscape:perspective
       sodipodi:type="inkscape:persp3d"
       inkscape:vp_x="0 : 526.18109 : 1"
       inkscape:vp_y="0 : 1000 : 0"
       inkscape:vp_z="744.09448 : 526.18109 : 1"
       inkscape:persp3d-origin="372.04724 : 350.78739 : 1"
       id="perspective10" />
    <linearGradient
       inkscape:collect="always"
       xlink:href="#linearGradient3592"
       id="linearGradient3598"
       x1="192.33765"
       y1="184.11545"
       x2="636.62335"
       y2="184.11545"
       gradientUnits="userSpaceOnUse"
       gradientTransform="matrix(0,-1,1,0,-304.31023,187.61329)" />
  </defs>
  <sodipodi:namedview
     id="base"
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1.0"
     inkscape:pageopacity="0.0"
     inkscape:pageshadow="2"
     inkscape:zoom="0.98994949"
     inkscape:cx="87.239069"
     inkscape:cy="250.3007"
     inkscape:document-units="px"
     inkscape:current-layer="layer1"
     showgrid="false"
     inkscape:window-width="1680"
     inkscape:window-height="973"
     inkscape:window-x="0"
     inkscape:window-y="26"
     inkscape:window-maximized="1">
    <inkscape:grid
       type="xygrid"
       id="grid2818"
       empspacing="5"
       visible="true"
       enabled="true"
       snapvisiblegridlinesonly="true" />
  </sodipodi:namedview>
  <metadata
     id="metadata7">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title></dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <g
     inkscape:label="Layer 1"
     inkscape:groupmode="layer"
     id="layer1"
     transform="translate(-49.999992,-2.3621826)">
    <rect
       style="fill:url(#linearGradient3598);fill-opacity:1;stroke:none"
       id="rect2816"
       width="40.38961"
       height="444.28571"
       x="-90.389603"
       y="-446.64789"
       transform="scale(-1,-1)" />
    <text
       xml:space="preserve"
       style="font-size:64px;font-style:normal;font-weight:normal;line-height:125%%;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;font-family:Helvetica;-inkscape-font-specification:Helvetica"
       x="96.045723"
       y="29.038183"
       id="text3600"
       sodipodi:linespacing="125%%"><tspan
         sodipodi:role="line"
         id="tspan3602"
         x="96.045723"
         y="29.038183"
         style="font-size:36px">%.4E</tspan></text>
    <text
       sodipodi:linespacing="125%%"
       id="text3604"
       y="445.81989"
       x="96.045723"
       style="font-size:64px;font-style:normal;font-weight:normal;line-height:125%%;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;font-family:Helvetica;-inkscape-font-specification:Helvetica"
       xml:space="preserve"><tspan
         style="font-size:36px"
         y="445.81989"
         x="96.045723"
         id="tspan3606"
         sodipodi:role="line">%.4E</tspan></text>
    <text
       sodipodi:linespacing="125%%"
       id="text3635"
       y="237.42903"
       x="96.045723"
       style="font-size:64px;font-style:normal;font-weight:normal;line-height:125%%;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;font-family:Helvetica;-inkscape-font-specification:Helvetica"
       xml:space="preserve"><tspan
         style="font-size:36px"
         y="237.42903"
         x="96.045723"
         id="tspan3637"
         sodipodi:role="line">0.0</tspan></text>
  </g>
</svg>"""
        % scale,
    )
    close(svgfile)


def SaveVisualizationState(obj, filename, prefix=""):
    """Save the current visualization state.

    This does not work when visualizing spheres.

    Args:
        obj: (dictionary) contains the data to be saved to disk
        filename: (string) the name of the file where the data shall be saved
        prefix: (string) path to the file (not including file)
    """
    oldstate = obj["firstrun"]
    obj["firstrun"] = True
    f = io.open(prefix + filename, "wb")
    # protocol version 2 stays compatible with Python 2 (but is slower than more recent versions)
    pickle.dump(obj, f, 2)
    f.close()
    obj["firstrun"] = oldstate
    print("Saved visualization state to file " + prefix + filename)


def LoadVisualization(filename):
    """Load the current visualization state.

    This does not work when visualizing spheres.

    Args:
        filename: (string) the name of the file from which the data shall be loaded
    """
    f = io.open(filename, "rb")
    obj = pickleload(f)
    f.close()
    for key in (
        "povray_data",
        "faces",
    ):
        if obj.get(key, None) is not None:
            if len(obj[key]) > 0:
                if not isinstance(obj[key][0], numpy.ndarray):
                    obj[key] = [numpy.array(a) for a in obj[key]]
    return obj


def _get_value_from_save(regex, dirs, key, fallback=None, warn=False):
    if len(regex) == 0 or len(key) == 0 or len(dirs) == 0:
        if warn:
            print(
                "WARNING: given key '%s' or regex '%s' or directory list %s is empty"
                % (key, regex, dirs),
                file=sys.stderr,
            )
        return fallback
    pattern = re.compile(regex)
    filenames = [
        d + os.sep + f
        for d in dirs.split("|")
        for f in os.listdir(d)
        if re.match(pattern, f)
    ]
    if len(filenames) == 0:
        if warn:
            print(
                "WARNING: no file found that matches the given regex: %s" % (regex),
                file=sys.stderr,
            )
        return fallback
    result = []
    for f in filenames:
        try:
            value = LoadVisualization(f)[key]
            result.append(value)
        except KeyError:
            if warn:
                print(
                    "WARNING: file %s matching pattern %s does not contain a dictionary with the necessary key."
                    % (f, regex),
                    file=sys.stderr,
                )
        except pickle.UnpicklingError:
            if warn:
                print(
                    "WARNING: file %s matching pattern %s could not be loaded (it probably was not saved by this programme)."
                    % (f, regex),
                    file=sys.stderr,
                )
    return result


def _TopLevelRenderFunction(gl_c, rendertrajectory):
    if re.match(".*,n(,d|,s|,p|,f)*$", rendertrajectory):
        snapping = False
    else:
        snapping = True
    if re.match(".*,d(,n|,s|,p|,f)*$", rendertrajectory):
        drop = True
    else:
        drop = False
    if re.match(".*,f(,n|,s|,p|,d)*$", rendertrajectory):
        first = False
    else:
        first = True
    povray_bool = False
    if gl_c["povray"] > 0:
        if re.match(".*,p(,d|,s|,n|,f)*$", rendertrajectory):
            povray_bool = True
    else:
        if re.match(".*,p(,d|,s|,n|,f)*$", rendertrajectory):
            print(
                "WARNING: PovRay support is not supported by this type of visualization.",
                file=sys.stderr,
            )
    save = False
    if re.match(".*,s(,n|,d|,p|,f)*$", rendertrajectory):
        if gl_c["savefile"] is not None:
            save = True
    while re.match(".*(,n|,d|,s|,p|,f)+$", rendertrajectory):
        rendertrajectory = rendertrajectory.rstrip(",n")
        rendertrajectory = rendertrajectory.rstrip(",d")
        rendertrajectory = rendertrajectory.rstrip(",s")
        rendertrajectory = rendertrajectory.rstrip(",p")
        rendertrajectory = rendertrajectory.rstrip(",f")
    actions = {
        "+": lambda a, b: a + b,
        "-": lambda a, b: a - b,
        "*": lambda a, b: a * b,
        "/": lambda a, b: a / b,
    }
    parsed = _parseTrajectory(rendertrajectory)
    digits = len(str(len(parsed) + 1))
    snapformat = "%+" + str(digits) + "d"
    if drop:
        while gl_c["running"]:
            glutMainLoopEvent()
            _main_control()
    _main_control()
    _main_control()
    if snapping:
        if first:
            snap(
                gl_c["resolution"],
                gl_c["snap_title"] + "_",
                "%3d",
                gl_c["snap_count"],
                "png",
            )
        gl_c["snap_count"] += 1
    if povray_bool:
        if first:
            povray(
                [gl_c["povray"] * i for i in gl_c["resolution"]],
                gl_c["snap_title"] + "_",
                "%3d",
                gl_c["povray_count"],
                gl_c["angles"],
                [t - t0 for t, t0 in zip(gl_c["translation"], org_translation)],
                gl_c["povray_data"],
                gl_c["face_colourscale"],
                globalscale=gl_c["globalscale"],
                colours=gl_c["colours"],
                borders=gl_c["borders"],
                arrow_transform=gl_c["povray_transform"],
            )
        gl_c["povray_count"] += 1
    if save:
        gl_c["savecount"] += 1
        if first:
            SaveVisualizationState(
                gl_c, gl_c["savefile"], prefix=str(gl_c["savecount"] - 1) + "_"
            )
    for ac in parsed:
        for com, val in zip(ac[0], ac[1]):
            if com[1] == 3:
                gl_c[com[0]] = actions[com[2]](gl_c[com[0]], val)
            else:
                gl_c[com[0]][com[1]] = actions[com[2]](gl_c[com[0]][com[1]], val)
        _main_control()
        if snapping:
            snap(
                gl_c["resolution"],
                gl_c["snap_title"] + "_",
                "%3d",
                gl_c["snap_count"],
                "png",
            )
            gl_c["snap_count"] += 1
        if save:
            gl_c["savecount"] += 1
            SaveVisualizationState(
                gl_c, gl_c["savefile"], prefix=str(gl_c["savecount"] - 1) + "_"
            )
        if povray_bool:
            povray(
                [gl_c["povray"] * i for i in gl_c["resolution"]],
                gl_c["snap_title"] + "_",
                "%3d",
                gl_c["povray_count"],
                gl_c["angles"],
                [t - t0 for t, t0 in zip(gl_c["translation"], org_translation)],
                gl_c["povray_data"],
                gl_c["face_colourscale"],
                globalscale=gl_c["globalscale"],
                colours=gl_c["colours"],
                borders=gl_c["borders"],
                arrow_transform=gl_c["povray_transform"],
            )
            gl_c["povray_count"] += 1


def RenderExtern(filename, agg=None, dictionary={}):
    """Restore a visualization state from a file. 

    It is also possible to update the visualization data using an aggregate or
    a dictionary.

    Args:
        filename: (string) the name of the file from which the data shall be loaded
        agg: (ManipulateAggregates.aggregate.agg) config data from this
            aggregate can be used to overwrite the loaded config
        dictionary: (dictionary) config data from this dictionary can be used
            to overwrite the loaded config. Takes precedence over agg.
    """
    global gl_c
    ext_gl_c = LoadVisualization(filename)
    gl_c.update(ext_gl_c)
    if agg is not None:
        vs = agg.get_vs
        agg.vs.update(dictionary)
    else:
        from ..aggregate import agg as agg_class

        heredict = copy.deepcopy(agg_class.default_vs)
        heredict.update(gl_c.get("vs", {}))
        heredict.update(dictionary)
        vs = lambda key: heredict.get(key, None)
    if vs("povray") > 0:
        if gl_c["povray_data"] is not None:
            gl_c["povray"] = vs("povray")
        else:
            gl_c["povray"] = 0
            print(
                "WARNING: no PovRay compatible data has been computed so PovRay support cannot be enabled.",
                file=sys.stderr,
            )
    if vs("high_contrast") is not None:
        if vs("high_contrast"):
            gl_c["colours"] = _set_high_contrast()
            gl_c["high_contrast"] = True
        else:
            gl_c["colours"] = _set_low_contrast()
            gl_c["high_contrast"] = False

    scales = None
    if (
        vs("colorscale") == "independent"
        or vs("colorscale") == ()
        or vs("colorscale") == ("", "")
    ):
        if gl_c["face_colourscale"][0] == gl_c["face_colourscale"][1]:
            print(
                "WARNING: the visualization saved in %s probably used a dependent color scale\n"
                "         (both borders of the scale are identical), hence, I cannot make it independent.\n"
                "         Will use the saved dependent scale instead.",
                file=sys.stderr,
            )
    elif vs("colorscale") == "dependent":
        abs_overall = max(
            [abs(gl_c["face_colourscale"][0]), abs(gl_c["face_colourscale"][1])]
        )
        gl_c["face_colourscale"] = (-abs_overall, abs_overall)
    elif isinstance(vs("colorscale"), (list, tuple)):
        if len(vs("colorscale")) >= 2:
            scale = vs("colorscale")
            scales = _get_value_from_save(
                scale[0], scale[1], "face_colourscale", warn=True
            )
    else:
        raise Exception("Unhandled internal exception: 'scale' has no matching value.")

    if scales is not None:
        gl_c["face_colourscale"] = (
            min(s[0] for s in scales),
            max(s[1] for s in scales),
        )
    if len(gl_c["colours"]) == 3:
        gl_c["borders"] = [
            0.0,
            -gl_c["face_colourscale"][0]
            / (gl_c["face_colourscale"][1] - gl_c["face_colourscale"][0]),
            1.0,
        ]
    elif len(gl_c["colours"]) == 5:
        zeroval = -gl_c["face_colourscale"][0] / (
            gl_c["face_colourscale"][1] - gl_c["face_colourscale"][0]
        )
        gl_c["borders"] = [0.0, zeroval / 2.0, zeroval, (zeroval + 1.0) / 2, 1.0]
    else:
        raise Exception("Unhandled internal exception: 'scales' has no matching value.")
    print("Colour scale: %.4E to %.4E" % gl_c["face_colourscale"])
    if vs("svgscale") is not None:
        _svg_colorscale(vs("svgscale"), gl_c["face_colourscale"], vs("high_contrast"))

    check = _TopLevelGlInitialization(
        gl_c,
        vs("zoom"),
        vs("resolution"),
        title=vs("title"),
        use_light=False,
        hide=vs("hide"),
    )

    if not check:
        print(
            "Cannot initialize OpenGL, will save visialization state.", file=sys.stderr
        )
        if vs("savefile") is not None:
            gl_c["savefile"] = vs("savefile")
        else:
            gl_c["savefile"] = "GLError"
        SaveVisualizationState(gl_c, "GLErrorDump_" + gl_c["savefile"])
        return
    if vs("savefile") is not None:
        gl_c["savefile"] = vs("savefile")
        if vs("savestart"):
            SaveVisualizationState(gl_c, "start_" + gl_c["savefile"])
    if vs("renderpath") is None:
        glutMainLoop()
    else:
        _TopLevelRenderFunction(gl_c, vs("renderpath"))
    if vs("savefile") is not None:
        if vs("saveend"):
            SaveVisualizationState(gl_c, "end_" + gl_c["savefile"])


def _PlotGL_Surface(agg, geom_manip_func=None, normal_manip_func=None):
    global gl_c
    vs = agg.get_vs

    corners, face_indices, normals = agg.get_surface()
    potential = agg.get_potential(corners)

    if vs("povray") > 0:
        povray_indices = numpy.array(face_indices)
        povray_vertices = corners
        povray_normals = numpy.array(normals)
        povray_potential = numpy.copy(potential)
        povray_potential.shape = (-1, 1)

    faces = numpy.array(_expand_surface_data(corners, face_indices))
    faces.shape = (-1, 3, 3)
    potential = numpy.array(_expand_surface_data(potential, face_indices))
    potential.shape = (-1, 3, 1)

    if geom_manip_func is not None:
        faces.shape = (-1, 3)
        faces = numpy.array([geom_manip_func(f) for f in faces])
        faces.shape = (-1, 3, 3)
        if vs("povray") > 0:
            povray_normals_end = numpy.array(
                [geom_manip_func(f) for f in povray_vertices + povray_normals]
            )
            povray_vertices = numpy.array([geom_manip_func(f) for f in povray_vertices])
            povray_normals = povray_normals_end - povray_vertices
            del povray_normals_end
            if normal_manip_func is not None:
                povray_normals = numpy.array(
                    [normal_manip_func(f) for f in povray_normals]
                )

    gl_c["faces"] = list(numpy.concatenate((faces, potential), axis=2))

    if vs("povray") > 0:
        gl_c["povray_data"] = [
            povray_indices,
            povray_vertices,
            povray_normals,
            povray_potential,
        ]

    scale = vs("colorscale")
    if scale == "independent" and (
        abs(numpy.min(potential)) <= 0.0 or abs(numpy.max(potential)) <= 0.0
    ):
        print(
            "WARNING: independent colour scaling won't work, will switch to dependent colour scaling.",
            file=sys.stderr,
        )
        scale = "dependent"

    if scale == "independent":
        gl_c["face_colourscale"] = (numpy.min(potential), numpy.max(potential))
    elif scale == "dependent":
        abs_overall = abs(max([abs(numpy.min(potential)), abs(numpy.max(potential))]))
        gl_c["face_colourscale"] = (-abs_overall, abs_overall)
    else:
        scales = _get_value_from_save(scale[0], scale[1], "face_colourscale", warn=True)
        if scales is not None:
            gl_c["face_colourscale"] = (
                min(s[0] for s in scales),
                max(s[1] for s in scales),
            )
        else:
            raise ValueError(
                "Scale must be either independent or dependent or an appropriate regex."
            )
    print("Colour scale: %.4E to %.4E" % gl_c["face_colourscale"])
    if vs("svgscale") is not None:
        _svg_colorscale(vs("svgscale"), gl_c["face_colourscale"], vs("high_contrast"))

    gl_c["povray"] = vs("povray")
    if vs("high_contrast"):
        gl_c["colours"] = _set_high_contrast()
        gl_c["high_contrast"] = True
    else:
        gl_c["colours"] = _set_low_contrast()
        gl_c["high_contrast"] = False
    # gl_c['colours']   =   [[0.0,0.0,1.0],[0.2,0.2,0.2],[1.0,0.0,0.0]]
    # gl_c['borders']   =   [0.0,-numpy.min(potential)/(numpy.max(potential)-numpy.min(potential)),1.0]
    # gl_c['borders']   =   [0.0,-gl_c['face_colourscale'][0]/(gl_c['face_colourscale'][1]-gl_c['face_colourscale'][0]),1.0]
    if len(gl_c["colours"]) == 3:
        gl_c["borders"] = [
            0.0,
            -gl_c["face_colourscale"][0]
            / (gl_c["face_colourscale"][1] - gl_c["face_colourscale"][0]),
            1.0,
        ]
    elif len(gl_c["colours"]) == 5:
        zeroval = -gl_c["face_colourscale"][0] / (
            gl_c["face_colourscale"][1] - gl_c["face_colourscale"][0]
        )
        gl_c["borders"] = [0.0, zeroval / 2.0, zeroval, (zeroval + 1.0) / 2, 1.0]
    else:
        raise Exception("Unhandled internal exception.")

    check = _TopLevelGlInitialization(
        gl_c,
        vs("zoom"),
        vs("resolution"),
        title=vs("title"),
        use_light=False,
        hide=vs("hide"),
    )

    if not check:
        print(
            "Cannot initialize OpenGL, will save visialization state.", file=sys.stderr
        )
        if vs("savefile") is not None:
            gl_c["savefile"] = vs("savefile")
        else:
            gl_c["savefile"] = "GLError"
        SaveVisualizationState(gl_c, "GLErrorDump_" + gl_c["savefile"])
        return
    if vs("savefile") is not None:
        gl_c["savefile"] = vs("savefile")
        if vs("savestart"):
            SaveVisualizationState(gl_c, "start_" + gl_c["savefile"])
    if vs("renderpath") is None:
        glutMainLoop()
    else:
        _TopLevelRenderFunction(gl_c, vs("renderpath"))
    if vs("savefile") is not None:
        if vs("saveend"):
            SaveVisualizationState(gl_c, "end_" + gl_c["savefile"])


def _PlotGL_Spheres(agg):
    global gl_c
    vs = agg.get_vs
    # get actual coordinates for indices
    coordinates = agg.get_coordinates()
    vdw_radii = agg.get_vdw_radii()
    gl_c["sphere_colours"] = agg.get_colours()
    spherescale = vs("vdw_scale")
    gl_c["spheres"] = [
        [c[0], c[1], c[2], r * spherescale] for c, r in zip(coordinates, vdw_radii)
    ]

    check = _TopLevelGlInitialization(
        gl_c,
        vs("zoom"),
        vs("resolution"),
        title=vs("title"),
        use_light=True,
        hide=vs("hide"),
    )

    if not check:
        print("Cannot initialize OpenGL, aborting.", file=sys.stderr)
        return
    if vs("renderpath") is None:
        glutMainLoop()
    else:
        _TopLevelRenderFunction(gl_c, vs("renderpath"))


def visualize(agg):
    """This function is a wrapper for visualizing the molecule using OpenGL.
    
    The aggregate's configuration dictionaries will be used.
    
    Args:
        agg: (ManipulateAggregates.aggregate.agg) the aggregate to be visualized

    Raises:
        ValueError.
    """
    vs = agg.get_vs
    gl_c["vs"] = copy.deepcopy(agg.vs)
    if vs("type").lower() == "vdw" or vs("type").lower() == "iso":
        if vs("align"):
            translate_before = -numpy.array(agg.get_center())
            translate_after = numpy.array(vs("align_center"))
            rotate = agg.get_align_matrix(vs("align_main3"), vs("align_main2"))
            rotmat = agg.get_povlight_matrix()
            geom_manip_func = (
                lambda e: numpy.dot(rotate, (numpy.array(e) + translate_before))
                + translate_after
            )
            normal_manip_func = lambda e: numpy.dot(rotmat, e)
            gl_c[
                "povray_transform"
            ] = """
      translate <%.10f,%.10f,%.10f>
      matrix <
          %.10f,%.10f,%.10f,
          %.10f,%.10f,%.10f,
          %.10f,%.10f,%.10f,
          %.10f,%.10f,%.10f
      >
      translate <%.10f,%.10f,%.10f>""" % (
                translate_before[0],
                translate_before[1],
                translate_before[2],
                rotate[0][0],
                rotate[1][0],
                rotate[2][0],
                rotate[0][1],
                rotate[1][1],
                rotate[2][1],
                rotate[0][2],
                rotate[1][2],
                rotate[2][2],
                0.0,
                0.0,
                0.0,
                translate_after[0],
                translate_after[1],
                translate_after[2],
            )
        else:
            geom_manip_func = None
            normal_manip_func = None
        _PlotGL_Surface(
            agg, geom_manip_func=geom_manip_func, normal_manip_func=normal_manip_func
        )
    elif vs("type").lower() == "simple":
        if vs("align"):
            newagg = agg.duplicate()
            newagg.align(vs("align_center"), vs("align_main3"), vs("align_main2"))
        else:
            newagg = agg
        _PlotGL_Spheres(newagg)
    else:
        raise ValueError("Selected method must be either complex or simple")
