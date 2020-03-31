"""Function definitions for drawing complex objects on the screen using OpenGL.

This is just a handy collection of wrapper functions for OpenGL. See
ManipulateAggregates.aggregate.visualize.PlotGL_Spheres and
ManipulateAggregates.aggregate.visualize.PlotGL_Surface as example
functions about how to use this module.

WARNING: Note that this module might be a bit slow. It was my first try at using OpenGL,
so please bear with me.
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
import logging

logger = logging.getLogger(__name__)

try:
    import numpy

    NPROTX = numpy.eye(4, dtype=float)
    NPROTY = numpy.eye(4, dtype=float)
    NPROTZ = numpy.eye(4, dtype=float)
except ImportError:
    logger.warning("Could not import numpy")
try:
    import FireDeamon.visualize as lvs

    use_cpp = True
except ImportError:
    logger.info("Could not import FireDeamon.visualize, will try outdated alternative")
    try:
        import libFDvisualize as lvs

        use_cpp = True
    except ImportError:
        logger.info(
            "Could also not import libFDvisualize, will use slower routines in pure Python."
        )
        use_cpp = False
try:
    import PIL.Image as Image
except ImportError:
    logger.info("Could not import PIL.Image.")
try:
    import re
except ImportError:
    logger.info("Could not import re.")

try:
    from ..p2p3IO import open, writeto, close
except ImportError:
    logger.warning("Could not import p2p3IO")

global COLOURS
# default RGB values for colourscale
COLOURS = [[0.0, 0.0, 0.0], [0.8, 0.3, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]]
global BORDERS
# default values on the scale from 0 to 1 where the colours in COLOURS shall be used
BORDERS = [0.0, 0.2, 0.7, 1.0]
global SHIFT
# default Cartesian displacement vector
SHIFT = [0.0, 0.0, 0.0]


def yield_values(
    values,
    minc=0,
    maxc=1,
    scale=1,
    maxextent_x=1,
    maxextent_y=1,
    xcol=0,
    ycol=1,
    zcol=2,
    ccol=2,
    shift=SHIFT,
    colours=COLOURS,
    borders=BORDERS,
    backcol=None,
    skip=None,
):
    """Scale objects that are to be drawn using OpenGL.
    
    Give back properly adjusted values to fit on the screen using a generator.

    Args:
        values: (list of lists) contains data in coloumn-major format
        minc: (float) minimum z value for colours (cbrange [minc:maxc] in gluplot)
        maxc: (float) maximum z value for colours (cbrange [minc:maxc] in gluplot)
        scale: (float) stretch z dimension by this factor
        xcol: (int) which coloumn of values contains x data (gnuplot: using
            ($xcol):($ycol):(scale*$zcol):($ccol))
        ycol: (int) which coloumn of values contains y data (gnuplot: using
            ($xcol):($ycol):(scale*$zcol):($ccol))
        zcol: (int) which coloumn of values contains z data (gnuplot: using
            ($xcol):($ycol):(scale*$zcol):($ccol))
        ccol: (int) which coloumn of values contains color data (gnuplot: using
            ($xcol):($ycol):(scale*$zcol):($ccol))
        shift: (list of 3 floats) Cartesian displacement vector
        maxextent_x: (float) if negative, desired maximum extent of structure
            in x direction, will be streched to fit perfectly. If positive,
            scale x direction by the absolute value of the given number.
        maxextent_y: (float) if negative, desired maximum extent of structure
            in y direction, will be streched to fit perfectly. If positive,
            scale y direction by the absolute value of the given number
        colours: (list of lists of 3 floats) colours for linear interpolation
            in rgb format
        borders: (list of floats) start and stop of colour intervals in
            fractions of 1 ([minc:maxc] will be mapped to [0:1])
        backcol: (int) if you want an additional coloumn returned, specify
            which one
        skip: (list of lists of 3 floats) should be a list of (m,n,o). Starting
            at point o, skip n data points every m data points. Is processed in
            the given order.

    Yields:
        tuples of 2 or 3 lists. If backcol is None, only 2 elements are
        present. Otherwise the third element is the coloumn indexed by
        backcol. The first list contains 3-element lists which are the
        Cartesian coordinates of the vertices to be drawn. The second list are
        the colour codes for OpenGL for these vertices.
    """
    # transform to numpy arrays
    colours = numpy.array(colours)
    borders = numpy.array(borders)
    # c_array contains all the values which we will map on a colour (here: 3rd coloumn)
    # Values are rescaled to fit in [0,1]
    # reshaping is necessary because otherwise numpy does not fill the shape tuple which
    # would give errors otherwise
    c_array = (numpy.array([p[ccol] for p in values]).reshape((-1, 1)) - minc) / (
        1.0 * (maxc - minc)
    )
    # m stands for matrix
    m = numpy.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    if maxextent_x < 0:
        # find maximum extent
        mx = max([abs(p[xcol]) for p in values])
        # rescale position of points to maxextent_x and maxextent_y and stretch
        # z-direction
        m[0][0] = -maxextent_x / (2.0 * mx)
    else:
        m[0][0] = maxextent_x
    if maxextent_y <= 0:
        # find maximum extent
        my = max([abs(p[ycol]) for p in values])
        # rescale position of points to maxextent_x and maxextent_y and stretch
        # z-direction
        m[1][1] = -maxextent_y / (2.0 * my)
    else:
        m[1][1] = maxextent_y
    m[2][2] = scale
    # v_array contains the points properly scaled and shifted
    # reshaping is necessary because otherwise numpy does not fill the shape tuple which
    # would give errors otherwise
    v_array = numpy.array([[v[xcol], v[ycol], v[zcol]] for v in values]).reshape(
        (-1, 3)
    ).dot(m) + numpy.array(shift)
    # b_array contains either nothing (i.e. Nones) or the values from the coloumn that
    # shall also be passed back
    if not backcol is None:
        b_array = numpy.array([b[backcol] for b in values]).reshape((-1, 1))
    else:
        b_array = numpy.array([None] * len(values)).reshape((-1, 1))
    nr_borders = len(borders)
    # process skipping
    if not skip is None:
        for m, n, o in skip:
            startindex = o
            new_c_array = numpy.array([]).reshape(0, c_array.shape[1])
            new_b_array = numpy.array([]).reshape(0, b_array.shape[1])
            new_v_array = numpy.array([]).reshape(0, v_array.shape[1])
            if m > 0 and n > 0:
                while startindex < len(c_array):
                    endindex = startindex + m
                    new_c_array = numpy.concatenate(
                        (new_c_array, c_array[startindex:endindex]), axis=0
                    )
                    new_b_array = numpy.concatenate(
                        (new_b_array, b_array[startindex:endindex]), axis=0
                    )
                    new_v_array = numpy.concatenate(
                        (new_v_array, v_array[startindex:endindex]), axis=0
                    )
                    startindex = endindex + n
            c_array = new_c_array
            b_array = new_b_array
            v_array = new_v_array
    # c: will contain colour code, now it contains the value that we want to assign a
    # colour code
    # v: cotains x,y and z coordinate of point as it will be drawn on the screen
    # b: contains value from coloumn that was requested to be passed back, if so desired
    # without the [] around c and b, a list of length 1 would be returned
    for [c], v, [b] in zip(c_array, v_array, b_array):
        j = -1
        # make sure we're never out of the interval [0,1]
        if c <= 0.0:
            c = 0
            j = 1
        elif c >= 1.0:
            c = 1.0
            j = nr_borders - 1
        else:
            # find the range defined by the borders that we are in
            for i in range(nr_borders):
                if c < borders[i]:
                    j = i
                    break
        # perform linear interpolation
        c1 = colours[j - 1]
        c2 = colours[j]
        b1 = borders[j - 1]
        b2 = borders[j]
        fraction = 1.0 * (c - b1) / (b2 - b1)
        c = fraction * c2 + (1.0 - fraction) * c1
        if b is None:
            yield (list(c), list(v))
        else:
            yield (list(c), list(v), b)


# opengl stuff is blatantly taken from the tutorial at
# http://pydoc.net/Python/PyOpenGL-Demo/3.0.0/PyOpenGL-Demo.NeHe.lesson5/
try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
except ImportError:
    logger.warning("Could not import OpenGL")

# Number of the glut window.
# A general OpenGL initialization function.  Sets all of the initial parameters.
def InitGL(
    Width, Height, use_light=False
):  # We call this right after our OpenGL window is created.
    """Initialize OpenGL.

    Args:
        Width: (int) width of window in pixels
        Height: (int) height of window in pixels
        use_light: (bool) whether or not to use lighting.

    Bug: 
        Nothing is displayed if use_light is True and a trimesh shall be
        displayed using ManipulateAggregates.collection.opengl.WritePovrayTrimesh
    """
    glClearColor(1.0, 1.0, 1.0, 0.0)  # This Will Clear The Background Color To White
    glClearDepth(1.0)  # Enables Clearing Of The Depth Buffer
    glDepthFunc(GL_LESS)  # The Type Of Depth Test To Do
    glEnable(GL_DEPTH_TEST)  # Enables Depth Testing
    glShadeModel(GL_SMOOTH)  # Enables Smooth Color Shading
    if use_light:
        specular = [0.1, 0.1, 0.1, 0.1]
        shininess = [30.0]
        position = [1.0, 0.0, 5.0, 0.0]
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular)
        glMaterialfv(GL_FRONT, GL_SHININESS, shininess)
        glEnable(GL_LIGHT0)
        glLightfv(GL_LIGHT0, GL_POSITION, position)
        glEnable(GL_COLOR_MATERIAL)
        glEnable(GL_LIGHTING)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()  # Reset The Projection Matrix
    # Calculate The Aspect Ratio Of The Window
    gluPerspective(60.0, float(Width) / float(Height), 0.01, 10000.0)
    glMatrixMode(GL_MODELVIEW)
    return


def GetCameraMatrix(angles, invert=True):
    """Get the proper transformation matrix for the camera.

    Args:
        angles: (tuple of 3 floats) angles around the x, y and z axes
            in degrees
        invert: (bool) if True, the matrices Mx, My and Mz (which are the
            rotation matrices around the x, y and z axes, respectively) will be
            multiplied in the order Mx*My*Mz. If False, the order will be inversed.

    Returns:
        A 4x4 numpy matrix that can be used by, say, PoVRay to correctly
        manipulate the camera.
    """
    # rotate the view properly
    x, y, z = angles
    x *= 0.017453292519943295
    y *= 0.017453292519943295
    z *= 0.017453292519943295
    if invert:
        x = -x
        y = -y
        z = -z
    NPROTX[1, 1] = numpy.cos(x)
    NPROTX[1, 2] = -numpy.sin(x)
    NPROTX[2, 1] = numpy.sin(x)
    NPROTX[2, 2] = numpy.cos(x)
    NPROTY[0, 0] = numpy.cos(y)
    NPROTY[0, 2] = numpy.sin(y)
    NPROTY[2, 0] = -numpy.sin(y)
    NPROTY[2, 2] = numpy.cos(y)
    NPROTZ[0, 0] = numpy.cos(z)
    NPROTZ[0, 1] = -numpy.sin(z)
    NPROTZ[1, 0] = numpy.sin(z)
    NPROTZ[1, 1] = numpy.cos(z)
    if invert:
        return numpy.dot(numpy.dot(NPROTX, NPROTY), NPROTZ)
    else:
        return numpy.dot(numpy.dot(NPROTZ, NPROTY), NPROTX)


# This function is called to properly adjust the relative positions of plot and camers
def GLAdjustCamera(angles, translation):
    """Properly adjust relative positions of plot and camera.

    Args:
        angles: (tuple of 3 floats) angles around the x, y and z axes
            in degrees
        translation: (tuple of 3 floats) camera displacement vector
    """
    glClear(
        GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT
    )  # Clear The Screen And The Depth Buffer
    glLoadIdentity()  # Reset The View
    glTranslatef(*translation)  # move the camera to the correct position
    glMultMatrixd(numpy.ndarray.flatten(GetCameraMatrix(angles)))
    return


global TRIMESHES
# static list that stores the trimeshes that were already drawn
TRIMESHES = {}

# This function is called to display the surface on the screen
def DrawGLTrimesh(
    index,
    faces=None,
    colourscale=None,
    globalscale=None,
    globalskip=0,
    elements_per_line=None,
    ccol=2,
    colours=COLOURS,
    borders=BORDERS,
):
    """This function tells OpenGL to draw a trimesh.

    The function uses a static variable to remember the trimeshes that it was
    used to draw in order to save a lot of time. The only input that is allowed
    to change for a trimesh that has already been drawn is globalscale. In
    fact, everything else will be ignored for subsequent calls. The parameters
    faces and colourscale are mandatory when calling this function with a
    new trimesh for the first time. A new trimesh is defined as one whose
    index has not yet been passed to this function.

    Args:
        index: (int) the index of the trimesh to be drawn. Can also be an new index.
        faces: (list of lists) should look like [A,B,C,...] where A,B and C are
            faces and have 3 elements each p1,p2 and p3, which all are Cartesian
            vectors in 3 dimensions
        colourscale: (tuple of 2 floats) the z-values that should be associated
            with the "lowest" colour and the "highest" colour
        globalscale: (float) scale the whole plot in all dimensions by this much
        globalskip: (int) every 1 triangle, do not draw the next this many
            triangles and every "line" of triangles, do not plot the next this
            many lines
        elements_per_line: (int) assuming a quadratic mesh (in x and y
            dimensions), how many points are there per line. Ignored if not set
        ccol: (int) which entry in each sublist of faces contains color data
            (a float between 0 and 1 that is mapped to the colour scale using
            the values in colours and borders.
        colours: (list of N lists of 3 floats) a colour scale for the plot in RGB format
        borders: (list of N floats) the borders for the new colour scale
    """
    global TRIMESHES

    if globalskip > 0:
        if elements_per_line is None:
            skip = [(3, 3 * globalskip, 0)]
        else:
            el = elements_per_line - 1
            # per space between 2 data points (there are el such spaces per line), there
            # are 6 points in the below array, 3 per triangle
            skip = [(6 * el, 6 * el * globalskip, 0), (3, 3 * globalskip, 0)]
    else:
        skip = None

    if not index in TRIMESHES:
        new_mesh = list(
            yield_values(
                [point for triangle in faces for point in triangle],
                minc=colourscale[0],
                maxc=colourscale[1],
                scale=1.0 * globalscale,
                skip=skip,
                maxextent_x=1.0 * globalscale,
                maxextent_y=1.0 * globalscale,
                ccol=ccol,
                colours=colours,
                borders=borders,
            )
        )
        if use_cpp:
            cpp_mesh = lvs.add_mesh(new_mesh)
            TRIMESHES[index] = {
                "mesh": cpp_mesh,
                "scale": globalscale,
            }
        else:
            TRIMESHES[index] = {
                "mesh": new_mesh,
                "scale": globalscale,
            }

    orgscale = TRIMESHES[index]["scale"]
    rescale = globalscale / orgscale

    if use_cpp:
        lvs.visualize_mesh(TRIMESHES[index]["mesh"], rescale)
    else:
        glBegin(GL_TRIANGLES)
        # get variables via a generator
        # c stands for colour and p stands for point, i.e. position
        for c, p in TRIMESHES[index]["mesh"]:
            glColor3f(c[0], c[1], c[2])
            glVertex3f(rescale * p[0], rescale * p[1], rescale * p[2])
        glEnd()
    return


# This function is called to translate the OpenGL data to povray format and write it to
# a file handle
def WritePovrayTrimesh(
    handle,
    matrix,
    translation,
    indices,
    points,
    normals,
    colorvalues,
    colourscale,
    globalscale=1,
    globalskip=0,
    elements_per_line=None,
    ccol=2,
    colours=COLOURS,
    borders=BORDERS,
):
    """This function writes a trimesh in PovRay format to a file handle.

    Args:
        handle: (file descriptor) the handle whose .write method will be used
        matrix: (4x4 matrix) the matrix that rotates the mesh to be properly visualized
        translation: (list of 3 floats) a Cartesian displacement vector that
            will be applied to all vertices
        indices: (list of lists of 3 ints) each sublist contains the indices of
            those vertives (given in points) that make up a face of the trimesh
        points: (numpy array of shape 1,3 and dtype float) the vertices
            describing the actual trimesh per element
        normals: (list of [float,float,float]) one normal vector per vertex
        colorvalues: (numpy array of shape 1,3 and dtype float) the RGB color
            values associated with each vertex (as given in points)
        ccol: (int) which entry in each sublist of faces contains color data
            (a float between 0 and 1 that is mapped to the colour scale using
            the values in colours and borders.
        colourscale: (2 element tuple) the z-values that should be associated
            with the "lowest" colour and the "highest" colour
        globalscale: (float) scale the whole plot in all dimensions by this much
        globalskip: (int) every 1 triangle, do not draw the next this many
            triangles and every "line" of triangles, do not plot the next this
            many lines
        elements_per_line: (int) assuming a quadratic mesh (in x and y
            dimensions), how many points are there per line. Ignored if not set
        colours: (list of N lists of 3 floats) a colour scale for the plot in RGB format
        borders: (list of N floats) the borders for the new colour scale
    """
    transparency = 0.0
    if len(points) != len(normals) or len(points) != len(colorvalues):
        raise ValueError(
            "Length of 'points', 'normals' and 'colorvalues' are not identical."
        )
    scale = 0.008
    scalemat = numpy.dot(
        numpy.array(
            [[scale, 0, 0, 0], [0, scale, 0, 0], [0, 0, scale, 0], [0, 0, 0, scale]],
            dtype=float,
        ),
        matrix,
    )
    scalemat = numpy.ndarray.flatten(scalemat)
    for v, c, n in zip(
        scalemat,
        range(1, 17),
        (
            "M11",
            "M12",
            "M13",
            "M14",
            "M21",
            "M22",
            "M23",
            "M24",
            "M31",
            "M32",
            "M33",
            "M34",
            "M41",
            "M42",
            "M43",
            "M44",
        ),
    ):
        if c % 4 != 0:
            writeto(handle, "#declare %s=%.10f;\n" % (n, v))
    points_colors = numpy.concatenate((points, colorvalues), axis=1)
    tab = "    "
    tabcount = 0
    writeto(handle, tab * tabcount + "mesh2 {\n")
    tabcount += 1
    writeto(handle, tab * tabcount + "vertex_vectors {\n")
    tabcount += 1
    writeto(handle, tab * tabcount + "%d,\n" % (len(points)))
    it = list(
        yield_values(
            points_colors,
            minc=colourscale[0],
            maxc=colourscale[1],
            scale=1.0 * globalscale,
            maxextent_x=1.0 * globalscale,
            maxextent_y=1.0 * globalscale,
            ccol=ccol,
            colours=colours,
            borders=borders,
        )
    )
    for c, p in it:
        writeto(handle, tab * tabcount + "<%.10f,%.10f,%.10f>,\n" % tuple(p))
    tabcount -= 1
    writeto(handle, tab * tabcount + "}\n")
    writeto(handle, tab * tabcount + "normal_vectors {\n")
    tabcount += 1
    writeto(handle, tab * tabcount + "%d,\n" % (len(normals)))
    for n in normals:
        writeto(handle, tab * tabcount + "<%.10f,%.10f,%.10f>,\n" % tuple(n))
    tabcount -= 1
    writeto(handle, tab * tabcount + "}\n")
    writeto(handle, tab * tabcount + "texture_list {\n")
    tabcount += 1
    writeto(handle, tab * tabcount + "%d,\n" % (len(colorvalues)))
    for c, p in it:
        writeto(handle, tab * tabcount + "RGBTVERT(<%.6f,%.6f,%.6f" % tuple(c))
        # writeto(handle,",1.0-((1.0-%.6f)*OPAQUE)>),\n"%(transparency))
        writeto(handle, ",1.0-OPAQUE>),\n")
    tabcount -= 1
    writeto(handle, tab * tabcount + "}\n")
    writeto(handle, tab * tabcount + "face_indices {\n")
    tabcount += 1
    writeto(handle, tab * tabcount + "%d\n" % (len(indices)))
    for i in indices:
        writeto(handle, tab * tabcount + "<%d,%d,%d>," % tuple(i))
        writeto(handle, "%d,%d,%d\n" % tuple(i))
    tabcount -= 1
    writeto(handle, tab * tabcount + "}\n")
    writeto(handle, tab * tabcount + "inside_vector <0, 0, 1>\n")
    writeto(handle, tab * tabcount + "no_shadow\n")
    writeto(handle, tab * tabcount + "matrix <\n")
    tabcount += 1
    writeto(handle, tab * tabcount + "M11,M12,M13,\n")
    writeto(handle, tab * tabcount + "M21,M22,M23,\n")
    writeto(handle, tab * tabcount + "M31,M32,M33,\n")
    writeto(handle, tab * tabcount + "M41,M42,M43\n")
    tabcount -= 1
    writeto(handle, tab * tabcount + ">\n")
    # writeto(handle,tab*tabcount+"translate <%.10f*M11, %.10f*M22, %.10f*M33>\n"%(translation[0],translation[1],-translation[2]))
    tabcount -= 1
    writeto(handle, tab * tabcount + "}\n")
    return


# This function is called to display the spheres on the screen
def DrawGLSpheres(
    spheres,
    colourscale,
    globalscale=1,
    globalskip=0,
    elements_per_line=None,
    sphere_elements=50,
    colour_list=None,
):
    """This function tells OpenGL to draw spheres via GLUT.

    Args:
        spheres: (list of lists of 4 floats x,y,z,r) x,y,z are sphere's Cartesian
            center coordinates and r is the radius
        colourscale: (2 element tuple) the z-values that should be associated
            with the "lowest" colour and the "highest" colour
        globalscale: (float) scale the whole plot in all dimensions by this much
        globalskip: (int) every 1 triangle, do not draw the next this many
            triangles and every "line" of triangles, do not plot the next this
            many lines
        elements_per_line: (int) assuming a quadratic mesh (in x and y
            dimensions), how many points are there per line. Ignored if not set
        sphere_elements: (int) controls how many azimutal and longitudinal
            elements are drawn per sphere
        colour_list: (list of lists of 3 floats) RGB values for sphere colors.
            If None is given, the spheres are coloured according to the
            z-coordinate of their centers.
    """

    if globalskip > 0:
        if elements_per_line is None:
            skip = [(1, globalskip, 0)]
        else:
            # per space between 2 data points (there are el such spaces per line), there
            # are 6 points in the below array, 3 per triangle
            skip = [
                (elements_per_line, elements_per_line * globalskip, 0),
                (1, globalskip, 0),
            ]
        skip = [(1, globalskip, 0)]
    else:
        skip = None
    use_gen = not (colour_list is None)
    if use_gen:
        gen = (c for c in colour_list)
    # c stands for colour and p stands for point, i.e. position and r stands for radius
    for c, p, r in yield_values(
        spheres,
        minc=colourscale[0],
        maxc=colourscale[1],
        scale=1.0 * globalscale,
        colours=[[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        borders=[0.0, 1.0],
        backcol=3,
        skip=skip,
        maxextent_x=1.0 * globalscale,
        maxextent_y=1.0 * globalscale,
    ):
        p = numpy.array(p)
        glTranslatef(*p)
        if use_gen:
            glColor3f(*next(gen))
        else:
            glColor3f(*c)
        glutSolidSphere(globalscale * r, sphere_elements, sphere_elements)
        glTranslatef(*(-p))
    return


# This function glues together all the other displaying functions
def GLMainDisplay(
    angles,
    translation,
    faces,
    face_colourscale,
    draw_faces,
    spheres,
    sphere_colourscale,
    draw_spheres,
    globalscale,
    elements_per_line_faces,
    elements_per_line_spheres,
    globalskip,
):

    # adjust camera
    GLAdjustCamera(angles, translation)
    # draw what is desired
    if draw_faces:
        DrawGLTrimesh(
            faces,
            face_colourscale,
            globalscale=globalscale,
            globalskip=globalskip,
            elements_per_line=elements_per_line_faces,
        )
    if draw_spheres:
        DrawGLSpheres(
            spheres,
            sphere_colourscale,
            globalscale=globalscale,
            globalskip=globalskip,
            elements_per_line=elements_per_line_spheres,
        )
    #  since this is double buffered, swap the buffers to display what just got drawn.
    glutSwapBuffers()


def snap(size, basename, format, count, extension):
    """Save an OpenGL screenshot to disk.

    Args:
        size: (tuple of 2 ints) size of the image in pixels in x and y directions
        basename: (string) path plus first part of filename
        count: (int) if the images shall be numbered, say what number this
            image shall have
        format: (string) to have fixed width numbering, declare the printf-type
            format string (like %6d)
        extension: (string) filetype, "png" recommended
    """
    filename = re.sub("\s", "0", basename + format % (count) + "." + extension)
    screenshot = glReadPixels(0, 0, size[0], size[1], GL_RGBA, GL_UNSIGNED_BYTE)
    snapshot = Image.frombuffer("RGBA", size, screenshot, "raw", "RGBA", 0, 0)
    snapshot.save(filename)
    print(filename)
    return count + 1


def povray(
    size,
    basename,
    format,
    count,
    angles,
    translation,
    povray_data,
    colourscale,
    globalscale=1,
    colours=COLOURS,
    borders=BORDERS,
    arrow_transform="",
):
    """Render the current OpenGL-view of a trimesh using PoVRay into a PNG.

    The PNGs that are generated all be called basename followed by an
    integer number (the index).

    Args:
        size: (tuple of 2 ints) the x and y sizes of the plot in pixels
        basename: (string) path to the file without extension
        format: (string) printf-type format string converting an integer to the
            count in string form (i.e., the index)
        count: (int) index of the image
        angles: (tuple of 3 floats) angles around the x, y and z axes
            in degrees
        translation: (tuple of 3 floats) camera displacement vector
        povray_data: (list of lists) the arguments will be passed as "indices",
            "points", "normals" and "colorvalues" to the
            function ManipulateAggregates.collection.opengl.WritePovrayTrimesh
            in that order
        colourscale: (tuple of 2 floats) the z-values that should be associated
            with the "lowest" colour and the "highest" colour
        globalscale: (float) scale the whole plot in all dimensions by this much
        colours: (list of lists of 3 floats) colours for linear interpolation
            in rgb format
        borders: (list of floats) start and stop of colour intervals in
            fractions of 1 ([minc:maxc] will be mapped to [0:1])
        arrow_transform: (string) valid PoVRay code that will be used to
            transform any arrows that the user might want to add using the
            provided "arrow" macro by editing the generated .pov-file
    """
    LEFTMAT = numpy.array(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]], dtype=float
    )
    extension = "pov"
    filename = re.sub(
        "\s",
        "0",
        basename + "%dx%d_" % (size[0], size[1]) + format % (count) + "." + extension,
    )
    handle = open(filename, "w")
    viewmat = GetCameraMatrix(angles, invert=False)
    writeto(handle, "//The command used to render this thing:\n")
    writeto(
        handle,
        "//povray +W%d +H%d -I%s -O%s +UA +D +X +A +FN\n"
        % (size[0], size[1], filename, filename + ".png"),
    )
    writeto(
        handle,
        """
#version 3.5;
#if (version < 3.5)
#error "This programme was designed to work with povray 3.5 or above."
#end
#macro RGBTVERT ( C1 )
  texture { pigment { rgbt C1 }}
#end
#declare OPAQUE=1.0;
//If you want to draw an arrow starting at (0,0,100) pointing in the direction (0.651156,-0.985225,1.302909)
//with width 1 for the line and 3 times that for the base of the cone of colour blue and the length scaled by 20 and
//the arrowhead starting 20%% before the end of the arrow, do it like so
//below the declarations of M11 through M43
//arrow(<0.0,0.0,100.0>,<0.651156,-0.985225,1.302909>,1.0,rgbt<0.000,0.000,1.000,0.000>,20,3,0.2,M11,M12,M13,M21,M22,M23,M31,M32,M33)
#default { texture {
    finish { ambient 0.800 diffuse 0.200 phong 0.2 phong_size 4.0 specular 0.05 roughness 0.10 }
} }
""",
    )
    # WritePovrayTrimesh(handle, numpy.dot(viewmat,LEFTMAT).T, povray_data[0], povray_data[1], povray_data[2], povray_data[3],
    #        colourscale, globalscale=globalscale, ccol=3, colours=colours, borders=borders)
    WritePovrayTrimesh(
        handle,
        numpy.dot(LEFTMAT, viewmat).T,
        translation,
        povray_data[0],
        povray_data[1],
        povray_data[2],
        povray_data[3],
        colourscale,
        globalscale=globalscale,
        ccol=3,
        colours=colours,
        borders=borders,
    )
    writeto(
        handle,
        """
camera {
  up <0, %.10f, 0>
  right <%.10f, 0, 0>
  location <0.0000, 0.0000, -2.0000>
  look_at <0.0000, 0.0000, 0.0000>
  direction <0.0000, 0.0000, 1.0000>
  translate <-(%.10f)*M11, -(%.10f)*M22, -(%.10f)*M33>
}
union {
    light_source {
      <0.0000, 0.0000, -2.0000>
      color rgb<1.000, 1.000, 1.000>
    }
    light_source {
      <10.0000, 10.0000, -2.0000>
        color rgb<0.500, 0.500, 0.500>
    }
    light_source {
        <10.0000, -10.0000, -2.0000>
        color rgb<0.500, 0.500, 0.500>
    }
    light_source {
        <-10.0000, 10.0000, -2.0000>
        color rgb<0.500, 0.500, 0.500>
    }
    light_source {
        <-10.0000, -10.0000, -2.0000>
        color rgb<0.500, 0.500, 0.500>
    }
    light_source {
      <5.0000, 5.0000, -2.0000>
        color rgb<0.500, 0.500, 0.500>
    }
    light_source {
        <5.0000, -5.0000, -2.0000>
        color rgb<0.500, 0.500, 0.500>
    }
    light_source {
        <-5.0000, 5.0000, -2.0000>
        color rgb<0.500, 0.500, 0.500>
    }
    light_source {
        <-5.0000, -5.0000, -2.0000>
        color rgb<0.500, 0.500, 0.500>
    }
    translate <-(%.10f)*M11, -(%.10f)*M22, -(%.10f)*M33>
}
background {
    color rgb<1.000, 1.000, 1.000>
}
#macro arrow (P, D, R, C, L, CF, CL, M11,M12,M13,M21,M22,M23,M31,M32,M33)
  #local T = texture { pigment { C } finish { ambient 0.800 diffuse 0.200 phong 0.3 phong_size 2.0 specular 0.05 roughness 0.10 } }
  #local Ecyl = P+((1.0-CL)*L*D);
  #local Econ = P+(L*D);
  union {
        object { cylinder {P, Ecyl, R open texture {T} no_shadow} }
        object { cone{Ecyl, CF*R, Econ, 0.0 open texture {T} no_shadow} }
        object { disc{ P, D, R texture {T} no_shadow} }
        object { disc{ Ecyl, D, CF*R texture {T} no_shadow} }
        %s
        matrix <
        M11, M12, M13,
        M21, M22, M23,
        M31, M32, M33
        0.0, 0.0, 0.0
        >
  }
#end
"""
        % (
            1.0,
            1.0 * size[0] / size[1],
            translation[0],
            translation[1],
            translation[2],
            translation[0],
            translation[1],
            translation[2],
            arrow_transform,
        ),
    )
    close(handle)
    from subprocess import Popen

    f = Popen(
        [
            "povray",
            "+W%d" % (size[0]),
            "+H%d" % (size[1]),
            "-I%s" % (filename),
            "-O%s" % (filename + ".png"),
            "+UA",
            "+D",
            "+X",
            "+A",
            "+FN",
        ]
    )
    f.wait()
    print(filename + ".png")
    return count + 1
