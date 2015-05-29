"""
function definitions for drawing stuff on the screen using OpenGL
This is just a handy collection of wrapper functions for opengl
"""

import numpy as np

def yield_values(values,minc=0,maxc=1,scale=1,maxextent_x=1,maxextent_y=1,xcol=0,ycol=1,zcol=2,ccol=2,shift=[0.0,0.0,0.0],colours=[[0.0,0.0,0.0],[0.8,0.3,0.0],[1.0,1.0,0.0],[1.0,1.0,1.0]],borders=[0.0,0.2,0.7,1.0],backcol=None,skip=None):
    """
    Give back properly adjusted values to fit on the screen using a generator
    values: contains the data in coloumn major format
    minc: minimum z value for colours (cbrange [minc:maxc] in gluplot)
    maxc: maximum z value for colours (cbrange [minc:maxc] in gluplot)
    scale: stretch z dimension by this factor (to change apparent height of plot, see next line)
    xcol, ycol, zcol, ccol: which coloumn of values contains x,y and z data (gnuplot: using ($xcol):($ycol):(scale*$zcol):($ccol))
    maxextent_{x,y}: if negative: desired maximum extent of structure in {x,y} direction, will be streched to fit perfectly
                     if positive: scale x and y direction by the absolute value of the given number
    colours: colours for linear interpolation in rgb format
    borders: start and stop of colour intervals in fractions of 1 ([minc:maxc] will be mapped to [0:1])
    backcol: if you want an additional coloumn returned, specify which one (only scalars are possible)
    skip: should be a list of (m,n,o): starting at point o, skip n data points every m data points. Is processed in the given order.
    """
    #transform to numpy arrays
    colours=np.array(colours)
    borders=np.array(borders)
    #c_array contains all the values which we will map on a colour (here: 3rd coloumn)
    #Values are rescaled to fit in [0,1]
    #reshaping is necessary because otherwise numpy does not fill the shape tuple which would give errors otherwise
    c_array=(np.array([p[ccol] for p in values]).reshape((-1,1))-minc)/(1.0*(maxc-minc))
    #m stands for matrix
    m=np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
    if maxextent_x<0:
        #find maximum extent
        mx=max([abs(p[xcol]) for p in values])
        #rescale position of points to maxextent_x and maxextent_y and stretch z-direction
        m[0][0]=-maxextent_x/(2.0*mx)
    else:
        m[0][0]=maxextent_x
    if maxextent_y<=0:
        #find maximum extent
        my=max([abs(p[ycol]) for p in values])
        #rescale position of points to maxextent_x and maxextent_y and stretch z-direction
        m[1][1]=-maxextent_y/(2.0*my)
    else:
        m[1][1]=maxextent_y
    m[2][2]=scale
    #v_array contains the points properly scaled and shifted
    #reshaping is necessary because otherwise numpy does not fill the shape tuple which would give errors otherwise
    v_array=np.array([[v[xcol],v[ycol],v[zcol]] for v in values]).reshape((-1,3)).dot(m)+np.array(shift)
    #b_array contains either nothing (i.e. Nones) or the values from the coloumn that shall also be passed back
    if not backcol==None:
        b_array=np.array([b[backcol] for b in values]).reshape((-1,1))
    else:
        b_array=np.array([None]*len(values)).reshape((-1,1))
    nr_borders=len(borders)
    #process skipping
    if not skip==None:
        for m,n,o in skip:
            startindex=o
            new_c_array=np.array([]).reshape(0,c_array.shape[1])
            new_b_array=np.array([]).reshape(0,b_array.shape[1])
            new_v_array=np.array([]).reshape(0,v_array.shape[1])
            if m>0 and n>0:
                while startindex<len(c_array):
                    endindex=startindex+m
                    new_c_array=np.concatenate((new_c_array,c_array[startindex:endindex]),axis=0)
                    new_b_array=np.concatenate((new_b_array,b_array[startindex:endindex]),axis=0)
                    new_v_array=np.concatenate((new_v_array,v_array[startindex:endindex]),axis=0)
                    startindex=endindex+n
            c_array=new_c_array
            b_array=new_b_array
            v_array=new_v_array
    #c: will contain colour code, now it contains the value that we want to assign a colour code
    #v: cotains x,y and z coordinate of point as it will be drawn on the screen
    #b: contains value from coloumn that was requested to be passed back, if so desired
    #without the [] around c and b, a list of length 1 would be returned
    for [c],v,[b] in zip(c_array,v_array,b_array):
        j=-1
        #make sure we're never out of the interval [0,1]
        if c<=0.0:
            c=0
            j=1
        elif c>=1.0:
            c=1.0
            j=nr_borders-1
        else:
            #find the range defined by the borders that we are in
            for i in xrange(nr_borders):
                if c<borders[i]:
                    j=i
                    break
        #perform linear interpolation
        c1=colours[j-1]
        c2=colours[j]
        b1=borders[j-1]
        b2=borders[j]
        fraction=1.0*(c-b1)/(b2-b1)
        c=fraction*c2+(1.0-fraction)*c1
        if b==None:
            yield (list(c),list(v))
        else:
            yield (list(c),list(v),b)


#opengl stuff is blatantly taken from the tutorial at
#http://pydoc.net/Python/PyOpenGL-Demo/3.0.0/PyOpenGL-Demo.NeHe.lesson5/
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
# Number of the glut window.
# A general OpenGL initialization function.  Sets all of the initial parameters. 
def InitGL(Width, Height):              # We call this right after our OpenGL window is created.
    """
    Initialize OpenGL
    """
    glClearColor(0.0, 0.0, 0.0, 0.0)    # This Will Clear The Background Color To Black
    glClearDepth(1.0)                   # Enables Clearing Of The Depth Buffer
    glDepthFunc(GL_LESS)                # The Type Of Depth Test To Do
    glEnable(GL_DEPTH_TEST)             # Enables Depth Testing
    glShadeModel(GL_SMOOTH)             # Enables Smooth Color Shading
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()                    # Reset The Projection Matrix
                                        # Calculate The Aspect Ratio Of The Window
    gluPerspective(60.0, float(Width)/float(Height), 0.01, 10000.0)
    glMatrixMode(GL_MODELVIEW)
    return

# This function is called to properly adjust the relative positions of plot and camers
def GLAdjustCamera(axes, angles, translation):
    """
    Properly adjust relative positions of plot and camera
    """
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) # Clear The Screen And The Depth Buffer
    glLoadIdentity()                   # Reset The View
    glTranslatef(*translation)     #move the camera to the correct position
    for an,ax in zip(angles,axes):  #rotate the view properly
        glRotatef(an,*ax)
    return

# This function is called to display the surface on the screen
def DrawGLTrimesh(faces, colourscale, globalscale=1, globalskip=0, elements_per_line=None, ccol=2, colours=[[0.0,0.0,0.0],[0.8,0.3,0.0],[1.0,1.0,0.0],[1.0,1.0,1.0]], borders=[0.0,0.2,0.7,1.0]):
    """
    This function tell OpenGL to draw a mesh.
    faces: should look like [A,B,C,...] where A,B and C are faces and have
           3 elements each p1,p2 and p3, which all are Cartesian vectors in
           3 dimensions
    colourscale: 2 element list or tuple with the z-value that should be
                 associated with the "lowest" colour and the "highest" colour
    globalscale: scale the whole plot in all dimensions by this much
    globalskip: every 1 triangle, do not draw the next this many triangles
                and every "line" of triangles, do not plot the next this many
                lines
    elements_per_line: assuming a quadratic mesh (in x and y dimensions), how
                       many points are there per line
                       ignored if not set
    colours: give a new colour scale for the plot
    borders: give the borders for the new colour scale
    """
    glBegin(GL_TRIANGLES)
    #get variables via a generator
    #c stands for colour and p stands for point, i.e. position
    if globalskip>0:
        if elements_per_line==None:
            skip=[(3,3*globalskip,0)]
        else:
            el=elements_per_line-1
            #per space between 2 data points (there are el such spaces per line), there are 6 points in the below array, 3 per triangle
            skip=[(6*el,6*el*globalskip,0),(3,3*globalskip,0)]
    else:
        skip=None

    for c,p in yield_values([point for triangle in faces for point in triangle],minc=colourscale[0],maxc=colourscale[1],scale=1.0*globalscale,skip=skip,maxextent_x=1.0*globalscale,maxextent_y=1.0*globalscale, ccol=ccol, colours=colours, borders=borders):
        glColor3f(*c)
        glVertex3f(*p)
    glEnd()
    return

# This function is called to display the spheres on the screen
def DrawGLSpheres(spheres, colourscale, globalscale=1, globalskip=0, elements_per_line=None, sphere_elements=50):
    """
    This function tell OpenGL to draw spheres.
    spheres: should look like [A,B,C,...] where A,B and C are spheres consisting of
             x,y,z,r where x,y,z are Cartesian coordinates and r is the radius
    colourscale: 2 element list or tuple with the z-value that should be
                 associated with the "lowest" colour and the "highest" colour
    globalscale: scale the whole plot in all dimensions by this much
    globalskip: every 1 sphere, do not draw the next this many spheres
                and every "line" of spheres, do not plot the next this many
                lines
    elements_per_line: assuming a quadratic mesh (in x and y dimensions), how
                       many spheres are there per line
                       ignored if not set
    sphere_elements: controls how many azimutal and longitudinal elements are drawn per
                     sphere
    """

    if globalskip>0:
        if elements_per_line==None:
            skip=[(1,globalskip,0)]
        else:
            #per space between 2 data points (there are el such spaces per line), there are 6 points in the below array, 3 per triangle
            skip=[(elements_per_line,elements_per_line*globalskip,0),(1,globalskip,0)]
        skip=[(1,globalskip,0)]
    else:
        skip=None
    #c stands for colour and p stands for point, i.e. position and r stands for radius 
    for c,p,r in yield_values(spheres,minc=colourscale[0],maxc=colourscale[1],scale=1.0*globalscale,colours=[[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]],borders=[0.0,1.0],backcol=3,skip=skip,maxextent_x=1.0*globalscale,maxextent_y=1.0*globalscale):
        p=np.array(p)
        glTranslatef(*p)
        glColor3f(*c)
        glutSolidSphere(globalscale*r, sphere_elements, sphere_elements)
        glTranslatef(*(-p))
    return

# This function glues together all the other displaying functions
def GLMainDisplay(axes, angles, translation, faces, face_colourscale, draw_faces, spheres, sphere_colourscale, draw_spheres, globalscale, elements_per_line_faces, elements_per_line_spheres, globalskip):

    #adjust camera
    GLAdjustCamera(axes, angles, translation)
    #draw what is desired
    if draw_faces:
        DrawGLTrimesh(faces, face_colourscale, globalscale=globalscale, globalskip=globalskip, elements_per_line=elements_per_line_faces)
    if draw_spheres:
        DrawGLSpheres(spheres, sphere_colourscale, globalscale=globalscale, globalskip=globalskip,  elements_per_line=elements_per_line_spheres)
    #  since this is double buffered, swap the buffers to display what just got drawn. 
    glutSwapBuffers()

import Image
import re
def snap(size,basename,format,count,extension):
    """
    Save a OpenGL screenshot to disk.
    size: (xres,yres): size of the image
    basename: path plus first part of filename
    count: if the images shall be numbered, say what number this image shall have
    format: to have fixed width numbering, declare the printf-type format string (like %6d)
    extension: filetype, png recommended
    """
    pixels=[]
    filename=re.sub('\s', '0', basename+format%(count)+"."+extension)
    screenshot = glReadPixels(0,0,size[0],size[1],GL_RGBA,GL_UNSIGNED_BYTE)
    snapshot = Image.frombuffer("RGBA",size,screenshot,"raw","RGBA",0,0)
    snapshot.save(filename)
    print filename
    return count+1


