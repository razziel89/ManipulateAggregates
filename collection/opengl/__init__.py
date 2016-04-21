"""
function definitions for drawing stuff on the screen using OpenGL
This is just a handy collection of wrapper functions for opengl
"""

import numpy as np
NPROTX = np.eye(4,dtype=float)
NPROTY = np.eye(4,dtype=float)
NPROTZ = np.eye(4,dtype=float)
ROTMAT = np.eye(4,dtype=float)

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
def InitGL(Width, Height, use_light=False):              # We call this right after our OpenGL window is created.
    """
    Initialize OpenGL
    """
    glClearColor(1.0, 1.0, 1.0, 0.0)    # This Will Clear The Background Color To White
    glClearDepth(1.0)                   # Enables Clearing Of The Depth Buffer
    glDepthFunc(GL_LESS)                # The Type Of Depth Test To Do
    glEnable(GL_DEPTH_TEST)             # Enables Depth Testing
    glShadeModel(GL_SMOOTH)             # Enables Smooth Color Shading
    if use_light:
        specular = [ 0.1, 0.1, 0.1, 0.1 ];
        shininess = [ 30.0 ];
        position = [ 1.0, 0.0, 5.0, 0.0 ];
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
        glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0, GL_POSITION, position);
        glEnable(GL_COLOR_MATERIAL)
        glEnable(GL_LIGHTING);
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()                    # Reset The Projection Matrix
                                        # Calculate The Aspect Ratio Of The Window
    gluPerspective(60.0, float(Width)/float(Height), 0.01, 10000.0)
    glMatrixMode(GL_MODELVIEW)
    return

def GetCameraMatrix(angles):
    """
    Get the proper transformation matrix for the camera.
    """
    #rotate the view properly
    x,y,z = angles
    x *= -0.017453292519943295
    y *= -0.017453292519943295
    z *= -0.017453292519943295
    NPROTX[1,1] = np.cos(x)
    NPROTX[1,2] = -np.sin(x)
    NPROTX[2,1] = np.sin(x)
    NPROTX[2,2] = np.cos(x)
    NPROTY[0,0] = np.cos(y)
    NPROTY[0,2] = np.sin(y)
    NPROTY[2,0] = -np.sin(y)
    NPROTY[2,2] = np.cos(y)
    NPROTZ[0,0] = np.cos(z)
    NPROTZ[0,1] = -np.sin(z)
    NPROTZ[1,0] = np.sin(z)
    NPROTZ[1,1] = np.cos(z)
    return np.dot(np.dot(NPROTX,NPROTY),NPROTZ)

# This function is called to properly adjust the relative positions of plot and camers
def GLAdjustCamera(angles, translation):
    """
    Properly adjust relative positions of plot and camera
    """
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) # Clear The Screen And The Depth Buffer
    glLoadIdentity()                   # Reset The View
    glTranslatef(*translation)     #move the camera to the correct position
    glMultMatrixd(np.ndarray.flatten(GetCameraMatrix(angles)))
    return

def GLGetViewMatrix():
    mat = [0.0]*16
    return [i for i in glGetDoublev(GL_MODELVIEW_MATRIX,mat)]

# This function is called to display the surface on the screen
def DrawGLTrimesh(faces, colourscale, globalscale=1, globalskip=0, elements_per_line=None, ccol=2, colours=[[0.0,0.0,0.0],[0.8,0.3,0.0],[1.0,1.0,0.0],[1.0,1.0,1.0]], borders=[0.0,0.2,0.7,1.0]):
    """
    This function tells OpenGL to draw a mesh.
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

# This function is called to translate the OpenGL data to povray format and write it to a file handle
def WritePovrayTrimesh(handle, indices, points, normals, colorvalues, colourscale, globalscale=1, globalskip=0, elements_per_line=None, ccol=2, colours=[[0.0,0.0,0.0],[0.8,0.3,0.0],[1.0,1.0,0.0],[1.0,1.0,1.0]], borders=[0.0,0.2,0.7,1.0]):
    """
    This function writes a trimesh in PovRay format to a file handle
    handle: file descriptor (or anything with a write method, really)
    indices: should look like [A,B,C,...] where A,B and C are faces and contain
             3 indices each, i.e., i1,i2 and i3
    points:  list of [float,float,float]
             The vertices describing the actual trimesh.
    normals: list of [float,float,float]
             One normal vector per vertex.
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
    transparency=0.0
    if len(points) != len(normals) or len(points) != len(colorvalues):
        raise ValueError("Length of 'points', 'normals' and 'colorvalues' are not identical.")
    points_colors=np.concatenate((points,colorvalues),axis=1)
    tab="    "
    tabcount=0
    handle.write(tab*tabcount+"mesh2 {\n")
    tabcount+=1
    handle.write(tab*tabcount+"vertex_vectors {\n")
    tabcount+=1
    handle.write(tab*tabcount+"%d,\n"%(len(points)))
    for c,p in yield_values(points_colors,
            minc=colourscale[0],maxc=colourscale[1],
            scale=1.0*globalscale,maxextent_x=1.0*globalscale,maxextent_y=1.0*globalscale,
            ccol=ccol, colours=colours, borders=borders):
        handle.write(tab*tabcount+"<%.4f,%.4f,%.4f>,\n"%tuple(p))
    tabcount-=1
    handle.write(tab*tabcount+"}\n")
    handle.write(tab*tabcount+"normal_vectors {\n")
    tabcount+=1
    handle.write(tab*tabcount+"%d,\n"%(len(normals)))
    for n in normals:
        handle.write(tab*tabcount+"<%.4f,%.4f,%.4f>,\n"%tuple(n))
    tabcount-=1
    handle.write(tab*tabcount+"}\n")
    handle.write(tab*tabcount+"texture_list {\n")
    tabcount+=1
    handle.write(tab*tabcount+"%d,\n"%(len(colorvalues)))
    for c,p in yield_values(points_colors,
            minc=colourscale[0],maxc=colourscale[1],
            scale=1.0*globalscale,maxextent_x=1.0*globalscale,maxextent_y=1.0*globalscale,
            ccol=ccol, colours=colours, borders=borders):
        handle.write(tab*tabcount+"RGBTVERT(<%.3f,%.3f,%.3f"%tuple(c))
        handle.write(",%.3f>),\n"%(transparency))
    tabcount-=1
    handle.write(tab*tabcount+"}\n")
    handle.write(tab*tabcount+"face_indices {\n")
    tabcount+=1
    handle.write(tab*tabcount+"%d\n"%(len(indices)))
    for i in indices:
        handle.write(tab*tabcount+"<%d,%d,%d>,"%tuple(i))
        handle.write("%d,%d,%d\n"%tuple(i))
    tabcount-=1
    handle.write(tab*tabcount+"}\n")
    handle.write(tab*tabcount+"inside_vector <0, 0, 1>\n")
    handle.write(tab*tabcount+"no_shadow\n")
    handle.write(tab*tabcount+"matrix <\n")
    tabcount+=1
    handle.write(tab*tabcount+"0.500000, 0.000000, 0.000000,\n")
    handle.write(tab*tabcount+"0.000000, 0.500000, 0.000000,\n")
    handle.write(tab*tabcount+"0.000000, 0.000000, 0.500000,\n")
    handle.write(tab*tabcount+"0.000000, 0.000000, 0.000000\n")
    tabcount-=1
    handle.write(tab*tabcount+">\n")
    tabcount-=1
    handle.write(tab*tabcount+"}\n")
    return

# This function is called to display the spheres on the screen
def DrawGLSpheres(spheres, colourscale, globalscale=1, globalskip=0, elements_per_line=None, sphere_elements=50, colour_list=None):
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
    use_gen=not(colour_list==None)
    if use_gen:
        gen=(c for c in colour_list)
    #c stands for colour and p stands for point, i.e. position and r stands for radius 
    for c,p,r in yield_values(spheres,minc=colourscale[0],maxc=colourscale[1],scale=1.0*globalscale,colours=[[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]],borders=[0.0,1.0],backcol=3,skip=skip,maxextent_x=1.0*globalscale,maxextent_y=1.0*globalscale):
        p=np.array(p)
        glTranslatef(*p)
        if use_gen:
            glColor3f(*gen.next())
        else:
            glColor3f(*c)
        glutSolidSphere(globalscale*r, sphere_elements, sphere_elements)
        glTranslatef(*(-p))
    return

# This function glues together all the other displaying functions
def GLMainDisplay(angles, translation, faces, face_colourscale, draw_faces, spheres, sphere_colourscale, draw_spheres, globalscale, elements_per_line_faces, elements_per_line_spheres, globalskip):

    #adjust camera
    GLAdjustCamera(angles, translation)
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
    filename=re.sub('\s', '0', basename+format%(count)+"."+extension)
    screenshot = glReadPixels(0,0,size[0],size[1],GL_RGBA,GL_UNSIGNED_BYTE)
    snapshot = Image.frombuffer("RGBA",size,screenshot,"raw","RGBA",0,0)
    snapshot.save(filename)
    print filename
    return count+1


def povray(size,basename,format,count,
        povray_data,colourscale,globalscale=1,
        colours=[[0.0,0.0,0.0],[0.8,0.3,0.0],[1.0,1.0,0.0],[1.0,1.0,1.0]],borders=[0.0,0.2,0.7,1.0]):
    """
    """
    extension="pov"
    filename=re.sub('\s', '0', basename+"%dx%d_"%(size[0],size[1])+format%(count)+"."+extension)
    handle = open(filename,"wb")
    handle.write("#macro RGBTVERT ( C1 )\n")
    handle.write("  texture { pigment { rgbt C1 }}\n")
    handle.write("#end\n")
    handle.write("#if (version < 3.5)\n")
    handle.write("#error \"This input is only compatible with PovRay 3.5 and above.\"\n")
    handle.write("#end\n")
    viewmat = [i for i in GLGetViewMatrix()]  #get the model view matrix from OpenGL
    camerapos = np.array(viewmat[12:15])
    #write out camera information and take it from OpenGL
    handle.write("camera {\n")
    handle.write("perspective\n")
    handle.write("direction <0,0,-1>\n")#set to OpenGL's default value
    handle.write("angle 45.000000\n")   #PovRay uses horizontal FOV and OpenGL vertical, so that needs to be adjuster
    handle.write("transform {\n")       #transform the camera to OpenGL's position
    handle.write("matrix <\n")
    handle.write("%.4f,%.4f,%.4f,\n"%tuple(viewmat[0:3]))  #first rotation
    handle.write("%.4f,%.4f,%.4f,\n"%tuple(viewmat[4:7]))  #second rotation
    handle.write("%.4f,%.4f,%.4f,\n"%tuple(viewmat[8:11])) #third rotation
    handle.write("%.4f,%.4f,%.4f\n"%tuple(camerapos))      #position
    handle.write(">\n")
    handle.write("inverse }\n") #needed because PovRay sets the camera position but OpenGL sets the models position
    handle.write("}\n")
    handle.write("light_source {\n")
    handle.write("  <%.4f,%.4f,%.4f>\n"%tuple(-camerapos-np.array([-50.0000, 50.0000, -50.0000])))
    handle.write("  color rgb<1.000, 1.000, 1.000>\n")
    handle.write("  parallel\n")
    handle.write("  point_at <0.0, 0.0, 0.0>\n")
    handle.write("}\n")
    handle.write("light_source {\n")
    handle.write("  <%.4f,%.4f,%.4f>\n"%tuple(-camerapos-np.array([50.0000, 100.0000, -50.0000])))
    handle.write("  color rgb<1.000, 1.000, 1.000>\n")
    handle.write("  parallel\n")
    handle.write("  point_at <0.0, 0.0, 0.0>\n")
    handle.write("}\n")
    handle.write("background {\n")
    handle.write("  color rgb<1.000, 1.000, 1.000>\n")
    handle.write("}\n")
    handle.write("#default { texture {\n")
    handle.write(" finish { ambient 0.000 diffuse 0.650 phong 0.1 phong_size 95.499 specular 0.220 }\n")
    handle.write("} }\n")
    WritePovrayTrimesh(handle, povray_data[0], povray_data[1], povray_data[2], povray_data[3],
            colourscale, globalscale=globalscale, ccol=3, colours=colours, borders=borders)
    handle.close()
    from subprocess import Popen
    f = Popen(["povray","+W%d"%(size[0]),"+H%d"%(size[1]),"-I%s"%(filename),"-O%s"%(filename+".png"),"+UA","+D","+X","+A","+FN"])
    f.wait()
    print filename+".png"
    return count+1
