import electric_potential as ep
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from collection.opengl import *
import sys
import re
#_c stands for collection
global gl_c        #contains all variables that need to be global for OpenGL to be able to display stuff
gl_c = {}          #initialize as empty dictionary

#default values for OpenGL related stuff
gl_c['fullscreen']         =   False               #whether fullscreen shall be activated right from the start
gl_c['resolution']         =   (1024,768)          #default resolution in pixels
gl_c['maxextent']          =   -1.0                #how to scale the view by default (default -1.0 means do not scale)
gl_c['translation']        =   [-0.5,0.0,-260.0]   #position of camera with respect to drawn structures
gl_c['angles']             =   [0.0,-10.0,0.0]     #rotational angles around axes defined below
gl_c['axes']               =   [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] #x,y and z axes (could be arbitrary)
gl_c['globalscale']        =   10                  #global scale for whole plot, i.e., zooming in or out
gl_c['faces']              =   []                  #will contain all the triangles that make up the surface and shall be drawn
gl_c['face_colourscale']   =   []                  #will contain all the minimum and maximum z value that is associated with
gl_c['spheres']            =   []                  #will contain all the spheres
gl_c['sphere_colours']     =   []                  #will contain colours for all spheres
gl_c['snap_count']         =   0                   #the counter for snapped images
gl_c['snap_title']         =   'snap'              #the title for the snapped images
                                                   
gl_c['colours']   =   []
gl_c['borders']   =   []
gl_c['graphical_output']   =   True                #whether graphical output is desired or not
gl_c['window']             =   1                   #the number of the window used for the main display

class WrongScalingTypeError(Exception):
    pass

class ParseError(Exception):
    pass

def _ReSizeGLScene(Width, Height):
    global gl_c
    if Height == 0:                     # Prevent A Divide By Zero If The Window Is Too Small 
        Height = 1
    glViewport(0, 0, Width, Height)     # Reset The Current Viewport And Perspective Transformation
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60.0, float(Width)/float(Height), 0.01, 10000.0)
    glMatrixMode(GL_MODELVIEW)
    gl_c['resolution']=(Width,Height)

#this function is repeatedly being called by OpenGL 
def _main_control():
    global gl_c
    #draw everything
    GLAdjustCamera(gl_c['axes'], gl_c['angles'], gl_c['translation'])
    if len(gl_c['faces'])>0:
        DrawGLTrimesh(gl_c['faces'], gl_c['face_colourscale'],globalscale=gl_c['globalscale'],ccol=3, colours=gl_c['colours'], borders=gl_c['borders'])
    if len(gl_c['spheres'])>0:
        DrawGLSpheres(gl_c['spheres'], (0,1), globalscale=gl_c['globalscale'], sphere_elements=50, colour_list=gl_c['sphere_colours'])
    glutSwapBuffers()

# The function called by the OpenGL main loop whenever a key is pressed
def _keyPressed(*args):
    global gl_c
    # If escape is pressed, kill everything.
    if args[0] == '\033': #this is the escape sequence for the ESC key
        glutLeaveMainLoop()
        gl_c['running']=False
    if args[0] == "+":
        gl_c['globalscale']=gl_c['globalscale']+0.1
    if args[0] == "=":
        gl_c['globalscale']=gl_c['globalscale']+0.1
    if args[0] == "-":
        gl_c['globalscale']=gl_c['globalscale']-0.1
    if args[0] == "w":
        gl_c['translation'][1]+=1.0
    if args[0] == "s":
        gl_c['translation'][1]-=1.0
    if args[0] == "a":
        gl_c['translation'][0]-=1.0
    if args[0] == "d":
        gl_c['translation'][0]+=1.0
    if args[0] == "q":
        gl_c['translation'][2]-=1.0
    if args[0] == "e":
        gl_c['translation'][2]+=1.0
    if args[0] == "i":
        gl_c['angles'][0]-=1.5
    if args[0] == "k":
        gl_c['angles'][0]+=1.5
    if args[0] == "j":
        gl_c['angles'][1]-=1.5
    if args[0] == "l":
        gl_c['angles'][1]+=1.5
    if args[0] == "u":
        gl_c['angles'][2]+=1.5
    if args[0] == "o":
        gl_c['angles'][2]-=1.5
    if args[0] == ".":
        snap(gl_c['resolution'],gl_c['snap_title']+"_","%3d",gl_c['snap_count'],"png")
        gl_c['snap_count']+=1

def TopLevelGlInitialization(gl_c,zoom,resolution,title="Molecule Visualization"):
    gl_c['running']=True
    gl_c['globalscale'] *= zoom
    gl_c['resolution'] = resolution
    if not title=="Molecule Visualization":
        gl_c['snap_title']=title
    glutInit()
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
    glutInitWindowSize(*gl_c['resolution'])
    glutInitWindowPosition(0, 0)
    gl_c['window'] = glutCreateWindow(title)
    #glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_CONTINUE_EXECUTION)
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS)
    InitGL(*gl_c['resolution'])
    glutDisplayFunc(_main_control)
    glutIdleFunc(_main_control)
    glutReshapeFunc(_ReSizeGLScene)
    glutKeyboardFunc(_keyPressed)

#def _subcommandRange(sc):
#    range=float(sc[0])
#    nsteps=int(sc[1])
#    return [[1.0*range/nsteps]]*nsteps

def _subcommandRangeNew(sc):
    s=sc[0].split("-")
    ranges=map(float,s)
    nsteps=int(sc[1])
    return [[1.0*r/nsteps for r in ranges]]*nsteps

def _parseTrajectoryNew(trajectory):
    actions=trajectory.split(",")
    commands=[a.split("|") for a in actions]
    parsed=[]
    for c in commands:
        p=[]
        if re.match("(r[123][\+-]|t[123][\+-]|z[\+-])+",c[0]):
            nr_commands=sum([c[0].count(s) for s in ["r","t","z"]])
            nr_ranges=c[1].count("-")+1
            if not nr_ranges==nr_commands:
                raise ParseError("Could not parse "+c[0]+" and "+c[1]+" for an interlaced command.")
        for match in re.findall("r[123][\+-]|t[123][\+-]|z[\+-]",c[0]):
            subp=[]
            if re.match("r",match):
                subp.append("angles")
            elif re.match("t",match):
                subp.append("translation")
            elif re.match("z",match):
                subp.append("globalscale")
            if re.match(".1",match):
                subp.append(0)
            elif re.match(".2",match):
                subp.append(1)
            elif re.match(".3",match):
                subp.append(2)
            else:
                subp.append(3)
            if re.match("(.|..)\+",match):
                subp.append("+")
            elif re.match("(.|..)-",match):
                subp.append("-")
            p.append(subp)
        for sc in _subcommandRangeNew(c[1].split("/")):
            parsed.append([p,sc])
    return parsed

#def _parseTrajectory(trajectory):
#    actions=trajectory.split(",")
#    commands=[a.split("|") for a in actions]
#    parsed=[]
#    for c in commands:
#        p=[]
#        if re.match("(r|t)[123][+-]",c[0]):
#            if re.match("r",c[0]):
#                p.append("angles")
#            elif re.match("t",c[0]):
#                p.append("translation")
#            else:
#                raise ParseError("Could not parse "+c[0]+" for a command.1")
#            if re.match(".1",c[0]):
#                p.append(0)
#            elif re.match(".2",c[0]):
#                p.append(1)
#            elif re.match(".3",c[0]):
#                p.append(2)
#            else:
#                raise ParseError("Could not parse "+c[0]+" for a command.2")
#            if re.match("..\+",c[0]):
#                p.append("+")
#            elif re.match("..-",c[0]):
#                p.append("-")
#            else:
#                raise ParseError("Could not parse "+c[0]+" for a command.3")
#        elif re.match("z[+-]",c[0]):
#            p.append("globalscale")
#            if re.match(".\+",c[0]):
#                p.append("+")
#            elif re.match(".-",c[0]):
#                p.append("-")
#            else:
#                raise ParseError("Could not parse "+c[0]+" for a command.4")
#        else:
#            raise ParseError("Could not parse "+c[0]+" for a command.5")
#        for sc in _subcommandRange(c[1].split("/")):
#            parsed.append(p+sc)
#    return parsed
    
def TopLevelRenderFunction(gl_c,rendertrajectory,title):
    if re.match(".*,n(,d)*$",rendertrajectory):
        snapping=False
    else:
        snapping=True
    if re.match(".*,d(,n)*$",rendertrajectory):
        drop=True
    else:
        drop=False
    while re.match(".*(,n|,d)+$",rendertrajectory):
        rendertrajectory=rendertrajectory.rstrip(",n")
        rendertrajectory=rendertrajectory.rstrip(",d")
    actions={"+":lambda a,b:a+b, "-":lambda a,b:a-b, "*":lambda a,b:a*b, "/":lambda a,b:a/b}
    parsed=_parseTrajectoryNew(rendertrajectory)
    digits=len(str(len(parsed)+1))
    snapformat="%+"+str(digits)+"d"
    if drop:
        while gl_c['running']:
            glutMainLoopEvent()
            _main_control()
    _main_control()
    if snapping:
        snap(gl_c['resolution'],gl_c['snap_title']+"_","%3d",gl_c['snap_count'],"png")
    gl_c['snap_count']+=1
    for ac in parsed:
        for com,val in zip(ac[0],ac[1]):
            if com[1]==3:
                gl_c[com[0]]=actions[com[2]](gl_c[com[0]],val)
            else:
                gl_c[com[0]][com[1]]=actions[com[2]](gl_c[com[0]][com[1]],val)
        _main_control()
        if snapping:
            snap(gl_c['resolution'],gl_c['snap_title']+"_","%3d",gl_c['snap_count'],"png")
            gl_c['snap_count']+=1

def PlotGL_Surface(mol,zoom,nr_refinements=1,title="Molecule Visualization",resolution=(1024,768),scale='independent',high_contrast=False,rendertrajectory=None):
    global gl_c
    faces=[]
    corners, potential = mol.get_vdw_surface_potential(vertex='corners', triangulation=faces,nr_refinements=nr_refinements)
    #get actual coordinates for indices
    coordinates=mol.get_coordinates()
    charges=mol.get_partial_charges()
    if scale == 'independent':
        gl_c['face_colourscale']=(min(potential),max(potential))
    elif scale == 'dependent':
        abs_overall=abs(max([abs(min(potential)),abs(max(potential))]))
        gl_c['face_colourscale']=(-abs_overall,abs_overall)
    else:
        raise WrongScalingTypeError("Scale must be either independent or dependent")
    gl_c['faces']=[[[f[0],f[1],f[2],ep.potential_at_points([f], charges, coordinates)[0]] for f in face] for face in faces]
    if high_contrast:
        middlecolour=[0.0,0.0,0.0]
        sidecolours=([0.4,0.4,1.0],[1.0,0.4,0.4])
    else:
        middlecolour=[0.2,0.2,0.2]
        sidecolours=([0.0,0.0,1.0],[1.0,0.0,0.0])
    #gl_c['colours']   =   [[0.0,0.0,1.0],[0.2,0.2,0.2],[1.0,0.0,0.0]]
    gl_c['colours']   =   [sidecolours[0],middlecolour,sidecolours[1]]
    gl_c['borders']   =   [0.0,-min(potential)/(max(potential)-min(potential)),1.0]

    TopLevelGlInitialization(gl_c,zoom,resolution,title=title)
    if rendertrajectory==None:
        glutMainLoop()
    else:
        TopLevelRenderFunction(gl_c,rendertrajectory,title)

def PlotGL_Spheres(mol,zoom,title="Molecule Visualization",resolution=(1024,768),spherescale=1,rendertrajectory=None):
    global gl_c

    #get actual coordinates for indices
    coordinates=mol.get_coordinates()
    vdw_radii=mol.get_vdw_radii()
    gl_c['sphere_colours']=mol.get_colours()
    gl_c['spheres']=[[c[0],c[1],c[2],r*spherescale] for c,r in zip(coordinates,vdw_radii)]

    TopLevelGlInitialization(gl_c,zoom,resolution,title=title)
    if rendertrajectory==None:
        glutMainLoop()
    else:
        TopLevelRenderFunction(gl_c,rendertrajectory,title)
