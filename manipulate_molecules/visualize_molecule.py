try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
    from collection.opengl import *
    gl_imported=True
except ImportError:
    gl_imported=False
import sys
import re
import os
#_c stands for collection
global gl_c        #contains all variables that need to be global for OpenGL to be able to display stuff
gl_c = {}          #initialize as empty dictionary

#default values for OpenGL related stuff
gl_c['fullscreen']         =   False               #whether fullscreen shall be activated right from the start
gl_c['resolution']         =   (1024,768)          #default resolution in pixels
gl_c['maxextent']          =   -1.0                #how to scale the view by default (default -1.0 means do not scale)
gl_c['translation']        =   [0.0,0.0,-260.0]   #position of camera with respect to drawn structures
gl_c['look_at']            =   [0.0,0.0,0.0]   
gl_c['angles']             =   [0.0,-10.0,0.0]     #rotational angles around x,y and z axes
gl_c['globalscale']        =   10                  #global scale for whole plot, i.e., zooming in or out
gl_c['faces']              =   []                  #will contain all the triangles that make up the surface and shall be drawn
gl_c['face_colourscale']   =   []                  #will contain all the minimum and maximum z value that is associated with
gl_c['spheres']            =   []                  #will contain all the spheres
gl_c['sphere_colours']     =   []                  #will contain colours for all spheres
gl_c['snap_count']         =   0                   #the counter for snapped images
gl_c['snap_title']         =   'snap'              #the title for the snapped images
gl_c['keys']               =   {}                  #will contain all the keys pressed
                                                   
gl_c['colours']   =   []
gl_c['borders']   =   []
gl_c['high_contrast']      =   False
gl_c['povray']             =   False
gl_c['povray_data']        =   []
gl_c['graphical_output']   =   True                #whether graphical output is desired or not
gl_c['window']             =   1                   #the number of the window used for the main display
gl_c['savefile']           =   None                #the name of the file where visualization shall be saved to
gl_c['savecount']          =   0                   #the number of visualization states already saved
gl_c['additional']         =   []                  #additional data. Only used to save the visualization

class WrongScalingTypeError(Exception):
    pass

class ParseError(Exception):
    pass

class ArbitraryInputError(Exception):
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
    GLAdjustCamera(gl_c['angles'], gl_c['translation'])
    if len(gl_c['faces'])>0:
        DrawGLTrimesh(gl_c['faces'], gl_c['face_colourscale'],globalscale=gl_c['globalscale'],ccol=3, colours=gl_c['colours'], borders=gl_c['borders'])
    if len(gl_c['spheres'])>0:
        DrawGLSpheres(gl_c['spheres'], (0,1), globalscale=gl_c['globalscale'], sphere_elements=50, colour_list=gl_c['sphere_colours'])
    glutSwapBuffers()

def _keyPressed(*args):
    global gl_c
    keys=gl_c['keys']
    if args[0] == '\033': #this is the escape sequence for the ESC key
        keys["quit"]=True
    if args[0] == "+":
        keys["zoom+"]=True
    if args[0] == "=":
        keys["zoom+"]=True
    if args[0] == "-":
        keys["zoom-"]=True
    if args[0] == "w":
        keys["up"]=True
    if args[0] == "s":
        keys["down"]=True
    if args[0] == "a":
        keys["left"]=True
    if args[0] == "d":
        keys["right"]=True
    if args[0] == "q":
        keys["front"]=True
    if args[0] == "e":
        keys["back"]=True
    if args[0] == "i":
        keys["rot1+"]=True
    if args[0] == "k":
        keys["rot1-"]=True
    if args[0] == "j":
        keys["rot2+"]=True
    if args[0] == "l":
        keys["rot2-"]=True
    if args[0] == "u":
        keys["rot3+"]=True
    if args[0] == "o":
        keys["rot3-"]=True
    if args[0] == ".":
        keys["snap"]=True
    if args[0] == "p":
        if gl_c["povray"]:
            keys["povray"]=True
        else:
            print >>sys.stderr,"WARNING: PovRay support has either not been activated or is not supported by this type of visualization."
    if args[0] == ",":
        if gl_c['savefile'] is not None:
            gl_c['savecount']+=1
            SaveVisualizationState(gl_c,gl_c['savefile'],prefix=str(gl_c['savecount']-1)+"_")
    _evaluateKeyPressed()

def _evaluateKeyPressed():
    global gl_c
    keys=gl_c['keys']
    # If escape is pressed, kill everything.
    if keys["quit"]: #this is the escape sequence for the ESC key
        glutLeaveMainLoop()
        gl_c['running']=False
    if keys["zoom+"]:
        gl_c['globalscale']=gl_c['globalscale']+0.1
    if keys["zoom-"]:
        gl_c['globalscale']=gl_c['globalscale']-0.1
    if keys["up"]:
        gl_c['translation'][1]+=1.0
    if keys["down"]:
        gl_c['translation'][1]-=1.0
    if keys["left"]:
        gl_c['translation'][0]-=1.0
    if keys["right"]:
        gl_c['translation'][0]+=1.0
    if keys["front"]:
        gl_c['translation'][2]-=1.0
    if keys["back"]:
        gl_c['translation'][2]+=1.0
    if keys["rot1+"]:
        gl_c['angles'][0]-=1.5
    if keys["rot1-"]:
        gl_c['angles'][0]+=1.5
    if keys["rot2+"]:
        gl_c['angles'][1]-=1.5
    if keys["rot2-"]:
        gl_c['angles'][1]+=1.5
    if keys["rot3+"]:
        gl_c['angles'][2]+=1.5
    if keys["rot3-"]:
        gl_c['angles'][2]-=1.5
    if keys["snap"]:
        snap(gl_c['resolution'],gl_c['snap_title']+"_","%3d",gl_c['snap_count'],"png")
        gl_c['snap_count']+=1
    if keys["povray"]:
        povray(gl_c['resolution'],gl_c['snap_title']+"_","%3d",gl_c['snap_count'],
                gl_c['povray_data'],gl_c['face_colourscale'],globalscale=gl_c['globalscale'],
                colours=gl_c['colours'], borders=gl_c['borders'])
        gl_c['snap_count']+=1

def _keyReleased(*args):
    global gl_c
    keys=gl_c['keys']
    if args[0] == "+":
        keys["zoom+"]=False
    if args[0] == "=":
        keys["zoom+"]=False
    if args[0] == "-":
        keys["zoom-"]=False
    if args[0] == "w":
        keys["up"]=False
    if args[0] == "s":
        keys["down"]=False
    if args[0] == "a":
        keys["left"]=False
    if args[0] == "d":
        keys["right"]=False
    if args[0] == "q":
        keys["front"]=False
    if args[0] == "e":
        keys["back"]=False
    if args[0] == "i":
        keys["rot1+"]=False
    if args[0] == "k":
        keys["rot1-"]=False
    if args[0] == "j":
        keys["rot2+"]=False
    if args[0] == "l":
        keys["rot2-"]=False
    if args[0] == "u":
        keys["rot3+"]=False
    if args[0] == "o":
        keys["rot3-"]=False
    if args[0] == ".":
        keys["snap"]=False
    if args[0] == "p":
        keys["povray"]=False

def _initializeKeys(keys):
    keys["quit"]=False
    keys["zoom+"]=False
    keys["zoom-"]=False
    keys["zoom-"]=False
    keys["up"]=False
    keys["down"]=False
    keys["left"]=False
    keys["right"]=False
    keys["front"]=False
    keys["back"]=False
    keys["rot1+"]=False
    keys["rot1-"]=False
    keys["rot2+"]=False
    keys["rot2-"]=False
    keys["rot3+"]=False
    keys["rot3-"]=False
    keys["snap"]=False
    keys["povray"]=False

## The function called by the OpenGL main loop whenever a key is pressed
#def _keyPressed(*args):
#    global gl_c
#    # If escape is pressed, kill everything.
#    if args[0] == '\033': #this is the escape sequence for the ESC key
#        glutLeaveMainLoop()
#        gl_c['running']=False
#    if args[0] == "+":
#        gl_c['globalscale']=gl_c['globalscale']+0.1
#    if args[0] == "=":
#        gl_c['globalscale']=gl_c['globalscale']+0.1
#    if args[0] == "-":
#        gl_c['globalscale']=gl_c['globalscale']-0.1
#    if args[0] == "w":
#        gl_c['translation'][1]+=1.0
#    if args[0] == "s":
#        gl_c['translation'][1]-=1.0
#    if args[0] == "a":
#        gl_c['translation'][0]-=1.0
#    if args[0] == "d":
#        gl_c['translation'][0]+=1.0
#    if args[0] == "q":
#        gl_c['translation'][2]-=1.0
#    if args[0] == "e":
#        gl_c['translation'][2]+=1.0
#    if args[0] == "i":
#        gl_c['angles'][0]-=1.5
#    if args[0] == "k":
#        gl_c['angles'][0]+=1.5
#    if args[0] == "j":
#        gl_c['angles'][1]-=1.5
#    if args[0] == "l":
#        gl_c['angles'][1]+=1.5
#    if args[0] == "u":
#        gl_c['angles'][2]+=1.5
#    if args[0] == "o":
#        gl_c['angles'][2]-=1.5
#    if args[0] == ".":
#        snap(gl_c['resolution'],gl_c['snap_title']+"_","%3d",gl_c['snap_count'],"png")
#        gl_c['snap_count']+=1

def TopLevelGlInitialization(gl_c,zoom,resolution,title="Molecule Visualization",use_light=False):
    _initializeKeys(gl_c['keys'])
    gl_c['running']=True
    gl_c['globalscale'] *= zoom
    gl_c['resolution'] = resolution
    if not title=="Molecule Visualization":
        gl_c['snap_title']=title
    if not gl_imported:
        return False
    try:
        os.environ['DISPLAY']
    except KeyError:
        return False
    if bool(glutInit):
        glutInit()
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
        #glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST) #this should enable antialiasing
        glutInitWindowSize(*gl_c['resolution'])
        glutInitWindowPosition(0, 0)
        gl_c['window'] = glutCreateWindow(title)
        #glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_CONTINUE_EXECUTION)
        glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS)
        InitGL(*gl_c['resolution'],use_light=use_light)
        glutDisplayFunc(_main_control)
        glutIdleFunc(_main_control)
        glutReshapeFunc(_ReSizeGLScene)
        glutKeyboardUpFunc(_keyReleased)
        glutKeyboardFunc(_keyPressed)
        return True
    else:
        return False

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

def _set_high_contrast():
        middlecolour=[0.0,0.0,0.0]
        #sidecolours=([0.4,0.4,1.0],[1.0,0.4,0.4])
        sidecolours=([0.3,0.6,1.0],[1.0,0.4,0.4])
        return [sidecolours[0],middlecolour,sidecolours[1]]

def _set_low_contrast():
        middlecolour=[0.2,0.2,0.2]
        sidecolours=([0.0,0.2,1.0],[1.0,0.0,0.0])
        return [sidecolours[0],middlecolour,sidecolours[1]]

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

def SaveVisualizationState(obj,filename,prefix=""):
    from cPickle import dump as save
    f=open(prefix+filename,'w')
    save(obj,f,-1)
    f.close()
    print 'Saved visualization state to file '+prefix+filename

def LoadVisualization(filename):
    from cPickle import load as load
    f=open(filename,'r')
    obj=load(f)
    f.close()
    return obj
    
def TopLevelRenderFunction(gl_c,rendertrajectory):
    if re.match(".*,n(,d|,s)*$",rendertrajectory):
        snapping=False
    else:
        snapping=True
    if re.match(".*,d(,n|,s)*$",rendertrajectory):
        drop=True
    else:
        drop=False
    save=False
    if re.match(".*,s(,n|,d)*$",rendertrajectory):
        if gl_c['savefile'] is not None:
            save=True
    while re.match(".*(,n|,d|,s)+$",rendertrajectory):
        rendertrajectory=rendertrajectory.rstrip(",n")
        rendertrajectory=rendertrajectory.rstrip(",d")
        rendertrajectory=rendertrajectory.rstrip(",s")
    actions={"+":lambda a,b:a+b, "-":lambda a,b:a-b, "*":lambda a,b:a*b, "/":lambda a,b:a/b}
    parsed=_parseTrajectoryNew(rendertrajectory)
    digits=len(str(len(parsed)+1))
    snapformat="%+"+str(digits)+"d"
    if drop:
        while gl_c['running']:
            glutMainLoopEvent()
            _main_control()
    _main_control()
    _main_control()
    if snapping:
        gl_c['snap_count']+=1
        snap(gl_c['resolution'],gl_c['snap_title']+"_","%3d",gl_c['snap_count']-1,"png")
    if save: 
        gl_c['savecount']+=1
        SaveVisualizationState(gl_c,gl_c['savefile'],prefix=str(gl_c['savecount']-1)+"_")
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
        if save: 
            gl_c['savecount']+=1
            SaveVisualizationState(gl_c,gl_c['savefile'],prefix=str(gl_c['savecount']-1)+"_")

def RenderExtern(ext_gl_c,resolution=(1024,768),rendertrajectory=None,title="Molecule Visualization",savefile=None,high_contrast=None,invert_potential=False):
    global gl_c
    for key in ext_gl_c:
        try:
            gl_c[key] = ext_gl_c[key]
        except KeyError:
            print >>sys.stderr, 'Key '+key+' not found when loading, skipping key.'
    if high_contrast is not None:
        if high_contrast:
            gl_c['colours'] = _set_high_contrast()
            gl_c['high_contrast'] = True
        else:
            gl_c['colours'] = _set_low_contrast()
            gl_c['high_contrast'] = False
    check=TopLevelGlInitialization(gl_c,1,resolution,title=title)
    if not check:
        print >>sys.stderr, "Cannot initialize OpenGL, will save visialization state."
        if savefile is not None:
            gl_c['savefile']=savefile['file']
        else:
            gl_c['savefile']="GLError"
        SaveVisualizationState(gl_c,"GLErrorDump_"+gl_c['savefile'])
        return
    if savefile is not None:
        gl_c['savefile']=savefile['file']
        if savefile['start']:
            SaveVisualizationState(gl_c,"start_"+gl_c['savefile'])
    if rendertrajectory==None:
        glutMainLoop()
    else:
        TopLevelRenderFunction(gl_c,rendertrajectory)
    if savefile is not None:
        if savefile['end']:
            SaveVisualizationState(gl_c,"end_"+gl_c['savefile'])

def PlotGL_Surface(mol,zoom,nr_refinements=1,title="Molecule Visualization",resolution=(1024,768),scale='independent',high_contrast=False,rendertrajectory=None,charges=None,ext_potential=None,manip_func=None,invert_potential=False,config=None,savefile=None,povray=False):

    if ext_potential is not None and charges is not None:
        raise ArbitraryInputError("Cannot use external charges and external potential at the same time.")

    global gl_c
    import numpy as np

    temp = np.array(mol.get_vdw_surface(get='faces',nr_refinements=nr_refinements,povray=povray))
    if povray:
        faces = np.array(temp[0])
        povray_indices  = np.array(temp[1])
        povray_vertices = np.array(temp[2])
        povray_normals  = np.array(temp[3])
    else:
        faces = np.array(temp)

    if charges is None:
        #get actual coordinates for indices
        charges=mol.get_partial_charges()
        coordinates=mol.get_coordinates()
    else:
        coordinates,charges=charges
        charges=-np.array(charges)

    try:
        if os.environ['PROGRESS']=='0':
            prog_report=False
        else:
            prog_report=True
    except KeyError:
        prog_report=True

    faces.shape=(-1,3)
    if ext_potential is not None:
        try:
            from FireDeamon import InterpolationPy as interpol
        except ImportError as e:
            raise ImportError("Error importing FireDeamon.InterpolationPy which is needed to use an external potential.",e)
        potential=np.array(interpol(ext_potential[0],ext_potential[1],faces,prog_report=prog_report,config=config))
        if povray:
            povray_potential=np.array(interpol(ext_potential[0],ext_potential[1],povray_vertices,prog_report=prog_report,config=config))
    else:
        try:
            from FireDeamon import ElectrostaticPotentialPy as potential_at_points
        except ImportError as e:
            raise ImportError("Error importing FireDeamon.ElectrostaticPotentialPy which is needed to visiualize an empirical potential.",e)
        potential=np.array(potential_at_points(faces, charges, coordinates, prog_report=prog_report))
        if povray:
            povray_potential=np.array(potential_at_points(povray_vertices, charges, coordinates, prog_report=prog_report))
    potential.shape=(-1,3,1)
    faces.shape=(-1,3,3)
    if povray:
        povray_potential.shape=(-1,1)

    if invert_potential:
        potential*=-1
        if povray:
            povray_potential*=-1
    
    if manip_func is not None:
        faces.shape=(-1,3)
        faces=np.array([ manip_func(f) for f in faces ])
        faces.shape=(-1,3,3)
        if povray:
            povray_vertices=np.array([ manip_func(f) for f in povray_vertices ])

    gl_c['faces']=list(np.concatenate((faces,potential),axis=2))
    if povray:
        gl_c['povray_data']=[povray_indices,povray_vertices,povray_normals,povray_potential]
        np.concatenate((povray_vertices,povray_potential),axis=1)

    if scale == 'independent' and ( abs(np.min(potential))<=0.0 or abs(np.max(potential))<=0.0 ):
        print >>sys.stderr, "WARNING: independent colour scaling won't work, will switch to dependent colour scaling."
        scale='dependent'

    if scale == 'independent':
        gl_c['face_colourscale']=(np.min(potential),np.max(potential))
    elif scale == 'dependent':
        abs_overall=abs(max([abs(np.min(potential)),abs(np.max(potential))]))
        gl_c['face_colourscale']=(-abs_overall,abs_overall)
    else:
        raise WrongScalingTypeError("Scale must be either independent or dependent")
    

    gl_c['povray'] = povray
    if high_contrast:
        gl_c['colours'] = _set_high_contrast()
        gl_c['high_contrast'] = True
    else:
        gl_c['colours'] = _set_low_contrast()
        gl_c['high_contrast'] = False
    #gl_c['colours']   =   [[0.0,0.0,1.0],[0.2,0.2,0.2],[1.0,0.0,0.0]]
    gl_c['borders']   =   [0.0,-np.min(potential)/(np.max(potential)-np.min(potential)),1.0]

    check=TopLevelGlInitialization(gl_c,zoom,resolution,title=title)
    if not check:
        print >>sys.stderr, "Cannot initialize OpenGL, will save visialization state."
        if savefile is not None:
            gl_c['savefile']=savefile['file']
        else:
            gl_c['savefile']="GLError"
        SaveVisualizationState(gl_c,"GLErrorDump_"+gl_c['savefile'])
        return
    if savefile is not None:
        gl_c['savefile']=savefile['file']
        if savefile['start']:
            SaveVisualizationState(gl_c,"start_"+gl_c['savefile'])
    if rendertrajectory==None:
        glutMainLoop()
    else:
        TopLevelRenderFunction(gl_c,rendertrajectory)
    if savefile is not None:
        if savefile['end']:
            SaveVisualizationState(gl_c,"end_"+gl_c['savefile'])

def PlotGL_Spheres(mol,zoom,title="Molecule Visualization",resolution=(1024,768),spherescale=1,rendertrajectory=None):
    global gl_c

    #get actual coordinates for indices
    coordinates=mol.get_coordinates()
    vdw_radii=mol.get_vdw_radii()
    gl_c['sphere_colours']=mol.get_colours()
    gl_c['spheres']=[[c[0],c[1],c[2],r*spherescale] for c,r in zip(coordinates,vdw_radii)]

    check=TopLevelGlInitialization(gl_c,zoom,resolution,title=title,use_light=True)
    if not check:
        print >>sys.stderr, "Cannot initialize OpenGL, aborting."
        return
    if rendertrajectory==None:
        glutMainLoop()
    else:
        TopLevelRenderFunction(gl_c,rendertrajectory)
