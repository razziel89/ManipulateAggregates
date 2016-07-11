import sys

import openbabel as op
import pybel as p

from energyscan.scan import _prepare_molecules, _double_array
from manipulate_molecules import read_from_file



def similarityscreening_main(parser):
    gets   = parser.get_str
    geti   = parser.get_int
    getf   = parser.get_float
    getb   = parser.get_boolean
    do_calculate = not(getb("config_check"))
    
    #value for progress reports
    if geti("progress") not in [0,1,2]:
        raise ValueError('Wrong value for parameter "progress" given. Must be 0,1 or 2.')
    progress = geti("progress")

    mol1 = read_from_file(gets("geometry1"),ff=None)
    mol2 = read_from_file(gets("geometry2"),ff=None)

    obmol = _prepare_molecules(mol1,mol2,align=getb("prealign"))

    std_map = op.StdMapStringString()
    #add the appropriate configuration paramters to the std::map<std::string,std::string>
    std_map['rcutoff'] = str(getf("rmsd_min"))      #this way I can be sure it's actually a floating point number
    std_map['ecutoff'] = str(getf("energy_cutoff"))
    std_map['ffname']  =     gets("forcefield")
    
    #try to find the chosen force field
    if gets("forcefield").lower() not in ["uff", "mmff94", "gaff", "ghemical"]:
        raise ValueError('Wrong foce field given. Only "uff", "mmff94", "gaff" and "ghemical" will be accepted.')
    temp_ff = op.OBForceField.FindType(gets("forcefield").lower())
    if temp_ff is None:
        raise ValueError("Somehow there was an error loading the forcefield %s (although it should be known to OpenBabel)."%(gets("forcefield").lower()))
    del temp_ff

    if not do_calculate:
        return

    #to avoid segfaults, define some bogus input parameters that would normally be given via the command-line
    in_out_options = op.OBConversion()
    in_out_options.SetInFormat('nul')
    in_out_options.SetOutFormat('nul')

    #create a new OBAggregate instance that can contain a single aggregate and
    # that will walk through all the minima that were found. Each of these geometries
    # will be added to obmol as a new conformer so that the OBOp SimSearch can
    # perform its screening duty
    saveobmol = op.OBAggregate(obmol)   #copy constructor
    obmol.DeleteConformer(0)            #clean all conformer information

    tempmol = None
    emin    = float("inf")

    old_angles = (None,None,None)
    ang        = [0.0,0.0,0.0] #current angles
    disp       = [0.0,0.0,0.0] #current displacement
    if progress>0:
        print "Adding minima geometries to data structure ..."
    with open(gets("minima_file_load")) as f:
        for line in f:
            if not line.startswith("#"):
                linevals = line.rstrip().split()
                disp     = list(map(float,linevals[1:4]))
                pos_disp = _double_array(disp)
                neg_disp = _double_array([-v for v in disp])
                ang      = tuple(map(float,linevals[4:7]))
                energy   = float(linevals[7])
                if energy < emin:
                    emin = energy
                if ang != old_angles:
                    if progress>0:
                        print "... re-creating aggregate with new angles: (%.2f,%.2f,%.2f) ..."%ang
                    del tempmol
                    tempmol   = op.OBAggregate(saveobmol) #copy constructor
                    #since python needs quite some time to access an objects member, saving a member saves time
                    transfunc   = tempmol.TranslatePart
                    rotfunc     = tempmol.RotatePart
                    coordfunc   = tempmol.GetCoordinates
                    a1,a2,a3    = ang
                    old_angles  = ang
                    rotfunc(0,1,a1)
                    rotfunc(0,2,a2)
                    rotfunc(0,3,a3)
                transfunc(0,pos_disp)
                #actually deep-copy the new coordinates to avoid segfaults
                obmol.AddConformer(coordfunc(),True) 
                transfunc(0,neg_disp)

    if emin != float("inf"):
        #add global minimum energy to list of config parameters
        std_map['emin'] = str(emin)
  
    #align all aggregates with their centers to (0,0,0) to improve screening success
    obmol.Center()

    print "... %d aggregates have been processed ..."%(obmol.NumConformers())
    if obmol.NumConformers() <= 0:
        print "\n... not a single conformer was processed, hence we're done ...\n"
        return

    #force openbabel to be verbose if detailed progress reports were requested
    if progress==1:
        std_map['verbose'] = "true"

    simscreen = op.OBOp.FindType('simscreen')

    if getf("energy_cutoff")>0:

        print "\n... starting RMSD-pre-screening by energy ...\n"

        #First, only sort out those aggregates that do not pass the energy filter.
        #No checks for symmetry are being performed
        std_map['rcutoff'] = str(0.0)
        #perform the pre-screening
        simscreen.Do(obmol,'', std_map, in_out_options)

        if progress>0:
            print "... %d aggregates passed energy filter ...\n\n"%(obmol.NumConformers())

    else:
        print "\n... no screening by energy requested so none will be performed ...\n\n"

    step = 0
    #energy screening has already been performed if it was desired so do not do it again
    std_map['ecutoff'] = str(-100.0)
    #screen until fewer than nr_geometries agregates are left
    rmsd     = getf("rmsd_min")
    rmsdstep = getf("rmsd_step")
    maxagg   = geti("nr_geometries")
    aggfunc  = obmol.NumConformers
    while aggfunc() > maxagg:
        step += 1
        std_map['rcutoff'] = str(rmsd)
        simscreen.Do(obmol,'', std_map, in_out_options)
        if progress>0:
            print "... %d aggregates passed screening step %d at rmsd %f ...\n\n"%(aggfunc(),step,rmsd)
        rmsd += rmsdstep

    if progress>0:
        print "... %d aggregates passed screening ..."%(aggfunc())

    #write all conformers that passed the filter to file
    writefile = p.Outputfile("xyz",gets("screened_xyz"),overwrite=True)
    pybelmol  = p.Molecule(obmol)
    nr_conformers = obmol.NumConformers()
    commentfunc   = obmol.SetTitle
    setconffunc   = obmol.SetConformer
    for conf in xrange(nr_conformers):
        commentfunc("Conformer %d/%d"%(conf+1,nr_conformers))
        setconffunc(conf)
        writefile.write(pybelmol)
    writefile.close()
