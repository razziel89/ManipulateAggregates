import sys

import openbabel as op
import pybel as p

from energyscan.scan import _prepare_molecules, _double_array, CURSOR_UP_ONE, ERASE_LINE
from manipulate_molecules import read_from_file

#conversion factors from meV to the given unit
E_UNIT_CONVERSION = {
        "kJ/mol"        : 0.09648500,
        "kcal/mol"      : 0.02306035,
        }

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
    #check whether partitioning over nodes was switched on
    if not gets("partition") == "1/1":
        raise ValueError("Parallelization unequal 1/1 not suported for similarity screening.")

    mol1 = read_from_file(gets("geometry1"),ff=None)
    mol2 = read_from_file(gets("geometry2"),ff=None)

    obmol = _prepare_molecules(mol1,mol2,align=getb("prealign"))

    std_map = op.StdMapStringString()
    #add the appropriate configuration paramters to the std::map<std::string,std::string>
    std_map['rcutoff'] = str(getf("rmsd_min"))      #this way I can be sure it's actually a floating point number
    std_map['ffname']  =     gets("forcefield")
    
    #try to find the chosen force field
    if gets("forcefield").lower() not in ["uff", "mmff94", "gaff", "ghemical"]:
        raise ValueError('Wrong force field given. Only "uff", "mmff94", "gaff" and "ghemical" will be accepted.')
    temp_ff = op.OBForceField.FindType(gets("forcefield").lower())
    if temp_ff is None:
        raise RuntimeError("Somehow there was an error loading the forcefield %s (although it should be known to OpenBabel)."%(gets("forcefield").lower()))
    try:
        if getb("use_ff_units"):
            print "Using given energy in force field units: %.6f %s (equals %.6f meV)"%(
                    getf("energy_cutoff"),temp_ff.GetUnit(),getf("energy_cutoff")/E_UNIT_CONVERSION[temp_ff.GetUnit()])
            std_map['ecutoff'] = str(getf("energy_cutoff"))
        else:
            std_map['ecutoff'] = str(getf("energy_cutoff")*E_UNIT_CONVERSION[temp_ff.GetUnit()])
            print "Converting energy cutoff to force field units: %s meV -> %.6f %s"%(
                    gets("energy_cutoff"),getf("energy_cutoff")*E_UNIT_CONVERSION[temp_ff.GetUnit()],temp_ff.GetUnit())
    except KeyError as e:
        raise RuntimeError("Unknown unit type '%s' of the chosen force field '%s', cannot convert the energy cutoff in meV to that unit. KeyError was: %s. Known units are: %s"%(
            temp_ff.GetUnit(),gets("forcefield").lower(),e,", ".join([t for t in E_UNIT_CONVERSION])))
    finally:
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
    sameff  = -1

    is_sorted  = True
    old_angles = (-float("inf"),-float("inf"),-float("inf"))
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
                if sameff!=0 and energy < emin:
                    emin = energy
                if ang != old_angles:
                    if progress>0:
                        print ERASE_LINE+"... re-creating aggregate with new angles: (%8.2f,%8.2f,%8.2f) ..."%ang+CURSOR_UP_ONE
                    del tempmol
                    tempmol   = op.OBAggregate(saveobmol) #copy constructor
                    #since python needs quite some time to access an objects member, saving a member saves time
                    transfunc   = tempmol.TranslatePart
                    rotfunc     = tempmol.RotatePart
                    coordfunc   = tempmol.GetCoordinates
                    a1,a2,a3    = ang
                    if is_sorted and ( a1<old_angles[0] or a2<old_angles[1] or a3<old_angles[2] ):
                        #angles should be in a monotonically nondecreasing order
                        is_sorted = False
                    old_angles  = ang
                    rotfunc(0,1,a1)
                    rotfunc(0,2,a2)
                    rotfunc(0,3,a3)
                transfunc(0,pos_disp)
                #actually deep-copy the new coordinates to avoid segfaults
                obmol.AddConformer(coordfunc(),True) 
                transfunc(0,neg_disp)
            else:
                l = line.split()
                if len(l) >= 3 and l[1] == "FF:":
                    print "\n... determining force field used to create the minima file %s ..."%(
                            gets("minima_file_load"))
                    if l[2].lower() != gets("forcefield").lower():
                        print "... old force field '%s' is not the same as the current one '%s' ..."%(
                                l[2].lower(),gets("forcefield").lower())
                        sameff = 0
                        emin = float("inf")
                    else:
                        print "... minima file was created using the current force field ..."
                        sameff = 1
        if progress>0:
            print

    if sameff == -1:
        print >>sys.stderr,"WARNING: could not determine force field used to create the minima file %s from the comments line."%(
                gets("minima_file_load"))
        print >>sys.stderr,"         Assuming it's not the current one."
        emin = float("inf")

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
    if progress>0:
        print "... writing %d aggregates to file %s ..."%(aggfunc(),gets("screened_xyz"))
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

    if not is_sorted:
        print "WARNING: minima file was not in sorted order with respect to the angles."
        print "         Beware that results might change slightly if the order is changed."
