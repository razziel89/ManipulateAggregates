"""Screen aggregate geometries with respect to similarity.

This subsubmodule is part of ManipulateAggregates.energyscan. It implements the
third step of the 3-step procedure that creates low energy aggregate geometries.

Parallelization is only supported for the determination of pointgroups.

@package ManipulateAggregates.energyscan.similarityscreening
"""

#This file is part of ManipulateAggregates.
#
#Copyright (C) 2016 by Torsten Sachse
#
#ManipulateAggregates is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ManipulateAggregates is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with ManipulateAggregates.  If not, see <http://www.gnu.org/licenses/>.
import sys
import os
import re
from multiprocessing import Pool, Event

import logging
logger = logging.getLogger(__name__)
try:
    import openbabel as openbabel
except ImportError:
    logger.warning("Could not import openbabel")
try:
    import pybel
except ImportError:
    logger.warning("Could not import pybel")

try:
    from ..energyscan.scan import _prepare_molecules, _double_array, CURSOR_UP_ONE, ERASE_LINE, _double_array
except ImportError:
    logger.warning("Could not import  _prepare_molecules, _double_array, CURSOR_UP_ONE, ERASE_LINE or _double_array from ..energyscan.scan")
try:
    from ..aggregate import read_from_file, FILETYPEDICT
except ImportError:
    logger.warning("Could not import read_from_file or FILETYPEDICT from ..aggregate")

## conversion factors from meV to the declared force field units
E_UNIT_CONVERSION = {
        "kJ/mol"        : 0.09648500,
        "kcal/mol"      : 0.02306035,
        }

## associate the name of each pointgroup with its subgroups
# C1 is implicit
SUBGROUPS = {
        'C4v' :  ('C2v', 'C2', 'C4', 'Cs'),
        'D4h' :  ('C4v', 'C2v', 'D2h', 'C2h', 'C4h', 'D2d', 'C2', 'C4', 'S4', 'Cs', 'Ci', 'D4', 'D2'),
        'C2v' :  ('C2', 'Cs'),
        'D6h' :  ('C2v', 'D3h', 'C6v', 'S6', 'D2h', 'C6h', 'C2h', 'C3', 'C2', 'C6', 'C3h', 'Cs', 'D3d', 'C3v', 'Ci', 'D6', 'D2', 'D3'),
        'C6v' :  ('C2v', 'C3', 'C2', 'C6', 'Cs', 'C3v'),
        'D6d' :  ('C2v', 'C6v', 'D2d', 'C3', 'C2', 'C6', 'S4', 'Cs', 'C3v', 'D6', 'D2', 'D3'),
        'S6' :  ('C3', 'Ci'),
        'D4d' :  ('C4v', 'C2v', 'C2', 'C4', 'Cs', 'S8', 'D4', 'D2'),
        'D2h' :  ('C2v', 'C2h', 'C2', 'Cs', 'Ci', 'D2'),
        'S8' :  ('C2', 'C4'),
        'C6h' :  ('S6', 'C2h', 'C3', 'C2', 'C6', 'C3h', 'Cs', 'Ci'),
        'C2h' :  ('C2', 'Cs', 'Ci'),
        'C4h' :  ('C2h', 'C2', 'C4', 'S4', 'Cs', 'Ci'),
        'D2d' :  ('C2v', 'C2', 'S4', 'Cs', 'D2'),
        'D3d' :  ('C2v', 'S6', 'D2h', 'C2h', 'C3', 'C2', 'Cs', 'C3v', 'Ci', 'D2', 'D3'),
        'C3v' :  ('C3', 'Cs'),
        'C8v' :  ('C4v', 'C2v', 'C8', 'C2', 'C4', 'Cs'),
        'C8' :  ('C2', 'C4'),
        'O' :  ('C3', 'C2', 'C4', 'T', 'D4', 'D2', 'D3'),
        'D8h' :  ('C4v', 'C8h', 'C2v', 'D4d', 'D2h', 'D4h', 'C2h', 'C4h', 'D2d', 'C8v', 'C8', 'C2', 'C4', 'S4', 'Cs', 'Ci', 'S8', 'D8', 'D4', 'D2'),
        'C8h' :  ('C2h', 'C4h', 'C8', 'C2', 'C4', 'S4', 'Cs', 'Ci', 'S8'),
        'C2' :  (),
        'D8d' :  ('C4v', 'C2v', 'C8v', 'C8', 'C2', 'C4', 'Cs', 'D8', 'D4', 'D2'),
        'C7' :  (),
        'C6' :  ('C3', 'C2'),
        'C5' :  (),
        'C4' :  ('C2',),
        'D5h' :  (),
        'C5v' :  ( 'C5', 'Cs'),
        'C3h' :  ('C3', 'Cs'),
        'Oh' :  ('C4v', 'C2v', 'S6', 'D2h', 'D4h', 'C2h', 'C4h', 'D2d', 'C3', 'C2', 'C4', 'Td', 'S4', 'Cs', 'T', 'D3d', 'C3v', 'O', 'Th', 'Ci', 'D4', 'D2', 'D3'),
        'I' :  ('C3', 'C2', 'C5', 'T', 'D5', 'D2', 'D3'),
        'K' :  (),
        'Kh' :  ( 'K', 'Cs', 'Cinfv', 'Ci'),
        'D5d' :  ('C2v', 'D2h', 'C2h', 'C2', 'C5', 'C5v', 'Cs', 'Ci', 'D5', 'D2'),
        'Td' :  ('C2v', 'D2d', 'C3', 'C2', 'S4', 'Cs', 'T', 'C3v', 'D2', 'D3'),
        'C7v' :  ( 'C7', 'Cs'),
        'C3' :  (),
        'T' :  ('C3', 'C2', 'D2', 'D3'),
        'D3h' :  ('C2v', 'C3', 'C2', 'C3h', 'Cs', 'C3v', 'D2', 'D3'),
        'Cinfv' :  ( 'K', 'Cs'),
        'D7h' :  ('C2v', 'C2', 'C7', 'C7v', 'Cs', 'C7h', 'D7', 'D2'),
        'Dinfh' :  ('C2v', 'C2h', 'C2', 'K', 'Kh', 'Cs', 'Cinfv', 'Ci'),
        'S4' :  ('C2',),
        'C5h' :  ( 'C5', 'Cs'),
        'D7d' :  ('C2v', 'D2h', 'C2h', 'C2', 'C7', 'C7v', 'Cs', 'Ci', 'D7', 'D2'),
        'Th' :  ('C2v', 'S6', 'D2h', 'C2h', 'C3', 'C2', 'Cs', 'T', 'D3d', 'C3v', 'Ci', 'D2', 'D3'),
        'Ci' :  (),
        'C7h' :  ( 'C7', 'Cs'),
        'Ih' :  ('C2v', 'S6', 'D2h', 'C2h', 'C3', 'C2', 'C5', 'C5v', 'I', 'D5d', 'Cs', 'T', 'D3d', 'C3v', 'Th', 'Ci', 'D5', 'D2', 'D3'),
        'D6' :  ('C3', 'C2', 'C6', 'D2', 'D3'),
        'D8' :  (),
        'Cs' :  (),
        'D7' :  ('C2', 'C7', 'D2'),
        'D4' :  ('C2', 'C4', 'D2'),
        'D5' :  ('C2', 'C5', 'D2'),
        'D2' :  ('C2',),
        'D3' :  ('C3', 'C2', 'D2'),
        'C1' :  (),
        }

global data_ss

def _pointgroup_parallel_init(initobmol, terminating):
    """Allow easy data sharing between processes without pickling.

    Args:
        initobmol: (OBMol) contains the conformers whose pointgroups are to
            be determined
        terminating: (return value of multiprocessing.Event()) specifies 
            whether or not a parallel computation has been terminated
            prematurely or not
    """
    global data_ss
    data_ss = (initobmol, terminating)

def _get_pg_thread(args):
    global data_ss

    threadobmol, terminating = data_ss

    try:
        if not terminating.is_set():
            i,tolerance = args
            sym = openbabel.OBPointGroup()
            sym.Setup(threadobmol,i)
            pg = sym.IdentifyPointGroup(tolerance)
            del sym
    except KeyboardInterrupt:
        print >>sys.stderr, "Terminating worker process "+str(os.getpid())+" prematurely."
    return i,pg

def _get_pg(obmol, defaultobmol, subgroups, c1, filename, progress, postalign, screenregex):
    #if this function is executed, at least one of "screening by pointgroup" or "saving all conformers
    #belonging to a pointgroup" is to be performed.
    if progress>0:
        print "...determining pointgroups..."
        if subgroups:
            print "...high-symmetry conformers are also included in subgroups..."
    do_write = (filename is not None)
    if do_write:
        extension = filename.split(".")[-1]
        if FILETYPEDICT.has_key(extension) and not extension==filename:
            filename = filename[0:-(len(extension)+1)]
        getname = lambda pg: filename+"_"+pg+".xyz"
        tot_nr_mols = obmol.NumConformers()
    do_screen = (screenregex is not None)
    #in this special case, all C1 geometries have to be screened
    #despite the fact that no regex was given
    if not c1 and not do_screen:
        #an empty regex matches everything
        screenregex = re.compile("")
        do_screen = True

    if not do_screen and not do_write:
        raise RuntimeError("Unhandled internal error - the programme should never get here.")

    try:
        nr_threads = int(os.environ["OMP_NUM_THREADS"])
    except KeyError:
        nr_threads = 1
    except ValueError:
        nr_threads = 1

    terminating = Event()

    pool = Pool(nr_threads, initializer=_pointgroup_parallel_init, initargs=(obmol, terminating))   #NODEBUG
    #global data_ss                  #DEBUG
    #data_ss = (obmol, terminating)  #DEBUG

    tolerance = 0.01    #hard-coded so far
    args = [[i,tolerance] for i in xrange(0,obmol.NumConformers())]

    if do_write:
        writemols = {}
    if do_screen:
        screenmol = openbabel.OBAggregate(defaultobmol)
        screenmol.DeleteConformers(0,screenmol.NumConformers()-1)
        if screenmol.NumConformers()!=0:
            raise RuntimeError("Could not clear conformer information, still %d left."%(screenmol.NumConformers()))
        screenpgs = {}

    count = 0
    if progress>0:
        print "...analysed pointgroup of conformer: %d/%d..."%(count,tot_nr_mols)+CURSOR_UP_ONE
    try:
        for temp in pool.imap(_get_pg_thread, args):    #NODEBUG
        #for arg in args:                               #DEBUG
        #    temp = _get_pg_thread(arg)                 #DEBUG
            i,thread_pg = temp
            count += 1
            if progress>0:
                print ERASE_LINE+"...analysed pointgroup of conformer: %d/%d..."%(count,tot_nr_mols)+CURSOR_UP_ONE
            addgroups = (thread_pg,)
            if thread_pg.lower() != 'c1':
                addgroups += ('C1',)
            if subgroups:
                addgroups += SUBGROUPS[thread_pg]
            if do_write:
                for pg in addgroups:
                    if (pg.lower() != "c1") or (pg.lower() == "c1" and c1):
                        if writemols.has_key(pg):
                            writemols[pg].AddConformer(obmol.GetConformer(i),True)
                        else:
                            tempmol = openbabel.OBAggregate(defaultobmol)
                            tempmol.DeleteConformers(0,tempmol.NumConformers()-1) #clear all conformer information
                            if tempmol.NumConformers()!=0:
                                raise RuntimeError("Could not clear conformer information, still %d left."%(tempmol.NumConformers()))
                            tempmol.AddConformer(obmol.GetConformer(i),True)
                            writemols[pg] = tempmol
            if do_screen:
                for pg in addgroups:
                    if (pg.lower() != "c1") or (pg.lower() == "c1" and c1):
                        #add if at least one of the subgroups matches the regex
                        if re.match(screenregex,pg) is not None:
                            screenmol.AddConformer(obmol.GetConformer(i),True)
                            screenpgs[thread_pg] = screenpgs.get(thread_pg,0)+1
                            break
        pool.close()    #NODEBUG
        pool.join()     #NODEBUG
    except KeyboardInterrupt as e:
        print >>sys.stderr,"Caught keyboard interrupt."
        pool.terminate()    #NODEBUG
        pool.join()         #NODEBUG
        print >>sys.stderr,"Terminating main routine prematurely."
        raise e

    if progress>0:
        print

    if do_write:
        for pg,writeobmol in writemols.iteritems():
            pgfilename = getname(pg)
            if progress>0:
                print "...writing %4d aggregates of pointgroup %s to file %s..."%(writeobmol.NumConformers(),pg,pgfilename)
            writefile = pybel.Outputfile("xyz",pgfilename,overwrite=True)
            pybelmol  = pybel.Molecule(writeobmol)
            nr_conformers = writeobmol.NumConformers()
            commentfunc   = writeobmol.SetTitle
            setconffunc   = writeobmol.SetConformer
            if postalign:
                alignfunc     = writeobmol.Align
                aligncenter   = _double_array([0.0,0.0,0.0])
                alignaxis1    = _double_array([1.0,0.0,0.0])
                alignaxis2    = _double_array([0.0,1.0,0.0])
            for conf in xrange(nr_conformers):
                commentfunc("Conformer %d/%d"%(conf+1,nr_conformers))
                setconffunc(conf)
                if postalign:
                    alignfunc(aligncenter, alignaxis1, alignaxis2)
                writefile.write(pybelmol)
            writefile.close()
        del writemols
        if progress>0:
            print

    if do_screen:
        if progress>0:
            print "...reporting pointgroups of conformers that were not screened..."
            if subgroups:
                print "...the highest pointgroups are reported even if only a subgroup matches the regex..."
            for pg,count in screenpgs.iteritems():
                print "...keeping %4d aggregates of pointgroup %s..."%(count,pg)
                if subgroups:
                    tempstring = ", ".join(SUBGROUPS[thread_pg])
                    if len(tempstring)!=0:
                        print "...they also belong to subgroups: %s..."%(tempstring)
        del obmol
        if progress>0:
            print
        return screenmol
    else:
        if progress>0:
            print
        return obmol

def similarityscreening_main(parser):
    """Main control function for the similarity screening procedure.

    Args:
        parser: (of class ManipulateAggregates.collection.read.SectionlessConfigParser)
            contains information about the config file. Defines the methods
            "get_str", "get_int", "get_float" and "get_boolean" to get the
            appropriate data type.
    """
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
    nr_ats1 = mol1.obmol.NumAtoms()
    nr_ats2 = mol2.obmol.NumAtoms()

    obmol = _prepare_molecules(mol1,mol2,align=getb("prealign"))

    std_map = openbabel.StdMapStringString()
    #add the appropriate configuration paramters to the std::map<std::string,std::string>
    std_map['ffname']  =     gets("forcefield")

    #try to find the chosen force field
    if gets("forcefield").lower() not in ["uff", "mmff94", "gaff", "ghemical"]:
        raise ValueError('Wrong force field given. Only "uff", "mmff94", "gaff" and "ghemical" will be accepted.')
    temp_ff = openbabel.OBForceField.FindType(gets("forcefield").lower())
    if temp_ff is None:
        raise RuntimeError("Somehow there was an error loading the forcefield %s (although it should be known to OpenBabel)."%(gets("forcefield").lower()))
    try:
        if getf("energy_cutoff")>=0:
            if getb("use_ff_units"):
                print "...using given energy in force field units: %.6f %s (equals %.6f meV)"%(
                        getf("energy_cutoff"),temp_ff.GetUnit(),getf("energy_cutoff")/E_UNIT_CONVERSION[temp_ff.GetUnit()])
                std_map['ecutoff'] = str(getf("energy_cutoff"))
            else:
                std_map['ecutoff'] = str(getf("energy_cutoff")*E_UNIT_CONVERSION[temp_ff.GetUnit()])
                print "...converting given energy cutoff to force field units: %s meV -> %.6f %s"%(
                        gets("energy_cutoff"),getf("energy_cutoff")*E_UNIT_CONVERSION[temp_ff.GetUnit()],temp_ff.GetUnit())
        else:
            std_map['ecutoff'] = str(-100)
    except KeyError as e:
        raise RuntimeError("Unknown unit type '%s' of the chosen force field '%s', cannot convert the energy cutoff in meV to that unit. KeyError was: %s. Known units are: %s"%(
            temp_ff.GetUnit(),gets("forcefield").lower(),e,", ".join([t for t in E_UNIT_CONVERSION])))
    finally:
        del temp_ff

    postalign = getb("postalign")
    geti('symprec')
    geti("maxscreensteps")

    if not(gets("consider_h1")=="" and (gets("consider_h2") in ("","SAME"))):
        try:
            tmp_h1 = list(map(int,gets("consider_h1").split(",")))
        except ValueError as e:
            raise ValueError("Could not parse consider_h1, must be comma-separated ints.")
        if len(tmp_h1)!=0:
            if min(tmp_h1)<1 or max(tmp_h1)>nr_ats1:
                raise ValueError("Indices for consider_h1 must be >=%d and <=%d"%(1,nr_ats1))
        if gets("consider_h2")=="SAME":
            if gets("geometry1")!=gets("geometry2") and len(tmp_h1)!=0:
                raise ValueError("Can only use 'SAME' for consider_h2 when geometry1 and geometry2 are identical.")
            else:
                tmp_h2 = [i+nr_ats1 for i in tmp_h1]
        else:
            try:
                tmp_h2 = list(map(lambda s: int(s)+nr_ats1,gets("consider_h2").split(",")))
            except ValueError as e:
                raise ValueError("Could not parse 'consider_h2', must be comma-separated ints.")
        if len(tmp_h2)!=0:
            if min(tmp_h2)<nr_ats1+1 or max(tmp_h2)>nr_ats1+nr_ats2:
                raise ValueError("Indices for consider_h1 must be >=%d and <=%d"%(1,nr_ats2))
        important_hs = ",".join(map(str,tmp_h1 + tmp_h2))
        std_map["--imp-H"] = important_hs

    pgstep = -1
    if getb("pointgroups"):
        if gets("pgstep") == "last":
            pgstep = gets("pgstep")
        elif gets("pgstep") == "first":
            pgstep = 1
        else:
            pgstep = geti("pgstep")
            if pgstep<0:
                raise ValueError("The given pgstep must be >=0.")
    getb("subgroups")
    getb("exclude_c1")
    if getb("pgwrite"):
        pgfilename = gets("pgfile")
    else:
        pgfilename = None
    if len(gets("pgregex")) != 0:
        pgregex = re.compile(gets("pgregex"))
    else:
        pgregex = None
    if pgstep!=-1 and (len(gets("pgregex")) == 0 and not getb("pgwrite")):
        pgstep = -1
        print >>sys.stderr,"WARNING: pgwrite is False and no pgregex given -> pointgroups will not be determined in step "+gets("pgstep")

    if not do_calculate:
        return

    #to avoid segfaults, define some bogus input parameters that would normally be given via the command-line
    in_out_options = openbabel.OBConversion()
    in_out_options.SetInFormat('nul')
    in_out_options.SetOutFormat('nul')

    #create a new OBAggregate instance that can contain a single aggregate and
    # that will walk through all the minima that were found. Each of these geometries
    # will be added to obmol as a new conformer so that the OBOp SimSearch can
    # perform its screening duty
    saveobmol = openbabel.OBAggregate(obmol)   #copy constructor
    obmol.DeleteConformers(0,obmol.NumConformers()-1) #clean all conformer information
    if obmol.NumConformers()!=0:
        raise RuntimeError("Could not clear conformer information, still %d left."%(obmol.NumConformers()))

    tempmol = openbabel.OBAggregate(saveobmol)
    sameff  = True

    with open(gets("minima_file_load")) as f:
        #angles should be in a monotonically nondecreasing order
        angles = [tuple(map(float,line.split()[4:7])) for line in f if not line.startswith("#")]
        if not angles == list(sorted(angles)):
            print >>sys.stderr,"WARNING: minima file was not in sorted order with respect to the angles."
            print >>sys.stderr,"         Beware that results might change slightly if the order is changed."
        del angles

    old_angles = (-float("inf"),-float("inf"),-float("inf"))
    ang        = [0.0,0.0,0.0] #current angles
    disp       = [0.0,0.0,0.0] #current displacement
    if progress>0:
        print "...adding minima geometries to data structure..."
    printcount=0
    with open(gets("minima_file_load")) as f:
        transfunc   = tempmol.TranslatePart
        rotfunc     = tempmol.RotatePart
        coordfunc   = tempmol.GetCoordinates
        for line in f:
            if not line.startswith("#"):
                linevals = line.rstrip().split()
                disp     = list(map(float,linevals[1:4]))
                pos_disp = _double_array(disp)
                neg_disp = _double_array([-v for v in disp])
                ang      = tuple(map(float,linevals[4:7]))
                if ang != old_angles:
                    if progress>0 and printcount%10==0:
                        print ERASE_LINE+"...re-creating aggregate with new angles: (%8.2f,%8.2f,%8.2f)..."%ang+CURSOR_UP_ONE
                        printcount=0
                    printcount += 1
                    tempmol.Assign(saveobmol)
                    #since python needs quite some time to access an objects member, saving a member saves time
                    a1,a2,a3    = ang
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
                    print "\n...determining force field used to create the minima file %s..."%(
                            gets("minima_file_load"))
                    if l[2].lower() != gets("forcefield").lower():
                        print "...old force field '%s' is not the same as the current one '%s'..."%(
                                l[2].lower(),gets("forcefield").lower())
                        sameff = False
                    else:
                        print "...minima file was created using the current force field..."
                        sameff = True
    if progress>0:
        print
  
    print "...%d aggregates have been processed..."%(obmol.NumConformers())
    if obmol.NumConformers() <= 0:
        print "\n...not a single conformer was processed, hence we're done...\n"
        return

    #force openbabel to be verbose if detailed progress reports were requested
    if progress==1:
        std_map['verbose'] = "true"

    simscreen = openbabel.OBOp.FindType('simscreen')

    prescreen = False
    screenstring = ""
    if geti("symprec")>=0:
        if prescreen:
            screenstring += "and "
        prescreen = True
        screenstring += "symmetry "
        std_map['prec'] = str(geti('symprec'))
        #align all aggregates with their centers to (0,0,0)
        #and their third and second main axes to the x axis and the y axis, respectively,
        #to improve symmetry screening success
        std_map['ssalign'] = 'b'
    if getf("energy_cutoff")>0:
        if prescreen:
            screenstring += "and "
        screenstring += "energy "
        prescreen = True
    else:
        std_map.erase('ecutoff')

    step = 0
    if pgstep == step:
        nr_aggs = obmol.NumConformers()
        obmol = _get_pg(obmol, saveobmol, getb("subgroups"), not(getb("exclude_c1")), pgfilename, progress, postalign, pgregex)
        if obmol.NumConformers()>nr_aggs:
            raise RuntimeError("Number of conformers increased (%d -> %d) during pointgroup screening."%(nr_aggs,obmol.NumConformers()))
    if prescreen:
        nr_aggs = obmol.NumConformers()
        step += 1
        #First, only sort out those aggregates that do not pass the energy and symmetry filter.
        if progress>0:
            print "\n...starting "+screenstring+"pre-screening...\n"
        #perform the pre-screening
        success = simscreen.Do(obmol,'', std_map, in_out_options)
        if obmol.NumConformers()>nr_aggs:
            raise RuntimeError("Number of conformers increased (%d -> %d) during symmetry screening."%(nr_aggs,obmol.NumConformers()))
        if not success:
            raise RuntimeError("Error executing the SimScreen OBOp in OpenBabel.")
        if progress>0:
            print "...%d aggregates passed %sfilter...\n\n"%(obmol.NumConformers(),screenstring)
        #energy and symmetry screening have already been performed if they were desired so do not do that again
        std_map.erase('ecutoff')
        std_map.erase('ssalign')
        std_map.erase('prec')
    else:
        print "\n...skipping energy and symmetry pre-screening...\n"
    if prescreen and pgstep==step:
        nr_aggs = obmol.NumConformers()
        obmol = _get_pg(obmol, saveobmol, getb("subgroups"), not(getb("exclude_c1")), pgfilename, progress, postalign, pgregex)
        if obmol.NumConformers()>nr_aggs:
            raise RuntimeError("Number of conformers increased (%d -> %d) during pointgroup screening."%(nr_aggs,obmol.NumConformers()))

    success = True
    maxstep = geti("maxscreensteps")
    #screen until fewer than nr_geometries agregates are left
    rmsd     = getf("rmsd_min")
    rmsdstep = getf("rmsd_step")
    maxagg   = geti("nr_geometries")
    while success and obmol.NumConformers() > maxagg and step < maxstep:
        step += 1
        nr_aggs = obmol.NumConformers()
        std_map['rcutoff'] = str(rmsd)
        success = simscreen.Do(obmol,'', std_map, in_out_options)
        if obmol.NumConformers()>nr_aggs and success:
            raise RuntimeError("Number of conformers increased (%d -> %d) during screening step %d."%(nr_aggs,obmol.NumConformers(),step))
        if progress>0:
            print "...%d aggregates passed screening step %d at rmsd %f...\n\n"%(obmol.NumConformers(),step,rmsd)
        if pgstep == "last":
            nr_aggs = obmol.NumConformers()
            obmol = _get_pg(obmol, saveobmol, getb("subgroups"), not(getb("exclude_c1")), pgfilename, progress, postalign, pgregex)
            if obmol.NumConformers()>nr_aggs:
                raise RuntimeError("Number of conformers increased (%d -> %d) during pointgroup screening."%(nr_aggs,obmol.NumConformers()))
        rmsd += rmsdstep
        if pgstep == step:
            nr_aggs = obmol.NumConformers()
            obmol = _get_pg(obmol, saveobmol, getb("subgroups"), not(getb("exclude_c1")), pgfilename, progress, postalign, pgregex)
            if obmol.NumConformers()>nr_aggs:
                raise RuntimeError("Number of conformers increased (%d -> %d) during pointgroup screening."%(nr_aggs,obmol.NumConformers()))

    if not success:
        raise RuntimeError("Error executing the SimScreen OBOp in OpenBabel.")
    if step >= maxstep:
        print >>sys.stderr,"WARNING: maximum number of similarity screening steps exceeded"
    if success:
        if progress>0:
            print "...%d aggregates passed screening..."%(obmol.NumConformers())

        #write all conformers that passed the filter to file
        if progress>0:
            print "...writing %d aggregates to file %s..."%(obmol.NumConformers(),gets("screened_xyz"))
        writefile = pybel.Outputfile("xyz",gets("screened_xyz"),overwrite=True)
        pybelmol  = pybel.Molecule(obmol)
        nr_conformers = obmol.NumConformers()
        commentfunc   = obmol.SetTitle
        setconffunc   = obmol.SetConformer
        if postalign:
            alignfunc     = obmol.Align
            aligncenter   = _double_array([0.0,0.0,0.0])
            alignaxis1    = _double_array([1.0,0.0,0.0])
            alignaxis2    = _double_array([0.0,1.0,0.0])
        for conf in xrange(nr_conformers):
            commentfunc("Conformer %d/%d"%(conf+1,nr_conformers))
            setconffunc(conf)
            if postalign:
                alignfunc(aligncenter, alignaxis1, alignaxis2)
            writefile.write(pybelmol)
        writefile.close()
        if pgstep == "last":
            nr_aggs = obmol.NumConformers()
            obmol = _get_pg(obmol, saveobmol, getb("subgroups"), not(getb("exclude_c1")), pgfilename, progress, postalign, pgregex)
            if obmol.NumConformers()>nr_aggs:
                raise RuntimeError("Number of conformers increased (%d -> %d) during pointgroup screening."%(nr_aggs,obmol.NumConformers()))
