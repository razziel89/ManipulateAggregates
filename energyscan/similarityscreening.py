import sys

import openbabel as op

from energyscan.scan import _prepare_molecules as prepare_molecules
from manipulate_molecules import read_from_file

def similarityscreening_main(parser):
    gets   = parser.get_str
    geti   = parser.get_int
    getf   = parser.get_float
    getb   = parser.get_boolean
    do_calculate = not(getb("config_check"))

    mol1 = read_from_file(gets("geometry1"),ff=None)
    mol2 = read_from_file(gets("geometry2"),ff=None)

    obmol = prepare_molecules(mol1,mol2)

    std_map = op.StdMapStringString()

    #add the appropriate configuration paramters to the std::map<std::string,std::string>
    std_map['rcutoff'] = str(getf("rmsd_min"))      #this way I can be sure it's actually a floating point number
    std_map['ecutoff'] = str(geti("energy_cutoff")) #this way I can be sure it's actually an integer
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

    #to avoid segfaults, define some input parameters that would normally be given via the command-line
    in_out_options = op.OBConversion()
    in_out_options.SetInFormat('nul')
    in_out_options.SetOutFormat('nul')
  
    #create own instance of the "SimScreen"-object
    simscreen = op.OBOp.FindType('simscreen')
    simscreen.Do(obmol,'', std_map, in_out_options)

    print >>sys.stderr, "WARNING: the jobtype 'smimilarityscreening' has not yet been fully implemented."
