"""
resp2.py
A template to create RESP2 charges

Handles the primary functions
"""
import sys, os
from openeye import oechem
from openeye import oeomega
import openforcefield
import subprocess
import logging as log
from create_mol2_pdb import run_create_mol2_pdb
import openbabel
import shutil

"""
Project functions
create_data.csv
create_target -- uses create data // create mol2 pdb.py
create_fb_input
DONE 1 create_conformers
create_respyte_input (includes folder structure and all files)
optimize_structures
run_respyte
make_resp1
make_resp2 both using run_respyte
functions from correct_atomtypes
"""



def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


def create_respyte(type='RESP1', name ='' ,resname ='MOL',number_of_conformers=1):
    """
    :param type:
        Allowed values RESP1 RESP2GAS RESP2LIQUID
    :return:
    """
    #1 Create folder structure for respyte

    foldername = name+'-'+type
    try:
        os.mkdir(foldername)
    except:
        log.warning('folder {} already exists'.format(foldername))

    input_folder=os.path.join(foldername,'input')
    molecule_folder=os.path.join(input_folder,'molecule')
    mol_folder=os.path.join(molecule_folder,'mol1')
    try:
        os.mkdir(input_folder)
    except:
        log.warning('folder {} already exists'.format(input_folder))

    try:
        os.mkdir(molecule_folder)
    except:
        log.warning('folder {} already exists'.format(molecule_folder))


    try:
        os.mkdir(mol_folder)
    except:
        log.warning('folder {} already exists'.format(mol_folder))

    for i in range(number_of_conformers):
        conf_folder=os.path.join(mol_folder,'conf'+str(i))
        try:
            os.mkdir(conf_folder)
        except:
            log.warning('folder {} already exists'.format(conf_folder))

    log.info('Create folder structure for {} with {} conformers'.format(name,number_of_conformers))
    create_respyte_input_files(type=type, name =name ,resname =resname ,number_of_conformers=number_of_conformers)


    """
    Outsourced in optimization.yml
    #2 Convert mol2 files to xyz files and put them in the corresponding folder
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "xyz")

    foldername = os.path.dirname(name)
    for i in range(number_of_conformers):
        inputfile = os.path.join(foldername, resname + '-conformers_'+str(i)+'.mol2')
        outputfile = os.path.join(foldername, resname + '-conformers_'+str(i)+'.xyz')

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, inputfile)
        obConversion.WriteFile(mol,outputfile)
    """




    #3 Copy optimized files
    for i in range(number_of_conformers):
        foldername = os.path.dirname(name)
        shutil.copyfile(os.path.join(foldername,resname+'-confermers_opt_'+str(i)+'.xyz'), os.path.join('{}-{}/input/molecule/mol1/conf{}/mol1_conf{}.xyz'.format(name,type,i,i)))


    #4 Generate input files

    create_respyte_input_files(type=type, name =name ,resname =resname,number_of_conformers=number_of_conformers)

    return 0


def create_respyte_input_files(type='RESP1', name ='' ,resname ='MOL',number_of_conformers=1):

    if type == 'RESP1':
        method = 'HF'
        basis = '6-31G*'
        pcm = 'N'
    elif type =='RESP2GAS':

        method = 'PW6B95'
        basis = 'aug-cc-pV(D+d)Z'
        pcm = 'N'
    elif type =='RESP2LIQUID':
        method = 'PW6B95'
        basis = 'aug-cc-pV(D+d)Z'
        pcm = 'Y\n    solvent   : water'
    else:
        log.error('Charge type not recognized')
        sys.exit()

    #input.yml
    input_file = open('{}-{}/input/input.yml'.format(name,type),'w')
    input_file.write("""# This file is a sample input.yml, input file for esp_generator.py
# 'input.yml' should be included in 'input' directory with 'molecules' directory to run esp_generator.

# In molecules, the number of molecules and the number of conformers for each molecule are assigned. 
# For example, if the user wants to generate espf data for 2 different molecules with 5 conformers each, 
# molecules is like below. (Each coordinate file should be either PDB or Mol2 file format and should be 
# located in molecules/mol(i)/mol(j)/ ((i), (j) is integer) with its filename, mol(i)_mol(j).pdb(or .mol2)  
# molecules : 
#     mol1 : 5
#     mol2 : 5
molecules:
    mol1 : {}
# In charges, the user should provide charge of each molecule. If not, it consider all molecules as neutral.
# if mol1 is neutral and mol2 has -1 charge,  
# charges :
#     mol1 :  0
#     mol2 : -1
charges :
    mol1 : 0
# In Cheminformatics, user can choose 'openeye' or 'rdkit' or 'None' 
# But for now,  'None' option is not implemented yet, please consider use open-source 'rdkit' instead'
cheminformatics : openeye

# In grid_Setting, (explanation) 
grid_setting :
    forcegen  : Y
    type      : msk # msk(default)/ extendedmsk/ fcc/ newfcc/ vdwfactors/ vdwconstants
    radii     : bondi # bondi(default)/ modbondi
    method    : {}
    basis     : {}
    pcm       : {}
    space     : 0.4
    innner    : 1.6
    outer     : 2.1

    
    
    """.format(number_of_conformers,method,basis,pcm))


    input_file.close()
    #respyte.yml
    respyte_file = open('{}-{}/input/respyte.yml'.format(name,type),'w')

    respyte_file.write("""## This file is a sample respyte.yml, input file for resp_optimzer.py
    ## 'respyte.yml' should be included in 'input' directory', generated from running esp_generator.

    ## Data structure (generated from running esp_generator.py or can manually be generated with an appropriate data structure):
    ## input/
    ## |----respyte.yml
    ## |----molecules/
    ##      |----mol1/
    ##           |----conf1/
    ##                |----mol1_conf1.pdb(or .mol2) , mol1_conf1.espf
    ##           |----conf2/
    ##                |----mol1_conf2.pdb(or .mol2) , mol1_conf2.espf
    ## if xyz format has been used for esp generator, the data structure from esp_generator.py is containing pdb format converted from input xyz.

    ## In molecules, the number of molecules and the number of conformers for each molecule are assigned.
    ## For example, if the user wants to generate espf data for 2 different molecules with 5 conformers each,
    ## molecules is like below. (Each coordinate file should be either PDB or Mol2 file format and should be
    ## located in molecules/mol(i)/mol(j)/ ((i), (j) is integer) with its filename, mol(i)_mol(j).pdb(or .mol2)
    ## molecules :
    ##     mol1 : 5
    ##     mol2 : 5
    molecules :
        mol1 : {}
    ## In charges, user should specify a total charge of each residue. 
    ## If the input molecule is small molecule, should specify net charge of each species.
    charges :
        mol1 : 0

    ## In Cheminformatics, user can choose 'openeye' or 'rdkit' or 'None'
    ## For now, rdkit option is not implemented yet.
    cheminformatics : openeye

    ## In boundary_select, user can select subset of grid.dat
    ## If boundary_select is not set, it uses all grid points in grid.dat for fitting.
    boundary_select:
        radii    : bondi
        inner    : 1.3
        outer    : 2.1

    ## In restraint, user can specify which restraint(model2, model3 and two-stage fit) to use. 
    restraint :
        # penalty : model2
        # a       : 0.005

        # penalty  : model3
        # matrices :
        #     - esp
        # #    - ef
        # a        : 0.0005
        # b        : 0.1

        penalty : 2-stg-fit
        matrices :
            - esp
        a1      : 0.0005
        a2      : 0.001
        b       : 0.1

    ## User can set user-defined constraints on specific charges 



        """.format(number_of_conformers))
    respyte_file.close()

    return 0


def create_conformers(infile = None, outfile= None):

    ifs = oechem.oemolistream()
    if not ifs.open(infile):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % infile)

    ofs = oechem.oemolostream()
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)

    if not oechem.OEIs3DFormat(ofs.GetFormat()):
        oechem.OEThrow.Fatal("Invalid output file format for 3D coordinates!")

    omegaOpts = oeomega.OEOmegaOptions()
    omega = oeomega.OEOmega(omegaOpts)
    omega.SetCommentEnergy(True)
    omega.SetEnumNitrogen(True)
    omega.SetSampleHydrogens(True)
    omega.SetEnergyWindow(10.0)
    omega.SetMaxConfs(5)
    omega.SetRangeIncrement(3)
    omega.SetRMSRange([0.5,1.0,1.5,2.0,2.5,3.0,3.5])
    folder=os.path.dirname(outfile)
    filename= os.path.basename(outfile).split('.')[0]
    for mol in ifs.GetOEMols():
        ret_code = omega.Build(mol)
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            oechem.OEWriteMolecule(ofs, mol)
            for k,conf in enumerate(mol.GetConfs()):
                ofs2 = oechem.oemolostream()
                if not ofs2.open(os.path.join(folder,filename+'_'+str(k)+'.mol2')):
                    oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)
                oechem.OEWriteMolecule(ofs2, conf)
                nconf = k+1
            log.info('Created conformations for {} and saved them to {}'.format(infile, outfile))

        else:
            oechem.OEThrow.Warning("%s: %s" % (mol.GetTitle(), oeomega.OEGetOmegaError(ret_code)))

    return nconf

def create_std_target_file(name='',density=None,hov=None,dielectric=None):
    header_csv = '''# This is documentation for the ForceBalance condensed phase reference data file,,,,,,
    "# Lines beginning with octothorpe are comments, empty lines are ignored",,,,,,
    "# A line can either be a comment, a global parameter, a single line with column headings, or a data line (after the column heading line)",,,,,,
    # This file should be saved as .xlsx (to preserve formatting) but exported to .csv for ForceBalance to use.,,,,,,
    ,,,,,,
    "# Global parameters are defined here; they have ""Global"" in the first column.",,,,,,
    # rho_denom : least squares denominator for the density objective function,,,,,,
    "# Note: w_rho, w_hvap etc. is set in the input file.  This is because it's considered an adjustable option rather than a property of the data set.",,,,,,
    # The overall prefactor for an observable is w_obs / obs_denom^2 ,,,,,,
    ,,,,,,
    Global,rho_denom,5,,,,
    Global,hvap_denom,0.5,,,,
    Global,alpha_denom,1,,,,
    Global,kappa_denom,5,,,,
    Global,cp_denom,2,,,,
    Global,eps0_denom,2,,,,
    Global,use_cvib_intra,FALSE,,,,
    Global,use_cvib_inter,FALSE,,,,
    Global,use_cni,FALSE,,,,
    ,,,,,,
    "# Before entering the data, there must be one line with column headings.",,,,,,
    # Allowed column headings are: (not case sensitive),,,,,,
    # T : Temperature in Kelvin,,,,,,
    "# P : Pressure in specified unit (append atm or bar), default is atm",,,,,,
    # MBAR : Include this state in multistate Bennett acceptance ratio calculation (more extreme states are not included),,,,,,
    # Rho : Density in kg m^-3,,,,,,
    # Hvap : Enthalpy of vaporization in kJ/mol.  Note: simulating liquid at 1 atm is an acceptable approximation for the vapor pressure.,,,,,,
    # Alpha : Thermal expansion coefficient in 10^-4 K^-1.,,,,,,
    # Kappa: Isothermal compressibility in 10^-6 bar^-1.,,,,,,
    # Cp: Isobaric heat capacity in kJ/mol K^-1.,,,,,,
    # Eps0: Static dielectric constant.,,,,,,
    # Rho_wt : Weight to use for this phase point. Analogous for Hvap etc. If column is missing then all weights are 1.,,,,,,
    # Cvib_intra : To be ADDED to the calculated enthalpy of vaporization ; energy difference due to intramolecular vibrational frequency shifts in going from gas to liquid,,,,,,
    # Cvib_inter : To be ADDED to the calculated enthalpy of vaporization ; energy difference due to intermolecular vibrational frequency shifts in going from gas to liquid,,,,,,
    # Cni : To be ADDED to the calculated enthalpy of vaporization ; accounts for the nonideality of the gas phase (more relevant close to critical point),,,,,,
    # dEvib : To be ADDED to the calculated heat capacity.  Note that these values are slightly different from the TIP4P-Ew paper for some reason.,,,,,,
    ,,,,,,
    '''
    if name == '':
        log.error('You did not specify a name for the target folder. Please use create_std_target_file(name=targetname,density=targetdensity,hov=target_heat,dielectric=target_eps0)')
    f=open(name+'-liquid/data.csv','w')
    f.write(header_csv)
    if dielectric == None:
        f.write('T,P,MBAR,Rho,Rho_wt,Hvap,Hvap_wt\n')
        f.write('298.0,1.0 atm, FALSE,{},1,{},1\n'.format(density,hov))
    else:
        f.write('T,P,MBAR,Rho,Rho_wt,Hvap,Hvap_wt,eps0,eps0_wt\n')
        f.write('298.0,1.0 atm, FALSE,{},1,{},1,{},1\n'.format(density,hov,dielectric))

    f.close()
    log.info('''Created target file {}-liquid/data.csv
    Density: {} g/l
    Heat of Vaporization: {} kcal/mol
    Dielectric Constant: {}'''.format(name,density,hov,dielectric))
    return 0

def create_target(smiles='',name='',density=None,hov=None,dielectric=None,resname ='MOL',nmol =700, tries =1000):
    try:
        os.mkdir(name+'-liquid')
    except:
        log.warning('folder {}-liquid already exists'.format(name))
    foldername = name+'-liquid/'
    create_std_target_file(name=name,density=density,hov=hov,dielectric=dielectric)
    create_smifile_from_string(smiles=smiles, filename=foldername+resname+'.smi',)
    run_create_mol2_pdb(nmol = nmol, density = density - 100, tries = tries, input=foldername+resname+'.smi', resname=resname )
    #os.chdir(name+'-liquid')
    return 0


def create_smifile_from_string(smiles='',filename=''):
    f=open(filename,'w')
    f.write(smiles)
    f.close()

    return 0

def optimize_conformers(opt=True, name ='' ,resname ='MOL',number_of_conformers=1):
    header="""memory 12 gb
molecule mol {
noreorient
nocom
    """
    tail_m1="""
}
set basis 6-31G*
optimize('HF')
set basis cc-pV(D+d)Z
optimize('HF')
set basis cc-pV(D+d)Z
optimize('PW6B95')

"""

    #2 Convert mol2 files to xyz files and put them in the corresponding folder
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "xyz")

    foldername = os.path.dirname(name)
    filename = os.path.basename(name)
    for i in range(number_of_conformers):
        inputfile = os.path.join(foldername, resname + '-conformers_'+str(i)+'.mol2')
        outputfile = os.path.join(foldername, resname + '-conformers_'+str(i)+'.xyz')

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, inputfile)
        obConversion.WriteFile(mol,outputfile)

    if opt == True:

        for i in range(number_of_conformers):
            xyz_file = os.path.join(foldername, resname + '-conformers_'+str(i)+'.xyz')
            psi4_input_file =os.path.join(foldername, resname + '-conformers_'+str(i)+'.in')
            psi4_output_file =os.path.join(foldername, resname + '-conformers_'+str(i)+'.out')
            f=open(xyz_file,'r')
            coordinates=f.readlines()[2:]
            f.close()
            f=open(psi4_input_file,'w')
            f.write(header)
            f.write('0 1\n')
            for line in coordinates:
                f.write(line)
            f.write(tail_m1)
            f.write("mol.save_xyz_file('{}',True)".format(os.path.join(foldername,resname+'-confermers_opt_'+str(i)+'.xyz')))

            f.close()

            os.system('psi4 {} -n 4'.format(psi4_input_file))
            if 'beer' in open(psi4_output_file).read():
                log.info('Optimization of {} and conformer {} succesful'.format(filename,i))

    else:
        for i in range(number_of_conformers):
            shutil.copy(os.path.join(foldername,resname+'-confermers_'+str(i)+'.xyz'),os.path.join(foldername,resname+'-confermers_opt_'+str(i)+'.xyz'))


if __name__ == "__main__":
    log.getLogger().setLevel(log.INFO)
    create_target(name ='data/ethanol',density=789.3,hov = 42.3, dielectric=32.7,smiles='CCO',resname ='ETH')
    number_of_conformers=create_conformers(infile='data/ETH.mol2', outfile='data/ETH-conformers.mol2')
    optimize_conformers(name ='data/ethanol',resname ='ETH', opt = True,number_of_conformers=number_of_conformers)
    create_respyte(name ='data/ethanol',resname ='ETH', type='RESP1',number_of_conformers=number_of_conformers)


