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
import resp2.create_mol2_pdb
import openbabel
import shutil


"""
Project functions
XXXcreate_data.csv
XXXcreate_target -- uses create data // create mol2 pdb.py
create_fb_input
XXXXDONE 1 create_conformers
XXXXcreate_respyte_input (includes folder structure and all files)
XXXXXoptimize_structures
XXXXXrun_respyte
XXXXXXmake_resp1
XXXXXXXmake_resp2 both using run_respyte
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
    molecule_folder=os.path.join(input_folder,'molecules')
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

    for i in range(1,number_of_conformers+1):
        conf_folder=os.path.join(mol_folder,'conf'+str(i))
        try:
            os.mkdir(conf_folder)
        except:
            log.warning('folder {} already exists'.format(conf_folder))

    log.info('Create folder structure for {} with {} conformers'.format(name,number_of_conformers))
    create_respyte_input_files(type=type, name =name ,resname =resname ,number_of_conformers=number_of_conformers)




    #3 Copy optimized files
    for i in range(1,number_of_conformers+1):
        foldername = os.path.join(os.path.dirname(name),name+'-liquid')
        shutil.copyfile(os.path.join(foldername,resname+'-confermers_opt_'+str(i)+'.xyz'), os.path.join('{}-{}/input/molecules/mol1/conf{}/mol1_conf{}.xyz'.format(name,type,i,i)))


    #4 Generate input files

    create_respyte_input_files(type=type, name =name ,resname =resname,number_of_conformers=number_of_conformers)

    cwdir = os.getcwd()
    #5 Run RESPyt and PSI4
    #os.chdir(foldername)
    calculate_respyte(type=type, name =name ,resname =resname,number_of_conformers=number_of_conformers)
    #os.chdir(cwdir)
    return 0


def calculate_respyte(type='RESP1', name ='' ,resname ='MOL',number_of_conformers=1):
    foldername = name + '-' + type
    cwdr = os.getcwd()
    os.chdir(foldername)
    mol_folder = 'input/molecules/mol1/'
    for i in range(1,number_of_conformers+1):
        conf_folder=os.path.join(mol_folder,'conf'+str(i))
        tmp_folder= os.path.join(conf_folder,'tmp/')
        try:
            shutil.rmtree(tmp_folder)
        except:
            pass
    os.system('python ~/programs/respyte/respyte/esp_generator.py')
    for i in range(1,number_of_conformers+1):
        conf_folder=os.path.join(mol_folder,'conf'+str(i))
        psi4_output_file = os.path.join(conf_folder,'tmp/output.dat')
        if 'beer' in open(psi4_output_file).read():
            log.info('ESP calculation for {} and conformer {} succesful'.format(name, i))
        else:
            log.error('ESP calculation for {} and conformer {} FAILED!!!!!!'.format(name, i))

    os.system('python ~/programs/respyte/respyte/resp_optimizer.py')
    os.chdir(cwdr)
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
    input_file.write("""molecules:
    mol1 : {}
charges :
    mol1 : 0
cheminformatics : openeye

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

    respyte_file.write("""
    molecules :
        mol1 : {}
    charges :
        mol1 : 0

    cheminformatics : openeye

    boundary_select:
        radii    : bondi
        inner    : 1.3
        outer    : 2.1

    restraint :
        penalty : 2-stg-fit
        matrices :
            - esp
        a1      : 0.0005
        a2      : 0.001
        b       : 0.1

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
                if not ofs2.open(os.path.join(folder,filename+'_'+str(k+1)+'.mol2')):
                    oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)
                oechem.OEWriteMolecule(ofs2, conf)
                nconf = k+1
            log.info('Created conformations for {} and saved them to {}'.format(infile, outfile))

        else:
            oechem.OEThrow.Warning("%s: %s" % (mol.GetTitle(), oeomega.OEGetOmegaError(ret_code)))

    return nconf

def create_std_target_file(name='',density=None,hov=None,dielectric=None):
    header_csv = '''# This is documentation for the ForceBalance condensed phase reference data file
,,,,,,
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
    resp2.create_mol2_pdb.run_create_mol2_pdb(nmol = nmol, density = density - 100, tries = tries, input=foldername+resname+'.smi', resname=resname )
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

    foldername = os.path.join(os.path.dirname(name),name+'-liquid')
    filename = os.path.basename(name)
    for i in range(1,number_of_conformers+1):
        inputfile = os.path.join(foldername, resname + '-conformers_'+str(i)+'.mol2')
        outputfile = os.path.join(foldername, resname + '-conformers_'+str(i)+'.xyz')

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, inputfile)
        obConversion.WriteFile(mol,outputfile)

    if opt == True:

        for i in range(1,number_of_conformers+1):
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
                log.error('Optimization of {} and conformer {} FAILED!!!!!!'.format(filename,i))

    else:
        for i in range(1,number_of_conformers+1):
            if not os.path.exists(os.path.join(foldername,resname+'-confermers_opt_'+str(i)+'.xyz')):
                shutil.copy(os.path.join(foldername,resname+'-confermers_'+str(i)+'.xyz'),os.path.join(foldername,resname+'-confermers_opt_'+str(i)+'.xyz'))

def create_charge_file(name = '', resname= 'MOL', delta = 0.0, type = 'RESP1'):
    if type == 'RESP1':
        mol2_resp1 = name+'-RESP1/resp_output/mol1_conf1.mol2'
        f = open(mol2_resp1)
        output_file= os.path.join(name+'-liquid',resname +'_R1_'+str(int(delta*100))+'.mol2')
        output = open(output_file,'w')

        # Read in RESP1 charges
        v = 0
        resp1charges = []
        lines = f.readlines()
        for line in lines:
            if '@<TRIPOS>ATOM' in line:
                v = 1
            elif '@<TRIPOS>BOND' in line:
                v = 2
            elif v == 1:
                entry = line.split()
                resp1charges.append(float(entry[8]))
        f.close()
        charges = []
        for i in range(len(resp1charges)):
            charges.append(resp1charges[i]*delta)


    elif type == 'RESP2':
        mol2_gas = name+'-RESP2GAS/resp_output/mol1_conf1.mol2'
        mol2_liquid = name+'-RESP2LIQUID/resp_output/mol1_conf1.mol2'
        output_file= os.path.join(name+'-liquid',resname +'_R2_'+str(int(delta*100))+'.mol2')
        f = open(mol2_gas, 'r')
        f2 = open(mol2_liquid, 'r')
        output = open(output_file,'w')

        # Read in gas phase charges (gpc)
        v = 0
        gpc = []
        lines = f.readlines()
        for line in lines:
            if '@<TRIPOS>ATOM' in line:
                v = 1
            elif '@<TRIPOS>BOND' in line:
                v = 2
            elif v == 1:
                entry = line.split()
                gpc.append(float(entry[8]))
        f.close()

        # Read in implicit solvent charges (isc)
        v = 0
        isc = []
        lines2 = f2.readlines()
        for line in lines2:
            if '@<TRIPOS>ATOM' in line:
                v = 1
            elif '@<TRIPOS>BOND' in line:
                v = 2
            elif v == 1:
                entry = line.split()
                isc.append(float(entry[8]))
        f2.close()

        charges = []
        for i in range(len(isc)):
            charges.append(gpc[i]*(1.0-delta)+isc[i]*delta)
        log.info('Created charges {} type charges with a delta value of {}'.format(type,delta))

    else:
        log.error('The type you defined is not recognized. Up to now only RESP1 and RESP2 are valid options')
        sys.exit(1)


    v = 0
    num = 0
    if lines[1].startswith('***') or lines[1].startswith('resp_gas') or lines[1].startswith('mol1_conf1'):
        lines[1] = '{}\n'.format(resname)
    for i, line in enumerate(lines):
        if '@<TRIPOS>ATOM' in line:
            v = 1
            output.write(line)
        elif '@<TRIPOS>BOND' in line:
            v = 2
            output.write(line)
        elif v == 1:
            entry = line.split()
            output.write(
            "{:>7} {:<3}{:>15}{:>10}{:>10} {:<3}{:>8}{:>5}{:>14.4f} \n".format(
                        entry[0], entry[1], entry[2], entry[3],
                        entry[4], entry[5], entry[6], resname,
                        charges[num]))
            num += 1
        else:
            output.write(line)
    f.close()
    output.close()


    return 0

if __name__ == "__main__":
    log.getLogger().setLevel(log.INFO)
    #create_target(name ='data/ethanol',density=789.3,hov = 42.3, dielectric=32.7,smiles='CCO',resname ='ETH')
    number_of_conformers=create_conformers(infile='data/ETH.mol2', outfile='data/ETH-conformers.mol2')
    #optimize_conformers(name ='data/ethanol',resname ='ETH', opt = False,number_of_conformers=number_of_conformers)
    #create_respyte(name ='data/ethanol',resname ='ETH', type='RESP2LIQUID',number_of_conformers=number_of_conformers)
    #create_respyte(name ='data/ethanol',resname ='ETH', type='RESP2GAS',number_of_conformers=number_of_conformers)
    #create_respyte(name ='data/ethanol',resname ='ETH', type='RESP1',number_of_conformers=number_of_conformers)
    #create_charge_file(name='data/ethanol',resname ='ETH', type = 'RESP1', delta= 2.0)




