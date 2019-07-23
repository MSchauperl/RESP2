import sys



def charge_scaling(mol2_with_charges=None, resname = 'MOL1'):
    f = open(mol2_with_charges, 'r')
    v = 0
    num = 0
    lines = f.readlines()
    if lines[1].startswith('***') or lines[1].startswith('UNL'):
        lines[1]='{}\n'.format(name)
    for line in lines:
        if '@<TRIPOS>ATOM' in line:
            v = 1
            output.write(line)
        elif '@<TRIPOS>BOND' in line:
            v = 2
            output.write(line)
        elif v == 1:
            entry = line.split()
            if float(entry[8])<0.0:
                output.write("{:>7} {:<3}{:>15}{:>10}{:>10} {:<3}{:>8}{:>5}{:>14} # EVAL 8 0.000*(1.0-PRM['parameters.txt:1.0']){}*PRM['parameters.txt:1.0']\n".format(entry[0], entry[1], entry[2], entry[3],
                                                                            entry[4], entry[5], entry[6], entry[7],
                                                                          entry[8],entry[8]))
            else:
                output.write("{:>7} {:<3}{:>15}{:>10}{:>10} {:<3}{:>8}{:>5}{:>14} # EVAL 8 0.000*(1.0-PRM['parameters.txt:1.0'])+{}*PRM['parameters.txt:1.0']\n".format(entry[0], entry[1], entry[2], entry[3],
                                                                            entry[4], entry[5], entry[6], entry[7],
                                                                          entry[8],entry[8]))
            num+=1
        else:
            output.write(line)
    f.close()

def delta_resp2(mol2_gas = None, mol2_liquid = None, resname = 'MOL1'):
    f = open(mol2_gas, 'r')
    f2 = open(mol2_liquid, 'r')
    v = 0

    v = 0
    # Read in implicit solvent charges (isc)
    isc = []
    lines = f2.readlines()
    for line in lines:
        if '@<TRIPOS>ATOM' in line:
            v = 1
        elif '@<TRIPOS>BOND' in line:
            v = 2
        elif v == 1:
            entry = line.split()
            isc.append(float(entry[8]))
    f2.close()


    v = 0
    num = 0
    lines = f.readlines()
    if lines[1].startswith('***') or lines[1].startswith('resp_gas') or lines[1].startswith('mol1_conf1'):
        lines[1]='{}\n'.format(resname)
    for i,line in enumerate(lines):
        if '@<TRIPOS>ATOM' in line:
            v = 1
            output.write(line)
        elif '@<TRIPOS>BOND' in line:
            v = 2
            output.write(line)
        elif v == 1:
            entry = line.split()
            if float(isc[num])<0.0:
                output.write("{:>7} {:<3}{:>15}{:>10}{:>10} {:<3}{:>8}{:>5}{:>14} # EVAL 8 {}*(1.0-PRM['parameters.txt:1.0']){}*PRM['parameters.txt:1.0']\n".format(entry[0], entry[1], entry[2], entry[3],
                                                                            entry[4], entry[5], entry[6], resname,
                                                                          entry[8],entry[8],isc[num]))
            else:
                output.write("{:>7} {:<3}{:>15}{:>10}{:>10} {:<3}{:>8}{:>5}{:>14} # EVAL 8 {}*(1.0-PRM['parameters.txt:1.0'])+{}*PRM['parameters.txt:1.0']\n".format(entry[0], entry[1], entry[2], entry[3],
                                                                            entry[4], entry[5], entry[6], resname,
                                                                          entry[8],entry[8],isc[num]))
            num+=1
        else:
            output.write(line)
    f.close()



def calc_parm(par_file,mol2_file):
    try:
        para=float(par_file)
    except:
        f=open(par_file,'r')
        para=float(f.readline().split()[0])
        f.close()
    print(para)
    f=open(mol2_file,'r')
    lines=f.readlines()
    f.close()
    f=open(mol2_file,'w')
    for line in lines:
        if "EVAL" in line:
            entry=line.split()
            c1=float(entry[12].split("*")[0])
            c2=float(entry[12].split("*")[1].split(')')[1])
            entry[8]=c1*(1.0-para)+c2*para
            f.write('{0:>7} {1:<3} {2:>14} {3:>9} {4:>9} {5:<4} {6:>6} {7:>4} {8:>13.6f}\n'.format(*entry))
        else:
            f.write(line)

    f.close()

opt='normal'

if __name__=='__main__':
    if len(sys.argv)<3:
        print("Usage: python charge_calculation.py charges.mol2 atomtyp.mol2 -out <output.mol2>--opt <normal/scale/resp2> <implicit_solvent_charges.mol2>")
        exit(0)
    else:
        m1 = sys.argv[1]
        m2 = sys.argv[2]
        for i in range(len(sys.argv)):
            if sys.argv[i]=='--opt':
                 opt=sys.argv[i+1]
                 if sys.argv[i+1]=='resp2':
                     m3=sys.argv[i+2]
            if sys.argv[i]=='-out':
                 out=sys.argv[i+1]
                 output=open(out,'w')
        name=out.split('/')[-1].split('-')[0]

    if opt=='normal':
        combine2mol2(m1, m2)
    if opt=='scale':
        charge_scaling(m1, m2)
    if opt=='resp2':
        name=out.split('/')[-1].split('-')[0]
        delta_resp2(m1, m2, m3)
    output.close()
