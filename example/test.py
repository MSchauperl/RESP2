import resp2
import logging as log

log.getLogger().setLevel(log.INFO)
resp2.create_target(name ='ethanol',density=789.3,hov = 42.3, dielectric=32.7,smiles='CCO',resname ='ETH')
number_of_conformers = resp2.create_conformers(infile='ethanol-liquid/ETH.mol2', outfile='ethanol-liquid/ETH-conformers.mol2')
resp2.optimize_conformers(name ='ethanol',resname ='ETH', opt = True,number_of_conformers=number_of_conformers)
resp2.create_respyte(name ='ethanol',resname ='ETH', type='RESP2LIQUID',number_of_conformers=number_of_conformers)
resp2.create_respyte(name ='ethanol',resname ='ETH', type='RESP2GAS',number_of_conformers=number_of_conformers)
resp2.create_respyte(name ='ethanol',resname ='ETH', type='RESP1',number_of_conformers=number_of_conformers)
resp2.create_charge_file(name='ethanol',resname ='ETH', type = 'RESP1', delta= 2.0)
resp2.create_charge_file(name='ethanol',resname ='ETH', type = 'RESP2', delta= 0.5)
