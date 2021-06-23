import numpy as np
import matplotlib.pyplot as plt
import os

class molecule():
    def __init__(self, atoms = [], totchg = 0., mocoeff = []):
        self.atoms = atoms
        self.totchg = totchg
        self.mocoeff = mocoeff

    def add_atoms(self,atom):
        self.atoms.append(atom)


class basis():
    def __init__(self, listecons, l):
        self.listecons = listecons
        self.l = l


class atom():
    def __init__(self,   basis = [], pos=[0.,0.,0.], chg=1.0, type='H'):
        self.type = type
        self.chg = chg
        self.pos = pos
        self.basis = basis

    def add_basis(self,basis):
        self.basis.append(basis)
REP = "../input/CO2/800nm_I5E13_NX_200_1MO_2e/BASIS/CO2/1/"
xmax=6
zmax=6

x=np.arange(-xmax,xmax,0.01)
z=np.arange(-zmax,zmax,0.01)

mocoeff=np.zeros((27,27))
fname=REP+'mocoef_scf.sp'
with open(fname) as f:
    content = f.readlines()
content = [x.strip() for x in content]
norb=0
morb=0
ncont=10
while morb<27:
 while norb<27:
  for number in content[ncont].split():
    mocoeff[morb][norb]=number
    norb=norb+1
  ncont=ncont+1
 morb=morb+1
 norb=0

def get_orb(argosls, line):
    if "orbitals" in argosls[line]:
        orb = []
        TYPE = argosls[line].split()[0]
        print("Orbitals Type = ",TYPE)
        line = line + 3
        while not ("symmetry orbital" in argosls[line+1]):
            print(argosls[line])
            orb.append([float(argosls[line].split()[0]),float(argosls[line].split()[1])])
            line = line + 1
            #print("Line + 2 = ",argosls[line+2])
        else: 
            line = line + 2
            if 's' in TYPE:
                ao = [argosls[line].split()[1]]
                l = 0
            if 'p' in TYPE:
                ao = [argosls[line].split()[1],argosls[line+1].split()[1],argosls[line+2].split()[1]]
                line = line + 2
                l = 1
            line = line + 2
            return line, orb, TYPE, ao, l
        
f = open(REP+'argosls.sp', "r")
argosls = f.readlines()
# title = argosls[2]
#print(str(mocoeff[1]))
#print("!!!!!!!!!!!!! TTTTESSSSSSSSSSSSTTTTTTT")
#print(argosls[3].split())
line=3
NGEN, NS, NAORDS, NCONS, NGCS, ITOL, ICUT, AOINTS, ONLY1E, \
 INRM, NCRS, L1REC, L2REC, AOINT2, FSPLIT = argosls[line].split()
 
 
while not ("nuclear repulsion energy" in argosls[line]):
    line = line + 1

#while not ("atoms" in argosls[line]):
#    line = line + 1


mol = molecule()
atome = atom()

def searchorb(atome,argosls,line,mol):
    if "orbitals" in argosls[line]:
        line,orb,TYPE,ao, l = get_orb(argosls, line)
        atome.add_basis(basis(listecons = orb, l=l))
        print(orb)
        print("line = ", line)
        print(argosls[line])
        if "orbitals" in argosls[line]:
            searchorb(atome,argosls,line,mol)
        else :
            searchatomorb(argosls,line,atome,mol)

def searchatomorb(argosls,line,atome,mol):
    while line < len(argosls)-1 and (("orbitals" in argosls[line]) or not ("atoms" in argosls[line])):
        line = line + 1
    else : 
        if ("atoms" in argosls[line]):
            
            ATYPE = argosls[line].split()[0]
            line = line + 2
            CHG = argosls[line].split()[2]
            line = line + 3

            ATOMNUM = argosls[line].split()[0]
            X = float(argosls[line].split()[1]) 
            Y = float(argosls[line].split()[2]) 
            Z = float(argosls[line].split()[3]) 
            print(ATYPE)
            print(CHG)
            print(ATOMNUM)
            print(X,Y,Z)
            
            atome = atom(type=ATYPE, pos = [X,Y,Z], chg = CHG)
            
            line = line + 2
            searchorb(atome,argosls,line,mol)
            mol.add_atoms(atome)
        else:
            #mol.add_atoms(atome)
            print("What else???")
            exit 
       
         

searchatomorb(argosls,line,atome,mol)


for atome in mol.atoms:
    print(atome.type)
    print(atome.pos)
    for base in atome.basis:
        print(base.listecons)