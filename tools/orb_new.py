import numpy as np
import matplotlib.pyplot as plt
import os

class molecule():
    def __init__(self, atoms = [], totchg = 0., mocoeff = []):
        self.atoms = atoms
        self.totchg = totchg
        self.mocoeff = mocoeff

    def add_atoms(self,atome):
        self.atoms.append(atome)


class basis():
    def __init__(self, listecons, l, ao=[100]):
        self.listecons = listecons
        self.l = l
        self.ao = ao


class atom():
    def __init__(self,   bases = [], pos=[0.,0.,0.], chg=1.0, atype='H'):
        self.atype = atype
        self.chg = chg
        self.pos = pos
        self.bases = bases

    def add_basis(self,basis):
        self.bases.append(basis)
REP = "../input/CO2/800nm_I5E13_NX_200_1MO_2e/BASIS/CO2/1/"
xmax=6
zmax=6

x=np.arange(-xmax,xmax,0.01)
z=np.arange(-zmax,zmax,0.01)

xv, zv = np.meshgrid(x, z, sparse=False, indexing='ij')

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
                line = line + 4
                print("test line ", argosls[line])
                ao = [argosls[line].split()[1]]
                l = 0
            if 'p' in TYPE:
                line = line + 6
                ao = [argosls[line].split()[1],argosls[line+1].split()[1],argosls[line+2].split()[1]]
                line = line + 2
                l = 1
            line = line + 2
            return line, orb, TYPE, ao, l
        
f = open(REP+'argosls.sp', "r")
argosls = f.readlines()
f.close()
# title = argosls[2]
#print(str(mocoeff[1]))
#print("!!!!!!!!!!!!! TTTTESSSSSSSSSSSSTTTTTTT")
#print(argosls[3].split())
line=3
NGEN, NS, NAORDS, NCONS, NGCS, ITOL, ICUT, AOINTS, ONLY1E, \
 INRM, NCRS, L1REC, L2REC, AOINT2, FSPLIT = argosls[line].split()

NAORDS = int(NAORDS)

fname=REP+'mocoef_scf.sp'
with open(fname) as f:
    content = f.readlines()
content = [x.strip() for x in content]
nmocoeffs = int(content[6].split()[0])

print(nmocoeffs)
mocoeff = np.zeros((nmocoeffs, nmocoeffs))
norb=0
morb=0
ncont=10
while morb < nmocoeffs:
#while morb<27:
    while norb < nmocoeffs:
        for number in content[ncont].split():
            print(morb,norb,number)
            mocoeff[morb][norb] = number
            norb = norb + 1
        ncont = ncont + 1
    morb=morb+1
    norb=0

f.close()
while not ("nuclear repulsion energy" in argosls[line]):
    line = line + 1

#while not ("atoms" in argosls[line]):
#    line = line + 1


mol = molecule()
atome = atom()

def searchorb(atome,argosls,line,mol):
    if "orbitals" in argosls[line]:
        line,orb,TYPE,ao, l = get_orb(argosls, line)
        atome.add_basis(basis(listecons=orb, l=l, ao=ao))
        print(orb)
        print("line = ", line)
        print(argosls[line])
        if "orbitals" in argosls[line]:
            searchorb(atome,argosls,line,mol)
        else :
            mol.add_atoms(atome)
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
            
            atome = atom(atype=ATYPE, pos = [X,Y,Z], chg = CHG, bases=[])
            print("test new basis", atome.bases )
            line = line + 2
            searchorb(atome,argosls,line,mol)
            
            print("testing atom add in molecule")
            for atome in mol.atoms:
                print(atome.atype)
        else:
            #mol.add_atoms(atome)
            print("What else???")
            exit 
       
         

searchatomorb(argosls,line,atome,mol)
aoorbs = []
print("Testing Molecule")
m = 0
p = 0
for n in range(len(mol.atoms)):
    atome = mol.atoms[n]
    print(atome.atype)
    print(atome.pos)
    for base in atome.bases:

        print(base.listecons)
        if base.l == 0:
            orb = np.zeros(xv.shape)
            for con in base.listecons:
                orb = orb + con[1] * np.exp(-con[0] * ((xv-atome.pos[0]) ** 2.0 + (zv-atome.pos[2]) ** 2.0))
            aoorbs.append(orb)

            plt.contourf(xv, zv, orb)
            plt.colorbar()
            plt.savefig('aoorb_' + str(m) + '.png')
            plt.close()
            m = m + 1
        if base.l == 1:
            for iao in base.ao:
                if iao == "100":
                    mult = xv
                elif iao == "010":
                    mult = 0
                elif iao == "001":
                    mult = zv
                orb = np.zeros(xv.shape)
                for con in base.listecons:
                    orb = orb + con[1] * np.exp(-con[0] * ((xv - atome.pos[0]) ** 2.0 + (zv - atome.pos[2]) ** 2.0))
                orb = orb * mult
                plt.contourf(xv, zv, orb)
                plt.colorbar()
                plt.savefig('aoorb_' + str(m) + '.png')
                plt.close()
                aoorbs.append(orb)
                m = m + 1

        #print('Min, Max')
        #print(min(orb))
        #print(max(orb))

moorb = []
#print("printing moceoeoof")
m = 0
for mocoefs in mocoeff:
    orb = np.zeros(xv.shape)
    print("MOCEOF (",m,") = ",mocoefs)
    print("Sum**2 = ",np.sum(mocoefs)**2)
    norm = 0
    for mocoef in mocoefs:
        orb = orb + float(mocoef)*aoorbs[m]
        norm = norm + mocoef**2
    print("norm =" , norm)
    plt.contourf(xv, zv, orb)
    plt.colorbar()
    plt.savefig('moorb_' + str(m) + '.png')
    plt.close()
    m = m + 1
    #for mocoef in mocoeff:


#        print(atome.atype)
#        print(atome.pos)
#        n = 0
#        for base in atome.basis:
#            print(base.listecons)#
#
#            for con in base.listecons:
#                if base.l == 0:
#                    orb = orb + con[1] * np.exp(-con[0] * ((xv-atome.pos[0]) ** 2.0 + (zv-atome.pos[2]) ** 2.0))
#
#            plt.contourf(xv, zv, orb)
#            plt.savefig('orb_'+str(n)+'.png')
#        n = n + 1