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
# IDEE : comp[iler le read argos avec f2py
#TODO réécrire sous forme matricielle
x=np.arange(-xmax,xmax,0.01)
z=np.arange(-zmax,zmax,0.01)

orbc1s_mod = np.zeros((len(x),len(z)))
orbc1s=np.zeros((len(x),len(z)))
orbc2s=np.zeros((len(x),len(z)))
orbc3s=np.zeros((len(x),len(z)))
orbc2px=np.zeros((len(x),len(z)))
orbc2py=np.zeros((len(x),len(z)))
orbc2pz=np.zeros((len(x),len(z)))
orbc3px=np.zeros((len(x),len(z)))
orbc3py=np.zeros((len(x),len(z)))
orbc3pz=np.zeros((len(x),len(z)))
orbc3s=np.zeros((len(x),len(z)))
orbO1_1s=np.zeros((len(x),len(z)))
orbO1_2s=np.zeros((len(x),len(z)))
orbO1_3s=np.zeros((len(x),len(z)))
orbO1_2px=np.zeros((len(x),len(z)))
orbO1_2py=np.zeros((len(x),len(z)))
orbO1_2pz=np.zeros((len(x),len(z)))
orbO1_3px=np.zeros((len(x),len(z)))
orbO1_3py=np.zeros((len(x),len(z)))
orbO1_3pz=np.zeros((len(x),len(z)))


orbO2_1s=np.zeros((len(x),len(z)))
orbO2_2s=np.zeros((len(x),len(z)))
orbO2_3s=np.zeros((len(x),len(z)))
orbO2_2px=np.zeros((len(x),len(z)))
orbO2_2py=np.zeros((len(x),len(z)))
orbO2_2pz=np.zeros((len(x),len(z)))
orbO2_3px=np.zeros((len(x),len(z)))
orbO2_3py=np.zeros((len(x),len(z)))
orbO2_3pz=np.zeros((len(x),len(z)))

MO1=np.zeros((len(x),len(z)))
MO2=np.zeros((len(x),len(z)))
MO=np.zeros((27,len(x),len(z)))
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


# TODO read the argosls as defined here :
#https://www.univie.ac.at/columbus/docs_COL70/documentation_main.html
# *******************************************************************************
#
#                                    INPUT DATA
#
# *******************************************************************************
#
#   The input is list directed except where a FORMAT statement is given.  It is
#   recommended that a / be put at the end of each input record, as is done on
#   line 2 and many other lines in the sample input data sets, in order to
#   allow for future extensions.
#
#   1)     TITLE(1)
#          FORMAT(A80)
#
#   2)     NGEN, NS, NAORDS, NCONS, NGCS, ITOL, ICUT, AOINTS, ONLY1E,
#          INRM, NCRS, L1REC, L2REC, AOINT2, FSPLIT
#
#   3)     NST, (ND(I), ITYP(I), I = 1, NST)
#          FORMAT(I3,12(I3,A3))
#
#   4)     NDPT
#
#   5)     DO I = 1, NDPT
#            P1(I), P2(I), P3(I)
#          ENDDO
#
#   6)     DO I = 1, NAORDS
#            NREP(I), (IREP(J), J = 1, NREP(I))
#          ENDDO
#
#   7)     DO I = 1, NGCS
#            ICSU(I), ICTU(I), IAORDS(I)
#            DO J = 1, ICSU(I)
#              (ISOCOEF(K,J), K = 1, ICTU(I))
#            ENDDO
#          ENDDO
#
#   8)     DO I = 1, NCONS
#            ICONU(I), LMNP1(I), NRCR(I)
#            DO J = 1, ICONU(I)
#              ZET(J,I), (ETA(K,J,I), K = 1, NRCR(I))
#            ENDDO
#          ENDDO
#
#   9)     IF(NCRS .NE. 0) THEN
#            DO ICRS = 1, NCRS
#              LCRU, LLSU
#              IF(LCRU .GE. 0) THEN
#                DO L = 0, LCRU
#                  NBFCR
#                  DO K = 1, NBFCR
#                    NCR(K), ZCR(K), CCR(K)
#                  ENDDO
#                ENDDO
#              ENDIF
#              IF(LLSU .GE. 1) THEN
#                DO L = 1, LLSU
#                  NBFCR
#                  DO K = 1, NBFCR
#                    NCR(K), ZCR(K), CCR(K)
#                  ENDDO
#                ENDDO
#              ENDIF
#            ENDDO
#          ENDIF
#
#  10)     DO IS = 1, NS
#            MTYPE(IS), NF(IS), NC(IS), CHG(IS)
#            FORMAT(A3,2I3,F3.0)
#            DO J = 1, NC
#              X(J), Y(J), Z(J)
#            ENDDO
#            IF(NC(IS) .NE. 1) THEN
#              DO J = 1, NGEN
#                (IGEN(K,J), K = 1, NC)
#              ENDDO
#            ENDIF
#            DO J = 1, NF
#              MCONS(J), IGCS(J)
#            ENDDO
#            IF(NCRS .NE. 0) THEN
#              MCRS(IS)
#            ENDIF
#          ENDDO
#
#  11)     DO I = 2, 4
#            TITLE(I)
#            IF(TITLE(I) .EQ. ' ') EXIT
#          ENDDO
f = open(REP+'argosls.sp', "r")
argosls = f.readlines()
# title = argosls[2]
#print(str(mocoeff[1]))
#print("!!!!!!!!!!!!! TTTTESSSSSSSSSSSSTTTTTTT")
#print(argosls[3].split())
NGEN, NS, NAORDS, NCONS, NGCS, ITOL, ICUT, AOINTS, ONLY1E, \
 INRM, NCRS, L1REC, L2REC, AOINT2, FSPLIT = argosls[3].split()
#   2)     NGEN, NS, NAORDS, NCONS, NGCS, ITOL, ICUT, AOINTS, ONLY1E,
#          INRM, NCRS, L1REC, L2REC, AOINT2, FSPLIT
#print(NGEN)
#   3)     NST, (ND(I), ITYP(I), I = 1, NST)
#          FORMAT(I3,12(I3,A3))
#
#   4)     NDPT
NST = int(argosls[4].split()[0])
j=1
#print("NST = ", NST)
ND = []
ITYPE = []
for i in range(NST):
    print("test",i)
    ND.append(argosls[4].split()[j])
    ITYPE.append(argosls[4].split()[j+1])
NDPT = int(argosls[5].split()[0])
# TODO : this if NPDT > 0
#   5)     DO I = 1, NDPT
#            P1(I), P2(I), P3(I)
#          ENDDO


#
#   6)     DO I = 1, NAORDS
#            NREP(I), (IREP(J), J = 1, NREP(I))
#          ENDDO
#
NAORDS = int(NAORDS)
line=6
#print("NAORDS = ",NAORDS)
NREP=[]
IREP=[]
k=0
for i in range(NAORDS):
    NREP.append(int(argosls[6+i].split()[0]))
    for j in range(NREP[k]):
        IREP.append(argosls[6+i].split()[j+1])
    k=k+1
line = line+NAORDS-1
#   7)     DO I = 1, NGCS
#            ICSU(I), ICTU(I), IAORDS(I)
#            DO J = 1, ICSU(I)
#              (ISOCOEF(K,J), K = 1, ICTU(I))
#            ENDDO
#          ENDDO
NGCS = int(NGCS)
ICSU = []
ICTU = []
IAORDS = []
ISOCOEF = []
k = 0

for i in range(NGCS):
    ICSU.append(int(argosls[line+k].split()[0]))
    ICTU.append(int(argosls[line+k].split()[1]))
    IAORDS.append(int(argosls[line+k].split()[2]))

    for j in range(ICSU[k]):
        for l in range(ICTU[k]):
            ISOCOEF.append(argosls[line+k+1].split())
            print(argosls[line+k+1].split())


    print(line)
print(NGCS+ICSU[0])

line=line+NGCS+sum(ICSU)
print("test line")
print(line)
print(NGCS)
line = line - 1
print(argosls[line].split())
print(NCONS)
ZET = []
ETA = []
ICONU = []
LMNP1 = []
NRCR = []
NCONS=int(NCONS)

liste_base = []





for i in range(NCONS):
    a,b,c = argosls[line].split()
    ICONU.append(int(a))
    LMNP1.append(int(b))
    NRCR.append(int(c))
    line = line + 1
    print(a,b,c)
    une_base1 = []
    for j in range(int(a)):
        ZET.append(argosls[line].split()[0])
        ETA.append(argosls[line].split()[1])
        print()
        une_base1.append([float(argosls[line].split()[0]),float(argosls[line].split()[1])])
        print("test base test base")
        line=line+1

    liste_base.append(une_base1)
print("testi")
print(liste_base)

print(liste_base[0])
#print(ICSU)
#line=line+ICSU[0]+ICSU[1]
#print(argosls[line].split())
#   8)     DO I = 1, NCONS
#            ICONU(I), LMNP1(I), NRCR(I)
#            DO J = 1, ICONU(I)
#              ZET(J,I), (ETA(K,J,I), K = 1, NRCR(I))
#            ENDDO
#          ENDDO
#
#   9)     IF(NCRS .NE. 0) THEN
#            DO ICRS = 1, NCRS
#              LCRU, LLSU
#              IF(LCRU .GE. 0) THEN
#                DO L = 0, LCRU
#                  NBFCR
#                  DO K = 1, NBFCR
#                    NCR(K), ZCR(K), CCR(K)
#                  ENDDO
#                ENDDO
#              ENDIF
#              IF(LLSU .GE. 1) THEN
#                DO L = 1, LLSU
#                  NBFCR
#                  DO K = 1, NBFCR
#                    NCR(K), ZCR(K), CCR(K)
#                  ENDDO
#                ENDDO
#              ENDIF
#            ENDDO
#          ENDIF
#
NS = int(NS)
#line=line+1
MTYPE = []
NF  = []
NC = []
CHG = np.zeros(NS)
pos=[]
print("testestestset")
print(type(argosls[line].split()))
MCONS = []
IGCS = []
molecule = molecule()
atoms = []
for i in range(NS):
    temp = argosls[line].split()
    print("test3")
    print(temp)
    MTYPE.append(temp[0])
    NF.append(int(temp[1]))
    print("testtestsetsettemp2")
    print(temp[2])
    NC.append(int(temp[2]))
    CHG[i] = temp[3]
    line = line+1
    print(NC[i])

    for j in range(NC[i]):

        postemp = np.array(argosls[line].split()).transpose()
        pos.append(postemp)
        atoms.append(atom(type=temp[0],chg=CHG[i],pos=postemp))
        #unatom =
        print('testi',i)
        molecule.add_atoms(atoms[i])
        line=line+1
        if (NC[i] != 1):
            print("TODO: NC[1] != 1")
        print("NF[",i,"]=",NF[i])
        for k in range(NF[i]):
            MCONS.append(argosls[line].split()[0])
            IGCS.append(argosls[line].split()[1])
            print("test base")
            print(argosls[line].split()[0], argosls[line].split()[1])
            print("Another Test Testing TestT")
            print(int(argosls[line].split()[0])-1)
            a = int(argosls[line].split()[0])-1
            print(liste_base[a])
            l = int(argosls[line].split()[1])
            basenumber = int(argosls[line].split()[0])-1
            print("l = ",'basenumber = ',basenumber)
            unebase = basis(liste_base[basenumber] ,l)
            print("Addinf  a base to atom["+str(i)+"]")
            molecule.atoms[i].add_basis(unebase)
            line = line +1

    #molecule.add_atoms(unatom)
        #            ENDDO
        #            IF(NCRS .NE. 0) THEN
        #              MCRS(IS)
        #            ENDIF
        #          ENDDO
print("testsetst")
print(pos[0])
#  10)     DO IS = 1, NS
#            MTYPE(IS), NF(IS), NC(IS), CHG(IS)
#            FORMAT(A3,2I3,F3.0)
#            DO J = 1, NC
#              X(J), Y(J), Z(J)
#            ENDDO
#            IF(NC(IS) .NE. 1) THEN
#              DO J = 1, NGEN
#                (IGEN(K,J), K = 1, NC)
#              ENDDO
#            ENDIF
#            DO J = 1, NF
#              MCONS(J), IGCS(J)
#            ENDDO
#            IF(NCRS .NE. 0) THEN
#              MCRS(IS)
#            ENDIF
#          ENDDO
#
#  11)     DO I = 2, 4
#            TITLE(I)
#            IF(TITLE(I) .EQ. ' ') EXIT
#          ENDDO




gridX=np.zeros((len(x),len(z)))
gridZ=np.zeros((len(x),len(z)))
r=0


x2 = np.zeros((len(x),len(z)))
z2 = np.zeros((len(x),len(z)))
# TODO : use this meshgrid and delete the loop
xv, zv = np.meshgrid(x, z, sparse=False, indexing='ij')

print("test molecule")
print(molecule.atoms)
n = 0
for atom in molecule.atoms:
    print(atom.type)
    print(atom.pos)
    print("Len atom.basis " + str(len(atom.basis)))
    for unebase in atom.basis:
        print('listecons = ',unebase.listecons)
        print('l = ',unebase.l)
        if unebase.l == 1:
            print("fonction 1s")
            orb = np.zeros(xv.shape)
            f = open('orb_' + str(n) + '.txt', 'w')
            for con in unebase.listecons:
                f.write(str(con[0])+" "+str(con[1])+"\n")
                orb = orb + con[1] * np.exp(-con[0] * (xv ** 2.0 + zv ** 2.0))
            f.close()
            plt.contourf(xv, zv, orb)
            plt.savefig('orb_'+str(n)+'.png')
            n = n + 1
#orbc1s_mod2 = 1.8347002E-03 * np.exp(-3047.525 * (x2 ** 2.0 + z2 ** 2.0)) \
#                  + 1.4037301E-02 * np.exp(-457.3695 * (x2 ** 2.0 + z2 ** 2.0)) \
#                  + 6.8842606E-02 * np.exp(-103.9487 * (x2 ** 2.0 + z2 ** 2.0)) \
#                  + 0.2321844 * np.exp(-29.21016 * (x2 ** 2.0 + z2 ** 2.0)) \
#                  + 0.4679413 * np.exp(-9.286663 * (x2 ** 2.0 + z2 ** 2.0)) \
#                  + 0.3623120 * np.exp(-3.163927 * (x2 ** 2.0 + z2 ** 2.0))

    #print(atom.basis)


for nx in range(len(x)):
    for nz in range(len(z)):
        x2[nx, nz] = x[nx]
        z2[nx, nz] = z[nz]
        r=np.sqrt(x[nx]**2.0+z[nz]**2.0)
        gridX[nx,nz]=x[nx]
        gridZ[nx,nz]=z[nz]
        #orbc1s[nx,nz]=1.8347002E-03*np.exp(-3047.525*r**2)+1.4037301E-02*np.exp(-457.3695*r**2)+6.8842606E-02*np.exp(-103.9487*r**2)+0.2321844*np.exp(-29.21016*r**2)+0.4679413*np.exp(-9.286663*r**2)+0.3623120*np.exp(-3.163927*r**2)
        #orbc1s_mod[nx, nz] = 1.8347002E-03 * np.exp(-3047.525 * (x[nx]**2.0+z[nz]**2.0)) \
                             #+ 1.4037301E-02 * np.exp(-457.3695 * (x[nx]**2.0+z[nz]**2.0)) \
                             #+ 6.8842606E-02 * np.exp(-103.9487 * (x[nx]**2.0+z[nz]**2.0)) \
                             #+ 0.2321844 * np.exp(-29.21016 * (x[nx]**2.0+z[nz]**2.0)) \
                             #+ 0.4679413 * np.exp(-9.286663 * (x[nx]**2.0+z[nz]**2.0)) \
                             #+ 0.3623120 * np.exp(-3.163927 * (x[nx]**2.0+z[nz]**2.0))

        orbc2s[nx,nz]=-0.1193324*np.exp(-7.868272*r**2)-0.1608542*np.exp(-1.881288*r**2)+1.143456*np.exp(-0.5442493*r**2)
        orbc3s[nx,nz]=np.exp(-0.1687144*r**2)
        orbc2px[nx,nz]=6.8999096E-02*np.exp(-7.868272*r**2)+0.3164240*np.exp(-1.881288*r**2)+0.7443083*np.exp(-0.5442493*r**2)*x[nx]
        orbc2pz[nx,nz]=6.8999096E-02*np.exp(-7.868272*r**2)+0.3164240*np.exp(-1.881288*r**2)+0.7443083*np.exp(-0.5442493*r**2)*z[nz]
        orbc3px[nx,nz]=np.exp(-0.1687144*r**2)*x[nx]
        orbc3pz[nx,nz]=np.exp(-0.1687144*r**2)*z[nz]
        orbc3s[nx,nz]=1.8310998E-03*np.exp(-5484.672*r**2)+1.3950099E-02*np.exp(-825.2350*r**2)+6.8445093E-0*np.exp(-188.0470*r**2)+0.2327143*np.exp(-52.96450*r**2)+0.4701930*np.exp(-16.89757*r**2)+0.3585209*np.exp(-5.799635*r**2)
        rO1=np.sqrt(x[nx]**2.0+(z[nz]+2.19208232)**2.0)
        orbO1_1s[nx,nz]=1.8310998E-03*np.exp(-5484.672*rO1**2)+1.3950099E-02*np.exp(-825.2350*rO1**2)+6.8445093E-02*np.exp(-188.0470*rO1**2)+0.2327143*np.exp(-52.96450*rO1**2)+0.4701930*np.exp(-16.89757*rO1**2)+0.3585209*np.exp(-5.799635*rO1**2)
        orbO1_2s[nx,nz]=-0.1107775*np.exp(-15.53962*rO1**2)-0.1480263*np.exp(-3.599934*rO1**2)+1.130767*np.exp(-1.013762*rO1**2)
        orbO1_3s[nx,nz]=np.exp(-0.2700058*rO1**2)
        orbO1_2px[nx,nz]=7.0874300E-02*np.exp(-15.53962*rO1**2)+0.3397528*np.exp(-3.599934*rO1**2)+0.7271586*np.exp(-1.013762*rO1**2)*x[nx]
        orbO1_2pz[nx,nz]=7.0874300E-02*np.exp(-15.53962*rO1**2)+0.3397528*np.exp(-3.599934*rO1**2)+0.7271586*np.exp(-1.013762*rO1**2)*(z[nz]+2.19208232)
        orbO1_3px[nx,nz]=np.exp(-0.2700058*rO1**2)*x[nx]
        orbO1_3pz[nx,nz]=np.exp(-0.2700058*rO1**2)*(z[nz]+2.19208232)

        rO2=np.sqrt(x[nx]**2.0+(z[nz]-2.19208232)**2.0)
        orbO2_1s[nx,nz]=1.8310998E-03*np.exp(-5484.672*rO2**2)+1.3950099E-02*np.exp(-825.2350*rO2**2)+6.8445093E-02*np.exp(-188.0470*rO2**2)+0.2327143*np.exp(-52.96450*rO2**2)+0.4701930*np.exp(-16.89757*rO2**2)+0.3585209*np.exp(-5.799635*rO2**2)
        #*np.exp(-*rO2**2)
        orbO2_2s[nx,nz]=-0.1107775*np.exp(-15.53962*rO2**2)-0.1480263*np.exp(-3.599934*rO2**2)+1.130767*np.exp(-1.013762*rO2**2)
        orbO2_3s[nx,nz]=np.exp(-0.2700058*rO2**2)
        #*np.exp(-*rO1**2)
        orbO2_2px[nx,nz]=7.0874300E-02*np.exp(-15.53962*rO2**2)+0.3397528*np.exp(-3.599934*rO2**2)+0.7271586*np.exp(-1.013762*rO2**2)*x[nx]
        orbO2_2pz[nx,nz]=7.0874300E-02*np.exp(-15.53962*rO2**2)+0.3397528*np.exp(-3.599934*rO2**2)+0.7271586*np.exp(-1.013762*rO2**2)*(z[nz]-2.19208232)
        orbO2_3px[nx,nz]=np.exp(-0.2700058*rO2**2)*x[nx]
        orbO2_3pz[nx,nz]=np.exp(-0.2700058*rO2**2)*(z[nz]-2.19208232)
        MO1[nx,nz]=mocoeff[0][0]*orbc1s[nx,nz]+mocoeff[0][1]*orbc2s[nx,nz]+mocoeff[0][2]*orbc3s[nx,nz]+mocoeff[0][3]*orbc2px[nx,nz]+mocoeff[0][4]*orbc2py[nx,nz]+mocoeff[0][5]*orbc2pz[nx,nz]+mocoeff[0][6]*orbc3px[nx,nz]+\
        mocoeff[0][7]*orbc3py[nx,nz]+\
        mocoeff[0][8]*orbc3pz[nx,nz]+\
        mocoeff[0][9]*orbO1_1s[nx,nz]+\
        mocoeff[0][10]*orbO1_2s[nx,nz]+\
        mocoeff[0][11]*orbO1_3s[nx,nz]+\
        mocoeff[0][12]*orbO1_2px[nx,nz]+\
        mocoeff[0][13]*orbO1_2py[nx,nz]+\
        mocoeff[0][14]*orbO1_2pz[nx,nz]+\
        mocoeff[0][15]*orbO1_3px[nx,nz]+\
        mocoeff[0][16]*orbO1_3py[nx,nz]+\
        mocoeff[0][17]*orbO1_3pz[nx,nz]+\
        mocoeff[0][18]*orbO2_1s[nx,nz]+\
        mocoeff[0][19]*orbO2_2s[nx,nz]+\
        mocoeff[0][20]*orbO2_3s[nx,nz]+\
        mocoeff[0][21]*orbO2_2px[nx,nz]+\
        mocoeff[0][22]*orbO2_2py[nx,nz]+\
        mocoeff[0][23]*orbO2_2pz[nx,nz]+\
        mocoeff[0][24]*orbO2_3px[nx,nz]+\
        mocoeff[0][25]*orbO2_3py[nx,nz]+\
        mocoeff[0][26]*orbO2_3pz[nx,nz]


        for n in range(27):
          MO[n,nx,nz]=mocoeff[n][0]*orbc1s[nx,nz]+mocoeff[n][1]*orbc2s[nx,nz]+mocoeff[n][2]*orbc3s[nx,nz]+mocoeff[n][3]*orbc2px[nx,nz]+mocoeff[n][4]*orbc2py[nx,nz]+mocoeff[n][5]*orbc2pz[nx,nz]+mocoeff[n][6]*orbc3px[nx,nz]+\
          mocoeff[n][7]*orbc3py[nx,nz]+\
          mocoeff[n][8]*orbc3pz[nx,nz]+\
          mocoeff[n][9]*orbO1_1s[nx,nz]+\
          mocoeff[n][10]*orbO1_2s[nx,nz]+\
          mocoeff[n][11]*orbO1_3s[nx,nz]+\
          mocoeff[n][12]*orbO1_2px[nx,nz]+\
          mocoeff[n][13]*orbO1_2py[nx,nz]+\
          mocoeff[n][14]*orbO1_2pz[nx,nz]+\
          mocoeff[n][15]*orbO1_3px[nx,nz]+\
          mocoeff[n][16]*orbO1_3py[nx,nz]+\
          mocoeff[n][17]*orbO1_3pz[nx,nz]+\
          mocoeff[n][18]*orbO2_1s[nx,nz]+\
          mocoeff[n][19]*orbO2_2s[nx,nz]+\
          mocoeff[n][20]*orbO2_3s[nx,nz]+\
          mocoeff[n][21]*orbO2_2px[nx,nz]+\
          mocoeff[n][22]*orbO2_2py[nx,nz]+\
          mocoeff[n][23]*orbO2_2pz[nx,nz]+\
          mocoeff[n][24]*orbO2_3px[nx,nz]+\
          mocoeff[n][25]*orbO2_3py[nx,nz]+\
          mocoeff[n][26]*orbO2_3pz[nx,nz]

orbc1s_mod2 = 1.8347002E-03 * np.exp(-3047.525 * (x2**2.0+z2**2.0)) \
                             + 1.4037301E-02 * np.exp(-457.3695 * (x2**2.0+z2**2.0)) \
                             + 6.8842606E-02 * np.exp(-103.9487 * (x2**2.0+z2**2.0)) \
                             + 0.2321844 * np.exp(-29.21016 * (x2**2.0+z2**2.0)) \
                             + 0.4679413 * np.exp(-9.286663 * (x2**2.0+z2**2.0)) \
                             + 0.3623120 * np.exp(-3.163927 * (x2**2.0+z2**2.0))

#*np.exp**(*r**2)
plt.contourf(gridX,gridZ ,orbc1s)
plt.savefig('orbc1s.png')
plt.contourf(gridX,gridZ ,orbc1s_mod)
plt.savefig('orbc1s_mod.png')
plt.contourf(gridX,gridZ ,orbc1s_mod2)
plt.savefig('orbc1s_mod2.png')
plt.contourf(gridX,gridZ ,orbc2s)
plt.savefig('orbc2s.png')
plt.contourf(gridX,gridZ ,orbc3s)
plt.savefig('orbc3s.png')
plt.contourf(gridX,gridZ ,orbc2px)
plt.savefig('orbc2px.png')
plt.contourf(gridX,gridZ ,orbc2pz)
plt.savefig('orbc2pz.png')
plt.contourf(gridX,gridZ ,orbc3px)
plt.savefig('orbc3px.png')
plt.contourf(gridX,gridZ ,orbc3pz)
plt.savefig('orbc3pz.png')
plt.contourf(gridX,gridZ ,orbO1_1s)
plt.savefig('orbO1_1s.png')
plt.contourf(gridX,gridZ ,orbO1_2s)
plt.savefig('orbO1_2s.png')
plt.contourf(gridX,gridZ ,orbO1_3s)
plt.savefig('orbO1_3s.png')
plt.contourf(gridX,gridZ ,orbO1_2px)
plt.savefig('orbO1_2px.png')
plt.contourf(gridX,gridZ ,orbO1_2pz)
plt.savefig('orbO1_2pz.png')
plt.contourf(gridX,gridZ ,orbO1_3px)
plt.savefig('orbO1_3px.png')
plt.contourf(gridX,gridZ ,orbO1_3pz)
plt.savefig('orbO1_3pz.png')
plt.contourf(gridX,gridZ ,orbO2_1s)
plt.savefig('orbO2_1s.png')
plt.contourf(gridX,gridZ ,orbO2_2s)
plt.savefig('orbO2_2s.png')
plt.contourf(gridX,gridZ ,orbO2_3s)
plt.savefig('orbO2_3s.png')
plt.contourf(gridX,gridZ ,orbO2_2px)
plt.savefig('orbO2_2px.png')
plt.contourf(gridX,gridZ ,orbO2_2pz)
plt.savefig('orbO2_2pz.png')
plt.contourf(gridX,gridZ ,orbO2_3px)
plt.savefig('orbO2_3px.png')
plt.contourf(gridX,gridZ ,orbO2_3pz)
plt.savefig('orbO2_3pz.png')
plt.contourf(gridX,gridZ ,MO1)
plt.savefig('orbMO1.png')
plt.contourf(gridX,gridZ ,MO2)
plt.savefig('orbMO2.png')
for n in range(27):
    plt.contourf(gridX,gridZ ,MO[n,:,:])
    plt.savefig('orbMO'+str(n)+'.png')

