echo of the argos input file:
 ------------------------------------------------------------------------
 BH3
    0   4   2  13   2  20   9   0   0   0   0   0   0   0   0
   1  1  a
   0
     1    1
     3    1    1    1
     1    1    1
     1
     3    3    2
     0    0    1
     1    0    0
     0    1    0
     3    1    1
     116.4340000        0.0629605
      17.4314000        0.3633040
       3.6801600        0.6972550
     2    1    1
       2.2818700       -0.3686620
       0.4652480        1.1994400
     1    1    1
       0.1243280        1.0000000
     1    1    1
       0.0315000        1.0000000
     2    2    1
       2.2818700        0.2311520
       0.4652480        0.8667640
     1    2    1
       0.1243280        1.0000000
     1    2    1
       0.0315000        1.0000000
     2    1    1
       5.4471780        0.1562850
       0.8245470        0.9046910
     1    1    1
       0.1831920        1.0000000
     2    1    1
       5.4471780        0.1562850
       0.8245470        0.9046910
     1    1    1
       0.1831920        1.0000000
     2    1    1
       5.4471780        0.1562850
       0.8245470        0.9046910
     1    1    1
       0.1831920        1.0000000
 B    7  1 5.
     0.00000000    0.00000000    0.00000000
    1   1
    2   1
    3   1
    4   1
    5   2
    6   2
    7   2
 H    2  1 1.
     0.00000000    2.24106969    0.00000000
    8   1
    9   1
 H    2  1 1.
     0.00000000   -1.12053390    1.94082244
   10   1
   11   1
 H    2  1 1.
     0.00000000   -1.12053390   -1.94082244
   12   1
   13   1
 ------------------------------------------------------------------------
                              program "argos" 5.9
                            columbus program system

             this program computes integrals over symmetry orbitals
               of generally contracted gaussian atomic orbitals.
                        programmed by russell m. pitzer

                           version date: 20-aug-2001

references:

    symmetry analysis (equal contributions):
    r. m. pitzer, j. chem. phys. 58, 3111 (1973).

    ao integral evaluation (hondo):
    m. dupuis, j. rys, and h. f. king, j. chem. phys. 65, 111 (1976).

    general contraction of gaussian orbitals:
    r. c. raffenetti, j. chem. phys. 58, 4452 (1973).

    core potential ao integrals (meldps):
    l. e. mcmurchie and e. r. davidson, j. comput. phys. 44, 289 (1981)

    spin-orbit and core potential integrals:
    r. m. pitzer and n. w. winter, int. j. quantum chem. 40, 773 (1991)

 This Version of Program ARGOS is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de

*********************** File revision status: ***********************
* argos1.f  Revision: 2.15.12.1         Date: 2013/04/11 14:37:29   * 
* argos2.f  Revision: 2.4.12.1          Date: 2013/04/11 14:37:29   * 
* argos3.f  Revision: 2.3.12.1           Date: 2013/04/11 14:37:29  * 
* argos4.f  Revision: 2.4.6.1           Date: 2013/04/11 14:37:29   * 
* argos5.f  Revision: 2.2.20.1          Date: 2013/04/11 14:37:29   * 
* argos6.f  Revision: 2.22.6.1          Date: 2013/04/11 14:37:29   * 
********************************************************************

 workspace allocation parameters: lcore= 131072000 mem1=         0 ifirst=         1

 filenames and unit numbers:
 unit description                   filename
 ---- -----------                   ----------
   6  listing file:                 argosls                                                     
   4  1-e integral file:            aoints                                                      
   8  2-e integral file [fsplit=2]: aoints2                                                     
   5  input file:                   argosin                                                     

 argos input parameters and titles:
 ngen   =  0 ns     =  4 naords =  2 ncons  = 13 ngcs   =  2 itol   = 20
 icut   =  9 aoints =  4 only1e =  0 inrm   =  0 ncrs   =  0
 l1rec  =         0      l2rec  =         0      aoint2 =  8 fsplit =  2

BH3                                                                             
aoints SIFS file created by argos.      colosse2          16:37:21.494 14-Jan-18


irrep            1
degeneracy       1
label             a


direct product table
   (  a) (  a) =   a


                     nuclear repulsion energy    7.46610285


primitive ao integrals neglected if exponential factor below 10**(-20)
contracted ao and so integrals neglected if value below 10**(- 9)
symmetry orbital integrals written on units  4  8


                                    B   atoms

                              nuclear charge   5.00

           center            x               y               z
             1             0.00000000      0.00000000      0.00000000

                    1s orbitals

 orbital exponents  contraction coefficients
    116.4340       6.2960466E-02
    17.43140       0.3633038    
    3.680160       0.6972546    

                     symmetry orbital labels
                     1  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
    2.281870      -0.3686635    
   0.4652480        1.199445    

                     symmetry orbital labels
                     2  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.1243280        1.000000    

                     symmetry orbital labels
                     3  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   3.1500000E-02    1.000000    

                     symmetry orbital labels
                     4  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
    2.281870       0.2311519    
   0.4652480       0.8667636    

                     symmetry orbital labels
                     5  a1
                     6  a1
                     7  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  0.000  1.000  0.000
  1, 010  0.000  0.000  1.000
  1, 001  1.000  0.000  0.000

                    2p orbitals

 orbital exponents  contraction coefficients
   0.1243280        1.000000    

                     symmetry orbital labels
                     8  a1
                     9  a1
                    10  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  0.000  1.000  0.000
  1, 010  0.000  0.000  1.000
  1, 001  1.000  0.000  0.000

                    2p orbitals

 orbital exponents  contraction coefficients
   3.1500000E-02    1.000000    

                     symmetry orbital labels
                    11  a1
                    12  a1
                    13  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  0.000  1.000  0.000
  1, 010  0.000  0.000  1.000
  1, 001  1.000  0.000  0.000


                                    H   atoms

                              nuclear charge   1.00

           center            x               y               z
             1             0.00000000      2.24106969      0.00000000

                    1s orbitals

 orbital exponents  contraction coefficients
    5.447178       0.1562850    
   0.8245470       0.9046909    

                     symmetry orbital labels
                    14  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.1831920        1.000000    

                     symmetry orbital labels
                    15  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000


                                    H   atoms

                              nuclear charge   1.00

           center            x               y               z
             1             0.00000000     -1.12053390      1.94082244

                    1s orbitals

 orbital exponents  contraction coefficients
    5.447178       0.1562850    
   0.8245470       0.9046909    

                     symmetry orbital labels
                    16  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.1831920        1.000000    

                     symmetry orbital labels
                    17  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000


                                    H   atoms

                              nuclear charge   1.00

           center            x               y               z
             1             0.00000000     -1.12053390     -1.94082244

                    1s orbitals

 orbital exponents  contraction coefficients
    5.447178       0.1562850    
   0.8245470       0.9046909    

                     symmetry orbital labels
                    18  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.1831920        1.000000    

                     symmetry orbital labels
                    19  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

lx:   a          ly:   a          lz:   a

output SIFS file header information:
BH3                                                                             
aoints SIFS file created by argos.      colosse2          16:37:21.494 14-Jan-18

output energy(*) values:
 energy( 1)=  7.466102847527E+00, ietype=   -1,    core energy of type: Nuc.Rep.

total core energy =  7.466102847527E+00

nsym = 1 nbft=  19

symmetry  =    1
slabel(*) =    a
nbpsy(*)  =   19

info(*) =         2      4096      3272      4096      2700         0

output orbital labels, i:bfnlab(i)=
   1:  1B__1s   2:  2B__1s   3:  3B__1s   4:  4B__1s   5:  5B__2p   6:  6B__2p
   7:  7B__2p   8:  8B__2p   9:  9B__2p  10: 10B__2p  11: 11B__2p  12: 12B__2p
  13: 13B__2p  14: 14H__1s  15: 15H__1s  16: 16H__1s  17: 17H__1s  18: 18H__1s
  19: 19H__1s

bfn_to_center map(*), i:map(i)
   1:  1   2:  1   3:  1   4:  1   5:  1   6:  1   7:  1   8:  1   9:  1  10:  1
  11:  1  12:  1  13:  1  14:  2  15:  2  16:  3  17:  3  18:  4  19:  4

bfn_to_orbital_type map(*), i:map(i)
   1:  1   2:  1   3:  1   4:  1   5:  2   6:  2   7:  2   8:  2   9:  2  10:  2
  11:  2  12:  2  13:  2  14:  1  15:  1  16:  1  17:  1  18:  1  19:  1


       19 symmetry orbitals,        a:  19

 socfpd: mcxu=      174 mcxu2=      137 left=131071826
 
oneint:   103 S1(*)    integrals were written in  1 records.
oneint:   103 T1(*)    integrals were written in  1 records.
oneint:   114 V1(*)    integrals were written in  1 records.
oneint:    30 X(*)     integrals were written in  1 records.
oneint:    87 Y(*)     integrals were written in  1 records.
oneint:    74 Z(*)     integrals were written in  1 records.
oneint:    30 Im(px)   integrals were written in  1 records.
oneint:    74 Im(py)   integrals were written in  1 records.
oneint:    70 Im(pz)   integrals were written in  1 records.
oneint:    51 Im(lx)   integrals were written in  1 records.
oneint:    21 Im(ly)   integrals were written in  1 records.
oneint:    27 Im(lz)   integrals were written in  1 records.
oneint:   103 XX(*)    integrals were written in  1 records.
oneint:    27 XY(*)    integrals were written in  1 records.
oneint:    21 XZ(*)    integrals were written in  1 records.
oneint:   103 YY(*)    integrals were written in  1 records.
oneint:    71 YZ(*)    integrals were written in  1 records.
oneint:   103 ZZ(*)    integrals were written in  1 records.
 

twoint:        8941 1/r12    integrals and     1720 pk flags
                                 were written in     4 records.

 twoint: maximum mblu needed =       657
 
driver: 1-e  integral workspace high-water mark =     21176
driver: 2-e  integral workspace high-water mark =     21381
driver: overall argos workspace high-water mark =     21381
