echo of the argos input file:
 ------------------------------------------------------------------------
 test
    0   2   2  10   2  20   9   0   0   0   0   0   0   0   0
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
     6    1    1
    4173.5110000        0.0018348
     627.4579000        0.0139950
     142.9021000        0.0685870
      40.2343300        0.2322410
      12.8202100        0.4690700
       4.3904370        0.3604550
     3    1    1
      11.6263580       -0.1149610
       2.7162800       -0.1691180
       0.7722180        1.1458520
     1    1    1
       0.2120313        1.0000000
     3    2    1
      11.6263580        0.0675800
       2.7162800        0.3239070
       0.7722180        0.7408950
     1    2    1
       0.2120313        1.0000000
     6    1    1
    4173.5110000        0.0018348
     627.4579000        0.0139950
     142.9021000        0.0685870
      40.2343300        0.2322410
      12.8202100        0.4690700
       4.3904370        0.3604550
     3    1    1
      11.6263580       -0.1149610
       2.7162800       -0.1691180
       0.7722180        1.1458520
     1    1    1
       0.2120313        1.0000000
     3    2    1
      11.6263580        0.0675800
       2.7162800        0.3239070
       0.7722180        0.7408950
     1    2    1
       0.2120313        1.0000000
 N    5  1 7.
     0.00000000    0.00000000   -0.94500000
    1   1
    2   1
    3   1
    4   2
    5   2
 N    5  1 7.
     0.00000000    0.00000000    0.94500000
    6   1
    7   1
    8   1
    9   2
   10   2
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
 ngen   =  0 ns     =  2 naords =  2 ncons  = 10 ngcs   =  2 itol   = 20
 icut   =  9 aoints =  4 only1e =  0 inrm   =  0 ncrs   =  0
 l1rec  =         0      l2rec  =         0      aoint2 =  8 fsplit =  2

test                                                                         
aoints SIFS file created by argos.      localhost.localdo 11:40:47.526 13-Jan-21


irrep            1
degeneracy       1
label             a


direct product table
   (  a) (  a) =   a


                     nuclear repulsion energy   25.92592593


primitive ao integrals neglected if exponential factor below 10**(-20)
contracted ao and so integrals neglected if value below 10**(- 9)
symmetry orbital integrals written on units  4  8


                                    N   atoms

                              nuclear charge   7.00

           center            x               y               z
             1             0.00000000      0.00000000     -0.94500000

                    1s orbitals

 orbital exponents  contraction coefficients
    4173.511       1.8347994E-03
    627.4579       1.3994996E-02
    142.9021       6.8586979E-02
    40.23433       0.2322409    
    12.82021       0.4690699    
    4.390437       0.3604549    

                     symmetry orbital labels
                     1  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
    11.62636      -0.1149610    
    2.716280      -0.1691180    
   0.7722180        1.145852    

                     symmetry orbital labels
                     2  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.2120313        1.000000    

                     symmetry orbital labels
                     3  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
    11.62636       6.7580023E-02
    2.716280       0.3239071    
   0.7722180       0.7408953    

                     symmetry orbital labels
                     4  a1
                     5  a1
                     6  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  0.000  1.000  0.000
  1, 010  0.000  0.000  1.000
  1, 001  1.000  0.000  0.000

                    2p orbitals

 orbital exponents  contraction coefficients
   0.2120313        1.000000    

                     symmetry orbital labels
                     7  a1
                     8  a1
                     9  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  0.000  1.000  0.000
  1, 010  0.000  0.000  1.000
  1, 001  1.000  0.000  0.000


                                    N   atoms

                              nuclear charge   7.00

           center            x               y               z
             1             0.00000000      0.00000000      0.94500000

                    1s orbitals

 orbital exponents  contraction coefficients
    4173.511       1.8347994E-03
    627.4579       1.3994996E-02
    142.9021       6.8586979E-02
    40.23433       0.2322409    
    12.82021       0.4690699    
    4.390437       0.3604549    

                     symmetry orbital labels
                    10  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
    11.62636      -0.1149610    
    2.716280      -0.1691180    
   0.7722180        1.145852    

                     symmetry orbital labels
                    11  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.2120313        1.000000    

                     symmetry orbital labels
                    12  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
    11.62636       6.7580023E-02
    2.716280       0.3239071    
   0.7722180       0.7408953    

                     symmetry orbital labels
                    13  a1
                    14  a1
                    15  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  0.000  1.000  0.000
  1, 010  0.000  0.000  1.000
  1, 001  1.000  0.000  0.000

                    2p orbitals

 orbital exponents  contraction coefficients
   0.2120313        1.000000    

                     symmetry orbital labels
                    16  a1
                    17  a1
                    18  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  0.000  1.000  0.000
  1, 010  0.000  0.000  1.000
  1, 001  1.000  0.000  0.000

lx:   a          ly:   a          lz:   a

output SIFS file header information:
test                                                                         
aoints SIFS file created by argos.      localhost.localdo 11:40:47.526 13-Jan-21

output energy(*) values:
 energy( 1)=  2.592592592593E+01, ietype=   -1,    core energy of type: Nuc.Rep.

total core energy =  2.592592592593E+01

nsym = 1 nbft=  18

symmetry  =    1
slabel(*) =    a
nbpsy(*)  =   18

info(*) =         2      4096      3272      4096      2700         0

output orbital labels, i:bfnlab(i)=
   1:  1N__1s   2:  2N__1s   3:  3N__1s   4:  4N__2p   5:  5N__2p   6:  6N__2p
   7:  7N__2p   8:  8N__2p   9:  9N__2p  10: 10N__1s  11: 11N__1s  12: 12N__1s
  13: 13N__2p  14: 14N__2p  15: 15N__2p  16: 16N__2p  17: 17N__2p  18: 18N__2p

bfn_to_center map(*), i:map(i)
   1:  1   2:  1   3:  1   4:  1   5:  1   6:  1   7:  1   8:  1   9:  1  10:  2
  11:  2  12:  2  13:  2  14:  2  15:  2  16:  2  17:  2  18:  2

bfn_to_orbital_type map(*), i:map(i)
   1:  1   2:  1   3:  1   4:  2   5:  2   6:  2   7:  2   8:  2   9:  2  10:  1
  11:  1  12:  1  13:  2  14:  2  15:  2  16:  2  17:  2  18:  2


       18 symmetry orbitals,        a:  18

 socfpd: mcxu=      174 mcxu2=      137 left=131071826
 
oneint:    63 S1(*)    integrals were written in  1 records.
oneint:    63 T1(*)    integrals were written in  1 records.
oneint:    75 V1(*)    integrals were written in  1 records.
oneint:    32 X(*)     integrals were written in  1 records.
oneint:    32 Y(*)     integrals were written in  1 records.
oneint:    66 Z(*)     integrals were written in  1 records.
oneint:    32 Im(px)   integrals were written in  1 records.
oneint:    32 Im(py)   integrals were written in  1 records.
oneint:    45 Im(pz)   integrals were written in  1 records.
oneint:    40 Im(lx)   integrals were written in  1 records.
oneint:    40 Im(ly)   integrals were written in  1 records.
oneint:    16 Im(lz)   integrals were written in  1 records.
oneint:    63 XX(*)    integrals were written in  1 records.
oneint:    16 XY(*)    integrals were written in  1 records.
oneint:    38 XZ(*)    integrals were written in  1 records.
oneint:    63 YY(*)    integrals were written in  1 records.
oneint:    38 YZ(*)    integrals were written in  1 records.
oneint:    75 ZZ(*)    integrals were written in  1 records.
 

twoint:        4350 1/r12    integrals and      651 pk flags
                                 were written in     2 records.

 twoint: maximum mblu needed =       690
 
driver: 1-e  integral workspace high-water mark =     21157
driver: 2-e  integral workspace high-water mark =     21395
driver: overall argos workspace high-water mark =     21395
