echo of the argos input file:
 ------------------------------------------------------------------------
 O2
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
    5484.6717000        0.0018311
     825.2349500        0.0139501
     188.0469600        0.0684451
      52.9645000        0.2327143
      16.8975700        0.4701930
       5.7996353        0.3585209
     3    1    1
      15.5396160       -0.1107775
       3.5999336       -0.1480263
       1.0137618        1.1307670
     1    1    1
       0.2700058        1.0000000
     3    2    1
      15.5396160        0.0708743
       3.5999336        0.3397528
       1.0137618        0.7271586
     1    2    1
       0.2700058        1.0000000
     6    1    1
    5484.6717000        0.0018311
     825.2349500        0.0139501
     188.0469600        0.0684451
      52.9645000        0.2327143
      16.8975700        0.4701930
       5.7996353        0.3585209
     3    1    1
      15.5396160       -0.1107775
       3.5999336       -0.1480263
       1.0137618        1.1307670
     1    1    1
       0.2700058        1.0000000
     3    2    1
      15.5396160        0.0708743
       3.5999336        0.3397528
       1.0137618        0.7271586
     1    2    1
       0.2700058        1.0000000
 O    5  1 8.
     0.00000000    0.00000000   -1.10100000
    1   1
    2   1
    3   1
    4   2
    5   2
 O    5  1 8.
     0.00000000    0.00000000    1.10100000
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

O2                                                                              
aoints SIFS file created by argos.      localhost.localdo 13:12:07.714 13-Jan-21


irrep            1
degeneracy       1
label             a


direct product table
   (  a) (  a) =   a


                     nuclear repulsion energy   29.06448683


primitive ao integrals neglected if exponential factor below 10**(-20)
contracted ao and so integrals neglected if value below 10**(- 9)
symmetry orbital integrals written on units  4  8


                                    O   atoms

                              nuclear charge   8.00

           center            x               y               z
             1             0.00000000      0.00000000     -1.10100000

                    1s orbitals

 orbital exponents  contraction coefficients
    5484.672       1.8310998E-03
    825.2350       1.3950099E-02
    188.0470       6.8445093E-02
    52.96450       0.2327143    
    16.89757       0.4701930    
    5.799635       0.3585209    

                     symmetry orbital labels
                     1  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
    15.53962      -0.1107775    
    3.599934      -0.1480263    
    1.013762        1.130767    

                     symmetry orbital labels
                     2  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.2700058        1.000000    

                     symmetry orbital labels
                     3  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
    15.53962       7.0874300E-02
    3.599934       0.3397528    
    1.013762       0.7271586    

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
   0.2700058        1.000000    

                     symmetry orbital labels
                     7  a1
                     8  a1
                     9  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  0.000  1.000  0.000
  1, 010  0.000  0.000  1.000
  1, 001  1.000  0.000  0.000


                                    O   atoms

                              nuclear charge   8.00

           center            x               y               z
             1             0.00000000      0.00000000      1.10100000

                    1s orbitals

 orbital exponents  contraction coefficients
    5484.672       1.8310998E-03
    825.2350       1.3950099E-02
    188.0470       6.8445093E-02
    52.96450       0.2327143    
    16.89757       0.4701930    
    5.799635       0.3585209    

                     symmetry orbital labels
                    10  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
    15.53962      -0.1107775    
    3.599934      -0.1480263    
    1.013762        1.130767    

                     symmetry orbital labels
                    11  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.2700058        1.000000    

                     symmetry orbital labels
                    12  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
    15.53962       7.0874300E-02
    3.599934       0.3397528    
    1.013762       0.7271586    

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
   0.2700058        1.000000    

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
O2                                                                              
aoints SIFS file created by argos.      localhost.localdo 13:12:07.714 13-Jan-21

output energy(*) values:
 energy( 1)=  2.906448683015E+01, ietype=   -1,    core energy of type: Nuc.Rep.

total core energy =  2.906448683015E+01

nsym = 1 nbft=  18

symmetry  =    1
slabel(*) =    a
nbpsy(*)  =   18

info(*) =         2      4096      3272      4096      2700         0

output orbital labels, i:bfnlab(i)=
   1:  1O__1s   2:  2O__1s   3:  3O__1s   4:  4O__2p   5:  5O__2p   6:  6O__2p
   7:  7O__2p   8:  8O__2p   9:  9O__2p  10: 10O__1s  11: 11O__1s  12: 12O__1s
  13: 13O__2p  14: 14O__2p  15: 15O__2p  16: 16O__2p  17: 17O__2p  18: 18O__2p

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
 

twoint:        4345 1/r12    integrals and      651 pk flags
                                 were written in     2 records.

 twoint: maximum mblu needed =       690
 
driver: 1-e  integral workspace high-water mark =     21157
driver: 2-e  integral workspace high-water mark =     21395
driver: overall argos workspace high-water mark =     21395
