program testGS
use variables
use basics
use string_module,only:string
use read_input
use grille
use qp_integrals
use read_columbus
use nouveaudrt
use observable
!use Transformation
use Tests
use dynamique
!use matrix_temp_highlvl
use matrix_rearrangement_module

 Complex(kind=comp_16), Allocatable,Dimension(:,:) :: Psi
 Integer(kind=int_4), Allocatable,Dimension(:)     :: Ni
 Real(kind=real_8),allocatable,dimension(:)        :: X0
 Complex(kind=comp_16), Allocatable,Dimension(:,:) :: coeffGS

allocate(Ni(2))
Ni(1)=2
Ni(2)=1000
allocate(Psi(Ni(1),Ni(2)))
allocate(X0(Ni(1)))
X0(1)=400
X0(2)=600
do i=1,Ni(1)
  do j=1,Ni(2)
    psi(i,j)=exp((-(j-X0(i))**(2.d0)/(125.d0**2.d0)))
  end do
end do
  do j=1,Ni(2)
    write(99,*) j,(abs(psi(i,j))**2.d0,i=1,2)
  end do
allocate(coeffGS(Ni(1),Ni(1))
allocate(psipsi(2,2))
psipsi(1,1)=(1.d0,0.d0)
psipsi(1,2)=(0.25d0,0.d0)
psipsi(2,1)=(0.25d0,0.d0)
psipsi(2,2)=(1.d0,0.d0)
coeffGS=0.d0

end program testGS
