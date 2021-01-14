program test
double precision :: wnm,cmm1,ua,E0wcm2,E0ua
wnm=2100
cmm1=1.d7/wnm
ua=cmm1/2.19475d5
E0wcm2=1.d14
E0ua=(E0wcm2/3.51d16)**(0.5d0)
write(*,*) "w(nm)=", wnm, "w(cm-1) = " , cmm1, "w(u.a.) = ", ua
write(*,*) "E(w/cm2)=", E0wcm2, "E0(u.a.) = " , E0ua

end program test

