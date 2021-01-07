Module precision_module

Implicit none
Public
Save

Integer,parameter :: int_4=selected_int_kind(9) !9

Integer,parameter :: int_8=selected_int_kind(18)

Integer,parameter :: real_4=selected_real_kind(p=6)

Integer,parameter :: real_8=selected_real_kind(p=15) !15

Integer,parameter :: comp_16=selected_real_kind(p=15) !15

Real(kind=real_8),parameter :: zero=0.0_real_8

Real(kind=real_8),parameter :: quarter=0.25_real_8

Real(kind=real_8),parameter :: half=0.5_real_8

Real(kind=real_8),parameter :: one=1.0_real_8

Real(kind=real_8),parameter :: smallest=tiny(one)

Real(kind=real_8),parameter :: largest=huge(one)

Real(kind=real_8),parameter :: onehalf=1.5_real_8

Real(kind=real_8),parameter :: two=2.0_real_8

Real(kind=real_8),parameter :: three=3.0_real_8

Real(kind=real_8),parameter :: four=4.0_real_8

Real(kind=real_8),parameter :: five=5.0_real_8

Real(kind=real_8),parameter :: six=6.0_real_8

Real(kind=real_8),parameter :: seven=7.0_real_8

Real(kind=real_8),parameter :: eight=8.0_real_8

Real(kind=real_8),parameter :: nine=9.0_real_8

Real(kind=real_8),parameter :: ten=10.0_real_8

Real(kind=real_8),parameter :: eleven=11.0_real_8

Real(kind=real_8),parameter :: twelve=12.0_real_8

Real(kind=real_8),parameter :: thirteen=13.0_real_8

Real(kind=real_8),parameter :: fourteen=14.0_real_8

Real(kind=real_8),parameter :: fifteen=15.0_real_8

Real(kind=real_8),parameter :: sixteen=16.0_real_8

Real(kind=real_8),parameter :: seventeen=17.0_real_8

Real(kind=real_8),parameter :: eighteen=18.0_real_8

Real(kind=real_8),parameter :: nineteen=19.0_real_8

Real(kind=real_8),parameter :: twenty=20.0_real_8

Real(kind=real_8),parameter :: twentyone=21.0_real_8

Real(kind=real_8),parameter :: twentytwo=22.0_real_8

Real(kind=real_8),parameter :: twentythree=23.0_real_8

Real(kind=real_8),parameter :: twentyfour=24.0_real_8

Real(kind=real_8),parameter :: twentyfive=25.0_real_8

Real(kind=real_8),parameter :: twentysix=26.0_real_8

Real(kind=real_8),parameter :: twentyseven=27.0_real_8

Real(kind=real_8),parameter :: twentyeight=28.0_real_8

Real(kind=real_8),parameter :: twentynine=29.0_real_8

Real(kind=real_8),parameter :: thirty=30.0_real_8

Real(kind=real_8),parameter :: thirtyone=31.0_real_8

Real(kind=real_8),parameter :: thirtytwo=32.0_real_8

Real(kind=real_8),parameter :: thirtythree=33.0_real_8

Real(kind=real_8),parameter :: thirtyfour=34.0_real_8

Real(kind=real_8),parameter :: thirtyfive=35.0_real_8

Real(kind=real_8),parameter :: thirtysix=36.0_real_8

Real(kind=real_8),parameter :: thirtyseven=37.0_real_8

Real(kind=real_8),parameter :: thirtyeight=38.0_real_8

Real(kind=real_8),parameter :: thirtynine=39.0_real_8

Real(kind=real_8),parameter :: forty=40.0_real_8

Real(kind=real_8),parameter :: fifty=50.0_real_8

Real(kind=real_8),parameter :: sixty=60.0_real_8

Real(kind=real_8),parameter :: seventy=70.0_real_8

Real(kind=real_8),parameter :: eighty=80.0_real_8

Real(kind=real_8),parameter :: ninety=90.0_real_8

Real(kind=real_8),parameter :: hundred=100.0_real_8

Real(kind=real_8),parameter :: twohundred=200.0_real_8

Real(kind=real_8),parameter :: threehundred=300.0_real_8

Real(kind=real_8),parameter :: fourhundred=400.0_real_8

Real(kind=real_8),parameter :: fivehundred=500.0_real_8

Real(kind=real_8),parameter :: sixhundred=600.0_real_8

Real(kind=real_8),parameter :: sevenhundred=700.0_real_8

Real(kind=real_8),parameter :: eighthundred=800.0_real_8

Real(kind=real_8),parameter :: ninehundred=900.0_real_8

Real(kind=real_8),parameter :: thousand=1000.0_real_8

end module
