! ############################
! ## SUBROUTINE LIST 
! ## -- MakeBodyFixedCoordinate 
! ## -- UnitCoord 
! ## -- Write_Unit 
! ## -- Read_OPLS_Topology 
! ## -- Read_OPLS_Parameter 
! ## -- Read_Charmm_Parameter 
! ## -- Read_Charmm_Topology 
! ############################


module MKBFXparam
   character(len=10) :: ForceField
   character(len=80) :: Topology_file     ! topology file
   character(len=80) :: Parameter_file    ! parameter file
   integer :: NumAtomTypeParam
   integer :: NumUnit
   character(len=4), dimension(300) :: AtomTypeParam
   character(len=6), dimension(300) :: UNITName
   real(8), dimension(300) :: MassParam
   integer, dimension(300) :: NumAtom_inUNITParam, NumAtom_inUNIT
   character(len=4), dimension(300,300) :: AtomNameParam, AtomType
   integer :: NumBondParam
   integer :: NumAngleParam
   character(len=4), dimension(2,500) :: BondPairAtoms
   real(8), dimension(500) :: kBondParam, rBondParam
   character(len=4), dimension(3,800) :: AnglePairAtoms
   real(8), dimension(800) :: kThetaParam, Theta0Param
   real(8), parameter :: pi = 3.14159265358979d0
end module MKBFXparam


program MakeBodyFixedCoordinate

use MKBFXparam

implicit none

character(len=80) :: String1, String

open(11,file='./param/RigidBodyModel.prm')

   print *, 'Which force field do you use ?'
   read(5,*) ForceField

   print *, 'Name of PARAMETER file ? if you use a default file, write "D"'
   read(5,'(a80)') String1
   String = adjustl(String1)

   if(String(1:1)=='D') then
     if(ForceField(1:5)=='CHARM') then
     write(Parameter_file,*) 'param/par_charmm27.prm'
     else if(ForceField(1:4)=='OPLS') then
     write(Parameter_file,*) 'param/par_oplsaa.prm'
     else if(ForceField(1:3)=='BKS') then
     write(Parameter_file,*) 'param/par_BKS.prm'
     end if
   else
     read(String,*) Parameter_file
   end if

   print *, 'Name of TOPOLOGY file ? if you use a default file, write "D"'
   read(5,'(a80)') String1
   String = adjustl(String1)

   if(String(1:1)=='D') then
     if(ForceField(1:5)=='CHARM') then
     write(Topology_file,*) 'param/top_charmm27.prm'
     else if(ForceField(1:4)=='OPLS') then
     write(Topology_file,*) 'param/top_oplsaa.prm'
     end if
   else
     read(String,*) Topology_file
   end if

   if(ForceField(1:4) == 'OPLS') then
     call Read_OPLS_Parameter
     call Read_OPLS_Topology
   else if(ForceField(1:5) == 'CHARM') then
     call Read_Charmm_Parameter
     call Read_CHARMM_Topology
   end if

   call UnitCoord

   write(11,'(a)') '!'
   write(11,'(a)') '<end>'

close(11)

end program MakeBodyFixedCoordinate


!######################################################################
!######################################################################


subroutine UnitCoord

use MKBFXparam

implicit NONE

real(8), dimension(3,100) :: R
real(8), dimension(100) :: Mass
logical, dimension(100) :: Duplicate
character(len=4) :: AName
real(8) :: cst, snt, RCC, RHH, RCH, ROH, theta, X
real(8) :: RCHx, RCHy
integer :: i , j, k, Nc

   Duplicate(1) = .False.
   do i = 2 , NumUnit
     Duplicate(i) = .False.
     do j = 1, i-1
       if(UnitName(j)==UnitName(i)) then
         Duplicate(i) = .True.
         exit
       end if
     end do
   end do

   do i = 1 , NumUnit
     if(UnitName(i) == 'SNGL') then
       Duplicate(i) = .True.
     end if
   end do

   do i = 1 , NumUnit
     if(NumAtom_inUNITParam(i)/=NumAtom_inUNIT(i)) then
       write(*,*) 'ERROR : number of atom'
       write(*,*) 'unit number = ',i
       write(*,*) NumAtom_inUNITParam(i), NumAtom_inUNIT(i)
       stop
     end if
   end do

   do i = 1 , NumUnit

     if(Duplicate(i)) cycle
     Nc = NumAtom_inUNIT(i)

     do j = 1 , Nc

       AName = AtomType(j,i)

       do k = 1 , NumAtomTypeParam

         if(AtomTypeParam(k) == AName) then

           Mass(j) = MassParam(k)
           exit

         end if

       end do

     end do

     if(UnitName(i) == 'PPEBZ') then

       snt = sin(pi/3.)
       cst = cos(pi/3.)

       call BondIdent(AtomType(1,i),AtomType(2,i),RCC)
       call BondIdent(AtomType(2,i),AtomType(3,i),RCH)

       R = 0.d0

       R(1, 1) = -       RCC

       R(1, 2) = - cst * RCC
       R(2, 2) =   snt * RCC

       R(1, 3) = - cst * (RCC+RCH)
       R(2, 3) =   snt * (RCC+RCH)

       R(1, 4) =   cst * RCC
       R(2, 4) =   snt * RCC

       R(1, 5) =         RCC

       R(1, 6) =   cst * RCC
       R(2, 6) = - snt * RCC

       R(1, 7) = - cst * RCC
       R(2, 7) = - snt * RCC

       R(1, 8) = - cst * (RCC+RCH)
       R(2, 8) = - snt * (RCC+RCH)

       call adjustpos(Mass,R,Nc)
       call Write_UNIT(UnitName(i),R,Mass,Nc,i)

     else if(UnitName(i) == 'METL') then

       call BondIdent(AtomType(1,i),AtomType(2,i),RCH)
       call AnglIdent(AtomType(2,i),AtomType(1,i),AtomType(3,i),theta)
       RHH = 2.d0 * RCH * sin( theta*0.5d0 )

       R = 0.d0
       snt = sin(pi/3.)

       X = RHH * sin(pi/3.) * 2. / 3.
       R(3,1) = sqrt( RCH*RCH - X*X )
       R(1,2) = - X
       R(1,3) =   RHH * sin(pi/3.) * 1. / 3.
       R(1,4) =   RHH * sin(pi/3.) * 1. / 3.
       R(2,3) = + RHH * 0.5
       R(2,4) = - RHH * 0.5
       if(ForceField(1:5)=='CHARM') then
         R(2,3) = - RHH * 0.5
         R(2,4) = + RHH * 0.5
       end if

       call adjustpos(Mass,R,Nc)
       call Write_UNIT(UnitName(i),R,Mass,Nc,i)

     else if(UnitName(i) == 'PPEBI') then

       snt = sin(pi/3.)
       cst = cos(pi/3.)

       call BondIdent(AtomType(1,i),AtomType(2,i),RCH)
       call BondIdent(AtomType(2,i),AtomType(3,i),RCC)

       R = 0.d0

       R(1, 1) = -       (RCC+RCH)

       R(1, 2) = -       RCC

       R(1, 3) = - cst * RCC
       R(2, 3) =   snt * RCC

       R(1, 4) = - cst * (RCC+RCH)
       R(2, 4) =   snt * (RCC+RCH)

       R(1, 5) =   cst * RCC
       R(2, 5) =   snt * RCC

       R(1, 6) =         RCC

       R(1, 7) =   cst * RCC
       R(2, 7) = - snt * RCC

       R(1, 8) = - cst * RCC
       R(2, 8) = - snt * RCC

       R(1, 9) = - cst * (RCC+RCH)
       R(2, 9) = - snt * (RCC+RCH)

       call adjustpos(Mass,R,Nc)
       call Write_UNIT(UnitName(i),R,Mass,Nc,i)

     else if(UnitName(i) == 'OH') then

       call BondIdent(AtomType(1,i),AtomType(2,i),ROH)

       R = 0.d0
       R(1, 2) = ROH

       call adjustpos(Mass,R,Nc)
       call Write_UNIT(UnitName(i),R,Mass,Nc,i)

     else if((UnitName(i) == 'TIP3').or.(UnitName(i) == 'SPCE')) then

       call BondIdent(AtomType(1,i),AtomType(2,i),ROH)
       call AnglIdent(AtomType(2,i),AtomType(1,i),AtomType(3,i),theta)
       RHH = 2.d0 * ROH * sin( theta*0.5d0 )

       R = 0.d0

       R(1,2) = - RHH / 2.
       R(1,3) =   RHH / 2.
       R(2,2) = - ROH * cos( theta*0.5d0 )
       R(2,3) = - ROH * cos( theta*0.5d0 )

       call adjustpos(Mass,R,Nc)
       call Write_UNIT(UnitName(i),R,Mass,Nc,i)

     else if(UnitName(i) == 'H3O') then ! NOTE urata's numbering H, O, H, H

       call BondIdent(AtomType(1,i),AtomType(2,i),ROH)
       call AnglIdent(AtomType(1,i),AtomType(2,i),AtomType(3,i),theta) ! HOH 
       RHH = 2.d0 * ROH * sin( theta*0.5d0 )

       R = 0.d0
       snt = sin(pi/3.)

       X = RHH * sin(pi/3.) * 2. / 3.
       R(3,2) = sqrt( ROH*ROH - X*X )
       R(1,1) = - X
       R(1,3) =   RHH * sin(pi/3.) * 1. / 3.
       R(1,4) =   RHH * sin(pi/3.) * 1. / 3.
       R(2,3) = + RHH * 0.5
       R(2,4) = - RHH * 0.5

       call adjustpos(Mass,R,Nc)
       call Write_UNIT(UnitName(i),R,Mass,Nc,i)

     else if(UnitName(i) == 'BENZ') then

       snt = sin(pi/3.)
       cst = cos(pi/3.)

       call BondIdent(AtomType(1,i),AtomType(2,i),RCH)
       call BondIdent(AtomType(1,i),AtomType(3,i),RCC)

       R = 0.d0

       R(1, 1) = -       RCC

       R(1, 2) = -      (RCC+RCH)

       R(1, 3) = - cst * RCC
       R(2, 3) =   snt * RCC

       R(1, 4) = - cst * (RCC+RCH)
       R(2, 4) =   snt * (RCC+RCH)

       R(1, 5) =   cst * RCC
       R(2, 5) =   snt * RCC

       R(1, 6) =   cst * (RCC+RCH)
       R(2, 6) =   snt * (RCC+RCH)

       R(1, 7) =         RCC

       R(1, 8) =         (RCC+RCH)

       R(1, 9) =   cst * RCC
       R(2, 9) = - snt * RCC

       R(1,10) =   cst * (RCC+RCH)
       R(2,10) = - snt * (RCC+RCH)

       R(1,11) = - cst * RCC
       R(2,11) = - snt * RCC

       R(1,12) = - cst * (RCC+RCH)
       R(2,12) = - snt * (RCC+RCH)

       call adjustpos(Mass,R,Nc)
       call Write_UNIT(UnitName(i),R,Mass,Nc,i)

     else if(UnitName(i) == 'MET2') then

       call BondIdent(AtomType(1,i),AtomType(2,i),RCH)
       call AnglIdent(AtomType(2,i),AtomType(1,i),AtomType(3,i),theta)
       RCHx = RCH * sin( theta*0.5d0 )
       RCHy = RCH * cos( theta*0.5d0 )

       R = 0.d0

       R(1,2) = - RCHx
       R(2,2) = - RCHy
       R(1,3) =   RCHx
       R(2,3) = - RCHy

       call adjustpos(Mass,R,Nc)
       call Write_UNIT(UnitName(i),R,Mass,Nc,i)

     else if(UnitName(i) == 'METH') then

       call BondIdent(AtomType(1,i),AtomType(2,i),RCH)

       R = 0.d0

       R(1,2) = RCH

       call adjustpos(Mass,R,Nc)
       call Write_UNIT(UnitName(i),R,Mass,Nc,i)

!     else if(UnitName(i) == '') then
!     else if(UnitName(i) == '') then
!     else if(UnitName(i) == '') then
!     else if(UnitName(i) == '') then
     end if

   end do

contains

   subroutine BondIdent(IName,JName,RR)

   integer :: l
   character(len=4) :: IName, JName
   real(8) :: RR

     do l = 1 , NumBondParam
         if(((BondPairAtoms(1,l)==IName).and.(BondPairAtoms(2,l)==JName))&
      & .or.((BondPairAtoms(2,l)==IName).and.(BondPairAtoms(1,l)==JName))) then
         RR = rBondParam(l)
         exit
         end if
         if(l==NumBondParam) then
           write(*,*) 'ERROR : missing bond parameters'
           write(*,*) IName,JName
           stop
         end if
     end do

   end subroutine BondIdent

   subroutine AnglIdent(IName,JName,KName,RR)

   integer :: l
   character(len=4) :: IName, JName, KName
   real(8) :: RR

     do l = 1 , NumAngleParam

       if( AnglePairAtoms(2,l) == JName ) then

         if(((AnglePairAtoms(1,l) == IName)  &
     &  .and.(AnglePairAtoms(3,l) == KName)) &
     &  .or.((AnglePairAtoms(1,l) == KName)  &
     &  .and.(AnglePairAtoms(3,l) == IName))) then

         RR = Theta0Param(l)
         exit
         end if
       end if
       if(l==NumAngleParam) then
         write(*,*) 'ERROR : missing angle parameters'
         write(*,*) IName,JName,KName
         stop
       end if
     end do

   end subroutine AnglIdent


   subroutine adjustpos(Mass,R,Nc)

   real(8), dimension(3,3) :: Wi
   real(8), dimension(3) :: Rg
   real(8) :: SumM
   real(8), dimension(3,100) :: R
   real(8), dimension(100) :: Mass
   integer :: Nc, j

       Rg = 0.
       SumM = 0.
       do j = 1 , Nc
         Rg = Rg + Mass(j) * R(:,j)
         SumM = SumM + Mass(j)
       end do
       Rg = Rg / SumM

       do j = 1 , Nc
         R(:,j) = R(:,j) - Rg
       end do

       do j = 1 , Nc

         Wi(1,1) = Wi(1,1) + Mass(j) * (R(2,j) * R(2,j) + R(3,j) * R(3,j))
         Wi(2,2) = Wi(2,2) + Mass(j) * (R(3,j) * R(3,j) + R(1,j) * R(1,j))
         Wi(3,3) = Wi(3,3) + Mass(j) * (R(1,j) * R(1,j) + R(2,j) * R(2,j))
         Wi(1,2) = Wi(1,2) - Mass(j) * R(1,j) * R(2,j)
         Wi(1,3) = Wi(1,3) - Mass(j) * R(1,j) * R(3,j)
         Wi(2,3) = Wi(2,3) - Mass(j) * R(2,j) * R(3,j)

       end do

       Wi(2,1) = Wi(1,2)
       Wi(3,1) = Wi(1,3)
       Wi(3,2) = Wi(2,3)

       print * , UnitName(i)
       print * , 'Wi_diag = ',Wi(1,1),Wi(2,2),Wi(3,3)
       print * , 'Wi_offd = ',Wi(2,1),Wi(3,1),Wi(3,2)

   end subroutine adjustpos

end subroutine UnitCoord


!######################################################################
!######################################################################


subroutine Write_Unit(UName,R,Mass,Nc,ii)

use MKBFXparam

implicit none

integer :: i, ii, Nc
character(len=6) :: UName
character(len=4) :: AName
real(8), dimension(3,100) :: R
real(8), dimension(100) :: Mass
character(len=20) :: FileName

   write(11,'(a)') '!'
   write(11,'(a,3x,a,i6)') 'UNIT',UName,Nc
   do i = 1 , Nc
     write(11,'(i5,f9.3,3f15.10)') i,Mass(i),R(:,i)
   end do

   write(FileName,'(a,a,a)') './check/',trim(UName),'.xyz'

   open(2,file=trim(FileName))
   write(2,'(i5/)') Nc
   do i = 1 , Nc
     AName= AtomType(i,ii)
     write(2,'(a,x,3f8.3)') AName(1:1),R(:,i)
   end do
   close(2)

end subroutine Write_Unit


!######################################################################
!######################################################################


! **************************
! ** OPLS Topology File   **
! **************************

subroutine Read_OPLS_Topology

use MKBFXparam

implicit NONE

character(len=4) :: Dummy1

integer :: i
character(len=80) :: String
character(len=1) :: Wk1
character(len=1) :: Wk2
character(len=1) :: Wk3
character(len=1) :: Wk4

! Topology

  open(1,file=trim(Topology_file),status='old')

   Wk1 = '*'
   Wk2 = '!'
   Wk3 = ' '
   Wk4 = '	'

! ---------------------------------------------------------------------

   NumAtomTypeParam=0
   NumUnit=0

   do

     read(1,'(80a)') String

     if(String(1:3)=='END') exit

     if(( String(1:5) == 'GROUP' ) .or. ( String(1:1) == Wk3 ) ) cycle

     if( String(1:4) == 'MASS' ) then

       NumAtomTypeParam = NumAtomTypeParam + 1
       read(String,*) Dummy1,i,AtomTypeParam(NumAtomTypeParam),&
       &          MassParam(NumAtomTypeParam)

     end if

! ---------------------------------------------------------------------

     if( String(1:4) == 'UNIT' ) then

       NumUnit = NumUnit + 1

       read(String,*) Dummy1, &
       &   UNITName(NumUnit), &
       &   NumAtom_inUNITParam(NumUnit)

       NumAtom_inUNIT(NumUnit) = 0

       cycle

     end if

     if( String(1:4) == 'ATOM' ) then

       NumAtom_inUNIT(NumUnit) = NumAtom_inUNIT(NumUnit) + 1

       read(String,*) Dummy1,                                              &
       &          AtomNameParam(NumAtom_inUNIT(NumUnit),NumUnit),&
       &          AtomType(NumAtom_inUNIT(NumUnit),NumUnit)

       cycle

     end if

   end do

! ---------------------------------------------------------------------

close(1)

   write(*,*) 'NumAtomTypeParam = ',NumAtomTypeParam
   write(*,*) 'NumUNIT  = ',NumUNIT

end subroutine Read_OPLS_Topology


!######################################################################
!######################################################################


! ***************************
! ** OPLS Parameter File   **
! ***************************

subroutine Read_OPLS_Parameter

use MKBFXparam

implicit NONE

integer :: i, j
character(len=72) :: String, String1
character(len=1) :: Wk1
character(len=1) :: Wk2
character(len=1) :: Wk3
character(len=1) :: Wk4

! Parameters

  open(1,file=trim(Parameter_file),status='old')

   Wk1 = '*'
   Wk2 = '!'
   Wk3 = ' '
   Wk4 = '	'

   do

     read(1,'(72a)') String

     if(String(1:5)=='BONDS') exit

   end do

! ----------------------------------------------------------------------
! ## Bond parameters
!
!    V(bond) = Kb(b - b0)**2
!
!    Kb: kcal/mole/A**2
!    b0: A
!
!    atom type    Kb          b0

   NumBondParam=0

   do

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:5)=='ANGLE') exit !until "ANGLE" is found

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumBondParam = NumBondParam + 1

     read(String,*) BondPairAtoms(1,NumBondParam), &
     &              BondPairAtoms(2,NumBondParam), &
     &              kBondParam(NumBondParam),      &
     &              rBondParam(NumBondParam)

   end do

   write(*,*) 'NumBondParam=',NumBondParam
!   do i = 1 , NumBondParam
!     write(*,*) BondPairAtoms(1,i),BondPairAtoms(2,i),rBondParam(i)
!   end do
! ----------------------------------------------------------------------
! ## Angle paramters
!
!    V(angle) = Ktheta(Theta - Theta0)**2
!
!    Ktheta: kcal/mole/rad**2
!    Theta0: degrees
!
!  atom types      Ktheta    Theta0

   NumAngleParam=0

   do

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:9)=='DIHEDRALS') exit !until "DIHEDRALS" is found

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumAngleParam = NumAngleParam + 1

     read(String,*) (AnglePairAtoms(j,NumAngleParam),j=1,3),  &
     &               kThetaParam(NumAngleParam),              &
     &               Theta0Param(NumAngleParam)

   end do

   write(*,*) 'NumAngleParam=',NumAngleParam

! ----
! Unit
! ----

   do i = 1, NumAngleParam

     Theta0Param(i) = Theta0Param(i) / 180.d0 * pi

   end do

close(1)

end subroutine Read_OPLS_Parameter


!######################################################################
!######################################################################


! *****************************
! ** Charmm Parameter File   **
! *****************************

subroutine Read_Charmm_Parameter

use MKBFXparam

implicit NONE

integer :: i, j
character(len=72) :: String, String1
character(len=1) :: Wk1
character(len=1) :: Wk2
character(len=1) :: Wk3
character(len=1) :: Wk4
character(len=72) :: Dummy

! Parameters

open(1,file=trim(Parameter_file),status='old')

   Wk1 = '*'
   Wk2 = '!'
   Wk3 = ' '
   Wk4 = '	'

   do

     read(1,'(72a)') String
     if(String(1:5)=='BONDS') exit

   end do

! ----------------------------------------------------------------------
! ## Bond parameters
!
!    V(bond) = Kb(b - b0)**2
!
!    Kb: kcal/mole/A**2
!    b0: A
!
!    atom type    Kb          b0

   NumBondParam=0

   do

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:5)=='ANGLE') exit !until "ANGLE" is found

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumBondParam = NumBondParam + 1

     read(String,*) BondPairAtoms(1,NumBondParam),  &
     &              BondPairAtoms(2,NumBondParam),  &
     &              kBondParam(NumBondParam),       &
     &              rBondParam(NumBondParam)

   end do

   write(*,*) 'NumBondParam=',NumBondParam
! ----------------------------------------------------------------------
! ## Angle paramters
!
!    V(angle) = Ktheta(Theta - Theta0)**2
!
!    V(Urey-Bradley) = Kub(S - S0)**2
!
!    Ktheta: kcal/mole/rad**2
!    Theta0: degrees
!    Kub: kcal/mole/A**2 (Urey-Bradley)
!    S0: A
!
!  atom types      Ktheta    Theta0    Kub      S0

   NumAngleParam=0

   do

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:9)=='DIHEDRALS') exit !until "DIHEDRALS" is found

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumAngleParam = NumAngleParam + 1

     read(String,*) (AnglePairAtoms(j,NumAngleParam),j=1,3), &
     &              kThetaParam(NumAngleParam),              &
     &              Theta0Param(NumAngleParam), Dummy

   end do

   write(*,*) 'NumAngleParam=',NumAngleParam

! ----
! Unit
! ----

   do i = 1, NumAngleParam

     Theta0Param(i) = Theta0Param(i) / 180.d0 * pi

   end do

close(1)

end subroutine Read_Charmm_Parameter


!######################################################################
!######################################################################


! ****************************
! ** Charmm Topology File   **
! ****************************

subroutine Read_Charmm_Topology

use MKBFXparam

implicit NONE

character(len=4) :: Dummy1

integer :: i
character(len=72) :: Ch
character(len=1) :: Wk1
character(len=1) :: Wk2
character(len=1) :: Wk3
character(len=1) :: Wk4

! Topology

  open(1,file=trim(Topology_file),status='old')

   Wk1 = '*'
   Wk2 = '!'
   Wk3 = ' '
   Wk4 = '	'

! ---------------------------------------------------------------------

   NumAtomTypeParam=0
   NumUnit=0

   do

     read(1,'(72a)') Ch

     if(Ch(1:3)=='END') exit

     if(( Ch(1:5) == 'GROUP'    ) .or. ( Ch(1:5) == 'DONOR' ) .or.  &
     &  ( Ch(1:8) == 'ACCEPTOR' ) .or. ( Ch(1:1) == Wk3     ) .or.  &
     &  ( Ch(1:2) == 'IC'       ) ) cycle

     if( Ch(1:4) == 'MASS' ) then

       NumAtomTypeParam = NumAtomTypeParam + 1
       read(Ch,*) Dummy1,i,AtomTypeParam(NumAtomTypeParam),&
       &          MassParam(NumAtomTypeParam)

     end if

! ---------------------------------------------------------------------

     if( Ch(1:4) == 'UNIT' ) then

       NumUnit = NumUnit + 1

       read(Ch,*) Dummy1, &
       &   UNITName(NumUnit), &
       &   NumAtom_inUNITParam(NumUnit)

       NumAtom_inUNIT(NumUnit) = 0

       cycle

     end if


     if( Ch(1:4) == 'ATOM' ) then

       NumAtom_inUNIT(NumUnit) = NumAtom_inUNIT(NumUnit) + 1

       read(Ch,*) Dummy1,                                        &
       &          AtomNameParam(NumAtom_inUNIT(NumUnit),NumUnit),&
       &          AtomType(NumAtom_inUNIT(NumUnit),NumUnit)

       cycle

     end if

   end do

! ---------------------------------------------------------------------

close(1)

   write(*,*) 'NumAtomTypeParam = ',NumAtomTypeParam
   write(*,*) 'NumUNIT  = ',NumUNIT

end subroutine Read_Charmm_Topology
