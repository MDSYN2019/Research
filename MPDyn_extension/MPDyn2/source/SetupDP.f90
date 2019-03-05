! ############################
! ## SUBROUTINE LIST 
! ## -- SetupDP 
! ## -- Read_ColloidModel 
! ## -- Read_MolArt 
! ## -- SetupDP_IO 
! ## -- Gene_ConfigDP 
! ## -- Read_ConfigDP 
! ############################


!######################################################################
!######################################################################


subroutine SetupDP

use Numbers, only : N, Nf, NumSpec, NumMol, NumAtm
use CommonBlocks, only : QPBC, QRigidBody, ForceField
use Configuration
use CommonDPD
use CellListMethod
use RBparam
use IOparam, only : Parameter_file
use BondedParam, only : NumBond, BondI, BondJ, kBond, rBond, Frc_Bond
use AtomParam, only : MolName

implicit none

integer :: Nc, NumAllType, NumAllType1, Num
integer :: i, j, k, l
integer :: MaxNc, ib, jj, ia, i1, i2
integer, dimension(:), allocatable :: NumPair


   N = 0
   do i = 1 , NumSpec

     N = N + NumMol(i) * NumAtm(i)

   end do

   allocate( TypeNum(N) )
   allocate( Ncal(N) )

   allocate( R  (3,N) )
   allocate( Vel(3,N) )

   allocate( FrcDP (3,N) )
   allocate( FrcDPt(3,N) )
   allocate( Frc_Bond(3,N) )

   MaxNc = NumAtm(1)

   if(NumSpec >= 2) then

     do i = 2, NumSpec

       if(NumAtm(i) > MaxNc) then
         MaxNc = NumAtm(i)
       end if

     end do

   end if

   allocate( NumBondSpec(NumSpec) )
   allocate( BondPairS(2,MaxNc,NumSpec) )
   allocate( kBondSpec(MaxNc,NumSpec) )
   allocate( rBondSpec(MaxNc,NumSpec) )

   allocate( ColloidFlag(NumSpec) )

   ColloidFlag = .False.
   QColloid    = .False.
   NumColloid  = 0

   k = 0
   l = 0
   Nf = 0
   Num = 0

   do i = 1 , NumSpec

     Nc = NumMol(i)

     if(MolName(i)(1:4) == 'Coll') then

       l = l + 1
       Num = Num + 1

       ColloidFlag(i) = .True.
       QColloid = .True.

       NumColloid = Nc
       NumCollAtm = NumAtm(i)

       call Read_ColloidModel(i)

       allocate( R_RB(3,NumColloid) )
       allocate( V_RB(3,NumColloid) )
       allocate( Quaternion(4,NumColloid) )
       allocate( Rotation(3,3,NumColloid) )
       allocate( Lmoment(3,NumColloid) )
       allocate( Omega(3,NumColloid) )
       allocate( Rmolec(3,NumCollAtm,NumColloid) )

       allocate( FrcCo(3,NumColloid) )
       allocate( FrcCot(3,NumColloid) )
       allocate( Lmomentt(3,NumColloid) )
       allocate( Omegat(3,NumColloid) )
       allocate( V_RBt(3,NumColloid) )
       allocate( Torque(3,NumColloid) )

       do j = 1, NumColloid * NumCollAtm

         k = k + 1
         TypeNum(k) = l

         Ncal(k) = Num

       end do

       NumBondSpec(i) = 0

       Nf = Nf + NumColloid * 6

     else if(NumAtm(i) /= 1) then

       call Read_MolArt(i,k,l,Num)

       Nf = Nf + NumMol(i) * NumAtm(i) * 3

     else

       l = l + 1

       do j = 1, NumMol(i) * NumAtm(i)

         k = k + 1
         TypeNum(k) = l

         Num = Num + 1
         Ncal(k) = Num

       end do

       NumBondSpec(i) = 0

       Nf = Nf + NumMol(i) * NumAtm(i) * 3

     end if

   end do

! #####################################

   NumBond = 0

   do i = 1 , NumSpec

     NumBond = NumBond + NumBondSpec(i) * NumMol(i)

   end do

   if(NumBond /= 0) then

     allocate( BondI(NumBond) )
     allocate( BondJ(NumBond) )
     allocate( kBond(NumBond) )
     allocate( rBond(NumBond) )

     ib = 0
     jj = 0
     do i = 1 , NumSpec

       if(NumBondSpec(i)/=0) then

         do j = 1 , NumMol(i)

           ia = (j-1)*NumAtm(i) + jj

           do k = 1 , NumBondSpec(i)

             ib = ib + 1

             BondI(ib) = BondPairS(1,k,i) + ia
             BondJ(ib) = BondPairS(2,k,i) + ia

             kBond(ib) = kBondSpec(k,i)
             rBond(ib) = rBondSpec(k,i)

           end do

         end do

       end if

       jj = jj + NumMol(i)*NumAtm(i)

     end do

   end if

   allocate( NpBond(4,N) )
   allocate( NumPair(N) )

   NpBond = 0
   NumPair = 0

   if(NumBond /= 0) then

     do j = 1 , NumBond

       i1 = BondI(j)
       i2 = BondJ(j)

       NumPair(i1)=NumPair(i1)+1
       NumPair(i2)=NumPair(i2)+1

       if((NumPair(i1)>4).or.(NumPair(i2)>4)) then

         write( 6,*) 'ERROR : Too many hands (illegal bond definition) '
         stop

       end if

       NpBond(NumPair(i1),i1)=i2
       NpBond(NumPair(i2),i2)=i1

     end do

   end if

   deallocate( NumPair )

! #####################################

   call Read_Condition(3)

   NumAllType = TypeNum(N)

   open(3,file=trim(Parameter_file),status='old')

   if(ForceField(1:1)=='F') then

     allocate( a( NumAllType, NumAllType) )

     read(3,*) NumAllType1
     read(3,*)

!     if(ForceField(1:6)=='Smooth') then
!       read(ForceField(7:7),*) nfact
!       read(ForceField(8:8),*) mfact
!     end if
!
!     write(*,*) 'nfact = ',nfact
!     write(*,*) 'mfact = ',mfact

     if(NumAllType1 /= NumAllType) then

       write(*,*) 'ERROR : integeraction.param'
       stop

     end if

     do i = 1, NumAllType

       read(3,*) a(i,:)

     end do

   else if(ForceField=='Morse') then

     allocate( e_morse( NumAllType, NumAllType ) )

     read(3,*) NumAllType1
     read(3,*)

     if(NumAllType1 /= NumAllType) then

       write(*,*) 'ERROR : integeraction.param'
       stop

     end if

     read(3,*) a_morse

     do i = 1, NumAllType

       read(3,*) e_morse(i,:)

     end do

   end if

   close(3)

! ## set as default
   QRigidBody = .False.
   QPBC = .True.

end subroutine SetupDP


!######################################################################
!######################################################################


subroutine Read_ColloidModel(is)

use Numbers, only : NumAtm
use CommonDPD
use RBparam

implicit none

integer :: i, is, Natm

   open(3,file='./param/ColloidModel.data',status='old')

   read(3,'(i5/)') Natm

   if(Natm /= NumAtm(is)) then

     write(*,*) 'ERROR : defined model of colloid particle'
     stop

   end if

   allocate( R_onMol(3,Natm,1) )
   allocate( InertiaRB(3,1) )
   allocate( InvInertiaRB(3,1) )
   allocate( MassRB(1) )
   allocate( InvMassRB(1) )

   do i = 1 , Natm

     read(3,'(3f10.5)') R_onMol(:,i,1)

   end do

   read(3,'( f10.5)') MassRB(1)
   read(3,'(3f10.5)') InertiaRB(:,1)

   InvMassRB    = 1.d0 / MassRB
   InvInertiaRB = 1.d0 / InertiaRB

   close(3)

end subroutine Read_ColloidModel


!######################################################################
!######################################################################


subroutine Read_MolArt(is,k,l,Num)

use Numbers, only : NumMol, NumAtm
use CommonDPD
use AtomParam, only : MolName

implicit none

integer :: i, j, is, k, l, Num
character(len=80) :: String
integer :: NumSType, NumAtoms, NumBond0
integer, dimension(:), allocatable :: AtomSpec
integer, dimension(:,:), allocatable :: BondPair
real(8), dimension(:), allocatable :: kBond0, rBond0

   write(String,'(a,a,a)') './param/',trim(adjustl(MolName(is))),'Model.data'

   open(3,file=trim(String),status='old')

! ####### Read from MolArt File ######

   read(3,'(3i5/)') NumAtoms,NumSType,NumBondSpec(is)

   NumBond0 = NumBondSpec(is)

   if(NumAtoms /= NumAtm(is)) then

     write(*,*) 'ERROR : defined model'
     stop

   end if

   allocate( AtomSpec(NumAtoms) )
   allocate( BondPair(2,NumBond0) )
   allocate( kBond0(NumBond0) )
   allocate( rBond0(NumBond0) )

   do i = 1, NumAtm(is)

     read(3,'(i5)') AtomSpec(i)

   end do

   read(3,*)

   do i = 1 , NumBondSpec(is)

     read(3,'(2i5,2f8.2)') BondPairS(1,i,is),BondPairS(2,i,is),&
     &                     kBondSpec(i,is),rBondSpec(i,is)

   end do
! #####################################

   do i = 1 , NumMol(is)

     do j = 1 , NumAtm(is)

       k = k + 1
       TypeNum(k) = l + AtomSpec(j)

       Num = Num + 1
       Ncal(k) = Num

     end do

   end do

   l = l + NumSType

end subroutine Read_MolArt


!######################################################################
!######################################################################


subroutine SetupDP_IO

use Numbers, only : N
use CommonBlocks, only : QInitial
use CommonDPD
use CellListMethod

implicit none

! # create map of cell list
! -------------------------------------------
     call Mapping
     call Mapping_TopLine
! -------------------------------------------

! # initial configuration
! -------------------------------------------
   SlideGap = 0.d0

   if((QInitial).and.(TypeNum(N)==1)) then

     call Gene_ConfigDP

   else

     call Read_ConfigDP

   end if
! -------------------------------------------

end subroutine SetupDP_IO


!#####################################################################
!#####################################################################


subroutine Gene_ConfigDP

use Numbers, only : N
use Configuration
use CommonDPD
use BathParam, only : Temp_o
use CellParam, only : CellL, InvCL
use TimeParam, only : Timeps

implicit none

integer, parameter :: maxtry= 50000
real(8), parameter :: rcrt = 0.5d0
real(8), dimension(3) :: CellLh
real(8), dimension(3) :: x, Ri, Sumv
real(8) :: rcrt2, R2
integer :: mtry, i, j
real(8) :: Gauss
real(8), dimension(3) :: Rij, Ori
real(8) :: ranf, Sc
external ranf
integer, dimension(3) :: Ibox,IL
integer :: Nbox, Ubox, TryBox, Ixy, count
real(8), dimension(3) :: CL
integer, dimension(:), allocatable :: OcIndex

   Ibox = CellL / rcrt

   CL   = CellL / dble(Ibox)

   Nbox = Ibox(1) * Ibox(2) * Ibox(3)

   allocate( OcIndex(Nbox) )

   OcIndex = 0

   Ubox = Nbox

   CellLh = CellL * 0.5

! random configuration
! generate initial configuration in cell unit

   write( 6,*) 'initial conf. ----> random'

   rcrt2 = rcrt * rcrt

   do i = 1 , N

     mtry=0
  60 mtry=mtry+1

     if(mtry >= maxtry) then

       write( 6,*) 'ERROR:mtry exceeds maxtry!'
#ifndef BMONI
       write(11,*) 'ERROR:mtry exceeds maxtry!'
#endif

       print *, i
       stop

     end if

!---------------------------------
     TryBox = ranf()*Ubox + 1

     count = 0

     do j = 1 , Nbox

       if(OcIndex(j)==0) then

         count = count + 1

         if(count==TryBox) then

           TryBox = j
           exit

         end if

       end if

     end do

     IL(3) = (TryBox-1) / (Ibox(1)*Ibox(2)) + 1
     Ixy   = mod(TryBox-1,Ibox(1)*Ibox(2)) + 1
     IL(2) = (Ixy-1) / Ibox(1) + 1
     IL(1) = mod(Ixy-1,Ibox(1)) + 1

     Ori = (IL - 0.5d0) * CL - CellLh

     call Randp(x(1),x(2),x(3))
!---------------------------------
!
     Ri = x * CL * 0.5d0 + Ori

     do j = 1 , i-1

       Rij = Ri - R(:,j)
       Rij = Rij - dnint(Rij*InvCL) * CellL
       R2  = dot_product( Rij, Rij )

       if(R2<=rcrt2) go to 60

     end do

     R(:,i) = Ri
     OcIndex(TryBox) = 1
     Ubox = Ubox - 1

   end do

   write(6,*) 'initial configuration has just been generated'

!
! assign random velocity
!
   Sc = sqrt(Temp_o)
!
   do i = 1 , N

     Vel(1,i) = Sc * Gauss()
     Vel(2,i) = Sc * Gauss()
     Vel(3,i) = Sc * Gauss()

   end do

!    ** remove net momentum **

   Sumv = 0.d0

   do i = 1 , N

     Sumv = Sumv + Vel(:,i)

   end do

   Sumv = Sumv / dble(N)

   do i = 1 , N

     Vel(:,i) = Vel(:,i) - Sumv

   end do

   Timeps = 0.d0

   write(6,*) 'Initialization of the configuration has been completed!'

   call Write_ConfigDP

end subroutine Gene_ConfigDP


!######################################################################
!######################################################################


subroutine Read_ConfigDP

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration
use CommonDPD
use RBparam
use TimeParam, only : Timeps

implicit none

integer :: i, j, k, l, Num

! #####################################
! #      read config from disk        #
! #####################################

open(1,file='restartDP.dat',form='formatted',status='old')

   l = 0
   Num = 0
   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       do j = 1 , NumMol(i)

         Num = Num + 1
         read(1,*) R_RB(:,j),V_RB(:,j)
         read(1,*) Quaternion(:,j),Lmoment(:,j)

         do k = 1 , NumAtm(i)

           l = l + 1
           Ncal(l) = Num

         end do

       end do

     else

       do j = 1, NumMol(i)

         do k = 1 , NumAtm(i)

           l = l + 1
           Num = Num + 1
           read(1,*) R(:,l),Vel(:,l)
           Ncal(l) = Num

         end do

       end do

     end if

   end do

!   do i = 1, N
!     R(:,i) = R(:,i) - CellL(:) * 0.5d0
!   end do

!   read(1,'(f10.6)') SlideGap
   read(1,*) SlideGap

!   read(1,'(f15.3)') timeps
   read(1,*) Timeps

   write(11,'(a20,f15.3)') 'Restart at time t=',timeps

close(1)

end subroutine Read_ConfigDP
