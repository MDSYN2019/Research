! ############################
! ## SUBROUTINE LIST 
! ## -- RadiusOfGyrationProtein 
! ## -- AvConformProtein 
! ## -- ChangeAtomName 
! ## -- MSD_protein 
! ## -- MSD_proteinD 
! ## -- MSD_proteinR 
! ## -- Rot_to_standard 
! ## -- MSD_proteinR_detail 
! ## -- MSD_proteinMinim 
! ## -- MSD_proteinMinim_detail 
! ## -- MinimEuler 
! ## -- MSD_MinimP_MS 
! ############################


!######################################################################
!######################################################################


subroutine RadiusOfGyrationProtein

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : Mass
use TimeParam, only : Timeps

implicit none

integer :: i, j, k
real(8), dimension(3) :: Rg
real(8) :: Rg2, Rg3, R2 ,R3
real(8) :: SumM
integer :: Numa, Numb

   open(1,file='./Analy/RofG_Protein.dat',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   R = 0.d0

   SumM = 0.d0

   do i = Numa, Numb

     SumM = SumM + Mass(i)

   end do

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

       call CellTransform

       call COM_SpecComponent(Rg,Kcomp)

       do k = Numa, Numb

         R(:,k) = R(:,k) - Rg(:)

       end do

       Rg2 = 0.
       Rg3 = 0.

       do k = Numa, Numb

         R2 = sqrt( R(1,k)**2+R(2,k)**2 )
         R3 = sqrt( dot_product(R(:,k),R(:,k)) )

         Rg2 = Rg2 + Mass(k) * R2
         Rg3 = Rg3 + Mass(k) * R3

       end do

       Rg2 = Rg2 / SumM
       Rg3 = Rg3 / SumM

       write(1,'(f9.3,2f10.3)') Timeps, Rg2, Rg3

     end do

   end do

close(1)

end subroutine RadiusOfGyrationProtein


!######################################################################
!######################################################################


! ***********************************
! **  Average Structure of Protein **
! ***********************************

subroutine AvConformProtein(Rorig)

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : DefaultAtomName, ResidName, ResidNum
use TimeParam, only : Timeps

implicit none

integer :: i, j, k, Numa, Numb, StepCount
real(8), dimension(3,N) :: Rorig
real(8), dimension(3) :: Rcom

character(len=4) :: NameA
character(len=1), dimension(10), parameter :: &
&   Flpr = (/'A','B','C','D','E','F','G','H','I','J'/)
real(8), parameter :: zero=0.d0

open(1,file='./Analy/ProteinAvConf.data',status='unknown')
open(2,file='./Analy/TrajProtCOM.dat',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   call COM_SpecComponent(Rcom,Kcomp)

   do i = Numa , Numb
     R(:,i) = R(:,i) - Rcom
   end do

   Rorig = R

   allocate( Rav(3,Numb) )

! ---------------------------------------------------------------------------

   Rav = 0.d0
   StepCount = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepCount = StepCount + 1

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

       call CellTransform

       call COM_SpecComponent(Rcom,Kcomp)

       write(2,'(f7.2,3f9.4)') Timeps, Rcom

       do k = Numa , Numb

         R(:,k) = R(:,k) - Rcom

       end do

       call MinimEuler(Rorig,Numa,Numb)

       do k = Numa , Numb

         Rav(:,k) = Rav(:,k) + R(:,k)

       end do

     end do

   end do

   Rav = Rav / dble(StepCount)

   do i = Numa , Numb

     NameA = DefaultAtomName(i)

     call ChangeAtomName(NameA)

     if(ResidName(k)=='HSD') ResidName(k)='HIS'

       write(1,'(a4,2x,i5,2x,a3,x,a4,a1,i4,4x,3f8.3,2f6.2,10x,1a)') &
       & 'ATOM',i,NameA(1:3),ResidName(i),Flpr(i),ResidNum(i),      &
       &  Rav(:,i),zero,zero,NameA(1:1)

   end do

   write(1,'(a)') 'TER'

close(1)
close(2)

end subroutine AvConformProtein


!######################################################################
!######################################################################


subroutine ChangeAtomName(NameA)

implicit none

character(len=4) :: NameA

   if(NameA == 'HD11') NameA = 'HD1'
   if(NameA == 'HD12') NameA = 'HD2'
   if(NameA == 'HD13') NameA = 'HD3'
   if(NameA == 'HD21') NameA = 'HD4'
   if(NameA == 'HD22') NameA = 'HD5'
   if(NameA == 'HD23') NameA = 'HD6'

   if(NameA == 'HH11') NameA = 'HH1'
   if(NameA == 'HH12') NameA = 'HH2'
   if(NameA == 'HH21') NameA = 'HH3'
   if(NameA == 'HH22') NameA = 'HH4'

   if(NameA == 'HG11') NameA = 'HG1'
   if(NameA == 'HG12') NameA = 'HG2'
   if(NameA == 'HG13') NameA = 'HG3'
   if(NameA == 'HG21') NameA = 'HG4'
   if(NameA == 'HG22') NameA = 'HG5'
   if(NameA == 'HG23') NameA = 'HG6'

   if(NameA == 'HE11') NameA = 'HE1'
   if(NameA == 'HE12') NameA = 'H2'

end subroutine ChangeAtomName


!######################################################################
!######################################################################


! *******************************
! **  Mean square displacement **
! *******************************

subroutine MSD_protein

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : AtomName, ResidNum
use TimeParam, only : Timeps

implicit none

integer :: i , j, k
real(8), dimension(3,N) :: Rorig
real(8), dimension(3) :: Rcom, DR
real(8) :: MSD, MSD_TM
integer :: NumCA, NumCA_TM
logical, dimension(N) :: Flag_MT
integer :: Numa , Numb

   open(1,file='./Analy/MSD_protein.data',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   Rorig = R

   call COM_SpecComponent(Rcom,Kcomp)

   do i = Numa, Numb

     Rorig(:,i) = Rorig(:,i) - Rcom

   end do

   Flag_MT = .False.

   do i = Numa, Numb

     j = ResidNum(i)

     if(ResidFlag(j) == 1) then
       Flag_MT(i) = .True.
     end if

   end do

! ---------------------------------------------------------------------------

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

       call CellTransform

       call COM_SpecComponent(Rcom,Kcomp)

       do k = Numa, Numb

         R(:,k) = R(:,k) - Rcom

       end do

       NumCA = 0
       NumCA_TM = 0
       MSD = 0.d0
       MSD_TM = 0.d0

       do k = Numa, Numb

         if(AtomName(k) == 'CA') then

           NumCA = NumCA + 1
           DR    = R(:,k) - Rorig(:,k)
           MSD   = MSD + dot_product( DR, DR )

           if(Flag_MT(k)) then

             NumCA_TM = NumCA_TM + 1
             MSD_TM   = MSD_TM   + dot_product( DR, DR )

           end if

         end if

       end do

       MSD    = MSD    / NumCA
       MSD_TM = MSD_TM / NumCA_TM

       MSD    = sqrt(MSD)
       MSD_TM = sqrt(MSD_TM)

       write(1,'(f12.4,2f10.3)') Timeps, MSD_TM, MSD

     end do

   end do

   close(1)

end subroutine MSD_protein


!######################################################################
!######################################################################


! *******************************
! **  Mean square displacement **
! *******************************

subroutine MSD_proteinD

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : AtomName, ResidName, ResidNum

implicit none

integer, parameter :: Steps=10
integer :: i , j, k, l, StepCount, Nsteps
real(8), dimension(3,N) :: Rorig
real(8), dimension(3) :: Rcom, DR
real(8), dimension(:), allocatable :: MSDex, MSDav
logical, dimension(N) :: Flag_MT
real(8), dimension(:,:,:), allocatable :: RCa
real(8), dimension(:,:), allocatable :: RCaAv
character(len=4), dimension(:), allocatable :: RName
integer :: Numa , Numb

   open(1,file='./Analy/MSD_proteinD.data',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   Nsteps = Nsnap/Steps

   Rorig = R

   call COM_SpecComponent(Rcom,Kcomp)

   do i = Numa, Numb

     Rorig(:,i) = Rorig(:,i) - Rcom

   end do

   Flag_MT = .False.

   j = ResidNum(Numb)

   allocate( RCa(3,j,0:Nsteps) )
   allocate( RCaAv(3,j) )
   allocate( MSDex(j) )
   allocate( MSDav(j) )
   allocate( RName(j) )

   do i = Numa, Numb

     if(AtomName(i) == 'CA') then

       j = ResidNum(i)
       RCa(:,j,0) = Rorig(:,i)

       if(ResidFlag(j) == 1) then

         Flag_MT(j) = .True.

       end if

       RName(j) = ResidName(i)

     end if

   end do

! ---------------------------------------------------------------------------

   StepCount = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepCount = StepCount + 1

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

       if(mod(StepCount,Steps)/=0) cycle

       call CellTransform

       call COM_SpecComponent(Rcom,Kcomp)

       do k = Numa, Numb

         R(:,k) = R(:,k) - Rcom

       end do

       do k = Numa, Numb

         if(AtomName(k) == 'CA') then

           l = ResidNum(k)
           RCa(:,l,StepCount/Steps) = R(:,k)

         end if

       end do

     end do

   end do

   RCaAv = 0.
   MSDav = 0.
   MSDex = 0.

   do i = 1 , Nsteps

     do j = ResidNum(Numa) , ResidNum(Numb)

       DR = RCa(:,j,i) - RCa(:,j,0)
       MSDex(j) = MSDex(j) + dot_product( DR, DR )
       RCaAv(:,j) = RCaAv(:,j) + RCa(:,j,i)

     end do

   end do

   MSDex = MSDex / dble(Nsteps)
   MSDex = sqrt( MSDex )

   RCaAv = RCaAv / dble(Nsteps)

   do i = 1 , Nsteps

     do j = ResidNum(Numa) , ResidNum(Numb)

       DR = RCa(:,j,i) - RCaAv(:,j)
       MSDav(j) = MSDav(j) + dot_product( DR, DR )

     end do

   end do

   MSDav = MSDav / dble(Nsteps)
   MSDav = sqrt( MSDav )

   write(1,'(a/)') 'Root square displacements of Ca atoms of the protein'
   write(1,'(a)')  '## Number, Name, MSD from Ex. , MSD from Av., Transmembrane=1'
   write(1,'(a/)') '## data format(i4,x,a4,x,2f10.3,i2)'

   do i = ResidNum(Numa) , ResidNum(Numb)

     j = 0
     if(Flag_MT(i)) j = 1

     write(1,'(i4,x,a4,x,2f10.3,i5)') i, RName(i), MSDex(i), MSDav(i), j

   end do

   close(1)

end subroutine MSD_proteinD


!######################################################################
!######################################################################


! *******************************
! **  Mean square displacement **
! *******************************

subroutine MSD_proteinR

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : AtomName, ResidNum
use TimeParam, only : Timeps

implicit none

integer :: i , j, k
real(8), dimension(3,N) :: Rorig
real(8), dimension(3) :: Rcom, DR, dia
real(8) :: MSD, MSD_TM
integer :: NumCA, NumCA_TM, Numa, Numb
logical, dimension(N) :: Flag_MT

open(1,file='./Analy/MSD_proteinR.data',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   call COM_SpecComponent(Rcom,Kcomp)

   do i = Numa, Numb
     R(:,i) = R(:,i) - Rcom
   end do

   call Rot_to_standard(dia,Numa,Numb)

   write(1,'(a,3e14.4)') 'Original Config.', dia

   Rorig = R

   Flag_MT = .False.

   do i = Numa, Numb

     j = ResidNum(i)

     if(ResidFlag(j) == 1) then
       Flag_MT(i) = .True.
     end if

   end do

! ---------------------------------------------------------------------------

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

       call CellTransform

       call COM_SpecComponent(Rcom,Kcomp)

       do k = Numa, Numb
         R(:,k) = R(:,k) - Rcom
       end do

       call Rot_to_standard(dia,Numa,Numb)

       NumCA = 0
       NumCA_TM = 0
       MSD = 0.d0
       MSD_TM = 0.d0

       do k = Numa, Numb

         if(AtomName(k) == 'CA') then

           NumCA = NumCA + 1
           DR    = R(:,k) - Rorig(:,k)
           MSD   = MSD + dot_product( DR, DR )

           if(Flag_MT(k)) then
             NumCA_TM = NumCA_TM + 1
             MSD_TM   = MSD_TM   + dot_product( DR, DR )
           end if

         end if

       end do

       MSD    = MSD    / NumCA
       MSD_TM = MSD_TM / NumCA_TM

       MSD    = sqrt(MSD)
       MSD_TM = sqrt(MSD_TM)

       write(1,'(f12.4,2f10.3,3e14.4)') Timeps, MSD_TM, MSD, dia

     end do

   end do

close(1)

end subroutine MSD_proteinR


!######################################################################
!######################################################################


subroutine Rot_to_standard(dia,Numa,Numb)

use Configuration, only : R
use AtomParam, only : Mass
use Numbers, only : NumAtm

implicit none

real(8), dimension(3,3) :: Inertia, urot, Inv_urot
real(8), dimension(3,NumAtm(1)) :: Rtmp
real(8), dimension(3) :: dia
integer :: Numa, Numb, i

   Inertia = 0.d0

   do i = Numa, Numb
     Inertia(1,1) = Inertia(1,1) + Mass(i) * ( R(2,i) * R(2,i) + R(3,i) * R(3,i) )
     Inertia(2,2) = Inertia(2,2) + Mass(i) * ( R(3,i) * R(3,i) + R(1,i) * R(1,i) )
     Inertia(3,3) = Inertia(3,3) + Mass(i) * ( R(1,i) * R(1,i) + R(2,i) * R(2,i) )
     Inertia(1,2) = Inertia(1,2) - Mass(i) * R(1,i) * R(2,i)
     Inertia(1,3) = Inertia(1,3) - Mass(i) * R(1,i) * R(3,i)
     Inertia(2,3) = Inertia(2,3) - Mass(i) * R(2,i) * R(3,i)
   end do

   Inertia(2,1) = Inertia(1,2)
   Inertia(3,1) = Inertia(1,3)
   Inertia(3,2) = Inertia(2,3)

   call jacobi(Inertia,urot,dia)

   Inv_urot = Transpose( urot )

   do i = Numa, Numb
     Rtmp(:,i) = matmul( Inv_urot, R(:,i) )
   end do

   do i = Numa, Numb
     R(:,i) = Rtmp(:,i)
   end do

end subroutine Rot_to_standard


!######################################################################
!######################################################################


! *******************************
! **  Mean square displacement **
! *******************************

subroutine MSD_proteinR_detail

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : AtomName, ResidName, ResidNum

implicit none

integer :: i , j, k, l, StepCount
real(8), dimension(3,N) :: Rorig
real(8), dimension(3) :: Rcom, DR, dia
real(8), dimension(:), allocatable :: MSDex, MSDav
integer :: Numa, Numb
logical, dimension(N) :: Flag_MT
real(8), dimension(:,:,:), allocatable :: RCa
real(8), dimension(:,:), allocatable :: RCaAv
character(len=4), dimension(:), allocatable :: RName

   open(1,file='./Analy/MSD_proteinRD.data',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   call COM_SpecComponent(Rcom,Kcomp)

   do i = Numa, Numb
     R(:,i) = R(:,i) - Rcom
   end do

   call Rot_to_standard(dia,Numa,Numb)

   Rorig = R

   Flag_MT = .False.

   j = ResidNum(Numb)

   allocate( RCa(3,j,0:Nsnap) )
   allocate( RCaAv(3,j) )
   allocate( MSDex(j) )
   allocate( MSDav(j) )
   allocate( RName(j) )

   do i = Numa, Numb

     if(AtomName(i) == 'CA') then

       j = ResidNum(i)
       RCa(:,j,0) = Rorig(:,i)

       if(ResidFlag(j) == 1) then
         Flag_MT(j) = .True.
       end if

       RName(j) = ResidName(i)

     end if

   end do

! ---------------------------------------------------------------------------

   StepCount = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepCount = StepCount + 1

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

       call CellTransform

       call COM_SpecComponent(Rcom,Kcomp)

       do k = Numa, Numb

         R(:,k) = R(:,k) - Rcom

       end do

       call Rot_to_standard(dia,Numa,Numb)

       do k = Numa, Numb

         if(AtomName(k) == 'CA') then

           l = ResidNum(k)
           RCa(:,l,StepCount) = R(:,k)

         end if

       end do

     end do

   end do

   RCaAv = 0.
   MSDav = 0.
   MSDex = 0.

   do i = 1 , Nsnap

     do j = ResidNum(Numa) , ResidNum(Numb)

       DR = RCa(:,j,i) - RCa(:,j,0)
       MSDex(j) = MSDex(j) + dot_product( DR, DR )
       RCaAv(:,j) = RCaAv(:,j) + RCa(:,j,i)

     end do

   end do

   MSDex = MSDex / dble(Nsnap)
   MSDex = sqrt( MSDex )

   RCaAv = RCaAv / dble(Nsnap)

   do i = 1 , Nsnap

     do j = ResidNum(Numa) , ResidNum(Numb)

       DR = RCa(:,j,i) - RCaAv(:,j)
       MSDav(j) = MSDav(j) + dot_product( DR, DR )

     end do

   end do

   MSDav = MSDav / dble(Nsnap)
   MSDav = sqrt( MSDav )

   write(1,'(a/)') 'Root square displacements of Ca atoms of the protein'
   write(1,'(a)')  '## Number, Name, MSD from Ex. , MSD from Av., Transmembrane=1'
   write(1,'(a/)') '## data format(i4,x,a4,x,2f10.3,i2)'

   do i = ResidNum(Numa) , ResidNum(Numb)
     write(1,'(i4,x,a4,x,2f10.3,i2)') i, RName(i), MSDex, MSDav, Flag_MT(i)
   end do

   close(1)

end subroutine MSD_proteinR_detail


!######################################################################
!######################################################################


! *******************************
! **  Mean square displacement **
! *******************************

subroutine MSD_proteinMinim

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : AtomName, ResidNum
use TimeParam, only : Timeps

implicit none

integer :: i , j, k, Numa, Numb
real(8), dimension(3,N) :: Rorig
real(8), dimension(3) :: Rcom, DR
real(8) :: MSD, MSD_TM
integer :: NumCA, NumCA_TM
logical, dimension(N) :: Flag_MT

   open(1,file='./Analy/MSD_Minim.data',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   call COM_SpecComponent(Rcom,Kcomp)

   do i = Numa , Numb
     R(:,i) = R(:,i) - Rcom
   end do

   Rorig = R

   Flag_MT = .False.

   do i = Numa , Numb

     j = ResidNum(i)

     if(ResidFlag(j) == 1) then
       Flag_MT(i) = .True.
     end if

   end do

! ---------------------------------------------------------------------------

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

       call CellTransform

       call COM_SpecComponent(Rcom,Kcomp)

       do k = Numa , Numb
         R(:,k) = R(:,k) - Rcom
       end do

       call MinimEuler(Rorig,Numa,Numb)

       NumCA = 0
       NumCA_TM = 0
       MSD = 0.d0
       MSD_TM = 0.d0

       do k = Numa , Numb

         if(AtomName(k) == 'CA') then

           NumCA = NumCA + 1
           DR    = R(:,k) - Rorig(:,k)
           MSD   = MSD + dot_product( DR, DR )

           if(Flag_MT(k)) then

             NumCA_TM = NumCA_TM + 1
             MSD_TM   = MSD_TM   + dot_product( DR, DR )

           end if

         end if

       end do

       MSD    = MSD    / NumCA
       MSD_TM = MSD_TM / NumCA_TM

       MSD    = sqrt(MSD)
       MSD_TM = sqrt(MSD_TM)

       write(1,'(f12.4,2f10.3,3e14.4)') Timeps, MSD_TM, MSD

     end do

   end do

   close(1)

end subroutine MSD_proteinMinim



!######################################################################
!######################################################################


! *******************************
! **  Mean square displacement **
! *******************************

subroutine MSD_proteinMinim_detail

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : AtomName, ResidName, ResidNum

implicit none

integer, parameter :: Steps=10
integer :: i , j, k, l, Numa, Numb, StepCount, Nsteps
real(8), dimension(3,N) :: Rorig
real(8), dimension(3) :: Rcom, DR
real(8), dimension(:), allocatable :: MSDex, MSDav
logical, dimension(N) :: Flag_MT
real(8), dimension(:,:,:), allocatable :: RCa
real(8), dimension(:,:), allocatable :: RCaAv
character(len=4), dimension(:), allocatable :: RName
integer :: flg

   open(1,file='./Analy/MSD_Minim_detail.data',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   Nsteps = Nsnap/Steps

   call COM_SpecComponent(Rcom,Kcomp)

   do i = Numa , Numb
     R(:,i) = R(:,i) - Rcom
   end do

   Rorig = R

   Flag_MT = .False.

   j = ResidNum(Numb)

   allocate( RCa(3,j,0:Nsteps) )
   allocate( RCaAv(3,j) )
   allocate( MSDex(j) )
   allocate( MSDav(j) )
   allocate( RName(j) )

   do i = Numa , Numb

     if(AtomName(i) == 'CA') then

       j = ResidNum(i)
       RCa(:,j,0) = Rorig(:,i)

       if(ResidFlag(j) == 1) then
         Flag_MT(j) = .True.
       end if

       RName(j) = ResidName(i)

     end if

   end do

! ---------------------------------------------------------------------------

   StepCount = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepCount = StepCount + 1

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

       if(mod(StepCount,Steps)/=0) cycle

       call CellTransform

       call COM_SpecComponent(Rcom,Kcomp)

       do k = Numa , Numb
         R(:,k) = R(:,k) - Rcom
       end do

       call MinimEuler(Rorig,Numa,Numb)

       do k = Numa , Numb

         if(AtomName(k) == 'CA') then

           l = ResidNum(k)
           RCa(:,l,StepCount/Steps) = R(:,k)

         end if

       end do

     end do

   end do

   RCaAv = 0.
   MSDav = 0.
   MSDex = 0.

   do i = 1 , Nsteps

     do j = ResidNum(Numa) , ResidNum(Numb)

       DR = RCa(:,j,i) - RCa(:,j,0)
       MSDex(j) = MSDex(j) + dot_product( DR, DR )
       RCaAv(:,j) = RCaAv(:,j) + RCa(:,j,i)

     end do

   end do

   MSDex = MSDex / dble(Nsteps)
   MSDex = sqrt( MSDex )

   RCaAv = RCaAv / dble(Nsteps)

   do i = 1 , Nsteps

     do j = ResidNum(Numa) , ResidNum(Numb)

       DR = RCa(:,j,i) - RCaAv(:,j)
       MSDav(j) = MSDav(j) + dot_product( DR, DR )

     end do

   end do

   MSDav = MSDav / dble(Nsteps)
   MSDav = sqrt( MSDav )

   write(1,'(a)') 'Root square displacements of Ca atoms of the protein'
   write(1,'(a/6x,a/)') 'Note: Modes of translation as well as rotation of',&
   &                          'the whole protein were eliminated'
   write(1,'(a)')  '## Number, Name, MSD from Ex. , MSD from Av., Transmembrane=1'
   write(1,'(a/)') '## data format(i4,x,a4,x,2f10.3,i2)'

   do i = ResidNum(Numa) , ResidNum(Numb)

     flg=0
     if(Flag_MT(i)) flg=1
     write(1,'(i4,x,a4,x,2f10.3,i2)') i+7, RName(i), MSDex(i), MSDav(i), flg

   end do

close(1)

end subroutine MSD_proteinMinim_detail


!######################################################################
!######################################################################


subroutine MinimEuler(Rorig,Numa,Numb)

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze

implicit none

integer :: i, Numa, Numb
real(8) :: fb, fc
real(8) :: phi,theta,psi
real(8) :: csph,snph,tnph
real(8) :: csth,snth,tnth
real(8) :: csps,snps,tnps
real(8), dimension(3,N) :: Rorig
real(8), dimension(3,3) :: Rot
real(8), dimension(3,Numb) :: Rtmp
real(8) :: Z1, Z2

! ## phi

   fb = 0.d0
   fc = 0.d0

   do i = Numa, Numb

     fb = fb + (   Rorig(2,i)*R(1,i) - Rorig(1,i)*R(2,i) )
     fc = fc + ( - Rorig(1,i)*R(1,i) - Rorig(2,i)*R(2,i) )

   end do

   tnph = fb / fc
   phi  = atan( tnph )

   snph = sin( phi )
   csph = cos( phi )

   Z1 =   fb * snph + fc * csph
   Z2 = - fb * snph - fc * csph

   if(Z2 < Z1) then

     snph = - snph
     csph = - csph

   end if

   Rot = 0.d0

   Rot(1,1) =   csph
   Rot(1,2) =   snph
   Rot(2,1) = - snph
   Rot(2,2) =   csph
   Rot(3,3) =   1.d0

   do i = Numa, Numb

     Rtmp(:,i) = matmul( Rot, R(:,i) )

   end do

   do i = Numa, Numb

     R(:,i) = Rtmp(:,i)

   end do

! ## theta

   fb = 0.d0
   fc = 0.d0

   do i = Numa, Numb

     fb = fb + (   Rorig(3,i)*R(2,i) - Rorig(2,i)*R(3,i) )
     fc = fc + ( - Rorig(2,i)*R(2,i) - Rorig(3,i)*R(3,i) )

   end do

   tnth  = fb / fc
   theta = atan( tnth )

   snth = sin( theta )
   csth = cos( theta )

   Z1 =   fb * snth + fc * csth
   Z2 = - fb * snth - fc * csth

   if(Z2 < Z1) then

     snth = - snth
     csth = - csth

   end if

   Rot = 0.d0

   Rot(1,1) =   1.d0
   Rot(2,2) =   csth
   Rot(2,3) =   snth
   Rot(3,2) = - snth
   Rot(3,3) =   csth

   do i = Numa, Numb

     Rtmp(:,i) = matmul( Rot, R(:,i) )

   end do

   do i = Numa, Numb

     R(:,i) = Rtmp(:,i)

   end do

! ## psi

   fb = 0.d0
   fc = 0.d0

   do i = Numa, Numb

     fb = fb + (   Rorig(2,i)*R(1,i) - Rorig(1,i)*R(2,i) )
     fc = fc + ( - Rorig(1,i)*R(1,i) - Rorig(2,i)*R(2,i) )

   end do

   tnps = fb / fc
   psi  = atan( tnps )

   snps = sin( psi )
   csps = cos( psi )

   Z1 =   fb * snps + fc * csps
   Z2 = - fb * snps - fc * csps


   if(Z2 < Z1) then

     snps = - snps
     csps = - csps

   end if

   Rot = 0.d0

   Rot(1,1) =   csps
   Rot(1,2) =   snps
   Rot(2,1) = - snps
   Rot(2,2) =   csps
   Rot(3,3) =   1.d0

   do i = Numa, Numb

     Rtmp(:,i) = matmul( Rot, R(:,i) )

   end do

   do i = Numa, Numb

     R(:,i) = Rtmp(:,i)

   end do

end subroutine MinimEuler


!######################################################################
!######################################################################


! *******************************
! **  Mean square displacement **
! *******************************

subroutine MSD_MinimP_MS(Rorig)

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : AtomName, ResidName, ResidNum

implicit none

integer, parameter :: Steps=10
integer :: i , j, k, l, Numa, Numb, StepCount
real(8), dimension(3,N) :: Rorig
real(8), dimension(3) :: Rcom, DR
real(8), dimension(:), allocatable :: MSDexMain, MSDavMain
real(8), dimension(:), allocatable :: MSDexSide, MSDavSide
character(len=4), dimension(:), allocatable :: RName
character(len=4) :: NameA
integer, dimension(N) :: ChainSp
integer, dimension(N) :: NumMainAtom, NumSideAtom

open(1,file='./Analy/RMSD_Min_MainChain.data',status='unknown')
open(2,file='./Analy/RMSD_Min_SideChain.data',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   j = ResidNum(Numb)

   ChainSp = 0

   allocate( MSDexMain(j) )
   allocate( MSDavMain(j) )
   allocate( MSDexSide(j) )
   allocate( MSDavSide(j) )
   allocate( RName(j) )

   NumMainAtom = 0
   NumSideAtom = 0

   do i = Numa, Numb

     j = ResidNum(i)

     if(AtomName(i) == 'CA') then
       RName(j) = ResidName(i)
     end if

     NameA = AtomName(k)

     if((NameA(1:1) /= 'H').or.(NameA /= 'O').or.(NameA(1:2) /= 'OT')) then

       if((NameA == 'CA').or.(NameA == 'N').or.(NameA == 'C' )) then

         ChainSp(i) = 1
         NumMainAtom(j) = NumMainAtom(j) + 1

       else

         ChainSp(i) = 2
         NumSideAtom(j) = NumSideAtom(j) + 1

       end if

     end if

   end do

! ---------------------------------------------------------------------------

   StepCount = 0
   MSDavMain = 0.d0
   MSDexMain = 0.d0
   MSDavSide = 0.d0
   MSDexSide = 0.d0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepCount = StepCount + 1

!     --------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     --------------------

       call CellTransform

       call COM_SpecComponent(Rcom,Kcomp)

       do k = Numa, Numb

         R(:,k) = R(:,k) - Rcom

       end do

       call MinimEuler(Rorig,Numa,Numb)

       do k = Numa, Numb

         l = ResidNum(k)

         if(ChainSp(k) == 1) then

           DR = R(:,k) - Rorig(:,k)
           MSDexMain(l) = MSDexMain(l) + dot_product( DR, DR )
           DR = R(:,k) - Rav(:,k)
           MSDavMain(l) = MSDavMain(l) + dot_product( DR, DR )

         else if(ChainSp(k) == 2) then

           DR = R(:,k) - Rorig(:,k)
           MSDexSide(l) = MSDexSide(l) + dot_product( DR, DR )
           DR = R(:,k) - Rav(:,k)
           MSDavSide(l) = MSDavSide(l) + dot_product( DR, DR )

         end if

       end do

     end do

   end do

   do i = ResidNum(Numa) , ResidNum(Numb)

     MSDexMain(i) = MSDexMain(i) / dble(StepCount*NumMainAtom(i))
     MSDavMain(i) = MSDavMain(i) / dble(StepCount*NumMainAtom(i))
     MSDexSide(i) = MSDexSide(i) / dble(StepCount*NumSideAtom(i))
     MSDavSide(i) = MSDavSide(i) / dble(StepCount*NumSideAtom(i))

   end do

   MSDexMain = sqrt( MSDexMain )
   MSDavMain = sqrt( MSDavMain )
   MSDexSide = sqrt( MSDexSide )
   MSDavSide = sqrt( MSDavSide )

! ########

   write(1,'(a)') 'Root square displacements of MAIN chain atoms of the protein'
   write(1,'(a)') ' only atoms, N, CA, and C are considered not including H, O or OT '
   write(1,'(a/6x,a/)') 'Note: Modes of translation as well as rotation of',&
   &                          'the whole protein were eliminated'
   write(1,'(a)')  '## Number, Name, MSD from Ex. , MSD from Av.'
   write(1,'(a/)') '## data format(i4,x,a4,x,2f10.3)'

   do i = ResidNum(Numa) , ResidNum(Numb)

     write(1,'(i4,x,a4,x,2f10.3)') i, RName(i), MSDexMain(i), MSDavMain(i)

   end do

   write(2,'(a)') 'Root square displacements of SIDE chain atoms of the protein'
   write(2,'(a)') ' only heavy atoms are considered not including H '
   write(2,'(a/6x,a/)') 'Note: Modes of translation as well as rotation of',&
   &                          'the whole protein were eliminated'
   write(2,'(a)')  '## Number, Name, MSD from Ex. , MSD from Av.'
   write(2,'(a/)') '## data format(i4,x,a4,x,2f10.3)'

   do i = ResidNum(Numa) , ResidNum(Numb)

     write(2,'(i4,x,a4,x,2f10.3)') i, RName(i), MSDexSide(i), MSDavSide(i)

   end do

close(1)
close(2)

end subroutine MSD_MinimP_MS
