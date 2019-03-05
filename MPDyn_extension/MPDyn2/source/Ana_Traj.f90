! ############################
! ## SUBROUTINE LIST 
! ## -- TraceG 
! ## -- ResidueTrajectroy 
! ## -- MakeCRD 
! ## -- Write_CRDa 
! ## -- MakePDB 
! ## -- Write_PDBan 
! ## -- Write_PDBsp 
! ## -- Write_XYZ 
! ## -- Write_DCD 
! ## -- Write_ARC 
! ## -- MakePDBmin 
! ## -- MinimEulerSelected 
! ############################


!######################################################################
!######################################################################


subroutine TraceG

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use UnitExParam, only : Avogadro
use CellParam, only : H, InvH
use AtomParam, only : Mass

implicit none

integer :: i, j, k, l, ll, ii, kk, Nm, Nc, Numa
integer :: TotalStepNumber, Ni, Ns
real(8), dimension(:,:,:), allocatable :: Rcom, Scom, HH
real(8), dimension(3) :: Rg, Sij, Rij
real(8) :: MassM

   open(1,file='./Analy/TraceG.dat',status='unknown')

   Numa = 0

   if(Kcomp/=1) then
     do i = 1 , Kcomp - 1
       Numa = Numa + NumMol(i)*NumAtm(i)
     end do
   end if

   Nm = NumMol(Kcomp)
   Nc = NumAtm(Kcomp)
   Ni = Nsnap

   MassM = 0.d0
   do i = 1 , Nc
     MassM = MassM + Mass(Numa+i)
   end do

   print *, 'Mass(Mol)=',MassM*Avogadro*1.d+3

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)
   end do

   Ns = TotalStepNumber / Ni

   allocate( Rcom(3,Nm,Ns) )
   allocate( Scom(3,Nm,Ns) )
   allocate( HH(3,3,Ns) )

   l  = 0
   ll = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ------------------
       ll = ll + 1

       if( mod(ll,Ni) == 0 ) then

         l = l + 1

!!         call CellTransform
         call Transform

         ii = Numa
         do k = 1 , Nm
           Rg = 0.d0
           do kk = 1 , Nc
             ii = ii + 1
             Rg(:) = Rg(:) + Mass(Numa+kk) * R(:,ii)
           end do
           Rg = Rg / MassM
           Scom(:,k,l) = matmul( InvH, Rg )
         end do

         HH(1,:,l) = H(1,:)
         HH(2,:,l) = H(2,:)
         HH(3,:,l) = H(3,:)

       end if

     end do

   end do

   H(:,1) = HH(:,1,1)
   H(:,2) = HH(:,2,1)
   H(:,3) = HH(:,3,1)

   do i = 1 , Nm
     Rcom(:,i,1) = matmul( H, Scom(:,i,1) )
   end do

   do i = 1 , Ns-1

     j = i + 1
     H(:,1) = HH(:,1,j)
     H(:,2) = HH(:,2,j)
     H(:,3) = HH(:,3,j)

     call InversMatrix(H,InvH)

     do k = 1 , Nm
       Sij = Scom(:,k,j) - Scom(:,k,i)
       Sij = Sij - nint( Sij )
       Rij = matmul( H, Sij )
       Rcom(:,k,j) = Rcom(:,k,i) + Rij
     end do

   end do

   do k = 1 , Nm
     do l = 1 , Ns
       write(1,'(3f10.5)') Rcom(:,k,l)
     end do
     write(1,'(a)') ' = = = '
   end do

   close(1)

end subroutine TraceG


!######################################################################
!######################################################################


subroutine ResidueTrajectroy

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam , only : ResidNum
use AtomParam, only : Mass

implicit none

integer, parameter :: Steps = 50 ! every 1 ps
integer :: i , j, k, l, TotalSteps, Istep
character(len=72) :: Ch
real(8), dimension(3) :: Rg
real(8), dimension(:,:,:), allocatable :: RgRes
real(8) :: SumM
integer :: Numa, Numb

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   TotalSteps = 0
   do i = 1 , NJobs
     TotalSteps = TotalSteps + NTrjStep(i)
   end do

   allocate( RgRes(3,NumResRZ,TotalSteps/Steps) )

   Istep = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       Istep = Istep + 1

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

     if( mod(Istep,Steps) == 0 ) then

       call CellTransform

! ## CHANGE the origin of axes
       call COM_SpecComponent(Rg,Kcomp)

       do k = Numa, Numb
         R(:,k) = R(:,k) - Rg(:)
       end do
! ####

       do l = 1 , NumResRZ
         Rg   = 0.
         SumM = 0.
         do k = Numa, Numb
           if(ResidNum(k) == ResRZ(l)) then
             Rg   = Rg + Mass(k) * R(:,k)
             SumM = SumM + Mass(k)
           end if
         end do
         Rg = Rg / SumM
         RgRes(:,l,Istep/Steps) = Rg
       end do

     end if

     end do

   end do

   do i = 1 , NumResRZ

     if(ResRZ(i)<10) then
       write(Ch,'(a,i1,a)') './Analy/ResTrj000',ResRZ(i),'.dat'
     else if(ResRZ(i)<100) then
       write(Ch,'(a,i2,a)') './Analy/ResTrj00' ,ResRZ(i),'.dat'
     else if(ResRZ(i)<1000) then
       write(Ch,'(a,i3,a)') './Analy/ResTrj0'  ,ResRZ(i),'.dat'
     else if(ResRZ(i)<10000) then
       write(Ch,'(a,i4,a)') './Analy/ResTrj'   ,ResRZ(i),'.dat'
     else
       write(*,*) 'too large Residue number'
       call Finalize
     end if

     open(10,file=Ch,status='unknown')

     do j = 1, TotalSteps/Steps
       write(10,'(3f8.3)') RgRes(:,i,j)
     end do

     close(10)

   end do

end subroutine ResidueTrajectroy


!######################################################################
!######################################################################


! ***********************
! **  Make a CRD File  **
! ***********************

subroutine MakeCRD

use Configuration, only : R
use ParamAnalyze

implicit none

integer :: i , j, NumF

   NumF = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ------------------

       if(mod(j,Interval(i))==0) then
         NumF = NumF + 1
         call Write_CRDa( NumF )
       end if

     end do

   end do

end subroutine MakeCRD


!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine Write_CRDa( NumF )

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : MolName, ResidNum, ResidName, DefaultAtomName
use TimeParam, only : Timeps

implicit NONE

integer :: i, NumF
character(len=4) :: ModelName

   ModelName = MolName(1)(1:4)

  open(1,file=ResultFile(NumF),status='unknown')

   write(1,'(a/)') ModelName
   write(1,*) N

   do i = 1 , N
     write(1,'(2i5,2(x,a4),3f10.5,x,a4)')  &
     &     i , ResidNum(i) , ResidName(i) ,  &
     &     DefaultAtomName(i) , R(:,i) , ModelName
   end do
   write(1,'(f12.4)') Timeps

close(1)

end subroutine Write_CRDa


!######################################################################
!######################################################################


! ***********************
! **  Make a PDB File  **
! ***********************

subroutine MakePDB

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : AtomName, DefaultAtomName

implicit none

integer :: i , j, k, NumF
character(len=4) :: TmpA
real(8), dimension(3) :: Rg
integer :: Numa, Numb, nstr
integer, dimension(20) :: icntrl

   if(cWhat=='XYZsp') then
     open(2,file='./Analy/XMolselected.xyz')
     open(3,file='./Analy/XMol.cell.dat')
   end if
   if(cWhat=='DCD') open(2,file='./Analy/traj.dcd',form='unformatted')
   if(cWhat=='ARC_P')  open(2,file='./Analy/TrajProt.arc')
   if(cWhat=='ARC')    open(2,file='./Analy/Traj.arc')

   if((cWhat=='ARC').or.(cWhat=='ARC_P')) then
     write(2,'(a)') '!BIOSYM archive 3'
     write(2,'(a)') 'PBC=OFF'
     Numa = 0
     Numb = NumMol(1) * NumAtm(1)
     if(Kcomp /= 1) then
       do i = 2, Kcomp
         Numa = Numb
         Numb = Numb + NumMol(i)*NumAtm(i)
       end do
     end if
     Numa = Numa + 1
   end if

   if(cWhat=='XYZsp') then
     do i = 1, N
       TmpA = DefaultAtomName(i)
       if(TmpA(1:3)=='CLA') then
         AtomName(i) = 'Cl'
       else if(TmpA(1:1)=='C') then
         AtomName(i) = 'C'
       else if(TmpA(1:1)=='H') then
         AtomName(i) = 'H'
       else if(TmpA(1:1)=='O') then
         AtomName(i) = 'O'
       else if(TmpA(1:1)=='N') then
         AtomName(i) = 'N'
       else if(TmpA(1:2)=='Si') then
         AtomName(i) = 'Si'
       else if(TmpA(1:2)=='Al') then
         AtomName(i) = 'Al'
       else if(TmpA(1:1)=='P') then
         AtomName(i) = 'P'
       else if(TmpA(1:3)=='SOD') then
         AtomName(i) = 'Na'
       end if
     end do
   end if

   if(cWhat=='DCD') then
     Numa = 0
     do i = 1, Ncomp
       Numa = Numa + NumMol(NumComp(i))*NumAtm(NumComp(i))
     end do
     TmpA='CORD'
     icntrl = 0
     nstr = 0
     icntrl(1) = Nsnap  ! number of frames
     icntrl(2) = 1      ! number of steps in previous run
     icntrl(3) = 1      ! frequency of saving
     icntrl(4) = Nsnap  ! total number of steps
     icntrl(8) = Numa*3 ! number of degrees of freedm
     icntrl(10) = 981668463 ! coded time step
     icntrl(11) = 1         ! coded crystallographic group (or zro)
     icntrl(20) = 27        ! CHARMM version number
     write(2) TmpA,icntrl
     write(2) nstr
     write(2) Numa
   end if

   NumF = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ------------------

       call CellTransform

       if( mod(j,Interval(i)) == 0 ) then

         NumF = NumF + 1

         if(cWhat=='PDB')    call Write_PDBan( NumF )

         if(cWhat=='PDBsp')  call Write_PDBsp( NumF )

         if(cWhat=='PDBex')  call Write_PDBex( NumF )

         if(cWhat=='XYZsp')  call Write_XYZ

         if(cWhat=='DCD')    call Write_DCD( Numa )

         if(cWhat=='ARC')    call Write_ARC(Numa,Numb)

         if(cWhat=='ARC_P') then
           call COM_SpecComponent(Rg,Kcomp)
           do k = Numa, Numb
             R(:,k) = R(:,k) - Rg
           end do
           call Write_ARC(Numa,Numb)
         end if

       end if

     end do

   end do

   if(cWhat=='XYZsp') then
     close(2)
     close(3)
   end if

   if((cWhat=='ARC').or.(cWhat=='ARC_P')) then
     write(2,'(a)') 'end'
     close(2)
   end if

   if(cWhat=='DCD') then
     close(2)
   end if

end subroutine MakePDB


!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine Write_PDBan( NumF )

use Numbers, only : N, NumSpec, NumMol
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : MolName, DefaultAtomName, DefaultResidName, ResidNum
use TimeParam, only : Timeps
use CellParam, only : H
use UnitExParam, only : pi
use CommonBlocks, only : QPBC

implicit NONE

integer :: i, NumF
character(len=4) :: NameA, RName
character(len=1), dimension(10), parameter :: &
&   Flpr = (/'A','B','C','D','E','F','G','H','I','J'/)
real(8), parameter :: zero=0.d0
real(8), dimension(3) :: va, vb, vc
real(8) :: LLa, LLb, LLc, Aab, Abc, Aca

open(2,file=ResultFile(NumF))

   write(2,'(a)') &
   & 'TITLE     An Instantaneous Structure by Molecular Dynamics Simulation'

   write(2,'(a)') &
   & 'REMARK   1 The Simulated System is composed with'

   do i = 1 , NumSpec

     write(2,'(a21,i1,a3,a10,i5)') &
     &     'REMARK   1 Component ',i,' : ',MolName(i),NumMol(i)

   end do

   write(2,'(a)')         'REMARK   2 '
   write(2,'(a,f11.3,a)') 'REMARK   2  Time = ',Timeps,' ps'
   write(2,'(a)')         'REMARK   2 '

   if(QPBC) then

     va = H(:,1)
     vb = H(:,2)
     vc = H(:,3)

     LLa = sqrt( dot_product(va,va) )
     LLb = sqrt( dot_product(vb,vb) )
     LLc = sqrt( dot_product(vc,vc) )
     Aab = acos( dot_product(va,vb) / LLa / LLb ) * 180. / pi
     Abc = acos( dot_product(vb,vc) / LLb / LLc ) * 180. / pi
     Aca = acos( dot_product(vc,va) / LLc / LLa ) * 180. / pi

     write(2,'(a,3f9.3,3f7.2)') 'CRYST1',LLa,LLb,LLc,Abc,Aca,Aab

   end if

   do i = 1, N

     NameA = DefaultAtomName(i)
     RName = DefaultResidName(i)

     if(DefaultResidName(i)=='HSD') RName='HIS'

     if(((DefaultResidName(i)=='DPPC').and.(NameA(4:4)/=' ')) .or. &
     &  ((DefaultResidName(i)=='DMPC').and.(NameA(4:4)/=' ')) .or. &
     &  ((DefaultResidName(i)=='PtPC').and.(NameA(4:4)/=' ')) .or. &
     &  ((DefaultResidName(i)=='TEPC').and.(NameA(4:4)/=' ')) .or. &
     &  ((DefaultResidName(i)=='TBPC').and.(NameA(4:4)/=' ')) .or. &
     &  ((DefaultResidName(i)=='PhPC').and.(NameA(4:4)/=' '))) then

       write(2,'(a4,2x,i5,x,a4,x,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
       & 'ATOM',i,DefaultAtomName(i),RName,ResidNum(i),  &
       &  R(:,i),zero,zero,NameA(1:1)

     else

       write(2,'(a4,2x,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
       & 'ATOM',i,DefaultAtomName(i),RName,ResidNum(i), &
       & R(:,i),zero,zero,NameA(1:1)

     end if

   end do

   write(2,'(a)') 'END'

close(2)

end subroutine Write_PDBan


!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine Write_PDBsp( NumF )

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : MolName, DefaultAtomName, ResidName, ResidNum, Mass
use TimeParam, only : Timeps
use CellParam, only : H
use UnitExParam, only : pi
use CommonBlocks, only : QPBC

implicit NONE

integer :: i, j, k, l, ii, NumF
character(len=4) :: NameA
character(len=1), dimension(10), parameter :: &
&   Flpr = (/'A','B','C','D','E','F','G','H','I','J'/)
real(8), parameter :: zero=0.d0
integer, dimension(Ncomp) :: Numa, Numb
real(8), dimension(3) :: Rcom
real(8) :: Mcom
real(8), dimension(3) :: va, vb, vc
real(8) :: LLa, LLb, LLc, Aab, Abc, Aca

open(2,file=ResultFile(NumF))

   do i = 1 , Ncomp
     Numa(i) = 0
     Numb(i) = NumMol(1) * NumAtm(1)
     if( NumComp(i) /= 1 ) then
       do j = 2, NumComp(i)
         Numa(i) = Numb(i)
         Numb(i) = Numb(i) + NumMol(j) * NumAtm(j)
       end do
     end if
   end do

   write(2,'(a)') &
   & 'TITLE     An Instantaneous Structure by Molecular Dynamics Simulation'

   write(2,'(a)') &
   & 'REMARK   1 The Simulated System is composed with'

   do i = 1 , NumSpec

     write(2,'(a21,i1,a3,a10,i5)') &
     &     'REMARK   1 Component ',i,' : ',MolName(i),NumMol(i)

   end do

   write(2,'(a)') &
   & 'REMARK   2 Here selected components are '

   do i = 1 , Ncomp
     write(2,'(a,i4)') &
     &    'REMARK   2 ',NumComp(i)
   end do

   write(2,'(a)')         'REMARK   3 '
   write(2,'(a,f11.3,a)') 'REMARK   3  Time = ',Timeps,' ps'
   write(2,'(a)')         'REMARK   3 '
!
!   Nall = NumMol(1) * NumAtm(1)
!
!   Rcom = 0.d0
!   Mcom = 0.d0
!
!   do i = 1, Nall
!
!     Rcom = Rcom + Mass(i) * R(:,i)
!     Mcom = Mcom + Mass(i)
!
!   end do
!
!   Rcom = Rcom / Mcom
!
!   do i = 1 , N
!
!     R(:,i) = R(:,i) - Rcom
!
!   end do

   if(QPBC) then

     va = H(:,1)
     vb = H(:,2)
     vc = H(:,3)

     LLa = sqrt( dot_product(va,va) )
     LLb = sqrt( dot_product(vb,vb) )
     LLc = sqrt( dot_product(vc,vc) )
     Aab = acos( dot_product(va,vb) / LLa / LLb ) * 180. / pi
     Abc = acos( dot_product(vb,vc) / LLb / LLc ) * 180. / pi
     Aca = acos( dot_product(vc,va) / LLc / LLa ) * 180. / pi

     write(2,'(a,3f9.3,3f7.2)') 'CRYST1',LLa,LLb,LLc,Abc,Aca,Aab

   end if


   if(Qregion) then

     do i = 1 , Ncomp

       do j = 1 , NumMol(NumComp(i))

         Rcom = 0.d0
         Mcom = 0.d0
         l = Numa(i) + (NumMol(NumComp(i))-1)*NumAtm(NumComp(i))

         do k = 1 , NumAtm(NumComp(i))

           ii = l + k
           Rcom = Rcom + Mass(ii) * R(:,ii)
           Mcom = Mcom + Mass(ii)

         end do

         Rcom = Rcom / Mcom

         if((Rcom(1)<Xmax).and.(Rcom(1)>Xmin).and. &
         &  (Rcom(2)<Ymax).and.(Rcom(2)>Ymin).and. &
         &  (Rcom(3)<Zmax).and.(Rcom(3)>Zmin)) then

         do k = 1, NumAtm(NumComp(i))

           ii = l + k
           NameA = DefaultAtomName(k)
           if(ResidName(k)=='HSD') ResidName(k)='HIS'

           if(((ResidName(i)=='DPPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='PtPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='TEPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='TBPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='DMPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='PhPC').and.(NameA(4:4)/=' '))) then

             write(2,'(a4,2x,i5,x,a4,x,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
             & 'ATOM',ii,DefaultAtomName(ii),ResidName(ii),ResidNum(ii),  &
             &  R(:,ii),zero,zero,NameA(1:1)

           else

             write(2,'(a4,2x,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
             & 'ATOM',ii,DefaultAtomName(ii),ResidName(ii),ResidNum(ii), &
             & R(:,ii),zero,zero,NameA(1:1)

           end if

         end do

         end if

       end do

     end do

   else

     do i = 1 , Ncomp

       do j = 1 , NumMol(NumComp(i))

         l = Numa(i) + (NumMol(NumComp(i))-1)*NumAtm(NumComp(i))

         do k = 1, NumAtm(NumComp(i))

           ii = l + k
           NameA = DefaultAtomName(k)
           if(ResidName(k)=='HSD') ResidName(k)='HIS'

           if(((ResidName(i)=='DPPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='DMPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='PtPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='TEPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='TBPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='PhPC').and.(NameA(4:4)/=' '))) then

             write(2,'(a4,2x,i5,x,a4,x,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
             & 'ATOM',ii,DefaultAtomName(ii),ResidName(ii),ResidNum(ii),  &
             &  R(:,ii),zero,zero,NameA(1:1)

           else

             write(2,'(a4,2x,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
             & 'ATOM',ii,DefaultAtomName(ii),ResidName(ii),ResidNum(ii), &
             & R(:,ii),zero,zero,NameA(1:1)

           end if

         end do

       end do

     end do

   end if

   write(2,'(8(a4,2x,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a/))') &
   & 'ATOM',30001,'Xe  ','Xe  ',4366,-40.0,-40.0,-40.0,0.,0.,'X',&
   & 'ATOM',30002,'Xe  ','Xe  ',4367,-40.0,-40.0, 40.0,0.,0.,'X',&
   & 'ATOM',30003,'Xe  ','Xe  ',4368,-40.0, 40.0,-40.0,0.,0.,'X',&
   & 'ATOM',30004,'Xe  ','Xe  ',4369, 40.0,-40.0,-40.0,0.,0.,'X',&
   & 'ATOM',30005,'Xe  ','Xe  ',4370, 40.0, 40.0,-40.0,0.,0.,'X',&
   & 'ATOM',30006,'Xe  ','Xe  ',4371, 40.0,-40.0, 40.0,0.,0.,'X',&
   & 'ATOM',30007,'Xe  ','Xe  ',4372,-40.0, 40.0, 40.0,0.,0.,'X',&
   & 'ATOM',30008,'Xe  ','Xe  ',4373, 40.0, 40.0, 40.0,0.,0.,'X'

close(2)

end subroutine Write_PDBsp


!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine Write_PDBex( NumF )

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H
use AtomParam, only : MolName, DefaultAtomName, ResidName, ResidNum, Mass
use TimeParam, only : Timeps

implicit NONE

integer :: i, j, k, l, ii, ll, NumF, ix, iy, iz
character(len=4) :: NameA
character(len=1), dimension(10), parameter :: &
&   Flpr = (/'A','B','C','D','E','F','G','H','I','J'/)
real(8), parameter :: zero=0.d0
real(8), dimension(3) :: Rcom, Scom, ddr
real(8) :: Mcom
real(8) :: LLa, LLb, LLc, Aab, Abc, Aca


open(2,file=ResultFile(NumF))

   write(2,'(a)') &
   & 'TITLE     An Instantaneous Structure by Molecular Dynamics Simulation'

   write(2,'(a)') &
   & 'REMARK   1 The Simulated System is composed with'

   do i = 1 , NumSpec

     write(2,'(a21,i1,a3,a10,i5)') &
     &     'REMARK   1 Component ',i,' : ',MolName(i),NumMol(i)

   end do

   write(2,'(a)') &
   & 'REMARK   2 Here selected components are '

   do i = 1 , Ncomp
     write(2,'(a,i4)') &
     &    'REMARK   2 ',NumComp(i)
   end do

   write(2,'(a)')         'REMARK   3 '
   write(2,'(a,f11.3,a)') 'REMARK   3  Time = ',Timeps,' ps'
   write(2,'(a)')         'REMARK   3 '
!
   LLa = Xmax - Xmin
   LLb = Ymax - Ymin
   LLc = Zmax - Zmin
   Aab = 90.
   Abc = 90.
   Aca = 90.

   write(2,'(a,3f9.3,3f7.2)') 'CRYST1',LLa,LLb,LLc,Aab,Abc,Aca

!   Nall = NumMol(1) * NumAtm(1)
!
!   Rcom = 0.d0
!   Mcom = 0.d0
!
!   do i = 1, Nall
!
!     Rcom = Rcom + Mass(i) * R(:,i)
!     Mcom = Mcom + Mass(i)
!
!   end do
!
!   Rcom = Rcom / Mcom
!
!   do i = 1 , N
!
!     R(:,i) = R(:,i) - Rcom
!
!   end do

    do i = 1 , NumSpec

      if(i==1) then
        l = 0
      else
        l = l + NumMol(i-1)*NumAtm(i-1)
      end if

      do j = 1 , NumMol(i)

        Rcom = 0.d0
        Mcom = 0.d0

        ll = l + (j-1) * NumAtm(i)

        do k = 1 , NumAtm(i)

          ii = ll + k
          Rcom = Rcom + Mass(ii) * R(:,ii)
          Mcom = Mcom + Mass(ii)

        end do

        Scom = Rcom / Mcom

        do ix = -1, 1
        do iy = -1, 1
        do iz = -1, 1

          ddr(1) = H(1,1) * ix + H(1,2) * iy + H(1,3) * iz
          ddr(2) = H(2,1) * ix + H(2,2) * iy + H(2,3) * iz
          ddr(3) = H(3,1) * ix + H(3,2) * iy + H(3,3) * iz

          Rcom = Scom + ddr

        if((Rcom(1)<Xmax).and.(Rcom(1)>Xmin).and. &
        &  (Rcom(2)<Ymax).and.(Rcom(2)>Ymin).and. &
        &  (Rcom(3)<Zmax).and.(Rcom(3)>Zmin)) then

         ddr(1) = ddr(1) + Xmax
         ddr(2) = ddr(2) + Ymax
         ddr(3) = ddr(3) + Zmax

         do k = 1, NumAtm(i)

           ii = ll + k
           NameA = DefaultAtomName(ii)
           if(ResidName(k)=='HSD') ResidName(k)='HIS'

           if(((ResidName(i)=='DPPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='DMPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='PtPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='TEPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='TBPC').and.(NameA(4:4)/=' ')) .or. &
           &  ((ResidName(i)=='PhPC').and.(NameA(4:4)/=' '))) then

             write(2,'(a4,2x,i5,2x,a3,x,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
             & 'ATOM',ii,DefaultAtomName(ii)(1:3),ResidName(ii),ResidNum(ii),  &
             &  R(:,ii)+ddr(:),zero,zero,NameA(1:1)

           else

             write(2,'(a4,2x,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
             & 'ATOM',ii,DefaultAtomName(ii),ResidName(ii),ResidNum(ii), &
             & R(:,ii)+ddr(:),zero,zero,NameA(1:1)

           end if

         end do

         end if

         end do
         end do
         end do

       end do

     end do

     write(2,'(a6,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)')         &
           & 'HETATM',N+1,'DUM ','DUM ',N+1, -15.0, -15.0, -15.0, &
           & zero,zero,'D'
     write(2,'(a6,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)')         &
           & 'HETATM',N+2,'DUM ','DUM ',N+1, LLa+15.0, -15.0, -15.0, &
           & zero,zero,'D'
     write(2,'(a6,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)')         &
           & 'HETATM',N+3,'DUM ','DUM ',N+1, -15.0, LLb+15.0, -15.0, &
           & zero,zero,'D'
     write(2,'(a6,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)')         &
           & 'HETATM',N+4,'DUM ','DUM ',N+1, -15.0, -15.0, LLc+15.0, &
           & zero,zero,'D'
     write(2,'(a6,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)')         &
           & 'HETATM',N+5,'DUM ','DUM ',N+1, LLa+15.0, LLb+15.0, -15.0, &
           & zero,zero,'D'
     write(2,'(a6,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)')         &
           & 'HETATM',N+6,'DUM ','DUM ',N+1, -15.0, LLb+15.0, LLc+15.0, &
           & zero,zero,'D'
     write(2,'(a6,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)')         &
           & 'HETATM',N+7,'DUM ','DUM ',N+1, LLa+15.0, -15.0, LLc+15.0, &
           & zero,zero,'D'
     write(2,'(a6,i5,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)')         &
           & 'HETATM',N+8,'DUM ','DUM ',N+1, LLa+15.0, LLb+15.0, LLc+15.0, &
           & zero,zero,'D'

close(2)

end subroutine Write_PDBex



!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine Write_XYZ

use CommonBlocks, only : QPBC
use Configuration, only : R
use CellParam, only : H
use ParamAnalyze
use AtomParam, only : AtomName, Mass
use TimeParam, only : Timeps
use Numbers, only : N, NumMol, NumAtm
use UnitExParam, only : pi

implicit NONE

integer :: i, j, k, l, ii
integer :: Nall
integer, dimension(Ncomp) :: Numa, Numb
logical, dimension(N) :: wflag
real(8), dimension(3) :: Rcom
real(8) :: Mcom
real(8), dimension(3) :: va, vb, vc
real(8) :: LLa, LLb, LLc, Aab, Abc, Aca

   do i = 1 , Ncomp
     if( NumComp(i) == 1 ) then
       Numa(i) = 0
       Numb(i) = NumMol(1) * NumAtm(1)
     else
       Numb(i) = 0
       do j = 1, NumComp(i)
         Numb(i) = Numb(i) + NumMol(j) * NumAtm(j)
       end do
       Numa(i) = Numb(i) - NumMol(NumComp(i)) * NumAtm(NumComp(i))
     end if
   end do

   wflag = .False.

   if(Qregion) then

     do i = 1 , Ncomp

       do j = 1 , NumMol(NumComp(i))

         Rcom = 0.d0
         Mcom = 0.d0
         l = Numa(i) + (j-1)*NumAtm(NumComp(i))

         do k = 1 , NumAtm(NumComp(i))

           ii = l + k
           Rcom = Rcom + Mass(ii) * R(:,ii)
           Mcom = Mcom + Mass(ii)

         end do

         Rcom = Rcom / Mcom

         if((Rcom(1)<Xmax).and.(Rcom(1)>Xmin).and. &
         &  (Rcom(2)<Ymax).and.(Rcom(2)>Ymin).and. &
         &  (Rcom(3)<Zmax).and.(Rcom(3)>Zmin)) then

         do k = 1, NumAtm(NumComp(i))

           ii = l + k
           wflag(ii) = .True.

         end do

         end if

       end do

     end do

   else

     do i = 1 , Ncomp

       do j = Numa(i)+1, Numb(i)

         wflag(j) = .True.

       end do

     end do

   end if

   Nall = 0

   do i = 1, N

     if(wflag(i)) then
       Nall = Nall + 1
     end if

   end do

   write(2,'(i6)') Nall
   write(2,'(a,f11.3,a)') 'Time = ',Timeps,' ps'

   if(N>100000) then
     do i = 1 , N
       if(wflag(i)) then
         write(2,'(a,3f11.4)') AtomName(i),R(:,i)
       end if
     end do
   else
     do i = 1 , N
       if(wflag(i)) then
         write(2,'(a,3f10.4)') AtomName(i),R(:,i)
       end if
     end do
   end if

   if(QPBC) then

     va = H(:,1)
     vb = H(:,2)
     vc = H(:,3)

     LLa = sqrt( dot_product(va,va) )
     LLb = sqrt( dot_product(vb,vb) )
     LLc = sqrt( dot_product(vc,vc) )
     Aab = acos( dot_product(va,vb) / LLa / LLb ) * 180. / pi
     Abc = acos( dot_product(vb,vc) / LLb / LLc ) * 180. / pi
     Aca = acos( dot_product(vc,va) / LLc / LLa ) * 180. / pi

     write(3,'(3f11.5,3f7.2)') LLa,LLb,LLc,Abc,Aca,Aab

   end if

end subroutine Write_XYZ



!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine Write_DCD( Natom )

use CommonBlocks, only : QPBC
use Configuration, only : R
use CellParam, only : H
use ParamAnalyze
use Numbers, only : N, NumMol, NumAtm

implicit NONE

integer :: i, j, k, ii, nn, Natom
real(8), dimension(6) :: xcell
real, dimension(Natom) :: X, Y, Z
integer, dimension(Ncomp) :: Numa, Numb

   xcell(1) = H(1,1)
   xcell(2) = 0.d0
   xcell(3) = H(2,2)
   xcell(4) = 0.d0
   xcell(5) = 0.d0
   xcell(6) = H(3,3)

   do i = 1 , Ncomp
     if( NumComp(i) == 1 ) then
       Numa(i) = 0
       Numb(i) = NumMol(1) * NumAtm(1)
     else
       Numb(i) = 0
       do j = 1, NumComp(i)
         Numb(i) = Numb(i) + NumMol(j) * NumAtm(j)
       end do
       Numa(i) = Numb(i) - NumMol(NumComp(i)) * NumAtm(NumComp(i))
     end if
   end do

   nn = 0
   do i = 1, Ncomp
     ii = Numa(i)
     do j = 1 , NumMol(NumComp(i))
       do k = 1 , NumAtm(NumComp(i))
         ii = ii + 1
         nn = nn + 1
         X(nn) = sngl(R(1,ii))
         Y(nn) = sngl(R(2,ii))
         Z(nn) = sngl(R(3,ii))
       end do
     end do
   end do

   write(2) xcell
   write(2) (X(i),i=1,Natom)
   write(2) (Y(i),i=1,Natom)
   write(2) (Z(i),i=1,Natom)

end subroutine Write_DCD


!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine Write_ARC(Numa,Numb)

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use UnitExParam, only : ec
use NonbondParam, only : Charge
use AtomParam, only : DefaultAtomName, ResidName, ResidNum
use TimeParam, only : Timeps

implicit NONE

integer :: i, numtmp
character(len=4) :: NameA
character(len=5) :: NameB
integer :: Numa, Numb

   write(2,'(a,46x,f11.4)') 'input file for discover',Timeps
   write(2,'(a)') '!DATE Wed Jul 18 23:05:14 2001'

   do i = Numa , Numb

     NameA = DefaultAtomName(i)

     call ChangeAtomName(NameA)

     numtmp = ResidNum(i)

     if(numtmp<10) then
       write(NameB,'(i1)') numtmp
     else if(numtmp<100) then
       write(NameB,'(i2)') numtmp
     else if(numtmp<1000) then
       write(NameB,'(i3)') numtmp
     else if(numtmp<10000) then
       write(NameB,'(i4)') numtmp
     end if

     if(ResidName(i)=='HSD') ResidName(i)='HIS'

     write(2,'(a3,2x,3f15.9,x,a4,x,a5,2x,a1,7x,a1,x,f7.3)')       &
     & NameA(1:3),R(:,i),ResidName(i),NameB,NameA(1:1),&
     & NameA(1:1),Charge(i)/sqrt(ec)

     write(2,'(a)') 'end'

   end do

   if(cWhat=='ARC_P') Return

   do i = 1, Numa - 1

     NameA = DefaultAtomName(i)
     numtmp = ResidNum(i)

     if(numtmp<10) then
       write(NameB,'(i1,3x)') numtmp
     else if(numtmp<100) then
       write(NameB,'(i2,2x)') numtmp
     else if(numtmp<1000) then
       write(NameB,'(i3,1x)') numtmp
     else if(numtmp<10000) then
       write(NameB,'(i4)') numtmp
     end if

     write(2,'(a3,2x,3f15.9,x,a4,x,a4,3x,a1,7x,a1,x,f7.3)')       &
     & NameA(1:3),R(:,i),ResidName(i),NameB,NameA(1:1),&
     & NameA(1:1),Charge(i)/sqrt(ec)

     write(2,'(a)') 'end'

   end do

   do i = Numb + 1, N

     NameA = DefaultAtomName(i)
     numtmp = ResidNum(i)

     if(numtmp<10) then
       write(NameB,'(i1,3x)') numtmp
     else if(numtmp<100) then
       write(NameB,'(i2,2x)') numtmp
     else if(numtmp<1000) then
       write(NameB,'(i3,1x)') numtmp
     else if(numtmp<10000) then
       write(NameB,'(i4)') numtmp
     end if

     write(2,'(a3,2x,3f15.9,x,a4,x,a4,3x,a1,7x,a1,x,f7.3)')       &
     & NameA(1:3),R(:,i),ResidName(i),NameB,NameA(1:1),&
     & NameA(1:1),Charge(i)/sqrt(ec)

     write(2,'(a)') 'end'

   end do

end subroutine Write_ARC


!######################################################################
!######################################################################


! ***********************
! **  Make a PDB File  **
! ***********************

subroutine MakePDBmin

use Numbers, only : N, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use SHAKEparam, only : R_o
use AtomParam, only : AtomName, ResidNum
use TimeParam, only : Timeps

implicit none

integer :: i , j, k, NumF, Nall, NumLS
integer :: flg, itme

character(len=4) :: TmpA

real(8), dimension(3) :: Rg
real(8) :: timett
real :: MSDex, MSDav
integer :: Numa, Numb
logical, dimension(N) :: Flag_Atom,Flag_Resid

   allocate( R_o(3,N) )

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   Flag_Resid = .False.
   Flag_Atom  = .False.

   call COM_SpecComponent(Rg,Kcomp)

   do i = 1 , N

     R(:,i) = R(:,i) - Rg

   end do

   R_o = R

open(1,file='./Analy/MSD_Minim_detail.data',status='old')

   read(1,'(//////)')

   NumLS = 0

   do i = 1, ResidNum(Nall)

     read(1,'(i4,x,a4,x,2f10.3,i2)') j, TmpA, MSDex, MSDav, flg

     if((flg==1).and.(MSDex<2.0)) then

       Flag_Resid(i) = .True.
       NumLS = NumLS + 1

     end if

   end do

   write(*,*) 'NumLS=',NumLS

close(1)

   NumLS = 0

   do i = Numa, Numb

     j = ResidNum(i)

     if((Flag_Resid(j)).and.(AtomName(i)=='CA')) then

       Flag_Atom(i) = .True.
       NumLS = NumLS + 1

     end if

   end do

   write(*,*) 'NumLS=',NumLS

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ------------------

       timett = Timeps/100.
       timett = timett - int(timett)

       if(timett==0.) then

       call CellTransform

       call COM_SpecComponent(Rg,Kcomp)

       do k = 1 , N

         R(:,k) = R(:,k) - Rg

       end do

       call MinimEulerSelected(Flag_Atom,Nall)

       NumF = 1
       itme = int(Timeps)

       if(itme<1000) then

         write(ResultFile(1),'(a,i3,a)') './Analy/snap0',itme,'.pdb'

       else

         write(ResultFile(1),'(a,i4,a)') './Analy/snap',itme,'.pdb'

       end if

       call Write_PDBsp(NumF)

       end if

     end do

   end do

end subroutine MakePDBmin


!######################################################################
!######################################################################


subroutine MinimEulerSelected(Flag_Atom,Nall)

use Numbers, only : N
use Configuration, only : R
use SHAKEparam, only : R_o

implicit none

integer :: i, Nall
real(8) :: fb, fc
real(8) :: phi,theta,psi
real(8) :: csph,snph,tnph
real(8) :: csth,snth,tnth
real(8) :: csps,snps,tnps
real(8), dimension(3,3) :: Rot
real(8), dimension(3,N) :: Rtmp
real(8) :: Z1, Z2
logical, dimension(N) :: Flag_Atom

! ## phi

   fb = 0.d0
   fc = 0.d0
   do i = 1, Nall
     if(Flag_Atom(i)) then
       fb = fb + (   R_o(2,i)*R(1,i) - R_o(1,i)*R(2,i) )
       fc = fc + ( - R_o(1,i)*R(1,i) - R_o(2,i)*R(2,i) )
     end if
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

   do i = 1 , N
     Rtmp(:,i) = matmul( Rot, R(:,i) )
   end do

   R = Rtmp

! ## theta

   fb = 0.d0
   fc = 0.d0
   do i = 1, Nall
     if(Flag_Atom(i)) then
       fb = fb + (   R_o(3,i)*R(2,i) - R_o(2,i)*R(3,i) )
       fc = fc + ( - R_o(2,i)*R(2,i) - R_o(3,i)*R(3,i) )
     end if
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

   do i = 1 , N
     Rtmp(:,i) = matmul( Rot, R(:,i) )
   end do

   R = Rtmp

! ## psi

   fb = 0.d0
   fc = 0.d0
   do i = 1, Nall
     if(Flag_Atom(i)) then
       fb = fb + (   R_o(2,i)*R(1,i) - R_o(1,i)*R(2,i) )
       fc = fc + ( - R_o(1,i)*R(1,i) - R_o(2,i)*R(2,i) )
     end if
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

   do i = 1 , N
     Rtmp(:,i) = matmul( Rot, R(:,i) )
   end do

   R = Rtmp

end subroutine MinimEulerSelected
