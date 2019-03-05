! ############################
! ## SUBROUTINE LIST 
! ## -- dPP 
! ## -- OrientChain 
! ## -- OrientLipidSegments 
! ## -- RotLipid 
! ## -- CoAngleChain 
! ## -- LengthChain 
! ## -- FracGauche_Lipid 
! ## -- CorrTG_LipidChain 
! ## -- SCD 
! ## -- DihedralAngle 
! ## -- AreaOccupL 
! ## -- VoronoiLipid 
! ## -- GrGG_Lipid 
! ## -- Szz_CG
! ############################


!######################################################################
!######################################################################


subroutine dPP

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H
use AtomParam, only : MolName, AtomName
use TimeParam, only : Timeps

implicit none

integer :: i, j, k, ii
integer :: Numa, Numb
real(8), dimension(3) :: a, b, ab
real(8) :: Pupp, Plow, dp, DisPP
integer :: CountUp, CountDn

   open(1,file='./Analy/dPP.dat',status='unknown')

   write(1,'(a)') '# Averaged phosphorus positions along the bilayer normal'
   write(1,'(a)') '# Instantaneous average at t'
   write(1,'(a)') '# z1(t) : the average for the upper monolayer'
   write(1,'(a)') '# z2(t) : the average for the lower monolayer'
   write(1,'(a)') '# d_PP(t) : the distance between above two / membrane thickness'
   write(1,'(a)') '# time[ps]  z1(t)[A]  z2(t)[A]  d_PP(t)[A]'

   ii = 0

   do i = 1, NumSpec

     if((MolName(i)(1:4) == 'DPPC') .or. &
     &  (MolName(i)(1:4) == 'DMPC') .or. &
     &  (MolName(i)(1:4) == 'PtPC') .or. &
     &  (MolName(i)(1:4) == 'TEPC') .or. &
     &  (MolName(i)(1:4) == 'TBPC') .or. &
     &  (MolName(i)(1:4) == 'POPE') .or. &
     &  (MolName(i)(1:4) == 'PhPC') ) then
       ii = ii + 1
       Nlipid = i
     end if

   end do

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

       Pupp = 0.d0
       Plow = 0.d0
       CountUp = 0
       CountDn = 0

       a = H(:,1)
       b = H(:,2)

       call VecProduct(a,b,ab)

       do k = Numa, Numb

         if(AtomName(k)(1:1)=='P') then

           dp = dot_product(ab,R(:,k))

           if( R(3,k) > 0. ) then
             Pupp = Pupp + dp
             CountUp = CountUP + 1
           else
             Plow = Plow - dp
             CountDn = CountDn + 1
           end if

         end if

       end do

       Pupp = Pupp / dble(CountUp)
       Plow = Plow / dble(CountDn)

       disPP = Pupp + Plow

       write(1,'(f10.3,3f8.4)') Timeps, Pupp, Plow, disPP

     end do

   end do

   close(1)

end subroutine dPP


!######################################################################
!######################################################################

! ******************************************************************
! **  DPPC & DPhPC                                                **
! **  C21 -> C215                                                 **
! **  C31 -> C315                                                 **
! **  C21 -> C23,   C23 -> C27,   C27 -> C211,   C211 -> C215     **
! **  C31 -> C33,   C33 -> C37,   C37 -> C311,   C311 -> C315     **
! ******************************************************************

subroutine OrientChain

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use UnitExParam, only : pi
use AtomParam, only : MolName, AtomName

implicit none

integer :: i , j, k, l, ii, jj, ib
integer :: Numa, Numb, nm, nmh, ns, nn, ll, ia
real(8) :: war, t
integer, dimension(:,:,:), allocatable :: ChPair
real(8), dimension(:,:,:,:,:), allocatable :: Chain
integer :: TotalStepNumber, StepNumber
real(8), dimension(3) :: Sij
integer, dimension(180,2,5) :: ic
real(8), dimension(180,2,5) :: cb,cbn
real(8) :: cst, sn, cc, cc2, cn, sumbn
character(len=50), dimension(3,5) :: PairName
character(len=50) :: FileName
integer, parameter :: ni = 10
real(8), dimension(:,:,:), allocatable :: Corr, Corr2
real(8) :: deg

   ii = 0

   do i = 1, NumSpec

     if((MolName(i)(1:4) == 'DPPC') .or. &
     &  (MolName(i)(1:4) == 'PhPC') ) then
       ii = ii + 1
       Nlipid = i
     end if

   end do

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   l = Numa

   nm = NumMol(Nlipid)
   nmh = nm / 2

   allocate( ChPair(2,5,nm) )

   ChPair = 0

   do i = 1 , NumMol(Nlipid)

     do j = 1 , NumAtm(Nlipid)

       l = l + 1

       if(AtomName(l)(1:1) /= 'C') cycle

       if(AtomName(l) == 'C21') then
         ChPair(1,1,i) = l
       else if(AtomName(l) == 'C23') then
         ChPair(1,2,i) = l
       else if(AtomName(l) == 'C27') then
         ChPair(1,3,i) = l
       else if(AtomName(l) == 'C211') then
         ChPair(1,4,i) = l
       else if(AtomName(l) == 'C215') then
         ChPair(1,5,i) = l
       else if(AtomName(l) == 'C31') then
         ChPair(2,1,i) = l
       else if(AtomName(l) == 'C33') then
         ChPair(2,2,i) = l
       else if(AtomName(l) == 'C37') then
         ChPair(2,3,i) = l
       else if(AtomName(l) == 'C311') then
         ChPair(2,4,i) = l
       else if(AtomName(l) == 'C315') then
         ChPair(2,5,i) = l
       end if

     end do

   end do

   PairName(1,1) = 'C21_C215'
   PairName(1,2) = 'C21_C23'
   PairName(1,3) = 'C23_C27'
   PairName(1,4) = 'C27_C211'
   PairName(1,5) = 'C211_C215'
   PairName(2,1) = 'C31_C315'
   PairName(2,2) = 'C31_C33'
   PairName(2,3) = 'C33_C37'
   PairName(2,4) = 'C37_C311'
   PairName(2,5) = 'C311_C315'
   PairName(3,1) = 'Av_C1_C15'
   PairName(3,2) = 'Av_C1_C3'
   PairName(3,3) = 'Av_C3_C7'
   PairName(3,4) = 'Av_C7_C11'
   PairName(3,5) = 'Av_C11_C15'

   do i = 1 , NumMol(Nlipid)
     do j = 1 , 2
       do k = 1 , 5
         if(ChPair(j,k,i)==0) then
           write(*,*) 'ERROR : atom assign'
           call Finalize
         end if
       end do
     end do
   end do

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)
   end do

   Ns = TotalStepNumber / ni

   allocate( Chain(3,2,5,NumMol(Nlipid),TotalStepNumber) )
   allocate( Corr(Ns,2,5) )
   allocate( Corr2(Ns,2,5) )

   Corr  = 0.d0
   Corr2 = 0.d0

   StepNumber = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepNumber = StepNumber + 1

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

!!       call CellTransform
       call Transform

       do k = 1 , NumMol(Nlipid)

         call UnitVect(Sij,ChPair(1,1,k),ChPair(1,5,k))
         Chain(:,1,1,k,StepNumber) = Sij

         call UnitVect(Sij,ChPair(1,1,k),ChPair(1,2,k))
         Chain(:,1,2,k,StepNumber) = Sij

         call UnitVect(Sij,ChPair(1,2,k),ChPair(1,3,k))
         Chain(:,1,3,k,StepNumber) = Sij

         call UnitVect(Sij,ChPair(1,3,k),ChPair(1,4,k))
         Chain(:,1,4,k,StepNumber) = Sij

         call UnitVect(Sij,ChPair(1,4,k),ChPair(1,5,k))
         Chain(:,1,5,k,StepNumber) = Sij

         call UnitVect(Sij,ChPair(2,1,k),ChPair(2,5,k))
         Chain(:,2,1,k,StepNumber) = Sij

         call UnitVect(Sij,ChPair(2,1,k),ChPair(2,2,k))
         Chain(:,2,2,k,StepNumber) = Sij

         call UnitVect(Sij,ChPair(2,2,k),ChPair(2,3,k))
         Chain(:,2,3,k,StepNumber) = Sij

         call UnitVect(Sij,ChPair(2,3,k),ChPair(2,4,k))
         Chain(:,2,4,k,StepNumber) = Sij

         call UnitVect(Sij,ChPair(2,4,k),ChPair(2,5,k))
         Chain(:,2,5,k,StepNumber) = Sij

       end do

     end do

   end do

   ic = 0
   cbn = 0.d0

   do i = 1, TotalStepNumber

     do j = 1 , nmh

       k = j + nmh

       do ii = 1 , 2
         do jj = 1 , 5

           deg = acos(Chain(3,ii,jj,j,i)) * 180. / pi + 1
           ib  = int(deg)
           ic(ib,ii,jj) = ic(ib,ii,jj) + 1
           deg = ib - 0.5
           sn  = sin(deg * pi / 180.d0)
           cbn(ib,ii,jj) = cbn(ib,ii,jj) + 1.d0/sn

           deg = acos(-Chain(3,ii,jj,k,i)) * 180. / pi + 1
           ib  = int(deg)
           ic(ib,ii,jj) = ic(ib,ii,jj) + 1
           deg = ib - 0.5
           sn  = sin(deg * pi / 180.d0)
           cbn(ib,ii,jj) = cbn(ib,ii,jj) + 1.d0/sn

         end do
       end do

     end do

   end do

   do ii = 1 , 2
     do jj = 1 , 5

       write(FileName,'(a,a,a)') './Analy/DisChain_',trim(PairName(ii,jj)),'.dat'

       open(51,file = trim(FileName), status='unknown')

       sumbn = 0.d0
       do i = 1 , 180
         sumbn = sumbn + cbn(i,ii,jj)
       end do

       do i = 1 , 180

         cst = i - 0.5
         cb(i,ii,jj) = dble(ic(i,ii,jj))/dble(nm*TotalStepNumber) * 1.d+02
         cbn(i,ii,jj) = cbn(i,ii,jj)/sumbn * 1.d+02
         write(51,'(f8.4,1x,2f20.10)') cst, cb(i,ii,jj), cbn(i,ii,jj)

       end do

       close(51)

     end do
   end do

   do jj = 1 , 5

     write(FileName,'(a,a,a)') './Analy/DisChain_',trim(PairName(3,jj)),'.dat'

     open(51,file = trim(FileName), status='unknown')

     do i = 1 , 180

       cst = i - 0.5
       cc = (cb(i,1,jj)+cb(i,2,jj)) * 0.5d0
       cn = (cbn(i,1,jj)+cbn(i,2,jj)) * 0.5d0
       write(51,'(f8.4,1x,2f20.10)') cst, cc, cn

     end do

     close(51)

   end do

   do ii = 1, 2
     do jj = 1, 5

       do k = 1 , nm
         do l = 1 , TotalStepNumber - ni

           Sij = Chain(:,ii,jj,k,l)
           nn  = (TotalStepNumber - l) / ni

           do ll = 1, nn
             ia = l + ll*ni
             cc = dot_product( Chain(:,ii,jj,k,ia), Sij )
             Corr(ll,ii,jj)  = Corr(ll,ii,jj)  + cc
             Corr2(ll,ii,jj) = Corr2(ll,ii,jj) + cc * cc
           end do

         end do
       end do

       do l = 1, Ns - 1
         war = 1.d0 / dble(nm * (Ns - l) * ni)
         Corr(l,ii,jj)  = Corr(l,ii,jj) * war
         Corr2(l,ii,jj) = 0.5d0 * ( 3.d0 * Corr2(l,ii,jj) * war - 1.d0 )
       end do

       write(FileName,'(a,a,a)') './Analy/CorrChain_',trim(PairName(ii,jj)),'.dat'

       open(51,file = trim(FileName), status='unknown')

       do l = 2 , Ns - 2, 2
         t = l * dtime * ni
         write(51,'(f7.2,2x,f12.7,2x,f12.7)') t, Corr(l,ii,jj), Corr2(l,ii,jj)
       end do

       close(51)

     end do
   end do

   do jj = 1, 5

     write(FileName,'(a,a,a)') './Analy/CorrChain_',trim(PairName(3,jj)),'.dat'

     open(51,file = trim(FileName), status='unknown')

     do l = 2 , Ns - 2, 2
       t = l * dtime * ni
       cc  = ( Corr(l,1,jj)  + Corr(l,2,jj)  ) * 0.5d0
       cc2 = ( Corr2(l,1,jj) + Corr2(l,2,jj) ) * 0.5d0
       write(51,'(f7.2,2x,f12.7,2x,f12.7)') t, cc, cc2
     end do

     close(51)

   end do

contains

   subroutine UnitVect(Sij,ii,jj)

   integer :: ii , jj
   real(8), dimension(3) :: Rij, Sij
   real(8) :: R2

      Rij = R(:,jj) - R(:,ii)
      R2  = dot_product(Rij,Rij)
      Sij = Rij / sqrt(R2)

   end subroutine UnitVect

end subroutine OrientChain


!######################################################################
!######################################################################

! ******************************************************************
! **  DPPC & DPhPC                                                **
! **  PN,  CC,  CO                                                **
! ******************************************************************

subroutine OrientLipidSegments

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : MolName, AtomName

implicit none

integer :: i , j, k, l, ii
integer :: Numa, Numb, nm
integer, dimension(:,:), allocatable :: ChPair
integer :: StepNumber
real(8), dimension(3) :: Sij
real(8), dimension(:,:), allocatable :: VecPN, VecCCB, VecCCG
real(8), dimension(:,:), allocatable :: VecCOB, VecCOG
real(8), dimension(:), allocatable :: R_PN, R_CCB, R_CCG
real(8) :: R1

   open(1,file='./Analy/PN.dat',form='unformatted',status='unknown')
   open(2,file='./Analy/CC.dat',form='unformatted',status='unknown')
   open(3,file='./Analy/CO.dat',form='unformatted',status='unknown')

   ii = 0

   do i = 1, NumSpec

     if((MolName(i)(1:4) == 'DPPC') .or. &
     &  (MolName(i)(1:4) == 'PhPC') .or. &
     &  (MolName(i)(1:4) == 'DMUG') .or. &
     &  (MolName(i)(1:4) == 'POUG') ) then
       ii = ii + 1
       Nlipid = i
     end if

   end do

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   l = Numa

   nm = NumMol(Nlipid)

   allocate( ChPair(8,nm) )

   allocate( VecPN(3,nm) )
   allocate( VecCCB(3,nm) )
   allocate( VecCCG(3,nm) )
   allocate( VecCOB(3,nm) )
   allocate( VecCOG(3,nm) )
   allocate( R_PN(nm) )
   allocate( R_CCB(nm) )
   allocate( R_CCG(nm) )


   ChPair = 0

   do i = 1 , NumMol(Nlipid)

     do j = 1 , NumAtm(Nlipid)

       l = l + 1

       if((MolName(Nlipid)(1:4) == 'DMUG').or. &
       &  (MolName(Nlipid)(1:4) == 'POUG'))  then
         if(AtomName(l) == 'P') then
           ChPair(1,i) = l
         else if(AtomName(l) == 'O13B') then
           ChPair(2,i) = l
         else if(AtomName(l) == 'C21') then
           ChPair(3,i) = l
         else if(AtomName(l) == 'C214') then
           ChPair(4,i) = l
         else if(AtomName(l) == 'O22') then
           ChPair(5,i) = l
         else if(AtomName(l) == 'C31') then
           ChPair(6,i) = l
         else if(AtomName(l) == 'C314') then
           ChPair(7,i) = l
         else if(AtomName(l) == 'O32') then
           ChPair(8,i) = l
         end if
       else
         if(AtomName(l) == 'P') then
           ChPair(1,i) = l
         else if(AtomName(l) == 'N') then
           ChPair(2,i) = l
         else if(AtomName(l) == 'C21') then
           ChPair(3,i) = l
         else if(AtomName(l) == 'C215') then
           ChPair(4,i) = l
         else if(AtomName(l) == 'O22') then
           ChPair(5,i) = l
         else if(AtomName(l) == 'C31') then
           ChPair(6,i) = l
         else if(AtomName(l) == 'C315') then
           ChPair(7,i) = l
         else if(AtomName(l) == 'O32') then
           ChPair(8,i) = l
         end if
       end if

     end do

   end do

   do i = 1 , NumMol(Nlipid)
     do j = 1 , 8
         if(ChPair(j,i)==0) then
           write(*,*) 'ERROR : atom assign'
           call Finalize
         end if
     end do
   end do

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepNumber = StepNumber + 1

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

!!       call CellTransform
       call Transform

       do k = 1 , NumMol(Nlipid)

         call UnitVect(Sij,R1,ChPair(1,k),ChPair(2,k))
         VecPN(:,k) = Sij
         R_PN(k)    = R1

         call UnitVect(Sij,R1,ChPair(3,k),ChPair(4,k))
         VecCCB(:,k) = Sij
         R_CCB(k)    = R1

         call UnitVect(Sij,R1,ChPair(3,k),ChPair(5,k))
         VecCOB(:,k) = Sij

         call UnitVect(Sij,R1,ChPair(6,k),ChPair(7,k))
         VecCCG(:,k) = Sij
         R_CCG(k)    = R1

         call UnitVect(Sij,R1,ChPair(6,k),ChPair(8,k))
         VecCOG(:,k) = Sij

       end do

       write(1) StepNumber,sngl(VecPN),sngl(R_PN)
       write(2) StepNumber,sngl(VecCCB),sngl(VecCCG),sngl(R_CCB),sngl(R_CCG)
       write(3) StepNumber,sngl(VecCOB),sngl(VecCOG)

     end do

   end do

   close(1)
   close(2)
   close(3)

contains

   subroutine UnitVect(Sij,R1,ii,jj)

   integer :: ii , jj
   real(8), dimension(3) :: Rij, Sij
   real(8) :: R2, R1

      Rij = R(:,jj) - R(:,ii)
      R2  = dot_product(Rij,Rij)
      R1  = sqrt( R2 )
      Sij = Rij / R1

   end subroutine UnitVect

end subroutine OrientLipidSegments


!######################################################################
!######################################################################

! ******************************************************************
! **  DPPC & DPhPC                                                **
! **  Vector from CO to CO in a lipid molecule                    **
! **  Vector between two centers of the chain vector              **
! ******************************************************************

subroutine RotLipid

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : MolName, AtomName

implicit none

integer :: i , j, k, l, ii
integer :: Numa, Numb, nm
integer, dimension(:,:), allocatable :: ChPair
integer :: StepNumber
real(8), dimension(3) :: Sij, VecI, VecJ
real(8), dimension(:,:), allocatable :: VecCCCC, VecCOCO
real(8), dimension(:), allocatable :: R_CCCC, R_COCO
real(8) :: R1

   open(1,file='./Analy/COCO.dat',form='unformatted',status='unknown')
   open(2,file='./Analy/CCCC.dat',form='unformatted',status='unknown')

   ii = 0

   do i = 1, NumSpec

     if((MolName(i)(1:4) == 'DPPC') .or. &
     &  (MolName(i)(1:4) == 'PhPC') ) then
       ii = ii + 1
       Nlipid = i
     end if

   end do

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   l = Numa

   nm = NumMol(Nlipid)

   allocate( ChPair(4,nm) )

   allocate( VecCOCO(3,nm) )
   allocate( VecCCCC(3,nm) )
   allocate( R_COCO(nm) )
   allocate( R_CCCC(nm) )


   ChPair = 0

   do i = 1 , NumMol(Nlipid)

     do j = 1 , NumAtm(Nlipid)

       l = l + 1

       if(AtomName(l) == 'C21') then
         ChPair(1,i) = l
       else if(AtomName(l) == 'C31') then
         ChPair(2,i) = l
       else if(AtomName(l) == 'C215') then
         ChPair(3,i) = l
       else if(AtomName(l) == 'C315') then
         ChPair(4,i) = l
       end if

     end do

   end do

   do i = 1 , NumMol(Nlipid)
     do j = 1 , 4
         if(ChPair(j,i)==0) then
           write(*,*) 'ERROR : atom assign'
           call Finalize
         end if
     end do
   end do

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepNumber = StepNumber + 1

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

!!       call CellTransform
       call Transform

       do k = 1 , NumMol(Nlipid)

         call UnitVect(Sij,R1,ChPair(1,k),ChPair(2,k))
         VecCOCO(:,k) = Sij(:)
         R_COCO(k)    = R1

         call CentVect(VecI,ChPair(1,k),ChPair(3,k))
         call CentVect(VecJ,ChPair(2,k),ChPair(4,k))

         call UnitVec1(Sij,R1,VecI,VecJ)
         VecCCCC(:,k) = Sij(:)
         R_CCCC(k)    = R1

       end do

       write(1) StepNumber,sngl(VecCOCO),sngl(R_COCO)
       write(2) StepNumber,sngl(VecCCCC),sngl(R_CCCC)

     end do

   end do

   close(1)
   close(2)

contains

   subroutine UnitVect(Sij,R1,ii,jj)

   integer :: ii , jj
   real(8), dimension(3) :: Rij, Sij
   real(8) :: R2, R1

      Rij = R(:,jj) - R(:,ii)
      R2  = dot_product(Rij,Rij)
      R1  = sqrt( R2 )
      Sij = Rij / R1

   end subroutine UnitVect

   subroutine UnitVec1(Sij,R1,V1,V2)

   real(8), dimension(3) :: Rij, Sij, V1, V2
   real(8) :: R2, R1

      Rij(:) = V1(:) - V2(:)
      R2  = dot_product(Rij,Rij)
      R1  = sqrt( R2 )
      Sij = Rij / R1

   end subroutine UnitVec1

   subroutine CentVect(Sij,ii,jj)

   integer :: ii , jj
   real(8), dimension(3) :: Rij, Sij

      Rij = R(:,jj) - R(:,ii)
      Rij = Rij * 0.5d0
      Sij = Rij + R(:,ii)

   end subroutine CentVect

end subroutine RotLipid


!######################################################################
!######################################################################

! ******************************************************************
! **  DPPC & DPhPC                                                **
! **  C21 -> C215                                                 **
! **  C31 -> C315                                                 **
! ******************************************************************

subroutine CoAngleChain

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use UnitExParam, only : pi
use AtomParam, only : MolName, AtomName

implicit none

integer :: i , j, k, l, ii
integer :: Numa, Numb, nm, ics
integer, dimension(:,:,:), allocatable :: ChPair
integer :: StepNumber
real(8), dimension(3) :: Sij1, Sij2
integer, dimension(180) :: ic
real(8), dimension(180) :: cb, cbn
real(8) :: cst, sn, th, sumbn
character(len=50) :: FileName

   ii = 0

   do i = 1, NumSpec

     if((MolName(i)(1:4) == 'DPPC') .or. &
     &  (MolName(i)(1:4) == 'PhPC') ) then
       ii = ii + 1
       Nlipid = i
     end if

   end do

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   l = Numa

   nm = NumMol(Nlipid)

   allocate( ChPair(2,2,nm) )

   ChPair = 0

   do i = 1 , NumMol(Nlipid)

     do j = 1 , NumAtm(Nlipid)

       l = l + 1

       if(AtomName(l)(1:1) /= 'C') cycle

       if(AtomName(l) == 'C21') then
         ChPair(1,1,i) = l
       else if(AtomName(l) == 'C215') then
         ChPair(1,2,i) = l
       else if(AtomName(l) == 'C31') then
         ChPair(2,1,i) = l
       else if(AtomName(l) == 'C315') then
         ChPair(2,2,i) = l
       end if

     end do

   end do

   do i = 1 , NumMol(Nlipid)
     do j = 1 , 2
       do k = 1 , 2
         if(ChPair(j,k,i)==0) then
           write(*,*) 'ERROR : atom assign'
           call Finalize
         end if
       end do
     end do
   end do

   StepNumber = 0
   ic = 0
   cbn = 0.d0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepNumber = StepNumber + 1

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

       do k = 1 , NumMol(Nlipid)

         call UnitVect(Sij1,ChPair(1,1,k),ChPair(1,2,k))
         call UnitVect(Sij2,ChPair(2,1,k),ChPair(2,2,k))

         cst = dot_product(Sij1,Sij2)
         th  = acos(cst) * 180.d0 / pi + 1
         ics = int( th )
         ic(ics) = ic(ics) + 1
         th  = ics - 0.5
         sn  = sin(th * pi / 180.d0)
         cbn(ics) = cbn(ics) + 1.d0/sn

       end do

     end do

   end do

   write(FileName,'(a)') './Analy/CoAngChain.dat'

   open(51,file = trim(FileName), status='unknown')

   sumbn = 0.d0
   do i = 1 , 180
     sumbn = sumbn + cbn(i)
   end do

   do i = 1 , 180

     cst = i - 0.5
     cb(i) = dble(ic(i))/dble(nm*StepNumber) * 1.d+02
     cbn(i) = cbn(i) / sumbn * 1.d+02
     write(51,'(f8.4,1x,2f20.10)') cst, cb(i), cbn(i)

   end do

   close(51)

contains

   subroutine UnitVect(Sij,ii,jj)

   integer :: ii , jj
   real(8), dimension(3) :: Rij, Sij
   real(8) :: R2

      Rij = R(:,jj) - R(:,ii)
      R2  = dot_product(Rij,Rij)
      Sij = Rij / sqrt(R2)

   end subroutine UnitVect

end subroutine CoAngleChain


!######################################################################
!######################################################################

! ******************************************************************
! **  DPPC & DPhPC                                                **
! **  C21 -> C215                                                 **
! **  C31 -> C315                                                 **
! ******************************************************************

subroutine LengthChain

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : MolName, AtomName

implicit none

integer :: i , j, k, l, ii, ib
integer :: Numa, Numb, nm, ns, nn, ll, ia
real(8) :: war, t
integer, dimension(:,:,:), allocatable :: ChPair
real(8), dimension(:,:,:), allocatable :: Chain
integer :: TotalStepNumber, StepNumber
real(8) :: R1
integer, dimension(250,2) :: ic
real(8), dimension(250,2) :: cb
real(8) :: cc, cc2
character(len=50), dimension(3) :: PairName
character(len=50) :: FileName
integer, parameter :: ni = 10
real(8), dimension(:,:), allocatable :: Corr, Corr2
real(8), dimension(3) :: AvL, AvL2

   ii = 0

   do i = 1, NumSpec

     if((MolName(i)(1:4) == 'DPPC') .or. &
     &  (MolName(i)(1:4) == 'PhPC') ) then
       ii = ii + 1
       Nlipid = i
     end if

   end do

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   l = Numa

   nm = NumMol(Nlipid)

   allocate( ChPair(2,2,nm) )

   ChPair = 0

   do i = 1 , NumMol(Nlipid)

     do j = 1 , NumAtm(Nlipid)

       l = l + 1

       if(AtomName(l)(1:1) /= 'C') cycle

       if(AtomName(l) == 'C21') then
         ChPair(1,1,i) = l
       else if(AtomName(l) == 'C215') then
         ChPair(1,2,i) = l
       else if(AtomName(l) == 'C31') then
         ChPair(2,1,i) = l
       else if(AtomName(l) == 'C315') then
         ChPair(2,2,i) = l
       end if

     end do

   end do

   PairName(1) = 'C21_C215'
   PairName(2) = 'C31_C315'
   PairName(3) = 'Av_C1_C15'

   do i = 1 , NumMol(Nlipid)
     do j = 1 , 2
       do k = 1 , 2
         if(ChPair(j,k,i)==0) then
           write(*,*) 'ERROR : atom assign'
           call Finalize
         end if
       end do
     end do
   end do

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)
   end do

   Ns = TotalStepNumber / ni

   allocate( Chain(2,NumMol(Nlipid),TotalStepNumber) )
   allocate( Corr(Ns,2) )
   allocate( Corr2(Ns,2) )

   Corr  = 0.d0
   Corr2 = 0.d0

   StepNumber = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepNumber = StepNumber + 1

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

       call CellTransform

       do k = 1 , NumMol(Nlipid)

         call ChainLength(R1,ChPair(1,1,k),ChPair(1,2,k))
         Chain(1,k,StepNumber) = R1

         call ChainLength(R1,ChPair(2,1,k),ChPair(2,2,k))
         Chain(2,k,StepNumber) = R1

       end do

     end do

   end do

   ic = 0
   AvL  = 0.d0
   AvL2 = 0.d0

   do i = 1, TotalStepNumber

     do j = 1 , nm

       do ii = 1 , 2

           ib  = int(Chain(ii,j,i)*10.d0)
           ic(ib,ii) = ic(ib,ii) + 1
           AvL(ii)  = AvL(ii)  + Chain(ii,j,i)
           AvL2(ii) = AvL2(ii) + Chain(ii,j,i) * Chain(ii,j,i)

       end do

     end do

   end do

   do ii = 1, 2

     AvL(ii) = AvL(ii) / dble(nm*TotalStepNumber)

   end do

   AvL(3)  = ( AvL(1) + AvL(2) ) * 0.5d0
   AvL2(3) = ( AvL2(1) + AvL2(2) ) * 0.5d0

   do ii = 1 , 2

     write(FileName,'(a,a,a)') './Analy/DisChainLength_',trim(PairName(ii)),'.dat'

     open(51,file = trim(FileName), status='unknown')

     write(51,'(2(a,f8.2)/)') 'Average Length = ',AvL(ii),' : dispersion = ',AvL2(ii)

     do i = 1 , 250

       cb(i,ii) = dble(ic(i,ii))/dble(nm*TotalStepNumber) * 1.d+02
       write(51,'(f8.4,1x,f20.10)') real(i), cb(i,ii)

     end do

     close(51)

   end do

   write(FileName,'(a,a,a)') './Analy/DisChainLength_',trim(PairName(3)),'.dat'

   open(51,file = trim(FileName), status='unknown')

   write(51,'(2(a,f8.2)/)') 'Average Length = ',AvL(3),' : dispersion = ',AvL2(3)

   do i = 1 , 250

     cc = (cb(i,1)+cb(i,2)) * 0.5d0
     write(51,'(f8.4,1x,f20.10)') real(i), cc

   end do

   close(51)

   do ii = 1, 2

     do k = 1 , nm

       do l = 1 , TotalStepNumber - ni

         R1 = Chain(ii,k,l)
         nn  = (TotalStepNumber - l) / ni

         do ll = 1, nn
           ia = l + ll*ni
           cc = Chain(ii,k,ia) * R1
           Corr(ll,ii)  = Corr(ll,ii)  + cc
         end do

       end do

     end do

     do l = 1, Ns - 1
       war = 1.d0 / dble(nm * (Ns - l) * ni)
       Corr(l,ii)  = Corr(l,ii) * war
     end do

     write(FileName,'(a,a,a)') './Analy/CorrChainLength_',trim(PairName(ii)),'.dat'

     open(51,file = trim(FileName), status='unknown')

     write(51,'(a)') '   0.00   1.0000000'
     do l = 2 , Ns - 2, 2
       t  = l * dtime * ni
       cc = ( Corr(l,ii) - AvL(ii)*AvL(ii) ) / ( AvL2(ii) - AvL(ii) * AvL(ii) )
       write(51,'(f7.2,2x,f12.7)') t, cc
     end do

     close(51)

   end do

   write(FileName,'(a,a,a)') './Analy/CorrChainLength_',trim(PairName(3)),'.dat'

   open(51,file = trim(FileName), status='unknown')

   write(51,'(a)') '   0.00   1.0000000'

   do l = 2 , Ns - 2, 2
     t = l * dtime * ni
     cc  = ( Corr(l,1)  + Corr(l,2)  ) * 0.5d0
     cc2 = ( cc - AvL(3)*AvL(3) ) / ( AvL2(3) - AvL(3) * AvL(3) )
     write(51,'(f7.2,2x,f12.7)') t, cc2
   end do

   close(51)

contains

   subroutine ChainLength(R3,ii,jj)

   integer :: ii , jj
   real(8), dimension(3) :: Rij
   real(8) :: R2, R3

      Rij = R(:,jj) - R(:,ii)
      R2  = dot_product(Rij,Rij)
      R3  = sqrt(R2)

   end subroutine ChainLength

end subroutine LengthChain


!######################################################################
!######################################################################


subroutine FracGauche_Lipid

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : MolName, AtomName

implicit none

integer :: i , j, k, l, id, m
integer :: Numb
integer, dimension(:), allocatable :: Numa
integer, dimension(:), allocatable :: Nslipid
integer, dimension(:), allocatable :: Num_Chain
real(8), parameter :: bp= 120.d0, bm=-120.d0
integer :: NumLip, TotalStep
integer :: MaxChain
integer, dimension(:,:,:,:,:), allocatable :: Dihedrals
integer, dimension(:,:), allocatable :: NumDih_chain
integer, dimension(:,:,:), allocatable :: GaucheNumb
integer, dimension(:,:), allocatable :: Gauche_Chain
real(8), dimension(:,:), allocatable :: FracG
real(8) :: Deg, FracGAv
real(8) :: avt
integer :: ii, jj, kk, ll, istep
character(len=4), dimension(:), allocatable :: LipidName
character(len=72) :: Filename

   do j = 1, 2

     ii = 0

     do i = 1, NumSpec

       if((MolName(i)(1:4) == 'DPPC') .or. &
       &  (MolName(i)(1:4) == 'DMPC') .or. &
       &  (MolName(i)(1:4) == 'PtPC') .or. &
       &  (MolName(i)(1:4) == 'TEPC') .or. &
       &  (MolName(i)(1:4) == 'TBPC') .or. &
       &  (MolName(i)(1:4) == 'POPC') .or. &
       &  (MolName(i)(1:4) == 'DOPC') .or. &
       &  (MolName(i)(1:4) == 'POPE') .or. &
       &  (MolName(i)(1:4) == 'PhPC') .or. &
       &  (MolName(i)(1:3) == 'SDS' ) .or. &
       &  (MolName(i)(1:4) == 'CTAB') .or. &
       &  (MolName(i)(1:4) == 'DTDA') ) then
         ii = ii + 1
         if(j==2) then
           Nslipid(ii) = i
           LipidName(ii) = MolName(i)(1:4)
         end if
       end if

     end do

     if(j==1) then
       if( ii == 0 ) then
         write(*,*) 'ERROR : there is no lipid molecule'
         call Finalize
       end if
       allocate( Numa(ii), Nslipid(ii), LipidName(ii) )
       allocate( Num_Chain(ii) )
       NumLip = ii
     end if

   end do

   do j = 1, NumLip
     write(Filename,'(3a)') './Analy/GaucheFrac_',trim(LipidName(j)),'.dat'
     open(j,file=trim(adjustl(Filename)),status='unknown')
     write(j,'(a)') '# Instantaneous average of the raio of gauche conformers'
     write(j,'(a)') '# in the lipid hydrophobic chains'
     write(j,'(a)') '#'
     write(j,'(a)') '# time[ps] gauche fraction for each chain[%]'
   end do

   Numa(:) = 0

   do j = 1, NumLip
     do i = 1 , Nslipid(j) - 1
       Numa(j) = Numa(j) + NumMol(i) * NumAtm(i)
     end do
   end do

   do j = 1, NumLip
     Num_Chain(j) = 2
     if( LipidName(j) == 'TBPC' ) Num_Chain(j) = 3
     if( LipidName(j) == 'SDS'  ) Num_Chain(j) = 1
     if( LipidName(j) == 'CTAB' ) Num_Chain(j) = 1
   end do

   MaxChain = 0
   do j = 1, NumLip
     if(Num_Chain(j)>MaxChain) MaxChain = Num_Chain(j)
   end do

   allocate( NumDih_Chain(MaxChain,NumLip) )
   allocate( Gauche_Chain(MaxChain,NumLip) )
   allocate( FracG(MaxChain,NumLip) )

   call Dihed_LipidChain

   GaucheNumb = 0
   TotalStep = 0
   istep = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     TotalStep = TotalStep + NTrjStep(i)

     do j = 1 , NTrjStep(i)

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

       istep = istep + 1

       Gauche_Chain = 0

       do m = 1, NumLip

         do k = 1 , NumMol(Nslipid(m))

           do l = 1 , Num_Chain(m)

             do id = 1, NumDih_Chain(l,m)

               ii = Dihedrals(1,id,l,k,m)
               jj = Dihedrals(2,id,l,k,m)
               kk = Dihedrals(3,id,l,k,m)
               ll = Dihedrals(4,id,l,k,m)

               call DihedralAngle( Deg, ii, jj, kk, ll )

               if((Deg > bm).and.(Deg < bp)) then
                 Gauche_Chain(l,m) = Gauche_Chain(l,m)   + 1
                 GaucheNumb(id,l,m) = GaucheNumb(id,l,m) + 1
               end if

             end do

           end do

         end do

         do l = 1 , Num_Chain(m)

           FracG(l,m) = dble(Gauche_Chain(l,m)) / (NumMol(Nslipid(m)) &
           & * NumDih_Chain(l,m)) * 1.d+2

         end do

         FracGAv = 0.
         do l = 1, Num_Chain(m)
         FracGAv = FracGAv + FracG(l,m)
         end do

         FracGAv = FracGAv / dble(Num_Chain(m))

         write(m,'(i10,4f8.4)') istep, (FracG(l,m), l = 1, Num_Chain(m)), FracGAv

       end do

     end do

   end do

   do m = 1, NumLip
     close(m)
   end do

   do m = 1, NumLip
   write(Filename,'(3a)') './Analy/Ave_GaucheFrac_',trim(LipidName(m)),'.dat'
   open(m,file=trim(adjustl(Filename)),status='unknown')

   do i = 1 , Num_Chain(m)

     write(m,'(a,i1,a)') '#### Chain ',i,' ####'

     do j = 1 , NumDih_Chain(i,m)

       avt = dble(GaucheNumb(j,i,m)) / dble(NumMol(Nslipid(m)) * TotalStep) * 1.d+02

       write(m,'(i3,f8.4)') j,avt

     end do

   end do

   close(m)
   end do

!################

contains

   subroutine Dihed_LipidChain

   implicit none

   character(len=4), dimension(:,:,:), allocatable :: AtomLineup
   integer , dimension(:,:,:), allocatable :: NumbLineup
   integer, dimension(:,:), allocatable :: NumCarb
   integer :: Numb1, MaxDih, MDih, MCarb, MaxCarb
   integer :: MaxChain, Mchain, MaxMOL
   integer :: m, NN

    MaxDih = 0
    MaxCarb = 0
    MaxChain = 0
    MaxMOL = 0

    do m = 1, NumLip

      if((LipidName(m) == 'DPPC').or. &
      &  (LipidName(m) == 'PhPC').or. &
      &  (LipidName(m) == 'PtPC')) then

        NumDih_Chain(1,m) = 13
        NumDih_Chain(2,m) = 13
        MDih = 13

      else if(LipidName(m) == 'TEPC') then

        NumDih_Chain(1,m) = 29
        NumDih_Chain(2,m) = 29
        MDih = 29 ! max(NumDih_Chain)

      else if(LipidName(m) == 'TBPC') then

        NumDih_Chain(1,m) = 13
        NumDih_Chain(2,m) = 13
        NumDih_Chain(3,m) = 29
        MDih = 29 ! max(NumDih_Chain)

      else if(LipidName(m) == 'DMPC') then

        NumDih_Chain(1,m) = 11
        NumDih_Chain(2,m) = 11
        MDih = 11 ! max(NumDih_Chain)

      else if((LipidName(m) == 'POPC').or.(LipidName(m) == 'POPE')) then

        NumDih_Chain(1,m) = 14
        NumDih_Chain(2,m) = 13
        MDih = 14 ! max(NumDih_Chain)

      else if(LipidName(m) == 'DOPC') then

        NumDih_Chain(1,m) = 14
        NumDih_Chain(2,m) = 14
        MDih = 14 ! max(NumDih_Chain)

      else if(LipidName(m) == 'SDS') then

        NumDih_Chain(1,m) = 9
        MDih = 9

      else if(LipidName(m) == 'CTAB') then

        NumDih_Chain(1,m) = 12
        MDih = 12

      else if(LipidName(m) == 'DTDA') then

        NumDih_Chain(1,m) = 10
        NumDih_Chain(2,m) = 10
        MDih = 10

      end if

      if(MDih>MaxDih) MaxDih = MDih
      Mchain = Num_Chain(m)
      if(Mchain>MaxChain) MaxChain = Mchain
      MCarb = MDih + 3
      if(MCarb>MaxCarb) MaxCarb = MCarb

    end do

    allocate( NumCarb(MaxChain,NumLip) )
    allocate( AtomLineup(MaxCarb,MaxChain,NumLip) )
    allocate( NumbLineup(MaxCarb,MaxChain,NumLip) )

    do m = 1, NumLip

      if((LipidName(m) == 'DPPC').or. &
      &  (LipidName(m) == 'PhPC').or. &
      &  (LipidName(m) == 'PtPC')) then

        NumCarb(:,m) = NumDih_Chain(:,m) + 3

        AtomLineup( 1,1,m) = 'C21'
        AtomLineup( 2,1,m) = 'C22'
        AtomLineup( 3,1,m) = 'C23'
        AtomLineup( 4,1,m) = 'C24'
        AtomLineup( 5,1,m) = 'C25'
        AtomLineup( 6,1,m) = 'C26'
        AtomLineup( 7,1,m) = 'C27'
        AtomLineup( 8,1,m) = 'C28'
        AtomLineup( 9,1,m) = 'C29'
        AtomLineup(10,1,m) = 'C210'
        AtomLineup(11,1,m) = 'C211'
        AtomLineup(12,1,m) = 'C212'
        AtomLineup(13,1,m) = 'C213'
        AtomLineup(14,1,m) = 'C214'
        AtomLineup(15,1,m) = 'C215'
        AtomLineup(16,1,m) = 'C216'

        AtomLineup( 1,2,m) = 'C31'
        AtomLineup( 2,2,m) = 'C32'
        AtomLineup( 3,2,m) = 'C33'
        AtomLineup( 4,2,m) = 'C34'
        AtomLineup( 5,2,m) = 'C35'
        AtomLineup( 6,2,m) = 'C36'
        AtomLineup( 7,2,m) = 'C37'
        AtomLineup( 8,2,m) = 'C38'
        AtomLineup( 9,2,m) = 'C39'
        AtomLineup(10,2,m) = 'C310'
        AtomLineup(11,2,m) = 'C311'
        AtomLineup(12,2,m) = 'C312'
        AtomLineup(13,2,m) = 'C313'
        AtomLineup(14,2,m) = 'C314'
        AtomLineup(15,2,m) = 'C315'
        AtomLineup(16,2,m) = 'C316'

      else if(LipidName(m) == 'TEPC') then

        NumCarb(:,m) = NumDih_Chain(:,m) + 3

        AtomLineup( 1,1,m) = 'C21'
        AtomLineup( 2,1,m) = 'C22'
        AtomLineup( 3,1,m) = 'C23'
        AtomLineup( 4,1,m) = 'C24'
        AtomLineup( 5,1,m) = 'C25'
        AtomLineup( 6,1,m) = 'C26'
        AtomLineup( 7,1,m) = 'C27'
        AtomLineup( 8,1,m) = 'C28'
        AtomLineup( 9,1,m) = 'C29'
        AtomLineup(10,1,m) = 'C210'
        AtomLineup(11,1,m) = 'C211'
        AtomLineup(12,1,m) = 'C212'
        AtomLineup(13,1,m) = 'C213'
        AtomLineup(14,1,m) = 'C214'
        AtomLineup(15,1,m) = 'C215'
        AtomLineup(16,1,m) = 'C216'
        AtomLineup(17,1,m) = 'C616'
        AtomLineup(18,1,m) = 'C615'
        AtomLineup(19,1,m) = 'C614'
        AtomLineup(20,1,m) = 'C613'
        AtomLineup(21,1,m) = 'C612'
        AtomLineup(22,1,m) = 'C611'
        AtomLineup(23,1,m) = 'C610'
        AtomLineup(24,1,m) = 'C69'
        AtomLineup(25,1,m) = 'C68'
        AtomLineup(26,1,m) = 'C67'
        AtomLineup(27,1,m) = 'C66'
        AtomLineup(28,1,m) = 'C65'
        AtomLineup(29,1,m) = 'C64'
        AtomLineup(30,1,m) = 'C63'
        AtomLineup(31,1,m) = 'C62'
        AtomLineup(32,1,m) = 'C61'

        AtomLineup( 1,2,m) = 'C51'
        AtomLineup( 2,2,m) = 'C52'
        AtomLineup( 3,2,m) = 'C53'
        AtomLineup( 4,2,m) = 'C54'
        AtomLineup( 5,2,m) = 'C55'
        AtomLineup( 6,2,m) = 'C56'
        AtomLineup( 7,2,m) = 'C57'
        AtomLineup( 8,2,m) = 'C58'
        AtomLineup( 9,2,m) = 'C59'
        AtomLineup(10,2,m) = 'C510'
        AtomLineup(11,2,m) = 'C511'
        AtomLineup(12,2,m) = 'C512'
        AtomLineup(13,2,m) = 'C513'
        AtomLineup(14,2,m) = 'C514'
        AtomLineup(15,2,m) = 'C515'
        AtomLineup(16,2,m) = 'C516'
        AtomLineup(17,2,m) = 'C316'
        AtomLineup(18,2,m) = 'C315'
        AtomLineup(19,2,m) = 'C314'
        AtomLineup(20,2,m) = 'C313'
        AtomLineup(21,2,m) = 'C312'
        AtomLineup(22,2,m) = 'C311'
        AtomLineup(23,2,m) = 'C310'
        AtomLineup(24,2,m) = 'C39'
        AtomLineup(25,2,m) = 'C38'
        AtomLineup(26,2,m) = 'C37'
        AtomLineup(27,2,m) = 'C36'
        AtomLineup(28,2,m) = 'C35'
        AtomLineup(29,2,m) = 'C34'
        AtomLineup(30,2,m) = 'C33'
        AtomLineup(31,2,m) = 'C32'
        AtomLineup(32,2,m) = 'C31'

      else if(LipidName(m) == 'TBPC') then

        NumCarb(:,m) = NumDih_Chain(:,m) + 3

        AtomLineup( 1,1,m) = 'C21'
        AtomLineup( 2,1,m) = 'C22'
        AtomLineup( 3,1,m) = 'C23'
        AtomLineup( 4,1,m) = 'C24'
        AtomLineup( 5,1,m) = 'C25'
        AtomLineup( 6,1,m) = 'C26'
        AtomLineup( 7,1,m) = 'C27'
        AtomLineup( 8,1,m) = 'C28'
        AtomLineup( 9,1,m) = 'C29'
        AtomLineup(10,1,m) = 'C210'
        AtomLineup(11,1,m) = 'C211'
        AtomLineup(12,1,m) = 'C212'
        AtomLineup(13,1,m) = 'C213'
        AtomLineup(14,1,m) = 'C214'
        AtomLineup(15,1,m) = 'C215'
        AtomLineup(16,1,m) = 'C216'

        AtomLineup( 1,2,m) = 'C61'
        AtomLineup( 2,2,m) = 'C62'
        AtomLineup( 3,2,m) = 'C63'
        AtomLineup( 4,2,m) = 'C64'
        AtomLineup( 5,2,m) = 'C65'
        AtomLineup( 6,2,m) = 'C66'
        AtomLineup( 7,2,m) = 'C67'
        AtomLineup( 8,2,m) = 'C68'
        AtomLineup( 9,2,m) = 'C69'
        AtomLineup(10,2,m) = 'C610'
        AtomLineup(11,2,m) = 'C611'
        AtomLineup(12,2,m) = 'C612'
        AtomLineup(13,2,m) = 'C613'
        AtomLineup(14,2,m) = 'C614'
        AtomLineup(15,2,m) = 'C615'
        AtomLineup(16,2,m) = 'C616'

        AtomLineup( 1,3,m) = 'C51'
        AtomLineup( 2,3,m) = 'C52'
        AtomLineup( 3,3,m) = 'C53'
        AtomLineup( 4,3,m) = 'C54'
        AtomLineup( 5,3,m) = 'C55'
        AtomLineup( 6,3,m) = 'C56'
        AtomLineup( 7,3,m) = 'C57'
        AtomLineup( 8,3,m) = 'C58'
        AtomLineup( 9,3,m) = 'C59'
        AtomLineup(10,3,m) = 'C510'
        AtomLineup(11,3,m) = 'C511'
        AtomLineup(12,3,m) = 'C512'
        AtomLineup(13,3,m) = 'C513'
        AtomLineup(14,3,m) = 'C514'
        AtomLineup(15,3,m) = 'C515'
        AtomLineup(16,3,m) = 'C516'
        AtomLineup(17,3,m) = 'C316'
        AtomLineup(18,3,m) = 'C315'
        AtomLineup(19,3,m) = 'C314'
        AtomLineup(20,3,m) = 'C313'
        AtomLineup(21,3,m) = 'C312'
        AtomLineup(22,3,m) = 'C311'
        AtomLineup(23,3,m) = 'C310'
        AtomLineup(24,3,m) = 'C39'
        AtomLineup(25,3,m) = 'C38'
        AtomLineup(26,3,m) = 'C37'
        AtomLineup(27,3,m) = 'C36'
        AtomLineup(28,3,m) = 'C35'
        AtomLineup(29,3,m) = 'C34'
        AtomLineup(30,3,m) = 'C33'
        AtomLineup(31,3,m) = 'C32'
        AtomLineup(32,3,m) = 'C31'

      else if(LipidName(m) == 'DMPC') then

        NumCarb(:,m) = NumDih_Chain(:,m) + 3

        AtomLineup( 1,1,m) = 'C21'
        AtomLineup( 2,1,m) = 'C22'
        AtomLineup( 3,1,m) = 'C23'
        AtomLineup( 4,1,m) = 'C24'
        AtomLineup( 5,1,m) = 'C25'
        AtomLineup( 6,1,m) = 'C26'
        AtomLineup( 7,1,m) = 'C27'
        AtomLineup( 8,1,m) = 'C28'
        AtomLineup( 9,1,m) = 'C29'
        AtomLineup(10,1,m) = 'C210'
        AtomLineup(11,1,m) = 'C211'
        AtomLineup(12,1,m) = 'C212'
        AtomLineup(13,1,m) = 'C213'
        AtomLineup(14,1,m) = 'C214'

        AtomLineup( 1,2,m) = 'C31'
        AtomLineup( 2,2,m) = 'C32'
        AtomLineup( 3,2,m) = 'C33'
        AtomLineup( 4,2,m) = 'C34'
        AtomLineup( 5,2,m) = 'C35'
        AtomLineup( 6,2,m) = 'C36'
        AtomLineup( 7,2,m) = 'C37'
        AtomLineup( 8,2,m) = 'C38'
        AtomLineup( 9,2,m) = 'C39'
        AtomLineup(10,2,m) = 'C310'
        AtomLineup(11,2,m) = 'C311'
        AtomLineup(12,2,m) = 'C312'
        AtomLineup(13,2,m) = 'C313'
        AtomLineup(14,2,m) = 'C314'

      else if((LipidName(m) == 'POPC').or.(LipidName(m) == 'POPE')) then

        NumCarb(1,m) = NumDih_Chain(1,m) + 4
        NumCarb(2,m) = NumDih_Chain(2,m) + 3

        AtomLineup( 1,1,m) = 'C21'
        AtomLineup( 2,1,m) = 'C22'
        AtomLineup( 3,1,m) = 'C23'
        AtomLineup( 4,1,m) = 'C24'
        AtomLineup( 5,1,m) = 'C25'
        AtomLineup( 6,1,m) = 'C26'
        AtomLineup( 7,1,m) = 'C27'
        AtomLineup( 8,1,m) = 'C28'
        AtomLineup( 9,1,m) = 'C29'
        AtomLineup(10,1,m) = 'C210'
        AtomLineup(11,1,m) = 'C211'
        AtomLineup(12,1,m) = 'C212'
        AtomLineup(13,1,m) = 'C213'
        AtomLineup(14,1,m) = 'C214'
        AtomLineup(15,1,m) = 'C215'
        AtomLineup(16,1,m) = 'C216'
        AtomLineup(17,1,m) = 'C217'
        AtomLineup(18,1,m) = 'C218'

        AtomLineup( 1,2,m) = 'C31'
        AtomLineup( 2,2,m) = 'C32'
        AtomLineup( 3,2,m) = 'C33'
        AtomLineup( 4,2,m) = 'C34'
        AtomLineup( 5,2,m) = 'C35'
        AtomLineup( 6,2,m) = 'C36'
        AtomLineup( 7,2,m) = 'C37'
        AtomLineup( 8,2,m) = 'C38'
        AtomLineup( 9,2,m) = 'C39'
        AtomLineup(10,2,m) = 'C310'
        AtomLineup(11,2,m) = 'C311'
        AtomLineup(12,2,m) = 'C312'
        AtomLineup(13,2,m) = 'C313'
        AtomLineup(14,2,m) = 'C314'
        AtomLineup(15,2,m) = 'C315'
        AtomLineup(16,2,m) = 'C316'

      else if(LipidName(m) == 'DOPC') then

        NumCarb(1,m) = NumDih_Chain(1,m) + 4
        NumCarb(2,m) = NumDih_Chain(2,m) + 4

        AtomLineup( 1,1,m) = 'C21'
        AtomLineup( 2,1,m) = 'C22'
        AtomLineup( 3,1,m) = 'C23'
        AtomLineup( 4,1,m) = 'C24'
        AtomLineup( 5,1,m) = 'C25'
        AtomLineup( 6,1,m) = 'C26'
        AtomLineup( 7,1,m) = 'C27'
        AtomLineup( 8,1,m) = 'C28'
        AtomLineup( 9,1,m) = 'C29'
        AtomLineup(10,1,m) = 'C210'
        AtomLineup(11,1,m) = 'C211'
        AtomLineup(12,1,m) = 'C212'
        AtomLineup(13,1,m) = 'C213'
        AtomLineup(14,1,m) = 'C214'
        AtomLineup(15,1,m) = 'C215'
        AtomLineup(16,1,m) = 'C216'
        AtomLineup(17,1,m) = 'C217'
        AtomLineup(18,1,m) = 'C218'

        AtomLineup( 1,2,m) = 'C31'
        AtomLineup( 2,2,m) = 'C32'
        AtomLineup( 3,2,m) = 'C33'
        AtomLineup( 4,2,m) = 'C34'
        AtomLineup( 5,2,m) = 'C35'
        AtomLineup( 6,2,m) = 'C36'
        AtomLineup( 7,2,m) = 'C37'
        AtomLineup( 8,2,m) = 'C38'
        AtomLineup( 9,2,m) = 'C39'
        AtomLineup(10,2,m) = 'C310'
        AtomLineup(11,2,m) = 'C311'
        AtomLineup(12,2,m) = 'C312'
        AtomLineup(13,2,m) = 'C313'
        AtomLineup(14,2,m) = 'C314'
        AtomLineup(15,2,m) = 'C315'
        AtomLineup(16,2,m) = 'C316'
        AtomLineup(16,2,m) = 'C317'
        AtomLineup(16,2,m) = 'C318'

      else if(LipidName(m) == 'SDS') then

        NumCarb(:,m) = NumDih_Chain(:,m) + 3

        AtomLineup( 1,1,m) = 'C1'
        AtomLineup( 2,1,m) = 'C2'
        AtomLineup( 3,1,m) = 'C3'
        AtomLineup( 4,1,m) = 'C4'
        AtomLineup( 5,1,m) = 'C5'
        AtomLineup( 6,1,m) = 'C6'
        AtomLineup( 7,1,m) = 'C7'
        AtomLineup( 8,1,m) = 'C8'
        AtomLineup( 9,1,m) = 'C9'
        AtomLineup(10,1,m) = 'C10'
        AtomLineup(11,1,m) = 'C11'
        AtomLineup(12,1,m) = 'C12'

      else if(LipidName(m) == 'CTAB') then

        NumCarb(:,m) = NumDih_Chain(:,m) + 3

        AtomLineup( 1,1,m) = 'C5'
        AtomLineup( 2,1,m) = 'C6'
        AtomLineup( 3,1,m) = 'C7'
        AtomLineup( 4,1,m) = 'C8'
        AtomLineup( 5,1,m) = 'C9'
        AtomLineup( 6,1,m) = 'C10'
        AtomLineup( 7,1,m) = 'C11'
        AtomLineup( 8,1,m) = 'C12'
        AtomLineup( 9,1,m) = 'C13'
        AtomLineup(10,1,m) = 'C14'
        AtomLineup(11,1,m) = 'C15'
        AtomLineup(12,1,m) = 'C16'
        AtomLineup(13,1,m) = 'C17'
        AtomLineup(14,1,m) = 'C18'
        AtomLineup(15,1,m) = 'C19'

      else if(LipidName(m) == 'DTDA') then

        NumCarb(:,m) = NumDih_Chain(:,m) + 3

        AtomLineup( 1,1,m) = 'C5'
        AtomLineup( 2,1,m) = 'C6'
        AtomLineup( 3,1,m) = 'C7'
        AtomLineup( 4,1,m) = 'C8'
        AtomLineup( 5,1,m) = 'C9'
        AtomLineup( 6,1,m) = 'C10'
        AtomLineup( 7,1,m) = 'C11'
        AtomLineup( 8,1,m) = 'C12'
        AtomLineup( 9,1,m) = 'C13'
        AtomLineup(10,1,m) = 'C14'
        AtomLineup(11,1,m) = 'C15'
        AtomLineup(12,1,m) = 'C16'
        AtomLineup(13,1,m) = 'C17'

        AtomLineup( 1,2,m) = 'C18'
        AtomLineup( 2,2,m) = 'C19'
        AtomLineup( 3,2,m) = 'C20'
        AtomLineup( 4,2,m) = 'C21'
        AtomLineup( 5,2,m) = 'C22'
        AtomLineup( 6,2,m) = 'C23'
        AtomLineup( 7,2,m) = 'C24'
        AtomLineup( 8,2,m) = 'C25'
        AtomLineup( 9,2,m) = 'C26'
        AtomLineup(10,2,m) = 'C27'
        AtomLineup(11,2,m) = 'C28'
        AtomLineup(12,2,m) = 'C29'
        AtomLineup(13,2,m) = 'C30'

      end if

    end do

    do m = 1, NumLip
      NN = NumMol(Nslipid(m))
      if(NN>MaxMOL) MaxMOL = NN
    end do

    allocate( Dihedrals(4,MaxDih,MaxChain,MaxMOL,NumLip) )

    do m = 1, NumLip

      do i = 1 , NumMol(Nslipid(m))

        Numb1 = Numa(m) + (i-1)*NumAtm(Nslipid(m))

        do j = 1 , NumAtm(Nslipid(m))

          k = Numb1 + j

          do l = 1, Num_Chain(m)

            do ii = 1 , NumCarb(l,m)

              if(AtomName(k) == AtomLineup(ii,l,m)) then

                NumbLineup(ii,l,m) = k

              end if

            end do

          end do

        end do

        do l = 1, Num_Chain(m)

          do ii = 1, NumDih_Chain(l,m)

            if((LipidName(m) == 'POPC').or.(LipidName(m) == 'POPE')) then
              if((l==1).and.(ii>=8)) then
                jj = ii + 1
              else
                jj = ii
              end if
            else if(LipidName(m) == 'DOPC') then
              if(ii>=8) then
                jj = ii + 1
              else
                jj = ii
              end if
            else
              jj = ii
            end if

            Dihedrals(1,ii,l,i,m) = NumbLineup(jj  ,l,m)
            Dihedrals(2,ii,l,i,m) = NumbLineup(jj+1,l,m)
            Dihedrals(3,ii,l,i,m) = NumbLineup(jj+2,l,m)
            Dihedrals(4,ii,l,i,m) = NumbLineup(jj+3,l,m)

          end do

        end do

      end do

    end do

    allocate( GaucheNumb(MaxDih,MaxChain,NumLip) )

   end subroutine Dihed_LipidChain


end subroutine FracGauche_Lipid


!######################################################################
!######################################################################


subroutine CorrTG_LipidChain

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : MolName, AtomName
use TimeParam, only : Timeps

implicit none

integer :: i , j, k, l, NumDihChain
integer :: ii, jj, kk, ll
integer :: lbb, lrun
integer :: Numa, Numb
integer :: NumF
integer :: TotalStepNumber
integer, dimension(:,:,:), allocatable :: DihedBeta,DihedGamm
real(8) :: DegG, DegB, as, bb, x
character(len=4) :: LipidName
character(len=50) :: FileName
integer, dimension(:,:,:), allocatable :: itb, itg
integer, dimension(:,:), allocatable :: iftb, iftg
real(8), dimension(:,:), allocatable :: ftb, ftg
real(8), dimension(:), allocatable :: savB, savG
real(8), dimension(:), allocatable :: gfB, gfG
real(8), parameter :: zero = 0.d0
real(8), parameter :: one = 1.d0

   ii = 0

   do i = 1, NumSpec

     if((MolName(i)(1:4) == 'DPPC') .or. &
     &  (MolName(i)(1:4) == 'DMPC') .or. &
     &  (MolName(i)(1:4) == 'PhPC') ) then
       LipidName = MolName(i)(1:4)
       ii = ii + 1
       Nlipid = i
     end if

   end do

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   if( ( LipidName == 'DPPC' ) .or. ( LipidName == 'PhPC' ) ) then

     allocate( DihedBeta(4,13,NumMol(Nlipid)) )
     allocate( DihedGamm(4,13,NumMol(Nlipid)) )

   else if( LipidName == 'DMPC' ) then

     allocate( DihedBeta(4,11,NumMol(Nlipid)) )
     allocate( DihedGamm(4,11,NumMol(Nlipid)) )

   end if

   call Dihed_LipidChain

   if(mod(Nsnap,BlockStep) /= 0) then
     write(*,*) 'ERROR : BlockStep'
     write(*,*) Nsnap, BlockStep
     call Finalize
   end if

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)/Interval(i)
     if(i>1) then
       if(Interval(i)/=Interval(i-1)) then
         write(*,*) 'ERROR : Interval'
         write(*,*) 'Interval = ', Interval(i),Interval(i-1)
         call Finalize
       end if
     end if
   end do

   dtime = dtime * Interval(1)

   if(Nsnap /= TotalStepNumber) then
     write(*,*) 'ERROR : Nsnap'
     write(*,*) 'total step number = ',TotalStepNumber
     write(*,*) 'Nsnap             = ',Nsnap
   end if

!   nrun = Nsnap / BlockStep
   lrun = (BlockStep - ilength) / interv + 1
   lbb  = (BlockStep - ilength) + 1

   allocate( itb(NumDihChain,NumMol(Nlipid),BlockStep) )
   allocate( itg(NumDihChain,NumMol(Nlipid),BlockStep) )

   allocate( iftb(ilength,NumDihChain) )
   allocate( iftg(ilength,NumDihChain) )
   allocate( ftb(ilength,NumDihChain) )
   allocate( ftg(ilength,NumDihChain) )

   allocate( savB(NumDihChain) )
   allocate( savG(NumDihChain) )
   allocate( gfB(NumDihChain) )
   allocate( gfG(NumDihChain) )

   do j = 1 , NumDihChain

     if(j<10) then
       write(FileName,'(a,i1,a)') './Analy/Cotg_Beta_Dih_',j,'.dat'
     else
       write(FileName,'(a,i2,a)') './Analy/Cotg_Beta_Dih_',j,'.dat'
     end if

     open(50+j,file=trim(FileName), status='unknown')

     if(j<10) then
       write(FileName,'(a,i1,a)') './Analy/Cotg_Gamm_Dih_',j,'.dat'
     else
       write(FileName,'(a,i2,a)') './Analy/Cotg_Gamm_Dih_',j,'.dat'
     end if

     open(70+j,file=trim(FileName), status='unknown')

   end do

   open(7,file='./Analy/GaucheCheck.dat',status='unknown')

   NumF = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

       if(mod(j,Interval(i))/=0) cycle

       NumF = NumF + 1

       do k = 1 , NumMol(Nlipid)

         do l = 1 , NumDihChain

           ii = DihedBeta(1,l,k)
           jj = DihedBeta(2,l,k)
           kk = DihedBeta(3,l,k)
           ll = DihedBeta(4,l,k)

           call DihedralAngle( DegB, ii, jj, kk, ll )

           ii = DihedGamm(1,l,k)
           jj = DihedGamm(2,l,k)
           kk = DihedGamm(3,l,k)
           ll = DihedGamm(4,l,k)

           call DihedralAngle( DegG, ii, jj, kk, ll )

           if(abs(DegB) > 120.d0) then
             itb(l,k,NumF) =  1
           else
             itb(l,k,NumF) = -1
           end if

           if(abs(DegG) > 120.d0) then
             itg(l,k,NumF) =  1
           else
             itg(l,k,NumF) = -1
           end if

         end do

       end do

       if(NumF == BlockStep) then
         call COTG
         NumF = 0
       end if

     end do

   end do

!-----------------------------------------------------------

Contains

   subroutine COTG

   implicit none

   integer :: ini, ia, iig, iib, iob, iog, NumBeta, NumGamma
   real(8) :: pref

   iftb = 0
   iftg = 0

     do ii = 1 , NumDihChain

       NumBeta = 0
       NumGamma = 0

       do kk = 1, NumMol(Nlipid)

         do ll = 1, BlockStep

           NumBeta  = NumBeta  + itb(ii,kk,ll)
           NumGamma = NumGamma + itg(ii,kk,ll)

         end do

       end do

       savB(ii) = dble(NumBeta)  / dble(NumMol(Nlipid)*BlockStep)
       savG(ii) = dble(NumGamma) / dble(NumMol(Nlipid)*BlockStep)

       gfB(ii)  = 0.5 * (1. - savB(ii)) * 100.
       gfG(ii)  = 0.5 * (1. - savG(ii)) * 100.

     end do

     do jj = 1 , NumDihChain

       do ini = 1 , lbb, interv

         do kk = 1 , NumMol(Nlipid)

           iob = itb(jj,kk,ini)
           iog = itg(jj,kk,ini)

           do ii = 1, ilength-1

             ia  = ini + ii
             iib = iob * itb(jj,kk,ia)
             iig = iog * itg(jj,kk,ia)
             iftb(ii,jj) = iftb(ii,jj) + iib
             iftg(ii,jj) = iftg(ii,jj) + iig

           end do

         end do

       end do

     end do

     Timeps = 0.d0

     pref = 1.d0 / dble(lrun*NumMol(Nlipid))

     do jj = 1 , NumDihChain

       do ll = 1 , ilength - 1

         ftb(ll,jj) = dble(iftb(ll,jj)) * pref
         ftg(ll,jj) = dble(iftg(ll,jj)) * pref

       end do

     end do

     do jj = 1, NumDihChain

       as = savB(jj)*savB(jj)
       bb = 1.d0 - as

       write(jj+50,'(f7.1,f15.10)') zero,one

       do ll = 1 , ilength - 1

         timeps = dtime * ll
         x  = (ftb(ll,jj) - as) / bb
         write(jj+50,'(f7.1,f15.10)') timeps,x

       end do

       write(jj+50,'(a)') ' = = '

     end do

     do jj = 1, NumDihChain

       as = savG(jj)*savG(jj)
       bb = 1.d0 - as

       write(jj+70,'(f7.1,f15.10)') zero,one

       do ll = 1 , ilength - 1

         timeps = dtime * ll
         x  = (ftg(ll,jj) - as) / bb
         write(jj+70,'(f7.1,f15.10)') timeps,x

       end do

       write(jj+70,'(a)') ' = = '

     end do

     do jj = 1 , NumDihChain
       write(7,'(3f7.1)') real(jj),gfB(jj),gfG(jj)
     end do

     write(7,'(a)') ' = = = '

   end subroutine COTG


   subroutine Dihed_LipidChain

   implicit none

   character(len=4), dimension(:), allocatable :: AtomBeta,AtomGamm
   integer , dimension(:), allocatable :: NaBeta, NaGamm
   integer :: NumCarb, Numb1
   logical, dimension(:), allocatable :: FlagBeta,FlagGamm

      if((LipidName == 'DPPC').or.(LipidName == 'PhPC')) then

        NumCarb = 16
        allocate( AtomBeta(16) )
        allocate( AtomGamm(16) )
        allocate( NaBeta(16) )
        allocate( NaGamm(16) )
        allocate( FlagBeta(16) )
        allocate( FlagGamm(16) )

        AtomBeta( 1) = 'C21'
        AtomBeta( 2) = 'C22'
        AtomBeta( 3) = 'C23'
        AtomBeta( 4) = 'C24'
        AtomBeta( 5) = 'C25'
        AtomBeta( 6) = 'C26'
        AtomBeta( 7) = 'C27'
        AtomBeta( 8) = 'C28'
        AtomBeta( 9) = 'C29'
        AtomBeta(10) = 'C210'
        AtomBeta(11) = 'C211'
        AtomBeta(12) = 'C212'
        AtomBeta(13) = 'C213'
        AtomBeta(14) = 'C214'
        AtomBeta(15) = 'C215'
        AtomBeta(16) = 'C216'

        AtomGamm( 1) = 'C31'
        AtomGamm( 2) = 'C32'
        AtomGamm( 3) = 'C33'
        AtomGamm( 4) = 'C34'
        AtomGamm( 5) = 'C35'
        AtomGamm( 6) = 'C36'
        AtomGamm( 7) = 'C37'
        AtomGamm( 8) = 'C38'
        AtomGamm( 9) = 'C39'
        AtomGamm(10) = 'C310'
        AtomGamm(11) = 'C311'
        AtomGamm(12) = 'C312'
        AtomGamm(13) = 'C313'
        AtomGamm(14) = 'C314'
        AtomGamm(15) = 'C315'
        AtomGamm(16) = 'C316'

      else if(LipidName == 'DMPC') then

        NumCarb = 14
        allocate( AtomBeta(14) )
        allocate( AtomGamm(14) )
        allocate( NaBeta(14) )
        allocate( NaGamm(14) )
        allocate( FlagBeta(14) )
        allocate( FlagGamm(14) )

        AtomBeta( 1) = 'C21'
        AtomBeta( 2) = 'C22'
        AtomBeta( 3) = 'C23'
        AtomBeta( 4) = 'C24'
        AtomBeta( 5) = 'C25'
        AtomBeta( 6) = 'C26'
        AtomBeta( 7) = 'C27'
        AtomBeta( 8) = 'C28'
        AtomBeta( 9) = 'C29'
        AtomBeta(10) = 'C210'
        AtomBeta(11) = 'C211'
        AtomBeta(12) = 'C212'
        AtomBeta(13) = 'C213'
        AtomBeta(14) = 'C214'

        AtomGamm( 1) = 'C31'
        AtomGamm( 2) = 'C32'
        AtomGamm( 3) = 'C33'
        AtomGamm( 4) = 'C34'
        AtomGamm( 5) = 'C35'
        AtomGamm( 6) = 'C36'
        AtomGamm( 7) = 'C37'
        AtomGamm( 8) = 'C38'
        AtomGamm( 9) = 'C39'
        AtomGamm(10) = 'C310'
        AtomGamm(11) = 'C311'
        AtomGamm(12) = 'C312'
        AtomGamm(13) = 'C313'
        AtomGamm(14) = 'C314'

      end if

      NumDihChain = NumCarb - 3

      do i = 1 , NumMol(Nlipid)

        Numb1 = Numa + (i-1)*NumAtm(Nlipid)

        FlagBeta = .false.
        FlagGamm = .false.

        do j = 1 , NumAtm(Nlipid)

          k = Numb1 + j

          do ii = 1 , NumCarb

            if(AtomName(k) == AtomBeta(ii)) then

              NaBeta(ii) = k
              FlagBeta(ii) = .true.

            else if(AtomName(k) == AtomGamm(ii)) then

              NaGamm(ii) = k
              FlagGamm(ii) = .true.

            end if

          end do

        end do

        do ii = 1 , NumCarb
          if((.not.FlagBeta(ii)).or.(.not.FlagGamm(ii))) then
            write(*,*) 'ERROR : Carbon number has not been assigned sucessfully.'
            write(*,*) ii
            call Finalize
          end if
        end do

        do j = 1, NumDihChain

          DihedBeta(1,j,i) = NaBeta(j  )
          DihedBeta(2,j,i) = NaBeta(j+1)
          DihedBeta(3,j,i) = NaBeta(j+2)
          DihedBeta(4,j,i) = NaBeta(j+3)

          DihedGamm(1,j,i) = NaGamm(j  )
          DihedGamm(2,j,i) = NaGamm(j+1)
          DihedGamm(3,j,i) = NaGamm(j+2)
          DihedGamm(4,j,i) = NaGamm(j+3)

        end do

      end do

   end subroutine Dihed_LipidChain


end subroutine CorrTG_LipidChain


!######################################################################
!######################################################################


subroutine SCD

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H
use AtomParam, only : MolName, AtomName

implicit none

integer :: ip
integer, dimension(:), allocatable :: Numa
integer, dimension(:), allocatable :: Nslipid
real(8), dimension(3) :: a, b, ab
real(8), dimension(:,:,:), allocatable :: ScdData
integer, dimension(:), allocatable :: Num_Chain
integer :: MaxCarbon, MaxChain, MaxMOL
integer, dimension(:,:), allocatable :: NumCarbon
integer, dimension(:,:,:), allocatable :: NumProton
integer, dimension(:,:,:,:), allocatable :: CarbonNumber
integer, dimension(:,:,:,:,:), allocatable :: ProtonNumber
integer :: NumLip, totalstep
integer :: i, j, k, l, ii, jj, ll, i1, i2
real(8), dimension(3) :: R12
real(8) :: R_12, cs1
character(len=4), dimension(:), allocatable :: LipidName

   open(1,file='./Analy/SCD.dat',status='unknown')

   do j = 1, 2

     ii = 0

     do i = 1, NumSpec

       if((MolName(i)(1:4) == 'DPPC') .or. &
       &  (MolName(i)(1:4) == 'DMPC') .or. &
       &  (MolName(i)(1:4) == 'PtPC') .or. &
       &  (MolName(i)(1:4) == 'TEPC') .or. &
       &  (MolName(i)(1:4) == 'TBPC') .or. &
       &  (MolName(i)(1:4) == 'POPC') .or. &
       &  (MolName(i)(1:4) == 'POPE') .or. &
       &  (MolName(i)(1:4) == 'DOPC') .or. &
       &  (MolName(i)(1:4) == 'PhPC') .or. &
       &  (MolName(i)(1:3) == 'SDS' ) .or. &
       &  (MolName(i)(1:4) == 'CTAB') .or. &
       &  (MolName(i)(1:4) == 'DTDA') .or. &
       &  (MolName(i)(1:4) == 'SPM2') .or. &
       &  (MolName(i)(1:4) == 'SM18')) then
         ii = ii + 1
         if(j==2) then
           Nslipid(ii) = i
           LipidName(ii) = MolName(i)(1:4)
         end if
       end if

     end do

     if(j==1) then
       if( ii == 0 ) then
         write(*,*) 'ERROR : there is no lipid molecule'
         call Finalize
       end if
       allocate( Numa(ii), Nslipid(ii), LipidName(ii) )
       allocate( Num_Chain(ii) )
       NumLip = ii
     end if

   end do

   Numa(:) = 0

   do j = 1, NumLip
     do i = 1 , Nslipid(j) - 1
       Numa(j) = Numa(j) + NumMol(i) * NumAtm(i)
     end do
   end do

   do j = 1, NumLip
     Num_Chain(j) = 2
     if( LipidName(j) == 'TBPC' ) Num_Chain(j) = 3
     if( LipidName(j) == 'SDS'  ) Num_Chain(j) = 1
     if( LipidName(j) == 'CTAB' ) Num_Chain(j) = 1
   end do

   MaxChain = 0
   do j = 1, NumLip
     if(Num_Chain(j)>MaxChain) MaxChain = Num_Chain(j)
   end do

   allocate( NumCarbon(MaxChain,NumLip) )

   MaxCarbon = 0
   do j = 1, NumLip
     if( LipidName(j) == 'DPPC' ) then
       NumCarbon(:,j) = 15
     else if( LipidName(j) == 'PhPC' ) then
       NumCarbon(:,j) = 19
     else if( LipidName(j) == 'PtPC' ) then
       NumCarbon(:,j) = 20
     else if( LipidName(j) == 'DMPC' ) then
       NumCarbon(:,j) = 13
     else if(LipidName(j) == 'TEPC') then
       NumCarbon(:,j) = 40
     else if(LipidName(j) == 'TBPC') then
       NumCarbon(:,j) = 20
       NumCarbon(3,j) = 40
     else if((LipidName(j) == 'POPC').or.(LipidName(j) == 'POPE')) then
       NumCarbon(1,j) = 17
       NumCarbon(2,j) = 15
     else if(LipidName(j) == 'DOPC') then
       NumCarbon(:,j) = 17
     else if(LipidName(j) == 'SDS') then
       NumCarbon(1,j) = 12
     else if(LipidName(j) == 'CTAB') then
       NumCarbon(1,j) = 15
     else if(LipidName(j) == 'DTDA') then
       NumCarbon(:,j) = 13
     else if(LipidName(j) == 'SPM2') then
       NumCarbon(1,j) = 17
       NumCarbon(2,j) = 16
     else if(LipidName(j) == 'SM18') then
       NumCarbon(1,j) = 17
       NumCarbon(2,j) = 16
     end if
   end do

   MaxCarbon = 0
   do j = 1, NumLip
     do i = 1, Num_Chain(j)
       if(NumCarbon(i,j)>MaxCarbon) MaxCarbon = NumCarbon(i,j)
     end do
   end do

   MaxMOL = 0
   do j = 1, NumLip
     i = Nslipid(j)
     if(NumMol(i)>MaxMOL) MaxMOL = NumMol(i)
   end do

   allocate( NumProton(MaxCarbon,MaxChain,NumLip) )
   allocate( CarbonNumber(MaxCarbon,MaxChain,MaxMOL,NumLip) )
   allocate( ProtonNumber(3,MaxCarbon,MaxChain,MaxMOL,NumLip) )

   allocate( ScdData(MaxCarbon,MaxChain,NumLip) )

   call define_SCD

   ScdData = 0.d0
   totalstep = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     totalstep = totalstep + NTrjStep(i)

     do j = 1 , NTrjStep(i)

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

       a  = H(:,1)
       b  = H(:,2)
       call VecProduct(a,b,ab)

       do ip = 1, NumLip

       do k = 1 , NumMol(Nslipid(ip))

         do ii = 1 , Num_Chain(ip)

           do l = 1 , NumCarbon(ii,ip)

             do ll = 1, NumProton(l,ii,ip)

               i1 = CarbonNumber(l,ii,k,ip)
               i2 = ProtonNumber(ll,l,ii,k,ip)

               R12 = R(:,i2) - R(:,i1)

               R_12 = sqrt( dot_product( R12, R12 ) )

               R12 = R12 / R_12

               cs1 = dot_product( R12, ab )

               ScdData(l,ii,ip) = ScdData(l,ii,ip) + cs1*cs1

             end do

           end do

         end do

       end do

       end do

     end do

   end do

   do j = 1, NumLip

   do ii = 1, Num_Chain(j)

     do l = 1, NumCarbon(ii,j)

       ScdData(l,ii,j) = ScdData(l,ii,j) / dble(TotalStep*NumMol(Nslipid(j))*NumProton(l,ii,j))

     end do

   end do

   end do

   ScdData = 0.5d0 * ( 3.d0*ScdData - 1.d0 )

!   write(1,'(a)') 'Scd of the acyl chain of lipid molecules'

   do k = 1, NumLip

     write(1,'(3a)') '### Order Parameter Profile for ',LipidName(k),' ###'

   do i = 1, Num_Chain(k)

     write(1,'(a,i1,a)') '#### Chain ',i,' ####'

     do j = 1, NumCarbon(i,k)

       if( (LipidName(k) == 'DPPC').or.(LipidName(k) == 'DMPC').or.(LipidName(k) == 'PhPC') &
     & .or.(LipidName(k) == 'POPC').or.(LipidName(k) == 'POPE').or.(LipidName(k) == 'DOPC') &
     & .or.((LipidName(k) == 'SM18').and.(i==1)).or.((LipidName(k) == 'SPM2').and.(i==1))) then
         ii = j + 1
       else if(((LipidName(k) == 'SM18').and.(i==2)).or.((LipidName(k) == 'SPM2').and.(i==2))) then
         ii = j + 2
       else
         ii = j
       end if

       write(1,'(i8,3f13.5)') ii, ScdData(j,i,k)

     end do

   end do

     write(1,*)

   end do

   close(1)

!##############################

contains

   subroutine define_SCD

   implicit none

   character(len=4), dimension(MaxCarbon,MaxChain,NumLip) :: CarbonName
   character(len=4), dimension(3,MaxCarbon,MaxChain,NumLip) :: ProtonName
   integer :: Numb1, kk

   do k = 1, NumLip

   if( ( LipidName(k) == 'DPPC' ).or. &
   &   ( LipidName(k) == 'DMPC' ).or. &
   &   ( LipidName(k) == 'PhPC' ) ) then

      do i = 1, 8
        write(CarbonName(i,1,k),  '(a2,i1)')    'C2',i+1
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i+1,'R'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i+1,'S'
        write(ProtonName(3,i,1,k),'(a1,i1,a1)') 'H',i+1,'T'

        write(CarbonName(i,2,k),  '(a2,i1)')    'C3',i+1
        write(ProtonName(1,i,2,k),'(a1,i1,a1)') 'H',i+1,'X'
        write(ProtonName(2,i,2,k),'(a1,i1,a1)') 'H',i+1,'Y'
        write(ProtonName(3,i,2,k),'(a1,i1,a1)') 'H',i+1,'Z'
      end do  

      do i = 9, MaxCarbon
        write(CarbonName(i,1,k),  '(a2,i2)')    'C2',i+1
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i+1,'R'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i+1,'S'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i+1,'T'

        write(CarbonName(i,2,k),  '(a2,i2)')    'C3',i+1
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i+1,'X'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i+1,'Y'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i+1,'Z'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if( (LipidName(k) == 'PhPC') .and. &
          & ( (j == 2).or.(j == 6).or.(j == 10).or.(j == 14) ) ) then
            NumProton(j,i,k) = 1
          else if( ( (LipidName(k) == 'DPPC') .and. (j == 15) ) .or. &
          &        ( (LipidName(k) == 'PhPC') .and. (j >= 15) ) .or. &
          &        ( (LipidName(k) == 'DMPC') .and. (j == 13) ) ) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

    else if(LipidName(k) == 'PtPC') then

      do i = 1, 9
        write(CarbonName(i,1,k),  '(a2,i1)')    'C2',i
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i,'R'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i,'S'
        write(ProtonName(3,i,1,k),'(a1,i1,a1)') 'H',i,'T'

        write(CarbonName(i,2,k),  '(a2,i1)')    'C3',i
        write(ProtonName(1,i,2,k),'(a1,i1,a1)') 'H',i,'X'
        write(ProtonName(2,i,2,k),'(a1,i1,a1)') 'H',i,'Y'
        write(ProtonName(3,i,2,k),'(a1,i1,a1)') 'H',i,'Z'
      end do  

      do i = 10, MaxCarbon
        write(CarbonName(i,1,k),  '(a2,i2)')    'C2',i
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i,'R'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i,'S'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i,'T'

        write(CarbonName(i,2,k),  '(a2,i2)')    'C3',i
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i,'X'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i,'Y'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i,'Z'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if( (j == 3).or.(j == 7).or.(j == 11).or.(j == 15) ) then
            NumProton(j,i,k) = 1
          else if( j >= 16 ) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

   else if(LipidName(k) == 'TEPC') then

      do i = 1, 9
        write(CarbonName(i,1,k),  '(a2,i1)')    'C2',i
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i,'R'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i,'S'
        write(ProtonName(3,i,1,k),'(a1,i1,a1)') 'H',i,'T'

        write(CarbonName(i,2,k),  '(a2,i1)')    'C3',i
        write(ProtonName(1,i,2,k),'(a1,i1,a1)') 'H',i,'X'
        write(ProtonName(2,i,2,k),'(a1,i1,a1)') 'H',i,'Y'
        write(ProtonName(3,i,2,k),'(a1,i1,a1)') 'H',i,'Z'
      end do  

      do i = 10, 20
        write(CarbonName(i,1,k),  '(a2,i2)')    'C2',i
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i,'R'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i,'S'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i,'T'

        write(CarbonName(i,2,k),  '(a2,i2)')    'C3',i
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i,'X'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i,'Y'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i,'Z'
      end do

      ii = 20

      do i = 20, 10, -1
        ii = ii + 1
        write(CarbonName(ii,1,k),  '(a2,i2)')    'C6',i
        write(ProtonName(1,ii,1,k),'(a1,i2,a1)') 'H',i,'U'
        write(ProtonName(2,ii,1,k),'(a1,i2,a1)') 'H',i,'V'
        write(ProtonName(3,ii,1,k),'(a1,i2,a1)') 'H',i,'W'

        write(CarbonName(ii,2,k),  '(a2,i2)')    'C5',i
        write(ProtonName(1,ii,2,k),'(a1,i2,a1)') 'H',i,'O'
        write(ProtonName(2,ii,2,k),'(a1,i2,a1)') 'H',i,'P'
        write(ProtonName(3,ii,2,k),'(a1,i2,a1)') 'H',i,'Q'
      end do

      do i = 9, 1, -1
        ii = ii + 1
        write(CarbonName(ii,1,k),  '(a2,i1)')    'C6',i
        write(ProtonName(1,ii,1,k),'(a1,i1,a1)') 'H',i,'U'
        write(ProtonName(2,ii,1,k),'(a1,i1,a1)') 'H',i,'V'
        write(ProtonName(3,ii,1,k),'(a1,i1,a1)') 'H',i,'W'

        write(CarbonName(ii,2,k),  '(a2,i1)')    'C5',i
        write(ProtonName(1,ii,2,k),'(a1,i1,a1)') 'H',i,'O'
        write(ProtonName(2,ii,2,k),'(a1,i1,a1)') 'H',i,'P'
        write(ProtonName(3,ii,2,k),'(a1,i1,a1)') 'H',i,'Q'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if( (j ==  3).or.(j ==  7).or.(j == 11).or.(j == 15).or. &
          &   (j == 26).or.(j == 30).or.(j == 34).or.(j == 38) ) then
            NumProton(j,i,k) = 1
          else if( ( j >= 17 ).and.(j <= 24) ) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

   else if(LipidName(k) == 'TBPC') then

      do i = 1, 9
        write(CarbonName(i,1,k),  '(a2,i1)')    'C2',i
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i,'R'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i,'S'
        write(ProtonName(3,i,1,k),'(a1,i1,a1)') 'H',i,'T'

        write(CarbonName(i,2,k),  '(a2,i1)')    'C6',i
        write(ProtonName(1,i,2,k),'(a1,i1,a1)') 'H',i,'U'
        write(ProtonName(2,i,2,k),'(a1,i1,a1)') 'H',i,'V'
        write(ProtonName(3,i,2,k),'(a1,i1,a1)') 'H',i,'W'

        write(CarbonName(i,3,k),  '(a2,i1)')    'C3',i
        write(ProtonName(1,i,3,k),'(a1,i1,a1)') 'H',i,'X'
        write(ProtonName(2,i,3,k),'(a1,i1,a1)') 'H',i,'Y'
        write(ProtonName(3,i,3,k),'(a1,i1,a1)') 'H',i,'Z'
      end do  

      do i = 10, 20
        write(CarbonName(i,1,k),  '(a2,i2)')    'C2',i
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i,'R'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i,'S'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i,'T'

        write(CarbonName(i,2,k),  '(a2,i2)')    'C6',i
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i,'U'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i,'V'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i,'W'

        write(CarbonName(i,3,k),  '(a2,i2)')    'C3',i
        write(ProtonName(1,i,3,k),'(a1,i2,a1)') 'H',i,'X'
        write(ProtonName(2,i,3,k),'(a1,i2,a1)') 'H',i,'Y'
        write(ProtonName(3,i,3,k),'(a1,i2,a1)') 'H',i,'Z'
      end do

      ii = 20

      do i = 20, 10, -1
        ii = ii + 1
        write(CarbonName(ii,3,k),  '(a2,i2)')    'C5',i
        write(ProtonName(1,ii,3,k),'(a1,i2,a1)') 'H',i,'O'
        write(ProtonName(2,ii,3,k),'(a1,i2,a1)') 'H',i,'P'
        write(ProtonName(3,ii,3,k),'(a1,i2,a1)') 'H',i,'Q'
      end do

      do i = 9, 1, -1
        ii = ii + 1
        write(CarbonName(ii,3,k),  '(a2,i1)')    'C5',i
        write(ProtonName(1,ii,3,k),'(a1,i1,a1)') 'H',i,'O'
        write(ProtonName(2,ii,3,k),'(a1,i1,a1)') 'H',i,'P'
        write(ProtonName(3,ii,3,k),'(a1,i1,a1)') 'H',i,'Q'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if( (j ==  3).or.(j ==  7).or.(j == 11).or.(j == 15).or. &
          &   (j == 26).or.(j == 30).or.(j == 34).or.(j == 38) ) then
            NumProton(j,i,k) = 1
          else if( ( j >= 17 ).and.(j <= 24) ) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

   else if((LipidName(k) == 'POPC').or.(LipidName(k) == 'POPE')) then

      do i = 1, 8
        write(CarbonName(i,1,k),  '(a2,i1)')    'C2',i+1
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i+1,'R'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i+1,'S'
        write(ProtonName(3,i,1,k),'(a1,i1,a1)') 'H',i+1,'T'

        write(CarbonName(i,2,k),  '(a2,i1)')    'C3',i+1
        write(ProtonName(1,i,2,k),'(a1,i1,a1)') 'H',i+1,'X'
        write(ProtonName(2,i,2,k),'(a1,i1,a1)') 'H',i+1,'Y'
        write(ProtonName(3,i,2,k),'(a1,i1,a1)') 'H',i+1,'Z'
      end do  

      do i = 9, MaxCarbon
        write(CarbonName(i,1,k),  '(a2,i2)')    'C2',i+1
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i+1,'R'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i+1,'S'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i+1,'T'

        write(CarbonName(i,2,k),  '(a2,i2)')    'C3',i+1
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i+1,'X'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i+1,'Y'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i+1,'Z'
      end do

      write(ProtonName(1,8,1,k),'(a)') 'H91'
      write(ProtonName(1,9,1,k),'(a)') 'H101'

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if( (i == 1) .and. ( (j == 8).or.(j == 9) ) ) then
            NumProton(j,i,k) = 1
          else if( ( (i == 1) .and. (j == 17) ) .or. &
          &        ( (i == 2) .and. (j == 15) ) ) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

   else if(LipidName(k) == 'DOPC') then

      do i = 1, 8
        write(CarbonName(i,1,k),  '(a2,i1)')    'C2',i+1
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i+1,'R'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i+1,'S'
        write(ProtonName(3,i,1,k),'(a1,i1,a1)') 'H',i+1,'T'

        write(CarbonName(i,2,k),  '(a2,i1)')    'C3',i+1
        write(ProtonName(1,i,2,k),'(a1,i1,a1)') 'H',i+1,'X'
        write(ProtonName(2,i,2,k),'(a1,i1,a1)') 'H',i+1,'Y'
        write(ProtonName(3,i,2,k),'(a1,i1,a1)') 'H',i+1,'Z'
      end do  

      do i = 9, MaxCarbon
        write(CarbonName(i,1,k),  '(a2,i2)')    'C2',i+1
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i+1,'R'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i+1,'S'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i+1,'T'

        write(CarbonName(i,2,k),  '(a2,i2)')    'C3',i+1
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i+1,'X'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i+1,'Y'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i+1,'Z'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if( (j == 8).or.(j == 9) ) then
            NumProton(j,i,k) = 1
          else if(j == 17) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

   else if(LipidName(k) == 'SPM2') then

      do i = 1, 8
        write(CarbonName(i,1,k),  '(a1,i1,a1)')    'C',i+1,'F'
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i+1,'F'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i+1,'G'
        write(ProtonName(3,i,1,k),'(a1,i1,a1)') 'H',i+1,'H'
      end do
      do i = 9, MaxCarbon
        write(CarbonName(i,1,k),  '(a1,i2,a1)') 'C',i+1,'F'
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i+1,'F'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i+1,'G'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i+1,'H'
      end do

      do i = 1, 7
        write(CarbonName(i,2,k),  '(a1,i1,a1)') 'C',i+2,'S'
        write(ProtonName(1,i,2,k),'(a1,i1,a1)') 'H',i+2,'S'
        write(ProtonName(2,i,2,k),'(a1,i1,a1)') 'H',i+2,'T'
        write(ProtonName(3,i,2,k),'(a1,i1,a1)') 'H',i+2,'U'
      end do

      do i = 8, 16
        write(CarbonName(i,2,k),  '(a1,i2,a1)') 'C',i+2,'S'
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i+2,'S'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i+2,'T'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i+2,'U'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if( (i == 2) .and. ( (j == 1).or.(j == 2).or.(j == 3) ) ) then
            NumProton(j,i,k) = 1
          else if( ( (i == 1) .and. (j == 17) ) .or. &
          &        ( (i == 2) .and. (j == 16) ) ) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

   else if(LipidName(k) == 'SM18') then

      do i = 1, 8
        write(CarbonName(i,1,k),  '(a2,i1)')    'C2',i+1
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i+1,'R'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i+1,'S'
        write(ProtonName(3,i,1,k),'(a1,i1,a1)') 'H',i+1,'T'
      end do
      do i = 9, MaxCarbon
        write(CarbonName(i,1,k),  '(a2,i2)')    'C2',i+1
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i+1,'R'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i+1,'S'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i+1,'T'
      end do

      write(CarbonName(1,2,k),  '(a2)') 'C3'
      write(ProtonName(1,1,2,k),'(a2)') 'HX'
      do i = 2, 7
        write(CarbonName(i,2,k),  '(a2,i1)')    'C3',i+2
        write(ProtonName(1,i,2,k),'(a1,i1,a1)') 'H',i+2,'X'
        write(ProtonName(2,i,2,k),'(a1,i1,a1)') 'H',i+2,'Y'
        write(ProtonName(3,i,2,k),'(a1,i1,a1)') 'H',i+2,'Z'
      end do

      do i = 8, 16
        write(CarbonName(i,2,k),  '(a2,i2)')    'C3',i+2
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i+2,'X'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i+2,'Y'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i+2,'Z'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if( (i == 2) .and. ( (j == 1).or.(j == 2).or.(j == 3) ) ) then
            NumProton(j,i,k) = 1
          else if( ( (i == 1) .and. (j == 17) ) .or. &
          &        ( (i == 2) .and. (j == 16) ) ) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

   else if(LipidName(k) == 'SDS') then

      do i = 1, 9
        write(CarbonName(i,1,k),  '(a1,i1)')    'C',i
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i,'1'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i,'2'
      end do  

      do i = 10, 12
        write(CarbonName(i,1,k),  '(a1,i2)')    'C',i
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i,'1'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i,'2'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i,'3'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if(j == 12) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

   else if(LipidName(k) == 'CTAB') then

      do i = 1, 5
        write(CarbonName(i,1,k),  '(a1,i1)')    'C',i+4
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i+4,'1'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i+4,'2'
      end do  

      do i = 6, 15
        write(CarbonName(i,1,k),  '(a1,i2)')    'C',i+4
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i+4,'1'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i+4,'2'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i+4,'3'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if(j == 15) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

   else if(LipidName(k) == 'DTDA') then

      do i = 1, 5
        write(CarbonName(i,1,k),  '(a1,i1)')    'C',i+4
        write(ProtonName(1,i,1,k),'(a1,i1,a1)') 'H',i+4,'1'
        write(ProtonName(2,i,1,k),'(a1,i1,a1)') 'H',i+4,'2'
        write(ProtonName(3,i,1,k),'(a1,i1,a1)') 'H',i+4,'3'

        write(CarbonName(i,2,k),  '(a1,i2)')    'C',i+17
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i+17,'1'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i+17,'2'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i+17,'3'
      end do  

      do i = 6, 13
        write(CarbonName(i,1,k),  '(a1,i2)')    'C',i+4
        write(ProtonName(1,i,1,k),'(a1,i2,a1)') 'H',i+4,'1'
        write(ProtonName(2,i,1,k),'(a1,i2,a1)') 'H',i+4,'2'
        write(ProtonName(3,i,1,k),'(a1,i2,a1)') 'H',i+4,'3'

        write(CarbonName(i,2,k),  '(a1,i2)')    'C',i+17
        write(ProtonName(1,i,2,k),'(a1,i2,a1)') 'H',i+17,'1'
        write(ProtonName(2,i,2,k),'(a1,i2,a1)') 'H',i+17,'2'
        write(ProtonName(3,i,2,k),'(a1,i2,a1)') 'H',i+17,'3'
      end do

      do i = 1, Num_Chain(k)

        do j = 1, NumCarbon(i,k)

          if(j == 13) then
            NumProton(j,i,k) = 3
          else
            NumProton(j,i,k) = 2
          end if

        end do

      end do

    end if

    do i = 1 , NumMol(Nslipid(k))

      Numb1 = Numa(k) + (i-1)*NumAtm(Nslipid(k))

      do j = 1 , NumAtm(Nslipid(k))

        ll = Numb1 + j

        do ii = 1 , Num_Chain(k)

          do jj = 1 , NumCarbon(ii,k)

            if(AtomName(ll) == CarbonName(jj,ii,k)) then

              CarbonNumber(jj,ii,i,k) = ll

            end if

            do kk = 1, NumProton(jj,ii,k)

              if(AtomName(ll) == ProtonName(kk,jj,ii,k)) then

                ProtonNumber(kk,jj,ii,i,k) = ll

              end if

            end do

          end do

        end do

      end do

    end do

    end do

   end subroutine define_SCD


end subroutine SCD


!######################################################################
!######################################################################


subroutine DihedralAngle(Deg,i,j,k,l)

use Configuration, only : R
use UnitExParam, only : pi

implicit none

integer :: i, j, k, l
real(8), dimension(3) :: Rij, Rkj, Rlj
real(8), dimension(3) :: Pla, Plb, Ori
real(8) :: CsPhi, Phi, Ds, Dir, Deg
real(8) :: Ra2, Rb2, rRa2, rRb2, rRab

   Rij = R(:,i) - R(:,j)
   Rkj = R(:,k) - R(:,j)
   Rlj = R(:,l) - R(:,j)

   Pla = VecProd(Rij,Rkj)
   Plb = VecProd(Rlj,Rkj)

   Ra2 = dot_product(Pla,Pla)
   Rb2 = dot_product(Plb,Plb)

   rRa2 = 1.d0 / Ra2
   rRb2 = 1.d0 / Rb2
   rRab = sqrt(rRa2*rRb2)

   CsPhi = dot_product(Pla,Plb) * rRab

   if(CsPhi >  1.d0) CsPhi =  1.d0
   if(CsPhi < -1.d0) CsPhi = -1.d0

   Ori = VecProd(Pla,Plb)
   Ds  = dot_product(Ori,Rkj)
   Dir = 1.d0
   if(Ds < 0.d0) Dir = -1.d0

   Phi = acos(CsPhi)
   Deg = Dir * Phi * 180.d0 / pi

! #############

Contains

   function VecProd(x,y) Result(z)

   real(8), dimension(3) :: x, y, z

     z(1) = x(2) * y(3) - y(2) * x(3)
     z(2) = x(3) * y(1) - y(3) * x(1)
     z(3) = x(1) * y(2) - y(1) * x(2)

   end function VecProd

end subroutine DihedralAngle


!######################################################################
!######################################################################


subroutine AreaOccupL

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use UnitExParam, only : pi, Avogadro
use CellParam, only : H, InvH
use AtomParam, only : Mass
use TimeParam, only : Timeps

implicit none

real(8), parameter :: dd = 0.05d0

integer :: i, j, k, kk, l, ix, iy
integer :: nx, ny, TotalStepNumber, NumF
integer :: nm, nmh, Numa, Numb, Num
integer, dimension(2000000) :: Iocc
real(8), dimension(2,2000000) :: Sgrid
real(8) :: MassLipid
real(8), dimension(:,:), allocatable :: Rcom, Rg, Sg
real(8), dimension(3) :: a, b
real(8) :: TArea
real(8), dimension(2,2) :: SH
real(8) :: x0, y0, dx, dy, dxH, dyH, Rc2, R2
real(8), dimension(3) :: Ri, Si
real(8), dimension(2) :: Rij, Sij
real(8) :: SOccRatio, OccRatioLow, OccRatioUpp
integer :: Isum

   Rc2 = AreaL / pi

   open( 1,file='./Analy/OccArea.dat',form='formatted',status='unknown')

   nm  = NumMol(Kcomp)
   nmh = nm / 2

   Numa = 0
   Numb = NumMol(1)*NumAtm(1)

   if(Kcomp /= 1) then
     do i = 2, Kcomp
       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)
     end do
   end if

   MassLipid = 0.d0
   do i = Numa+1, Numa+NumAtm(Kcomp)
     MassLipid = MassLipid + Mass(i)
   end do

   print *, 'MassLipid=',MassLipid*Avogadro*1.d+3

   allocate( Rcom(3,nm) )
   allocate( Rg(2,nm) )
   allocate( Sg(2,nm) )

! ## clear

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)
   end do

   NumF = 0
   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       NumF = NumF + 1
!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

       call CellTransform

       a  = H(:,1)
       b  = H(:,2)

       SH(1,1) = H(1,1)
       SH(1,2) = H(1,2)
       SH(2,1) = H(2,1)
       SH(2,2) = H(2,2)

       call CArea(a,b,TArea)

       Num = Numa

       do k = 1, NumMol(Kcomp)
         Rcom(:,k) = 0.d0
         do l = 1, NumAtm(Kcomp)
           Num = Num + 1
           Rcom(:,k) = Rcom(:,k) + Mass(Num)*R(:,Num)
         end do
         Rcom(:,k) = Rcom(:,k) / MassLipid
       end do

       do k = 1 , nm
         Ri = Rcom(:,k)
         Si = matmul(InvH,Ri)
         Si = Si - nint(Si)
         Sg(:,k) = Si(:)
         Rg(1,k) = dot_product( H(1,:),Si(:) )
         Rg(2,k) = dot_product( H(2,:),Si(:) )
       end do

       nx = int(a(1)/dd) + 1
       ny = int(b(2)/dd) + 1

       dx  = 1.d0 / dble(nx)
       dy  = 1.d0 / dble(ny)

       dxH  = dx  * 0.5d0
       dyH  = dy  * 0.5d0

       Ngrid = nx * ny

       x0 = -0.5d0 + dxH
       y0 = -0.5d0 + dyH

       do ix = 1, nx
         k = (ix-1) * ny

         do iy = 1, ny
           kk = k + iy

           Sgrid(1,kk) = x0 + ix * dx
           Sgrid(2,kk) = y0 + iy * dy

         end do

       end do

!  ****************** upper layer ***************

       do k = 1 , Ngrid
         Iocc(k) = 0
       end do

       do k = 1 , nmh

         do kk = 1 , Ngrid

           if(Iocc(kk)==1) cycle

           Sij = Sg(:,k) - Sgrid(:,kk)
           Sij = Sij - nint(Sij)
           Rij = matmul( SH , Sij )

           R2  = dot_product( Rij, Rij )

           if(R2 <= Rc2) Iocc(kk)=1

         end do

       end do

       Isum = 0
       do k = 1 , Ngrid
         Isum = Isum + Iocc(k)
       end do

       OccRatioUpp = dble(Isum) / dble( Ngrid )

!  ****************** lower layer ***************

       do k = 1 , Ngrid
         Iocc(k) = 0
       end do

       do k = nmh+1 , nm

         do kk = 1 , Ngrid

           if(Iocc(kk)==1) cycle

           Sij = Sg(:,k) - Sgrid(:,kk)
           Sij = Sij - nint(Sij)
           Rij = matmul( SH , Sij )

           R2  = dot_product( Rij, Rij )

           if(R2 <= Rc2) Iocc(kk)=1

         end do

       end do

       Isum = 0
       do k = 1 , Ngrid
         Isum = Isum + Iocc(k)
       end do

       OccRatioLow = dble(Isum) / dble( Ngrid )

       SOccRatio = SOccRatio + OccRatioUpp + OccRatioLow

       write(1,'(f9.4,2f9.5)') Timeps, OccRatioUpp, OccRatioLow

     end do

   end do

   SOccRatio = SOccRatio / dble(TotalStepNumber * 2)

   write(1,'(a,f9.5)') 'average = ',SOccRatio

   close(1)

end subroutine AreaOccupL


!######################################################################
!######################################################################


subroutine VoronoiLipid

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use UnitExParam, only : Avogadro
use CellParam, only : H, InvH
use AtomParam, only : AtomName, Mass
use TimeParam, only : Timeps

implicit none

character(len=4) :: AName
integer :: i, j, TotalStepNumber, nm, nmh, nmh9
integer :: k, l, ix, iy, nx, ny, lx, ly, nn
integer :: ii, jj, kk, ll, j0, k0, l0, NumF, iq
integer :: Numa, Numb, Num, k1, k2
real(8), dimension(:,:), allocatable :: Rcom, Rg
real(8), dimension(:,:), allocatable :: Rupp, Rlow
real(8), dimension(:,:), allocatable :: RList
real(8), dimension(:), allocatable :: RR2
integer, dimension(:,:), allocatable :: NumSide
integer, dimension(:), allocatable :: Lpair
real(8), dimension(:,:), allocatable :: MArea
real(8), dimension(15) :: Side, Sr
real(8), dimension(15,2) :: Pr
integer, dimension(15,2) :: Ip
real(8), parameter :: RVcutoff = 20.d0
real(8), parameter :: RVcut2 = RVcutoff * RVcutoff
real(8) :: MassLipid
real(8), dimension(3) :: Ri, Si
real(8), dimension(2) :: Zr, Rgi, Rij, Rji, Rki, Rli, Rc, Sc
real(8) :: Rji2, Rki2, R1, R2, ab, ab2, aaa, area
real(8), dimension(3) :: a, b
real(8) :: TArea, SumAu, SumAd, s
integer :: ITArea, ISumAu, ISumAd
integer, dimension(:), allocatable :: AFlag
real(8) :: NumC2, NumC3


if(cWhat=='VoronoiL') then

   open( 1,file='./Analy/Vor_Area.dat',form='unformatted',status='unknown')
   open( 2,file='./Analy/Vor_Side.dat',form='unformatted',status='unknown')
   open( 3,file='./Analy/Dis_Side.dat',status='unknown')

   open(51,file='./Analy/PolygonLupp.dat',status='unknown')
   open(52,file='./Analy/PolygonLlow.dat',status='unknown')
   open(53,file='./Analy/Vor_Rg.dat',status='unknown')
   open(54,file='./Analy/Vor_Rupp.dat',status='unknown')
   open(55,file='./Analy/Vor_Rlow.dat',status='unknown')

else if(cWhat=='VoronoiC') then

   open( 1,file='./Analy/Vor_AreaC.dat',form='unformatted',status='unknown')
   open( 2,file='./Analy/Vor_SideC.dat',form='unformatted',status='unknown')
   open( 3,file='./Analy/Dis_SideC.dat',status='unknown')

   open(51,file='./Analy/PolygonLuppC.dat',status='unknown')
   open(52,file='./Analy/PolygonLlowC.dat',status='unknown')
   open(53,file='./Analy/Vor_RgC.dat',status='unknown')
   open(54,file='./Analy/Vor_RuppC.dat',status='unknown')
   open(55,file='./Analy/Vor_RlowC.dat',status='unknown')

end if

   nm  = NumMol(Kcomp)
   nmh = nm / 2

   nmh9 = nmh * 9

   Numa = 0
   Numb = NumMol(1)*NumAtm(1)

   if(Kcomp /= 1) then
     do i = 2, Kcomp
       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)
     end do
   end if

   MassLipid = 0.d0
   do i = Numa+1, Numa+NumAtm(Kcomp)
     MassLipid = MassLipid + Mass(i)
   end do

   print *, 'MassLipid=',MassLipid*Avogadro*1.d+3

if(cWhat=='VoronoiC') then
   nm   = nm * 2
   nmh  = nmh * 2
   nmh9 = nmh9 * 2

   allocate( AFlag(NumAtm(Kcomp)) )
   AFlag = 0

   do i = 1 , NumAtm(Kcomp)
     j = i + Numa
     AName = AtomName(j)
     if((AName(1:2) == 'C2').and.(AName(3:3) /= ' ')) then
       AFlag(i) = 1
     else if((AName(1:2) == 'C3').and.(AName(3:3) /= ' ')) then
       AFlag(i) = 2
     end if
   end do

   ii = 0
   jj = 0
   do i = 1 , NumAtm(Kcomp)
     if(AFlag(i) == 1) then
       ii = ii + 1
     else if(AFlag(i) == 2) then
       jj = jj + 1
     end if
   end do

   print * , 'Number of C in the sn-1 chain = ', jj
   print * , 'Number of C in the sn-2 chain = ', jj

   NumC2 = dble(ii)
   NumC3 = dble(jj)

end if

   allocate( Rcom(3,nm) )
   allocate( Rg(2,nm) )
   allocate( Rupp(2,nmh9) )
   allocate( Rlow(2,nmh9) )
   allocate( RList(2,nm*10) )
   allocate( RR2(nm*10) )
   allocate( Lpair(nm*10) )

! ## clear
   Side = 0.d0

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)
   end do

   allocate( MArea(nm,TotalStepNumber) )
   allocate( NumSide(nm,TotalStepNumber) )

   NumF = 0
   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       NumF = NumF + 1
!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

       call CellTransform

       a  = H(:,1)
       b  = H(:,2)

       call CArea(a,b,TArea)

       Num = Numa

if(cWhat=='VoronoiL') then

       do k = 1, NumMol(Kcomp)
         Rcom(:,k) = 0.d0
         do l = 1, NumAtm(Kcomp)
           Num = Num + 1
           Rcom(:,k) = Rcom(:,k) + Mass(Num)*R(:,Num)
         end do
         Rcom(:,k) = Rcom(:,k) / MassLipid
       end do

else if(cWhat=='VoronoiC') then

       do k = 1, NumMol(Kcomp)
         k1 = 2*k - 1
         k2 = 2*k
         Rcom(:,k1) = 0.d0
         Rcom(:,k2) = 0.d0
         do l = 1, NumAtm(Kcomp)
           Num = Num + 1
           if(AFlag(l) == 1) then
             Rcom(:,k1) = Rcom(:,k1) + R(:,Num)
           else if(AFlag(l) == 2) then
             Rcom(:,k2) = Rcom(:,k2) + R(:,Num)
           end if
         end do
         Rcom(:,k1) = Rcom(:,k1) / NumC2
         Rcom(:,k2) = Rcom(:,k2) / NumC3
       end do

end if

       do k = 1 , nm
         Ri = Rcom(:,k)
         Si = matmul(InvH,Ri)
         Si = Si - nint(Si)
         Rg(1,k) = dot_product( H(1,:),Si(:) )
         Rg(2,k) = dot_product( H(2,:),Si(:) )
       end do

       do ix = 1 , 3

         nx = ix - 2
         lx = (ix-1) * 3 * nmh

         do iy = 1 , 3
           ny = iy - 2
           ly = (iy-1) * nmh + lx

           Zr(1) = nx * H(1,1) + ny * H(1,2)
           Zr(2) = nx * H(2,1) + ny * H(2,2)

           do k = 1 , nmh
             nn = ly + k
             kk = k + nmh
             Rupp(:,nn) = Rg(:, k) + Zr(:)
             Rlow(:,nn) = Rg(:,kk) + Zr(:)
           end do

         end do

       end do

!-----------------------------------------------------------
!  ****************** upper layer ***************
!-----------------------------------------------------------

       do k = 1 , nmh

         nn   = 4 * nmh + k
         numb = 0
         Rgi = Rg(:,k)

         do l = 1 , nmh9

           if(l==nn) cycle

           Rij = Rupp(:,l) - Rgi
           R2  = dot_product(Rij,Rij)

           if(R2 < RVcut2) then

             numb = numb + 1

             Lpair(numb) = l      ! the number of tagged particle
             RList(:,numb)  = Rij    ! vector
             RR2(numb)   = R2     ! distance^2

           end if

         end do

         call Voronoi2d(iq,area,1)

         NumSide(k,NumF) = iq
         MArea(k,NumF) = area

         Side(iq) = Side(iq)+1

       end do

!-----------------------------------------------------------
!  ****************** lower layer ***************
!-----------------------------------------------------------

       do k = 1 , nmh

         nn   = 4 * nmh + k
         numb = 0
         ii   = nmh + k
         Rgi  = Rg(:,ii)

         do l = 1 , nmh9

           if(l==nn) cycle

           Rij = Rlow(:,l) - Rgi
           R2  = dot_product(Rij,Rij)

           if(R2 < RVcut2) then

             numb = numb + 1

             Lpair(numb) = l      ! the number of tagged particle
             RList(:,numb)  = Rij    ! vector
             RR2(numb)   = R2     ! distance^2

           end if

         end do

         call Voronoi2d(iq,area,2)

         NumSide(ii,NumF) = iq
         MArea(ii,NumF) = area

         Side(iq) = Side(iq)+1

       end do

! ## check the consistency

       SumAu=0.
       SumAd=0.

       do k = 1 , nmh

         l = k + nmh
         SumAu = SumAu + MArea(k,NumF)
         SumAd = SumAd + MArea(l,NumF)

       end do

       ITArea = int(TArea * 100.)
       ISumAu = int(SumAu * 100.)
       ISumAd = int(SumAd * 100.)

       if( (ISumAu /= ITArea) .or. (ISumAd /= ITArea) ) then
         write(*,*) 'sum_a not conserved'
         write(*,*) SumAu,SumAd,TArea
         write(*,*) 'step=',NumF
         write(*,*) 'time=',Timeps
         write(*,*) 'rc length may be insufficient?'
         write(*,*) 'area='
         write(*,'(10(8f9.2/))') (sngl(MArea(k,NumF)),k=1,nm)
         write(*,*) 'number of side='
         write(*,'(5(16i4/))') (NumSide(k,NumF),k=1,nm)
         write(*,*) 'Rg='
         write(*,'(20(8f9.2/))') (sngl(Rg(:,k)),k=1,nm)
         write(*,*) 'Rupp='
         write(*,'(180(8f9.2/))') (sngl(Rupp(:,k)),k=1,nmh9)
         write(*,*) 'Rlow='
         write(*,'(180(8f9.2/))') (sngl(Rlow(:,k)),k=1,nmh9)
         call Finalize
       end if

       write(1) Timeps , ( MArea(k,NumF) , k = 1 , nm )
       write(2) Timeps , ( NumSide(k,NumF) , k = 1 , nm )

       if(NumF == TotalStepNumber) then

         do k = 1 , nm
           write(53,'(2f8.3)') Rg(:,k)
         end do

         do k = 1 , nmh9
           write(54,'(2f8.3)') Rupp(:,k)
           write(55,'(2f8.3)') Rlow(:,k)
         end do

       end if

     end do

   end do

   do i = 1 , 15
     s = Side(i) / dble(TotalStepNumber*nm)
     write(3,'(i10,f15.6)') i, s*1.d+02
   end do

   close(1)
   close(2)
   close(3)
   close(51)
   close(52)
   close(53)
   close(54)
   close(55)

contains

   subroutine Voronoi2d(iq,area,iflag)

   implicit none

     integer :: iq, ifile, iflag
     real(8) :: area, bo
     logical :: Vflag

     ifile = 50 + iflag

     iq   = 0
     area = 0.d0

     do jj = 1 , numb - 1

       j0   = Lpair(jj)
       Rji  = RList(:,jj)
       Rji2 = RR2(jj)

       do kk = jj + 1 , numb

         k0   = Lpair(kk)
         Rki  = RList(:,kk)
         Rki2 = RR2(kk)

         bo = 2.d0 * (Rji(1) * Rki(2) - Rki(1) * Rji(2))
         Sc(1) = (Rji2*Rki(2) - Rki2*Rji(2)) / bo
         Sc(2) = (Rki2*Rji(1) - Rji2*Rki(1)) / bo

         Rc = Rgi + Sc

         R1 = dot_product(Sc,Sc)

         Vflag = .True.

         do ll = 1 , numb

           l0 = Lpair(ll)
           if( (l0 == j0) .or. (l0 == k0) ) cycle

           if(iflag == 1) Rli = Rupp(:,l0) - Rc
           if(iflag == 2) Rli = Rlow(:,l0) - Rc

           R2 = dot_product(Rli,Rli)

           if(R2 < R1) then
             Vflag = .False.
             exit
           end if

         end do

         if(Vflag) then

           iq = iq + 1

           Ip(iq,1) = j0
           Ip(iq,2) = k0

           Pr(iq,:) = Sc

           Sr(iq)   = R1

         end if

       end do

     end do

     do l = 1 , iq - 1

       j0 = Ip(l,1)
       k0 = Ip(l,2)

       Sc = Pr(l,:)
       R1 = Sr(l)

       do ll = l + 1 , iq

         if( (ip(ll,1)==j0) .or. (ip(ll,1)==k0) ) then
           ab   = dot_product(Pr(ll,:),Sc(:))
           ab2  = ab * ab
           aaa  = 0.5d0 * sqrt(R1 * Sr(ll) - ab2)
           area = area + aaa

           if(NumF==TotalStepNumber) then
             write(ifile,'(2f8.2)') Sc(1)+Rgi(1),Sc(2)+Rgi(2)
             write(ifile,'(2f8.2)') Pr(ll,1)+Rgi(1),Pr(ll,2)+Rgi(2)
             write(ifile,*) '= ='
           end if

         else if( (ip(ll,2) == j0) .or. (ip(ll,2) == k0) ) then

           ab   = dot_product(Pr(ll,:),Sc(:))
           ab2  = ab * ab
           aaa  = 0.5 * sqrt( R1 * Sr(ll) - ab2 )
           area = area + aaa

           if(NumF==TotalStepNumber) then
             write(ifile,'(2f8.2)') Sc(1)+Rgi(1),Sc(2)+Rgi(2)
             write(ifile,'(2f8.2)') Pr(ll,1)+Rgi(1),Pr(ll,2)+Rgi(2)
             write(ifile,*) '= ='
           end if

         end if

       end do

     end do

   end subroutine Voronoi2d

end subroutine VoronoiLipid


!######################################################################
!######################################################################


subroutine GrGG_Lipid

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use CommonBlocks, only : QMaster
use ParamAnalyze
use IOparam, only : trajectory_file
use UnitExParam, only : pi
use CellParam, only : H, InvH
use CutoffParam, only : Rcutoff2
use AtomParam, only : Mass

implicit none

integer :: i, j, k, l, ir, Num, NumF
integer :: nm, nmh, Numa, Numb, TotalStepNumber
real(8) :: MassLipid, R2, R1
real(8), dimension(3) :: Rij, Sij, a, b
real(8) :: AreaT, SumArea
integer, dimension(:), allocatable :: igg
real(8), dimension(:), allocatable :: ggg, xgg
real(8), dimension(:,:), allocatable :: Rcom, Scom
real(8) :: zero, dr, qkpp, zz

   IRcut = int(sqrt(Rcutoff2)*4.) + 1

   allocate( igg(IRcut) )
   allocate( ggg(IRcut) )
   allocate( xgg(IRcut) )

   igg = 0

   nm  = NumMol(Kcomp)
   nmh = nm / 2

   Numa = 0
   Numb = NumMol(1)*NumAtm(1)

   if(Kcomp /= 1) then
     do i = 2, Kcomp
       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)
     end do
   end if

   MassLipid = 0.d0
   do i = Numa+1, Numa+NumAtm(Kcomp)
     MassLipid = MassLipid + Mass(i)
   end do

   allocate( Rcom(3,nm) )
   allocate( Scom(3,nm) )

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)
   end do

   SumArea = 0.d0

   NumF = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

       NumF = NumF + 1

       call CellTransform

       a  = H(:,1)
       b  = H(:,2)

       call CArea(a,b,AreaT)

       SumArea = SumArea + AreaT

       Num = Numa

       do k = 1, NumMol(Kcomp)

         Rcom(:,k) = 0.d0

         do l = 1, NumAtm(Kcomp)

           Num = Num + 1
           Rcom(:,k) = Rcom(:,k) + Mass(Num)*R(:,Num)

         end do

         Rcom(:,k) = Rcom(:,k) / MassLipid

       end do

       do k = 1 , nm

         Scom(:,k) = matmul( InvH, Rcom(:,k) )

       end do

       do k = 1 , nmh-1

         do l = k+1, nmh

           Sij = Scom(:,k) - Scom(:,l)
           Sij = Sij - nint( Sij )
           Rij = matmul( H, Sij )
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2)

           if( R2 < Rcutoff2 ) then
             R1 = sqrt(R2)
             ir = nint(R1*4.d0)
             if(ir==0) then
               print *, 'ir=',ir
               print *, 'numf=',NumF
               print *, 'pair=',k,l
               print *, trajectory_file
               cycle
             end if
             igg(ir) = igg(ir) + 1
           end if

         end do

       end do

       do k = nmh+1 , nm-1

         do l = k+1, nm

           Sij = Scom(:,k) - Scom(:,l)
           Sij = Sij - nint( Sij )
           Rij = matmul( H, Sij )
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2)

           if( R2 < Rcutoff2 ) then
             R1 = sqrt(R2)
             ir = nint(R1*4.d0)
             if(ir==0) then
               print *, 'ir=',ir
               print *, 'numf=',NumF
               print *, 'pair=',k,l
               print *, trajectory_file
               cycle
             end if
             igg(ir) = igg(ir) + 1
           end if

         end do

       end do

     end do

     write(*,*) 'file, ',i,' was completed '

   end do

   AreaT = SumArea / dble(TotalStepNumber)

   do i = 1, IRcut-1
     dr   = i * 0.25d0
     qkpp = AreaT / ( nmh * 2.d0 * pi * dr * 0.25d0 )
     xgg(i) = dble(igg(i)) / (nmh*TotalStepNumber*2.d0) * 2.d0
     ggg(i) = xgg(i) * qkpp
   end do

   i  = IRcut
   dr = i * 0.25d0
   qkpp = AreaT / (nmh * 2.d0 * pi * dr * 0.25d0 )
   xgg(i) = dble(igg(i)) / (nmh*TotalStepNumber*2.d0) * 2.d0 * 2.d0
   ggg(i) = xgg(i) * qkpp

   zz=0.
   do i = 1 , IRcut
     zz     = zz + xgg(i)
     xgg(i) = zz
   end do

   open(1,file='./Analy/GR_LipidGG.dat',status='unknown')

   zero=0.
   write(1,'(1x,f7.2,1x,f12.7,1x,f12.7)') zero,zero,zero
   do i = 1 , IRcut
     write(1,'(1x,f7.2,1x,f12.7,1x,f12.7)') real(i)/4. , ggg(i) , xgg(i)
   end do

   close(1)

end subroutine GrGG_Lipid


!######################################################################
!######################################################################


subroutine TGChain

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : MolName, AtomName

implicit none

integer :: i , j, k, l, NumDihChain
integer :: ii, jj, kk, ll
integer :: Numa, Numb
integer :: NumF
integer :: TotalStepNumber
integer, dimension(:,:,:), allocatable :: DihedBeta,DihedGamm
real(8) :: DegG, DegB
character(len=4) :: LipidName
integer, dimension(:,:), allocatable :: curtb, curtg
integer, dimension(:,:), allocatable :: pretb, pretg
integer, dimension(:,:), allocatable :: tgb, tgg
real, dimension(:,:,:), allocatable :: rdihb, rdihg
real(8), dimension(3) :: RiB, RiG

   ii = 0

   do i = 1, NumSpec

     if((MolName(i)(1:4) == 'DPPC') .or. &
     &  (MolName(i)(1:4) == 'DMPC') .or. &
     &  (MolName(i)(1:4) == 'PhPC') ) then
       LipidName = MolName(i)(1:4)
       ii = ii + 1
       Nlipid = i
     end if

   end do

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   if( ( LipidName == 'DPPC' ) .or. ( LipidName == 'PhPC' ) ) then

     allocate( DihedBeta(4,13,NumMol(Nlipid)) )
     allocate( DihedGamm(4,13,NumMol(Nlipid)) )

   else if( LipidName == 'DMPC' ) then

     allocate( DihedBeta(4,11,NumMol(Nlipid)) )
     allocate( DihedGamm(4,11,NumMol(Nlipid)) )

   end if

   call Dihed_LipidChain

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)/Interval(i)
     if(i>1) then
       if(Interval(i)/=Interval(i-1)) then
         write(*,*) 'ERROR : Interval'
         write(*,*) 'Interval = ', Interval(i),Interval(i-1)
         call Finalize
       end if
     end if
   end do

   if(Nsnap /= TotalStepNumber) then
     write(*,*) 'ERROR : Nsnap'
     write(*,*) 'total step number = ',TotalStepNumber
     write(*,*) 'Nsnap             = ',Nsnap
   end if

   allocate( curtb(NumDihChain,NumMol(Nlipid)) )
   allocate( curtg(NumDihChain,NumMol(Nlipid)) )
   allocate( pretb(NumDihChain,NumMol(Nlipid)) )
   allocate( pretg(NumDihChain,NumMol(Nlipid)) )
   allocate( tgb(NumDihChain,NumMol(Nlipid)) )
   allocate( tgg(NumDihChain,NumMol(Nlipid)) )
   allocate( rdihb(3,NumDihChain,NumMol(Nlipid)) )
   allocate( rdihg(3,NumDihChain,NumMol(Nlipid)) )

   open(7,file='./Analy/TGdata.dat',form='unformatted',status='unknown')

   NumF = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

       if(mod(j,Interval(i))/=0) cycle

       NumF = NumF + 1

       do k = 1 , NumMol(Nlipid)

         do l = 1 , NumDihChain

           ii = DihedBeta(1,l,k)
           jj = DihedBeta(2,l,k)
           kk = DihedBeta(3,l,k)
           ll = DihedBeta(4,l,k)

           call DihedralAngle( DegB, ii, jj, kk, ll )
           call DihPosition( RiB, jj, kk )

           ii = DihedGamm(1,l,k)
           jj = DihedGamm(2,l,k)
           kk = DihedGamm(3,l,k)
           ll = DihedGamm(4,l,k)

           call DihedralAngle( DegG, ii, jj, kk, ll )
           call DihPosition( RiG, jj, kk )

           if(abs(DegB) > 120.d0) then
             curtb(l,k) =  1
           else
             curtb(l,k) = -1
           end if

           if(abs(DegG) > 120.d0) then
             curtg(l,k) =  1
           else
             curtg(l,k) = -1
           end if

           rdihb(:,l,k) = sngl(RiB(:))
           rdihg(:,l,k) = sngl(RiG(:))

         end do

       end do

       if(NumF /= 1) then

         do k = 1 , NumMol(Nlipid)

           do l = 1 , NumDihChain

             ii = curtb(l,k) * pretb(l,k)
             if(ii==1) then
               tgb(l,k) = 0
             else
               tgb(l,k) = 1
             end if

             ii = curtg(l,k) * pretg(l,k)
             if(ii==1) then
               tgg(l,k) = 0
             else
               tgg(l,k) = 1
             end if

           end do

         end do

         write(7) tgb,tgg,rdihb,rdihg

       end if

       pretb = curtb
       pretg = curtg

     end do

   end do

!-----------------------------------------------------------

Contains

   subroutine DihPosition( Ri, i, j )

   implicit none

   integer :: i, j
   real(8), dimension(3) :: Ri

      Ri(:) = 0.5d0 * ( R(:,i) + R(:,j) )

   end subroutine DihPosition


   subroutine Dihed_LipidChain

   implicit none

   character(len=4), dimension(:), allocatable :: AtomBeta,AtomGamm
   integer , dimension(:), allocatable :: NaBeta, NaGamm
   integer :: NumCarb, Numb1
   logical, dimension(:), allocatable :: FlagBeta,FlagGamm

      if((LipidName == 'DPPC').or.(LipidName == 'PhPC')) then

        NumCarb = 16
        allocate( AtomBeta(16) )
        allocate( AtomGamm(16) )
        allocate( NaBeta(16) )
        allocate( NaGamm(16) )
        allocate( FlagBeta(16) )
        allocate( FlagGamm(16) )

        AtomBeta( 1) = 'C21'
        AtomBeta( 2) = 'C22'
        AtomBeta( 3) = 'C23'
        AtomBeta( 4) = 'C24'
        AtomBeta( 5) = 'C25'
        AtomBeta( 6) = 'C26'
        AtomBeta( 7) = 'C27'
        AtomBeta( 8) = 'C28'
        AtomBeta( 9) = 'C29'
        AtomBeta(10) = 'C210'
        AtomBeta(11) = 'C211'
        AtomBeta(12) = 'C212'
        AtomBeta(13) = 'C213'
        AtomBeta(14) = 'C214'
        AtomBeta(15) = 'C215'
        AtomBeta(16) = 'C216'

        AtomGamm( 1) = 'C31'
        AtomGamm( 2) = 'C32'
        AtomGamm( 3) = 'C33'
        AtomGamm( 4) = 'C34'
        AtomGamm( 5) = 'C35'
        AtomGamm( 6) = 'C36'
        AtomGamm( 7) = 'C37'
        AtomGamm( 8) = 'C38'
        AtomGamm( 9) = 'C39'
        AtomGamm(10) = 'C310'
        AtomGamm(11) = 'C311'
        AtomGamm(12) = 'C312'
        AtomGamm(13) = 'C313'
        AtomGamm(14) = 'C314'
        AtomGamm(15) = 'C315'
        AtomGamm(16) = 'C316'

      else if(LipidName == 'DMPC') then

        NumCarb = 14
        allocate( AtomBeta(14) )
        allocate( AtomGamm(14) )
        allocate( NaBeta(14) )
        allocate( NaGamm(14) )
        allocate( FlagBeta(14) )
        allocate( FlagGamm(14) )

        AtomBeta( 1) = 'C21'
        AtomBeta( 2) = 'C22'
        AtomBeta( 3) = 'C23'
        AtomBeta( 4) = 'C24'
        AtomBeta( 5) = 'C25'
        AtomBeta( 6) = 'C26'
        AtomBeta( 7) = 'C27'
        AtomBeta( 8) = 'C28'
        AtomBeta( 9) = 'C29'
        AtomBeta(10) = 'C210'
        AtomBeta(11) = 'C211'
        AtomBeta(12) = 'C212'
        AtomBeta(13) = 'C213'
        AtomBeta(14) = 'C214'

        AtomGamm( 1) = 'C31'
        AtomGamm( 2) = 'C32'
        AtomGamm( 3) = 'C33'
        AtomGamm( 4) = 'C34'
        AtomGamm( 5) = 'C35'
        AtomGamm( 6) = 'C36'
        AtomGamm( 7) = 'C37'
        AtomGamm( 8) = 'C38'
        AtomGamm( 9) = 'C39'
        AtomGamm(10) = 'C310'
        AtomGamm(11) = 'C311'
        AtomGamm(12) = 'C312'
        AtomGamm(13) = 'C313'
        AtomGamm(14) = 'C314'

      end if

      NumDihChain = NumCarb - 3

      do i = 1 , NumMol(Nlipid)

        Numb1 = Numa + (i-1)*NumAtm(Nlipid)

        FlagBeta = .false.
        FlagGamm = .false.

        do j = 1 , NumAtm(Nlipid)

          k = Numb1 + j

          do ii = 1 , NumCarb

            if(AtomName(k) == AtomBeta(ii)) then

              NaBeta(ii) = k
              FlagBeta(ii) = .true.

            else if(AtomName(k) == AtomGamm(ii)) then

              NaGamm(ii) = k
              FlagGamm(ii) = .true.

            end if

          end do

        end do

        do ii = 1 , NumCarb
          if((.not.FlagBeta(ii)).or.(.not.FlagGamm(ii))) then
            write(*,*) 'ERROR : Carbon number has not been assigned sucessfully.'
            write(*,*) ii
            call Finalize
          end if
        end do

        do j = 1, NumDihChain

          DihedBeta(1,j,i) = NaBeta(j  )
          DihedBeta(2,j,i) = NaBeta(j+1)
          DihedBeta(3,j,i) = NaBeta(j+2)
          DihedBeta(4,j,i) = NaBeta(j+3)

          DihedGamm(1,j,i) = NaGamm(j  )
          DihedGamm(2,j,i) = NaGamm(j+1)
          DihedGamm(3,j,i) = NaGamm(j+2)
          DihedGamm(4,j,i) = NaGamm(j+3)

        end do

      end do

   end subroutine Dihed_LipidChain


end subroutine TGChain


!######################################################################
!######################################################################


subroutine NeighborChainCorr

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H, InvH
use CutoffParam, only : Rcutoff2
use AtomParam, only : MolName, AtomName, Mass

implicit none

integer :: i, j, k, l, ii, jj
integer :: Numa, Numb, NumF
character(len=4) :: LipidName
character(len=4), dimension(:,:), allocatable :: AtomChain
integer, dimension(:,:,:), allocatable :: IAtmChain
integer :: NumAtmChain, NumLipid
real(8), dimension(:,:,:), allocatable :: chainvector, Rgl
real(8) :: MassMol, Rc
! ## <SP>
real(8), dimension(:,:,:), allocatable :: Sgl
integer :: NumSample1
real(8), dimension(3) :: Rij, Sij, v1, v2
real(8) :: R2, d1, d2, cs, corr
integer :: kk, ll
! ## <SP>


   ii = 0

   do i = 1, NumSpec

     if((MolName(i)(1:4) == 'DPPC') .or. &
     &  (MolName(i)(1:4) == 'DMPC') .or. &
     &  (MolName(i)(1:4) == 'PhPC') ) then
       LipidName = MolName(i)(1:4)
       ii = ii + 1
       Nlipid = i
     end if

   end do

   print *, 'Lipid name =',LipidName

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   NumLipid = NumMol(Nlipid)

   if( LipidName == 'DPPC' ) then
     NumAtmChain = 16
! ## <SP>
     Rc = 8.7d0
     Rcutoff2 = Rc**2
! ## <SP>
   else if( LipidName == 'PhPC' ) then
     NumAtmChain = 20
! ## <SP>
     Rc = 9.6d0
     Rcutoff2 = Rc**2
! ## <SP>
   end if

   allocate( AtomChain(2,NumAtmChain) )
   allocate( IAtmChain(2,NumAtmChain,NumLipid) )
   allocate( chainvector(3,2,NumLipid) )
   allocate( Rgl(3,2,NumLipid) )
! ## <SP>
   allocate( Sgl(3,2,NumLipid) )
   NumSample1 = 0
   corr = 0.d0
! ## <SP>

   call NameOfChain

   do i = 1, NumLipid

     do j = 1, NumAtm(Nlipid)

       ii = Numa + (i-1)*NumAtm(Nlipid) + j

       do k = 1, NumAtmChain

         if(AtomChain(1,k) == AtomName(ii)) then
           IAtmChain(1,k,i) = ii
         else if(AtomChain(2,k) == AtomName(ii)) then
           IAtmChain(2,k,i) = ii
         end if

       end do

     end do

   end do

   MassMol = 0.d0

   do j = 1, NumAtmChain

inn: do i = Numa+1, Numa+NumAtm(Nlipid)

       if(AtomName(i) == AtomChain(1,j)) then
         MassMol = MassMol + Mass(i)
         exit inn
       end if

       if(i==Numa+NumAtm(Nlipid)) then
         write(*,*) 'ERROR : no atom'
         call Finalize
       end if

     end do inn

   end do

!   print *, 'MassMol = ', MassMol*Avogadro*1.d3

   NumF = 0

   open(7,file='Rg.dat',form='unformatted')
   open(8,file='Chainvec.dat',form='unformatted')

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

!       call CellTransform

       call InversMatrix(H,InvH)

       if(mod(j,Interval(i))/=0) cycle

       NumF = NumF + 1

       do k = 1, NumLipid

         Rgl(:,:,k) = 0.d0

         do l = 1, NumAtmChain

           ii = IAtmChain(1,l,k)
           jj = IAtmChain(2,l,k)

           Rgl(:,1,k) = Rgl(:,1,k) + Mass(ii)*R(:,ii)
           Rgl(:,2,k) = Rgl(:,2,k) + Mass(jj)*R(:,jj)

         end do

         Rgl(:,1,k) = Rgl(:,1,k) / MassMol
         Rgl(:,2,k) = Rgl(:,2,k) / MassMol
                          
         ii = IAtmChain(1, 1,k)
         jj = IAtmChain(1,16,k)
         chainvector(:,1,k) = R(:,jj) - R(:,ii)

         ii = IAtmChain(2, 1,k)
         jj = IAtmChain(2,16,k)
         chainvector(:,2,k) = R(:,jj) - R(:,ii)

       end do

! ## <SP>
       do k = 1, NumLipid

         do l = 1, 2

           Sgl(:,l,k) = matmul( InvH, Rgl(:,l,k) )

         end do

       end do

       do k = 1, NumLipid/2

         do l = 1, 2

           do kk = k, NumLipid/2

             do ll = 1, 2

               if((k==kk).and.(l==ll)) cycle

                 Sij = Sgl(:,ll,kk) - Sgl(:,l,k)
                 Sij = Sij - nint(Sij)
                 Rij = matmul( H, Sij )
                 R2 = dot_product( Rij , Rij )

                 if(R2 < Rcutoff2) then

                   v1 = chainvector(:,l,k)
                   v2 = chainvector(:,ll,kk)
                   d1 = dot_product( v1, v1 )
                   d2 = dot_product( v2, v2 )
                   cs = dot_product( v1, v2 ) / sqrt( d1 * d2 )

                   corr = corr + cs

                   NumSample1 = NumSample1 + 1

                 end if

             end do

           end do

         end do

       end do

       do k = NumLipid/2+1, NumLipid

         do l = 1, 2

           do kk = k, NumLipid

             do ll = 1, 2

               if((k==kk).and.(l==ll)) cycle

                 Sij = Sgl(:,ll,kk) - Sgl(:,l,k)
                 Sij = Sij - nint(Sij)
                 Rij = matmul( H, Sij )
                 R2 = dot_product( Rij , Rij )

                 if(R2 < Rcutoff2) then

                   v1 = chainvector(:,l,k)
                   v2 = chainvector(:,ll,kk)
                   d1 = dot_product( v1, v1 )
                   d2 = dot_product( v2, v2 )
                   cs = dot_product( v1, v2 ) / sqrt( d1 * d2 )

                   corr = corr + cs

                   NumSample1 = NumSample1 + 1

                 end if

             end do

           end do

         end do

       end do
! ## <SP>

       write(7) sngl(Rgl)
       write(8) sngl(chainvector)

     end do

   end do

   open(1,file='CScorr.dat',status='unknown')

   write(1,*) corr / dble(NumSample1)

   close(1)

Contains

   subroutine NameOfChain

      AtomChain(1, 1) = 'C21'
      AtomChain(1, 2) = 'C22'
      AtomChain(1, 3) = 'C23'
      AtomChain(1, 4) = 'C24'
      AtomChain(1, 5) = 'C25'
      AtomChain(1, 6) = 'C26'
      AtomChain(1, 7) = 'C27'
      AtomChain(1, 8) = 'C28'
      AtomChain(1, 9) = 'C29'
      AtomChain(1,10) = 'C210'
      AtomChain(1,11) = 'C211'
      AtomChain(1,12) = 'C212'
      AtomChain(1,13) = 'C213'
      AtomChain(1,14) = 'C214'
      AtomChain(1,15) = 'C215'
      AtomChain(1,16) = 'C216'

      AtomChain(2, 1) = 'C31'
      AtomChain(2, 2) = 'C32'
      AtomChain(2, 3) = 'C33'
      AtomChain(2, 4) = 'C34'
      AtomChain(2, 5) = 'C35'
      AtomChain(2, 6) = 'C36'
      AtomChain(2, 7) = 'C37'
      AtomChain(2, 8) = 'C38'
      AtomChain(2, 9) = 'C39'
      AtomChain(2,10) = 'C310'
      AtomChain(2,11) = 'C311'
      AtomChain(2,12) = 'C312'
      AtomChain(2,13) = 'C313'
      AtomChain(2,14) = 'C314'
      AtomChain(2,15) = 'C315'
      AtomChain(2,16) = 'C316'

      if(LipidName == 'PhPC') then

        AtomChain(1,17) = 'C217'
        AtomChain(1,18) = 'C218'
        AtomChain(1,19) = 'C219'
        AtomChain(1,20) = 'C220'

        AtomChain(2,17) = 'C317'
        AtomChain(2,18) = 'C318'
        AtomChain(2,19) = 'C319'
        AtomChain(2,20) = 'C320'

      end if

    end subroutine NameOfChain

end subroutine NeighborChainCorr


!######################################################################
!######################################################################


subroutine PNcorr

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H, InvH
use CutoffParam, only : Rcutoff2
use AtomParam, only : MolName, AtomName

implicit none

integer :: i , j, k, l, ii, kk
integer :: Numa, Numb, nm, nmh
integer, dimension(:,:), allocatable :: ChPair
integer :: StepNumber, NumSample
real(8), dimension(3) :: Zij, Xij
real(8), dimension(:,:), allocatable :: VecPN, Sp
real(8), dimension(:), allocatable :: R_PN
real(8) :: Z1, Z2, cs, corr, Rc


   ii = 0

   do i = 1, NumSpec

     if(MolName(i)(1:4) == 'DPPC') then
       ii = ii + 1
       Nlipid = i
       Rc = 11.6d0
     else if(MolName(i)(1:4) == 'PhPC') then
       ii = ii + 1
       Nlipid = i
       Rc = 12.9d0
     end if

   end do

   Rcutoff2 = Rc*Rc

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no lipid molecule or more than 2 kinds of lipids'
     call Finalize
   end if

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Nlipid /= 1) then

     do i = 2 , Nlipid

       Numa = Numb
       Numb = Numb + NumMol(i) * NumAtm(i)

     end do

   end if

   l = Numa

   nm = NumMol(Nlipid)
   nmh = nm / 2


   allocate( ChPair(2,nm) )
   allocate( VecPN(3,nm) )
   allocate( Sp(3,nm) )
   allocate( R_PN(nm) )

   ChPair = 0

   do i = 1 , NumMol(Nlipid)

     do j = 1 , NumAtm(Nlipid)

       l = l + 1

       if(AtomName(l) == 'P') then
         ChPair(1,i) = l
       else if(AtomName(l) == 'N') then
         ChPair(2,i) = l
       end if

     end do

   end do

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       StepNumber = StepNumber + 1

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

!!       call CellTransform
       call Transform

       call InversMatrix(H,InvH)

       do k = 1 , nm

         call UnitVect(Xij,Z1,ChPair(1,k),ChPair(2,k))
         VecPN(:,k) = Xij
         R_PN(k)    = Z1
         Sp(:,k)    = matmul(InvH, R(:,ChPair(1,k)))

       end do

       do k = 1, nmh-1

         do kk = k+1, nmh

           Xij = Sp(:,k) - Sp(:,kk)
           Xij = Xij - nint(Xij)
           Zij = matmul(H,Xij)
           Z2 = Zij(1)*Zij(1) + Zij(2)*Zij(2)

           if(Z2 < Rcutoff2) then

             NumSample = NumSample + 1
             cs = dot_product(VecPN(:,k),VecPN(:,kk))
             corr = corr + cs

           end if

         end do

       end do

       do k = nmh + 1, nm - 1

         do kk = k+1, nm

           Xij = Sp(:,k) - Sp(:,kk)
           Xij = Xij - nint(Xij)
           Zij = matmul(H,Xij)
           Z2 = Zij(1)*Zij(1) + Zij(2)*Zij(2)

           if(Z2 < Rcutoff2) then

             NumSample = NumSample + 1
             cs = dot_product(VecPN(:,k),VecPN(:,kk))
             corr = corr + cs

           end if

         end do

       end do

     end do

   end do

   corr = corr / dble(NumSample)

   open(1,file='PNcorr.dat',form='unformatted',status='unknown')

   write(1,*) corr

   close(1)

contains

   subroutine UnitVect(Sij,R1,ii,jj)

   integer :: ii , jj
   real(8), dimension(3) :: Rij, Sij
   real(8) :: R2, R1

      Rij = R(:,jj) - R(:,ii)
      R2  = dot_product(Rij,Rij)
      R1  = sqrt( R2 )
      Sij = Rij / R1

   end subroutine UnitVect

end subroutine PNcorr


!######################################################################
!######################################################################


subroutine Szz_CG

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H
use AtomParam, only : MolName, AtomName

implicit none

integer :: ip
integer, dimension(:), allocatable :: Numa
integer, dimension(:), allocatable :: Nslipid
real(8), dimension(:,:,:), allocatable :: SzzData
integer :: MaxCarbon, MaxChain, MaxMOL
integer, dimension(:), allocatable :: Num_Chain
integer, dimension(:,:), allocatable :: NumCarbon
integer, dimension(:,:,:,:), allocatable :: CarbonNumber
integer :: NumLip, totalstep
integer :: i, j, k, l, ii, jj, i1, i2
real(8), dimension(3) :: R12
real(8) :: R_12, cs1
character(len=4), dimension(:), allocatable :: LipidName

   open(1,file='./Analy/Szz_CG.dat',status='unknown')

   do j = 1, 2

     ii = 0

     do i = 1, NumSpec

       if((MolName(i)(1:4) == 'DPPC') .or. &
       &  (MolName(i)(1:4) == 'DMPC') .or. &
       &  (MolName(i)(1:4) == 'DOPC') .or. &
       &  (MolName(i)(1:4) == 'POPC') .or. &
       &  (MolName(i)(1:4) == 'POPE') .or. &
       &  (MolName(i)(1:3) == 'SDS' ) .or. &
       &  (MolName(i)(1:4) == 'CTAB') .or. &
       &  (MolName(i)(1:4) == 'DHDA')) then
         ii = ii + 1
         if(j==2) then
           Nslipid(ii) = i
           LipidName(ii) = MolName(i)(1:4)
         end if
       end if

     end do

     if(j==1) then
       if( ii == 0 ) then
         write(*,*) 'ERROR : there is no lipid molecule'
         call Finalize
       end if
       allocate( Numa(ii), Nslipid(ii), LipidName(ii) )
       allocate( Num_Chain(ii) )
       NumLip = ii
     end if

   end do

   Numa(:) = 0

   do j = 1, NumLip
     do i = 1 , Nslipid(j) - 1
       Numa(j) = Numa(j) + NumMol(i) * NumAtm(i)
     end do
   end do

   do j = 1, NumLip
     Num_Chain(j) = 2
     if( LipidName(j) == 'TBPC' ) Num_Chain(j) = 3
     if( LipidName(j) == 'SDS'  ) Num_Chain(j) = 1
     if( LipidName(j) == 'CTAB' ) Num_Chain(j) = 1
   end do

   MaxChain = 0
   do j = 1, NumLip
     if(Num_Chain(j)>MaxChain) MaxChain = Num_Chain(j)
   end do

   allocate( NumCarbon(MaxChain,NumLip) )

   do j = 1, NumLip
     if( LipidName(j) == 'DPPC' ) then
       NumCarbon(:,j) = 6
     else if( LipidName(j) == 'DMPC' ) then
       NumCarbon(:,j) = 5
     else if((LipidName(j) == 'POPC').or.(LipidName(j) == 'POPE')) then
       NumCarbon(1,j) = 7
       NumCarbon(2,j) = 6
     else if( LipidName(j) == 'DOPC' ) then
       NumCarbon(:,j) = 7
     else if( LipidName(j) == 'SDS' ) then
       NumCarbon(:,j) = 5
     else if((LipidName(j) == 'CTAB').or.(LipidName(j) == 'DHDA')) then
       NumCarbon(:,j) = 6
     end if
   end do

   MaxCarbon = 0
   do j = 1, NumLip
     do i = 1, Num_Chain(j)
       if(NumCarbon(i,j)>MaxCarbon) MaxCarbon = NumCarbon(i,j)
     end do
   end do

   MaxMOL = 0
   do j = 1, NumLip
     i = Nslipid(j)
     if(NumMol(i)>MaxMOL) MaxMOL = NumMol(i)
   end do

   allocate( CarbonNumber(MaxCarbon,MaxChain,MaxMOL,NumLip) )

   allocate( SzzData(MaxCarbon,MaxChain,NumLip) )

   call define_Szz_CG

   SzzData = 0.d0
   totalstep = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     totalstep = totalstep + NTrjStep(i)

     do j = 1 , NTrjStep(i)

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

       do ip = 1, NumLip

       do k = 1 , NumMol(Nslipid(ip))
         do ii = 1 , Num_Chain(ip)
           do l = 1 , NumCarbon(ii,ip) - 1

             i1 = CarbonNumber(l,ii,k,ip)
             i2 = CarbonNumber(l+1,ii,k,ip)

             R12 = R(:,i2) - R(:,i1)
             R_12 = sqrt( dot_product( R12, R12 ) )
             cs1 = R12(3) / R_12
             SzzData(l,ii,ip) = SzzData(l,ii,ip) + cs1*cs1

           end do
         end do
       end do

       end do

     end do

   end do

   do j = 1, NumLip
   do ii = 1, Num_Chain(j)
     do l = 1, NumCarbon(ii,j) - 1
       SzzData(l,ii,j) = SzzData(l,ii,j) / dble(TotalStep*NumMol(Nslipid(j)))
     end do
   end do
   end do

   SzzData = 0.5d0 * ( 3.d0*SzzData - 1.d0 )

   do k = 1, NumLip

     write(1,'(3a)') '### Szz of the acyl chain of ',LipidName(k),' ###'

     do i = 1, Num_Chain(k)

       write(1,'(a,i1,a)') '#### Chain ',i,' ####'

       do j = 1, NumCarbon(i,k) - 1

         write(1,'(i8,3f13.5)') j, SzzData(j,i,k)

       end do

     end do

   end do

   close(1)

!##############################

contains

   subroutine define_Szz_CG

   implicit none

   character(len=4), dimension(MaxCarbon,MaxChain,NumLip) :: CarbonName
   integer :: Numb1, ll

   do k = 1, NumLip

     if( LipidName(k)=='SDS' ) then

       write(CarbonName(1,1,k),'(a)') 'SO4'

       do i = 2, MaxCarbon
         write(CarbonName(i,1,k),'(a1,i1)') 'C',i-1
       end do

     else if( LipidName(k)=='CTAB' ) then

       write(CarbonName(1,1,k),'(a)') 'NC4'

       do i = 2, MaxCarbon
         write(CarbonName(i,1,k),'(a1,i1)') 'C',i-1
       end do

     else if( LipidName(k)=='DHDA' ) then

       write(CarbonName(1,1,k),'(a)') 'NC4'
       write(CarbonName(1,2,k),'(a)') 'NC4'

       do i = 2, MaxCarbon
         write(CarbonName(i,1,k),'(a2,i1)') 'C1',i-1
         write(CarbonName(i,2,k),'(a2,i1)') 'C2',i-1
       end do

     else

       write(CarbonName(1,1,k),'(a)') 'EST1'
       write(CarbonName(1,2,k),'(a)') 'EST2'

       do i = 2, MaxCarbon
         write(CarbonName(i,1,k),'(a2,i1)') 'C1',i-1
         write(CarbonName(i,2,k),'(a2,i1)') 'C2',i-1
       end do

     end if

     do i = 1 , NumMol(Nslipid(k))
       Numb1 = Numa(k) + (i-1)*NumAtm(Nslipid(k))
       do j = 1 , NumAtm(Nslipid(k))
         ll = Numb1 + j

         do ii = 1 , Num_Chain(k)
           do jj = 1 , NumCarbon(ii,k)
             if(AtomName(ll) == CarbonName(jj,ii,k)) then
               CarbonNumber(jj,ii,i,k) = ll
             end if
           end do
         end do

       end do
     end do

   end do

   end subroutine define_Szz_CG


end subroutine Szz_CG
