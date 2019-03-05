! ############################
! ## SUBROUTINE LIST 
! ## -- CoarseningData 
! ## -- C_GC_def
! ## -- RadialC 
! ## -- Write_CGSites 
! ############################


module CGchW
   integer :: Kw, Kcomp, IRcut, NumTypes
   real(8), dimension(:,:), allocatable :: gr_cgb, EneNonBond, EneTotalNonBond
   real(8), dimension(:,:), allocatable :: cc_cgb
   real(8), dimension(:,:,:,:), allocatable :: Rsite_onGroup
   integer, dimension(:,:), allocatable :: NAtomGroup
   integer, dimension(:), allocatable :: NtypeGroup
   integer, dimension(:), allocatable :: NinTypeGroup
   integer :: NumGroup
   integer, dimension(:), allocatable :: GroupNumber, NumAtm_Group
   real(8), dimension(:), allocatable :: Mass_Group
   real(8), dimension(:,:,:), allocatable :: Rcom_Group
end module CGchW


!######################################################################
!######################################################################


subroutine CG_Chain_W

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use CellListMethod
use CGchW
use CGanalyW
use ParamAnalyze, only :NJobs, NTrjStep, Interval
use UnitExParam, only : pi, cvol
use CellParam, only : H, InvH, CellL, InvCL, Volume
use AtomParam, only : ResidName, Mass

implicit none

character(len=72) :: FileName
integer :: i, j, k, l, ii, Numa, NumF
integer :: imol
real(8), dimension(3) :: Rg

real(8) :: Vol, det, qk, gg, drr, gn
real(8) :: cr, en, ett
integer :: NumGR, Ncell_pre
real(8), parameter :: Rth = 5.d0
external det

   ii = 1
   Kw = 0
   do i = 1 , NumSpec
     if(ResidName(ii)=='TIP3') then
       Kw = i
       NumIni = ii - 1
     end if
     ii = ii + NumMol(i)*NumAtm(i)
   end do

   call C_GC_def(Numa)

   NinTypeGroup(:) = 0
   do i = 1, NumGroup
     j = NTypeGroup(i)
     NinTypeGroup(j) = NinTypeGroup(j) + 1
   end do

   NumGR = NumTypes

   MassW = 0.d0
   do i = NumIni+1, NumIni+3
     MassW = MassW + Mass(i)
   end do
   NumW = NumMol(Kw)
   allocate( Rw(3,NumW) )
   allocate( number_pair(NumW) )
   allocate( id_pair(NumW,50) )
   allocate( WCount(NumW) )
   allocate( RWCount(NumW,10) )
   allocate( NextP(NumW) )

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

       if( mod(j,Interval(i)) == 0 ) then

         NumF = NumF + 1

         print *, 'step=', NumF

         call CheckCell
         call InversMatrix(H,InvH)
         do k = 1, 3
           CellL(k) = H(k,k)
           InvCL(k) = InvH(k,k)
         end do
         call PBC

         call CalcCG_COM(Numa)

         do imol = 1, NumW
           Rg = 0.d0
           do l = 1, 3
             ii = (imol-1)*3 + l + NumIni
             Rg(:) = Rg(:) + Mass(ii) * R(:,ii)
           end do
           Rw(:,imol) = Rg(:) / MassW
         end do

         Ndiv(:) = int( CellL(:) / Rth )

         Ncell= Ndiv(1) * Ndiv(2) * Ndiv(3)

         if(NumF == 1) then
           Maps = 26*Ncell
           allocate( Head(Ncell) )
           allocate( Map(Maps) )
           call MappingAll ! Mapping neighboring cells
         else if(Ncell/=Ncell_pre) then
           Maps = 26*Ncell
           deallocate( Head, Map )
           allocate( Head(Ncell) )
           allocate( Map(Maps) )
           call MappingAll ! Mapping neighboring cells
         end if

         call NeighborPair(Rth) ! pick up neighboring water

         call RadialCG_blobW(Numa)

         Volume = det(H)
         Vol = Vol + Volume

         Ncell_pre = Ncell

       end if

     end do

   end do

   Vol = Vol / dble(NumF)

   do i = 1, NumTypes

     write(FileName,'(a,i1,a)') './Analy/Ene_Type',i,'--blobW.dat'

     open(1,file = trim(FileName), status='unknown')

     cr = 0.d0

     do j = 1 , IRcut

       drr = (j-0.5d0) * 0.1d0
       qk = Vol / (NumMol(Kw)/3 * 4.d0 * pi * drr * drr * 0.1d0)
       gn = gr_cgb(i,j) / dble(NumF*NinTypeGroup(i)*NumMol(Kcomp))
       gg = gn * qk
       cr = cr + gn
       en = EneNonBond(i,j) / cc_cgb(i,j) * cvol
       ett = EneTotalNonBond(i,j) / cc_cgb(i,j) * cvol + 6.1d0
!      - ( - 6.1 kcal/mol) subtrancion intra-blob energy

       write(1,'(f8.3,2f12.6,2d18.10)') drr, gg, cr, en, ett

     end do

     close(1)

   end do

Contains

   subroutine CheckCell

       if((H(1,2)/=0.).or.(H(1,3)/=0.).or. &
       &  (H(2,1)/=0.).or.(H(2,3)/=0.).or. &
       &  (H(3,1)/=0.).or.(H(3,2)/=0.)) then
         write(*,*) 'error : this analysis is just for the ortholombic cell'
         call Finalize
       end if

   end subroutine CheckCell

end subroutine CG_Chain_W


!######################################################################
!######################################################################


subroutine RadialCG_blobW(Numa)

use Numbers, only : NumMol
use CGchW, only : NumGroup, Rcom_Group, Kcomp

implicit none

integer :: Numa
integer :: i, j
real(8), dimension(3) :: Rg

   do i = 1, NumGroup
     do j = 1 , NumMol(Kcomp)
       Rg(:) = Rcom_Group(:,i,j)
       call find_blobW(Rg,i,j,Numa)
     end do
   end do

end subroutine RadialCG_blobW


!######################################################################
!######################################################################


subroutine find_blobW(R_CGi,Igrp,Imol,Numa)

use CellParam, only : CellL, InvCL
use CutoffParam, only : Rcutoff2
use CGchW, only : NtypeGroup, gr_cgb
use CGanalyW, only : NumW, Rw, number_pair, id_pair, WCount, RWCount

implicit none

integer :: l, i, j, k, itrial, itype, ll, ir
integer :: iq, ip, Igrp,Imol,Numa
real(8), dimension(3) :: R_CGi, Rij, R_CGw_j
integer, dimension(3) :: Nij
real(8) :: R2, ranf, subt, plus
logical, dimension(NumW) :: Qdone
logical :: Qgoodpair
external ranf

   Qdone(:) = .False.
   WCount(:) = 0
   RWCount(:,:) = 0

   itype = NTypeGroup(Igrp)

   do l = 1, NumW

! ## pick up the first mol

     if(Qdone(l)) cycle
     i = l
     Rij(:) = Rw(:,i) - R_CGi(:)
     call MImageO(Rij,CellL,InvCL,R2,Nij)
     if(R2 > Rcutoff2+100.d0) then
       Qdone(l)=.True.
       cycle
     end if

#ifdef DEBUG
     print *, 'pick up l=',l
#endif

     Qgoodpair = .False.
     itrial = 0

     do while(.not.Qgoodpair)

       ip = int(ranf()*number_pair(l)) + 1

       j = id_pair(l,ip)

       itrial = itrial + 1

       iq = ip
       do while (iq==ip)
         iq = int(ranf()*number_pair(l)) + 1
       end do

       k = id_pair(l,iq)

       if(itrial > 1000) then
         write(*,*) 'cannot find a reasonable pair for ',l
         exit
       end if

       call find_rg(R_CGw_j,i,j,k)
       call check_blob(Qgoodpair,R_CGw_j,i,j,k)

     end do

#ifdef DEBUG
     print *, 'pick up blob = ', i,j,k
#endif

     Qdone(i) = .True.
     Qdone(j) = .True.
     Qdone(k) = .True.

     call Ene_CG_blobW(R_CGi, R_CGw_j, Igrp, Imol, i, j, k, Numa)

   end do

   subt = - 1.d0 / 3.d0
   do l = 1, NumW
#ifdef DEBUG
     print *, l, WCount(l)
#endif
     if(WCount(l) <= 1) cycle
     plus = 1.d0/3.d0/dble(WCount(l)) + subt
     do ll = 1, WCount(l)
       ir = RWCount(l,ll)
       gr_cgb(itype,ir) = gr_cgb(itype,ir) + plus
     end do
   end do

end subroutine find_blobW


!######################################################################
!######################################################################


subroutine CalcCG_COM(Numa)

use Configuration, only : R
use Numbers, only : NumMol, NumAtm
use AtomParam , only : Mass
use CGchW, only : Rcom_Group, NumGroup, Rsite_onGroup, NAtomGroup, &
&   NumAtm_Group, Kcomp, GroupNumber, Mass_Group

implicit none

integer :: k, ii, ll, jj, l, Numa

   do k = 1 , NumMol(Kcomp)

     ii = Numa + (k-1)*NumAtm(Kcomp)

! ##  Center of Mass

     Rcom_Group(:,:,k) = 0.d0

     do l = 1 , NumAtm(Kcomp)
       ll = GroupNumber(l)
       jj = ii + l
       Rcom_Group(:,ll,k) = Rcom_Group(:,ll,k) + R(:,jj) * Mass(jj)
     end do

     do l = 1 , NumGroup
       Rcom_Group(:,l,k) = Rcom_Group(:,l,k) / Mass_Group(l)
     end do

! ## internal coordinates for each CG site

     do l = 1 , NumGroup
       do ll = 1, NumAtm_Group(l)
         jj = ii + NAtomGroup(ll,l)
         Rsite_onGroup(:,ll,l,k) = R(:,jj) - Rcom_Group(:,l,k)
       end do
     end do

   end do

end subroutine CalcCG_COM


!######################################################################
!######################################################################


subroutine C_GC_def(Numa)

use Numbers, only : NumSpec, NumMol, NumAtm
use AtomParam, only : Mass
use CutoffParam, only : Rcutoff2
use CGchW

implicit none

integer :: i, ii, Numa
integer :: j, k, MaxNumAtm_Group, nn
character(len=72) :: String1, String

   open(21,file='CGmodel.data',status='unknown')

! ## this number is used to identify the numbering of molecule to be analyzed 

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     select case(nn)

     case(1)
       read(String,*) Kcomp

! ## to get 'Numa'

       ii = 1

       do i = 1 , NumSpec

         if(i == Kcomp) then
           Numa = ii - 1
         end if

         ii = ii + NumMol(i) * NumAtm(i)

       end do

       allocate( GroupNumber(NumAtm(Kcomp)) ) ! Group Number ( atom number ) within a molecule 

! ## Number of Group in a molecule, Number of Group Types in a molecule 

     case(2)

       read(String,*) NumGroup, NumTypes

       allocate(NumAtm_Group(NumGroup)) ! ## Number of atoms in a defined group
       allocate(Mass_Group(NumGroup))   ! ## Mass of a defined group 
       allocate(Rcom_Group(3,NumGroup,NumMol(Kcomp))) ! ## Center of mass of a defined group

       allocate( NTypeGroup(NumGroup) )   ! ## Type number of a group in a molecule
       allocate( NinTypeGroup(NumTypes) ) ! ## Number of groups of a type

       exit

     end select

   end do

! ## read group number of an atom in a molecule

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     read(String,*) j, GroupNumber(j)

     if(nn == NumAtm(Kcomp)) then
       if(j /= NumAtm(Kcomp)) then
         write(*,*) 'ERROR : reading GroupNumber'
         call Finalize
       end if
       exit
     end if

   end do

! ## read type number of a group

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     read(String,*) j, NTypeGroup(j)

     if(nn == NumGroup) then
       if(j /= NumGroup) then
         write(*,*) 'ERROR : reading NTypeGroup'
         call Finalize
       end if
       exit
     end if

   end do

   close(21)

! ## Number of type pairs

   IRcut = int( sqrt(Rcutoff2) * 10.d0 ) + 1

   allocate( EneNonBond(NumTypes,IRcut) )
   allocate( EneTotalNonBond(NumTypes,IRcut) )
   allocate( gr_cgb(NumTypes,IRcut) )
   allocate( cc_cgb(NumTypes,IRcut) )

   EneNonBond = 0.d0
   EneTotalNonBond = 0.d0
   gr_cgb = 0.d0
   cc_cgb = 0.d0

! ## Number of atoms in a group, Mass of a group

   NumAtm_Group = 0
   Mass_Group = 0.d0

   do i = 1, NumAtm(Kcomp)

     j = Numa + i
     k = GroupNumber(i)

     NumAtm_Group(k) = NumAtm_Group(k) + 1
     Mass_Group(k)   = Mass_Group(k)   + Mass(j)

   end do

! ## Maximum value of the number of atoms in a group

   MaxNumAtm_Group = 0

   do i = 1, NumGroup

     if(NumAtm_Group(i) > MaxNumAtm_Group) then
       MaxNumAtm_Group = NumAtm_Group(i)
     end if

   end do

! ## atom number of i-th atom in a group 

   allocate( NAtomGroup(MaxNumAtm_Group,NumGroup) )

   NumAtm_Group = 0

   do i = 1, NumAtm(Kcomp)

     k = GroupNumber(i)

     NumAtm_Group(k) = NumAtm_Group(k) + 1
     NAtomGroup(NumAtm_Group(k),k) = i

   end do

   allocate( Rsite_onGroup(3,MaxNumAtm_Group,NumGroup,NumMol(Kcomp)) )

end subroutine C_GC_def


!######################################################################
!######################################################################


subroutine Ene_CG_blobW(RcomCGsite, RcomWblob, Igrp, Imol, iw1, iw2, iw3, Numa)

use Configuration, only : R
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellL, InvCL
use CommonBlocks, only : QSwitch
use Numbers, only : NumAtm
use NonbondParam, only : Charge, Rminh, EpsLJ
use CGchW
use CGanalyW, only : WCount, RWCount, NumIni

implicit none

real(8), dimension(3) :: RcomCGsite, RcomWblob
integer :: Igrp, Imol, iw1, iw2, iw3, Numa
integer :: i, j, n1, j1, k1, ii, jj
real(8), dimension(3) :: Rij, RGij
integer, dimension(3) :: Nij, Nwlist
real(8) :: R2, R1, Upot, fk, Uint00
integer :: nnn, itype
real(8), dimension(3,3,3) :: RwJ ! (dim,atm,mol)
real(8) :: cf, Sgm, Sgm2, Eps, InvR2, SR2, SR6, SR12, ek, fk1

   Nwlist(1) = iw1
   Nwlist(2) = iw2
   Nwlist(3) = iw3

   itype = NTypeGroup(Igrp)

   do n1 = 1, 3 !mol
     k1 = Nwlist(n1)
     do j1 = 1, 3 !atom
       j = (k1-1)*3+j1+NumIni
       RwJ(:,j1,n1) = R(:,j) - RcomWblob(:)
       RwJ(:,j1,n1) = RwJ(:,j1,n1) - CellL(:)*nint(RwJ(:,j1,n1)*InvCL(:))
!       R2 = dot_product(RwJ(:,j1,n1),RwJ(:,j1,n1))
!       if(R2>20.) write(*,*) 'WARNING:too large distance from com for blob J', sqrt(R2)
     end do
   end do

   Rij = RcomCGsite - RcomWblob
   call MImageO(Rij,CellL,InvCL,R2,Nij)
   if(R2 > Rcutoff2) Return

   RGij = Rij

   R1 = sqrt(R2)
   nnn = int(10.d0*R1) + 1
   gr_cgb(itype,nnn) = gr_cgb(itype,nnn) + 1.d0
   cc_cgb(itype,nnn) = cc_cgb(itype,nnn) + 1.d0

   do i = 1, 3
     j = Nwlist(i)
     WCount(j) = WCount(j) + 1
     RWCount(j,WCount(j)) = nnn
   end do

   Upot = 0.d0
   do ii = 1, NumAtm_Group(Igrp)
     jj = NAtomGroup(ii,Igrp)
     i  = jj + (Imol-1)*NumAtm(Kcomp) + Numa
     do n1 = 1, 3 !mol
       k1 = Nwlist(n1)
       do j1 = 1, 3 !atom
         j = (k1-1)*3+j1+NumIni
         Rij(:) = Rsite_onGroup(:,ii,Igrp,Imol) - RwJ(:,j1,n1) + RGij(:)
         R2 = dot_product(Rij,Rij)

         Sgm   = Rminh(i) + Rminh(j)
         Sgm2  = Sgm * Sgm
         Eps   = EpsLJ(i) * EpsLJ(j)
         InvR2 = 1.d0 / R2

         SR2  = Sgm2 * InvR2                    !(sigma/r)^2
         SR6  = SR2 * SR2 * SR2                 !         ^6
         SR12 = SR6 * SR6                       !         ^12
         ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
   if(QSwitch) call SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

         cf = Charge(i) * Charge(j)
         fk1 = cf * sqrt(InvR2)

         Upot = Upot + (ek + fk1)

       end do

     end do
   end do

   EneNonBond(itype,nnn) = EneNonBond(itype,nnn) + Upot

   call cal_internal(Uint00,iw1,iw2,iw3)

   EneTotalNonBond(itype,nnn) = EneTotalNonBond(itype,nnn) + Upot + Uint00

end subroutine Ene_CG_blobW
