! ############################
! ## SUBROUTINE LIST 
! ############################

module CGanalyW
   integer :: NumW, NumIni
   real(8) :: MassW
   real(8), dimension(:,:), allocatable :: Rw
   real(8), dimension(:), allocatable :: Uintra
   real(8), dimension(:), allocatable :: gra, pota
   integer, dimension(:), allocatable :: number_pair
   integer, dimension(:,:), allocatable :: id_pair
   real(8), dimension(:), allocatable :: gbb, grb, Ene_blob, Ene_total_blob
   integer, dimension(:), allocatable :: WCount
   integer, dimension(:,:), allocatable :: RWCount
   real(8) :: subUin
end module CGanalyW


!######################################################################
!######################################################################


subroutine CoarseGrainW

use CGanalyW
use Configuration, only : R
use Numbers, only : NumSpec, NumMol, NumAtm
use AtomParam, only : MolName, Mass
use CommonBlocks, only : QMaster
use CellParam, only : H, InvH, CellL, InvCL
use CutoffParam, only : Rcutoff2
use ParamAnalyze, only : NJobs, NTrjStep, IRcut
use CellListMethod
use UnitExParam, only : cvol, pi

implicit none

integer :: i, ii, j, k, l
integer :: ix, iy, iz, imol, im, jm, km
integer :: totalstep, Nwater, Ncell_pre
real(8), dimension(3) :: Rg, Ro
real(8), parameter :: Rth = 5.d0
real(8) :: sUintra, Vol, det, drr, qk, gg, qkb, ggb, tot
logical :: Qgoodpair
real(8) :: Umax, Umin, searchmax, searchmin, UU
external det, searchmax, searchmin

   ii = 0
   do i = 1, NumSpec
     if(MolName(i)(1:3) == 'TIP') then
       ii = ii + 1
       Nwater = i
     end if
   end do

   if( ii /= 1 ) then
     write(*,*) 'ERROR : there is no TIP3 water molecule'
     call Finalize
   end if

   NumIni = 0

   if(Nwater /= 1) then
     do i = 2 , Nwater
       NumIni = NumIni + NumMol(i-1) * NumAtm(i-1)
     end do
   end if

   MassW = 0.d0
   do i = NumIni+1, NumIni+3
     MassW = MassW + Mass(i)
   end do

   NumW = NumMol(Nwater)

   allocate( Rw(3,NumW) )
   allocate( number_pair(NumW) )
   allocate( id_pair(NumW,50) )
   allocate( WCount(NumW) )
   allocate( RWCount(NumW,10) )

   IRcut = int(sqrt(Rcutoff2)*10.d0) + 1
   allocate( gra(IRcut) )
   allocate( pota(IRcut) )
   allocate( gbb(IRcut) )
   allocate( grb(IRcut) )
   allocate( Ene_blob(IRcut) )
   allocate( Ene_total_blob(IRcut) )

   sUintra = 0.d0
   subUin = 0.d0
   totalstep = 0
   Vol = 0
   gra = 0.d0
   pota = 0.d0
   grb = 0.d0
   gbb = 0.d0
   Ene_blob = 0.d0
   Ene_total_blob = 0.d0

   open(31,file='ene_blob_intra.dat',form='unformatted',status='unknown')
   open(32,file='ene_blob_intra_bound.dat',form='formatted',status='unknown')
   open(33,file='gr_blob_mol.dat',form='formatted',status='unknown')
   open(34,file='ene_blob_mol.dat',form='formatted',status='unknown')
   open(35,file='ene_blob_blob.dat',form='formatted',status='unknown')
   open(36,file='gr_blob_blob.dat',form='formatted',status='unknown')
   open(37,file='ene_total_blob_blob.dat',form='formatted',status='unknown')

   allocate( NextP(NumW) )

   do i = 1 , NJobs

     if(QMaster) call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       totalstep = totalstep + 1

       print *, 'step=',totalstep

!     -----------------------------
#ifdef MOLFILE
       if(QMaster) call Read_RTraj(i)
#else
       if(QMaster) call Read_RTraj
#endif
!     -----------------------------

       Vol = Vol + det(H)

       call CheckCell
       call BcastRH
       call InversMatrix(H,InvH)
       call PBC
       do k = 1, 3
         CellL(k) = H(k,k)
         InvCL(k) = InvH(k,k)
       end do

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

       if(totalstep == 1) then
         Maps = 26*Ncell
         allocate( Head(Ncell) )
         allocate( Map(Maps) )
         allocate( Uintra(Ncell) )
       else if(Ncell/=Ncell_pre) then
         Maps = 26*Ncell
         deallocate( Head, Map, Uintra )
         allocate( Head(Ncell) )
         allocate( Map(Maps) )
         allocate( Uintra(Ncell) )
       end if

       call MappingAll ! Mapping neighboring cells
       call NeighborPair(Rth) ! pick up neighboring water

! ## based upon the grid points

       do ix = 1, Ndiv(1)
       do iy = 1, Ndiv(2)
       do iz = 1, Ndiv(3)

         ii = (ix-1)*Ndiv(2)*Ndiv(3) + (iy-1)*Ndiv(3) + iz
         Ro(1) = ((ix-0.5d0)/Ndiv(1) - 0.5d0)*CellL(1)
         Ro(2) = ((iy-0.5d0)/Ndiv(2) - 0.5d0)*CellL(2)
         Ro(3) = ((iz-0.5d0)/Ndiv(3) - 0.5d0)*CellL(3)

#ifdef DEBUG
         print *, 'Ro'
         print *, 'Ro = ',Ro(:)
#endif

         call find_closest(Ro,im,jm,km) ! finding the closest three molecules i,j,k

#ifdef DEBUG
         print *, 'find_closest'
         print *, 'i = ',im
         print *, 'j = ',jm
         print *, 'k = ',km
#endif

         call find_rg(Rg,im,jm,km)      ! calculate the centers of mass

#ifdef DEBUG
         print *, 'find_rg'
         print *, 'Rg = ',Rg(:)
#endif

         call check_blob(Qgoodpair,Rg,im,jm,km)  ! check 

#ifdef DEBUG
         print *, 'check_blob'
         print *, 'good pair = ',Qgoodpair
         print *, 'i = ',im
         print *, 'j = ',jm
         print *, 'k = ',km
#endif

         if(.not.Qgoodpair) call search_possible_pair(Rg,im,jm,km) ! search better pair if needed 
#ifdef DEBUG
         print *, 'search_possible_pair'
         print *, 'good pair = ',Qgoodpair
         print *, 'i = ',im
         print *, 'j = ',jm
         print *, 'k = ',km
#endif

         call cal_internal(UU,im,jm,km)  ! calculate energies inside the blob
         Uintra(ii) = UU
         call cal_rdf_ene_a(Rg,im,jm,km)  ! calculate rdf and energy between water and the blob
         call cal_blob_blob(Rg,im,jm,km) !calculate energy between blobs

       end do
       end do
       end do

       Uintra(:) = Uintra(:) * cvol
       do ii = 1, Ncell
         sUintra = sUintra + Uintra(ii)
       end do
       Umax = searchmax(Uintra,Ncell)
       Umin = searchmin(Uintra,Ncell)
       write(31) Ncell
       write(31) Uintra(:)
       write(32,'(2e16.6)') Umax, Umin

       Ncell_pre = Ncell

     end do

   end do

   sUintra = sUintra / dble(Ncell*totalstep)

   write(32,'(d24.16)') sUintra

! ## rdf 

   Vol = Vol / dble(totalstep)

   do i = 1, IRcut

     drr = (i-0.5) * 0.1d0
     qk = Vol / (NumW * 4.d0 * pi * drr * drr * 0.1d0)

     gg = gra(i) * qk / dble(totalstep*Ncell)
     write(33,'(f8.3,f12.6)') drr, gg

     qkb = Vol / (NumW / 3.d0 * 4.d0 * pi * drr * drr * 0.1d0)

     ggb = grb(i) * qkb / dble(totalstep*Ncell)
     write(36,'(f8.3,f12.6)') drr, ggb

   end do

   pota(:) = pota(:) * cvol
   Ene_blob(:) = Ene_blob(:) * cvol
   Ene_total_blob(:) = Ene_total_blob(:) * cvol

   tot = 0.d0
   do i = 1, IRcut
     tot = tot + gbb(i)
   end do
   subUin = 0.5d0 * subUin * cvol / tot

   do i = 1, IRcut

     drr = (i-0.5) * 0.1d0

     if(gra(i) /= 0.) then
       pota(i) = pota(i) / gra(i)
       write(34,'(f8.3,d18.10)') drr, pota(i)
     else
       write(34,'(f8.3,d18.10)') drr, 0.
     end if

     if(gbb(i) /= 0.) then
       Ene_blob(i) = Ene_blob(i) / gbb(i)
       write(35,'(f8.3,d18.10)') drr, Ene_blob(i)
     else
       write(35,'(f8.3,d18.10)') drr, 0.
     end if

     if(gbb(i) /= 0.) then
       Ene_total_blob(i) = Ene_total_blob(i) / gbb(i) - subUin * 2.d0
       write(37,'(f8.3,d18.10)') drr, Ene_total_blob(i)
     else
       write(37,'(f8.3,d18.10)') drr, 0.
     end if

   end do

   write(32,'(d24.16)') subUin

   close(31)
   close(32)
   close(33)
   close(34)
   close(35)
   close(36)
   close(37)

Contains

   subroutine CheckCell

       if((H(1,2)/=0.).or.(H(1,3)/=0.).or. &
       &  (H(2,1)/=0.).or.(H(2,3)/=0.).or. &
       &  (H(3,1)/=0.).or.(H(3,2)/=0.)) then
         write(*,*) 'error : this analysis is just for the ortholombic cell'
         call Finalize
       end if

   end subroutine CheckCell

end subroutine CoarseGrainW


!######################################################################
!######################################################################


function searchmax(a,num)

   integer :: num, ll
   real(8), dimension(num) :: a
   real(8) :: searchmax

   searchmax = a(1)
   do ll = 2, num
     if(a(ll)>searchmax) searchmax = a(ll)
   end do

end function searchmax


!######################################################################
!######################################################################


function searchmin(a,num)

   integer :: num, ll
   real(8), dimension(num) :: a
   real(8) :: searchmin

   searchmin = a(1)
   do ll = 2, num
     if(a(ll)<searchmin) searchmin = a(ll)
   end do

end function searchmin


!######################################################################
!######################################################################


subroutine NeighborPair(Rth)

use CGanalyW
use CellParam, only : CellL, InvCL
use CommonMPI, only : NProcs, MyRank
use CellListMethod

implicit none

real(8) :: Rth, Rth2
integer :: Nas, icell, i, j, k, jcell0, jcell, nabor
real(8), dimension(3) :: Ri, Rij
integer, dimension(3) :: Nij
real(8) :: R2

   Nas = NProcs - MyRank

   Rth2 = Rth*Rth

   Head = 0

   number_pair = 0
   id_pair = 0

   do k = 1, NumW
     icell = 1 + int((Rw(1,k) * InvCL(1) + 0.5d0 ) * Ndiv(1) )                   &
     &         + int((Rw(2,k) * InvCL(2) + 0.5d0 ) * Ndiv(2) ) * Ndiv(1)         &
     &         + int((Rw(3,k) * InvCL(3) + 0.5d0 ) * Ndiv(3) ) * Ndiv(1) * Ndiv(2)
     NextP(k) = Head(icell)
     Head(icell) = k
   end do

   do icell = Nas, Ncell, NProcs

     i = Head(icell)

     do while( i /= 0 )

       Ri = Rw(:,i)
       j = Head(icell)

       do while( j /= 0 )
         if(i==j) then
           j = NextP(j)
           cycle
         end if

         Rij = Ri - Rw(:,j)
         R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
         if(R2 < Rth2) then
           number_pair(i) = number_pair(i) + 1
           id_pair(i,number_pair(i)) = j
         end if

         j = NextP(j)

       end do

       jcell0 = 26 * ( icell - 1 )

       do nabor = 1 , 26

         jcell = Map( jcell0 + nabor )

         if(jcell > 0) then

           j = Head(jcell)

           do while( j > 0 )

             Rij = Ri - Rw(:,j)
             call MImageO(Rij,CellL,InvCL,R2,Nij)
             if(R2 < Rth2) then
               number_pair(i) = number_pair(i) + 1
               id_pair(i,number_pair(i)) = j
             end if

             j = NextP(j)

           end do

         end if

       end do

       i = NextP(i)

     end do

   end do

   call SumList( NumW*51, number_pair, id_pair )

   do i = 1, NumW
     if(number_pair(i)<2) then
       write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       write(*,*) 'warning : no pair was found for water',i
       write(*,*) 'number =',number_pair(i)
       write(*,*) 'recalculing the nearest neighbor...'
       write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       call search_pair_all(i)
     else if(number_pair(i)>50) then
       write(*,*) 'error : too many pairs were found for water',i
       write(*,*) 'number =',number_pair(i)
       call Finalize
     end if
   end do

end subroutine NeighborPair


!######################################################################
!######################################################################


subroutine search_pair_all(ii)

use CGanalyW, only : Rw, NumW, number_pair, id_pair
use CellParam, only : CellL, InvCL

implicit none

integer :: ii, idup
integer :: j
real(8), dimension(3) :: Rij
real(8) :: R1st, R2nd, R2
integer :: N1st, N2nd

   if(number_pair(ii)==1) then
     idup = id_pair(ii,1)
   else
     idup = 0
   end if

   R1st = 400.d0
   R2nd = 400.d0

   do j = 1, NumW
     if(j==ii) cycle
     Rij = Rw(:,j) - Rw(:,ii)
     Rij = Rij - nint(InvCL*Rij) * CellL
     R2  = dot_product(Rij,Rij)
     if(R2 < R1st) then
       R2nd = R1st
       R1st = R2
       N2nd = N1st
       N1st = j
     else if(R2 < R2nd) then
       R2nd = R2
       N2nd = j
     end if
   end do

   if((idup/=0).and.(idup/=N1st)) then
     write(*,*) 'possible error for neighboring list'
     call Finalize
   end if

   number_pair(ii) = 2
   id_pair(ii,1) = N1st
   id_pair(ii,2) = N2nd

end subroutine search_pair_all


!######################################################################
!######################################################################


subroutine check_blob(Qgoodpair,Rg,i,j,k)

use CGanalyW
use CellParam, only : CellL, InvCL
use CellListMethod

implicit none

integer :: ii, jj, i, j, k
integer :: nabor, icell, jcell0, jcell
integer, dimension(3) :: mlist, Nij
real(8), dimension(3) :: Rg, Rij
real(8) :: R2, rmx2
logical :: Qgoodpair

   Qgoodpair = .True.

   mlist(1) = i
   mlist(2) = j
   mlist(3) = k

#ifdef DEBUG
   print *, '****In check_blob***'
   print '(a,3f8.4)', 'Rg=',Rg(:)
#endif
   rmx2 = 0.d0
   do ii = 1, 3
     jj = mlist(ii)
     Rij(:) = Rw(:,jj) - Rg(:)
     call MImageO(Rij,CellL,InvCL,R2,Nij)
#ifdef DEBUG
     print '(i5,4f8.4)', jj,sqrt(R2),Rw(:,jj)
#endif
     if(R2 > rmx2) rmx2 = R2
   end do

#ifdef DEBUG
   print *, 'rmax=',sqrt(rmx2)
#endif

   icell = 1 + int((Rg(1) * InvCL(1) + 0.5d0 ) * Ndiv(1) )                   &
   &         + int((Rg(2) * InvCL(2) + 0.5d0 ) * Ndiv(2) ) * Ndiv(1)         &
   &         + int((Rg(3) * InvCL(3) + 0.5d0 ) * Ndiv(3) ) * Ndiv(1) * Ndiv(2)

   jj = Head(icell)

#ifdef DEBUG
   print *, 'cell=',icell
#endif

   do while(jj /= 0)
     if((jj==i).or.(jj==j).or.(jj==k)) then
       jj = NextP(jj)
       cycle
     end if
     Rij = Rw(:,jj) - Rg(:)
     R2     = dot_product(Rij,Rij)
     if(R2 < rmx2) then
       Qgoodpair = .False.
       Return
     end if
     jj = NextP(jj)
   end do

   jcell0 = 26 * ( icell - 1 )
   do nabor = 1 , 26
     jcell = Map( jcell0 + nabor )
     if(jcell > 0) then
       jj = Head(jcell)
       do while ( jj /= 0 )
         if((jj==i).or.(jj==j).or.(jj==k)) then
           jj = NextP(jj)
           cycle
         end if
         Rij = Rg - Rw(:,jj)
         call MImageO(Rij,CellL,InvCL,R2,Nij)
         if(R2 < rmx2) then
           Qgoodpair = .False.
           Return
         end if
         jj = NextP(jj)
       end do
     end if
   end do

#ifdef DEBUG
   print *, 'Fin: check_blob'
#endif

end subroutine check_blob


!######################################################################
!######################################################################


subroutine find_closest(Ro,i,j,k)

use CGanalyW
use CellParam, only : CellL, InvCL
use CellListMethod

implicit none

integer :: i,j,k
integer :: icell, jj, jcell0, jcell, nabor
integer, dimension(3) :: Nminim, Nij
real(8), dimension(3) :: Ro, Rij
real(8), dimension(3) :: Rminim2
real(8) :: R2

   Rminim2(:) = 1000.d0
   Nminim(:)  = 0

   icell = 1 + int((Ro(1) * InvCL(1) + 0.5d0 ) * Ndiv(1) )                   &
   &         + int((Ro(2) * InvCL(2) + 0.5d0 ) * Ndiv(2) ) * Ndiv(1)         &
   &         + int((Ro(3) * InvCL(3) + 0.5d0 ) * Ndiv(3) ) * Ndiv(1) * Ndiv(2)

   jj = Head(icell)

   do while(jj > 0)
     Rij(:) = Rw(:,jj) - Ro(:)
     R2  = dot_product(Rij,Rij)
     call dis_judge
     jj = NextP(jj)
   end do

   jcell0 = 26 * ( icell - 1 )
   do nabor = 1 , 26
     jcell = Map( jcell0 + nabor )
     if(jcell > 0) then
       jj = Head(jcell)
       do while(jj > 0)
         Rij(:) = Ro(:) - Rw(:,jj)
         call MImageO(Rij,CellL,InvCL,R2,Nij)
         call dis_judge
         jj = NextP(jj)
       end do
     end if
   end do

#ifdef DEBUG
   print *, 'i=',Nminim(1),Rminim2(1)
   print *, 'j=',Nminim(2),Rminim2(2)
   print *, 'k=',Nminim(3),Rminim2(3)
#endif

   i = Nminim(1)
   j = Nminim(2)
   k = Nminim(3)

Contains

   subroutine dis_judge

     if(R2 < Rminim2(1)) then
       Rminim2(3) = Rminim2(2)
       Rminim2(2) = Rminim2(1)
       Rminim2(1) = R2
       Nminim(3)  = Nminim(2)
       Nminim(2)  = Nminim(1)
       Nminim(1)  = jj
     else if(R2 < Rminim2(2)) then
       Rminim2(3) = Rminim2(2)
       Rminim2(2) = R2
       Nminim(3)  = Nminim(2)
       Nminim(2)  = jj
     else if(R2 < Rminim2(3)) then
       Rminim2(3) = R2
       Nminim(3)  = jj
     end if

   end subroutine dis_judge

end subroutine find_closest


!######################################################################
!######################################################################


subroutine find_rg(Rg,i,j,k)

use CellParam, only : CellL, InvCL
use CGanalyW, only : Rw

implicit none

integer :: i, j, k, l
real(8), dimension(3) :: Rg, Rij
real(8), dimension(3,3) :: Rshift
integer, dimension(3) :: Nij
real(8) :: R2

   Rshift(:,1) = Rw(:,i)

   Rij(:) = Rw(:,j) - Rshift(:,1)
   call MImageO(Rij,CellL,InvCL,R2,Nij)
   Rshift(:,2) = Rij(:) + Rshift(:,1)

   Rij(:) = Rw(:,k) - Rshift(:,1)
   call MImageO(Rij,CellL,InvCL,R2,Nij)
   Rshift(:,3) = Rij(:) + Rshift(:,1)

   Rg(:) = 0.d0
   do l = 1, 3
     Rg(:) = Rg(:) + Rshift(:,l)
   end do
   Rg(:) = Rg(:) / 3.d0

   Rg(:) = Rg(:) - nint( InvCL(:) * Rg(:) )* CellL(:)

end subroutine find_rg


!######################################################################
!######################################################################


subroutine cal_internal(Upot,ii,jj,kk)

use CGanalyW
use Configuration, only : R
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellL, InvCL
use CommonBlocks, only : QSwitch
use NonbondParam, only : Rminh, EpsLJ, Charge

implicit none

integer :: i, j, l
integer :: ii, jj, kk, im1, im2, ia1, ia2
integer, dimension(2,3) :: kpair
real(8), dimension(3) :: Rij
integer, dimension(3) :: Nij
real(8) :: Upot, fk, ek, fk1, R2
real(8) :: Sgm, Sgm2, Eps, InvR2, SR2, SR6, SR12, cf

   kpair(1,1) = ii; kpair(2,1) = jj
   kpair(1,2) = jj; kpair(2,2) = kk
   kpair(1,3) = kk; kpair(2,3) = ii

   Upot = 0.d0

   fk = 0.d0

   do l = 1, 3
     im1 = kpair(1,l)
     im2 = kpair(2,l)
     do ia1 = 1, 3
       i = (im1-1)*3 + ia1 + NumIni
       do ia2 = 1, 3
         j = (im2-1)*3 + ia2 + NumIni
         Rij(:) = R(:,i) - R(:,j)
         call MImageO(Rij,CellL,InvCL,R2,Nij)

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


end subroutine cal_internal


!######################################################################
!######################################################################


subroutine cal_rdf_ene_a(Rg,n1,n2,n3)

use CGanalyW
use Configuration, only : R
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellL, InvCL
use CommonBlocks, only : QSwitch
use NonbondParam, only : Charge, Rminh, EpsLJ

implicit none

integer, dimension(3) :: molnumber
integer :: i, j, ir, n1, n2, n3
integer :: ii, im1, im2, ia1, ia2
real(8), dimension(3) :: Rij, Rg
integer, dimension(3) :: Nij
real(8) :: Upot, fk, ek, fk1
real(8) :: Sgm, Sgm2, Eps, SR2, SR6, SR12, cf, InvR2
real(8) :: R1, R2

   molnumber(1) = n1
   molnumber(2) = n2
   molnumber(3) = n3

   do im2 = 1, NumW

     if(im2==n1.or.im2==n2.or.im2==n3) cycle
     Rij(:) = Rw(:,im2) - Rg(:)
     call MImageO(Rij,CellL,InvCL,R2,Nij)

     if(R2 > Rcutoff2) cycle

     R1 = sqrt(R2)
     ir = int(R1*10.d0) + 1
     gra(ir) = gra(ir) + 1.d0

     Upot = 0.d0

     do ii = 1, 3
       im1 = molnumber(ii)
       do ia1 = 1, 3
         i = (im1-1)*3 + ia1 + NumIni
         do ia2 = 1, 3
           j = (im2-1)*3 + ia2 + NumIni
           Rij(:) = R(:,i) - R(:,j)
           call MImageO(Rij,CellL,InvCL,R2,Nij)

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

     pota(ir) = pota(ir) + Upot

   end do


end subroutine cal_rdf_ene_a


!######################################################################
!######################################################################


subroutine search_possible_pair(Rg,i,j,k)

use CGanalyW

implicit none

integer :: i, j, k, ip, iq
integer :: itrial
logical :: Qgoodpair
real(8), dimension(3) :: Rg
real(8) :: ranf
external ranf

#ifdef DEBUG
   print *, 'In search_possible_pair'
   print *, 'i =',i,number_pair(i)
   print *, 'j =',j,number_pair(j)
   print *, 'k =',k,number_pair(k)
#endif

   Qgoodpair = .False.

   if(number_pair(i)<2) then
     i = j
     if(number_pair(j)<2) then
       i = k
       if(number_pair(k)<2) then
         write(*,*) "no pairs for ",i,j,k
         call Finalize
       end if
     end if
   end if

   itrial = 0

   do while(.not.Qgoodpair)

     ip = int(ranf()*number_pair(i)) + 1

     j = id_pair(i,ip)

     iq = ip
     do while (iq==ip)
       iq = int(ranf()*number_pair(i)) + 1
     end do

     k = id_pair(i,iq)

     itrial = itrial + 1
     if(itrial > 1000) then
       write(*,*) 'cannot find a reasonable pair'
       call Finalize
     end if

     call find_rg(Rg,i,j,k)
     call check_blob(Qgoodpair,Rg,i,j,k)

   end do

end subroutine search_possible_pair


!######################################################################
!######################################################################


subroutine cal_blob_blob(Ri,im,jm,km)

use CGanalyW
use Configuration, only : R
use CellParam, only : CellL, InvCL
use CutoffParam, only : Rcutoff2

implicit none

integer :: im, jm, km, l, ll, ip, iq, ir, n1, j1, k1, i1
integer :: i, j, k, itrial, reducedpair, jj
integer, dimension(50) :: id_reduced_pair
real(8), dimension(3) :: Ri, Rj
real(8), dimension(3,3,3) :: RwI
logical, dimension(NumW) :: Qdone
real(8), dimension(3) :: Rij
integer, dimension(3) :: Nij, Ilist
real(8) :: R2, subt, plus, ranf
logical :: Qgoodpair
external ranf

#ifdef DEBUG
   print *, 'In cal_blob_blob'
#endif

   Qdone(:) = .False.
   Qdone(im) = .True.
   Qdone(jm) = .True.
   Qdone(km) = .True.
   WCount(:) = 0
   RWCount(:,:) = 0

   Ilist(1) = im
   Ilist(2) = jm
   Ilist(3) = km

   do n1 = 1, 3 !mol
     k1 = Ilist(n1)
     do j1 = 1, 3 !atom
       i1 = (k1-1)*3+j1+NumIni
       RwI(:,j1,n1) = R(:,i1) - Ri(:)
       RwI(:,j1,n1) = RwI(:,j1,n1) - CellL(:)*nint(RwI(:,j1,n1)*InvCL(:))
!       R2 = dot_product(RwI(:,j1,n1),RwI(:,j1,n1))
!       if(R2>16.) write(*,*) 'WARNING:too large distance from com for blob'
     end do
   end do

   do l = 1, NumW

! ## pick up the first mol

     if(Qdone(l)) cycle
     i = l
     Rij(:) = Rw(:,i) - Ri(:)
     call MImageO(Rij,CellL,InvCL,R2,Nij)
     if(R2 > Rcutoff2+100.d0) then
       Qdone(l)=.True.
       cycle
     end if

#ifdef DEBUG
     print *, 'pick up l=',l
#endif

     reducedpair = 0
     id_reduced_pair(:) = 0

     do ll = 1, number_pair(l)
       jj = id_pair(l,ll)
       if(jj==im.or.jj==jm.or.jj==km) cycle
       reducedpair = reducedpair + 1
       id_reduced_pair(reducedpair) = jj
     end do

     if(reducedpair<2) then
       write(*,*) 'error : too few pairs for ',i
       call Finalize
     end if

     Qgoodpair = .False.
     itrial = 0

     do while(.not.Qgoodpair)

       ip = int(ranf()*reducedpair) + 1

       j = id_reduced_pair(ip)

       itrial = itrial + 1

       if(j==im.or.j==jm.or.j==km) then
         write(*,*) 'ERROR:WRONG CODE'
         call Finalize
       end if

       iq = ip
       do while (iq==ip)
         iq = int(ranf()*reducedpair) + 1
       end do

       k = id_reduced_pair(iq)

       if(itrial > 1000) then
         write(*,*) 'cannot find a reasonable pair for ',l
         exit
       end if

       call find_rg(Rj,i,j,k)
       call check_blob(Qgoodpair,Rj,i,j,k)

     end do

     if(.not.Qgoodpair) then
       print *, im, jm, km, l
       cycle
     end if

#ifdef DEBUG
     print *, 'pick up blob = ', i,j,k
#endif

     Qdone(i) = .True.
     Qdone(j) = .True.
     Qdone(k) = .True.

     call Ene_blob_blob(Ri, Rj, RwI, im, jm, km, i, j, k)

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
       grb(ir) = grb(ir) + plus
     end do
   end do

end subroutine cal_blob_blob




subroutine Ene_blob_blob(Ri, Rj, RwI, ii, ij, ik, ji, jj, jk)

use CGanalyW
#ifdef DEBUG
use Configuration, only : R
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellL, InvCL
use CommonBlocks, only : QSwitch
use NonbondParam, only : Charge, Rminh, EpsLJ
use UnitExparam, only : cvol
integer :: ir
#else
use Configuration, only : R
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellL, InvCL
use CommonBlocks, only : QSwitch
use NonbondParam, only : Charge, Rminh, EpsLJ
#endif
integer :: ii, ij, ik, ji, jj, jk, i, j, nnn
real(8), dimension(3) :: Ri, Rj, Rij, RGij
integer, dimension(3) :: Nij
real(8) :: R2, R1, Uintra1, Uintra2
integer, dimension(3) :: Ilist, Jlist
real(8) :: cf, Upot, Sgm, Sgm2, Eps, InvR2, SR2, SR6, SR12, ek, fk1
real(8), dimension(3,3,3) :: RwI, RwJ ! (dim,atm,mol)

#ifdef DEBUG
   print *, 'In Ene_blob_blob'
#endif

   Ilist(1) = ii
   Ilist(2) = ij
   Ilist(3) = ik
   Jlist(1) = ji
   Jlist(2) = jj
   Jlist(3) = jk

   do n1 = 1, 3 !mol
     k2 = Jlist(n1)
     do j1 = 1, 3 !atom
       j = (k2-1)*3+j1+NumIni
       RwJ(:,j1,n1) = R(:,j) - Rj(:)
       RwJ(:,j1,n1) = RwJ(:,j1,n1) - CellL(:)*nint(RwJ(:,j1,n1)*InvCL(:))
!       R2 = dot_product(RwJ(:,j1,n1),RwJ(:,j1,n1))
!       if(R2>20.) write(*,*) 'WARNING:too large distance from com for blob J', sqrt(R2)
     end do
   end do

   Rij = Ri - Rj
   call MImageO(Rij,CellL,InvCL,R2,Nij)

   if(R2 > Rcutoff2) Return

   RGij(:) = Rij(:)
   R1 = sqrt(R2)
   nnn = int(10.d0*R1) + 1
   gbb(nnn) = gbb(nnn) + 1.d0
   grb(nnn) = grb(nnn) + 1.d0

#ifdef DEBUG
   print *, 'Distance blob_blob = ', R1
   print *, 'ir=',nnn
#endif

   WCount(ji) = WCount(ji) + 1
   WCount(jj) = WCount(jj) + 1
   WCount(jk) = WCount(jk) + 1
   RWCount(ji,WCount(ji)) = nnn
   RWCount(jj,WCount(jj)) = nnn
   RWCount(jk,WCount(jk)) = nnn

#ifdef DEBUG
   print *, 'ir=',nnn
#endif

   Upot = 0.d0
   do n1 = 1, 3
     k1 = Ilist(n1)
     do n2 = 1, 3
       k2 = Jlist(n2)
       do j1 = 1, 3
       do j2 = 1, 3

       i = (k1-1)*3+j1+NumIni
       j = (k2-1)*3+j2+NumIni

       Rij(:) = RwI(:,j1,n1) - RwJ(:,j2,n2) + RGij(:)
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
   end do

#ifdef DEBUG
   print *, 'Upot (blob-blob) = ',Upot*cvol
   print *, 'ir = ', nnn
#endif

   Ene_Blob(nnn) = Ene_Blob(nnn) + Upot

   call cal_internal(Uintra1,ii,ij,ik)
   call cal_internal(Uintra2,ji,jj,jk)

   Ene_total_blob(nnn) = Ene_total_blob(nnn) + Upot + Uintra1 + Uintra2
   subUin = subUin + Uintra1 + Uintra2

#ifdef DEBUG
   print *, 'Add Ene_blob'
#endif

end subroutine Ene_blob_blob
