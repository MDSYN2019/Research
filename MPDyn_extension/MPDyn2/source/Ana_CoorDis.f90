! ############################
! ## SUBROUTINE LIST 
! ## -- CoorDis
! ############################

!######################################################################
!######################################################################


subroutine CoorDis

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use CellParam, only : H, InvH, CellShft
use CutoffParam, only : Rcutoff2
use AtomParam, only : MolName, AtomName, Mass
use ParamAnalyze
use CommonMPI, only : NProcs, MyRank

implicit none

integer :: i, j, k, l, ii, kk, ll, ka, kb, lc
integer, dimension(iatom) :: iia
integer, dimension(jatom) :: jja
integer, dimension(iatom,imol) :: lista
integer, dimension(jatom,jmol) :: listb
real(8), dimension(3,iatom,imol) :: Rsi
real(8), dimension(3,jatom,jmol) :: Rsj
real(8), dimension(3,imol) :: Rcomi,ScRi
real(8), dimension(3,jmol) :: Rcomj,ScRj
real(8) :: iDRgrid, iDNgrid, weig1, weig2
real(8) :: Six, Siy, Siz, Rix, Riy, Riz
real(8) :: Sx, Sy, Sz, Rx, Ry, Rz
real(8) :: R2, co, dix, diy, diz, dx, dy, dz, d2, d1
real(8) :: bs, bb, ci, R1, xx
integer :: Nx, Ny, Nz, Nas
integer :: iatom_all, jatom_all, ini1, ini2
real(8) :: wmass_1, wmass_2, inv_mass1, inv_mass2
real(8), dimension(ngrid) :: HistC
integer, dimension(ngrid) :: Count
integer, dimension(ngrid,ncgrid) :: CrossCount

   open(51,file='./Analy/CoorDisAv.dat',status='unknown')
   open(52,file='./Analy/CoorDisP.dat',status='unknown')

   ii = 0
   do i = 1, NumSpec
     if(MolName(i)==MOLNAME1) then
       ini1 = ii
       iatom_all = NumAtm(i)
       do kk = 1, iatom
lpa:     do k = 1, iatom_all
           if(AtomName(ii+k)==NAtomI(kk)) then
             iia(kk) = k
             exit lpa
           end if
           if(k==NumAtm(i)) then
             write(*,*) 'ERROR: no atom found for ',NAtomI(kk)
             call Finalize
           end if
         end do lpa
       end do
       do j = 1, imol
         do kk = 1, iatom
           lista(kk,j) = ii + iia(kk)
         end do
         ii = ii + NumAtm(i)
       end do
     else if(MolName(i)==MOLNAME2) then
       ini2 = ii
       jatom_all = NumAtm(i)
       do kk = 1, jatom
lpb:     do k = 1, jatom_all
           if(AtomName(ii+k)==NAtomJ(kk)) then
             jja(kk) = k
             exit lpb
           end if
           if(k==NumAtm(i)) then
             write(*,*) 'ERROR: no atom found for ',NAtomJ(kk)
             call Finalize
           end if
         end do lpb
       end do
       do j = 1, jmol
         do kk = 1, jatom
           listb(kk,j) = ii + jja(kk)
         end do
         ii = ii + NumAtm(i)
       end do
     else
       ii = ii + NumMol(i)*NumAtm(i)
     end if
   end do

   weig1 = 0.d0
   do k = 1, iatom_all
     weig1 = weig1 + Mass(ini1+k)
   end do
   weig1 = 1.d0 / weig1

   weig2 = 0.d0
   do k = 1, jatom_all
     weig2 = weig2 + Mass(ini2+k)
   end do
   weig2 = 1.d0 / weig2

   dAB = 1.d0/(dAB*dAB)
   iDRgrid = 1.d0/DRgrid
   iDNgrid = 1.d0/DNgrid

   HistC(:) = 0.d0
   Count(:) = 0
   CrossCount(:,:) = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif

!     ---------------
       call BcastRH
!     ---------------

       call InversMatrix(H,InvH)
       call TransCellList

       ii = ini1
       do k = 1, imol
         Rcomi(:,k) = 0.d0
         do kk = 1, iatom_all
           ii = ii + 1
           Rcomi(:,k) = Rcomi(:,k) + Mass(ii)*R(:,ii)
         end do
         do kk = 1, iatom
           l = lista(kk,k)
           Rsi(:,kk,k) = R(:,l)
         end do
         Rcomi(:,k) = Rcomi(:,k) * weig1
       end do
       ii = ini2
       do k = 1, jmol
         Rcomj(:,k) = 0.d0
         do kk = 1, jatom_all
           ii = ii + 1
           Rcomj(:,k) = Rcomj(:,k) + Mass(ii)*R(:,ii)
         end do
         do kk = 1, jatom
           l = listb(kk,k)
           Rsj(:,kk,k) = R(:,l)
         end do
         Rcomj(:,k) = Rcomj(:,k) * weig2
       end do

       do k = 1 , imol
         ScRi(:,k) = matmul( InvH , Rcomi(:,k) )
       end do
       do k = 1 , jmol
         ScRj(:,k) = matmul( InvH , Rcomj(:,k) )
       end do

       Nas = NProcs - MyRank

       do k = Nas, imol, NProcs
         Six = ScRi(1,k)
         Siy = ScRi(2,k)
         Siz = ScRi(3,k)
         Rix = Rcomi(1,k)
         Riy = Rcomi(2,k)
         Riz = Rcomi(3,k)
         do kk = 1, jmol
           Sx = Six - ScRj(1,kk)
           Sy = Siy - ScRj(2,kk)
           Sz = Siz - ScRj(3,kk)
           Rx = Rix - Rcomj(1,kk)
           Ry = Riy - Rcomj(2,kk)
           Rz = Riz - Rcomj(3,kk)
           if(Sx>0.5) then
             Nx = -9
           else if(Sx<-0.5) then
             Nx =  9
           else 
             Nx =  0
           end if
           if(Sy>0.5) then
             Ny = -3
           else if(Sy<-0.5) then
             Ny =  3
           else
             Ny =  0
           end if
           if(Sz>0.5) then
             Nz = -1
           else if(Sz<-0.5) then
             Nz =  1
           else
             Nz =  0
           end if
           l = Nx + Ny + Nz
           Rx = Rx + CellShft(1,l)
           Ry = Ry + CellShft(2,l)
           Rz = Rz + CellShft(3,l)
           R2 = Rx*Rx + Ry*Ry + Rz*Rz

           if ( R2 <= Rcutoff2 ) then

             co = 0.d0
             do ka = 1, iatom
               dix = Rsi(1,ka,k)
               diy = Rsi(2,ka,k)
               diz = Rsi(3,ka,k)
               do kb = 1, jatom
                 dx = dix - Rsj(1,kb,kk) + CellShft(1,l)
                 dy = diy - Rsj(2,kb,kk) + CellShft(2,l)
                 dz = diz - Rsj(3,kb,kk) + CellShft(3,l)
                 d2 = dx*dx + dy*dy + dz*dz
                 d1 = d2*dAB
                 bs = d1
                 do ii = 1, mex-1
                   bs = bs * d1
                 end do
                 bb = bs
                 do ii = mex, nex-1
                   bb = bb * d1
                 end do
                 ci = (1.d0-bs)/(1.d0-bb)
                 co = co + ci
               end do
             end do
!             co = co / iatom

             R1 = sqrt(R2)
             ll = int(R1*iDRgrid) + 1
             lc = int(co*iDNgrid) + 1
             Count(ll) = Count(ll) + 1
             HistC(ll) = HistC(ll) + co
             CrossCount(ll,lc) = CrossCount(ll,lc) + 1

           end if

         end do

       end do

     end do

   end do

   do i = 1, ngrid
     if(Count(i)/=0) then
       HistC(i) = HistC(i) / Count(i)
       R1 = i*DRgrid
       write(51,'(f6.2,f10.5)') R1, HistC(i)
     end if
   end do

   do i = 1, ngrid
     R1 = (i-0.5)*DRgrid
     do j = 1, ncgrid
       co = (j-0.5)*DNgrid
       xx = dble(CrossCount(i,j))/dble(Count(i))*100.d0
       write(52,'(2f6.2,f10.5)') R1,co,xx
     end do
     write(52,*)
   end do

   close(51)
   close(52)

end subroutine CoorDis

