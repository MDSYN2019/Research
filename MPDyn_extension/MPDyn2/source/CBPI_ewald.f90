! ############################
! ## SUBROUTINE LIST 
! ## -- Cav_ParticleInsertion 
! ## -- Pre_CBPI 
! ## -- List_Insertion 
! ## -- Ene_Insertion 
! ############################

module CBPI_kspace
real(8), dimension(10000) ::CsSum_ol, SnSum_ol, Epkn_ol
real(8), dimension(3,10000) :: Wk_ol
real(8) :: Ene_k_ol, EneT_intra_ol
real(8), dimension(:), allocatable :: EneT_intra_PI, Ene_Eslf_PI, Ene_Eksp_PI
end module CBPI_kspace

!######################################################################
!######################################################################


! **************************************************
! **  Cavity based particle insertion technique   **
! **  Procedure                                   **
! **     1. map the system on the uniform grid    **
! **     2. find cavity grids                     **
! **     3. trial insertion                       **
! **************************************************

subroutine Cav_ParticleInsertion

use Numbers, only : N, NumSpec, NumMol
use AtomParam, only: MolName
use CommonBlocks, only : QMaster
use Configuration, only : R
use ParamAnalyze
use CommonMPI
use RBparam
use BathParam, only : Temp_o, Pressure_o, Beta
use CellParam, only : H, Volume
use CBPI_kspace

implicit none

character(len=72) :: String
integer :: i, j, Nas, k, kk, l, jj, NumF, Ncavz
integer :: lx, ly, lz, ix, iy, iz, Icav, TotalStepNumber
integer :: Nx, Ny, Nz, lzh
real(8) :: Rxmaxh, Rymaxh, Rzmaxh
integer, dimension(:,:,:), allocatable :: Flcav
integer, dimension(:), allocatable :: NumCavZ, NumCavZAcc
integer, dimension(:,:,:), allocatable :: NumCavXY
real(8), dimension(:,:), allocatable :: EneZ
real(8), dimension(:,:), allocatable :: EneZPT
real(8), dimension(:), allocatable :: ProbCavZ
! >> Ewald
real(8), dimension(:,:), allocatable :: AveExpEN, FreeEne
! <<
real(8), dimension(4) :: qt
real(8), dimension(10) :: qt2
real(8), dimension(3) :: Rg, Ri, Rext
real(8), dimension(3,3) :: Rot
real(8), dimension(3,N) :: dR
real(8) :: ranf, q2
real(8), parameter :: kbmol = 1.98720d-03 ! gas const. [kcal/(K mol)]
real(8) :: kTmol, det, PTm
real(8) :: www, ww
real(8) :: AveEneZ, AveEneZPT
real(8) :: CPot, CPotPT
external ranf, det

! ## kT [kcal / mol]
   kTmol = kbmol * Temp_o

   Nas = NProcs - MyRank

! ## Read the parameters of the molecules to be inserted.
   call Pre_CBPI

!   print *, 'Rzmax=',Rzmax
!   print *, 'R_cav=',R_cav

   lx    = int( Rxmax / R_cav )
   ly    = int( Rymax / R_cav )
   lz    = int( Rzmax / R_cav )

!   print *, 'lz=',lz

   Rxmax = dble(lx) * R_cav
   Rymax = dble(ly) * R_cav
   Rzmax = dble(lz) * R_cav

!   print *, 'Rzmax=',Rzmax

   Rxmaxh = Rxmax * 0.5d0
   Rymaxh = Rymax * 0.5d0
   Rzmaxh = Rzmax * 0.5d0

   TotalStepNumber = 0

   do i = 1 , NJobs

     TotalStepNumber = TotalStepNumber + NTrjStep(i)

   end do

   allocate( Flcav(lx,ly,lz) )
   allocate( NumCavZ(lz) )
   allocate( NumCavZAcc(lz) )
   allocate( EneZ(lz,NumSpecInsert) )
   allocate( EneZPT(lz,NumSpecInsert) )
   allocate( NumCavXY(2,300,lz) )

   NumCavZAcc = 0
   EneZ   = 0.d0
   EneZPT = 0.d0

   NumF = 0

   ww = 1.d0 / dble(NumInsertion * NumOriTry)


! ## Start the main routine

   do i = 1 , NJobs

     if(QMaster) call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ------------------------------
#ifdef MOLFILE
       if(QMaster) call Read_RTraj(i)
#else
       if(QMaster) call Read_RTraj
#endif
!     ------------------------------

       NumF = NumF + 1

       if(QMaster.and.mod(NumF,10)==0) print *, 'Step = ',NumF

       call BcastRH

! ## Cell rotation !
       call Transform

! >> Ewald
       if(i == 1 .and. j == 1) then
         call ErrorFuncList           !  make a table of error-function
         call RecLatticeList_CBPI          !  define reciprocal lattice

         call Ewald_SelfTerm_CBPI
       endif
!<<

! ## correction for the NPT ensemble
       Volume = det(H)

       PTm = Pressure_o * Volume * Beta

! >> Ewald
       if(mod(NumF,10)==1) then
         call Ene_kspace_preCBPI
       endif
!<<

! ## Search the cavity,  Flcav = 0 : cavity, 1 : occupied 

       Flcav = 0

       do jj = 1 , N

         Ri = R(:,jj)

         do ix = -1, 1
           do iy = -1, 1
             do iz = -1, 1

               Rext(1) = Ri(1) + H(1,1) * ix + H(1,2) * iy + H(1,3) * iz + Rxmaxh
               Rext(2) = Ri(2) + H(2,1) * ix + H(2,2) * iy + H(2,3) * iz + Rymaxh
               Rext(3) = Ri(3) + H(3,1) * ix + H(3,2) * iy + H(3,3) * iz + Rzmaxh

               if((Rext(1) < Rxmax).and.(Rext(1) >= 0.d0).and.&
               &  (Rext(2) < Rymax).and.(Rext(2) >= 0.d0).and.&
               &  (Rext(3) < Rzmax).and.(Rext(3) >= 0.d0)) then

                 Nx = int(Rext(1) / R_cav) + 1
                 Ny = int(Rext(2) / R_cav) + 1
                 Nz = int(Rext(3) / R_cav) + 1

                 Flcav(Nx,Ny,Nz)  = 1

               end if

             end do
           end do
         end do

       end do

! ## Number of cavity as a function of z 

       NumCavXY = 0

       do iz = 1 , lz

         NumCavZ(iz) = 0

         do ix = 1, lx
           do iy = 1, ly

             if(Flcav(ix,iy,iz) == 0) then

               NumCavZ(iz) = NumCavZ(iz) + 1

               if(NumCavZ(iz) > 300) then

                 write(*,*) 'ERROR : too large NumCavZ in (CBPI)'
                 call Finalize

               end if

               NumCavXY(1,NumCavZ(iz),iz) = ix
               NumCavXY(2,NumCavZ(iz),iz) = iy

             end if

           end do
         end do

         NumCavZAcc(iz) = NumCavZAcc(iz) + NumCavZ(iz)

       end do

! ## Cavity Insertion

       do iz = 1 , lz

         Ncavz = NumCavZ(iz)

         EneInsert = 0.d0
         EneInsertPT = 0.d0

         do k = Nas , NumInsertion, NProcs

           Icav = int( Ncavz * ranf() ) + 1
           ix = NumCavXY(1,Icav,iz)
           iy = NumCavXY(2,Icav,iz)

           Rg(1) = ( dble(ix) - 0.5d0 ) * R_cav - Rxmaxh
           Rg(2) = ( dble(iy) - 0.5d0 ) * R_cav - Rymaxh
           Rg(3) = ( dble(iz) - 0.5d0 ) * R_cav - Rzmaxh

! ## make a pair-list
           call List_Insertion(Rg,dR)

! ## NumOriTry is the number of insertion trials for one COM with different orientation of 
! ## the molecule
           do kk = 1 , NumOriTry

             qt(1) = 2.d0 * ranf() - 1.d0
             qt(2) = 2.d0 * ranf() - 1.d0
             qt(3) = 2.d0 * ranf() - 1.d0
             qt(4) = 2.d0 * ranf() - 1.d0

             q2 = dot_product( qt, qt )
             qt = qt / sqrt(q2)

             qt2( 1) = qt(1) * qt(1)
             qt2( 2) = qt(2) * qt(2)
             qt2( 3) = qt(3) * qt(3)
             qt2( 4) = qt(4) * qt(4)

             qt2( 5) = qt(1) * qt(2)
             qt2( 6) = qt(1) * qt(3)
             qt2( 7) = qt(1) * qt(4)
             qt2( 8) = qt(2) * qt(3)
             qt2( 9) = qt(2) * qt(4)
             qt2(10) = qt(3) * qt(4)

             Rot(1,1) = qt2(1) + qt2(2) - qt2(3) - qt2(4)
             Rot(2,2) = qt2(1) - qt2(2) + qt2(3) - qt2(4)
             Rot(3,3) = qt2(1) - qt2(2) - qt2(3) + qt2(4)

             Rot(1,2) = 2.d0 * ( qt2( 8) + qt2( 7) )
             Rot(2,1) = 2.d0 * ( qt2( 8) - qt2( 7) )
             Rot(1,3) = 2.d0 * ( qt2( 9) - qt2( 6) )
             Rot(3,1) = 2.d0 * ( qt2( 9) + qt2( 6) )
             Rot(2,3) = 2.d0 * ( qt2(10) + qt2( 5) )
             Rot(3,2) = 2.d0 * ( qt2(10) - qt2( 5) )

             do l = 1 , NumSpecInsert

               do jj = 1 , NumAtmInsert(l)

                 Rmolec(1,jj,l) = Rot(1,1) * R_onMol(1,jj,l) &
                 &              + Rot(2,1) * R_onMol(2,jj,l) & 
                 &              + Rot(3,1) * R_onMol(3,jj,l)
                 Rmolec(2,jj,l) = Rot(1,2) * R_onMol(1,jj,l) &
                 &              + Rot(2,2) * R_onMol(2,jj,l) & 
                 &              + Rot(3,2) * R_onMol(3,jj,l)
                 Rmolec(3,jj,l) = Rot(1,3) * R_onMol(1,jj,l) &
                 &              + Rot(2,3) * R_onMol(2,jj,l) &
                 &              + Rot(3,3) * R_onMol(3,jj,l)

               end do

             end do

! ##  calculation of exp(- beta * E_N+1)
! >> EWALD
             call Ene_kspace_CBPI(0,Rg)
! <<
             call Ene_Insertion(dR)

! ##  summation of exp(- beta * E_N+1)

             EneInsert   = EneInsert   + EneIns
             EneInsertPT = EneInsertPT + EneIns * PTm

           end do

         end do

         do k = 1, NumSpec
           if(MolName(k) == 'TIP3') then
            EneInsertPT(1) = EneInsertPT(1) /(NumMol(k)+1.0) 
           endif
         end do

         call SumEneInsertion

         if(QMaster) then

           do l = 1, NumSpecInsert

             EneZ(iz,l)   = EneZ(iz,l)   + EneInsert(l)   * ww
             EneZPT(iz,l) = EneZPT(iz,l) + EneInsertPT(l) * ww

           end do

         end if

       end do

     end do

   end do

! ## Averaging

   if(QMaster) then

   allocate( ProbCavZ(lz) )

   lzh = lz / 2

   open(51,file='./Analy/CavityMU',status='unknown')

   do iz = 1 , lz

     ProbCavZ(iz) = dble( NumCavZAcc(iz) ) / dble(lx * ly * TotalStepNumber)

     if(lzh*2==lz) then

       write(51,'(f7.2,e14.6)') (iz-lzh-0.5)*R_cav, -kTmol * log( ProbCavZ(iz) )

     else

       write(51,'(f7.2,e14.6)') (iz-lzh-1.0)*R_cav, -kTmol * log( ProbCavZ(iz) )

     end if

   end do

   close(51)

!   print *, 'NumOriTry=',NumOriTry
!   print *, 'NumInsertion=',NumInsertion
!   print *, 'TotalStepNumber=',TotalStepNumber
!   print *, 'lx,ly=',lx,ly
!   print *, 'ProbCavZ=',ProbCavz

!   print *, 'kT=',kT

   www = 1.d0 / dble(TotalStepNumber)

   do i = 1 , NumSpecInsert

     write(String,'(a,i1,a)') './Analy/CavityInsert_',i,'.data'
     open(51,file=trim(String),status='unknown')

     do iz = 1 , lz

       AveEneZ   = EneZ(iz,i)   * www
       AveEneZPT = EneZPT(iz,i) * www

       CPot   = - kTmol * ( log( AveEneZ   ) + log(ProbCavZ(iz)) )
       CPotPT = - kTmol * ( log( AveEneZPT ) + log(ProbCavZ(iz)) )

       if(lzh*2==lz) then

         write(51,'(f7.2,2e14.6)') (iz-lzh-0.5)*R_cav, CPot, CPotPT

       else

         write(51,'(f7.2,2e14.6)') (iz-lzh-1.0)*R_cav, CPot, CPotPT

       end if

     end do

     close(51)

   end do

   end if

end subroutine Cav_ParticleInsertion


!######################################################################
!######################################################################


! *********************************
! ** set up for cavity insertion **
! *********************************

subroutine Pre_CBPI

use CommonBlocks, only : QMaster
use ParamAnalyze
use RBparam
use BookParam, only : MaxPair
use UnitExParam, only : ExParam, ec
use CutoffParam, only : Rcutoff2
use CBPI_kspace

implicit none

integer :: i, j
integer , parameter :: MaxA = 5

   allocate( ListI(MaxPair) )

   open(51,file=trim(adjustl(CBPI_FILE)))

   read(51,*) NumSpecInsert

   allocate( MolNameCBPI(NumSpecInsert) )
   allocate( NumAtmInsert(NumSpecInsert) )
   allocate( Rmolec(3,MaxA,NumSpecInsert) )
   allocate( R_onMol(3,MaxA,NumSpecInsert) )

   allocate( EpsLJPI(MaxA,NumSpecInsert) )
   allocate( RminhPI(MaxA,NumSpecInsert) )
   allocate( chargePI(MaxA,NumSpecInsert) )

   allocate( EneIns(NumSpecInsert) )
   allocate( EneInsert(NumSpecInsert) )
   allocate( EneInsertPT(NumSpecInsert) )

   allocate( Ene_Eslf_PI(NumSpecInsert) )
   allocate( EneT_intra_PI(NumSpecInsert) )
   allocate( Ene_Eksp_PI(NumSpecInsert) )

   do i = 1 , NumSpecInsert

     read(51,*) MolNameCBPI(i)

     read(51,*) NumAtmInsert(i)

     if(NumAtmInsert(i)>MaxA) then
       write(6,*) 'ERROR : too large NumAtmInsert in (Pre_CBPI)'
       call Finalize
     end if

     do j = 1 , NumAtmInsert(i)

       read(51,*) EpsLJPI(j,i), RminhPI(j,i), chargePI(j,i)

     end do

     do j = 1 , NumAtmInsert(i)

       read(51,*) R_onMol(:,j,i)

     end do

   end do

   close(51)

   if(QMaster) then

   write(6,'(a/)')      '>>>  Cavity Insertion has just been started ! '
   write(6,'(a,f8.2/)') '>>>  Cut off length = ', sqrt(Rcutoff2)
   write(6,'(a,i2,a)')  '>>>  ',NumSpecInsert,' molecules to be inserted !! '
   write(6,'(a/)')      '>>>  parameters of the inserted molecules'

   do i = 1 , NumSpecInsert

     write(6,'(a)')    '##########################################################'
     write(6,'(a,i3)') '###  Molecule ',i
     write(6,'(a,i3)') '###  Number of Atom per Molecule',NumAtmInsert(i)

     do j = 1 , NumAtmInsert(i)

       write(6,'(a,i3)')      '---  Atom number : ',j
       write(6,'(a,f10.5,a)') '         Epsilon = ',EpsLJPI(j,i),'kcal/mol'
       write(6,'(a,f10.5,a)') '         Rminh   = ',RminhPI(j,i),'A'
       write(6,'(a,f10.5,a/)') '         charge  = ',chargePI(j,i),'e'

     end do

   end do

   end if

   do i = 1 , NumSpecInsert

     do j = 1 , NumAtmInsert(i)

       EpsLJPI(j,i) = EpsLJPI(j,i) * ExParam

       if( EpsLJPI(j,i) /= 0. ) then

         EpsLJPI(j,i) = sqrt( abs( EpsLJPI(j,i) ) )

       end if

       ChargePI(j,i) = ChargePI(j,i) * sqrt(ec)

     end do

   end do

end subroutine Pre_CBPI


!######################################################################
!######################################################################


! ***********************************
! ** list up the neighboring atoms **
! ***********************************

subroutine List_Insertion(Rg,dR)

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use BookParam, only : Npair
use CellParam, only : H, InvH
use CutoffParam, only : Rcutoff2

implicit none

integer :: i
real(8), dimension(3) :: Rg, Sg, Rij, Sij
real(8), dimension(3,N) :: ScR, dR
real(8) :: R2

   Sg = - matmul(InvH, Rg)

   do i = 1 , N

     ScR(:,i) = matmul( InvH, R(:,i) )

   end do

   Npair = 0

   do i = 1, N

     Sij = ScR(:,i) + Sg
     Sij = Sij - nint( Sij )
     Rij = matmul( H, Sij )

     R2  = dot_product( Rij, Rij )

     if( R2 <= Rcutoff2 ) then

       Npair = Npair + 1
       ListI(Npair) = i
       dR(:,Npair)  = Rij

     end if

   end do

end subroutine List_Insertion


!######################################################################
!######################################################################


! ************************
! ** energy calculation **
! ************************

subroutine Ene_Insertion(dR)

use Numbers, only : N
use CutoffParam, only : Rcutoff2, Ron2
use ParamAnalyze
use RBparam
use BookParam, only : Npair
use NonbondParam, only : Charge, Rminh, EpsLJ
use BathParam, only : Beta
use EwaldParam, only : Alpha
use CBPI_kspace

implicit none

integer :: i, j, k, l
real(8), dimension(3) :: Rij
real(8), dimension(3,N) :: dR
real(8) :: R2, Sgm, Sgm2, Eps, InvR2
real(8) :: SR2, SR6, SR12, ek, cf, fk1
real(8) :: EneT
real(8) :: Ben

!ewald 
! r-space switching func
real(8) :: swf1, sRoff2, sRon2
real(8) :: Xoff,Xon
real(8) :: Switch
real(8) :: ErrorFunc
real(8) :: x, R1
!temporary
real(8) :: Ene_check, ek_cutoff, E_cloumb_cutoff, Ben_check, dumcloumb,dum_Ersp
real(8), parameter :: kbmol = 1.98720d-03 ! gas const. [kcal/(K mol)]
real(8) :: Error_Function
external Error_Function

   swf1   = 1.d0 / (Rcutoff2 - Ron2)**3

   sRoff2 = -Rcutoff2
   sRon2  = -Ron2

   do i = 1, NumSpecInsert

     EneT = 0.d0

     do j = 1, NumAtmInsert(i)

       do k = 1, Npair

         l   = ListI(k)
         Rij = dR(:,k) - Rmolec(:,j,i)
         R2  = dot_product( Rij, Rij )

         Sgm   = Rminh(l) + RminhPI(j,i)
         Sgm2  = Sgm * Sgm

         Eps   = EpsLJ(l) * EpsLJPI(j,i)
         InvR2 = 1.d0 / R2

         SR2  = Sgm2 * InvR2                    !(sigma/r)^2
         SR6  = SR2 * SR2 * SR2                 !         ^6
         SR12 = SR6 * SR6                       !         ^12
         ek   = Eps * ( SR12 - 2.d0 * SR6 )

! ----------------------------------------
!           Switching Function
! ----------------------------------------
         if ( R2 > Ron2 ) then        ! Switching range!!
           Xoff = R2 + sRoff2         ! x-xoff
           Xon  = R2 + sRon2          ! x-xon
           Switch  = ( 3.d0*Xoff*Xon - Xoff*Xoff*Xoff ) * swf1
           ek = ek * Switch
         end if
! ----------------------------------------

         cf = Charge(l) * ChargePI(j,i)

! >> Ewald
         if(cf /= 0.) then
           R1  = sqrt( R2 )
           x   = Alpha * R1
           ErrorFunc = Error_Function(x)
           fk1 = cf * ErrorFunc * sqrt(InvR2)
         else
           fk1 = 0.d0
         end if
! <<
         EneT = EneT + ek + fk1

       end do

     end do

! >> Ewald
     EneT = EneT + Ene_Eksp_PI(i) + Ene_Eslf_PI(i) + EneT_intra_PI(i)
! <<
     Ben = - Beta * EneT

     if(Ben > 700.d0) then

       print *, 'too large value of -beta*ene_N+1'
       print *, Ben
       Ben = 700.d0
       EneIns(i) = exp( Ben )

     else if(Ben < -708.d0) then

       EneIns(i) = 0.d0

     else

       EneIns(i) = exp( Ben )

     end if

   end do

end subroutine Ene_Insertion


!>> EWALD

!#####################################################################
!#####################################################################

! ***************************************************
! ** subroutine for making reciprocal lattice list **
! ** NOTE: Cell Shape is not taken into account !! **
! ***************************************************

subroutine RecLatticeList_CBPI

use CommonBlocks, only : QMaster
use CommonMPI
use EwaldParam, only : ih2mx, Nh, TmpNh, ih, kmaxx, kmaxy, kmaxz
use CellParam, only : H, InvH

implicit NONE

integer :: NAddition, Count
integer :: i, ihb, k, ix, iy, iz
real(8) :: a2, b2, c2, l2, hr2mx, hr2
real(8), dimension(3) :: a, b, c, hr
real(8), dimension(3,3) :: InvHt
integer, dimension(3,10000) :: ih_temp

   a = H(:,1)
   b = H(:,2)
   c = H(:,3)
   a2 = dot_product( a , a )
   b2 = dot_product( b , b )
   c2 = dot_product( c , c )
   l2 = (a2*b2*c2)**(1./3.d0)

   ihb    = int( sqrt( dble( ih2mx ) ) ) / 2 * 3
   hr2mx  = dble(ih2mx) / l2
   InvHt = Transpose(InvH)

   k = 0

   do ix =    0 , ihb
     do iy = -ihb , ihb
       do iz = -ihb , ihb

         hr(1) = InvHt(1,1)*ix + InvHt(1,2)*iy + InvHt(1,3)*iz
         hr(2) = InvHt(2,1)*ix + InvHt(2,2)*iy + InvHt(2,3)*iz
         hr(3) = InvHt(3,1)*ix + InvHt(3,2)*iy + InvHt(3,3)*iz

         hr2 = dot_product( hr , hr )

         if( hr2 > hr2mx ) cycle
         if( ( ix == 0 ) .and. ( iy <  0 ) ) cycle
         if( ( ix == 0 ) .and. ( iy == 0 ) .and. ( iz <= 0 ) ) cycle

         k = k + 1

         ih_temp(1,k) = ix
         ih_temp(2,k) = iy
         ih_temp(3,k) = iz

       end do

     end do
   end do

   TmpNh = k

   if(TmpNh > 10000) then
     if(QMaster) then
       write(*,*)  'error : Nh exceeds dimension of ih_temp(3,*)'
       write(11,*) 'error : Nh exceeds dimension of ih_temp(3,*)'
     end if

     stop

   end if

   do ix = 1, ihb
     hr2 = (ix * InvHt(1,1))**2
     if( hr2 > hr2mx ) cycle
     kmaxx = ix
   end do

   do iy = 1, ihb
     hr2 = (iy * InvHt(2,2))**2
     if( hr2 > hr2mx ) cycle
     kmaxy = iy
   end do

   do iz = 1, ihb
     hr2 = (iz * InvHt(3,3))**2
     if( hr2 > hr2mx ) cycle
     kmaxz = iz
   end do

   if(QMaster) then
     write(*,'(a,3i5)') 'kmax = ', kmaxx, kmaxy, kmaxz
   end if

   allocate( ih(3,Nh) )  !for CPOL

   do i = 1, Nh
     ih(:,i) = ih_temp(:,i)
   end do

end subroutine RecLatticeList_CBPI


!#####################################################################
!#####################################################################


! ***********************************************
! **  Ewald Self Energy term for partial space **
! ***********************************************

subroutine Ewald_SelfTerm_CBPI

use CommonBlocks
use ParamAnalyze, only : ChargePI, NumSpecInsert, NumAtmInsert
use CommonMPI
use CBPI_kspace
use UnitExParam, only : pi
use EwaldParam, only : Alpha
use RBparam, only : R_onMol

implicit NONE

integer :: i, j, k,m
real(8) :: InvRootPi
real(8), dimension(3) :: Rij
real(8) :: R2
real(8) ::  cf, fk1, R1, x
real(8) :: ErrorFunc
real(8) :: Error_Function
external Error_Function

   InvRootPi = 1.d0 / sqrt( pi )
   Ene_Eslf_PI = 0.d0
   EneT_intra_PI = 0.d0

!----------------------------------------------------------------------

   do i = 1, NumSpecInsert
     do j = 1, NumAtmInsert(i)
       Ene_Eslf_PI(i) = Ene_Eslf_PI(i) + ChargePI(j,i) * ChargePI(j,i)
     end do
   end do

   Ene_Eslf_PI = - Ene_Eslf_PI * Alpha * InvRootPi

   do i = 1, NumSpecInsert
     do m = 1, NumAtmInsert(i)-1
       do j = m +1 , NumAtmInsert(i)

         Rij = R_onMol(:,m,i) - R_onMol(:,j,i)

         R2  = dot_product( Rij, Rij )
         cf  = ChargePI(m,i) * ChargePI(j,i)

         R1  = sqrt( R2 )
         x   = Alpha * R1
         ErrorFunc = Error_Function(x)
         fk1 = cf * (ErrorFunc-1.d0) / R1
         EneT_intra_PI(i) = EneT_intra_PI(i) + fk1 

       end do
     end do
   end do

end subroutine Ewald_SelfTerm_CBPI


!#####################################################################
!#####################################################################


subroutine Ene_kspace_preCBPI

use NonbondParam, only : Charge
use Numbers, only : N
use Configuration, only : R
use CBPI_kspace
use CellParam, only : InvH
use UnitExParam, only : sqpi, pi2, InvPi
use EwaldParam, only : alp2, Nh, ih

implicit NONE
integer :: i, j, k, l, m
! k-space
real(8), dimension(3,3) :: InvHt
real(8) :: tyo, pref
real(8) :: zz
real(8), dimension(3) :: kn
real(8) :: kn2,  Trp

   InvHt = Transpose( InvH )
   pref = -sqpi * alp2

!--------------------------------------------------------------------------

   CsSum_ol = 0.d0
   SnSum_ol = 0.d0
   Ene_k_ol = 0.d0

   do l = 1 , Nh

     kn(:) = matmul( InvHt, ih(:,l) )
     kn2   = dot_product( kn(:) , kn(:) )
     Epkn_ol(l)  = exp( pref * kn2 ) / kn2

     Wk_ol(:,l) = kn(:) * pi2

     do m = 1, N
       zz = Charge(m)
       tyo = dot_product( Wk_ol(:,l) , R(:,m) )

       CsSum_ol(l) = CsSum_ol(l) + zz * cos(tyo)
       SnSum_ol(l) = SnSum_ol(l) + zz * sin(tyo)
     end do                             

     Trp = CsSum_ol(l) * CsSum_ol(l) + SnSum_ol(l) * SnSum_ol(l)

     Ene_k_ol = Ene_k_ol + Epkn_ol(l) * Trp

   end do

end subroutine Ene_kspace_preCBPI

!#####################################################################
!#####################################################################

subroutine Ene_kspace_CBPI(k,Rg)

use CommonBlocks
use CommonMPI
use ParamAnalyze
use CBPI_kspace
use EwaldParam, only : Nh
use RBparam, only : Rmolec
use CellParam, only : Volume
use UnitExParam, only : InvPi

implicit NONE
integer :: i, j, k, l, m
integer, parameter :: Ncite = 10
! k-space
real(8), dimension(3) :: Rg
real(8) :: tyo, pref, CslwSum, SnlwSum
real(8) :: zz
real(8), dimension(Ncite) :: Cslw, Snlw
real(8) ::  Trp

!ewald
real(8), dimension(:,:), allocatable :: Rdum
real(8) :: tempCsSum, tempSnSum


   Ene_Eksp_PI = 0.d0

!--------------------------------------------------------------------------

   allocate(Rdum(3,NumAtmInsert(1)))       

   do i = 1, NumSpecInsert

     do l = 1 , Nh

       tempCsSum = CsSum_ol(l)
       tempSnSum = SnSum_ol(l)

       do j = 1, NumAtmInsert(i)

         Rdum(:,j)= Rg(:) + Rmolec(:,j,i)
         zz = ChargePI(j,i)

         if(zz /= 0.) then
           tyo = dot_product( Wk_ol(:,l) , Rdum(:,j) )

           Cslw(j) = zz * cos(tyo)
           Snlw(j) = zz * sin(tyo)
           tempCsSum = tempCsSum + Cslw(j)
           tempSnSum = tempSnSum + Snlw(j)
         end if
       end do
       CslwSum = tempCsSum * tempCsSum 
       SnlwSum = tempSnSum * tempSnSum 

       Trp = CslwSum + SnlwSum

       Ene_Eksp_PI(i) =  Ene_Eksp_PI(i) + Epkn_ol(l) * Trp

     end do

   end do

   Ene_Eksp_PI = Ene_Eksp_PI - Ene_k_ol

   deallocate(Rdum)

!----------------------------------------------------------------------

   Ene_Eksp_PI = Ene_Eksp_PI / Volume * InvPi


end subroutine Ene_kspace_CBPI

! <<
