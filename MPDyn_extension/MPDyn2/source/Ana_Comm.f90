! ############################
! ## SUBROUTINE LIST 
! ## -- CArea 
! ## -- CellTransform 
! ## -- Transform 
! ## -- COM_SpecComponent 
! ## -- MSD_Diffusion 
! ## -- MSD_DiffusionCom 
! ## -- GR 
! ## -- RadialDist 
! ## -- GRG 
! ## -- RadialDistGG
! ## -- IntraMolDist 
! ## -- CorrRotationW 
! ## -- CavityDist 
! ## -- RotRing5 
! ############################


!######################################################################
!######################################################################


subroutine CArea(x,y,dz)

implicit none

real(8), dimension(3) :: x, y
real(8) :: dz, LenA, LenB, csAB

   LenA = sqrt( dot_product(x,x) )
   LenB = sqrt( dot_product(y,y) )
   csAB = dot_product(x,y)/(LenA*LenB)
   dz   = LenA * LenB * sin( acos( csAB ) )

end subroutine CArea


!######################################################################
!######################################################################


subroutine CellTransform

use Numbers, only : N
use Configuration, only : R
use CellParam, only : H, InvH

implicit none

integer :: i
real(8), dimension(3,3) :: A1, A2, A3, Rott
real(8), dimension(3) :: a, b, ab
real(8) :: cst, snt, csp, snp
real(8), dimension(3,N) :: ScR
real :: cs

   call InversMatrix(H,InvH)

   if(H(1,2)==0..and.H(1,3)==0..and.H(2,3)==0.) Return

   do i = 1 , N
     ScR(:,i) = matmul( InvH , R(:,i) )
   end do

   a  = H(:,1)
   b  = H(:,2)
   call VecProduct(a,b,ab)

   cst = ab(3)
   cs  = real(cst)

   if(cs/=1.) then
     snt = sqrt( 1.d0 - cst * cst )
     csp = ab(1) / snt
     snp = ab(2) / snt

     Rott = 0.d0
     Rott(1,1) =  csp
     Rott(1,2) =  snp
     Rott(2,1) = -snp
     Rott(2,2) =  csp
     Rott(3,3) = 1.d0

     A1 = matmul( Rott, H )

     Rott = 0.d0
     Rott(1,1) =  cst
     Rott(1,3) = -snt
     Rott(2,2) = 1.d0
     Rott(3,1) =  snt
     Rott(3,3) =  cst

     A2 = matmul( Rott, A1 )

     Rott = 0.d0
     Rott(1,1) =  csp
     Rott(1,2) = -snp
     Rott(2,1) =  snp
     Rott(2,2) =  csp
     Rott(3,3) = 1.d0

     A3 = matmul( Rott, A2 )

     do i = 1 , N
       R(:,i) = matmul( A3 , ScR(:,i) )
     end do

     H = A3

     call InversMatrix(H,InvH)

   end if

end subroutine CellTransform

! #######

subroutine Transform

use Numbers, only : N
use Configuration, only : R
use CellParam, only : H, InvH

implicit none

integer :: i
real(8), dimension(3,N) :: ScR
real(8), dimension(3,3) :: Tfm, Ht, Ht1, Ht2, Ht3
real(8), dimension(3) :: Va
real(8) :: ct, st, cp, sp
real(8) :: ABx, ABy, ABz, Z, La
real :: cs

   Ht = H

   call InversMatrix(H,InvH)

   do i = 1 , N
     ScR(:,i) = matmul( InvH, R(:,i) )
   end do

   ABx = H(2,1)*H(3,2) - H(3,1)*H(2,2)
   ABy = H(3,1)*H(1,2) - H(1,1)*H(3,2)
   ABz = H(1,1)*H(2,2) - H(2,1)*H(1,2)
   Z = sqrt(ABx*ABx + ABy*ABy + ABz*ABz)
   ABx = ABx / Z
   ABy = ABy / Z
   ABz = ABz / Z

   ct = ABz

   cs = real(ct)

   if(cs/=1.) then

     st = sqrt(1.d0 - ABz*ABz)
     cp = ABx / sqrt(1.d0 - ABz*ABz)
     sp = ABy / sqrt(1.d0 - ABz*ABz)

! ##

     Tfm(1,1) = cp
     Tfm(1,2) = sp
     Tfm(1,3) = 0.d0
     Tfm(2,1) = -sp
     Tfm(2,2) = cp
     Tfm(2,3) = 0.d0
     Tfm(3,1) = 0.d0
     Tfm(3,2) = 0.d0
     Tfm(3,3) = 1.d0

     call matmat(Tfm,Ht,Ht1)

! ##

     Tfm(1,1) = ct
     Tfm(1,2) = 0.d0
     Tfm(1,3) = -st
     Tfm(2,1) = 0.d0
     Tfm(2,2) = 1.d0
     Tfm(2,3) = 0.d0
     Tfm(3,1) = st
     Tfm(3,2) = 0.d0
     Tfm(3,3) = ct

     call matmat(Tfm,Ht1,Ht2)

! ##

     Tfm(1,1) = cp
     Tfm(1,2) = -sp
     Tfm(1,3) = 0.d0
     Tfm(2,1) = sp
     Tfm(2,2) = cp
     Tfm(2,3) = 0.d0
     Tfm(3,1) = 0.d0
     Tfm(3,2) = 0.d0
     Tfm(3,3) = 1.d0

     call matmat(Tfm,Ht2,Ht3)

   else

     Ht3 = H

   end if

! ##

   Va = Ht3(:,1)

   La = sqrt( dot_product( Va, Va ) )
   cp = Va(1) / La
   sp = Va(2) / La

   Tfm(1,1) = cp
   Tfm(1,2) = sp
   Tfm(1,3) = 0.d0
   Tfm(2,1) = -sp
   Tfm(2,2) = cp
   Tfm(2,3) = 0.d0
   Tfm(3,1) = 0.d0
   Tfm(3,2) = 0.d0
   Tfm(3,3) = 1.d0

   call MatMat(Tfm,Ht3,H)

   do i = 1 , N
     R(:,i) = matmul( H, ScR(:,i) )
   end do

   call InversMatrix(H,InvH)

end subroutine Transform


!######################################################################
!######################################################################


subroutine COM_SpecComponent(Rg,j)

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use AtomParam, only : Mass

implicit none

real(8), dimension(3) :: Rg
real(8) :: SumM
integer :: i, j, Numa, Numb

   Numa = 0
   Numb = NumMol(1)*NumAtm(1)

   if(j/=1) then
     do i = 2 , j
       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)
     end do
   end if

   Numa = Numa + 1

   Rg   = 0.d0
   SumM = 0.d0

   do i = Numa , Numb

     Rg   = Rg   + Mass(i) * R(:,i)
     SumM = SumM + Mass(i)

   end do

   Rg = Rg / SumM

end subroutine COM_SpecComponent


!######################################################################
!######################################################################


subroutine MSD_Diffusion

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QMaster
use ParamAnalyze
use CellParam, only : H, InvH
use AtomParam, only : AtomName, ResidName
use TimeParam, only : Timeps

implicit none

integer :: i, j, k, l, ii
integer :: lbb, lrun, TotalStepNumber
character(len=4) :: RName, AName
integer, dimension(NumMSD) :: NumMSDPart
integer, dimension(:,:), allocatable :: MSDnumatom
real, dimension(:,:,:,:), allocatable :: Si
integer :: NumF, MaxMSDPart
real, dimension(:,:,:), allocatable :: Ht
real(8), dimension(3) :: Stmp
integer :: ix, iy

   NumMSDPart = 0

   if(MSDdim==2) then
     ix = MSDdir(1)
     iy = MSDdir(2)
   else if(MSDdim==1) then
     ix = MSDdir(1)
   end if

   do i = 1 , NumMSD

     RName = MSDResi(i)
     AName = MSDAtom(i)

     do j = 1 , N

       if((ResidName(j) == RName).and.(AtomName(j) == AName)) then

         NumMSDPart(i) = NumMSDPart(i) + 1

       end if

     end do

     print *, NumMSDPart(i)

   end do

   MaxMSDPart = NumMSDPart(1)

   do i = 2 , NumMSD

     if( MaxMSDPart < NumMSDPart(i)) then
       MaxMSDPart = NumMSDPart(i)
     end if

   end do

   print *, 'maxmsdpart=',MaxMSDPart

   allocate( MSDnumatom(MaxMSDPart,NumMSD) )

   do i = 1 , NumMSD

     RName = MSDResi(i)
     AName = MSDAtom(i)

     k = 0

     do j = 1 , N

       if((ResidName(j) == RName).and.(AtomName(j) == AName)) then

         k = k + 1
         MSDnumatom(k,i) = j

       end if

     end do

     if(k/=NumMSDPart(i)) then
       if(QMaster) write(*,*) 'FATAL ERROR'
       call Finalize
     end if

   end do

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

   if(mod(Nsnap,BlockStep) /= 0) then
     write(*,*) 'ERROR : BlockStep'
     write(*,*) Nsnap, BlockStep
     call Finalize
   end if

!   nrun = Nsnap / BlockStep
   lrun = (BlockStep - ilength) / interv + 1
   lbb  = (BlockStep - ilength) + 1

   allocate( Si(3,MaxMSDPart,NumMSD,BlockStep) )
   allocate( Ht(3,3,BlockStep) )

   do i = 1, NumMSD

     open(i+50, file=trim(MSDfile(i)),status='unknown')

   end do

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

       call InversMatrix(H,InvH)

       NumF = NumF + 1

       Ht(:,:,NumF)    = sngl(H(:,:))

       do k = 1, NumMSD

         do l = 1, NumMSDPart(k)

           ii = MSDnumatom(l,k)

           Stmp = matmul( InvH, R(:,ii) )

           Si(:,l,k,NumF) = sngl(Stmp)

         end do

       end do

       if(NumF == BlockStep) then
         call MSD_atom
         NumF = 0
       end if

     end do

   end do


Contains

   subroutine MSD_atom

   implicit none

   integer :: imsd, jmsd
   integer :: ini, jj
   real, dimension(3,ilength) :: Rw
   real, dimension(3) :: Ri, Sij
   real, dimension(3,3) :: Hs
   real, dimension(NumMSD,ilength) :: sd
   real :: R2
   real :: pref

     sd = 0.d0

     do imsd = 1 , NumMSD

       do jmsd = 1 , NumMSDPart(imsd)

         do ini = 1, lbb, interv

           Hs(:,:) = Ht(:,:,ini)

           Rw(:,1) = matmul( Hs, Si(:,jmsd,imsd,ini) )

           do jj = 1, ilength-1

             Sij = Si(:,jmsd,imsd,jj+ini) - Si(:,jmsd,imsd,jj+ini-1)
             Sij = Sij - nint( Sij )

             Ri(1) = Ht(1,1,jj+ini)*Sij(1) + Ht(1,2,jj+ini)*Sij(2) + Ht(1,3,jj+ini)*Sij(3)
             Ri(2) = Ht(2,1,jj+ini)*Sij(1) + Ht(2,2,jj+ini)*Sij(2) + Ht(2,3,jj+ini)*Sij(3)
             Ri(3) = Ht(3,1,jj+ini)*Sij(1) + Ht(3,2,jj+ini)*Sij(2) + Ht(3,3,jj+ini)*Sij(3)

             Rw(:,jj+1) = Rw(:,jj) + Ri(:)

           end do

           do jj = 2, ilength

             Ri(:) = Rw(:,jj) - Rw(:,1)
             if(MSDdim==3) then
               R2 = dot_product( Ri, Ri )
             else if(MSDdim==2) then
               R2 = Ri(ix)*Ri(ix) + Ri(iy)*Ri(iy)
             else if(MSDdim==1) then
               R2 = Ri(ix)*Ri(ix)
             end if

             sd(imsd,jj-1) = sd(imsd,jj-1) + R2

           end do

         end do

       end do

       timeps = 0.d0

       pref = 1.d0 / dble(lrun*NumMSDPart(imsd))

       write(imsd+50,'(a)') '# mean square displacement(MSD) '
       write(imsd+50,'(a)') '# time[ps]    MSD[A^2] '

       write(imsd+50,'(5x,a6,11x,a7)') '0.0000','0.00000'

       do jj = 1, ilength - 1

         timeps = dtime * jj

         write(imsd+50,'(f9.3,2f10.4)') timeps, sd(imsd,jj)*pref

       end do

       write(imsd+50,'(a)') ' = = = '

     end do

   end subroutine MSD_atom

end subroutine MSD_Diffusion


!######################################################################
!######################################################################


subroutine MSD_DiffusionCom

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H, InvH
use AtomParam, only : MolName, Mass
use TimeParam, only : Timeps

implicit none

integer :: i, j, k, l, ii, jj, kk
integer :: lbb, lrun
character(len=8) :: MName
integer, dimension(NumMSD) :: CompNumb, NumMolStart
real(8), dimension(:,:,:), allocatable :: Si, Sj, Sk, Sl, Sn
real(8), dimension(:,:,:), allocatable :: Ht
real(8), dimension(:,:), allocatable :: disXYZ
real(8), dimension(:), allocatable :: MassMol
integer :: NumF
real(8), dimension(3) :: Rg, Sg
integer :: ix, iy

   if(MSDdim==2) then
     ix = MSDdir(1)
     iy = MSDdir(2)
   else if(MSDdim==1) then
     ix = MSDdir(1)
   end if

   do i = 1 , NumMSD

     MName = MSDMol(i)

     ii = 0

     do j = 1 , NumSpec

       if(MolName(j) == MName) then

         CompNumb(i) = j
         NumMolStart(i) = ii

       end if

       ii = ii + NumMol(j)*NumAtm(j)

     end do

   end do

   allocate( MassMol(NumMSD) )

   do i = 1, NumMSD

     MassMol(i) = 0.d0

     do j = NumMolStart(i) + 1, NumMolStart(i) + NumAtm(CompNumb(i))

       MassMol(i) = MassMol(i) + Mass(j)

     end do

   end do

   if(mod(Nsnap,BlockStep) /= 0) then
     write(*,*) 'ERROR : BlockStep'
     write(*,*) Nsnap, BlockStep
     call Finalize
   end if

!   nrun = Nsnap / BlockStep
   lrun = (BlockStep - ilength) / interv + 1
   lbb  = (BlockStep - ilength) + 1

   allocate( Ht(3,3,BlockStep) )

   allocate( Si(3,NumMol(CompNumb(1)),BlockStep) )
   if(NumMSD>=2) then
     allocate( Sj(3,NumMol(CompNumb(2)),BlockStep) )
   end if
   if(NumMSD>=3) then
     allocate( Sk(3,NumMol(CompNumb(3)),BlockStep) )
   end if
   if(NumMSD>=4) then
     allocate( Sl(3,NumMol(CompNumb(4)),BlockStep) )
   end if
   if(NumMSD>=5) then
     allocate( Sn(3,NumMol(CompNumb(5)),BlockStep) )
   end if
   if(NumMSD>=6) then
     write(*,*) 'ERROR : too many particles to be analyzed at the same time'
     call Finalize
   end if

   allocate( disXYZ(NumMSD,ilength) )

   do i = 1, NumMSD

     open(i, file=trim(MSDfile(i)),status='unknown')

   end do

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

       if(mod(j,Interval(i))/=0) cycle

       call InversMatrix(H,InvH)

       NumF = NumF + 1

       Ht(:,:,NumF)    = H(:,:)

       do k = 1, NumMSD

         l = CompNumb(k)

         kk = NumMolStart(k)

         do ii = 1, NumMol(l)

           Rg = 0.d0

           do jj = 1, NumAtm(l)

             kk = kk + 1

             Rg = Rg + Mass(kk) * R(:,kk)

           end do

           Rg = Rg / MassMol(k)

           Sg = matmul( InvH, Rg )

           if(k==1) then

             Si(:,ii,NumF) = Sg

           else if(k==2) then

             Sj(:,ii,NumF) = Sg

           else if(k==3) then

             Sk(:,ii,NumF) = Sg

           else if(k==4) then

             Sl(:,ii,NumF) = Sg

           else if(k==5) then

             Sn(:,ii,NumF) = Sg

           end if

         end do

       end do


       if(NumF == BlockStep) then
         call MSD_com
         NumF = 0
       end if

     end do

   end do

Contains

   subroutine MSD_com

   implicit none

   integer :: ini, nn
   real(8), dimension(3,ilength) :: Rw
   real(8), dimension(3) :: Sij, Ri, Rf
   real(8) :: pref
   real(8) :: R2

     disXYZ = 0.d0

     do nn = 1, NumMSD

       do ini = 1, lbb, interv

         do k = 1, NumMol(CompNumb(nn))

           H(:,:) = Ht(:,:,ini)

           select case(nn)

           case(1)

             Rw(:,1) = matmul( H, Si(:,k,ini) )

             do ii = ini+1, ini+ilength-1

               Rf(:) = Rw(:,ii-ini)

               Sij = Si(:,k,ii) - Si(:,k,ii-1)
               Sij = Sij - nint( Sij )
               Ri(1) = Ht(1,1,ii)*Sij(1) + Ht(1,2,ii)*Sij(2) + Ht(1,3,ii)*Sij(3)
               Ri(2) = Ht(2,1,ii)*Sij(1) + Ht(2,2,ii)*Sij(2) + Ht(2,3,ii)*Sij(3)
               Ri(3) = Ht(3,1,ii)*Sij(1) + Ht(3,2,ii)*Sij(2) + Ht(3,3,ii)*Sij(3)

               Rw(:,ii-ini+1) = Rf(:) + Ri(:)

             end do

           case(2)

             Rw(:,1) = matmul( H, Sj(:,k,ini) )

             do ii = ini+1, ini+ilength-1

               Rf(:) = Rw(:,ii-ini)

               Sij = Sj(:,k,ii) - Sj(:,k,ii-1)
               Sij = Sij - nint( Sij )
               Ri(1) = Ht(1,1,ii)*Sij(1) + Ht(1,2,ii)*Sij(2) + Ht(1,3,ii)*Sij(3)
               Ri(2) = Ht(2,1,ii)*Sij(1) + Ht(2,2,ii)*Sij(2) + Ht(2,3,ii)*Sij(3)
               Ri(3) = Ht(3,1,ii)*Sij(1) + Ht(3,2,ii)*Sij(2) + Ht(3,3,ii)*Sij(3)

               Rw(:,ii-ini+1) = Rf(:) + Ri(:)

             end do

           case(3)

             Rw(:,1) = matmul( H, Sk(:,k,ini) )

             do ii = ini+1, ini+ilength-1

               Rf(:) = Rw(:,ii-ini)

               Sij = Sk(:,k,ii) - Sk(:,k,ii-1)
               Sij = Sij - nint( Sij )
               Ri(1) = Ht(1,1,ii)*Sij(1) + Ht(1,2,ii)*Sij(2) + Ht(1,3,ii)*Sij(3)
               Ri(2) = Ht(2,1,ii)*Sij(1) + Ht(2,2,ii)*Sij(2) + Ht(2,3,ii)*Sij(3)
               Ri(3) = Ht(3,1,ii)*Sij(1) + Ht(3,2,ii)*Sij(2) + Ht(3,3,ii)*Sij(3)

               Rw(:,ii-ini+1) = Rf(:) + Ri(:)

             end do

           case(4)

             Rw(:,1) = matmul( H, Sl(:,k,ini) )

             do ii = ini+1, ini+ilength-1

               Rf(:) = Rw(:,ii-ini)

               Sij = Sl(:,k,ii) - Sl(:,k,ii-1)
               Sij = Sij - nint( Sij )
               Ri(1) = Ht(1,1,ii)*Sij(1) + Ht(1,2,ii)*Sij(2) + Ht(1,3,ii)*Sij(3)
               Ri(2) = Ht(2,1,ii)*Sij(1) + Ht(2,2,ii)*Sij(2) + Ht(2,3,ii)*Sij(3)
               Ri(3) = Ht(3,1,ii)*Sij(1) + Ht(3,2,ii)*Sij(2) + Ht(3,3,ii)*Sij(3)

               Rw(:,ii-ini+1) = Rf(:) + Ri(:)

             end do

           case(5)

             Rw(:,1) = matmul( H, Sn(:,k,ini) )

             do ii = ini+1, ini+ilength-1

               Rf(:) = Rw(:,ii-ini)

               Sij = Sn(:,k,ii) - Sn(:,k,ii-1)
               Sij = Sij - nint( Sij )
               Ri(1) = Ht(1,1,ii)*Sij(1) + Ht(1,2,ii)*Sij(2) + Ht(1,3,ii)*Sij(3)
               Ri(2) = Ht(2,1,ii)*Sij(1) + Ht(2,2,ii)*Sij(2) + Ht(2,3,ii)*Sij(3)
               Ri(3) = Ht(3,1,ii)*Sij(1) + Ht(3,2,ii)*Sij(2) + Ht(3,3,ii)*Sij(3)

               Rw(:,ii-ini+1) = Rf(:) + Ri(:)

             end do

           end select

           do ii = 2, ilength

             Ri(:) = Rw(:,ii) - Rw(:,1)
             if(MSDdim==3) then
               R2 = dot_product( Ri, Ri )
             else if(MSDdim==2) then
               R2 = Ri(ix)*Ri(ix) + Ri(iy)*Ri(iy)
             else if(MSDdim==1) then
               R2 = Ri(ix)*Ri(ix)
             end if
             disXYZ(nn,ii-1) = disXYZ(nn,ii-1) + R2

           end do

         end do

       end do

       Timeps = 0.d0

       pref = 1.d0 / dble(lrun*NumMol(CompNumb(nn)))

       write(nn,'(a)') '# mean square displacement(MSD) '
       write(nn,'(a)') '# time[ps]    MSD[A^2] '

       write(nn,'(3x,a6,4x,a6)') ' 0.000','0.0000'

       do jj = 1, ilength - 1

         timeps = dtime * jj

         write(nn,'(f9.3,f10.4)') timeps, disXYZ(nn,jj)*pref

       end do

       write(nn,'(a)') ' = = = '

     end do

   end subroutine MSD_com

end subroutine MSD_DiffusionCom


!######################################################################
!######################################################################


subroutine GR

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : R
use CommonBlocks, only : QMaster
use ParamAnalyze
use UnitExParam, only : pi
use CellParam, only : H, InvH, Volume
use CutoffParam, only : Rcutoff2
use AtomParam, only : AtomName, ResidName

implicit none

integer :: i, j, k, l, ii, NumChar
integer :: TotalStep
character(len=4) :: AIName, RIName
character(len=4) :: AJName, RJName
character(len=4) :: AName
real(8) :: Vol, det, qk, gg, drr, crdn, xx
external det

   IRcut = int(sqrt(Rcutoff2)/DRgrid)

   allocate( igr(IRcut,NumGR) )
   allocate( NumGR_I(NumGR) )
   allocate( NumGR_J(NumGR) )
   allocate( ListGR_I(N,NumGR) )
   allocate( ListGR_J(N,NumGR) )
   if(QIntraSubt) then
     allocate( NumGR_Subt(NumGR) )
     allocate( MoleculeID(N) )

     ii = 0
     l  = 0
     do i = 1, NumSpec
       do j = 1, NumMol(i)
         l = l + 1
         do k = 1, NumAtm(i)
           ii = ii + 1
           MoleculeID(ii) = l
         end do
       end do
     end do

   end if

   igr = 0

!
! ## pick up the specified atom pairs
!
   NumGR_I = 0
   NumGR_J = 0

   do j = 1 , NumGR

     AIName = GRAtomI(j)
     RIName = GRResiI(j)

     NumChar = len( trim(AIName) )

     if( NumChar == 1 ) then

       do i = 1 , N

         if( ResidName(i) == RIName ) then

           Aname = AtomName(i)

           if( Aname(1:1) == trim(AIName) ) then

             NumGR_I(j) = NumGR_I(j) + 1
             ListGR_I(NumGR_I(j),j) = i

           end if

         end if

       end do

     else

       do i = 1 , N

         if( ResidName(i) == RIName ) then

           Aname = AtomName(i)

           if( Aname == AIName ) then

             NumGR_I(j) = NumGR_I(j) + 1
             ListGR_I(NumGR_I(j),j) = i

           end if

         end if

       end do

     end if

     AJName = GRAtomJ(j)
     RJName = GRResiJ(j)

     NumChar = len( trim(AJName) )

     if( NumChar == 1 ) then

       do i = 1 , N

         if( ResidName(i) == RJName ) then

           Aname = AtomName(i)

           if( Aname(1:1) == trim(AJName) ) then

             NumGR_J(j) = NumGR_J(j) + 1
             ListGR_J(NumGR_J(j),j) = i

           end if

         end if

       end do

     else

       do i = 1 , N

         if( ResidName(i) == RJName ) then

           Aname = AtomName(i)

           if( Aname == AJName ) then

             NumGR_J(j) = NumGR_J(j) + 1
             ListGR_J(NumGR_J(j),j) = i

           end if

         end if

       end do

     end if

     if(QIntraSubt) then
       l  = 0
       if(RIName==RJName) then
         ii = 0
         do i = 1, NumSpec
           do k = 1, NumAtm(i)
             if(ResidName(ii+k)==RJName) then
               NumChar = len( trim(AJName) )
               Aname = AtomName(ii+k)
               if(NumChar==1) then
                 if( Aname(1:1) == trim(AJName) ) then
                   l = l + 1
                 end if
               else
                 if( Aname == AJName ) then
                   l = l + 1
                 end if
               end if
             end if
           end do
           ii = ii + NumAtm(i)*NumMol(i)
         end do
       end if
       NumGR_Subt(j) = l
     end if

     if(NumGR_I(j)==0) then
       write(*,'(a,i5)') "ERROR: no atom (I) was found for PAIR",j
       call Finalize
     else if(NumGR_J(j)==0) then
       write(*,'(a,i5)') "ERROR: no atom (J) was found for PAIR",j
       call Finalize
     end if

     print *,"NUMGR=",NumGR_I(j),NumGR_J(j)

   end do

! -------------------------------------------------

   TotalStep = 0

   Vol = 0.d0

   do i = 1 , NJobs

     if(QMaster) call OpenTraj(i)

     TotalStep = TotalStep + NTrjStep(i)

     do j = 1 , NTrjStep(i)

!     ------------------------------
#ifdef MOLFILE
       if(QMaster) call Read_RTraj(i)
#else
       if(QMaster) call Read_RTraj
#endif
!     ------------------------------

!     ---------------
       call BcastRH
!     ---------------

       call InversMatrix(H,InvH)

!     -----------------
       call RadialDist
!     -----------------

       Volume = det(H)
       Vol = Vol + Volume

     end do

   end do

   call SumRDF( igr )

   if(QMaster) then

     Vol = Vol / dble(TotalStep)

     do i = 1, NumGR

       open(1,file = trim(GRfile(i)), status='unknown')

       write(1,'(a)') '# radial distribution function g(r) '
       write(1,'(a)') '# distance[A]  g(r)[-]  Coordination number[-]'

       crdn = 0.d0
       xx   = 1.d0 / dble(TotalStep*NumGR_I(i))

       do j = 1 , IRcut

         drr = (j-0.5d0) * DRgrid
         if(QIntraSubt) then
           qk = Vol / ((NumGR_J(i)-NumGR_Subt(i)) * 4.d0 * pi * drr * drr * DRgrid)
         else
           qk = Vol / (NumGR_J(i) * 4.d0 * pi * drr * drr * DRgrid)
         end if

         crdn = crdn + dble(igr(j,i)) * xx
         gg = dble( igr(j,i) ) * qk * xx

         write(1,'(f7.2,f12.6,e15.6)') drr, gg, crdn

       end do

       close(1)

     end do

   end if

end subroutine GR


!######################################################################
!######################################################################


subroutine RadialDist

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use CommonMPI
use CellParam, only : H, InvH
use CutoffParam, only : Rcutoff2

integer :: i , j, k, ii, jj
integer :: Nas, ir
real(8), dimension(3) :: Si, Rij, Sij
real(8), dimension(3,N) :: ScR
real(8) :: R2, R1, InvDR

   InvDR = 1.d0 / DRgrid

   Nas = NProcs - MyRank

   do i = 1 , N

     ScR(:,i) = matmul( InvH , R(:,i) )

   end do


   do k = 1 , NumGR

     do ii = Nas , NumGR_I(k), NProcs

       i  = ListGR_I(ii,k)
       Si = ScR(:,i)

       do jj = 1 , NumGR_J(k)

         j = ListGR_J(jj,k)

         if( i == j ) cycle
         if(QIntraSubt) then
           if( MoleculeID(i) == MoleculeID(j) ) cycle
         end if

         Sij = Si - ScR(:,j)
         Sij = Sij - nint(Sij)
         Rij = matmul( H, Sij )

         R2 = dot_product( Rij, Rij )

         if(R2 < Rcutoff2) then

           R1 = sqrt(R2)
           ir = int(R1*InvDR) + 1
           igr(ir,k) = igr(ir,k) + 1

         end if

       end do

     end do

   end do


end subroutine RadialDist


!######################################################################
!######################################################################


subroutine GRG

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : R
use CommonBlocks, only : QMaster
use ParamAnalyze
use UnitExParam, only : pi
use CellParam, only : H, InvH, Volume
use CutoffParam, only : Rcutoff2
use AtomParam, only : AtomName, ResidName

implicit none

integer :: i, j, ii, NmolMax
integer :: TotalStep
character(len=4) :: RIName, NameI, NameJ
logical :: MatchI, MatchJ
real(8) :: Vol, det, qk, gg, drr, crdn, xx
external det

   IRcut = int(sqrt(Rcutoff2)/DRgrid)

   allocate( igr(IRcut,NumGR) )
   allocate( NumGR_I(NumGR) )
   allocate( NumGR_J(NumGR) )
   allocate( ListGR_I(N,NumGR) )
   allocate( ListGR_J(N,NumGR) )

   igr = 0

!
! ## pick up the specified atom pairs
!
   NmolMax = NumMol(1)
   if(NumSpec>1) then
     do i = 2, NumSpec
       if(NmolMax<NumMol(i)) NmolMax = NumMol(i)
     end do
   end if

   do j = 1, NumGR
     NameI = GRResiI(j)
     NameJ = GRResiJ(j)
     MatchI = .True.
     MatchJ = .True.

     ii = 0
     do i = 1, NumSpec
       RIName = ResidName(ii+1)
       if( NameI == RIName ) then
         NumGR_I(j) = i
         MatchI= .False.
       end if
       if( NameJ == RIName ) then
         NumGR_J(j) = i
         MatchJ= .False.
       end if
       if((i==NumSpec).and.(MatchI.or.MatchJ)) then
         write(*,*) "error: no matching name is available for pair molecule"
         if(MatchI) write(*,*) NameI
         if(MatchJ) write(*,*) NameJ
         call Finalize
       end if
       ii = ii + NumMol(i)*NumAtm(i)
     end do

   end do

! -------------------------------------------------

   TotalStep = 0

   Vol = 0.d0

   do i = 1 , NJobs

     if(QMaster) call OpenTraj(i)

     TotalStep = TotalStep + NTrjStep(i)

     do j = 1 , NTrjStep(i)

!     ------------------------------
#ifdef MOLFILE
       if(QMaster) call Read_RTraj(i)
#else
       if(QMaster) call Read_RTraj
#endif
!     ------------------------------

!     ---------------
       call BcastRH
!     ---------------

       call InversMatrix(H,InvH)

!     ----------------------------
       call RadialDistGG(NmolMax)
!     ----------------------------

       Volume = det(H)
       Vol = Vol + Volume

     end do

   end do

   call SumRDF( igr )

   if(QMaster) then

     Vol = Vol / dble(TotalStep)

     do i = 1, NumGR

       open(1,file = trim(GRfile(i)), status='unknown')

       write(1,'(a)') '# radial distribution function g(r) '
       write(1,'(a)') '# distance[A]  g(r)[-]  Coordination number[-]'

       crdn = 0.d0
       xx   = 1.d0 / dble(TotalStep*NumMol(NumGR_I(i)))

       do j = 1 , IRcut

         drr = (j-0.5d0) * DRgrid
         qk = Vol / (NumMol(NumGR_J(i)) * 4.d0 * pi * drr * drr * DRgrid)

         crdn = crdn + dble(igr(j,i)) * xx
         gg = dble( igr(j,i) ) * qk * xx

         write(1,'(f7.2,f12.6,e15.6)') drr, gg, crdn

       end do

       close(1)

     end do

   end if

end subroutine GRG


!######################################################################
!######################################################################


subroutine RadialDistGG(NmolMax)

use Numbers, only : N, NumSpec, NumMol, NumAtm
use AtomParam, only : Mass
use Configuration, only : R
use ParamAnalyze
use CommonMPI
use CellParam, only : H, InvH
use CutoffParam, only : Rcutoff2

integer :: i , j, k, ii, NmolMax
integer :: Nas, ir
real(8), dimension(3) :: Si, Rij, Sij
real(8) :: R2, R1, InvDR
real(8), dimension(3,NmolMax,NumSpec) :: Rcom, ScR
real(8), dimension(NumSpec) :: M_Mol

   InvDR = 1.d0 / DRgrid

   Nas = NProcs - MyRank

   ii = 0
   do i = 1, NumSpec

     M_Mol = 0.d0

     do j = 1, NumAtm(i)
       M_Mol = M_Mol + Mass(ii+j)
     end do

     M_Mol = 1.d0 / M_Mol

     do j = 1, NumMol(i)
       Rcom(:,j,i) = 0.d0
       do k = 1, NumAtm(i)
         ii = ii + 1
         Rcom(:,j,i) = Rcom(:,j,i) + Mass(ii) * R(:,ii)
       end do
       Rcom(:,j,i) = Rcom(:,j,i)*M_Mol
     end do

   end do

   do i = 1 , NumSpec
     do j = 1, NumMol(i)

       ScR(:,j,i) = matmul( InvH , Rcom(:,j,i) )

     end do
   end do

!--------------

   do k = 1 , NumGR

     i1 = NumGR_I(k)
     i2 = NumGR_J(k)

     do i = Nas , NumMol(i1), NProcs

       Si = ScR(:,i,i1)

       do j = 1 , NumMol(i2)
         if((i1==i2).and.(i==j)) cycle

         Sij = Si - ScR(:,j,i2)
         Sij = Sij - nint(Sij)
         Rij = matmul( H, Sij )

         R2 = dot_product( Rij, Rij )

         if(R2 < Rcutoff2) then

           R1 = sqrt(R2)
           ir = int(R1*InvDR) + 1
           igr(ir,k) = igr(ir,k) + 1

         end if

       end do

     end do

   end do


end subroutine RadialDistGG


!######################################################################
!######################################################################


subroutine IntraMolDist

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use CommonBlocks, only : QMaster
use ParamAnalyze
use UnitExParam, only : InvPi
use AtomParam, only : MolName, AtomName

integer :: i, j, ii, ia, iia, ib, iib, id, iid, ic, imin, imax
integer :: it, ir
integer, dimension(20) :: MinNumbMol
integer :: MaxNumbMol
real(8) :: Ra2, Rb2
real(8) :: rRa2, rRb2, rRab
real(8) :: CsPhi, Phi
real(8) :: Ds, Dir
real(8) :: Deg, R1, R2, cst, th
real(8), dimension(3) :: Rij, Rkj, Rlj
real(8), dimension(3) :: Pla, Plb, Ori
character(len=90) :: Filename
integer, dimension(:,:,:), allocatable :: BDpair, ANpair, DHpair
integer, dimension(:), allocatable :: BDnum, ANnum, DHnum
real(8), dimension(:,:), allocatable :: BDdist, ANdist, DHdist
real(8), dimension(:), allocatable :: AveBD, Ave2BD
real(8), dimension(:), allocatable :: AveDH, Ave2DH
real(8), dimension(:), allocatable :: AveAN, Ave2AN

   ii = 0
   MaxNumbMol = 0
   do i = 1, NumSpec
     MinNumbMol(i) = ii
     ii = ii + NumMol(i)*NumAtm(i)
     if(NumMol(i)>MaxNumbMol) MaxNumbMol = NumMol(i)
   end do

   if(NBdAna/=0) then

     imin = nint(Rrangemin/DRgrid)
     imax = nint(Rrangemax/DRgrid)
     allocate(BDpair(2,MaxNumbMol,NBdAna))
     allocate(BDdist(imin:imax,NBdAna))
     allocate(BDnum(NBdAna))
     allocate(AveBD(NBdAna))
     allocate(Ave2BD(NBdAna))
     BDdist = 0.d0
     AveBD  = 0.d0
     Ave2BD = 0.d0

     do ib = 1, NBdAna
       do i = 1, NumSpec
         if(MolName(i)==BdMol(ib)) then
           ic = i
         end if
       end do
       do i = 1, NumMol(ic)
         do j = 1, NumAtm(ic)
           ii = MinNumbMol(ic) + (i-1)*NumAtm(ic) + j
           if(AtomName(ii)==BdAtomI(ib)) then
             BDpair(1,i,ib) = ii
           else if(AtomName(ii)==BdAtomJ(ib)) then
             BDpair(2,i,ib) = ii
           end if
         end do
       end do
       BDnum(ib) = NumMol(ic)
     end do

   end if

   if(NAnAna/=0) then

     allocate(ANpair(3,MaxNumbMol,NAnAna))
     allocate(ANdist(180,NAnAna))
     allocate(ANnum(NAnAna))
     allocate(AveAN(NAnAna))
     allocate(Ave2AN(NAnAna))
     ANdist = 0.d0
     AveAN  = 0.d0
     Ave2AN = 0.d0

     do ia = 1, NAnAna
       do i = 1, NumSpec
         if(MolName(i)==AnMol(ia)) then
           ic = i
         end if
       end do
       do i = 1, NumMol(ic)
         do j = 1, NumAtm(ic)
           ii = MinNumbMol(ic) + (i-1)*NumAtm(ic) + j
           if(AtomName(ii)==AnAtomI(ia)) then
             ANpair(1,i,ia) = ii
           else if(AtomName(ii)==AnAtomJ(ia)) then
             ANpair(2,i,ia) = ii
           else if(AtomName(ii)==AnAtomK(ia)) then
             ANpair(3,i,ia) = ii
           end if
         end do
       end do
       ANnum(ia) = NumMol(ic)
     end do

   end if

   if(NDhAna/=0) then

     allocate(DHpair(4,MaxNumbMol,NDhAna))
     allocate(DHdist(-180:180,NDhAna))
     allocate(DHnum(NDhAna))
     allocate(AveDH(NDhAna))
     allocate(Ave2DH(NDhAna))
     DHdist = 0.d0
     AveDH  = 0.d0
     Ave2DH = 0.d0

     do id = 1, NDhAna
       do i = 1, NumSpec
         if(MolName(i)==DhMol(id)) then
           ic = i
         end if
       end do
       do i = 1, NumMol(ic)
         do j = 1, NumAtm(ic)
           ii = MinNumbMol(ic) + (i-1)*NumAtm(ic) + j
           if(AtomName(ii)==DhAtomI(id)) then
             DHpair(1,i,id) = ii
           else if(AtomName(ii)==DhAtomJ(id)) then
             DHpair(2,i,id) = ii
           else if(AtomName(ii)==DhAtomK(id)) then
             DHpair(3,i,id) = ii
           else if(AtomName(ii)==DhAtomL(id)) then
             DHpair(4,i,id) = ii
           end if
         end do
       end do
       DHnum(id) = NumMol(ic)
     end do

   end if

! -------------------------------------------------
   TotalStep = 0

   do i = 1 , NJobs

     if(QMaster) call OpenTraj(i)

     TotalStep = TotalStep + NTrjStep(i)

     do j = 1 , NTrjStep(i)

!     ------------------------------
#ifdef MOLFILE
       if(QMaster) call Read_RTraj(i)
#else
       if(QMaster) call Read_RTraj
#endif
!     ------------------------------

       if(NBdAna/=0) then
         do ib = 1, NBdAna
           do iib = 1, BDnum(ib)
             i1 = BDpair(1,iib,ib)
             i2 = BDpair(2,iib,ib)
             Rij = R(:,i1) - R(:,i2)
             R2  = dot_product(Rij,Rij)
             R1  = sqrt(R2)
             ir  = nint(R1/DRgrid)
             if(ir<imin.or.ir>imax) then
               write(*,*) 'ERROR : out of range : Bond length', R1
               call Finalize
             end if
             BDdist(ir,ib) = BDdist(ir,ib) + 1.
             AveBD(ib)  = AveBD(ib)  + R1
             Ave2BD(ib) = Ave2BD(ib) + R2
           end do
         end do
       end if
       if(NAnAna/=0) then
         do ia = 1, NAnAna
           do iia = 1, ANnum(ia)
             i1 = ANpair(1,iia,ia)
             i2 = ANpair(2,iia,ia)
             i3 = ANpair(3,iia,ia)
             Rij = R(:,i1) - R(:,i2)
             Rkj = R(:,i3) - R(:,i2)
             R1  = dot_product(Rij,Rij)
             R2  = dot_product(Rkj,Rkj)
             cst = dot_product(Rij,Rkj)/sqrt(R1*R2)
             th  = acos(cst) * 180.d0 * InvPi
             it  = int(th) + 1
             ANdist(it,ia) = ANdist(it,ia) + 1
             AveAN(ia)  = AveAN(ia)  + th
             Ave2AN(ia) = Ave2AN(ia) + th*th
           end do
         end do
       end if
       if(NDhAna/=0) then
         do id = 1, NDhAna
           do iid = 1, DHnum(id)
             i1 = DHpair(1,iid,id)
             i2 = DHpair(2,iid,id)
             i3 = DHpair(3,iid,id)
             i4 = DHpair(4,iid,id)

             Rij = R(:,i1) - R(:,i2)
             Rkj = R(:,i3) - R(:,i2)
             Rlj = R(:,i4) - R(:,i2)

             Pla = VecProd(Rij,Rkj)
             Plb = VecProd(Rlj,Rkj)

             Ra2 = dot_product(Pla,Pla)
             Rb2 = dot_product(Plb,Plb)

             rRa2 = 1.d0 / Ra2
             rRb2 = 1.d0 / Rb2
             rRab = sqrt(rRa2*rRb2)
             CsPhi = dot_product(Pla,Plb) * rRab

             Ori = VecProd(Pla,Plb)
             Ds  = dot_product(Ori,Rkj)

             Dir = 1.d0
             if(Ds < 0.d0) Dir = -1.d0

             Phi = acos(CsPhi)
             Deg = Dir * Phi * 180.d0 * InvPi
             it  = nint(Deg)
             DHdist(it,id) = DHdist(it,id) + 1

             AveDH(id)  = AveDH(id)  + Deg
             Ave2DH(id) = Ave2DH(id) + Deg*Deg

           end do
         end do
       end if

     end do

   end do

   open(111,file='SummaryDISTinMOL.dat')

   do i = 1, NBdAna
     write(Filename,'(7a)') 'bondlength_',trim(BdAtomI(i)),'--',&
     &   trim(BdAtomJ(i)),'_',trim(BdMol(i)),'.dat'
     open(51,file=trim(Filename))

     AveBD(i)  = AveBD(i)  / dble(BDnum(i)*TotalStep)
     Ave2BD(i) = Ave2BD(i) / dble(BDnum(i)*TotalStep)
     Std = sqrt( Ave2BD(i) - AveBD(i)*AveBD(i) )

     write(111,'(/6a)') 'For Bond ',BdAtomI(i),':',BdAtomJ(i),' in ',BdMol(i)
     write(111,'(a,f15.8)') '    Ave. = ',AveBD(i)
     write(111,'(a,f15.8)') '    Std. = ',Std

     do j = imin, imax
       BDdist(j,i) = BDdist(j,i) / dble(BDnum(i)*TotalStep*DRgrid)!*100.
       write(51,'(f8.2,f15.10)') j*DRgrid, BDdist(j,i)
     end do

     close(51)
   end do

   do i = 1, NAnAna
     write(Filename,'(9a)') 'angle_',trim(AnAtomI(i)),'--',trim(AnAtomJ(i)),&
     &   '--',trim(AnAtomK(i)),'_',trim(AnMol(i)),'.dat'
     open(51,file=trim(Filename))

     AveAN(i)  = AveAN(i)  / dble(ANnum(i)*TotalStep)
     Ave2AN(i) = Ave2AN(i) / dble(ANnum(i)*TotalStep)
     Std = sqrt( Ave2AN(i) - AveAN(i)*AveAN(i) )

     write(111,'(/8a)') 'For Angle ',AnAtomI(i),':',AnAtomJ(i),':',AnAtomK(i),' in ',AnMol(i)
     write(111,'(a,f15.8)') '    Ave. = ',AveAN(i)
     write(111,'(a,f15.8)') '    Std. = ',Std

     do j = 1, 180
       ANdist(j,i) = ANdist(j,i) / dble(ANnum(i)*Totalstep)!*100.
       write(51,'(f8.2,f15.10)') j-0.5, ANdist(j,i)
     end do

     close(51)
   end do

   do i = 1, NDhAna
     write(Filename,'(11a)') 'dihedral_',trim(DhAtomI(i)),'--',trim(DhAtomJ(i)),&
     & '--',trim(DhAtomK(i)),'--',trim(DhAtomL(i)),'_',trim(DhMol(i)),'.dat'
     open(51,file=trim(Filename))

     AveDH(i)  = AveDH(i)  / dble(DHnum(i)*TotalStep)
     Ave2DH(i) = Ave2DH(i) / dble(DHnum(i)*TotalStep)
     Std = sqrt( Ave2DH(i) - AveDH(i)*AveDH(i) )

     write(111,'(/10a)') 'For Dihedral ',DhAtomI(i),':',DhAtomJ(i),':',&
     &                  DhAtomK(i),':',DhAtomL(i),' in ',DhMol(i)
     write(111,'(a,f15.8)') '    Ave. = ',AveDH(i)
     write(111,'(a,f15.8)') '    Std. = ',Std

     do j = -180, 180
       DHdist(j,i) = DHdist(j,i) / dble(DHnum(i)*Totalstep)!*100.
       write(51,'(f8.2,f15.10)') real(j), DHdist(j,i)
     end do

     close(51)
   end do

Contains

   function VecProd(x,y) Result(z)

   real(8), dimension(3) :: x, y, z

     z(1) = x(2) * y(3) - y(2) * x(3)
     z(2) = x(3) * y(1) - y(3) * x(1)
     z(3) = x(1) * y(2) - y(1) * x(2)

   end function VecProd


end subroutine IntraMolDist


!######################################################################
!######################################################################


subroutine CorrRotationW

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use CommonBlocks, only : QMaster
use ParamAnalyze
use CommonMPI
use CellParam, only : H, InvH
use AtomParam, only : ResidName

implicit none

integer :: i , j, k, l, ii, jj, Nw, nm, ita
integer :: i1, i2, i3, TotalStepNumber, NumF
real(8) :: cs, cst, cst2, R2
real(8), dimension(3) :: Rij
real(8), dimension(:,:,:), allocatable :: Si
real(8), dimension(:), allocatable :: Fcorr, Fcorr2

   ii = 1
   do i = 1 , NumSpec
     if(ResidName(ii)=='TIP3') then
       Nw = i
       nm = NumMol(i)
     end if
     ii = ii + NumAtm(i) * NumMol(i)
   end do

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)
   end do

   allocate( Si(3,nm,TotalStepNumber) )

   allocate( Fcorr(TotalStepNumber) )
   allocate( Fcorr2(TotalStepNumber) )

   NumF = 0

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

       call InversMatrix(H,InvH)

       NumF = NumF + 1

       jj = 0

       do k = 1, NumSpec

         if(k==Nw) then

           do ii = 1 , NumMol(k)
             i1 = jj + NumAtm(k)*(ii-1) + 1
             i2 = jj + NumAtm(k)*(ii-1) + 2
             i3 = jj + NumAtm(k)*(ii-1) + 3
             Rij = R(:,i2) + R(:,i3) - 2.d0 * R(:,i1)
             R2  = dot_product( Rij, Rij )
             Si(:,ii,NumF) = Rij / sqrt(R2)
           end do

           jj = jj + NumMol(k) * NumAtm(k)

         else

           jj = jj + NumMol(k) * NumAtm(k)

         end if

       end do

     end do

   end do

   do l = 1 , TotalStepNumber - 1

     ita = TotalStepNumber - l

     cst  = 0.d0
     cst2 = 0.d0

     do k = 1 , ita

       do j = 1 , nm

         cs   = dot_product( Si(:,j,k), Si(:,j,k+l) )
         cst  = cst + cs
         cst2 = cst2 + cs * cs

       end do

     end do

     Fcorr(l)  = cst  / dble(nm*ita)
     Fcorr2(l) = cst2 / dble(nm*ita)

   end do

   open(1,file= './Analy/CorrRotationW.dat')

   write(1,'(a)') '# Rotational correlation function '
   write(1,'(a)') '# Legendre polynomial 1 C1(t) = cos(theta)'
   write(1,'(a)') '# Legendre polynomial 2 C2(t) = 0.5*(3.d0*cos(theta)-1)'
   write(1,'(a)') '# cos(theta) = < u(t)u(0) >'
   write(1,'(a)') '# where u(t) is the dipole unit-vector at t'
   write(1,'(a)') '# Time[ps]    C1(t)    C2(t)'
   write(1,'(5x,a6,2(11x,a7))') '0.0000','0.00000','0.00000'

   do l = 1 , TotalStepNumber - 1
     write(1,'(1x,f10.4,2(1x,f17.5))') dtime*l,Fcorr(l),0.5d0*(3.d0*Fcorr2(l)-1.d0)
   end do

   close(1)

end subroutine CorrRotationW


!######################################################################
!######################################################################


subroutine CavityDist

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use CommonBlocks, only : QMaster
use ParamAnalyze
use CellParam, only : H

implicit none

character(len=72) :: FileName
integer :: i , j, k, l, ii, jj, NumF
integer :: lx, ly, lz, lzh, ix, iy, iz, Nx, Ny, Nz, TotalStepNumber
integer, dimension(:,:,:), allocatable :: Flcav, Flfil, FlcavAcc
!! real(8), dimension(:,:,:), allocatable :: Rhofil
real(8), dimension(3) :: Ri, Rext
real(8) :: xx, zz!, xi, yi, zi
real(8) :: Rxmaxh, Rymaxh, Rzmaxh

   lx    = int( Rxmax / R_cav )
   ly    = int( Rymax / R_cav )
   lz    = int( Rzmax / R_cav )

   Rxmax = dble(lx) * R_cav
   Rymax = dble(ly) * R_cav
   Rzmax = dble(lz) * R_cav

   Rxmaxh = Rxmax * 0.5d0
   Rymaxh = Rymax * 0.5d0
   Rzmaxh = Rzmax * 0.5d0

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)
   end do

   allocate( Flcav(lx,ly,lz) )
   allocate( FlcavAcc(lx,ly,lz) )
   allocate( Flfil(lx,ly,lz) )
!!   allocate( Rhofil(lx,ly,lz) )

   FlcavAcc = 0
   NumF = 0

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

!       call CellTransform
       call Transform

       Flcav = 0
       Flfil = 0
!!       Rhofil = 0.d0

       jj = 0
       do ii = 1 , NumSpec

         if(ii == Kcomp) then

           do k = 1, NumMol(ii)
             do l = 1 , NumAtm(ii)
               jj = jj + 1
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
!!                       Rhofil(Nx,Ny,Nz) = Mass(jj)

                     end if

                   end do
                 end do
               end do

             end do
           end do

         else

           do k = 1, NumMol(ii)
             do l = 1 , NumAtm(ii)
               jj = jj + 1
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

                       Flfil(Nx,Ny,Nz)  = 2
!!                       Rhofil(Nx,Ny,Nz) = Mass(jj)

                     end if

                   end do
                 end do
               end do

             end do
           end do

         end if

       end do

       Flfil = Flfil + Flcav

       Flcav = 0

       do Nx = 1, lx
         do Ny = 1, ly
           do Nz = 1, lz
             if(Flfil(Nx,Ny,Nz) == 0) then
               Flcav(Nx,Ny,Nz) = 1
             end if
           end do
         end do
       end do

       FlcavAcc = FlcavAcc + Flcav

       if(mod(j,Interval(i))==0) then

       NumF = NumF + 1

       if(NumF<10) then
         write(FileName,'(a,i1,a)') './Analy/distCav0000',NumF,'.dat'
       else if(NumF<100) then
         write(FileName,'(a,i2,a)') './Analy/distCav000', NumF,'.dat'
       else if(NumF<1000) then
         write(FileName,'(a,i3,a)') './Analy/distCav00',  NumF,'.dat'
       else if(NumF<10000) then
         write(FileName,'(a,i4,a)') './Analy/distCav0',   NumF,'.dat'
       else if(NumF<100000) then
         write(FileName,'(a,i5,a)') './Analy/distCav',    NumF,'.dat'
       else
         write(*,*) 'ERROR : too large number of files'
         call Finalize
       end if

       open(1,file=trim(FileName),status='unknown')

!       xx = 1.d+3 / (R_cav*1.d-8)**3

       do Nz = 1, lz
         do Ny = 1, ly
           do Nx = 1, lx

!             zz = Rhofil(Nx,Ny,Nz) * xx
!             xi = ( dble(Nx) - 0.5d0*lx ) * R_cav
!             yi = ( dble(Ny) - 0.5d0*ly ) * R_cav
!             zi = ( dble(Nz) - 0.5d0*lz ) * R_cav
!             write(1,'(3f7.1,2f5.0,f10.5)') xi, yi, zi, &
!             &           real(Flcav(Nx,Ny,Nz)), real(Flfil(Nx,Ny,Nz)), zz

             write(1,'(2i2)') Flcav(Nx,Ny,Nz),Flfil(Nx,Ny,Nz)

           end do
         end do
       end do

       close(1)

       end if

     end do
   end do

   open(1,file='./Analy/AveCavity.dat',status='unknown')
   open(2,file='./Analy/AveCavityZ.dat',status='unknown')

   lzh = lz / 2

   do Nz = 1 , lz
     zz = 0.d0
     do Ny = 1, ly
       do Nx = 1, lx
         xx = dble( FlcavAcc(Nx,Ny,Nz) ) / dble( TotalStepNumber )
         zz = zz + xx
         write(1,'(f15.10)') xx
       end do
     end do
     zz = zz / dble( lx * ly )
     if(lzh*2==lz) then
       write(2,'(f7.2,f15.10)') (Nz-lzh-0.5)*R_cav, zz
     else
       write(2,'(f7.2,f15.10)') (Nz-lzh-1.0)*R_cav, zz
     end if
   end do

   close(1)
   close(2)

   open(1,file='./Analy/cavity.fld',status='unknown')

   write(1,'(a)') '# AVS field file'
   write(1,'(a)') '#'
   write(1,'(a,i4)') 'nstep= ',NumF
   write(1,'(a)') 'ndim= 3'
   write(1,'(a,i6)') 'dim1= ',lx
   write(1,'(a,i6)') 'dim2= ',ly
   write(1,'(a,i6)') 'dim3= ',lz
   write(1,'(a)') 'nspace= 3'
   write(1,'(a)') 'veclen= 2'
   write(1,'(a)') 'data= integer'
   write(1,'(a)') 'field= uniform'
   write(1,'(a)') 'label= cavity  allfilled'
   write(1,*)

   do i = 1 , NumF

     if(i<10) then
       write(FileName,'(a,i1,a)') 'distCav000',i,'.dat'
       write(1,'(a,i1)') 'time value= ',i
     else if(i<100) then
       write(FileName,'(a,i2,a)') 'distCav00', i,'.dat'
       write(1,'(a,i2)') 'time value= ',i
     else if(i<1000) then
       write(FileName,'(a,i3,a)') 'distCav0',  i,'.dat'
       write(1,'(a,i3)') 'time value= ',i
     else if(i<10000) then
       write(FileName,'(a,i4,a)') 'distCav',   i,'.dat'
       write(1,'(a,i4)') 'time value= ',i
     else
       write(*,*) 'ERROR : too large number of files'
       call Finalize
     end if

   write(1,'(a,a,a)') 'variable 1 file=./',trim(FileName),&
   &                  ' filetype=ascii offset=0 stride=2'
   write(1,'(a,a,a)') 'variable 2 file=./',trim(FileName),&
   &                  ' filetype=ascii offset=1 stride=2'
   write(1,'(a)') 'EOT'
   write(1,'(a)') '#'

   end do

   close(1)

end subroutine CavityDist


!######################################################################
!######################################################################


subroutine RotRing5

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use CommonBlocks, only : QMaster
use ParamAnalyze
use CommonMPI
use CellParam, only : H, InvH
use AtomParam, only : AtomName, MolName, Mass

implicit none

integer :: i , j, k, l, ii, jj, Nw, nm, ita, Nha
integer :: TotalStepNumber, NumF
integer, dimension(:,:), allocatable :: IdRing
real(8) :: cs, cst, cst2, R2, MassM
real(8), dimension(3) :: Rij, Rg
real(8), dimension(:,:,:), allocatable :: Si
real(8), dimension(:), allocatable :: Fcorr, Fcorr2
character(len=72) :: Fname

   ii = 0
   do i = 1 , NumSpec
     if(trim(MolName(i))==trim(MSDMol(1))) then
       Nw = i
       nm = NumMol(i)
       exit
     end if
     ii = ii + NumAtm(i) * NumMol(i)
   end do

   if(cWhat=='Rot5Ring') then
     Nha = 5
   else if(cWhat=='RotVec') then
     Nha = 2
   end if

   allocate( IdRing(Nha,nm) )

   IdRing = 0
   k = ii
   do i = 1, nm
     do j = 1, NumAtm(Nw)
       k = k + 1
inn9:  do l = 1, Nha
         if(AtomName(k)==MSDAtom(l)) then
           IdRing(l,i) = k
           exit inn9
         end if
       end do inn9
     end do
   end do

   MassM = 0.d0
   do i = 1, Nha
     MassM = MassM + Mass(IdRing(i,1))
   end do
   MassM = 1.d0 / MassM

   do i = 1, nm
     do j = 1, Nha
       if(IdRing(j,i)==0) then
         write(*,*) "ERROR: no atom is found for ",MSDAtom(j)
         call Finalize
       end if
     end do
   end do

   TotalStepNumber = 0
   do i = 1 , NJobs
     TotalStepNumber = TotalStepNumber + NTrjStep(i)
   end do

   allocate( Si(3,nm,TotalStepNumber) )

   allocate( Fcorr(TotalStepNumber) )
   allocate( Fcorr2(TotalStepNumber) )

   NumF = 0

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

       call InversMatrix(H,InvH)

       NumF = NumF + 1

       jj = 0

       if(cWhat=='Rot5Ring') then

       do k = 1, nm

         Rg = 0.d0
         do ii = 1 , Nha
           l = IdRing(ii,k)
           Rg = Rg + Mass(l)*R(:,l)
         end do
         Rg = Rg * MassM

         Rij = R(:,IdRing(2,k)) - Rg(:)
         R2  = dot_product( Rij, Rij )
         Si(:,k,NumF) = Rij / sqrt(R2)
       end do

       else if(cWhat=='RotVec') then

       do k = 1, nm
         Rij = R(:,IdRing(2,k)) - R(:,IdRing(1,k))
         R2  = dot_product( Rij, Rij )
         Si(:,k,NumF) = Rij / sqrt(R2)
       end do

       end if

     end do

   end do

   do l = 1 , TotalStepNumber - 1

     ita = TotalStepNumber - l

     cst  = 0.d0
     cst2 = 0.d0

     do k = 1 , ita

       do j = 1 , nm

         cs   = dot_product( Si(:,j,k), Si(:,j,k+l) )
         cst  = cst + cs
         cst2 = cst2 + cs * cs

       end do

     end do

     Fcorr(l)  = cst  / dble(nm*ita)
     Fcorr2(l) = cst2 / dble(nm*ita)

   end do

   if(cWhat=='Rot5Ring') then
     write(Fname,'(a,a,a)') './Analy/Rot_Ring_',trim(MSDMol(1)),'.dat'
   else if(cWhat=='RotVec') then
     write(Fname,'(5a)') './Analy/Rot_Vec_',trim(MSDAtom(1)),'-',trim(MSDAtom(2)),'.dat'
   end if

   open(1,file= trim(Fname))

   write(1,'(a)') '# Rotational correlation function '
   write(1,'(a)') '# Legendre polynomial 1 C1(t) = cos(theta)'
   write(1,'(a)') '# Legendre polynomial 2 C2(t) = 0.5*(3.d0*cos(theta)-1)'
   write(1,'(a)') '# cos(theta) = < u(t)u(0) >'
   write(1,'(a)') '# where u(t) is the dipole unit-vector at t'
   write(1,'(a)') '# Time[ps]    C1(t)    C2(t)'
   write(1,'(5x,a6,2(11x,a7))') '0.0000','1.00000','1.00000'

   do l = 1 , TotalStepNumber - 1
     write(1,'(1x,f10.4,2(1x,f17.5))') dtime*l,Fcorr(l),0.5d0*(3.d0*Fcorr2(l)-1.d0)
   end do

   close(1)

end subroutine RotRing5


!######################################################################
!######################################################################


subroutine Anal_InerM

use ParamInertia
use ParamAnalyze
use Numbers, only : NumSpec, NumMol, NumAtm
use CommonBlocks, only : QMaster

integer :: i, ii
real(8), dimension(3) :: EValue, Rg
real(8), dimension(3,3) :: RMat, Rot

   ii = 0
   do i = 1, NumSpec
     if(i==Kcomp) exit
     ii = ii + NumMol(i)*NumAtm(i)
   end do
   NumInert = NumMol(Kcomp)*NumAtm(Kcomp)
   allocate( ListInert(NumInert) )
   allocate( Msel(NumInert) )
   allocate( Rsel(3,NumInert) )
   do i = 1, NumInert
     ii = ii + 1
     ListInert(i) = ii
   end do

   if(QMaster) open(1,file='Analy/MomentOfInertia.dat')

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

       if(QMaster) call CalcInertia(RMat,Rg)
       print *, RMat(1,:)
       print *, RMat(2,:)
       print *, RMat(3,:)
       if(QMaster) call Jacobi(RMat,Rot,Evalue)

       if(QMaster) write(1,'(3e16.8)') Evalue(:)

     end do

   end do

   close(1)

end subroutine Anal_InerM
