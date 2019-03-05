


!######################################################################
!######################################################################


subroutine NMAcalc

use Numbers, only : N
use CommonBlocks, only : QPBC
use Configuration, only : R
use EAM_param, only : Frc_EAM
use HeatCap
use EwaldParam, only : Frc_Eksp
use OptConstraintParam, only : Frc_OptC
use NonbondParam, only : Frc_Ersp, Frc_Elec
use BondedParam, only : Frc_Bond, Frc_Angle, Frc_UB, Frc_Dihed, Frc_Impro
use AtomParam, only : Mass

implicit none

integer :: i, j, ii, jj, k, l, Nall
real(8), dimension(3,N) :: FP, FM
real(8), dimension(:,:), allocatable :: Evector
real(8), dimension(:), allocatable :: Evalue
real(8), parameter :: deltah = 1.d-5 ! tunable parameter in A 
real(8) :: xx, yy, zz, rdh2

   Nall = 3 * N
   rdh2 = 1.d0 / (deltah * 2.d0)

   allocate( Hessian(Nall,Nall) )
   allocate( Evector(Nall,Nall) )
   allocate( Evalue(Nall) )

! ## Hessian

   do i = 1, N
     xx = R(1,i)
     R(1,i) = xx + deltah
     call Fcalc(FP)
     R(1,i) = xx - deltah
     call Fcalc(FM)
     R(1,i) = xx
     ii = (i-1) * 3 + 1
     do j = 1, N
       jj = (j-1) * 3
       xx = - ( FP(1,j) - FM(1,j) ) * rdh2
       Hessian(ii,jj+1) = xx
       yy = - ( FP(2,j) - FM(2,j) ) * rdh2
       Hessian(ii,jj+2) = yy
       zz = - ( FP(3,j) - FM(3,j) ) * rdh2
       Hessian(ii,jj+3) = zz
     end do
   end do

   do i = 1, N
     yy = R(2,i)
     R(2,i) = yy + deltah
     call Fcalc(FP)
     R(2,i) = yy - deltah
     call Fcalc(FM)
     R(2,i) = yy
     ii = (i-1) * 3 + 2
     do j = 1, N
       jj = (j-1) * 3
       xx = - ( FP(1,j) - FM(1,j) ) * rdh2
       Hessian(ii,jj+1) = xx
       yy = - ( FP(2,j) - FM(2,j) ) * rdh2
       Hessian(ii,jj+2) = yy
       zz = - ( FP(3,j) - FM(3,j) ) * rdh2
       Hessian(ii,jj+3) = zz
     end do
   end do

   do i = 1, N
     zz = R(3,i)
     R(3,i) = zz + deltah
     call Fcalc(FP)
     R(3,i) = zz - deltah
     call Fcalc(FM)
     R(3,i) = zz
     ii = (i-1) * 3 + 3
     do j = 1, N
       jj = (j-1) * 3
       xx = - ( FP(1,j) - FM(1,j) ) * rdh2
       Hessian(ii,jj+1) = xx
       yy = - ( FP(2,j) - FM(2,j) ) * rdh2
       Hessian(ii,jj+2) = yy
       zz = - ( FP(3,j) - FM(3,j) ) * rdh2
       Hessian(ii,jj+3) = zz
     end do
   end do

! ##

   do ii = 1, 3
   do jj = 1, 3
     do i = 1, N
     do j = 1, N
       k = i + (ii-1)*N
       l = j + (jj-1)*N
       Hessian(k,l) = Hessian(k,l)/sqrt(Mass(i)*Mass(j))*1.d+24  ! [ 1 / s^2 ]
     end do
     end do
   end do
   end do

! ## 

   call diag (Hessian, Evalue, Evector, Nall, Nall)

   call Print_NormalMode(Evalue,Nall)

Contains

   subroutine Fcalc(FF)

   integer :: kk
   real(8), dimension(3,N) :: FF

! -------------------------------
     if(QPBC) then
       call GetForce(0,0)
     else
       call GetForceIso(0)
     end if
! -------------------------------

     if(QPBC) then

       do kk = 1 , N

         FF(:,kk) = ( Frc_Bond (:,kk) &
         &          + Frc_Angle(:,kk) &
         &          + Frc_UB   (:,kk) &
         &          + Frc_Dihed(:,kk) &
         &          + Frc_Impro(:,kk) &
         &          + Frc_OptC (:,kk) &
         &          + Frc_EAM  (:,kk) &
         &          + Frc_Ersp (:,kk) &
         &          + Frc_Eksp (:,kk) )

       end do

     else

       do kk = 1, N
         FF(:,kk) = ( Frc_Bond (:,kk) &
         &          + Frc_Angle(:,kk) &
         &          + Frc_UB   (:,kk) &
         &          + Frc_OptC (:,kk) &
         &          + Frc_Dihed(:,kk) &
         &          + Frc_Impro(:,kk) &
         &          + Frc_Elec (:,kk) )
       end do

     end if

     call SumFrc( FF )

   end subroutine Fcalc

end subroutine NMAcalc


!######################################################################
!######################################################################


subroutine diag ( a, e, c, n, lda )

implicit real(8) (a-h,o-z)

integer :: n, lwork
real(8), dimension(n,n) :: a, c
real(8), dimension(n) :: e
real(8), dimension((n+2)*n) :: work
character(len=1) :: jobz, uplo

   lwork = (n+2) * n

   jobz = 'V'
   uplo = 'U'

   c(:,:) = a(:,:)

#ifdef LAPACK
   call DSYEV( jobz, uplo, n, c, lda, e, work, lwork, info )
#endif

end subroutine diag


!######################################################################
!######################################################################


subroutine Print_NormalMode(EigenVector, N)

implicit none

integer :: i, N
real(8), dimension(N) :: EigenVector
real(8), dimension(N) :: Hz, cm
real(8), parameter :: c_light = 2.99792458d+10 ! [cm / s]
real(8) :: pi2, xx

   pi2 = 1.d0 / 8.d0 / atan(1.d0)
   xx = pi2/c_light

   open(51,file='NormalMode_Hz.dat',status='unknown')
   open(52,file='NormalMode_cm.dat',status='unknown')

   do i = 1, N
     if( EigenVector(i) >= 0.) then
       Hz(i) =   sqrt( EigenVector(i))
     else
       Hz(i) = - sqrt(-EigenVector(i))
     end if
   end do

   do i = 1, N
     cm(i) = Hz(i) * xx
   end do

   do i = 1, N
     write(51,'(i4,d23.16)') i, Hz(i)
     write(52,'(i4,d23.16)') i, cm(i)
   end do

   close(51)
   close(52)

end subroutine Print_NormalMode
