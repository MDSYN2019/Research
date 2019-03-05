! ############################
! ## SUBROUTINE LIST 
! ## -- ElecPot 
! ############################


!######################################################################
!######################################################################


! *************************************************************
! **  Electric potential, polarization, charge distribution  **
! **  along the z axis                                       **
! **  heterogeneous system such as lipid bilayers            **
! *************************************************************

subroutine ElecPot

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use UnitExParam, only : eps_o, ec, e
use NonbondParam, only : Charge
use CellParam, only : H, InvH

implicit none

integer :: i, j, k, l, ii
real(8), dimension(3) :: a, b
real(8) :: Invrl2
integer, parameter :: ndz = 400
real(8), parameter :: drz = 0.1  ! [A]
real(8), dimension(-ndz:ndz) :: RhoT, PolT, PhiT
real(8), dimension(-ndz:ndz,NumSpec) :: RhoI, PolI, PhiI
real(8), dimension(-ndz:ndz) ::  RhoT_sym , PolT_sym, PhiT_sym
real(8), dimension(-ndz:ndz,NumSpec) :: RhoI_sym, PolI_sym,PhiI_sym
integer :: TotalStep
real(8) :: Area, zz, dz
real(8) :: PT0, PT0_sym
real(8), dimension(NumSpec) :: PI0, PI0_sym

integer, dimension(NumSpec) :: Numa, Numb
integer :: Mesh, iz, idum

real(8), dimension(3,N) :: Rtr

real(8) :: ek

integer, dimension(3) :: Ic
real(8), dimension(3) :: Si
real(8), dimension(3) :: Tc

   Invrl2 = 1.d-20
   ek = 1.d0 / eps_o

   Numa(1) = 0
   Numb(1) = NumMol(1)*NumAtm(1)

   if( NumSpec > 1 ) then

     do i = 2 , NumSpec

       Numa(i) = Numb(i-1)
       Numb(i) = Numb(i-1) + NumMol(i) * NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   Mesh  = nint( 1.d0 / drz )

   Charge = Charge * e / sqrt(ec)

   TotalStep = 0

   RhoT  = 0.d0
   RhoI  = 0.d0

   PolT  = 0.d0
   PolI  = 0.d0

   PhiT  = 0.d0
   PhiI  = 0.d0

   do i = 1 , NJobs

     call OpenTraj(i)

     TotalStep = TotalStep + NTrjStep(i)

     do j = 1 , NTrjStep(i)

!     -----------------------------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -----------------------------------------

       a  = H(:,1)
       b  = H(:,2)

       call CArea(a,b,Area)
       Area = Area * Invrl2

       call CellTransform

       Rtr = R

       do k = 1, N

         Si = matmul( InvH, Rtr(:,k) )
         Ic = -nint(Si)
         Tc = matmul( H, Ic )
         Rtr(:,k) = Rtr(:,k) + Tc

       end do

       do k = 1 , NumSpec

         do l = Numa(k), Numb(k)

           iz = nint( Rtr(3,l) * Mesh )
           if(abs(iz)>400) cycle
           RhoI(iz,k) = RhoI(iz,k) + Charge(l) / Area

         end do

       end do

     end do

   end do

   dz = drz * 1.d-10

   RhoI = RhoI / (TotalStep * dz)

   do ii = -ndz, ndz
     idum = (-1*ii)
     RhoI_sym(ii,:) = 0.5*(RhoI(ii,:)+RhoI(idum,:))
   end do

   do i = 1 , NumSpec

     RhoT_sym(:) = RhoT_sym(:) + RhoI_sym(:,i)
     RhoT(:) = RhoT(:) + RhoI(:,i)

   end do

   do i = 1 , ndz
     do j = 0 , i

       PolT_sym (i)   = PolT_sym (i)   + RhoT_sym (j)   * dz
       PolI_sym (i,:) = PolI_sym (i,:) + RhoI_sym(j,:) * dz

       PolI (i,:) = PolI (i,:) + RhoI (j,:) * dz
       PolT (i)   = PolT (i)   + RhoT (j)   * dz

     end do
   end do
   
   do i = 0, -ndz , -1
     do j = 0 , i, -1
     
       PolT_sym  (i)   = PolT_sym  (i)   + RhoT_sym (j)   * dz
       PolI_sym  (i,:) = PolI_sym  (i,:) + RhoI_sym(j,:) * dz
       PolI (i,:) = PolI (i,:) + RhoI (j,:) * dz
       PolT (i)   = PolT (i)   + RhoT (j)   * dz

     end do
   end do

   do i = 1 , ndz
     do j = 0 , i

       PhiT_sym (i)   = PhiT_sym (i)   - PolT_sym (j)   * dz
       PhiI_sym (i,:) = PhiI_sym (i,:) - PolI_sym (j,:) * dz
       PhiT (i)   = PhiT (i)   - PolT (j)   * dz
       PhiI (i,:) = PhiI (i,:) - PolI (j,:) * dz

     end do
   end do

   do i = 0, -ndz ,  -1
     do j = 0 , i, -1

       PhiT_sym (i)   = PhiT_sym (i)   - PolT_sym (j)   * dz
       PhiI_sym (i,:) = PhiI_sym (i,:) - PolI_sym (j,:) * dz
       PhiT (i)   = PhiT (i)   - PolT (j)   * dz
       PhiI (i,:) = PhiI (i,:) - PolI (j,:) * dz

     end do
   end do

   PT0_sym     = PhiT_sym  (0)
   PI0_sym (:) = PhiI_sym  (0,:)
   PT0    = PhiT (0)
   PI0(:) = PhiI (0,:)

   do i = -ndz , ndz

     PhiT_sym (i)   = ( PhiT_sym (i)   - PT0_sym     ) * ek
     PhiI_sym (i,:) = ( PhiI_sym (i,:) - PI0_sym (:) ) * ek
     PhiT(i)   = ( PhiT(i)   - PT0    ) * ek
     PhiI(i,:) = ( PhiI(i,:) - PI0(:) ) * ek

   end do

   open(81,file='./Analy/ElecRhoS_CLp.dat',status='unknown')
   open(82,file='./Analy/ElecPolS_CLp.dat',status='unknown')
   open(83,file='./Analy/ElecPhiS_CLp.dat',status='unknown')
 
   write(81,'(a)') "# i, zz (A),     RhoT,     RhoT_sym,    RhoI,       RhoI_sym"
   write(82,'(a)') "# i, zz (A),     PolT,     PolT_sym,    PolI,       PolI_sym"
   write(83,'(a)') "# i, zz (A),     PhiT,     PhiT_sym,    PhiI,       PhiI_sym"

   do i = -ndz , ndz

     zz = dble(i) * drz

     write(82,'(i5,f7.2,10(e15.5))') i,zz, PolT(i), PolT_sym(i), PolI(i,:), PolI_sym(i,:)
     write(83,'(i5,f7.2,10(e15.5))') i,zz, PhiT(i), PhiT_sym(i), PhiI(i,:), PhiI_sym(i,:)
     write(81,'(i5,f7.2,10(e15.5))') i,zz, RhoT(i), RhoT_sym(i), RhoI(i,:), RhoI_sym(i,:)

   end do

   close(81)
   close(82)
   close(83)

end subroutine ElecPot
