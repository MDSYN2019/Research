! ############################
! ## SUBROUTINE LIST 
! ## -- Force_EAM 
! ## -- Set_EAM 
! ## -- spline 
! ## -- splint 
! ############################


!**********************************************************************
!**********************************************************************


subroutine Force_EAM

use Numbers, only : N
use CommonBlocks, only : QPathInt, ForceField
use Configuration, only : R
use EAM_param
use CommonMPI
use CommonPI
use BookParam, only : Npair, ListIJ
use UnitExParam, only : reng, e
use CellParam, only : H, InvH
use CutoffParam, only : Rcutoff2
use AtomParam , only : ResidNum

implicit none

real(8), dimension(3,N) :: ScR
real(8), dimension(N) :: Phi, F
real(8) :: Rho, FRho, dRhoiRdR, dRhojRdR, PhiR, dPhiRdR
real(8), dimension(3,N) :: Frc_f, Frc_Phi
real(8), dimension(3) :: Sij, Rij
integer, dimension(3) :: Nij
integer :: i, j, l, kk
integer :: Numi, Numj, Num2ij
real(8), dimension(3) :: dummy
real(8) :: R1, R2, Raa, Rinv, SPhi, SF, SRhoi

#ifdef MEAM
real(8), dimension(N,NumEAM) :: SRho, dFdRho
integer :: k, Numii, Numjj
real(8) :: dFdRhoi, dFdRhoj
#else
real(8), dimension(N) :: SRho, dFdRho
#endif

integer :: NProcsTemp

real(8), dimension(3,N,-13:13) :: ForceD
real(8), dimension(3,-13:13) :: Gk

external PhiR

   if(QPathInt) then
     NProcsTemp = NumProcess
   else
     NProcsTemp = NProcs
   end if

!**********************************************************************

   Vir_EAM = 0.d0
   if(ForceField(1:3) /= 'EAM') then
     Frc_EAM = 0.d0
     Return
   end if

   SRho    = 0.d0
   Phi     = 0.d0
   Frc_F   = 0.d0
   Frc_Phi = 0.d0
   F       = 0.d0
   ForceD  = 0.d0

   do i = 1 , N

#ifdef PCC
     ScR(1,i) = InvH(1,1) * R(1,i) + InvH(1,2) * R(2,i) + InvH(1,3) * R(3,i)
     ScR(2,i) = InvH(2,1) * R(1,i) + InvH(2,2) * R(2,i) + InvH(2,3) * R(3,i)
     ScR(3,i) = InvH(3,1) * R(1,i) + InvH(3,2) * R(2,i) + InvH(3,3) * R(3,i)
#else
     ScR(:,i) = matmul( InvH , R(:,i) )
#endif

   end do


!**********************************************************************

   do l = 1 , Npair

     i = ListIJ(1,l)
     j = ListIJ(2,l)

     Sij = ScR(:,i) - ScR(:,j)
#ifdef PCC
     Nij(:) = 0
     if(Sij(1) >  0.5d0) Nij(1) = - 1
     if(Sij(1) < -0.5d0) Nij(1) =   1
     if(Sij(2) >  0.5d0) Nij(2) = - 1
     if(Sij(2) < -0.5d0) Nij(2) =   1
     if(Sij(3) >  0.5d0) Nij(3) = - 1
     if(Sij(3) < -0.5d0) Nij(3) =   1
     Sij(:) = Sij(:) + Nij(:)
     Rij(1) = H(1,1) * Sij(1) + H(1,2) * Sij(2) + H(1,3) * Sij(3)
     Rij(2) = H(2,1) * Sij(1) + H(2,2) * Sij(2) + H(2,3) * Sij(3)
     Rij(3) = H(3,1) * Sij(1) + H(3,2) * Sij(2) + H(3,3) * Sij(3)
     R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
     Nij = - nint( Sij )
     Sij = Sij + Nij
     Rij = matmul( H , Sij )
     R2  = dot_product( Rij, Rij )
#endif

     if(R2 < Rcutoff2) then

       R1 = sqrt(R2)

       Raa = R1

#ifdef MEAM

       Numii = ResidNum(i)
       Numjj = ResidNum(j)

       do k = 1, NumEAM

         Numi = (Numii - 1)*NumEAM + k
         Numj = (Numjj - 1)*NumEAM + k

         SRho(i,k) = SRho(i,k) + Rho(Raa,0,Numj)
         SRho(j,k) = SRho(j,k) + Rho(Raa,0,Numi)

       end do

#else

       Numi = ResidNum(i)
       Numj = ResidNum(j)

       SRho(i) = SRho(i) + Rho(Raa,0,Numj)
       SRho(j) = SRho(j) + Rho(Raa,0,Numi)

#endif

     end if

   end do


   call SumRho(SRho)

!**********************************************************************

   do i = 1, N

#ifdef MEAM
     do k = 1, NumEAM
       SRhoi = SRho(i,k)
       Numi  = (ResidNum(i)-1)*NumEAM + k
       F(i)  = F(i) + FRho(SRhoi,0,Numi) *e*reng
       dFdRho(i,k) = FRho(SRhoi,1,Numi) *e*reng
     end do
#else
     SRhoi = SRho(i)
     Numi  = ResidNum(i)
     F(i)  = FRho(SRhoi,0,Numi) *e*reng
     dFdRho(i) = FRho(SRhoi,1,Numi) *e*reng
#endif

   end do

!**********************************************************************

   do l = 1 , Npair

     i = ListIJ(1,l)
     j = ListIJ(2,l)

     Sij = ScR(:,i) - ScR(:,j)
#ifdef PCC
     Nij(:) = 0
     if(Sij(1) >  0.5d0) Nij(1) = - 1
     if(Sij(1) < -0.5d0) Nij(1) =   1
     if(Sij(2) >  0.5d0) Nij(2) = - 1
     if(Sij(2) < -0.5d0) Nij(2) =   1
     if(Sij(3) >  0.5d0) Nij(3) = - 1
     if(Sij(3) < -0.5d0) Nij(3) =   1
     Sij(:) = Sij(:) + Nij(:)
     Rij(1) = H(1,1) * Sij(1) + H(1,2) * Sij(2) + H(1,3) * Sij(3)
     Rij(2) = H(2,1) * Sij(1) + H(2,2) * Sij(2) + H(2,3) * Sij(3)
     Rij(3) = H(3,1) * Sij(1) + H(3,2) * Sij(2) + H(3,3) * Sij(3)
     R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
     Nij = - nint( Sij )
     Sij = Sij + Nij
     Rij = matmul( H , Sij )
     R2  = dot_product( Rij, Rij )
#endif

     R1 = sqrt(R2)

     Rinv = 1.d0/R1

     Raa = R1

     kk = 9 * Nij(1) + 3 * Nij(2) + Nij(3)

#ifdef MEAM

     do k = 1, NumEAM
        Numi = (ResidNum(i)-1)*NumEAM + k
        Numj = (ResidNum(j)-1)*NumEAM + k

        dRhoiRdR = Rho(Raa,1,Numi)
        dRhojRdR = Rho(Raa,1,Numj)

        dFdRhoi = dFdRho(i,k)
        dFdRhoj = dFdRho(j,k)

        dummy(:) =    + dFdRhoi * dRhojRdR * Rij(:) * Rinv &
        &             + dFdRhoj * dRhoiRdR * Rij(:) * Rinv

        ForceD(:,i,kk) = ForceD(:,i,kk) - dummy(:)
        Frc_F(:,j)     = Frc_F(:,j)     + dummy(:)

     end do

     Num2ij = ResidPair( ResidNum(i), ResidNum(j) )

     dPhiRdR = PhiR(Raa,1,Num2ij) *1.6021892d-23

     dummy(1) = PhiR(Raa,0,Num2ij) *1.6021892d-23

     Phi(i) = Phi(i) + dummy(1)
     Phi(j) = Phi(j) + dummy(1)

     dummy(:) = dPhiRdR * Rij(:) * Rinv

     ForceD(:,i,kk) = ForceD(:,i,kk) - dummy(:)
     Frc_Phi(:,j)   = Frc_Phi(:,j)   + dummy(:)

#else

     Numi = ResidNum(i)
     Numj = ResidNum(j)

     dRhoiRdR = Rho(Raa,1,Numi)
     dRhojRdR = Rho(Raa,1,Numj)

     dummy(:) = ( dFdRho(i) * dRhojRdR &
     &          + dFdRho(j) * dRhoiRdR ) * Rij(:) * Rinv

     ForceD(:,i,kk) = ForceD(:,i,kk) - dummy(:)
     Frc_F(:,j)     = Frc_F(:,j)     + dummy(:)

     Num2ij = ResidPair( ResidNum(i), ResidNum(j) )

     dPhiRdR = PhiR(Raa,1,Num2ij) *1.6021892d-23

     dummy(1) = PhiR(Raa,0,Num2ij) *1.6021892d-23

     Phi(i) = Phi(i) + dummy(1)
     Phi(j) = Phi(j) + dummy(1)

     dummy(:) = dPhiRdR * Rij(:) * Rinv

     ForceD(:,i,kk) = ForceD(:,i,kk) - dummy
     Frc_Phi(:,j)   = Frc_Phi(:,j)   + dummy

#endif

   end do

!**********************************************************************

   Phi(:) = Phi(:) * 0.5d0

   SPhi = sum( Phi(:) )
   SF   = sum( F(:) ) / dble( NProcsTemp )
   Ene_EAM = SPhi + SF

   Frc_EAM(:,:) = Frc_F(:,:) + Frc_Phi(:,:)

   Gk = 0.d0

   do kk = -13, 13
     do i = 1, N
       Gk(:,kk)     = Gk(:,kk)     + ForceD(:,i,kk)
       Frc_EAM(:,i) = Frc_EAM(:,i) + ForceD(:,i,kk)
     end do
   end do

   if(QPathInt) then
     call VirialBekkerPI(Frc_EAM,Vir_EAM,Gk)
   else
     call VirialBekker(Frc_EAM,Vir_EAM,Gk)
   end if


end subroutine Force_EAM


!**********************************************************************
!**********************************************************************


subroutine Set_EAM

use Numbers, only : N, NumSpec
use EAM_param
use UnitExParam, only : Avogadro
use OptConstraintParam, only : NumOptC, NumPLC
use BondedParam, only : NumBond, NumAngle, NumUB, NumDihedral, NumImproper
use AtomParam, only : MolName, AtomName, Mass, InvMass

implicit none

character(len=80) :: namf
integer :: i, j, ijSpec, iSpec, jSpec
real(8) :: R1, FRhoPlus, FRhoMinus, PhiPlus, PhiMinus, RhoPlus, RhoMinus
real(8) :: FRho, Rho, PhiR
external FRho, PhiR, Rho
#ifdef MEAM
integer :: k
character(len=1) :: alphabet
external alphabet
#endif

   NumBond = 0
   NumAngle = 0
   NumUB = 0
   NumDihedral = 0
   NumImproper = 0
   NumOptC = 0
   NumPLC = 0

   ijSpec = 0
   do i = 1, NumSpec
     do j = i, NumSpec
       ijSpec = ijSpec + 1
       ResidPair(i,j) = ijSpec
       ResidPair(j,i) = ijSpec
     end do
   end do

   allocate( Mass(N) )
   allocate( InvMass(N) )

   do i = 1 , N
     if(AtomName(i)(1:2) == 'Ni') then
       Mass(i) = 58.6934d-3/Avogadro
     else if(AtomName(i)(1:1) == 'H') then
       Mass(i) = 1.00794d-3/Avogadro
     else if(AtomName(i)(1:1) == 'Fe') then
       Mass(i) = 55.854d-3/Avogadro
     else if(AtomName(i)(1:1) == 'Au') then
       Mass(i) = 196.96655d-3/Avogadro
     else if(AtomName(i)(1:1) == 'Cu') then
       Mass(i) = 63.546d-3/Avogadro
     else if(AtomName(i)(1:1) == 'Ag') then
       Mass(i) = 107.8682d-3/Avogadro
     else if(AtomName(i)(1:1) == 'Pt') then
       Mass(i) = 195.08d-3/Avogadro
     else if(AtomName(i)(1:1) == 'Pd') then
       Mass(i) = 106.42d-3/Avogadro
     else
       write(*,*) 'ERROR : EAM is valid only for Ni, Fe, Au, Cu. Ag, Pt, Pd, and H'
       call Finalize
     end if
     InvMass(i) = 1.d0 / Mass(i)
   end do

#ifdef MEAM

   do jSpec = 1, NumSpec

     do k = 1, NumEAM

     iSpec = (jSpec - 1) * NumEAM + k

     write(namf,'(a,a,a)') trim(adjustl(MolName(jSpec))),'.1',alphabet(k)
     open (15,file=trim(namf),status='unknown' )

#else

   do iSpec = 1, NumSpec

     write(namf,'(a,a)') trim(adjustl(MolName(iSpec))),'.1'
     open (15, file = trim(namf),status = 'unknown' )

#endif

       read (15,'(/)')

       do i = 1, Nref
         read (15,*) Rrefa(i,iSpec), Rhorefa(i,iSpec)
       end do

     close (15)

     Rmina(iSpec) = Rrefa(1,iSpec)
     Rmaxa(iSpec) = Rrefa(Nref,iSpec)

#ifdef MEAM
     call spline ( Rrefa, Rhorefa, Rhoy2, Nref, iSpec, NumEAMSpec )
#else
     call spline ( Rrefa, Rhorefa, Rhoy2, Nref, iSpec, NumSpec )
#endif

     do i = 1, Nref

       R1 = Rrefa(i,iSpec) + dh

#ifdef MEAM
       call splint( Rrefa, Rhorefa, Rhoy2, Nref, R1, RhoPlus, iSpec, NumEAMSpec )
#else
       call splint( Rrefa, Rhorefa, Rhoy2, Nref, R1, RhoPlus, iSpec, NumSpec )
#endif

       R1 = Rrefa(i,iSpec) - dh

#ifdef MEAM
       call splint( Rrefa, Rhorefa, Rhoy2, Nref, R1, RhoMinus, iSpec, NumEAMSpec )
#else
       call splint( Rrefa, Rhorefa, Rhoy2, Nref, R1, RhoMinus, iSpec, NumSpec )
#endif

       dRhodRref(i,iSpec) = (RhoPlus-RhoMinus)/(2.d0*dh)

     end do

#ifdef MEAM
     call spline ( Rrefa, dRhodRref, dRhodRy2, Nref,iSpec, NumEAMSpec )

     end do
#else
     call spline ( Rrefa, dRhodRref, dRhodRy2, Nref,iSpec, NumSpec )
#endif

   end do

#ifdef MEAM
   do jSpec = 1, NumSpec

     do k = 1, NumEAM

      iSpec = (jSpec-1)*NumEAM + k

     write(namf,'(a,a,a)') trim(adjustl(MolName(jSpec))),'.2',alphabet(k)
     open (15, file = trim(namf), status = 'unknown' )
#else
   do iSpec = 1, NumSpec

     write(namf,'(a,a)') trim(adjustl(MolName(iSpec))),'.2'
     open (15, file = trim(namf),status = 'unknown' )
#endif

       read (15,'(/)')

       do i = 1, Nref

         read (15,*) Rhorefb(i,iSpec), FRhoref(i,iSpec)

       end do

     close (15)

     Rhomin(iSpec) = Rhorefb(1,iSpec)
     Rhomax(iSpec) = Rhorefb(Nref,iSpec)

#ifdef MEAM
     call spline ( Rhorefb, FRhoref, FRhoy2, Nref,iSpec, NumEAMSpec )
#else
     call spline ( Rhorefb, FRhoref, FRhoy2, Nref,iSpec, NumSpec )
#endif

     do i = 1, Nref

       RhoPlus = Rhorefb(i,iSpec) + dh

#ifdef MEAM
       call splint( Rhorefb, FRhoref, FRhoy2, Nref, RhoPlus, FRhoPlus,iSpec, NumEAMSpec )
#else
       call splint( Rhorefb, FRhoref, FRhoy2, Nref, RhoPlus, FRhoPlus,iSpec, NumSpec )
#endif

       RhoMinus = Rhorefb(i,iSpec) - dh

#ifdef MEAM
       call splint( Rhorefb, FRhoref, FRhoy2, Nref, RhoMinus, FRhoMinus,iSpec, NumEAMSpec )
#else
       call splint( Rhorefb, FRhoref, FRhoy2, Nref, RhoMinus, FRhoMinus,iSpec, NumSpec )
#endif

       dFdRhoref(i,iSpec) = (FRhoPlus-FRhoMinus)/(2.d0*dh)

     end do

#ifdef MEAM
     call spline ( Rhorefb, dFdRhoref, dFdRhoy2, Nref,iSpec, NumEAMSpec )
     end do
#else
     call spline ( Rhorefb, dFdRhoref, dFdRhoy2, Nref,iSpec, NumSpec )
#endif

   end do

   ijSpec = 0

   do iSpec = 1, NumSpec
     do jSpec = iSpec, NumSpec

       ijSpec = ijSpec + 1

       write(namf,'(a,a,a,a)') trim(adjustl(MolName(iSpec))),'-',  &
       &  trim(adjustl(MolName(jSpec))),'.3'

       open(15,file=trim(namf),status='unknown')

         read (15,'(/)')

         do i = 1, Nref
           read (15,*) Rrefb(i,ijSpec), PhiRef(i,ijSpec)
         end do

       close (15)

       Rminb(ijSpec) = Rrefb(1,ijSpec)
       Rmaxb(ijSpec) = Rrefb(Nref,ijSpec)

       call spline ( Rrefb, PhiRef, Phiy2, Nref,ijSpec, NumSpecPair )

       do i = 1, Nref

         R1 = Rrefb(i,ijSpec) + dh

         call splint( Rrefb, PhiRef, Phiy2, Nref, R1, PhiPlus,ijSpec, NumSpecPair )

         R1 = Rrefb(i,ijSpec) - dh

         call splint( Rrefb, PhiRef, Phiy2, Nref, R1, PhiMinus,ijSpec, NumSpecPair )

         dPhidRref(i,ijSpec) = (PhiPlus-PhiMinus)/(2.d0*dh)

       end do

       call spline ( Rrefb, dPhidRref, dPhidRy2, Nref,ijSpec, NumSpecPair )

     end do

   end do

end subroutine Set_EAM


!######################################################################
!######################################################################


Function Rho ( R1, iop, iSpec )

use Numbers, only : NumSpec
use EAM_param

implicit none

integer :: iop, iSpec
real(8) :: R1, Rho

   if ( iop == 0 ) then

     if ( ( R1 < Rmina(iSpec) ) .or. ( R1 > Rmaxa(iSpec) ) ) then
       Rho = 0.d0
       return
     end if

#ifdef MEAM
     call splint( Rrefa, Rhorefa, Rhoy2, Nref, R1, Rho, iSpec, NumEAMSpec )
#else
     call splint( Rrefa, Rhorefa, Rhoy2, Nref, R1, Rho, iSpec, NumSpec )
#endif

   else if ( iop == 1 ) then

     if ( ( R1 < Rmina(iSpec) ) .or. ( R1 > Rmaxa(iSpec) ) ) then
       Rho = 0.d0
       return
     end if

#ifdef MEAM
     call splint( Rrefa, dRhodRref, dRhodRy2, Nref, R1, Rho, iSpec, NumEAMSpec )
#else
     call splint( Rrefa, dRhodRref, dRhodRy2, Nref, R1, Rho, iSpec, NumSpec )
#endif

   end if

end Function Rho


!######################################################################
!######################################################################


Function FRho ( Rho, iop, iSpec )

use Numbers, only : NumSpec
use EAM_param

implicit none

integer :: iop, iSpec
real(8) :: Rho, FRho

   if ( iop == 0 ) then

     if ( ( Rho < Rhomin(iSpec) ) .or. ( Rho > Rhomax(iSpec) ) ) then
       FRho = 0.d0
       return
     end if

#ifdef MEAM
     call splint( Rhorefb, FRhoref, FRhoy2, Nref, Rho, FRho, iSpec, NumEAMSpec )
#else
     call splint( Rhorefb, FRhoref, FRhoy2, Nref, Rho, FRho, iSpec, NumSpec )
#endif

   else if ( iop == 1 ) then

     if ( ( Rho < Rhomin(iSpec) ) .or. ( Rho > Rhomax(iSpec) ) ) then
       FRho = 0.d0
       return
     end if

#ifdef MEAM
     call splint( Rhorefb, dFdRhoref, dFdRhoy2, Nref, Rho, FRho, iSpec, NumEAMSpec )
#else
     call splint( Rhorefb, dFdRhoref, dFdRhoy2, Nref, Rho, FRho, iSpec, NumSpec )
#endif

   end if

end Function FRho


!######################################################################
!######################################################################


Function PhiR ( R1, iop, iSpecPair )

use EAM_param

implicit none

integer :: iop, iSpecPair
real(8) :: R1, PhiR

  if ( iop == 0 ) then

    if ( ( R1 < Rminb(iSpecPair) ) .or. &
    &    ( R1 > Rmaxb(iSpecPair) ) ) then
      PhiR = 0.d0
      return
    end if

    call splint( Rrefb, Phiref, Phiy2, Nref, R1, PhiR, iSpecPair, NumSpecPair )

  else if ( iop == 1 ) then

    if ( ( R1 < Rminb(iSpecPair) ) .or. &
    &    ( R1 > Rmaxb(iSpecPair) ) ) then
      PhiR = 0.d0
      return
    end if

    call splint( Rrefb, dPhidRref, dPhidRy2, Nref, R1, PhiR, iSpecPair, NumSpecPair )

  end if

end Function PhiR


!######################################################################
!######################################################################


subroutine spline ( x, y, y2, n, idim, ndim )

implicit none

integer :: n, idim, i, ndim
real(8) :: yp1, ypn, sig, p, qn, un
real(8), dimension(n,ndim) :: x, y, y2
real(8), dimension(n,ndim) :: u

   yp1 = (y(2,idim)-y(1,idim)) / (x(2,idim)-x(1,idim))
   ypn = (y(n,idim)-y(n-1,idim)) / (x(n,idim)-x(n-1,idim))

   y2(1,idim) = -0.5d0
   u(1,idim) = (3.d0 / (x(2,idim)-x(1,idim)) )*((y(2,idim)-y(1,idim)) &
   &           / (x(2,idim)-x(1,idim)) - yp1)

   do i = 2, n-1
     sig = (x(i,idim)-x(i-1,idim)) / (x(i+1,idim)-x(i-1,idim))
     p = sig*y2(i-1,idim)+2.d0
     y2(i,idim) = (sig-1.d0) / p
     u(i,idim) = (6.d0*((y(i+1,idim)-y(i,idim))  &
     &              / (x(i+1,idim)-x(i,idim)) - (y(i,idim)-y(i-1,idim))  &
     &              / (x(i,idim)-x(i-1,idim))) / (x(i+1,idim)-x(i-1,idim))  &
     &              -sig*u(i-1,idim)) / p
   end do

   qn = 0.5d0
   un = (3.d0/(x(n,idim)-x(n-1,idim)))*(ypn-(y(n,idim)-y(n-1,idim)) &
   &    / (x(n,idim)-x(n-1,idim)))

   y2(n,idim) = (un-qn*u(n-1,idim)) / (qn*y2(n-1,idim)+1.d0)

   do i = n-1, 1, -1
     y2(i,idim) = y2(i,idim)*y2(i+1,idim) + u(i,idim)
   end do

end subroutine spline


!######################################################################
!######################################################################


subroutine splint( xa, ya, y2a, n, x, y, idim, ndim )

implicit none

integer :: klo, khi, k, n, idim, ndim
real(8) :: x, y, h, a, b
real(8), dimension(n,ndim) :: xa, ya, y2a

   klo = 1
   khi = n

   do

     if ( khi-klo <= 1 ) exit

     k = (khi+klo)/2

     if ( xa(k,idim) .gt. x ) then
       khi = k
     else
       klo = k
     end if

   end do

   h = xa(khi,idim) - xa(klo,idim)

   if ( h == 0.d0 ) then
     write(*,*) 'bad xa input in splint'
     call Finalize
   end if

   a = (xa(khi,idim)-x)/h
   b = (x-xa(klo,idim))/h

   y = a*ya(klo,idim) + b*ya(khi,idim) &
   &  + ((a*a*a-a)*y2a(klo,idim) &
   &    +(b*b*b-b)*y2a(khi,idim))*(h*h)/6.d0

end subroutine splint


!######################################################################
!######################################################################


#ifdef MEAM

function alphabet(k)

character(len=1) :: alphabet
integer :: k
character(len=1), dimension(26), parameter :: Ch = &
&      (/'a','b','c','d','e','f','g','h','i','j','k','l','m','n',&
&        'o','p','q','r','s','t','u','v','w','x','y','z'/)

   if(k<=0.or.k>=27) then
     alphabet = '0'
   else
     alphabet = Ch(k)
   end if

end function alphabet

#endif
