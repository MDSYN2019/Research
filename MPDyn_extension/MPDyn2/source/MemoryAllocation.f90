! ############################
! ## SUBROUTINE LIST 
! ## -- MemoryAllocation 
! ############################


!######################################################################
!######################################################################


subroutine MemoryAllocation

use Numbers, only : NfT, NumSpec, NumMol, NumAtm, NumMer, Number
use CommonBlocks, only : QPBC, PolyFlag, ForceField, QSHAKE, QPathInt
use Configuration
use FFParameters
use EAM_param
use CommonPI
use HeatCap
use EwaldParam, only : Frc_Eksp, Vir_Eksp, Ene_Eksp
use OptConstraintParam, only : Frc_OptC, Vir_OptC, Ene_OptC, HamR
use SHAKEparam, only : Frc_Const
use BookParam
use NonbondParam, only : Frc_Ersp, Frc_Elec, Frc_NBlong, Frc_NBshrt, &
&   Vir_Ersp, Vir_NBlong, Vir_NBshrt, Ene_Elec, Ene_Ersp,  &
&   Ene_NBlong, Ene_NBshrt
use SMDparam, only : Vir_ConstR, Vir_ConstV
use BondedParam, only : &
&   Frc_Bond, Frc_Angle, Frc_UB, Frc_Dihed, Frc_Impro, &
&   Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro, &
&   Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use AtomParam, only : MolName

implicit NONE

integer :: i, j, k, NumAdd, maxpsh
character(len=3) :: String
character(len=4) :: RName

   do i = 1 , NumSpec

     if(PolyFlag(i)) then

       if( ForceField(1:5) == 'CHARM' ) then
         write(*,*) 'ERROR : CHARMm for polymers has not been supported yet!'
         call Finalize
       end if

       String = MolName(i)(5:7)

       NumAtm(i) = 0

       do j = 1 , NumMer(i)

         if( j == 1 ) then
           write(RName,'(a,a)') trim(adjustl(String)),'I'
         else if( j == NumMer(i) ) then
           write(RName,'(a,a)') trim(adjustl(String)),'F'
         else
           write(RName,'(a)') trim(adjustl(String))
         end if

! ##################

         do k = 1 , NumResidueParam
           if(ResiNameParam(k) == RName) then
             NumAdd = NumAtom_inResi(k)
             exit
           end if
           if(k == NumResidueParam) then
             write(*,*) 'ERROR: polymer resid.'
             write(*,*) 'residue = ',RName
             call Finalize
           end if
         end do

         NumAtm(i) = NumAtm(i) + NumAdd

       end do

     end if

   end do

! ## Total Number of the Simulated System
   Number = 0

   do i = 1, NumSpec

     Number = Number + NumMol(i) * NumAtm(i)

   end do

   if(QPathInt) then

     allocate( Rpi(3,Number,Nbead) )
     allocate( Rnm(3,Number,Nbead) )
     allocate( Vnm(3,Number,Nbead) )

     allocate( Rcentroid(3,Number) )

     allocate( Fpi(3,Number,Nbead) )
     allocate( Fnm(3,Number,Nbead) )
     allocate( Fnm_ref(3,Number,Nbead) )

     allocate( NmMass(Number,Nbead) )
     allocate( FictMass(Number,Nbead) )
     allocate( InvFictMass(Number,Nbead) )

     allocate( Tnm(Nbead,Nbead) )
     allocate( InvTnm(Nbead,Nbead) )

   end if

   allocate( R(3,Number) )
   allocate( Vel(3,Number) )

   allocate( HamR(Number) )

   allocate( Frc_OptC(3,Number) )

   allocate( Frc_Bond (3,Number) )
   allocate( Frc_Angle(3,Number) )
   allocate( Frc_UB   (3,Number) )
   allocate( Frc_Dihed(3,Number) )
   allocate( Frc_Impro(3,Number) )

   allocate( Frc_Elec(3,Number) )
   allocate( Frc_Ersp(3,Number) )
   allocate( Frc_Eksp(3,Number) )

   allocate( Frc_EAM(3,Number) )

   allocate( Frc_NBlong(3,Number) )
   allocate( Frc_NBshrt(3,Number) )

   Frc_OptC   = 0.d0
   Frc_Bond   = 0.d0
   Frc_Angle  = 0.d0
   Frc_UB     = 0.d0
   Frc_Dihed  = 0.d0
   Frc_Impro  = 0.d0
   Frc_Elec   = 0.d0
   Frc_Ersp   = 0.d0
   Frc_Eksp   = 0.d0
   Frc_EAM    = 0.d0
   Frc_NBlong = 0.d0
   Frc_NBshrt = 0.d0

   Ene_OptC   = 0.d0
   Ene_Bond   = 0.d0
   Ene_Angle  = 0.d0
   Ene_UB     = 0.d0
   Ene_Dihed  = 0.d0
   Ene_Impro  = 0.d0
   Ene_Elec   = 0.d0
   Ene_Ersp   = 0.d0
   Ene_Eksp   = 0.d0
   Ene_EAM    = 0.d0
   Ene_NBlong = 0.d0
   Ene_NBshrt = 0.d0

   Vir_OptC   = 0.d0
   Vir_Bond   = 0.d0
   Vir_Angle  = 0.d0
   Vir_UB     = 0.d0
   Vir_Dihed  = 0.d0
   Vir_Impro  = 0.d0
   Vir_Ersp   = 0.d0
   Vir_Eksp   = 0.d0
   Vir_NBlong = 0.d0
   Vir_NBshrt = 0.d0

!## EAM >>

   if(ForceField(1:3) == 'EAM') then

     NumSpecPair = NumSpec*(NumSpec+1)/2

#ifdef MEAM

     NumEAMSpec = NumEAM*NumSpec

     allocate( Rrefa(Nref,NumEAMSpec))
     allocate( Rhorefa(Nref,NumEAMSpec) )
     allocate( Rhoy2(Nref,NumEAMSpec) )
     allocate( dRhodRref(Nref,NumEAMSpec) )
     allocate( dRhodRy2(Nref,NumEAMSpec) )
     allocate( Rmaxa(NumEAMSpec) )
     allocate( Rmina(NumEAMSpec) )

     allocate( Rhorefb(Nref,NumEAMSpec) )
     allocate( FRhoref(Nref,NumEAMSpec) )
     allocate( FRhoy2(Nref,NumEAMSpec) )
     allocate( dFdRhoref(Nref,NumEAMSpec) )
     allocate( dFdRhoy2(Nref,NumEAMSpec) )
     allocate( Rhomax(NumEAMSpec) )
     allocate( Rhomin(NumEAMSpec) )

#else

     allocate( Rrefa(Nref,NumSpec)) 
     allocate( Rhorefa(Nref,NumSpec) )
     allocate( Rhoy2(Nref,NumSpec) )
     allocate( dRhodRref(Nref,NumSpec) )
     allocate( dRhodRy2(Nref,NumSpec) )
     allocate( Rmaxa(NumSpec) )
     allocate( Rmina(NumSpec) )

     allocate( Rhorefb(Nref,NumSpec) )
     allocate( FRhoref(Nref,NumSpec) )
     allocate( FRhoy2(Nref,NumSpec) )
     allocate( dFdRhoref(Nref,NumSpec) )
     allocate( dFdRhoy2(Nref,NumSpec) )
     allocate( Rhomax(NumSpec) )
     allocate( Rhomin(NumSpec) )

#endif

     allocate( Rrefb(Nref,NumSpecPair) )
     allocate( PhiRef(Nref,NumSpecPair) )
     allocate( Phiy2(Nref,NumSpecPair) )
     allocate( dPhidRref(Nref,NumSpecPair) )
     allocate( dPhidRy2(Nref,NumSpecPair) )
     allocate( Rmaxb(NumSpecPair) )
     allocate( Rminb(NumSpecPair) )

     allocate( ResidPair( NumSpec,NumSpec) )

   end if

! ## << EAM

   if(QPBC) then

     allocate( ListIJ(3,MaxPair) )

     ListIJ = 0

   end if

   if(ForceField(1:2) == 'CG') then

     maxpsh = MaxPair / 3
     allocate( List_shortIJ(3,maxpsh) )

     List_shortIJ = 0

   end if

   if(QSHAKE) then

     allocate( Frc_Const(3,Number) )

     Frc_Const = 0.d0
     Vir_ConstR = 0.d0
     Vir_ConstV = 0.d0

   end if

end subroutine MemoryAllocation
