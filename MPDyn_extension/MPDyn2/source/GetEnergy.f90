! ############################
! ############################
! ## SUBROUTINE LIST 
! ## -- GetEnergy_NV 
! ## -- GetEnergy_NP 
! ############################


!######################################################################
!######################################################################


subroutine GetEnergy_NV(Potential)

use CommonBlocks, only : QMaster, QCorrectCutoff, ForceField
use EAM_param, only : Vir_EAM, Ene_EAM
use EwaldParam, only : Vir_Eksp, Ene_Eksp, Ene_Eslf
use OptConstraintParam, only : Vir_OptC, Ene_OptC
use NonbondParam, only : Vir_Ersp, Vir_NBshrt, Vir_NBlong, &
&   Ene_LJ, Ene_Elec, Ene_Ersp, Ene_ELshrt, Ene_ELlong,    &
&   Ene_NBshrt, Ene_NBlong
use BondedParam, only : &
&   Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro, &
&   Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use TailCorrect, only : Ene_LJ_co
use ThermoData, only : Virial

implicit none

real(8) :: Potential

   Potential = 0.d0

   Ene_LJ = Ene_LJ + Ene_EAM

   if(ForceField(1:2) == 'CG') then
     Ene_Ersp = Ene_ELshrt + Ene_ELlong
     Ene_LJ   = Ene_NBshrt + Ene_NBlong
   end if

!----------------------------------------------------------------------
   call SumEnergy(Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro, &
   &              Ene_LJ, Ene_Elec, Ene_Ersp, Ene_Eksp, Ene_OptC)
!----------------------------------------------------------------------

   Virial = Vir_Bond  + Vir_Angle + Vir_UB                &
   &      + Vir_Dihed + Vir_Impro + Vir_Ersp + Vir_NBshrt &
   &      + Vir_Eksp  + Vir_OptC  + Vir_EAM  + Vir_NBlong

!-------------------------
   call SumVir( Virial )
!-------------------------

   if(QMaster) then

     Potential = Ene_Bond + Ene_Angle + Ene_UB   + Ene_Dihed + Ene_Impro + Ene_OptC &
     &         + Ene_LJ   + Ene_Ersp  + Ene_Eksp + Ene_Eslf

     if(QCorrectCutoff) then

       Potential = Potential + Ene_LJ_co

     end if

   end if

end subroutine GetEnergy_NV


!######################################################################
!######################################################################


subroutine GetEnergy_NP(EneSystem)

use CommonBlocks, only : QMaster, QCorrectCutoff, ForceField, &
&   cBarostatMethod
use EAM_param, only : Ene_EAM
use BathParam, only: Pressure_o, SigmaS
use EwaldParam, only : Ene_Eksp, Ene_Eslf
use OptConstraintParam, only : Ene_OptC
use NonbondParam, only : Ene_LJ, Ene_Elec, Ene_Ersp, Ene_ELshrt, &
&   Ene_ELlong, Ene_NBshrt, Ene_NBlong
use BondedParam, only : &
&   Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use TailCorrect, only : Ene_LJ_co
use CellParam, only : H, Volume

implicit none

real(8) :: Pot_p, Potential, EneSystem
! ## Stress >>
real(8), dimension(3,3) :: ElasM, Htrans, G
real(8) :: Pot_elastic
! ## << Stress

   EneSystem = 0.d0

   Ene_LJ = Ene_LJ + Ene_EAM

   if(ForceField(1:2) == 'CG') then
     Ene_Ersp = Ene_ELshrt + Ene_ELlong
     Ene_LJ   = Ene_NBshrt + Ene_NBlong
   end if

!----------------------------------------------------------------------
   call SumEnergy(Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro, &
   &              Ene_LJ, Ene_Elec, Ene_Ersp, Ene_Eksp, Ene_OptC)
!----------------------------------------------------------------------

   if(QMaster) then

!-------------------------
     call CalcTemp
!-------------------------

     Pot_p = Pressure_o * Volume

! ## Stress >>
     if( cBarostatMethod == 'ST' ) then
       Htrans = transpose( H )
       G = matmul(Htrans,H)
       ElasM = matmul(SigmaS,G)
       Pot_elastic = 0.5 * ( ElasM(1,1) + ElasM(2,2) + ElasM(3,3) )
       Pot_p = Pot_p + Pot_elastic
     end if
! ## << Stress

     Potential = Ene_Bond + Ene_Angle + Ene_UB   + Ene_Dihed + Ene_Impro + Ene_OptC &
     &         + Ene_LJ   + Ene_Ersp  + Ene_Eksp + Ene_Eslf

     if(QCorrectCutoff) then

       Potential = Potential + Ene_LJ_co

     end if

     EneSystem = Pot_p + Potential

   end if

end subroutine GetEnergy_NP
