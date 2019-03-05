

!################################################################
!################################################################


! ***************************************************************
! **                                                           **
! **  MPDyn2                                  Version 2.00     **
! **                                                           **
! **  Author --- Wataru Shinoda                                **
! **                                                           **
! **  Department of Applied Chemistry, Nagoya University       **
! **  Furo-cho, Chikusa-ku, Nagoya, Aichi 464-8603 Japan       **
! **  e-mail : w.shinoda@apchem.nagoya-u.ac.jp                 **
! **                                         July 20, 2014     **
! **                                                           **
! ***************************************************************
!
!  This program is distributed WITHOUT ANY WARRANTY
!  See the GNU General Public License for more details.
!
! ***************************************************************
! **                                                           **
! **  Molecular simulations                                    **
! **  1. Molecular dynamics                                    **
! **  2. Hybrid Monte Carlo                                    **
! **  3. Path integral molecular dynamics                      **
! **  4. Centroid molecular dynamics                           **
! **     unit: Length[A] Time[ps]                              **
! **                                                           **
! **  Mesoscale simulation                                     **
! **  5. Dissipative particle dynamics                         **
! **     unit: reduced unit                                    **
! **                                                           **
! **  Analysis                                                 **
! **  a. Distribution functions (ex. g(r) etc. )               **
! **  b. Auto-correlation function (ex. MSD, velocity          **
! **     autocorrelation, etc. )                               **
! **  c. Cavity insertion (Chemical potential)                 **
! **  d. 2-dimensional Voronoi                                 **
! **                                                           **
! ***************************************************************
!
! ***************************************************************
! **                                                           **
! **  Develop History                                          **
! **  This program is an extension of the molecular simulation **
! **  package, MPDyn, which was developed at AIST.(-2004)      **
! **                                                           **
! ***************************************************************


!################################################################
!################################################################


program MPDyn2

use CommonBlocks, only : QMaster, Qdebug, QPathInt, SimMethod, QPDB

implicit NONE

character(len=17), save :: start_time, end_time

   call SetupMPI

   call Setup
   if(QMaster.and.Qdebug) write(*,*) 'OK:Setup'

   call PrepDyna

   if(QMaster) call Print_Ini                 !  print the initial conditions

   call InfoConf

   if(QMaster) call AcTime(start_time)

   if(SimMethod == 'MD') then

     call Select_Ensemble           ! Link to MD main part
     QPDB = .True.
     call SaveParam

   else if(QPathInt) then

     call Select_Ensemble_PI        ! Link to PIMD main part
     QPDB = .True.
     call SaveParam

   else if(SimMethod == 'HMC') then

     call Select_Ensemble_HMC       ! Link to HMC main part
     QPDB = .True.
     call SaveParam

   else if(SimMethod == 'Analysis') then

     call Select_Analyses           ! Link to main part

   else if(SimMethod == 'DPD') then

     call Select_DPDalgorithm

   else if(SimMethod == 'MM') then

     call MMcalc

   else if(SimMethod == 'NMA') then

     call NMAcalc

   else if(SimMethod == 'DynaLib') then

     call Select_Ensemble_DynaLib

   end if

   if(QMaster) call AcTime(end_time)

   call InfoConc(start_time,end_time)

   call FinMPI

end program MPDyn2
