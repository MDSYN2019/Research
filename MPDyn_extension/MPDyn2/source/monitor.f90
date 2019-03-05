! ############################
! ## SUBROUTINE LIST 
! ## -- Print_Ini 
! ## -- Print_Energy_NV 
! ## -- Print_Energy_NP 
! ## -- Print_Energy_iso 
! ## -- Print_Config 
! ## -- Print_Config_PI 
! ## -- Print_Velocity 
! ## -- Print_Velocity_PI 
! ## -- DPD_monitor 
! ## -- Print_Energy_DP 
! ## -- Print_Average_DP 
! ## -- RDF 
! ## -- Print_XYZ 
! ## -- Print_Energy_HMC_NV 
! ## -- Print_Energy_HMC_NP 
! ## -- Print_Energy_NV_PI 
! ## -- Print_Energy_NP_PI 
! ## -- Print_Energy_iso_PI 
! ## -- Monitor_Force 
! ## -- Print_Energy_MM 
! ## -- Print_Energy_MM_iso 
! ## -- Print_Energy_DynaLib_NV 
! ## -- Print_Energy_DynaLib_NP 
! ## -- Print_Energy_DynaLib_iso 
! ## -- Print_Eflux
! ############################


!#####################################################################
!#####################################################################


! ****************************************
! ** output initial condition           **
! ****************************************

subroutine Print_Ini

use Numbers, only : N, NumSpec, NumMol
use CommonBlocks, only : SimMethod, QAveTh, QPBC, QEps_r, Qstdout,        &
&   QPathInt, QThermostat, QBarostat, cBarostatMethod, cThermostatMethod, &
&   QATaup, QSHAKE, QRigidBody, QCoulomb, cCOULOMB, QOption, QREPWALL,    &
&   QCGWALL, QTICR, QVolScale, QOpFix, QJarzynski, ForceField
use CommonMPI
use ParamAnalyze
use CommonDPD
use CommonHMC
use PMEparam, only : Nfft, Bsp_order
use WallParam, only : kcap
use CommonPI, only : Nbead, deltat_ref, GammaC2
use BathParam
use SMDparam
use FEparam
use UnitExParam, only : rprs, pi, Avogadro, ExParam
use EwaldParam, only : Alpha, ih2mx, TmpNh, kmaxx, kmaxy, kmaxz
use OptConstraintParam, only : NHam, kHam, NumOptC, OptCI, OptCJ, rOptC
use CellParam, only : Vsc_Rate
use CutoffParam, only : Rcutoff2
use AtomParam, only : MolName
use TimeParam, only : Nstep, itgn, ixc, ixv, lp, lk, deltat
use ThermoData, only : Temp, CellMomentum

implicit NONE

integer :: i, j, k, l, ii
character(len=1) :: DIREC
integer, dimension(3), parameter :: IQ = (/11,12,6/)

   if(SimMethod == 'DynaLib') Return

   do k = 1, 3

#ifdef BMONI
     if(k==1.or.k==2) cycle
#endif

   if((.not.QAveTh).and.(k==2)) cycle

   ii = IQ(k)

   write(ii,'(a)') '########################################################################'
   write(ii,'(a)') '                                                                        '
   write(ii,'(a)') '                              MPDyn2                                    '
   write(ii,'(a)') '                                                                        '
   write(ii,'(a)') '########################################################################'

   if(SimMethod == 'MD') then
     write(ii,'(/a/)') '             ###  Molecular Dynamics Simulation ###                 '
   else if(SimMethod == 'MC') then
     write(ii,'(/a/)') '                ###  Monte Carlo Simulation ###                     '
   else if(SimMethod == 'HMC') then
     write(ii,'(/a/)') '             ###  Hybrid Monte Carlo Simulation ###                 '
   else if(SimMethod == 'Analysis') then
     write(ii,'(/a/)') '                     ### Data Analysis ###                          '
   else if(SimMethod == 'DPD') then
     write(ii,'(/a/)') '              ### Dissipative Particle Dynamics ###                 '
   else if(SimMethod == 'PIMD') then
     write(ii,'(/a/)') '        ###  Path Integral Molecular Dynamics Simulation ###        '
   else if(SimMethod == 'CMD') then
     write(ii,'(/a/)') '         ###  Centroid Molecular Dynamics Simulation ###            '
   else if(SimMethod == 'MM') then
     write(ii,'(/a/)') '                ###  MM Energy Calculation ###                      '
   else if(SimMethod == 'NMA') then
     write(ii,'(/a/)') '                ###  Normal Mode Analysis  ###                      '
   end if

   if(NProcs /= 1) then
     write(ii,'(a,i3,a/)') &
     &  '                Parallelized by using ',&
     &  NProcs,' CPUs                  '
   end if

   write(ii,'(15x,a)') '********* Simulated System *********'

   do i = 1 , NumSpec

     write(ii,'(16x,a9,i2,a3,a8,a1,i10)') 'Component',i,' : ',MolName(i),'=',NumMol(i)

   end do

   write(ii,'(15x,a/)') '************************************'

   write(ii,'(a )') '########################################################################'
   write(ii,'(a/)') '########################################################################'

   end do

   l = 2
   if(Qstdout) l = 3

   do k = 1, l

#ifdef BMONI
   if(k==1.or.k==2) cycle
#endif

   if((.not.QAveTh).and.(k==2)) cycle

   ii = IQ(k)

   if((SimMethod == 'MD').or.(QPathInt)) then

     write(ii,'(/5x,a)')  '-----------------------------'
     write(ii,'(5x,a,a)') ' Force Field = ',ForceField
     write(ii,'(5x,a/)')  '-----------------------------'

     if(QPBC) then

       write(ii,'(a/)') 'Three-dimensional Periodic Boundary Condition'

     else

       if(QEps_r) then
         write(ii,'(a/)') 'Isolated System  ::  Eps_r = 4R '
       else
         write(ii,'(a/)') 'Isolated System  ::  Eps_r = 1 '
       end if

     end if

     if(QThermostat) then

       if(QBarostat) then
         write(ii,'(a/)')            'Isothermal-Isobaric Ensemble'

         if((cBarostatMethod == 'PR').or.(cBarostatMethod == 'ST')) then

           write(ii,'(a)')           'Parrinello-Rahman  -  Fully Flexible Cell!'

         else if(cBarostatMethod == 'AN') then

           write(ii,'(a)')           'Andersen  -  Isotropic Cell Fluctuation!'

         else if(cBarostatMethod == 'A3') then

           write(ii,'(a)')           'Parrinello-Rahman  -  Freezed angles!'

         else if(cBarostatMethod == 'A2') then

           if(CoupleEdge(2)==2) then
             write(ii,'(a)')           'Parrinello-Rahman  -  Pzz, Pxx=Pyy, Freezed angles!'
           else if(CoupleEdge(1)==2) then
             write(ii,'(a)')           'Parrinello-Rahman  -  Pxx, Pyy=Pzz, Freezed angles!'
           else
             write(ii,'(a)')           'Parrinello-Rahman  -  Pyy, Pxx=Pzz, Freezed angles!'
           end if

         end if

         if(cThermostatMethod == 'NH') then

           write(ii,'(a/)') 'Nose-Hoover method'

         else if(cThermostatMethod == 'NHC') then

           write(ii,'(a,i2,a/)') 'Nose-Hoover chain method (Chain Length=',NHchain,')'

         else if(cThermostatMethod == 'MNHC') then

           write(ii,'(a,i2,a/)') 'Massive Nose-Hoover chain method (Chain Length=',NHchain,')'

         else if(cThermostatMethod == 'VSCALE') then

           write(ii,'(a/)') 'Velocity scaling method'

         end if

         write(ii,'(5x,a,f10.3,a)')  'P       = ', Pressure_o/rprs*1.d-6, '[MPa]'
         write(ii,'(5x,a,f10.1,a/)') 'T       = ', Temp_o,        '[K]'

         if(cBarostatMethod == 'ST') then
           write(ii,'(/5x,a)') ' *******  External Stress [GPa] ******* '
           write(ii,'(7x,3f12.3)')  Stress_ext(1,:)/rprs*1.d-9
           write(ii,'(7x,3f12.3)')  Stress_ext(2,:)/rprs*1.d-9
           write(ii,'(7x,3f12.3)')  Stress_ext(3,:)/rprs*1.d-9
           write(ii,'(5x,a/)') ' ************************************** '
         end if

         write(ii,'(5x,a)')          'Bath Parameters :: '

         if(QATaup) then

           write(ii,'(5x,a)') '## Anisotropic Tau_barostat ##'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,xx= ', Tau_p * prefTaup(1,1), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,yy= ', Tau_p * prefTaup(2,2), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,zz= ', Tau_p * prefTaup(3,3), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,xy= ', Tau_p * prefTaup(1,2), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,yz= ', Tau_p * prefTaup(2,3), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,zx= ', Tau_p * prefTaup(3,1), '[ps]'

         else

           write(ii,'(5x,a,f10.3,a)')  '     Tau_p   = ', Tau_p , '[ps]'

         end if

         if(cThermostatMethod == 'NH') then

           write(ii,'(5x,a,f10.3,a)')  '     Tau_s   = ', Tau_s0, '[ps]'

         else if((cThermostatMethod == 'NHC').or.(cThermostatMethod == 'MNHC')) then

           write(ii,'(5x,a,f10.3,a)')  '     Tau_s0  = ', Tau_s0, '[ps]'
           write(ii,'(5x,a,f10.3,a)')  '     Tau_s1  = ', Tau_s1, '[ps]'

         end if

         if(QATaup) then

           write(ii,'(5x,a,d14.5)')    '     M barostat,xx= ',Mp(1,1)
           write(ii,'(5x,a,d14.5)')    '     M barostat,yy= ',Mp(2,2)
           write(ii,'(5x,a,d14.5)')    '     M barostat,zz= ',Mp(3,3)
           write(ii,'(5x,a,d14.5)')    '     M barostat,xy= ',Mp(1,2)
           write(ii,'(5x,a,d14.5)')    '     M barostat,yz= ',Mp(2,3)
           write(ii,'(5x,a,d14.5)')    '     M barostat,zx= ',Mp(3,1)

         else

           write(ii,'(5x,a,d14.5)')    '     M barostat   = ',Mp(1,1)

         end if

         if(cThermostatMethod == 'NH') then

           write(ii,'(5x,a,d14.5/)')  '     M thermostat = ',Mts(1)

         else if(cThermostatMethod == 'NHC') then

           write(ii,'(5x,a,2d14.5/)')  '     M thermostat = ',Mts(1),Mts(2)

         else if(cThermostatMethod == 'MNHC') then

           write(ii,'(5x,a,2d14.5/)')  '     M thermostat = ',MMNHC(1,1),MMNHC(2,1)

         end if

         if(QSHAKE) then

           write(ii,'(a/)')  'SHAKE/ROLL, RATTLE/ROLL'

         else if(QRigidBody) then

           write(ii,'(a/)')  'Rigid-Body Model'

         else

           write(ii,'(a/)')  'Flexible Model'

         end if

       else

         write(ii,'(a/)')            'Canonical Ensemble'

         if(cThermostatMethod == 'NH') then

           write(ii,'(a/)') 'Nose-Hoover method'

         else if(cThermostatMethod == 'NHC') then

           write(ii,'(a,i2,a/)') 'Nose-Hoover chain method (Chain Length=',NHchain,')'

         else if(cThermostatMethod == 'MNHC') then

           write(ii,'(a,i2,a/)') 'Massive Nose-Hoover chain method (Chain Length=',NHchain,')'

         else if(cThermostatMethod == 'VSCALE') then

           write(ii,'(a/)') 'Velocity scaling method'

         end if

         write(ii,'(5x,a,f10.1,a/)') 'T       = ', Temp_o,        '[K]'
         write(ii,'(5x,a)')          'Bath Parameters :: '

         if(cThermostatMethod == 'NH') then

           write(ii,'(5x,a,f10.3,a)')  '     Tau_s0  = ', Tau_s0, '[ps]'
           write(ii,'(5x,a,d14.5/)')  '     M thermostat = ',Mts(1)

         else if(cThermostatMethod == 'NHC') then

           write(ii,'(5x,a,f10.3,a)')  '     Tau_s0  = ', Tau_s0, '[ps]'
           write(ii,'(5x,a,f10.3,a)')  '     Tau_s1  = ', Tau_s1, '[ps]'
           write(ii,'(5x,a,2d14.5/)')  '     M thermostat = ',Mts(1),Mts(2)

         else if(cThermostatMethod == 'MNHC') then

           write(ii,'(5x,a,f10.3,a)')  '     Tau_s0  = ', Tau_s0, '[ps]'
           write(ii,'(5x,a,f10.3,a)')  '     Tau_s1  = ', Tau_s1, '[ps]'
           write(ii,'(5x,a,2d14.5/)')  '     M thermostat = ',MMNHC(1,1),MMNHC(2,1)

         end if

         if(QSHAKE) then

           write(ii,'(a/)')  'SHAKE/RATTLE'

         else if(QRigidBody) then

           write(ii,'(a/)')  'Rigid-Body Model'

         else

           write(ii,'(a/)')  'Flexible Model'

         end if

       end if

     else

       if(QBarostat) then

         write(ii,'(a/)')           'NPH Ensemble'

         if((cBarostatMethod == 'PR').or.(cBarostatMethod == 'ST')) then

           write(ii,'(a)')           'Parrinello-Rahman  -  Fully Flexible Cell!'

         else if(cBarostatMethod == 'AN') then

           write(ii,'(a)')           'Andersen  -  Isotropic Cell Fluctuation!'

         else if(cBarostatMethod == 'A3') then

           write(ii,'(a)')           'Parrinello-Rahman  -  Freezed angles!'

         else if(cBarostatMethod == 'A2') then

           write(ii,'(a)')           'Parrinello-Rahman  -  Pzz, Pxx=Pyy, Freezed angles!'

         end if

         write(ii,'(5x,a,f10.3,a)') 'P       = ', Pressure_o/rprs*1.d-6, '[MPa]'

         if(cBarostatMethod == 'ST') then
           write(ii,'(/5x,a)') ' *******  External Stress [GPa] ******* '
           write(ii,'(7x,3f12.3)')  Stress_ext(1,:)/rprs*1.d-9
           write(ii,'(7x,3f12.3)')  Stress_ext(2,:)/rprs*1.d-9
           write(ii,'(7x,3f12.3)')  Stress_ext(3,:)/rprs*1.d-9
           write(ii,'(5x,a/)') ' ************************************** '
         end if

         write(ii,'(5x,a)')         'Bath Parameters :: '

         if(QATaup) then

           write(ii,'(5x,a)') '## Anisotropic Tau_barostat ##'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,xx= ', Tau_p * prefTaup(1,1), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,yy= ', Tau_p * prefTaup(2,2), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,zz= ', Tau_p * prefTaup(3,3), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,xy= ', Tau_p * prefTaup(1,2), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,yz= ', Tau_p * prefTaup(2,3), '[ps]'
           write(ii,'(5x,a,d10.3,a)')  '     Tau_p,zx= ', Tau_p * prefTaup(3,1), '[ps]'

         else

           write(ii,'(5x,a,f10.3,a)')  '     Tau_p   = ', Tau_p , '[ps]'

         end if

         if(QATaup) then

           write(ii,'(5x,a,d14.5)')    '     M barostat,xx= ',Mp(1,1)
           write(ii,'(5x,a,d14.5)')    '     M barostat,yy= ',Mp(2,2)
           write(ii,'(5x,a,d14.5)')    '     M barostat,zz= ',Mp(3,3)
           write(ii,'(5x,a,d14.5)')    '     M barostat,xy= ',Mp(1,2)
           write(ii,'(5x,a,d14.5)')    '     M barostat,yz= ',Mp(2,3)
           write(ii,'(5x,a,d14.5)')    '     M barostat,zx= ',Mp(3,1)

         else

           write(ii,'(5x,a,d14.5)')    '     M barostat   = ',Mp(1,1)

         end if

         if(QSHAKE) then

           write(ii,'(a/)')  'SHAKE/ROLL, RATTLE/ROLL'

         else if(QRigidBody) then

           write(ii,'(a/)')  'Rigid-Body Model'

         else

           write(ii,'(a/)')  'Flexible Model'

         end if

       else

         write(ii,'(a/)')  'Microcanonical Ensemble'

         if(QSHAKE) then

           write(ii,'(a/)')  'SHAKE/RATTLE'

         else if(QRigidBody) then

           write(ii,'(a/)')  'Rigid-Body Model'

         else

           write(ii,'(a/)')  'Flexible Model'

         end if

       end if

     end if

     write(ii,'(a)') 'Simulation Condition Details'

     if(QPBC) then

       write(ii,'(5x,a,f10.2,a)') 'Rcutoff = ', sqrt(Rcutoff2),       '[A]'

     end if

     write(ii,'(5x,a,f10.5,a)') 'deltat  = ', deltat   ,' [ps] (fast)'
     write(ii,'(15x, f10.5,a)')               deltat*lp,' [ps] (intermediate)'
     write(ii,'(15x, f10.5,a)')               deltat*lk,' [ps] (slow)'
     write(ii,'(5x,a,i10)')     'Steps   = ', Nstep
     write(ii,'(5x,a,i8)')      '   output freq. (Energy, H etc. ) = ', itgn
     write(ii,'(5x,a,i8)')      '                (Configuration  ) = ', ixc
     write(ii,'(5x,a,i8)')      '                (Velocity       ) = ', ixv

     if(QPBC.and.QCoulomb) then

       if( trim(cCOULOMB) == 'EWALD' ) then
         write(ii,'(/5x,a)')       'Ewald Parameters :: '
         write(ii,'(5x,a,f6.3)')   '     Alpha  = ', Alpha
         write(ii,'(5x,a,i6)')     '     ih2max = ', ih2mx
         write(ii,'(5x,a,i6)')     '     Nh     = ', TmpNh
         write(ii,'(5x,a,3i4)')    '     kMax   = ', kmaxx,kmaxy,kmaxz
       else if( trim(cCOULOMB) == 'PME' ) then
         write(ii,'(/5x,a)')       'PMEwald Parameters :: '
         write(ii,'(5x,a,f6.3)')   '     Alpha    = ', Alpha
         write(ii,'(5x,a,i6)')     '     Grid_x   = ', Nfft(1)
         write(ii,'(5x,a,i6)')     '     Grid_y   = ', Nfft(2)
         write(ii,'(5x,a,i6)')     '     Grid_z   = ', Nfft(3)
         write(ii,'(5x,a,i6)')     '     B-spline = ', Bsp_order
       end if

     end if

     if( NHam /= 0 ) then

       write(ii,'(/5x,a)') '%%%% Some atoms are fixed! %%%%'
       write(ii,'(5x,a,f5.1,a)')  '% Harmonic Constraint k=',kHam/ExParam,' %'
       write(ii,'(5x,a,i10,a)') '%      ',      NHam,' Atoms       %'
       write(ii,'(5x,a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

     end if

     if(QOption) then

       write(ii,'(/5x,a)') '%%%%% Optional constraints ! %%%%%'
       do i = 1, NumOptC
         write(ii,'(5x,a,i5,a,i5,3x,f10.3,a)') &
         &     '% Atom',OptCI(i),' - ',OptCJ(i),rOptC(i),' %'
       end do
       write(ii,'(5x,a)')  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

     end if

     if(QREPWALL) then

       write(ii,'(/5x,a)') '%%%%% Wall potential ! %%%%%%'
       write(ii,'(5x,a,f8.3,a)')  '%      k=', kcap/ExParam,'           %'
       write(ii,'(5x,a)')  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

     end if

     if(QCGWALL) then

       write(ii,'(/5x,a)')        '%%%%% CG Wall potential ! %%%%%%'
       write(ii,'(5x,a)')         '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

     end if

     if(QTICR) then

       write(ii,'(/5x,a)') '################## Thermodynamic Integration #################'
       write(ii,'( 5x,a)') '#                                                            #'
       write(ii,'( 5x,a)') '#    TI calculation for molecule creation or annihilation    #'
       write(ii,'( 5x,a)') '#                                                            #'
       if(QCreation) then
       write(ii,'( 5x,a)') '#    Action : CREATION                                       #'
       else
       write(ii,'( 5x,a)') '#    Action : ANNIHILATION                                   #'
       end if
       if(QLJ) then
       if(QCharge) then
       write(ii,'( 5x,a)') '#    Interaction : LJ & COULOMB                              #'
       else
       write(ii,'( 5x,a)') '#    Interaction : LJ                                        #'
       end if
       else
       write(ii,'( 5x,a)') '#    Interaction : COULOMB                                   #'
       end if
       write(ii,'( 5x,a)') '#                                                            #'
       write(ii,'( 5x,a,a,a,i4,a)') '#    TARGET MOLECULE : ',MolName(NTIspec), &
       &                            ' MOL NUMBER : ',NumTIMol,'            #'
       write(ii,'( 5x,a)') '#                                                            #'
       write(ii,'( 5x,a)') '#------------------------------------------------------------#'
       write(ii,'( 5x,a)') '#                                                            #'
       write(ii,'( 5x,a,i4,a)') '# Steps for lambda values = ',iTIstep,&
       &                        '                             #'
       write(ii,'( 5x,a)') '#                                                            #'
       do i = 1, iTIstep
         write(ii,'( 5x,a,i2,a,f7.3,a,i7,a,i7,a)')  &
         &  '# Step-',i,', lambda = ',TIlambda(i),  &
         &  ',  Equil.:',iequil_TICR(i),', Sampl.:',&
         &  isample_TICR(i),' #'
       end do
       write(ii,'( 5x,a)') '#                                                            #'
       write(ii,'( 5x,a)') '##############################################################'

     end if

     if(QVolScale) then
       write(ii,'(/5x,a)') '############################################################'
       write(ii,'( 5x,a)') '#                                                          #'
       write(ii,'( 5x,a)') '#    Cell-size scaling was selected : NVT ensemble         #'
       write(ii,'( 5x,a)') '#                                                          #'
       write(ii,'( 5x,a,f9.3,a)') '#       Scaling rate in x-axis =  ',Vsc_Rate(1), &
       &   '  [A/ps]        #'
       write(ii,'( 5x,a,f9.3,a)') '#       Scaling rate in y-axis =  ',Vsc_Rate(2), &
       &   '  [A/ps]        #'
       write(ii,'( 5x,a,f9.3,a)') '#       Scaling rate in z-axis =  ',Vsc_Rate(3), &
       &   '  [A/ps]        #'
       write(ii,'( 5x,a)') '#                                                          #'
       write(ii,'(5x,a/)') '############################################################'
     end if

     if(QOpFix) then
       write(ii,'(/5x,a)') '############################################################'
       write(ii,'( 5x,a)') '#                                                          #'
       write(ii,'( 5x,a)') '#    Optional Constraints was applied for :                #'
       write(ii,'( 5x,a)') '#                                                          #'
       do i = 1, NumOpFixS
         j=FixDir(i)
         select case(j)
         case(1)
           DIREC='X'
         case(2)
           DIREC='Y'
         case(3)
           DIREC='Z'
         end select
         write(ii,'(5x,a,i5,a,i5,a,a,a,f7.1,a)') '#    SNGL   (',NumFixF(i),'-',NumFixL(i),&
         &     ')  ',DIREC,'=',RGConst(i),'                       #'
       end do
       do i = 1, NumOpFixP
         write(ii,'(5x,a,i5,a,i5,a,i5,a,i5,a,f7.2,a)') '#    PAIR  (',NumFixFi(i),'-',&
         &     NumFixLi(i),')  (',NumFixFj(i),'-',NumFixLj(i),')   R0=',ConstDis(i),&
         &     '       #'
       end do
       write(ii,'( 5x,a)') '#                                                          #'
       write(ii,'( 5x,a,i9,a)') '#  Monitoring constraint force every ',NumFreqConst/lk,&
       &                        ' step        #'
       write(ii,'( 5x,a)') '#                                                          #'
       write(ii,'(5x,a/)') '############################################################'
     end if

     if(QJarzynski) then
       write(ii,'(/5x,a)') '############################################################'
       write(ii,'( 5x,a)') '#                                                          #'
       write(ii,'( 5x,a)') '#    Steering MD for a free energy calculation using       #'
       write(ii,'( 5x,a)') "#    Jarzynski's equality                                  #"
       write(ii,'( 5x,a)') '#                                                          #'
       do i = 1, NumOpFixS
       write(ii,'(5x,a,i5,a,i5,a,a,a,f7.1,a)') '#    SNGL   (',NumFixF(i),'-',NumFixL(i),&
       &     ')                                  #'
       end do
       do i = 1, NumOpFixP
       write(ii,'(5x,a,i5,a,i5,a,i5,a,i5,a,f7.2,a)') '#    PAIR  (',NumFixFi(i),'-',&
       &     NumFixLi(i),')  (',NumFixFj(i),'-',NumFixLj(i),')                    #'
       end do
       write(ii,'( 5x,a)') '#                                                          #'
       write(ii,'( 5x,a,f10.3,a)') '#    Spring constant  = ',k_steering/ExParam,&
       & ' [kcal/mol/A^2]          #'
       write(ii,'( 5x,a,f10.3,a)') '#       (unit convert)  ',&
       & k_steering/ExParam*4.184d+25/Avogadro,' [pN/A]                  #'
       write(ii,'( 5x,a,f10.3,a)') '#    Initial position = ',R0_steering,&
       & ' [A]                     #'
       if(NumOpFixS/=0) then
       write(ii,'( 5x,a,3f10.7,a)') '#    along the vector = (',Uvec_steering(:),')   #'
       end if
       write(ii,'( 5x,a,f10.3,a)') '#    Pulling rate = ',Vshift_steering/deltat*1.d+03,&
       &  ' [A/ns]                      #'
       write(ii,'( 5x,a,i9,a)') '#  Monitoring the work W(0->t) every ',NumFreqConst/lk,&
       &                        ' step        #'
       write(ii,'( 5x,a)') '#                                                          #'
       write(ii,'(5x,a/)') '############################################################'
     end if

     write(ii,'(/5x,a,3d8.1)')   'Cell Momentum = ',(CellMomentum(i),i=1,3)

     write(ii,'(/5x,a,f7.2/)')   'Initial Temperature = ', Temp

   else if(SimMethod == 'HMC') then

     write(ii,'(/5x,a)')  '-----------------------------'
     write(ii,'(5x,a,a)') ' Force Field = ',ForceField
     write(ii,'(5x,a/)')  '-----------------------------'

     if(QPBC) then

       write(ii,'(a/)') 'Three-dimensional Periodic Boundary Condition'

     else

       if(QEps_r) then
         write(ii,'(a/)') 'Isolated System  ::  Eps_r = 4R '
       else
         write(ii,'(a/)') 'Isolated System  ::  Eps_r = 1 '
       end if

     end if

     if(QBarostat) then
       write(ii,'(a/)')            'Isothermal-Isobaric Ensemble'

       if((cBarostatMethod == 'PR').or.(cBarostatMethod == 'ST')) then

         write(ii,'(a)')           'Parrinello-Rahman  -  Fully Flexible Cell!'

       else if(cBarostatMethod == 'AN') then

         write(ii,'(a)')           'Andersen  -  Isotropic Cell Fluctuation!'

       else if(cBarostatMethod == 'A3') then

         write(ii,'(a)')           'Parrinello-Rahman  -  Freezed angles!'

       else if(cBarostatMethod == 'A2') then

         write(ii,'(a)')           'Parrinello-Rahman  -  Pzz, Pxx=Pyy, Freezed angles!'

       end if

       write(ii,'(5x,a,f10.3,a)')  'P       = ', Pressure_o/rprs*1.d-6, '[MPa]'
       write(ii,'(5x,a,f10.1,a/)') 'T       = ', Temp_o,        '[K]'

       if(cBarostatMethod == 'ST') then
         write(ii,'(/5x,a)') ' *******  External Stress [GPa] ******* '
         write(ii,'(7x,3f12.3)')  Stress_ext(1,:)/rprs*1.d-9
         write(ii,'(7x,3f12.3)')  Stress_ext(2,:)/rprs*1.d-9
         write(ii,'(7x,3f12.3)')  Stress_ext(3,:)/rprs*1.d-9
         write(ii,'(5x,a/)') ' ************************************** '
       end if

       write(ii,'(5x,a)')          'Bath Parameters :: '

       if(QATaup) then

         write(ii,'(5x,a)') '## Anisotropic Tau_barostat ##'
         write(ii,'(5x,a,d10.3,a)')  '     Tau_p,xx= ', Tau_p * prefTaup(1,1), '[ps]'
         write(ii,'(5x,a,d10.3,a)')  '     Tau_p,yy= ', Tau_p * prefTaup(2,2), '[ps]'
         write(ii,'(5x,a,d10.3,a)')  '     Tau_p,zz= ', Tau_p * prefTaup(3,3), '[ps]'
         write(ii,'(5x,a,d10.3,a)')  '     Tau_p,xy= ', Tau_p * prefTaup(1,2), '[ps]'
         write(ii,'(5x,a,d10.3,a)')  '     Tau_p,yz= ', Tau_p * prefTaup(2,3), '[ps]'
         write(ii,'(5x,a,d10.3,a)')  '     Tau_p,zx= ', Tau_p * prefTaup(3,1), '[ps]'

       else

         write(ii,'(5x,a,f10.3,a)')  '     Tau_p   = ', Tau_p , '[ps]'

       end if

       if(QATaup) then

         write(ii,'(5x,a,d14.5)')    '     M barostat,xx= ',Mp(1,1)
         write(ii,'(5x,a,d14.5)')    '     M barostat,yy= ',Mp(2,2)
         write(ii,'(5x,a,d14.5)')    '     M barostat,zz= ',Mp(3,3)
         write(ii,'(5x,a,d14.5)')    '     M barostat,xy= ',Mp(1,2)
         write(ii,'(5x,a,d14.5)')    '     M barostat,yz= ',Mp(2,3)
         write(ii,'(5x,a,d14.5)')    '     M barostat,zx= ',Mp(3,1)

       else

         write(ii,'(5x,a,d14.5)')    '     M barostat   = ',Mp(1,1)

       end if

       if(QRigidBody) then

         write(ii,'(a/)')  'Rigid-Body Model'

       else

         write(ii,'(a/)')  'Flexible Model'

       end if

     else

       write(ii,'(a/)')            'Canonical Ensemble'

       write(ii,'(5x,a,f10.1,a/)') 'T       = ', Temp_o,        '[K]'

       if(QRigidBody) then

         write(ii,'(a/)')  'Rigid-Body Model'

       else

         write(ii,'(a/)')  'Flexible Model'

       end if

     end if

     write(ii,'(a)') 'Simulation Condition Details'

     if(QPBC) then

       write(ii,'(5x,a,f10.2,a)') 'Rcutoff = ', sqrt(Rcutoff2),       '[A]'

     end if

     write(ii,'(5x,a,f10.5,a)') 'deltat  = ', deltat   ,' [ps] (fast)'
     write(ii,'(15x, f10.5,a)')               deltat*lp,' [ps] (intermediate)'
     write(ii,'(15x, f10.5,a)')               deltat*lk,' [ps] (slow)'
     write(ii,'(5x,a,i10)')     'Steps   = ', Nstep
     write(ii,'(5x,a,i8)')      '   output freq. (Energy, H etc. ) = ', itgn
     write(ii,'(5x,a,i8)')      '                (Configuration  ) = ', ixc
     write(ii,'(5x,a,i8)')      '                (Velocity       ) = ', ixv
     write(ii,'(5x,a,i8/)')      '   MD steps per MC sweep          = ', MDstep
     if(QPartial) then
        write(ii,'(5x,a,f8.1, a)')      ' Partial Momentum Refreshment : theta = ', &
        &    ThetaPartial * 180.d0 / pi,' [degree] '
     end if

     if(QPBC.and.QCoulomb) then

       if( trim(cCOULOMB) == 'EWALD' ) then
         write(ii,'(/5x,a)')       'Ewald Parameters :: '
         write(ii,'(5x,a,f6.3)')   '     Alpha  = ', Alpha
         write(ii,'(5x,a,i6)')     '     ih2max = ', ih2mx
         write(ii,'(5x,a,i6)')     '     Nh     = ', TmpNh
         write(ii,'(5x,a,3i4)')    '     kMax   = ', kmaxx,kmaxy,kmaxz
       else if( trim(cCOULOMB) == 'PME' ) then
         write(ii,'(/5x,a)')       'PMEwald Parameters :: '
         write(ii,'(5x,a,f6.3)')   '     Alpha    = ', Alpha
         write(ii,'(5x,a,i6)')     '     Grid_x   = ', Nfft(1)
         write(ii,'(5x,a,i6)')     '     Grid_y   = ', Nfft(2)
         write(ii,'(5x,a,i6)')     '     Grid_z   = ', Nfft(3)
         write(ii,'(5x,a,i6)')     '     B-spline = ', Bsp_order
       end if

     end if

     if( NHam /= 0 ) then

       write(ii,'(/5x,a)') '%%%% Some atoms are fixed! %%%%'
       write(ii,'(5x,a,f5.1,a)')  '% Harmonic Constraint k=',kHam/ExParam,' %'
       write(ii,'(5x,a,i10,a)') '%      ',      NHam,' Atoms       %'
       write(ii,'(5x,a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

     end if

     if(QOption) then

       write(ii,'(/5x,a)') '%%%%% Optional constraints ! %%%%%'
       do i = 1, NumOptC
         write(ii,'(5x,a,i5,a,i5,3x,f10.3,a)') &
         &     '% Atom',OptCI(i),' - ',OptCJ(i),rOptC(i),' %'
       end do
       write(ii,'(5x,a)')  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

     end if

     if(QREPWALL) then

       write(ii,'(/5x,a)') '%%% repulsive wall! %%%%'
       write(ii,'(5x,a,f5.1,a)')  '%    k=', kcap/ExParam,'           %'
       write(ii,'(5x,a)')  '%%%%%%%%%%%%%%%%%%%%%%%%'

     end if

     if(QCGWALL) then

       write(ii,'(/5x,a)')        '%%%%% CG Wall potential ! %%%%%%'
       write(ii,'(5x,a)')         '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

     end if

   else if(SimMethod == 'DPD') then

     write(ii,'(/5x,a)')  '********************************************'
     write(ii,'(5x,a,a)') ' Integration of the EOM by ', IntegrMethod
     write(ii,'(5x,a)')  '********************************************'

     write(ii,'(/5x,a)')  '-----------------------------'
     write(ii,'(5x,a,a)') ' Force Type = ',ForceField
     write(ii,'(5x,a/)')  '-----------------------------'

     if( ForceField == 'Morse' ) then

       write(ii,'(5x,a,f5.1)') 'Morse function parameter alpha = ',a_morse
       do i = 1, TypeNum(N)
         write(ii,'(5x,7(f8.2,2x))') e_morse(i,:)
       end do

     else
       write(ii,'(5x,a)') ' %%% Interaction parameter matrix A %%% '
       do i = 1 , TypeNum(N)
         write(ii,'(5x,7(f8.1,2x))') a(i,:)
       end do

     end if

     write(ii,'(/a)') 'Simulation Condition Details'
     write(ii,'(5x,a,f10.1)') 'T       = ', Temp_o
     write(ii,'(5x,a,f10.2)') 'gamma   = ', gamma
     write(ii,'(5x,a,f10.2)') 'sigma   = ', sigma
     write(ii,'(5x,a,f10.5)') 'deltat  = ', deltat
     write(ii,'(5x,a,i10)')   'Steps   = ', Nstep
     write(ii,'(5x,a,i8)')    '   output freq. (Energy, H etc. ) = ', itgn
     write(ii,'(5x,a,i8)')    '                (Configuration  ) = ', ixc

     write(ii,'(/5x,a,f7.2/)')   'Initial Temperature = ', Temp

   else if(SimMethod == 'Analysis'.and.Qstdout) then

     if(cWhat=='ElecZ') then
       write(ii,'(/a)') '**************************************************'
       write(ii,'(a)')  '**  Electric potential along the bilayer normal **'
       write(ii,'(a/)') '**************************************************'
       write(ii,'(a )') '--> Output file are ./Analy/ElecRhoS.dat'
       write(ii,'(a )') '                    ./Analy/ElecPolS.dat'
       write(ii,'(a )') '                    ./Analy/ElecPhiS.dat'
     end if

     if(cWhat=='TraceG') then
       write(ii,'(/a)') '*********************************'
       write(ii,'(a)')  '**  Trace the COM of molecule  **'
       write(ii,'(a/)') '*********************************'
       write(ii,'(a/)') '--> Output file is ./Analy/TraceG.dat'
     end if

     if(cWhat=='Edens') then
       write(ii,'(/a)') '********************************************************'
       write(ii,'(a)')  '**  Electron density profile along the bilayer normal **'
       write(ii,'(a/)') '********************************************************'
       write(ii,'(a )') '--> Output file is ./Analy/ElectronDensityZ.dat'
     end if

     if(cWhat=='OriChain') then
       write(ii,'(/a)') '********************************************************'
       write(ii,'(a)')  '**  distribution and time correlation of orientation  **'
       write(ii,'(a)')  '**  of alkyl chain of the lipid molecules(DPPC,DPhPC) **'
       write(ii,'(a/)') '********************************************************'
       write(ii,'(a )') '--> Output file are ./Analy/DisChain_*_*.dat  '
       write(ii,'(a )') '                    ./Analy/CorrChain_*_*.dat '
     end if

     if(cWhat=='OriLipid') then
       write(ii,'(/a)') '********************************************************'
       write(ii,'(a)')  '**  distribution and time correlation of orientation  **'
       write(ii,'(a)')  '**  of alkyl chain of the lipid molecules(DPPC,DPhPC) **'
       write(ii,'(a/)') '********************************************************'
       write(ii,'(a )') '--> Output file are ./Analy/PN.dat '
       write(ii,'(a )') '                    ./Analy/CC.dat '
       write(ii,'(a )') '                    ./Analy/CO.dat '
     end if

     if(cWhat=='RotLipid') then
       write(ii,'(/a)') '********************************************************'
       write(ii,'(a)')  '**  distribution and time correlation of rotaion of   **'
       write(ii,'(a)')  '**  the lipid molecules(DPPC,DPhPC)                   **'
       write(ii,'(a/)') '********************************************************'
       write(ii,'(a )') '--> Output file are ./Analy/CCCC.dat '
       write(ii,'(a )') '                    ./Analy/COCO.dat '
     end if

     if(cWhat=='LenChain') then
       write(ii,'(/a)') '********************************************************'
       write(ii,'(a)')  '**  distribution and time correlation of length       **'
       write(ii,'(a)')  '**  of alkyl chain of the lipid molecules(DPPC,DPhPC) **'
       write(ii,'(a/)') '********************************************************'
       write(ii,'(a )') '--> Output file are ./Analy/DisChainLength_*_*.dat  '
       write(ii,'(a )') '                    ./Analy/CorrChainLength_*_*.dat '
     end if

     if(cWhat=='LipidChAng') then
       write(ii,'(/a)') '********************************************************'
       write(ii,'(a)')  '**  distribution of the angle between two alkyl chains**'
       write(ii,'(a)')  '**  in a lipid molecule (DPPC,DPhPC)                  **'
       write(ii,'(a/)') '********************************************************'
       write(ii,'(a/)') '--> Output file is ./Analy/CoAngChain.dat  '
     end if

     if(cWhat=='dPP') then
       write(ii,'(/a)') '********************************************'
       write(ii,'(a)')  '**  P-P Distance in the lipid bilayer     **'
       write(ii,'(a)')  '**  NOTE: only DPPC, DMPC, PhPC(DPhPC)!!  **'
       write(ii,'(a/)') '********************************************'
       write(ii,'(a/)') '--> Output file is ./Analy/dPP.dat '
     end if

     if(cWhat=='SCD') then
       write(ii,'(/a)') '**************************************************'
       write(ii,'(a)')  '**  order parameter of the acyl chains of LIPID **'
       write(ii,'(a)')  '**  NOTE: only DPPC, DMPC, PhPC(DPhPC)!!        **'
       write(ii,'(a/)') '**************************************************'
       write(ii,'(a )') '--> Output file is ./Analy/SCD.dat '
     end if

     if(cWhat=='SZZ') then
       write(ii,'(/a)') '**************************************************'
       write(ii,'(a)')  '**  order parameter of the acyl chains of LIPID **'
       write(ii,'(a)')  '**  NOTE: only DPPC, DMPC, POPC, POPE, DOPC     **'
       write(ii,'(a/)') '**************************************************'
       write(ii,'(a )') '--> Output file is ./Analy/Szz_CG.dat '
     end if

     if(cWhat=='FracGauche') then
       write(ii,'(/a)') '***********************************************'
       write(ii,'(a)')  '**  Gauche ratio in the acyl chains of LIPID **'
       write(ii,'(a)')  '**  NOTE: only DPPC, DMPC, PhPC(DPhPC)!!     **'
       write(ii,'(a/)') '***********************************************'
       write(ii,'(a )') '--> Output files are 1. ./Analy/GaucheFrac.dat           '
       write(ii,'(a )') '-->            ( instantaneous average over acyl chain ) '
       write(ii,'(a )') '-->                  2. ./Analy/Ave_GaucheFrac.dat       '
       write(ii,'(a/)') '-->            ( time average for individual dihedrals ) '
     end if

     if(cWhat=='CoTGLipidC') then
       write(ii,'(/a)') '**************************************************'
       write(ii,'(a)')  '**  Time correlation function of isomerizations **'
       write(ii,'(a)')  '** of dihedrals in the acyl chains of LIPID     **'
       write(ii,'(a)')  '**  NOTE: only DPPC, DMPC, PhPC(DPhPC)!!        **'
       write(ii,'(a/)') '**************************************************'
       write(ii,'(a )') '--> Output files are 1. ./Analy/Cotg_Gamm_Dih_**.dat'
       write(ii,'(a )') '-->                     ./Analy/Cotg_Beta_Dih_**.dat'
       write(ii,'(a/)') '-->                  2. ./Analy/GaucheCheck.dat     '
     end if

     if(cWhat=='PDB') then
       write(ii,'(/a)') '********************************************'
       write(ii,'(a)') '**    Make a PDB file from Trajectories   **'
       write(ii,'(a/)') '********************************************'
       write(ii,'(a/)') '--> Output files are named as ./Analy/snap****.pdb '
     end if

     if(cWhat=='XYZsp') then
       write(ii,'(/a)') '*******************************************************'
       write(ii,'(a)') '**    Make a X-Mol animation file from Trajectories   **'
       write(ii,'(a)') '**    for selected components                         **'
       write(ii,'(a/)') '*******************************************************'
       write(ii,'(a/)') '--> Output files are named as ./Analy/XMolProt.xyz '
     end if

     if(cWhat=='PDBsp') then
       write(ii,'(/a)') '********************************************'
       write(ii,'(a)') '**    Make a PDB file from Trajectories   **'
       write(ii,'(a)') '**    Only selected part                  **'
       write(ii,'(a/)') '********************************************'
       write(ii,'(a/)') '--> Output files are named as ./Analy/snap****.pdb '
     end if

     if(cWhat=='PDBminR') then
       write(ii,'(/a)') '********************************************'
       write(ii,'(a)') '**    Make a PDB file from Trajectories   **'
       write(ii,'(a)') '**    Only selected part                  **'
       write(ii,'(a)') '**    Rotation was determined by the LS   **'
       write(ii,'(a)') '**    fitting of selected amino residue   **'
       write(ii,'(a/)') '********************************************'
       write(ii,'(a/)') '--> Output files are named as ./Analy/snap****.pdb '
     end if

     if(cWhat=='CRD') then
       write(ii,'(/a)') '********************************************'
       write(ii,'(a)') '**    Make a CRD file from Trajectories   **'
       write(ii,'(a/)') '********************************************'
       write(ii,'(a/)') '--> Output files are named as ./Analy/snap****.crd '
     end if

     if(cWhat=='ARC') then
       write(ii,'(/a)') '********************************************'
       write(ii,'(a)') '**    Make a ARC file from Trajectories   **'
       write(ii,'(a)') '********************************************'
       write(ii,'(a/)') '--> Output file is ./Analy/Traj.arc '
     end if

     if(cWhat=='ARC_P') then
       write(ii,'(/a)') '********************************************'
       write(ii,'(a)') '**    Make a ARC file from Trajectories   **'
       write(ii,'(a)') '********************************************'
       write(ii,'(a/)') '--> Output file is ./Analy/TrajProt.arc '
     end if

     if(cWhat=='AMBER_PW') then
       write(ii,'(/a)') '**************************************'
       write(ii,'(a)') '**    Make a AMBER trajectory file  **'
       write(ii,'(a)') '**************************************'
       write(ii,'(a )') '--> Output file is ./Analy/TrajAMBER.trj'
       write(ii,'(a/)') '                   ./Analy/RefProtein.pdb'
     end if

     if(cWhat=='DiffMSD') then
       write(ii,'(/a)') '****************************************'
       write(ii,'(a)')  '** Diffusion coefficient calculation  **'
       write(ii,'(a)')  '**     Mean Square Displacement       **'
       write(ii,'(a/)') '****************************************'
       write(ii,'(a)')  '--> Output files are                           '
       write(ii,'(a)')  '-->     ./Analy/DMSD_(atomname)_(resiname).dat '
     end if

     if(cWhat=='DiffMSDcom') then
       write(ii,'(/a)') '****************************************'
       write(ii,'(a)')  '** Diffusion coefficient calculation  **'
       write(ii,'(a)')  '**     Mean Square Displacement       **'
       write(ii,'(a/)') '****************************************'
       write(ii,'(a)')  '--> Output files are               '
       write(ii,'(a)')  '-->     ./Analy/DMSD_(molname).dat '
     end if

     if(cWhat=='CorrRotW') then
       write(ii,'(/a)') '****************************************'
       write(ii,'(a)')  '** Rotational correlation function    **'
       write(ii,'(a)')  '** relaxation of the water dipole axis**'
       write(ii,'(a/)') '****************************************'
       write(ii,'(a)')  '--> Output files is              '
       write(ii,'(a)')  '-->    ./Analy/CorrRotationW.dat '
     end if

     if(cWhat=='RZ') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  distribution Atoms along Bilayer normal **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a )') '--> Output files are                          '
       write(ii,'(a/)') '-->      ./Analy/RZ_(AtomName)_(ResiName).dat '
     end if

     if(cWhat=='RD') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  Radial distribution of Atoms from the    *'
       write(ii,'(a)')  '**  COM of the simulation box                *'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a )') '--> Output files are                          '
       write(ii,'(a/)') '-->      ./Analy/RD_(AtomName)_(ResiName).dat '
     end if

     if(cWhat=='RDF') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  radial distribution function between    **'
       write(ii,'(a)')  '**          the specified atoms             **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a )') '--> Output files are                          '
       write(ii,'(a/)') &
       & '-->  ./Analy/GR_(AtomName)_(ResiName)--(AtomName)_(ResiName).dat '
       if(QIntraSubt) then
         write(ii,'(a/)') 'Intramolecular contribution is excluded!'
       end if
     end if

     if(cWhat=='GR_LipidGG') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  radial distribution function between    **'
       write(ii,'(a)')  '**  the COMs of lipid molecules (2dimension)**'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a/)') '--> Output file is ./Analy/GR_LipidGG.dat '
     end if

     if(cWhat=='Cavity') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  cavity distribution in the lipid bilayer**'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a )') '--> Output file are ./Analy/distCav****.dat '
       write(ii,'(a/)') '-->                 ./Analy/Cavity.dat      '
     end if

     if(cWhat=='RZCyl') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  distribution Atoms along Bilayer normal **'
       write(ii,'(a)')  '**  inside or outside of defined cylinder   **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a )') '--> Output files are                          '
       write(ii,'(a )') '-->   ./Analy/RZint_(AtomName)_(ResiName).dat '
       write(ii,'(a/)') '-->   ./Analy/RZout_(AtomName)_(ResiName).dat '
     end if

     if(cWhat=='RZRes') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  distribution Atoms along Bilayer normal **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a )') '--> Output files are ./Analy/RZ(Residue Number).dat '
     end if

     if(cWhat=='RZResG') then
       write(ii,'(/a)') '***********************************************'
       write(ii,'(a)')  '**  distribution Atoms along Bilayer normal  **'
       write(ii,'(a)')  '**  origin of the axes is the COM of protein **'
       write(ii,'(a/)') '***********************************************'
       write(ii,'(a )') '--> Output files are ./Analy/RZG(Residue Number).dat '
     end if

     if(cWhat=='RofG_Prot') then
       write(ii,'(/a)') '**************************************************'
       write(ii,'(a)')  '**  Radius of Gyration of the protein molecule  **'
       write(ii,'(a/)') '**************************************************'
       write(ii,'(a )') '--> Output file is ./Analy/RofG_Protein.dat '
     end if

     if(cWhat=='MSD_P') then
       write(ii,'(/a)') '*******************************'
       write(ii,'(a)')  '**  MSD of protein molecule  **'
       write(ii,'(a/)') '*******************************'
       write(ii,'(a )') '--> Output file is ./Analy/MSD_protein.dat '
     end if

     if(cWhat=='MSD_PD') then
       write(ii,'(/a)') '***********************************************'
       write(ii,'(a)')  '**  MSD of Ca atoms of the protein molecule  **'
       write(ii,'(a/)') '***********************************************'
       write(ii,'(a )') '--> Output file is ./Analy/MSD_proteinRD.dat '
     end if

     if(cWhat=='MSD_Pr') then
       write(ii,'(/a)') '*******************************'
       write(ii,'(a)')  '**  MSD of protein molecule  **'
       write(ii,'(a/)') '*******************************'
       write(ii,'(a )') '--> Output file is ./Analy/MSD_proteinR.dat '
     end if

     if(cWhat=='MSD_Pm') then
       write(ii,'(/a)') '*******************************'
       write(ii,'(a)')  '**  MSD of protein molecule  **'
       write(ii,'(a/)') '*******************************'
       write(ii,'(a )') '--> Output file is ./Analy/MSD_Minim.dat '
     end if

     if(cWhat=='MSD_PrD') then
       write(ii,'(/a)') '***********************************************'
       write(ii,'(a)')  '**  MSD of Ca atoms of the protein molecule  **'
       write(ii,'(a/)') '***********************************************'
       write(ii,'(a )') '--> Output file is ./Analy/MSD_proteinRD.dat '
     end if

     if(cWhat=='MSD_PmD') then
       write(ii,'(/a)') '***********************************************'
       write(ii,'(a)')  '**  MSD of Ca atoms of the protein molecule  **'
       write(ii,'(a/)') '***********************************************'
       write(ii,'(a )') '--> Output file is ./Analy/MSD_Minim_detail.dat '
     end if

     if(cWhat=='MSD_PmD_MS') then
       write(ii,'(/a)') '****************************************************************'
       write(ii,'(a)')  '**  Average structure of the protein molecule                 **'
       write(ii,'(a)')  '**  MSD of main and side chain atoms of the protein molecule  **'
       write(ii,'(a/)') '****************************************************************'
       write(ii,'(a )') '--> Output file is ./Analy/ProteinAvConf.dat     '
       write(ii,'(a )') '                   ./Analy/TrajProtCOM.dat       '
       write(ii,'(a )') '                   ./Analy/RMSD_Min_MainChain.dat'
       write(ii,'(a/)') '                   ./Analy/RMSD_Min_SideChain.dat'
     end if

     if(cWhat=='WatOrient') then
       write(ii,'(/a)') '*************************************************'
       write(ii,'(a)')  '**  Water Orientation along the bilayer normal **'
       write(ii,'(a/)') '*************************************************'
       write(ii,'(a )') '--> Output files are ./Analy/DensWater_inCyl.dat  '
       write(ii,'(a )') '-->                  ./Analy/DensWater_outCyl.dat '
       write(ii,'(a )') '-->                  ./Analy/OriWater_inCyl.dat '
       write(ii,'(a/)') '-->                  ./Analy/OriWater_outCyl.dat '
     end if

     if(cWhat=='WatOrientG') then
       write(ii,'(/a)') '*************************************************'
       write(ii,'(a)')  '**  Water Orientation along the bilayer normal **'
       write(ii,'(a/)') '*************************************************'
       write(ii,'(a )') '--> Output files are ./Analy/DensWater_inCylG.dat  '
       write(ii,'(a )') '-->                  ./Analy/DensWater_outCylG.dat '
       write(ii,'(a )') '-->                  ./Analy/OriWater_inCylG.dat '
       write(ii,'(a/)') '-->                  ./Analy/OriWater_outCylG.dat '
     end if

     if(cWhat=='PoreWatTrj') then
       write(ii,'(/a)') '******************************************************'
       write(ii,'(a)')  '**  Trajectories of Water molecules within the Pore **'
       write(ii,'(a/)') '******************************************************'
       write(ii,'(a )') '--> Output files are ./Analy/NumWatPore.dat'
       write(ii,'(a )') '-->                  ./Analy/WatTraj****.dat '
     end if

     if(cWhat=='PoreWater') then
       write(ii,'(/a)') '******************************************************************'
       write(ii,'(a)')  '**  Probability density of Water molecules near the Pore center **'
       write(ii,'(a/)') '******************************************************************'
       write(ii,'(a )') '--> Output files are ./Analy/GridDensityWater.dat'
       write(ii,'(a )') '-->                  ./Analy/WatTraj****.dat '
     end if

     if(cWhat=='ResidTrj') then
       write(ii,'(/a)') '*******************************************'
       write(ii,'(a)')  '**  Trajectories of the selected residue **'
       write(ii,'(a/)') '*******************************************'
       write(ii,'(a )') '--> Output files are ./Analy/ResTraj***.dat '
     end if

     if(cWhat=='AreaOccL') then
       write(ii,'(/a)') '*****************************************************'
       write(ii,'(a)')  '**  Measurement of the order of lipid arrangement  **'
       write(ii,'(a)')  '**  for center of mass of lipid molecule           **'
       write(ii,'(a/)') '*****************************************************'
       write(ii,'(a/)') '--> Output files are ./Analy/OccArea.dat    '
     end if

     if(cWhat=='VoronoiL') then
       write(ii,'(/a)') '*******************************************'
       write(ii,'(a)')  '**  two-dimensional voronoi tesselation  **'
       write(ii,'(a)')  '**  for center of mass of lipid molecule **'
       write(ii,'(a/)') '*******************************************'
       write(ii,'(a )') '--> Output files are ./Analy/Vor_Area.dat    '
       write(ii,'(a )') '                     ./Analy/Vor_Side.dat    '
       write(ii,'(a )') '                     ./Analy/Dis_Side.dat    '
       write(ii,'(a )') '                     ./Analy/PolygonLupp.dat '
       write(ii,'(a )') '                     ./Analy/PolygonLlow.dat '
       write(ii,'(a )') '                     ./Analy/Rg.dat          '
       write(ii,'(a )') '                     ./Analy/Rupp.dat        '
       write(ii,'(a/)') '                     ./Analy/Rlow.dat        '
     end if

     if(cWhat=='VoronoiC') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  two-dimensional voronoi tesselation     **'
       write(ii,'(a)')  '**  for center of mass of lipid alkyl chain **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a )') '--> Output files are ./Analy/Vor_AreaC.dat    '
       write(ii,'(a )') '                     ./Analy/Vor_SideC.dat    '
       write(ii,'(a )') '                     ./Analy/Dis_SideC.dat    '
       write(ii,'(a )') '                     ./Analy/PolygonLuppC.dat '
       write(ii,'(a )') '                     ./Analy/PolygonLlowC.dat '
       write(ii,'(a )') '                     ./Analy/RgC.dat          '
       write(ii,'(a )') '                     ./Analy/RuppC.dat        '
       write(ii,'(a/)') '                     ./Analy/RlowC.dat        '
     end if

     if(cWhat=='LS_anim') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  Making PDB files for animating lipids   **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(/a)') '--> Output files are ./Analy/snap******.pdb   '
       write(ii,'(a/)') '-->                  ./Analy/Rg******         '
     end if

     if(cWhat=='CGdata') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  Analysis for Coarse Graining Model      **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(/a)') '--> Output files are ./Analy/CGconfig.dat     '
       write(ii,'(a )') '-->                  ./Analy/Bond*--*.dat     '
       write(ii,'(a )') '-->                  ./Analy/Angle*--*--*.dat '
       write(ii,'(a/)') '-->                  ./Analy/GR*--*.dat       '
     end if

     if(cWhat=='CGalkane') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  Analysis for Coarse Graining Model      **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(/a)') '--> Output files are ./Analy/CGconfig.dat     '
       write(ii,'(a )') '-->                  ./Analy/Bond*--*.dat     '
       write(ii,'(a )') '-->                  ./Analy/Angle*--*--*.dat '
       write(ii,'(a/)') '-->                  ./Analy/NonBond*--*.dat  '
     end if

     if(cWhat=='DiffW_Z') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  MSD of water along the Z axis           **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a/)') '--> Output files are ./Analy/MSD_W_Z=**.dat'
     end if

     if(cWhat=='LifeTimeW') then
       write(ii,'(/a)') '**********************************************'
       write(ii,'(a)')  '**  Life Time of Coordinating Water         **'
       write(ii,'(a/)') '**********************************************'
       write(ii,'(a/)') '--> Output files are ./Analy/LifeTimeW.dat'
     end if

   else if((SimMethod == 'MM').or.(SimMethod == 'NMA')) then

     write(ii,'(/5x,a)')  '-----------------------------'
     write(ii,'(5x,a,a)') ' Force Field = ',ForceField
     write(ii,'(5x,a/)')  '-----------------------------'

     if(QPBC) then

       write(ii,'(a/)') 'Three-dimensional Periodic Boundary Condition'

     else

       if(QEps_r) then
         write(ii,'(a/)') 'Isolated System  ::  Eps_r = 4R '
       else
         write(ii,'(a/)') 'Isolated System  ::  Eps_r = 1 '
       end if

     end if

     write(ii,'(a)') 'Simulation Condition Details'

     if(QPBC) then

       write(ii,'(5x,a,f10.2,a)') 'Rcutoff = ', sqrt(Rcutoff2),       '[A]'

     end if

     if(QPBC.and.(ForceField/='EAM').and.QCoulomb) then

       if( trim(cCOULOMB) == 'EWALD' ) then
         write(ii,'(/5x,a)')       'Ewald Parameters :: '
         write(ii,'(5x,a,f6.3)')   '     Alpha  = ', Alpha
         write(ii,'(5x,a,i6)')     '     ih2max = ', ih2mx
         write(ii,'(5x,a,i6)')     '     Nh     = ', TmpNh
         write(ii,'(5x,a,3i4)')    '     kMax   = ', kmaxx,kmaxy,kmaxz
       else if( trim(cCOULOMB) == 'PME' ) then
         write(ii,'(/5x,a)')       'PMEwald Parameters :: '
         write(ii,'(5x,a,f6.3)')   '     Alpha    = ', Alpha
         write(ii,'(5x,a,i6)')     '     Grid_x   = ', Nfft(1)
         write(ii,'(5x,a,i6)')     '     Grid_y   = ', Nfft(2)
         write(ii,'(5x,a,i6)')     '     Grid_z   = ', Nfft(3)
         write(ii,'(5x,a,i6)')     '     B-spline = ', Bsp_order
       end if

     end if

   end if

   if(QPathInt) then

     write(ii,'(5x,a)') '************Path Integral*************'
     write(ii,'(5x,a,i6)') 'The number of bead elements, P = ',Nbead
     write(ii,'(5x,a,f6.3,a/)') 'Delta t for bead elements,  dt = ',deltat_ref*1.d+3,'[fs]'

     if(SimMethod == 'CMD') then
       write(ii,'(5x,a)') '************ Centroid MD *************'
       write(ii,'(5x,a,f6.3/)') 'Scaling factor, Gamma_C = ',sqrt(GammaC2)
     end if

   end if

   end do


end subroutine Print_Ini


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_NV(istep)

use CommonBlocks, only : QMaster, QSHAKE, QCorrectCutoff, QOpFix, QAveTh, &
&   ForceField, QThermostat, cThermostatMethod, QRigidBody, Qstdout
use EAM_param, only : Vir_EAM, Ene_EAM
use SHAKEparam, only : Vir_SHAKE, Vir_RATTL
use UnitExParam, only : rprs, cvol
use BathParam
use EwaldParam, only : Vir_Eksp, Ene_Eksp, Ene_Eslf
use OptConstraintParam, only : Vir_OptC, Ene_OptC
use NonbondParam, only : Vir_Ersp, Vir_NBshrt, Vir_NBlong, Ene_LJ, &
&   Ene_Elec, Ene_Ersp, Ene_ELshrt, Ene_ELlong, Ene_NBshrt, Ene_NBlong
use SMDparam, only : Vir_ConstR, Vir_ConstV
use BondedParam, only : &
&   Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro, &
&   Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use CellParam, only : Volume
use TailCorrect, only : Ene_LJ_co, Virial_co
use TimeParam, only : itgn, lk, Timeps
use ThermoData, only : Temp, Temp_Translation, Temp_Rotation, &
&   Ene_kin, Ene_kinT, Ene_kinR, Pkinp, Virial

implicit none

integer :: istep, i
real(8) :: Hamiltonian
real(8), dimension(3,3) :: Pressure
real(8) :: Kinetic, Potential, EneSystem
real(8) :: E_BondMol, E_AngleMol, E_UBMol, E_DihedMol, E_ImproMol
real(8) :: E_LJMol, E_ErMol, E_EkMol, E_EsMol, E_OptC
real(8) :: E_ElecMol, E_Int_Mol, E_Ext_Mol, E_LJcoMol
real(8), dimension(NHchain) :: Kin_ss
real(8), dimension(NHchain) :: Pot_ss
real(8) :: ThermoBath
!# average
real(8) :: fprint
real(8), save :: Av_Ene_Bond, Av_Ene_Angle, Av_Ene_UB, Av_Ene_Dihed
real(8), save :: Av_Ene_Impro, Av_Ene_LJ, Av_Ene_Elec, Av_Ene_Ersp
real(8), save :: Av_Ene_Eksp, Av_Ene_OptC, Av_Temp, Av_Ene_kin
real(8), dimension(3,3), save :: Av_Pkinp, Av_Vir
! ## Thermostat
real(8), save :: Av_ThermoBath
! ## RigidBody
real(8), save :: Av_Temp_Translation, Av_Temp_Rotation

   if((.not.QAveTh).and.(mod(istep,itgn)/=0)) Return

   Ene_LJ = Ene_LJ + Ene_EAM

   if(ForceField(1:2) == 'CG') then
     Ene_Ersp = Ene_ELshrt + Ene_ELlong
     Ene_LJ   = Ene_NBshrt + Ene_NBlong
   end if

!----------------------------------------------------------------------
   call SumEnergy(Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro, &
   &              Ene_LJ, Ene_Elec, Ene_Ersp, Ene_Eksp, Ene_OptC)
!----------------------------------------------------------------------

   if(.not.QSHAKE) then

     Vir_SHAKE = 0.d0
     Vir_RATTL = 0.d0

   end if

   Virial = Vir_Bond  + Vir_Angle + Vir_UB                &
   &      + Vir_Dihed + Vir_Impro + Vir_Ersp + Vir_NBshrt &
   &      + Vir_Eksp  + Vir_OptC  + Vir_EAM  + Vir_NBlong

!-------------------------
   call SumVir( Virial )
!-------------------------

   if(QMaster) then

     Virial = Virial - Vir_SHAKE + Vir_RATTL
     if(QCorrectCutoff) then
       Virial = Virial - Virial_co
     end if
     if(QOpFix) then
       Virial = Virial - Vir_ConstR + Vir_ConstV
     end if

!-------------------------
     call CalcTemp
!-------------------------

     if(QThermostat) then

       if((cThermostatMethod == 'NH').or.(cThermostatMethod == 'NHC')) then

         Kin_ss = 0.5d0 * Mts * Vss * Vss * cvol
         Pot_ss = kT * Rss                * cvol
         Pot_ss(1) = gkT * Rss(1) * cvol

       else if( cThermostatMethod == 'MNHC' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0
         do i = 1, NumMNHC
           Kin_ss(:) = Kin_ss(:) + MMNHC(:,i) * VMNHC(:,i) * VMNHC(:,i)
           Pot_ss(:) = Pot_ss(:) + RMNHC(:,i)
         end do
         Kin_ss = Kin_ss * 0.5d0 * cvol
         Pot_ss = Pot_ss * kT    * cvol

       else if( cThermostatMethod == 'VSCALE' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0

       end if

       ThermoBath = sum( Pot_ss ) + sum( Kin_ss )

     else

       ThermoBath = 0.d0

     end if

     if(QAveTh) then

       if(mod(istep-lk,itgn)==0) then

         Av_Ene_Bond  = 0.d0
         Av_Ene_Angle = 0.d0
         Av_Ene_UB    = 0.d0
         Av_Ene_Dihed = 0.d0
         Av_Ene_Impro = 0.d0
         Av_Ene_LJ    = 0.d0
         Av_Ene_Elec  = 0.d0
         Av_Ene_Ersp  = 0.d0
         Av_Ene_Eksp  = 0.d0
         Av_Ene_OptC  = 0.d0

         Av_Ene_kin   = 0.d0

         Av_Temp  = 0.d0
         Av_Pkinp = 0.d0
         Av_Vir   = 0.d0

         Av_ThermoBath = 0.d0

         Av_Temp_Translation = 0.d0
         Av_Temp_Rotation    = 0.d0

       end if

       Av_Ene_Bond  = Av_Ene_Bond  + Ene_Bond
       Av_Ene_Angle = Av_Ene_Angle + Ene_Angle
       Av_Ene_UB    = Av_Ene_UB    + Ene_UB
       Av_Ene_Dihed = Av_Ene_Dihed + Ene_Dihed
       Av_Ene_Impro = Av_Ene_Impro + Ene_Impro
       Av_Ene_LJ    = Av_Ene_LJ    + Ene_LJ
       Av_Ene_Elec  = Av_Ene_Elec  + Ene_Elec
       Av_Ene_Ersp  = Av_Ene_Ersp  + Ene_Ersp
       Av_Ene_Eksp  = Av_Ene_Eksp  + Ene_Eksp
       Av_Ene_OptC  = Av_Ene_OptC  + Ene_OptC

       Av_Ene_kin   = Av_Ene_kin   + Ene_kin

       Av_Temp  = Av_Temp  + Temp
       Av_Pkinp = Av_Pkinp + Pkinp
       Av_Vir   = Av_Vir   + Virial

       Av_ThermoBath = Av_ThermoBath + ThermoBath

       if(QRigidBody) then
         Av_Temp_Translation = Av_Temp_Translation + Temp_Translation
         Av_Temp_Rotation    = Av_Temp_Rotation    + Temp_Rotation
       end if

     end if

! ---------------------------------------------
     if(mod(istep,itgn)==0) then

       Kinetic    = Ene_kin   * cvol * 0.5d0

       E_BondMol  = Ene_Bond  * cvol
       E_AngleMol = Ene_Angle * cvol
       E_UBMol    = Ene_UB    * cvol
       E_DihedMol = Ene_Dihed * cvol
       E_ImproMol = Ene_Impro * cvol
       E_OptC     = Ene_OptC  * cvol

       E_LJMol    = Ene_LJ    * cvol

       E_ErMol    = Ene_Ersp  * cvol
       E_EkMol    = Ene_Eksp  * cvol
       E_EsMol    = Ene_Eslf  * cvol
       E_ElecMol  = E_ErMol + E_EkMol + E_EsMol

       E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
       &          + E_DihedMol + E_ImproMol + E_OptC
       E_Ext_Mol  = E_LJMol    + E_ElecMol

       if(QCorrectCutoff) then
         E_LJcoMol = Ene_LJ_co * cvol
         E_Ext_Mol = E_Ext_Mol + E_LJcoMol
       end if

       Potential  = E_Int_Mol + E_Ext_Mol

       EneSystem = Kinetic + Potential

       Pressure = ( Pkinp + Virial ) / Volume / rprs * 1.d-6 ! [MPa]

       Hamiltonian = Kinetic + Potential + ThermoBath

       if(QRigidBody) then
#ifdef BMONI
         write(11) Timeps, Temp, Temp_Translation, Temp_Rotation
#else
         write(11,'(f12.4,3f10.2)') Timeps, Temp, Temp_Translation, Temp_Rotation
#endif
         if(Qstdout) then
           write( 6,'(f12.4,3f10.2)') Timeps, Temp, Temp_Translation, Temp_Rotation
         end if
       else
#ifdef BMONI
         write(11) Timeps, Temp
#else
         write(11,'(f12.4,f10.2)') Timeps, Temp
#endif
         if(Qstdout) then
           write( 6,'(f12.4,f10.2)') Timeps, Temp
         end if
       end if

#ifdef BMONI
       write(11) E_BondMol, E_AngleMol, E_UBMol,&
       &         E_DihedMol,E_ImproMol, E_OptC, &
       &         E_LJMol, E_ErMol,              &
       &         E_EkMol, E_EsMol, E_ElecMol
       if(QThermostat) then
         write(11) E_Int_Mol, E_Ext_Mol, Potential, &
         &         Kinetic, EneSystem, ThermoBath
       else
         write(11) E_Int_Mol, E_Ext_Mol, Potential, &
         &         Kinetic, EneSystem
       end if
       write(11) Hamiltonian, Pressure
#else
       write(11,'(6d13.5/5d13.5)') &
       &          E_BondMol, E_AngleMol, E_UBMol,&
       &          E_DihedMol,E_ImproMol, E_OptC, &
       &          E_LJMol, E_ErMol,              &
       &          E_EkMol, E_EsMol, E_ElecMol
       if(QThermostat) then
         write(11,'(6d13.5)') &
         &          E_Int_Mol, E_Ext_Mol, Potential, &
         &          Kinetic, EneSystem, ThermoBath
       else
         write(11,'(5d13.5)') &
         &          E_Int_Mol, E_Ext_Mol, Potential, &
         &          Kinetic, EneSystem
       end if
       write(11,'(d16.8/3(3f12.4/))') &
       &          Hamiltonian,                  &
       &          ( Pressure(1,i) , i = 1 , 3 ),&
       &          ( Pressure(2,i) , i = 1 , 3 ),&
       &          ( Pressure(3,i) , i = 1 , 3 )
#endif
       if(Qstdout) then
         write( 6,'(6d13.5/5d13.5)') &
         &          E_BondMol, E_AngleMol, E_UBMol,&
         &          E_DihedMol,E_ImproMol, E_OptC, &
         &          E_LJMol, E_ErMol,              &
         &          E_EkMol, E_EsMol, E_ElecMol
         if(QThermostat) then
           write( 6,'(6d13.5)') &
           &          E_Int_Mol, E_Ext_Mol, Potential, &
           &          Kinetic, EneSystem, ThermoBath
         else
           write( 6,'(5d13.5)') &
           &          E_Int_Mol, E_Ext_Mol, Potential, &
           &          Kinetic, EneSystem
         end if
         write( 6,'(d16.8/3(3f12.4/))') &
         &          Hamiltonian,                   &
         &          ( Pressure(1,i) , i = 1 , 3 ), &
         &          ( Pressure(2,i) , i = 1 , 3 ), &
         &          ( Pressure(3,i) , i = 1 , 3 )
       end if

       if(QAveTh) then

         fprint = dble(lk) / dble(itgn)

         Av_Ene_Bond  = Av_Ene_Bond  * fprint
         Av_Ene_Angle = Av_Ene_Angle * fprint
         Av_Ene_UB    = Av_Ene_UB    * fprint
         Av_Ene_Dihed = Av_Ene_Dihed * fprint
         Av_Ene_Impro = Av_Ene_Impro * fprint
         Av_Ene_LJ    = Av_Ene_LJ    * fprint
         Av_Ene_Elec  = Av_Ene_Elec  * fprint
         Av_Ene_Ersp  = Av_Ene_Ersp  * fprint
         Av_Ene_Eksp  = Av_Ene_Eksp  * fprint
         Av_Ene_OptC  = Av_Ene_OptC  * fprint

         Av_Ene_kin   = Av_Ene_kin   * fprint

         Av_Temp  = Av_Temp  * fprint
         Av_Pkinp = Av_Pkinp * fprint
         Av_Vir   = Av_Vir   * fprint

         Av_ThermoBath = Av_ThermoBath * fprint

         if(QRigidBody) then
           Av_Temp_Translation = Av_Temp_Translation * fprint
           Av_Temp_Rotation    = Av_Temp_Rotation    * fprint
         end if

         Kinetic    = Av_Ene_kin   * cvol * 0.5d0

         E_BondMol  = Av_Ene_Bond  * cvol
         E_AngleMol = Av_Ene_Angle * cvol
         E_UBMol    = Av_Ene_UB    * cvol
         E_DihedMol = Av_Ene_Dihed * cvol
         E_ImproMol = Av_Ene_Impro * cvol
         E_OptC     = Av_Ene_OptC  * cvol

         E_LJMol    = Av_Ene_LJ    * cvol

         E_ErMol    = Av_Ene_Ersp  * cvol
         E_EkMol    = Av_Ene_Eksp  * cvol
         E_EsMol    = Ene_Eslf  * cvol
         E_ElecMol  = E_ErMol + E_EkMol + E_EsMol

         E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
         &          + E_DihedMol + E_ImproMol + E_OptC
         E_Ext_Mol  = E_LJMol    + E_ElecMol

         if(QCorrectCutoff) then
           E_LJcoMol = Ene_LJ_co * cvol
           E_Ext_Mol = E_Ext_Mol + E_LJcoMol
         end if

         Potential  = E_Int_Mol + E_Ext_Mol

         EneSystem = Kinetic + Potential

         Pressure = ( Av_Pkinp + Av_Vir ) / Volume / rprs * 1.d-6 ! [MPa]

         Hamiltonian = Kinetic + Potential + Av_ThermoBath

#ifdef BMONI
         if(QRigidBody) then
           write(12) Timeps, Av_Temp, Av_Temp_Translation, Av_Temp_Rotation
         else
           write(12) Timeps, Av_Temp
         end if

         write(12) E_BondMol, E_AngleMol, E_UBMol,     &
         &         E_DihedMol,E_ImproMol, E_OptC,      &
         &         E_LJMol, E_ErMol,                   &
         &         E_EkMol, E_EsMol,   E_ElecMol
         if(QThermostat) then
           write(12) E_Int_Mol, E_Ext_Mol, Potential,    &
           &         Kinetic, EneSystem, Av_ThermoBath
         else
           write(12) E_Int_Mol, E_Ext_Mol, Potential,    &
           &         Kinetic, EneSystem
         end if
         write(12) Hamiltonian, Pressure
#else
         if(QRigidBody) then
           write(12,'(f12.4,3f10.2)') &
           &    Timeps, Av_Temp, Av_Temp_Translation, Av_Temp_Rotation
         else
           write(12,'(f12.4,f10.2)') Timeps, Av_Temp
         end if

         write(12,'(6d13.5/5d13.5)') &
         &          E_BondMol, E_AngleMol, E_UBMol,     &
         &          E_DihedMol,E_ImproMol, E_OptC,      &
         &          E_LJMol, E_ErMol,                   &
         &          E_EkMol, E_EsMol,   E_ElecMol
         if(QThermostat) then
           write(12,'(6d13.5)') &
           &          E_Int_Mol, E_Ext_Mol, Potential,    &
           &          Kinetic, EneSystem, Av_ThermoBath
         else
           write(12,'(5d13.5)') &
           &          E_Int_Mol, E_Ext_Mol, Potential,    &
           &          Kinetic, EneSystem
         end if
         write(12,'(d16.8/3(3f12.4/))') &
         &          Hamiltonian,                        &
         &          ( Pressure(1,i) , i = 1 , 3 ),      &
         &          ( Pressure(2,i) , i = 1 , 3 ),      &
         &          ( Pressure(3,i) , i = 1 , 3 )
#endif

       end if

     end if

   end if


end subroutine Print_Energy_NV


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_NP(istep)

use CommonBlocks, only : QMaster, QAveTh, ForceField, QOpFix, QThermostat, &
&   cThermostatMethod, cBarostatMethod, QRigidBody, QCorrectCutoff, Qstdout
use EAM_param, only : Ene_EAM
use UnitExParam, only : rprs, pi, cvol
use BathParam
use EwaldParam, only : Ene_Eksp, Ene_Eslf
use OptConstraintParam, only : Ene_OptC
use NonbondParam, only : Ene_LJ, Ene_Ersp, Ene_ELshrt, Ene_ELlong, &
&   Ene_NBshrt, Ene_NBlong, Ene_Elec
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use SMDparam, only : Vir_ConstR, Vir_ConstV
use CellParam, only : H, Volume
use TailCorrect, only : Ene_LJ_co, Virial_co
use AtomParam, only : Mass
use TimeParam, only : itgn, lk, Timeps
use ThermoData, only : Temp, Temp_Translation, Temp_Rotation, &
&   Ene_kin, Ene_kinT, Ene_kinR, Pkinp, Virial

implicit none

integer :: istep, i
real(8) :: Hamiltonian
real(8), dimension(3,3) :: Pressure
real(8) :: Kinetic, Potential, EneSystem
real(8) :: E_BondMol, E_AngleMol, E_UBMol, E_DihedMol, E_ImproMol
real(8) :: E_LJMol, E_ErMol, E_EkMol, E_EsMol, E_OptC
real(8) :: E_ElecMol, E_Int_Mol, E_Ext_Mol, E_LJcoMol
real(8), dimension(NHchain) :: Kin_ss
real(8), dimension(NHchain) :: Pot_ss
real(8) :: ThermoBath
real(8) :: Pot_p, Kin_p
real(8) :: density, Area, TotalM
real(8), dimension(3) :: VecA, VecB, VecC
real(8) :: LenA, LenB, LenC, csAB, csBC, csCA
real(8) :: AngAB, AngBC, AngCA
! ## Stress >>
real(8), dimension(3,3) :: ElasM, Htrans, G
real(8) :: Pot_elastic
! ## << Stress
!# average
real(8) :: fprint
real(8), save :: Av_Ene_Bond, Av_Ene_Angle, Av_Ene_UB, Av_Ene_Dihed
real(8), save :: Av_Ene_Impro, Av_Ene_LJ, Av_Ene_Elec, Av_Ene_Ersp
real(8), save :: Av_Ene_Eksp, Av_Ene_OptC, Av_Temp, Av_Ene_kin
real(8), dimension(3,3), save :: Av_Pkinp, Av_Vir
! ## Thermostat
real(8), save :: Av_ThermoBath
! ## RigidBody
real(8), save :: Av_Temp_Translation, Av_Temp_Rotation
! ## NPT
real(8), save :: Av_Ene_LJ_co, Av_Volume, Av_Baro_kin
real(8), dimension(3,3), save :: Av_H
! ## ST
real(8), save :: Av_Pot_elastic

   if((.not.QAveTh).and.(mod(istep,itgn)/=0)) Return

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
     call BathTemp
!-------------------------

     if(QOpFix) then
       Virial = Virial - Vir_ConstR + Vir_ConstV
     end if

     if(QThermostat) then

       if((cThermostatMethod == 'NH').or.(cThermostatMethod == 'NHC')) then

         Kin_ss = 0.5d0 * Mts * Vss * Vss * cvol
         Pot_ss = kT * Rss                * cvol
         Pot_ss(1) = gkT * Rss(1) * cvol

       else if( cThermostatMethod == 'MNHC' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0
         do i = 1, NumMNHC
           Kin_ss(:) = Kin_ss(:) + MMNHC(:,i) * VMNHC(:,i) * VMNHC(:,i)
           Pot_ss(:) = Pot_ss(:) + RMNHC(:,i)
         end do
         Kin_ss = Kin_ss * 0.5d0 * cvol
         Pot_ss = Pot_ss * kT    * cvol

       else if( cThermostatMethod == 'VSCALE' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0

       end if

       ThermoBath = sum( Pot_ss ) + sum( Kin_ss )

     else

       ThermoBath = 0.d0

     end if

! ## Stress >>
     if( cBarostatMethod == 'ST' ) then
       Htrans = transpose( H )
       G = matmul(Htrans,H)
       ElasM = matmul(SigmaS,G)
       Pot_elastic = 0.5 * ( ElasM(1,1) + ElasM(2,2) + ElasM(3,3) ) * cvol
     end if
! ## << Stress

     if(QAveTh) then

       if(mod(istep-lk,itgn)==0) then

         Av_Ene_Bond  = 0.d0
         Av_Ene_Angle = 0.d0
         Av_Ene_UB    = 0.d0
         Av_Ene_Dihed = 0.d0
         Av_Ene_Impro = 0.d0
         Av_Ene_LJ    = 0.d0
         Av_Ene_Elec  = 0.d0
         Av_Ene_Ersp  = 0.d0
         Av_Ene_Eksp  = 0.d0
         Av_Ene_OptC  = 0.d0

         Av_Ene_kin   = 0.d0

         Av_Temp  = 0.d0
         Av_Pkinp = 0.d0
         Av_Vir   = 0.d0

         Av_ThermoBath = 0.d0

         Av_Temp_Translation = 0.d0
         Av_Temp_Rotation    = 0.d0

         Av_Ene_LJ_co = 0.d0
         Av_H         = 0.d0
         Av_Volume    = 0.d0
         Av_Baro_kin  = 0.d0

         Av_Pot_elastic = 0.d0

       end if

       Av_Ene_Bond  = Av_Ene_Bond  + Ene_Bond
       Av_Ene_Angle = Av_Ene_Angle + Ene_Angle
       Av_Ene_UB    = Av_Ene_UB    + Ene_UB
       Av_Ene_Dihed = Av_Ene_Dihed + Ene_Dihed
       Av_Ene_Impro = Av_Ene_Impro + Ene_Impro
       Av_Ene_LJ    = Av_Ene_LJ    + Ene_LJ
       Av_Ene_Elec  = Av_Ene_Elec  + Ene_Elec
       Av_Ene_Ersp  = Av_Ene_Ersp  + Ene_Ersp
       Av_Ene_Eksp  = Av_Ene_Eksp  + Ene_Eksp
       Av_Ene_OptC  = Av_Ene_OptC  + Ene_OptC

       Av_Ene_kin   = Av_Ene_kin   + Ene_kin

       Av_Temp  = Av_Temp  + Temp
       Av_Pkinp = Av_Pkinp + Pkinp
       Av_Vir   = Av_Vir   + Virial

       Av_ThermoBath = Av_ThermoBath + ThermoBath

       if(QRigidBody) then
         Av_Temp_Translation = Av_Temp_Translation + Temp_Translation
         Av_Temp_Rotation    = Av_Temp_Rotation    + Temp_Rotation
       end if

       Av_Ene_LJ_co = Av_Ene_LJ_co + Ene_LJ_co
       Av_H         = Av_H         + H
       Av_Volume    = Av_Volume    + Volume
       Av_Baro_kin  = Av_Baro_kin  + Baro_kin

       if( cBarostatMethod == 'ST' ) then
         Av_Pot_elastic = Av_Pot_elastic + Pot_elastic
       end if

     end if

! ---------------------------------------------

     if(mod(istep,itgn)==0) then

       Kinetic    = Ene_kin   * cvol * 0.5d0

       E_BondMol  = Ene_Bond  * cvol
       E_AngleMol = Ene_Angle * cvol
       E_UBMol    = Ene_UB    * cvol
       E_DihedMol = Ene_Dihed * cvol
       E_ImproMol = Ene_Impro * cvol
       E_OptC     = Ene_OptC  * cvol

       E_LJMol    = Ene_LJ    * cvol

       E_ErMol    = Ene_Ersp  * cvol
       E_EkMol    = Ene_Eksp  * cvol
       E_EsMol    = Ene_Eslf  * cvol

       E_ElecMol  = E_ErMol + E_EkMol + E_EsMol

       E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
       &          + E_DihedMol + E_ImproMol + E_OptC
       E_Ext_Mol  = E_LJMol    + E_ElecMol

       if(QCorrectCutoff) then
         E_LJcoMol = Ene_LJ_co * cvol
         E_Ext_Mol = E_Ext_Mol + E_LJcoMol
       end if

       Potential  = E_Int_Mol + E_Ext_Mol

       EneSystem = Kinetic + Potential

       Kin_p = 0.5d0 * Baro_kin    * cvol
       Pot_p = Pressure_o * Volume * cvol

! ## Stress >>
       if( cBarostatMethod == 'ST' ) then
         Pot_p = Pot_p + Pot_elastic
       end if
! ## << Stress

       Pressure = ( Pkinp + Virial ) / Volume / rprs * 1.d-6 ! [MPa]

       TotalM  = Sum(Mass)
       density = ( TotalM * 1.d3 ) / ( Volume * 1.d-24 ) ! [g/cm^3]

       VecA = H(:,1)
       VecB = H(:,2)
       VecC = H(:,3)

       LenA = sqrt( dot_product(VecA,VecA) )
       LenB = sqrt( dot_product(VecB,VecB) )
       LenC = sqrt( dot_product(VecC,VecC) )

       csAB = dot_product(VecA,VecB)/(LenA*LenB)
       csBC = dot_product(VecB,VecC)/(LenB*LenC)
       csCA = dot_product(VecC,VecA)/(LenC*LenA)
       AngAB = acos( csAB ) / pi * 180.d0
       AngBC = acos( csBC ) / pi * 180.d0
       AngCA = acos( csCA ) / pi * 180.d0

       Area = LenA * LenB * sin( acos( csAB ) )

       if(ForceField(1:3) == 'EAM') then
         write(13) sngl(Timeps),H
       end if

       if(QRigidBody) then
#ifdef BMONI
         write(11) Timeps, Temp, Temp_Translation, Temp_Rotation
#else
         write(11,'(f12.4,3f10.2)')    &
         &    Timeps, Temp, Temp_Translation, Temp_Rotation
#endif
         if(Qstdout) then
           write( 6,'(f12.4,3f10.2)')    &
           &    Timeps, Temp, Temp_Translation, Temp_Rotation
         end if
       else
#ifdef BMONI
         write(11) Timeps, Temp
#else
         write(11,'(f12.4,f10.2)') Timeps, Temp
#endif
         if(Qstdout) then
           write( 6,'(f12.4,f10.2)') Timeps, Temp
         end if
       end if

       Hamiltonian = Kinetic + Potential + ThermoBath + Kin_p + Pot_p

#ifdef BMONI
       write(11) E_BondMol,  E_AngleMol, E_UBMol,         &
       &         E_DihedMol, E_ImproMol, E_OptC,          &
       &         E_LJMol, E_ErMol,                        &
       &         E_EkMol, E_EsMol,   E_ElecMol
       if(QThermostat) then
         write(11) E_Int_Mol, E_Ext_Mol, Potential, &
         &         Kinetic, EneSystem, ThermoBath
       else
         write(11) E_Int_Mol, E_Ext_Mol, Potential, &
         &         Kinetic, EneSystem
       end if
       if( cBarostatMethod == 'ST' ) then
         write(11) Pot_p
       end if
       write(11) Hamiltonian, Area, Volume, density,      &
       &         LenA, LenB, LenC, AngBC, AngCA, AngAB
       write(11) H, Pressure
#else
       write(11,'(6d13.5/5d13.5)')                  &
       &          E_BondMol,  E_AngleMol, E_UBMol,         &
       &          E_DihedMol, E_ImproMol, E_OptC,          &
       &          E_LJMol, E_ErMol,                        &
       &          E_EkMol, E_EsMol,   E_ElecMol
       if(QThermostat) then
         write(11,'(6d13.5)')                        &
         &          E_Int_Mol, E_Ext_Mol, Potential, &
         &          Kinetic, EneSystem, ThermoBath
       else
         write(11,'(5d13.5)')                        &
         &          E_Int_Mol, E_Ext_Mol, Potential, &
         &          Kinetic, EneSystem
       end if
       if( cBarostatMethod == 'ST' ) then
         write(11,'(d13.5)') Pot_p
       end if
       write(11,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')      &
       &          Hamiltonian, Area, Volume, density,      &
       &          LenA, LenB, LenC, AngBC, AngCA, AngAB
       write(11,'(3(3f8.3,5x,3f12.4/))')                   &
       &          (H(1,i),i=1,3) , (Pressure(1,i),i=1,3),  &
       &          (H(2,i),i=1,3) , (Pressure(2,i),i=1,3),  &
       &          (H(3,i),i=1,3) , (Pressure(3,i),i=1,3)
#endif

       if(Qstdout) then
         write( 6,'(6d13.5/5d13.5)')                  &
         &          E_BondMol,  E_AngleMol, E_UBMol,         &
         &          E_DihedMol, E_ImproMol, E_OptC,          &
         &          E_LJMol, E_ErMol,                        &
         &          E_EkMol, E_EsMol,   E_ElecMol
         if(QThermostat) then
           write( 6,'(6d13.5)')                        &
           &          E_Int_Mol, E_Ext_Mol, Potential, &
           &          Kinetic, EneSystem, ThermoBath
         else
           write( 6,'(5d13.5)')                        &
           &          E_Int_Mol, E_Ext_Mol, Potential, &
           &          Kinetic, EneSystem
         end if
         if( cBarostatMethod == 'ST' ) then
           write( 6,'(d13.5)') Pot_p
         end if
         write( 6,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')      &
         &          Hamiltonian, Area, Volume, density,      &
         &          LenA, LenB, LenC, AngBC, AngCA, AngAB
         write( 6,'(3(3f8.3,5x,3f12.4/))')                   &
         &          (H(1,i),i=1,3) , (Pressure(1,i),i=1,3),  &
         &          (H(2,i),i=1,3) , (Pressure(2,i),i=1,3),  &
         &          (H(3,i),i=1,3) , (Pressure(3,i),i=1,3)
       end if

       if(QAveTh) then

         fprint = dble(lk) / dble(itgn)

         Av_Ene_Bond  = Av_Ene_Bond  * fprint
         Av_Ene_Angle = Av_Ene_Angle * fprint
         Av_Ene_UB    = Av_Ene_UB    * fprint
         Av_Ene_Dihed = Av_Ene_Dihed * fprint
         Av_Ene_Impro = Av_Ene_Impro * fprint
         Av_Ene_LJ    = Av_Ene_LJ    * fprint
         Av_Ene_Elec  = Av_Ene_Elec  * fprint
         Av_Ene_Ersp  = Av_Ene_Ersp  * fprint
         Av_Ene_Eksp  = Av_Ene_Eksp  * fprint
         Av_Ene_OptC  = Av_Ene_OptC  * fprint

         Av_Ene_kin   = Av_Ene_kin   * fprint

         Av_Temp  = Av_Temp  * fprint
         Av_Pkinp = Av_Pkinp * fprint
         Av_Vir   = Av_Vir   * fprint

         Av_ThermoBath = Av_ThermoBath * fprint

         if(QRigidBody) then
           Av_Temp_Translation = Av_Temp_Translation * fprint
           Av_Temp_Rotation    = Av_Temp_Rotation    * fprint
         end if

         Av_Ene_LJ_co = Av_Ene_LJ_co * fprint
         Av_H         = Av_H         * fprint
         Av_Volume    = Av_Volume    * fprint
         Av_Baro_kin  = Av_Baro_kin  * fprint

         Av_Pot_elastic = Av_Pot_elastic * fprint

         Kinetic    = Av_Ene_kin   * cvol * 0.5d0

         E_BondMol  = Av_Ene_Bond  * cvol
         E_AngleMol = Av_Ene_Angle * cvol
         E_UBMol    = Av_Ene_UB    * cvol
         E_DihedMol = Av_Ene_Dihed * cvol
         E_ImproMol = Av_Ene_Impro * cvol
         E_OptC     = Av_Ene_OptC  * cvol

         E_LJMol    = Av_Ene_LJ    * cvol

         E_ErMol    = Av_Ene_Ersp  * cvol
         E_EkMol    = Av_Ene_Eksp  * cvol
         E_EsMol    = Ene_Eslf  * cvol

         E_ElecMol  = E_ErMol + E_EkMol + E_EsMol

         E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
         &          + E_DihedMol + E_ImproMol + E_OptC
         E_Ext_Mol  = E_LJMol    + E_ElecMol

         if(QCorrectCutoff) then
           E_LJcoMol = Av_Ene_LJ_co * cvol
           E_Ext_Mol = E_Ext_Mol + E_LJcoMol
         end if

         Potential  = E_Int_Mol + E_Ext_Mol

         EneSystem = Kinetic + Potential

         Kin_p = 0.5d0 * Av_Baro_kin    * cvol
         Pot_p = Pressure_o * Av_Volume * cvol

! ## Stress >>
         if( cBarostatMethod == 'ST' ) then
           Pot_p = Pot_p + Av_Pot_elastic
         end if
! ## << Stress

         Pressure = ( Av_Pkinp + Av_Vir ) / Av_Volume / rprs * 1.d-6 ! [MPa]

         density = ( TotalM * 1.d3 ) / ( Av_Volume * 1.d-24 ) ! [g/cm^3]

         VecA = Av_H(:,1)
         VecB = Av_H(:,2)
         VecC = Av_H(:,3)

         LenA = sqrt( dot_product(VecA,VecA) )
         LenB = sqrt( dot_product(VecB,VecB) )
         LenC = sqrt( dot_product(VecC,VecC) )

         csAB = dot_product(VecA,VecB)/(LenA*LenB)
         csBC = dot_product(VecB,VecC)/(LenB*LenC)
         csCA = dot_product(VecC,VecA)/(LenC*LenA)
         AngAB = acos( csAB ) / pi * 180.d0
         AngBC = acos( csBC ) / pi * 180.d0
         AngCA = acos( csCA ) / pi * 180.d0

         Area = LenA * LenB * sin( acos( csAB ) )

#ifdef BMONI
         if(QRigidBody) then
           write(12) Timeps, Av_Temp, Av_Temp_Translation, Av_Temp_Rotation
         else
           write(12) Timeps, Av_Temp
         end if

         Hamiltonian = Kinetic + Potential + Av_ThermoBath + Kin_p + Pot_p

         write(12) E_BondMol,  E_AngleMol, E_UBMol, &
         &         E_DihedMol, E_ImproMol, E_OptC,  &
         &         E_LJMol, E_ErMol,                &
         &         E_EkMol, E_EsMol, E_ElecMol
         if(QThermostat) then
           write(12) E_Int_Mol, E_Ext_Mol, Potential, &
           &         Kinetic, EneSystem, Av_ThermoBath
         else
           write(12) E_Int_Mol, E_Ext_Mol, Potential, &
           &         Kinetic, EneSystem
         end if
         if( cBarostatMethod == 'ST' ) then
           write(12) Pot_p
         end if
         write(12) Hamiltonian, Area, Av_Volume, density,   &
         &         LenA, LenB, LenC, AngBC, AngCA, AngAB
         write(12) Av_H, Pressure
#else
         if(QRigidBody) then
           write(12,'(f12.4,3f10.2)')    &
           &    Timeps, Av_Temp, Av_Temp_Translation, Av_Temp_Rotation
         else
           write(12,'(f12.4,f10.2)') Timeps, Av_Temp
         end if

         Hamiltonian = Kinetic + Potential + Av_ThermoBath + Kin_p + Pot_p

         write(12,'(6d13.5/5d13.5)')                 &
         &          E_BondMol,  E_AngleMol, E_UBMol, &
         &          E_DihedMol, E_ImproMol, E_OptC,  &
         &          E_LJMol, E_ErMol,                &
         &          E_EkMol, E_EsMol, E_ElecMol
         if(QThermostat) then
           write(12,'(6d13.5)')                        &
           &          E_Int_Mol, E_Ext_Mol, Potential, &
           &          Kinetic, EneSystem, Av_ThermoBath
         else
           write(12,'(5d13.5)')                        &
           &          E_Int_Mol, E_Ext_Mol, Potential, &
           &          Kinetic, EneSystem
         end if
         if( cBarostatMethod == 'ST' ) then
           write(12,'(d13.5)') Pot_p
         end if
         write(12,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')      &
         &          Hamiltonian, Area, Av_Volume, density,   &
         &          LenA, LenB, LenC, AngBC, AngCA, AngAB
         write(12,'(3(3f8.3,5x,3f12.4/))')                      &
         &          (Av_H(1,i),i=1,3) , (Pressure(1,i),i=1,3),  &
         &          (Av_H(2,i),i=1,3) , (Pressure(2,i),i=1,3),  &
         &          (Av_H(3,i),i=1,3) , (Pressure(3,i),i=1,3)
#endif
       end if

     end if

   end if

end subroutine Print_Energy_NP


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_iso(istep)

use CommonBlocks, only : QMaster, QRigidBody, QThermostat, Qstdout, &
&   ForceField, QAveTh, cThermostatMethod
use UnitExParam, only : cvol
use BathParam, only : NHchain, NumMNHC, Mts, Rss, Vss, &
&   MMNHC, RMNHC, VMNHC, kT, gkT
use EwaldParam, only : Ene_Eksp
use OptConstraintParam, only : Ene_OptC
use NonbondParam, only : Ene_LJ, Ene_Ersp, Ene_Elec, Ene_ELshrt, Ene_ELlong, &
&   Ene_NBshrt, Ene_NBlong
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use TimeParam, only : itgn, lk, Timeps
use ThermoData, only : Temp, Temp_Translation, Temp_Rotation, &
&   Ene_kin, Ene_kinT, Ene_kinR

implicit none

integer :: istep, i
real(8) :: Hamiltonian
real(8) :: Kinetic, Potential, EneSystem
real(8) :: E_BondMol, E_AngleMol, E_UBMol, E_DihedMol, E_ImproMol
real(8) :: E_LJMol, E_OptC
real(8) :: E_ElecMol, E_Int_Mol, E_Ext_Mol
real(8), dimension(NHchain) :: Kin_ss
real(8), dimension(NHchain) :: Pot_ss
real(8) :: ThermoBath
! ## average
real(8) :: fprint
real(8), save :: Av_Ene_Bond, Av_Ene_Angle, Av_Ene_UB, Av_Ene_Dihed
real(8), save :: Av_Ene_Impro, Av_Ene_LJ, Av_Ene_Elec
real(8), save :: Av_Ene_OptC, Av_Temp, Av_Ene_kin
! ## Thermostat
real(8), save :: Av_ThermoBath
! ## RigidBody
real(8), save :: Av_Temp_Translation, Av_Temp_Rotation

   if((.not.QAveTh).and.(mod(istep,itgn)/=0)) Return

   if(ForceField(1:2) == 'CG') then
     Ene_Elec = Ene_ELshrt + Ene_ELlong
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

     if(QThermostat) then

       if((cThermostatMethod == 'NH').or.(cThermostatMethod == 'NHC')) then

         Kin_ss = 0.5d0 * Mts * Vss * Vss * cvol
         Pot_ss = kT * Rss                * cvol
         Pot_ss(1) = gkT * Rss(1) * cvol

       else if( cThermostatMethod == 'MNHC' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0
         do i = 1, NumMNHC
           Kin_ss(:) = Kin_ss(:) + MMNHC(:,i) * VMNHC(:,i) * VMNHC(:,i)
           Pot_ss(:) = Pot_ss(:) + RMNHC(:,i)
         end do
         Kin_ss = Kin_ss * 0.5d0 * cvol
         Pot_ss = Pot_ss * kT    * cvol

       else if( cThermostatMethod == 'VSCALE' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0

       end if

       ThermoBath = sum( Pot_ss ) + sum( Kin_ss )

     else

       ThermoBath = 0.d0

     end if

     if(QAveTh) then

       if(mod(istep-lk,itgn)==0) then

         Av_Ene_Bond  = 0.d0
         Av_Ene_Angle = 0.d0
         Av_Ene_UB    = 0.d0
         Av_Ene_Dihed = 0.d0
         Av_Ene_Impro = 0.d0
         Av_Ene_LJ    = 0.d0
         Av_Ene_Elec  = 0.d0
         Av_Ene_OptC  = 0.d0

         Av_Ene_kin   = 0.d0

         Av_Temp  = 0.d0

         Av_ThermoBath = 0.d0

         Av_Temp_Translation = 0.d0
         Av_Temp_Rotation    = 0.d0

       end if

       Av_Ene_Bond  = Av_Ene_Bond  + Ene_Bond
       Av_Ene_Angle = Av_Ene_Angle + Ene_Angle
       Av_Ene_UB    = Av_Ene_UB    + Ene_UB
       Av_Ene_Dihed = Av_Ene_Dihed + Ene_Dihed
       Av_Ene_Impro = Av_Ene_Impro + Ene_Impro
       Av_Ene_LJ    = Av_Ene_LJ    + Ene_LJ
       Av_Ene_Elec  = Av_Ene_Elec  + Ene_Elec
       Av_Ene_OptC  = Av_Ene_OptC  + Ene_OptC

       Av_Ene_kin   = Av_Ene_kin   + Ene_kin

       Av_Temp  = Av_Temp  + Temp

       Av_ThermoBath = Av_ThermoBath + ThermoBath

       if(QRigidBody) then
         Av_Temp_Translation = Av_Temp_Translation + Temp_Translation
         Av_Temp_Rotation    = Av_Temp_Rotation    + Temp_Rotation
       end if

     end if

! ---------------------------------------------
     if(mod(istep,itgn)==0) then

       Kinetic    = Ene_kin   * cvol * 0.5d0

       E_BondMol  = Ene_Bond  * cvol
       E_AngleMol = Ene_Angle * cvol
       E_UBMol    = Ene_UB    * cvol
       E_DihedMol = Ene_Dihed * cvol
       E_ImproMol = Ene_Impro * cvol
       E_OptC     = Ene_OptC  * cvol

       E_LJMol    = Ene_LJ    * cvol

       E_ElecMol  = Ene_Elec  * cvol

       E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
       &          + E_DihedMol + E_ImproMol + E_OptC
       E_Ext_Mol  = E_LJMol    + E_ElecMol

       Potential  = E_Int_Mol + E_Ext_Mol

       EneSystem = Kinetic + Potential

       if(QRigidBody) then
#ifdef BMONI
         write(11) Timeps, Temp, Temp_Translation, Temp_Rotation
#else
         write(11,'(f12.4,3f10.2)') Timeps, Temp, Temp_Translation, Temp_Rotation
#endif
         if(Qstdout) then
           write( 6,'(f12.4,3f10.2)') Timeps, Temp, Temp_Translation, Temp_Rotation
         end if
       else
#ifdef BMONI
         write(11) Timeps, Temp
#else
         write(11,'(f12.4,f10.2)') Timeps, Temp
#endif
         if(Qstdout) then
           write( 6,'(f12.4,f10.2)') Timeps, Temp
         end if
       end if

       Hamiltonian = Kinetic + Potential + ThermoBath

#ifdef BMONI
       write(11) E_BondMol,E_AngleMol,E_UBMol,    &
       &         E_DihedMol,E_ImproMol,           &
       &         E_LJMol,E_ElecMol,E_OptC,        &
       &         E_Int_Mol,E_Ext_Mol
       if(QThermostat) then
         write(11)  Potential,Kinetic,EneSystem,     &
         &          ThermoBath,Hamiltonian
       else
         write(11)  Potential,Kinetic,EneSystem,     &
         &          Hamiltonian
       end if
#else
       write(11,'(5d13.5/5d13.5)')   &
       &          E_BondMol,E_AngleMol,E_UBMol,    &
       &          E_DihedMol,E_ImproMol,           &
       &          E_LJMol,E_ElecMol,E_OptC,        &
       &          E_Int_Mol,E_Ext_Mol
       if(QThermostat) then
         write(11,'(4d13.5,d16.8/)')   &
         &          Potential,Kinetic,EneSystem,     &
         &          ThermoBath,Hamiltonian
       else
         write(11,'(3d13.5,d16.8/)')   &
         &          Potential,Kinetic,EneSystem,     &
         &          Hamiltonian
       end if
#endif
       if(Qstdout) then
         write( 6,'(5d13.5/5d13.5)')   &
         &          E_BondMol,E_AngleMol,E_UBMol,    &
         &          E_DihedMol,E_ImproMol,           &
         &          E_LJMol,E_ElecMol,E_OptC,        &
         &          E_Int_Mol,E_Ext_Mol
         if(QThermostat) then
           write( 6,'(4d13.5,d16.8/)')   &
           &          Potential,Kinetic,EneSystem,     &
           &          ThermoBath,Hamiltonian
         else
           write( 6,'(3d13.5,d16.8/)')   &
           &          Potential,Kinetic,EneSystem,     &
           &          Hamiltonian
         end if
       end if

       if(QAveTh) then

         fprint = dble(lk) / dble(itgn)

         Av_Ene_Bond  = Av_Ene_Bond  * fprint
         Av_Ene_Angle = Av_Ene_Angle * fprint
         Av_Ene_UB    = Av_Ene_UB    * fprint
         Av_Ene_Dihed = Av_Ene_Dihed * fprint
         Av_Ene_Impro = Av_Ene_Impro * fprint
         Av_Ene_LJ    = Av_Ene_LJ    * fprint
         Av_Ene_Elec  = Av_Ene_Elec  * fprint
         Av_Ene_OptC  = Av_Ene_OptC  * fprint

         Av_Ene_kin   = Av_Ene_kin   * fprint

         Av_Temp  = Av_Temp  * fprint

         Av_ThermoBath = Av_ThermoBath * fprint

         if(QRigidBody) then
           Av_Temp_Translation = Av_Temp_Translation * fprint
           Av_Temp_Rotation    = Av_Temp_Rotation    * fprint
         end if

         Kinetic    = Av_Ene_kin   * cvol * 0.5d0

         E_BondMol  = Av_Ene_Bond  * cvol
         E_AngleMol = Av_Ene_Angle * cvol
         E_UBMol    = Av_Ene_UB    * cvol
         E_DihedMol = Av_Ene_Dihed * cvol
         E_ImproMol = Av_Ene_Impro * cvol
         E_OptC     = Av_Ene_OptC  * cvol

         E_LJMol    = Av_Ene_LJ    * cvol

         E_ElecMol  = Av_Ene_Elec  * cvol

         E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
         &          + E_DihedMol + E_ImproMol + E_OptC
         E_Ext_Mol  = E_LJMol    + E_ElecMol

         Potential  = E_Int_Mol + E_Ext_Mol

         EneSystem = Kinetic + Potential

#ifdef BMONI
         if(QRigidBody) then
           write(12) Timeps, Av_Temp, Av_Temp_Translation, Av_Temp_Rotation
         else
           write(12) Timeps, Av_Temp
         end if

         Hamiltonian = Kinetic + Potential + Av_ThermoBath

         write(12)  E_BondMol,E_AngleMol,E_UBMol,    &
         &          E_DihedMol,E_ImproMol,           &
         &          E_LJMol,E_ElecMol,E_OptC,        &
         &          E_Int_Mol,E_Ext_Mol
         if(QThermostat) then
           write(12)  Potential,Kinetic,EneSystem,     &
           &          Av_ThermoBath,Hamiltonian
         else
           write(12)  Potential,Kinetic,EneSystem,     &
           &          Hamiltonian
         end if
#else
         if(QRigidBody) then
           write(12,'(f12.4,3f10.2)') Timeps, Av_Temp, Av_Temp_Translation, Av_Temp_Rotation
         else
           write(12,'(f12.4,f10.2)') Timeps, Av_Temp
         end if

         Hamiltonian = Kinetic + Potential + Av_ThermoBath

         write(12,'(5d13.5/5d13.5)')   &
         &          E_BondMol,E_AngleMol,E_UBMol,    &
         &          E_DihedMol,E_ImproMol,           &
         &          E_LJMol,E_ElecMol,E_OptC,        &
         &          E_Int_Mol,E_Ext_Mol
         if(QThermostat) then
           write(12,'(4d13.5,d16.8/)')   &
           &          Potential,Kinetic,EneSystem,     &
           &          Av_ThermoBath,Hamiltonian
         else
           write(12,'(3d13.5,d16.8/)')   &
           &          Potential,Kinetic,EneSystem,     &
           &          Hamiltonian
         end if
#endif
       end if

     end if

   end if

end subroutine Print_Energy_iso


!######################################################################
!######################################################################


! *******************************
! * Instantaneous Configuration *
! *******************************
!
subroutine Print_Config

use Numbers, only : N
use IOparam, only : DirectoryName, trajectory_file, ItrjF, NtrjF, StoreTrj
use Configuration, only : R
use CommonBlocks, only : Job_name, QMaster, iTrjForm, QPBC
use CellParam, only : H
use AtomParam, only : AtomName
use TimeParam, only : Timeps

implicit none

character(len=4) :: Aname
integer :: i, nstr
integer, dimension(20) :: icntrl
real(8), dimension(6) :: xcell
real, dimension(N) :: X, Y, Z

   ItrjF = ItrjF + 1

   if(iTrjForm==1) then

     if(QMaster) then
       if(ItrjF==1) then
         open(21,file=trim(DirectoryName)//trim(trajectory_file),&
         &    form='formatted',status='unknown')
       end if
       write(21,'(i9)') N
       if(QPBC) then
         write(21,'(10f9.4)') Timeps,H
       else
         write(21,'(f9.4)') Timeps
       end if
       do i = 1 , N
         Aname = AtomName(i)
         write(21,'(a,3f9.3)') Aname(1:1),R(:,i)
       end do
     end if

     if(ItrjF==StoreTrj) then

       call StoreRestart
       NtrjF = NtrjF + 1
       if(NtrjF<10) then
         write(trajectory_file,'(a,a,i1,a)') trim(adjustl(Job_name)),'.r00',NtrjF,'.xyz'
       else if(NtrjF<100) then
         write(trajectory_file,'(a,a,i2,a)') trim(adjustl(Job_name)),'.r0',NtrjF,'.xyz'
       else
         write(trajectory_file,'(a,a,i3,a)') trim(adjustl(Job_name)),'.r',NtrjF,'.xyz'
       end if

       if(QMaster) close(21)

       ItrjF = 0

     end if

   else if(iTrjForm==2) then

     if(QMaster) then
       if(ItrjF==1) then
         open(21,file=trim(DirectoryName)//trim(trajectory_file),form='unformatted',&
         &    status='unknown')
       end if
       write(21) sngl(Timeps), N
       if(QPBC) write(21) sngl(H)
       write(21) sngl(R)
     end if

     if(ItrjF==StoreTrj) then
       call StoreRestart
       NtrjF = NtrjF + 1
       if(NtrjF<10) then
         write(trajectory_file,'(a,a,i1)') trim(adjustl(Job_name)),'.r00',NtrjF
       else if(NtrjF<100) then
         write(trajectory_file,'(a,a,i2)') trim(adjustl(Job_name)),'.r0',NtrjF
       else
         write(trajectory_file,'(a,a,i3)') trim(adjustl(Job_name)),'.r',NtrjF
       end if
       if(QMaster) close(21)
       ItrjF = 0
     end if

   else if(iTrjForm==3) then

     if(QMaster) then
       if(ItrjF==1) then
         open(21,file=trim(DirectoryName)//trim(trajectory_file),form='unformatted',&
         &    status='unknown')
         Aname = 'CORD'
         icntrl = 0
         nstr = 0
         icntrl(1) = StoreTrj  ! number of frames
         icntrl(2) = 1         ! number of steps in previous run
         icntrl(3) = 1         ! frequency of saving
         icntrl(4) = StoreTrj  ! total number of steps
         icntrl(8) = N*3       ! number of degrees of freedm
         icntrl(10) = 981668463 ! coded time step
         icntrl(11) = 1         ! coded crystallographic group (or zro)
         icntrl(20) = 27        ! CHARMM version number
         write(21) Aname,icntrl
         write(21) nstr
         write(21) N
       end if

       xcell(1) = H(1,1)
       xcell(2) = 0.d0
       xcell(3) = H(2,2)
       xcell(4) = 0.d0
       xcell(5) = 0.d0
       xcell(6) = H(3,3)

       do i = 1, N
         X(i) = sngl(R(1,i))
         Y(i) = sngl(R(2,i))
         Z(i) = sngl(R(3,i))
       end do

       write(21) xcell
       write(21) (X(i),i=1,N)
       write(21) (Y(i),i=1,N)
       write(21) (Z(i),i=1,N)

     end if

     if(ItrjF==StoreTrj) then
       call StoreRestart
       NtrjF = NtrjF + 1
       if(NtrjF<10) then
         write(trajectory_file,'(a,a,i1,a)') trim(adjustl(Job_name)),'.r00',NtrjF,'.dcd'
       else if(NtrjF<100) then
         write(trajectory_file,'(a,a,i2,a)') trim(adjustl(Job_name)),'.r0',NtrjF,'.dcd'
       else
         write(trajectory_file,'(a,a,i3,a)') trim(adjustl(Job_name)),'.r',NtrjF,'.dcd'
       end if
       if(QMaster) close(21)
       ItrjF = 0
     end if

   end if

end subroutine Print_Config


!######################################################################
!######################################################################


! *******************************
! * Instantaneous Configuration *
! *******************************
!
subroutine Print_Config_PI

use Numbers, only : N
use IOparam, only : DirectoryName, trajectory_file, ItrjF, NtrjF, StoreTrj
use CommonBlocks, only : Job_name, QMaster, iTrjForm, QPBC
use CommonPI
use CellParam, only : H
use AtomParam, only : AtomName
use TimeParam, only : Timeps

implicit none

character(len=4) :: Aname
integer :: i, j

   ItrjF = ItrjF + 1

   if(iTrjForm==1) then

     if(QMaster) then

       if(ItrjF==1) then

         open(21,file=trim(DirectoryName)//trim(trajectory_file),&
         &    form='formatted',status='unknown')

       end if

       write(21,'(i9)') N * Nbead
       write(21,'(a,f10.4)') 'time = ',Timeps

       do i = 1 , N
       do j = 1 , Nbead

         Aname = AtomName(i)
         write(21,'(a,3f8.3)') Aname(1:1),Rpi(:,i,j)

       end do
       end do

     end if

     if(ItrjF==StoreTrj) then

       call StoreRestart

       NtrjF = NtrjF + 1

       if(NtrjF<10) then

         write(trajectory_file,'(a,a,i1,a)') trim(adjustl(Job_name)),'.r00',NtrjF,'.xyz'

       else if(NtrjF<100) then

         write(trajectory_file,'(a,a,i2,a)') trim(adjustl(Job_name)),'.r0',NtrjF,'.xyz'

       else

         write(trajectory_file,'(a,a,i3,a)') trim(adjustl(Job_name)),'.r',NtrjF,'.xyz'

       end if

       if(QMaster) close(21)

       ItrjF = 0

     end if

   else

     if(QMaster) then

       if(ItrjF==1) then

         open(21,file=trim(DirectoryName)//trim(trajectory_file),&
         &    form='unformatted',status='unknown')

       end if

       write(21) sngl(Timeps), N, Nbead
       if(QPBC) write(21) sngl(H)
       write(21) sngl(Rpi)

     end if

     if(ItrjF==StoreTrj) then

       call StoreRestart

       NtrjF = NtrjF + 1

       if(NtrjF<10) then

         write(trajectory_file,'(a,a,i1)') trim(adjustl(Job_name)),'.r00',NtrjF

       else if(NtrjF<100) then

         write(trajectory_file,'(a,a,i2)') trim(adjustl(Job_name)),'.r0',NtrjF

       else

         write(trajectory_file,'(a,a,i3)') trim(adjustl(Job_name)),'.r',NtrjF

       end if

       if(QMaster) close(21)

       ItrjF = 0

     end if

   end if

end subroutine Print_Config_PI


!######################################################################
!######################################################################


! *******************************
! *   Instantaneous Velocity    *
! *******************************
!
subroutine Print_Velocity

use Numbers, only : N
use IOparam, only : DirectoryName, velocity_file, NvelF, IvelF, StoreVel
use Configuration, only : Vel
use CommonBlocks, only : Job_name, QMaster, QPBC, QRigidBody
use RBparam
use CellParam, only : H
use TimeParam, only : Timeps

implicit none

   IvelF = IvelF + 1

   if(QMaster) then

     if(IvelF==1) then

       open(31,file=trim(DirectoryName)//trim(velocity_file),&
       &    form='unformatted',status='unknown')

     end if

     write(31) sngl(Timeps), N
     if(QPBC) write(31) sngl(H)

     if(QRigidBody) then
       write(31) sngl(V_RB)
       write(31) sngl(Lmoment)
     else
       write(31) sngl(Vel)
     end if

   end if

   if(IvelF==StoreVel) then

     NvelF = NvelF + 1

     if(NvelF<10) then

       write(velocity_file,'(a,a,i1)') trim(adjustl(Job_name)),'.v00',NvelF

     else if(NvelF<100) then

       write(velocity_file,'(a,a,i2)') trim(adjustl(Job_name)),'.v0',NvelF

     else

       write(velocity_file,'(a,a,i3)') trim(adjustl(Job_name)),'.v',NvelF

     end if

     if(QMaster) close(31)

     IvelF = 0

   end if

end subroutine Print_Velocity


!######################################################################
!######################################################################


! *******************************
! *   Instantaneous Velocity    *
! *******************************
!
subroutine Print_Velocity_PI

use Numbers, only : N
use IOparam, only : DirectoryName, velocity_file, NvelF, IvelF, StoreVel
use CommonBlocks, only : Job_name, QMaster, QPBC
use CommonPI
use CellParam, only : H
use TimeParam, only : Timeps

implicit none

! ## for parallel
   call VnmPI

   IvelF = IvelF + 1

   if(QMaster) then

     if(IvelF==1) then

       open(31,file=trim(DirectoryName)//trim(velocity_file),&
       &    form='unformatted',status='unknown')

     end if

     write(31) sngl(Timeps), N, Nbead
     if(QPBC) write(31) sngl(H)
     write(31) sngl(Vnm)

   end if

   if(IvelF==StoreVel) then

     NvelF = NvelF + 1

     if(NvelF<10) then

       write(velocity_file,'(a,a,i1)') trim(adjustl(Job_name)),'.v00',NvelF

     else if(NvelF<100) then

       write(velocity_file,'(a,a,i2)') trim(adjustl(Job_name)),'.v0',NvelF

     else

       write(velocity_file,'(a,a,i3)') trim(adjustl(Job_name)),'.v',NvelF

     end if

     if(QMaster) close(31)

     IvelF = 0

   end if

end subroutine Print_Velocity_PI


!######################################################################
!######################################################################


! ##################################################
! ##      calculate thermodynamic quantities      ##
! ##################################################

subroutine DPD_monitor(istep)

use Numbers, only : N, NumSpec, NumMol, NumAtm
use CommonBlocks, only : QHfunc, Qmixord
use Configuration
use CommonDPD
use RBparam
use CellParam, only : Volume, CellL
use TimeParam, only : itgn, ixc, ixv
use ThermoData, only : Temp, Temp_Translation, Temp_Rotation, &
&   Ene_kin, Pkinp, Virial

implicit none

integer :: istep, i, j, k, l
integer :: Num
real(8) :: TotalE, PressInst
real(8), save :: Pressure
real(8), save :: Vir
real(8), dimension(3), save :: Temperature
real(8), dimension(3), save :: Energy
real(8), save :: Pxy, Pr_xy  ! for sheared system 
real(8), dimension(3,3) :: Ptensor
integer :: iy
real(8) :: Rcom, Vcom
real(8), save :: VirC, VirD, VirR
integer :: iix, iiy, iiz

! ## for taking averages during the simulation

   if( istep == 1 ) then

     Pressure    = 0.d0
     Temperature = 0.d0
     Energy      = 0.d0
     Vir         = 0.d0

     VirC        = 0.d0
     VirD        = 0.d0
     VirR        = 0.d0

     if(QSheared) then
       Pxy = 0.d0
       Pr_xy = 0.d0
     end if

!##  for a calculation of a velocity profile along the y axis

     Nslab = nint(CellL(2))

     allocate( AveVslab( Nslab ) )

     AveVslab = 0.d0
     ThSlab = CellL(2) / Nslab

     ErrIter = 0.d0

     if(QHfunc) then
       allocate( Vlist(3,-50:50) )
       Vlist = 0
     end if

     if(Qmixord) sf1 = 0.d0

   end if

! ## Temperature
! ------------------
   call CalcTempDP
! ------------------

! ##  for colloids (rigid-body)
! ##### Virial correction due to Rigid-Body constrain (internal term)

   if(QColloid) then

     Num = 0

     do i = 1 , NumSpec

       if(ColloidFlag(i)) then

         do k = 1 , NumColloid

           do j = 1 , NumCollAtm

             Num = Num + 1

             Virial(1,1) = Virial(1,1) - 2.d0 * Rmolec(1,j,k) * FrcDP(1,Num)
             Virial(1,2) = Virial(1,2) - 2.d0 * Rmolec(1,j,k) * FrcDP(2,Num)
             Virial(1,3) = Virial(1,3) - 2.d0 * Rmolec(1,j,k) * FrcDP(3,Num)
             Virial(2,1) = Virial(2,1) - 2.d0 * Rmolec(2,j,k) * FrcDP(1,Num)
             Virial(2,2) = Virial(2,2) - 2.d0 * Rmolec(2,j,k) * FrcDP(2,Num)
             Virial(2,3) = Virial(2,3) - 2.d0 * Rmolec(2,j,k) * FrcDP(3,Num)
             Virial(3,1) = Virial(3,1) - 2.d0 * Rmolec(3,j,k) * FrcDP(1,Num)
             Virial(3,2) = Virial(3,2) - 2.d0 * Rmolec(3,j,k) * FrcDP(2,Num)
             Virial(3,3) = Virial(3,3) - 2.d0 * Rmolec(3,j,k) * FrcDP(3,Num)

           end do

         end do

       else

         Num = Num + NumMol(i) * NumAtm(i)

       end if

     end do

   end if

! ## for check the pressure tensor
   Ptensor(:,:) = ( Pkinp(:,:) + Virial(:,:) ) / Volume

#ifdef DPDcheck
   write(10) Ptensor
#endif

! ## pressure

   PressInst = 0.d0

   do i = 1 , 3

     PressInst = PressInst + Ptensor(i,i)

   end do

   PressInst = PressInst / 3.d0

! ## H-function

   if(QHfunc) then
     do i = 1, N
       iix = nint(Vel(1,i) * 10.d0)
       iiy = nint(Vel(2,i) * 10.d0)
       iiz = nint(Vel(3,i) * 10.d0)
       Vlist(1,iix) = Vlist(1,iix) + 1
       Vlist(2,iiy) = Vlist(2,iiy) + 1
       Vlist(3,iiz) = Vlist(3,iiz) + 1
     end do
   end if

! ## mixing order parameter
   if(Qmixord) then
     j = 0
     do i = 1, N/2
       if(R(1,i)<0.d0) j = j + 1
     end do
     sf1 = sf1 + (dble(j)-N/4.d0)*4.d0/dble(N)
   end if

!##  for a calculation of a velocity profile along the y axis

   if(istep > Iequil) then

   l = 0
   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       do j = 1, NumColloid

         iy = int(( R_RB(2,j) + CellL(2)*0.5d0 ) / ThSlab) + 1
         AveVslab(iy) = AveVslab(iy) + V_RB(1,j)

       end do

       l = l + NumColloid * NumCollAtm

     end if

     do j = 1 , NumMol(i)

       Rcom = 0.d0
       Vcom = 0.d0

       do k = 1 , NumAtm(i)

         l = l + 1
         Rcom = Rcom + R(2,l)
         Vcom = Vcom + Vel(1,l)

       end do

       Rcom = Rcom / dble(NumAtm(i))

       iy = int(( Rcom + CellL(2)*0.5d0 ) / ThSlab) + 1
       AveVslab(iy) = AveVslab(iy) + Vcom

     end do

   end do

   end if

! ## for averaging

   Pressure = Pressure + PressInst
   Vir      = Vir      + ( Virial(1,1) + Virial(2,2) + Virial(3,3) ) / (3.d0*Volume)

   VirC = VirC + (VirialC(1,1) + VirialC(2,2) + VirialC(3,3)) / (3.d0*Volume)
   VirD = VirD + (VirialD(1,1) + VirialD(2,2) + VirialD(3,3)) / (3.d0*Volume)
   VirR = VirR + (VirialR(1,1) + VirialR(2,2) + VirialR(3,3)) / (3.d0*Volume)

   Temperature(1) = Temperature(1) + Temp

   if(QSheared) then

     Pxy = Pxy + (Pkinp(1,2) + Virial(1,2)) / Volume
! ## check
     Pr_xy = Pr_xy + VirialR(1,2) / Volume

   end if

   if(QColloid) then

     Temperature(2) = Temperature(2) + Temp_Translation
     Temperature(3) = Temperature(3) + Temp_Rotation

   end if

   Ene_kin = Ene_kin * 0.5

   TotalE = Ene_kin + PotDP

   Energy(1) = Energy(1) + TotalE  ! / dble(N)
   Energy(2) = Energy(2) + PotDP   ! / dble(N)
   Energy(3) = Energy(3) + Ene_kin ! / dble(N)

   SdR = SdR + dR
   if(QColloid) then
     SdRc = SdRc + dRc
   end if

   if(mod(istep,itgn)==0) then

     call Print_Energy_DP(Pressure,Vir,VirC,VirD,VirR,Temperature,Energy,Pxy,Pr_xy,istep)

     if((istep>Iequil).and.(NumSpec==1).and.(NumAtm(1)==1)) call RDF

   end if

   if(mod(istep,ixc)==0) then

     call Print_Config

   end if

   if(mod(istep,ixv)==0) then

     call Print_Velocity

   end if

   if(((IntegrMethod=='VVerletSC').or.(IntegrMethod=='LeapFrogSC').or.&
       (IntegrMethod=='PetersSC')).and.(.not.QColloid)) then
     do i = 1 , N
       ErrIter = ErrIter + dabs( VelP(1,i) - Vel(1,i) ) &
       &                 + dabs( VelP(2,i) - Vel(2,i) ) &
       &                 + dabs( VelP(3,i) - Vel(3,i) )
     end do
   end if

end subroutine DPD_monitor


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_DP(Pressure,Vir,VirC,VirD,VirR,Temperature,Energy,Pxy,Pr_xy,istep)

use Numbers, only : N, NumSpec, NumMol, NumAtm
use CommonBlocks, only : QHfunc, Qstdout, Qmixord
use CommonDPD
use TimeParam, only : itgn, Timeps

implicit none

real(8) :: Pressure
real(8) :: Vir
real(8), dimension(3) :: Temperature
real(8), dimension(3) :: Energy
real(8), dimension(NumSpec) :: MSD
real(8) :: fprint, Pxy, Pr_xy
integer :: i, j, k, l, istep
real(8) :: VirC, VirD, VirR
real(8) :: ffx, ffy, ffz, Hx, Hy, Hz

   fprint = 1.d0 / dble(itgn)

   if(istep == itgn) then

     Pressure_av    = 0.d0
     Temperature_av = 0.d0
     Energy_av      = 0.d0
     Vir_av         = 0.d0

     VirC_av = 0.d0
     VirR_av = 0.d0
     VirD_av = 0.d0

     if(QSheared) then
       Pxy_av = 0.d0
       Pr_xy_av = 0.d0
     end if

   end if

   if(istep > Iequil) then

     Pressure_av    = Pressure_av    + Pressure
     Vir_av         = Vir_av         + Vir
     Energy_av      = Energy_av      + Energy
     Temperature_av = Temperature_av + Temperature

     VirC_av = VirC_av + VirC
     VirR_av = VirR_av + VirR
     VirD_av = VirD_av + VirD

     if(QSheared) then
       Pxy_av = Pxy_av + Pxy
       Pr_xy_av = Pr_xy_av + Pr_xy
     end if

   end if

   Pressure    = Pressure    * fprint
   Temperature = Temperature * fprint
   Energy      = Energy      * fprint
   Vir         = Vir         * fprint

   VirC = VirC * fprint
   VirR = VirR * fprint
   VirD = VirD * fprint

   if(QSheared) then
     Pxy  = Pxy  * fprint
     Pr_xy  = Pr_xy  * fprint
   end if

! ## H-function 
   if(QHfunc) then

     Hx = 0.d0
     Hy = 0.d0
     Hz = 0.d0

     do i = -50, 50
       if(Vlist(1,i)/=0) then
         ffx = Vlist(1,i) * fprint / dble(N)
         Hx  = Hx + ffx * log(ffx) * 0.1d0
       end if
       if(Vlist(2,i)/=0) then
         ffy = Vlist(2,i) * fprint / dble(N)
         Hy  = Hy + ffy * log(ffy) * 0.1d0
       end if
       if(Vlist(3,i)/=0) then
         ffz = Vlist(3,i) * fprint / dble(N)
         Hz  = Hz + ffz * log(ffz) * 0.1d0
       end if
     end do

   end if

   if(Qmixord) then

     sf1 = sf1 * fprint

   end if

! ## IO 

   if(Qstdout) then
   write(6,'(f15.3,3f10.3/4d15.7)') timeps, Temperature, Energy,Pressure
   end if

   write(28,'(f15.3,3f15.7)') timeps, Temperature
   write(22,'(f15.3,3e15.7)') timeps, Energy
   if(QSheared) then
     write(23,'(f15.3,4e15.7)') timeps, Pressure, Vir, Pxy, Pr_xy
   else
     write(23,'(f15.3,2e15.7)') timeps, Pressure, Vir
   end if

   write(25,'(f15.3,3e15.7)') timeps, VirC, VirR, VirD

   if(QHfunc) then
     write(29,'(f15.3,3e15.7)') timeps, Hx, Hy, Hz
     print *, Hx, Hy, Hz
   end if

   if(Qmixord) then
     write(30,'(f15.3,e15.7)') timeps, sf1
     print *, sf1
   end if

! ## memory refreshment for a new average

   Pressure    = 0.d0
   Temperature = 0.d0
   Energy      = 0.d0
   Vir         = 0.d0

   VirC = 0.d0
   VirR = 0.d0
   VirD = 0.d0

   MSD(:) = 0.d0

   Pxy = 0.d0
   Pr_xy = 0.d0

   if(QHfunc) Vlist = 0

! ## Meansquare displacement

   l = 0

   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       do j = 1 , NumColloid

         MSD(i) = MSD(i) + dot_product(SdRc(:,j),SdRc(:,j))

       end do

       l = l + NumColloid*NumCollAtm

     else

       do j = 1 , NumMol(i)

         do k = 1 , NumAtm(i)

           l = l + 1
           MSD(i) = MSD(i) + dot_product(SdR(:,l),SdR(:,l))

         end do

       end do

     end if

   end do

   do i = 1, NumSpec

     if(ColloidFlag(i)) then
       MSD(i) = MSD(i) / dble(NumColloid)
     else
       MSD(i) = MSD(i) / dble(NumMol(i)*NumAtm(i))
     end if

   end do

   write(24,'(f15.3,10e15.7)') timeps, (MSD(i), i=1,NumSpec)

end subroutine Print_Energy_DP


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Average_DP

use Numbers, only : N, NumSpec, NumMol, NumAtm
use IOparam, only : DirectoryName
use CommonBlocks, only : Job_name
use CommonDPD
use UnitExParam, only : pi
use CellParam, only : Volume, CellL
use TimeParam, only : Nstep, itgn

implicit none

character(len=80) :: String
integer :: i, j, Nfree
real(8) :: yy
real(8) :: drr, qk, gr

  if(QColloid) then

    Nfree = 0
    do i = 1 , NumSpec

      if(ColloidFlag(i)) then
        Nfree = Nfree + NumColloid
      else
        Nfree = Nfree + NumMol(i)*NumAtm(i)
      end if

    end do

   else

     Nfree = N

   end if

   write(String,'(2a)') trim(Job_Name),'.Average'

   open(7,file=trim(DirectoryName)//trim(String),status='unknown')

   if(Nstep > Iequil) then

     Pressure_av    = Pressure_av    / dble(Nstep-Iequil)
     Vir_av         = Vir_av         / dble(Nstep-Iequil)
     Energy_av      = Energy_av      / dble(Nstep-Iequil)
     Temperature_av = Temperature_av / dble(Nstep-Iequil)

     VirC_av = VirC_av / dble(Nstep-Iequil)
     VirR_av = VirR_av / dble(Nstep-Iequil)
     VirD_av = VirD_av / dble(Nstep-Iequil)

     if(QSheared) then
       Pxy_av = Pxy_av / dble(Nstep-Iequil)
       Pr_xy_av = Pr_xy_av / dble(Nstep-Iequil)
     end if

     do i = 6, 7

       write(i,'(a,e15.7)') ' Pressure     : ',Pressure_av
       write(i,'(a,e15.7)') ' Virial       : ',Vir_av

       write(i,'(a,e15.7)') ' Virial_C     : ',VirC_av
       write(i,'(a,e15.7)') ' Virial_R     : ',VirR_av
       write(i,'(a,e15.7)') ' Virial_D     : ',VirD_av

       write(i,'(a,e15.7)') ' Temperature  : ',Temperature_av(1)
       write(i,'(a,e15.7)') '    (trans.)  : ',Temperature_av(2)
       write(i,'(a,e15.7)') '    (rotat.)  : ',Temperature_av(3)
       write(i,'(a      )') ' Energy '
       write(i,'(a,e15.7)') '    (Total )  : ',Energy_av(1)
       write(i,'(a,e15.7)') '    (Poten.)  : ',Energy_av(2)
       write(i,'(a,e15.7)') '    (Kinet.)  : ',Energy_av(3)

       if(QSheared) then
         write(i,'(/a,e15.7)') ' Pxy          : ',Pxy_av
#ifdef DPDcheck
         write(i,'( a,e15.7)') ' Pxy (random) : ',Pr_xy_av
#endif
         write(i,'(/a)') ' Velocity Profile along the y-axis '

         do j = 1 , Nslab

           yy = (j - 0.5) * ThSlab - CellL(2)*0.5d0
           write(i,*) yy, AveVslab(j) / dble((Nstep-Iequil)*Nfree)

         end do

       end if

     end do

   end if

   if((IntegrMethod=='VVerletSC').or.(IntegrMethod=='LeapFrogSC')) then
     do i = 6, 7
       write(i,'(a,e15.7)') '  Iteration error  :  ',ErrIter/dble(Nstep * N)
     end do
   end if

   close(7)

   if(Nstep > Iequil) then

     if((NumSpec==1).and.(NumAtm(1)==1)) then

       write(String,'(2a)') trim(Job_Name),'.gr'

       open(7,file=trim(DirectoryName)//trim(String),status='unknown')

       do i = 1, irc

         drr = i * 0.05d0
         qk = Volume / (N*4.d0*pi*drr*drr*0.05)
         gr = irdf(i) * 2.d0 * qk / dble(N*(Nstep-Iequil)/itgn)

         write(7,'(f7.2,f12.4)') drr,gr

       end do

       close(7)

     end if

   end if

! ----------------------
   call Write_ConfigDP
! ----------------------

   call Print_XYZ

end subroutine Print_Average_DP


!#####################################################################
!#####################################################################


subroutine RDF

use Numbers, only : N
use Configuration, only : R
use CommonDPD
use CellParam, only : CellL,InvCL

implicit none

integer :: i , j, ir
real(8) :: Cellh, R2, R1
real(8), dimension(3) :: Rij

   Cellh = 0.5 * CellL(1)

   do i = 1 , N-1

     do j = i+1 , N

       Rij = R(:,i) - R(:,j)
       Rij = Rij - dnint( Rij * InvCL ) * CellL
       R2  = dot_product( Rij , Rij )

       R1  = sqrt(R2)

       if(R1 < Cellh) then

         ir = nint( R1 * 20.d0 )
         irdf(ir) = irdf(ir) + 1

       end if

     end do

   end do

end subroutine RDF


!######################################################################
!######################################################################


! ************************************************
! **  write a configuration in the XMOL format  **
! ************************************************

subroutine Print_XYZ

use Numbers, only : N
use Configuration, only : R
use CommonDPD
use TimeParam, only : Timeps

implicit none

character(len=3), dimension(6), parameter :: &
&          AtomT = (/'C','O','N','P','S','F'/)
integer :: i, k
real, dimension(3,N) :: Rr
real, parameter :: scale = 3.

   open( 8 , file = 'final.xyz', status = 'unknown')

   Rr = sngl(R) * scale

   write(8,*) N
   write(8,'(a,f15.3)') 'time step =',timeps

   do i = 1 , N

     k = TypeNum(i)

     write(8,'(a3,x,3f10.5)') AtomT(k),Rr(:,i)

   end do

   close(8)

end subroutine Print_XYZ


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_HMC_NV(istep)

use CommonBlocks, only : QRigidBody, QCorrectCutoff, Qstdout
use CommonHMC
use EAM_param, only : Ene_EAM
use UnitExParam, only : rprs, cvol
use EwaldParam, only : Ene_Eksp, Ene_Eslf
use OptConstraintParam, only : Ene_OptC
use NonbondParam, only : Ene_LJ, Ene_Ersp
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use TailCorrect, only : Ene_LJ_co, Virial_co
use CellParam, only : Volume
use ThermoData, only : Temp, Temp_Translation, Temp_Rotation, &
&   Ene_kin, Pkinp, Virial

implicit none

integer :: i, istep
real(8), dimension(3,3) :: Pressure
real(8) :: Kinetic, Potential, EneSystem
real(8) :: E_BondMol, E_AngleMol, E_UBMol, E_DihedMol, E_ImproMol
real(8) :: E_LJMol, E_ErMol, E_EkMol, E_EsMol, E_OptC
real(8) :: E_ElecMol, E_Int_Mol, E_Ext_Mol, E_LJcoMol
real(8) :: Accept_ratio

   Ene_LJ = Ene_LJ + Ene_EAM

   if(QCorrectCutoff) then

     Virial = Virial - Virial_co

   end if

!-------------------------
   call CalcTemp
!-------------------------

   Accept_ratio = dble(NumAccept) / dble(istep) * 100.d0

! ---------------------------------------------

   Kinetic    = Ene_kin   * cvol * 0.5d0

   E_BondMol  = Ene_Bond  * cvol
   E_AngleMol = Ene_Angle * cvol
   E_UBMol    = Ene_UB    * cvol
   E_DihedMol = Ene_Dihed * cvol
   E_ImproMol = Ene_Impro * cvol
   E_OptC     = Ene_OptC  * cvol

   E_LJMol    = Ene_LJ    * cvol

   E_ErMol    = Ene_Ersp  * cvol
   E_EkMol    = Ene_Eksp  * cvol
   E_EsMol    = Ene_Eslf  * cvol
   E_ElecMol  = E_ErMol + E_EkMol + E_EsMol

   E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
   &          + E_DihedMol + E_ImproMol + E_OptC
   E_Ext_Mol  = E_LJMol + E_ElecMol

   if(QCorrectCutoff) then

     E_LJcoMol = Ene_LJ_co * cvol
     E_Ext_Mol = E_Ext_Mol + E_LJcoMol

   end if

   Potential  = E_Int_Mol + E_Ext_Mol

   EneSystem = Kinetic + Potential

   Pressure = ( Pkinp + Virial ) / Volume / rprs * 1.d-6 ! [MPa]

   if(QRigidBody) then

#ifdef BMONI
     write(11) TimeMC, Temp, Temp_Translation, &
     &         Temp_Rotation, Accept_Ratio
#else
     write(11,'(i12,3f10.2,f10.4)')    &
     &    TimeMC, Temp, Temp_Translation, &
     &    Temp_Rotation, Accept_Ratio
#endif
     if(Qstdout) then
     write( 6,'(i12,3f10.2,f10.4)')    &
     &    TimeMC, Temp, Temp_Translation, &
     &    Temp_Rotation, Accept_Ratio
     end if

   else

#ifdef BMONI
     write(11) TimeMC, Temp, Accept_Ratio
#else
     write(11,'(i12,f10.2,f10.6)')    &
     &    TimeMC, Temp, Accept_Ratio
#endif
     if(Qstdout) then
     write( 6,'(i12,f10.2,f10.6)')    &
     &    TimeMC, Temp, Accept_Ratio
     end if

   end if

#ifdef BMONI
   write(11) E_BondMol, E_AngleMol, E_UBMol,     &
   &         E_DihedMol,E_ImproMol, E_OptC,      &
   &         E_LJMol, E_ErMol,                   &
   &         E_EkMol, E_EsMol,   E_ElecMol,      &
   &         E_Int_Mol, E_Ext_Mol, Potential,    &
   &         Kinetic, EneSystem, Pressure
#else
   write(11,'(6d13.5/5d13.5/5d13.5/3(3f12.4/))') &
   &          E_BondMol, E_AngleMol, E_UBMol,     &
   &          E_DihedMol,E_ImproMol, E_OptC,      &
   &          E_LJMol, E_ErMol,                   &
   &          E_EkMol, E_EsMol,   E_ElecMol,      &
   &          E_Int_Mol, E_Ext_Mol, Potential,    &
   &          Kinetic, EneSystem,                 &
   &          ( Pressure(1,i) , i = 1 , 3 ),      &
   &          ( Pressure(2,i) , i = 1 , 3 ),      &
   &          ( Pressure(3,i) , i = 1 , 3 )
#endif
   if(Qstdout) then
   write( 6,'(6d13.5/5d13.5/5d13.5/3(3f12.4/))') &
   &          E_BondMol, E_AngleMol, E_UBMol,     &
   &          E_DihedMol,E_ImproMol, E_OptC,      &
   &          E_LJMol, E_ErMol,                   &
   &          E_EkMol, E_EsMol,   E_ElecMol,      &
   &          E_Int_Mol, E_Ext_Mol, Potential,    &
   &          Kinetic, EneSystem,                 &
   &          ( Pressure(1,i) , i = 1 , 3 ),      &
   &          ( Pressure(2,i) , i = 1 , 3 ),      &
   &          ( Pressure(3,i) , i = 1 , 3 )
   end if

end subroutine Print_Energy_HMC_NV


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_HMC_NP(istep)

use CommonBlocks, only : QCorrectCutoff, cBarostatMethod, ForceField, &
&   QRigidBody, Qstdout
use CommonHMC
use EAM_param, only : Ene_EAM
use UnitExParam, only : rprs, pi, cvol
use BathParam, only : Baro_kin, Pressure_o, SigmaS
use EwaldParam, only : Ene_Eksp, Ene_Eslf
use OptConstraintParam, only : Ene_OptC
use NonbondParam, only : Ene_LJ, Ene_Ersp
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use CellParam, only : H, Volume
use TailCorrect, only : Ene_LJ_co, Virial_co
use AtomParam, only : Mass
use ThermoData, only : Temp, Temp_Translation, Temp_Rotation, &
&   Ene_kin, Pkinp, Virial

implicit none

integer :: i, istep
real(8) :: Hamiltonian
real(8), dimension(3,3) :: Pressure
real(8) :: Kinetic, Potential, EneSystem
real(8) :: E_BondMol, E_AngleMol, E_UBMol, E_DihedMol, E_ImproMol
real(8) :: E_LJMol, E_ErMol, E_EkMol, E_EsMol, E_OptC
real(8) :: E_ElecMol, E_Int_Mol, E_Ext_Mol, E_LJcoMol
real(8) :: Pot_p, Kin_p
real(8) :: density, Area, TotalM
real(8), dimension(3) :: VecA, VecB, VecC
real(8) :: LenA, LenB, LenC, csAB, csBC, csCA
real(8) :: AngAB, AngBC, AngCA
real(8) :: Accept_ratio, det
! ## Stress >>
real(8), dimension(3,3) :: ElasM, Htrans, G
real(8) :: Pot_elastic
! ## << Stress
external det

   Ene_LJ = Ene_LJ + Ene_EAM

   Accept_ratio = dble(NumAccept) / dble(istep) * 100.d0

!-------------------------
   call CalcTemp
   call BathTemp
!-------------------------

! ---------------------------------------------

   Kinetic    = Ene_kin   * cvol * 0.5d0

   E_BondMol  = Ene_Bond  * cvol
   E_AngleMol = Ene_Angle * cvol
   E_UBMol    = Ene_UB    * cvol
   E_DihedMol = Ene_Dihed * cvol
   E_ImproMol = Ene_Impro * cvol
   E_OptC     = Ene_OptC  * cvol

   E_LJMol    = Ene_LJ    * cvol

   E_ErMol    = Ene_Ersp  * cvol
   E_EkMol    = Ene_Eksp  * cvol
   E_EsMol    = Ene_Eslf  * cvol

   E_ElecMol  = E_ErMol + E_EkMol + E_EsMol

   E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
   &          + E_DihedMol + E_ImproMol + E_OptC
   E_Ext_Mol  = E_LJMol    + E_ElecMol

   if(QCorrectCutoff) then

     E_LJcoMol = Ene_LJ_co * cvol
     E_Ext_Mol = E_Ext_Mol + E_LJcoMol

   end if

   Potential  = E_Int_Mol + E_Ext_Mol

   EneSystem = Kinetic + Potential

   Volume = det( H )

   Kin_p = 0.5d0 * Baro_kin    * cvol
   Pot_p = Pressure_o * Volume * cvol

! ## Stress >>
   if( cBarostatMethod == 'ST' ) then
     Htrans = transpose( H )
     G = matmul(Htrans,H)
     ElasM = matmul(SigmaS,G)
     Pot_elastic = 0.5 * ( ElasM(1,1) + ElasM(2,2) + ElasM(3,3) ) * cvol
     Pot_p = Pot_p + Pot_elastic
   end if
! ## << Stress

   Pressure = ( Pkinp + Virial ) / Volume / rprs * 1.d-6 ! [MPa]

   TotalM  = Sum(Mass)
   density = ( TotalM * 1.d3 ) / ( Volume * 1.d-24 ) ! [g/cm^3]

   VecA = H(:,1)
   VecB = H(:,2)
   VecC = H(:,3)

   LenA = sqrt( dot_product(VecA,VecA) )
   LenB = sqrt( dot_product(VecB,VecB) )
   LenC = sqrt( dot_product(VecC,VecC) )

   csAB = dot_product(VecA,VecB)/(LenA*LenB)
   csBC = dot_product(VecB,VecC)/(LenB*LenC)
   csCA = dot_product(VecC,VecA)/(LenC*LenA)
   AngAB = acos( csAB ) / pi * 180.d0
   AngBC = acos( csBC ) / pi * 180.d0
   AngCA = acos( csCA ) / pi * 180.d0

   Area = LenA * LenB * sin( acos( csAB ) )

   if(ForceField(1:3) == 'EAM') then
! ## cell matrix
     write(13) TimeMC,H
   end if


   if(QRigidBody) then

#ifdef BMONI
     write(11) TimeMC, Temp, Temp_Translation, &
     &         Temp_Rotation, Accept_Ratio
#else
     write(11,'(i12,3f10.2,f10.4)')     &
     &    TimeMC, Temp, Temp_Translation, &
     &    Temp_Rotation, Accept_Ratio
#endif
     if(Qstdout) then
     write( 6,'(i12,3f10.2,f10.4)')     &
     &    TimeMC, Temp, Temp_Translation, &
     &    Temp_Rotation, Accept_Ratio
     end if

   else

#ifdef BMONI
     write(11) TimeMC, Temp, Accept_Ratio
#else
     write(11,'(i12,f10.2,f10.4)')    &
     &    TimeMC, Temp, Accept_Ratio
#endif
     if(Qstdout) then
     write( 6,'(i12,f10.2,f10.4)')    &
     &    TimeMC, Temp, Accept_Ratio
     end if

   end if

   Hamiltonian = Kinetic + Potential + Kin_p + Pot_p

#ifdef BMONI
   write(11)  E_BondMol,  E_AngleMol, E_UBMol,         &
   &          E_DihedMol, E_ImproMol, E_OptC,          &
   &          E_LJMol, E_ErMol,                        &
   &          E_EkMol, E_EsMol,   E_ElecMol,           &
   &          E_Int_Mol, E_Ext_Mol, Potential,         &
   &          Kinetic, EneSystem
   if( cBarostatMethod == 'ST' ) then
     write(11) Pot_p
   end if
   write(11)  Hamiltonian, Area, Volume, density,      &
   &          LenA, LenB, LenC, AngBC, AngCA, AngAB
   write(11)  H, Pressure
#else
   write(11,'(6d13.5/5d13.5/5d13.5)')                  &
   &          E_BondMol,  E_AngleMol, E_UBMol,         &
   &          E_DihedMol, E_ImproMol, E_OptC,          &
   &          E_LJMol, E_ErMol,                        &
   &          E_EkMol, E_EsMol,   E_ElecMol,           &
   &          E_Int_Mol, E_Ext_Mol, Potential,         &
   &          Kinetic, EneSystem
   if( cBarostatMethod == 'ST' ) then
     write(11,'(d13.5)') Pot_p
   end if
   write(11,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')      &
   &          Hamiltonian, Area, Volume, density,      &
   &          LenA, LenB, LenC, AngBC, AngCA, AngAB
   write(11,'(3(3f8.3,5x,3f12.4/))')                   &
   &          (H(1,i),i=1,3) , (Pressure(1,i),i=1,3),  &
   &          (H(2,i),i=1,3) , (Pressure(2,i),i=1,3),  &
   &          (H(3,i),i=1,3) , (Pressure(3,i),i=1,3)
#endif
   if(Qstdout) then
   write( 6,'(6d13.5/5d13.5/5d13.5)')                  &
   &          E_BondMol,  E_AngleMol, E_UBMol,         &
   &          E_DihedMol, E_ImproMol, E_OptC,          &
   &          E_LJMol, E_ErMol,                        &
   &          E_EkMol, E_EsMol,   E_ElecMol,           &
   &          E_Int_Mol, E_Ext_Mol, Potential,         &
   &          Kinetic, EneSystem
   if( cBarostatMethod == 'ST' ) then
     write(6,'(d13.5)') Pot_p
   end if
   write( 6,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')      &
   &          Hamiltonian, Area, Volume, density,      &
   &          LenA, LenB, LenC, AngBC, AngCA, AngAB
   write( 6,'(3(3f8.3,5x,3f12.4/))')                   &
   &          (H(1,i),i=1,3) , (Pressure(1,i),i=1,3),  &
   &          (H(2,i),i=1,3) , (Pressure(2,i),i=1,3),  &
   &          (H(3,i),i=1,3) , (Pressure(3,i),i=1,3)
   end if

end subroutine Print_Energy_HMC_NP


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_NV_PI(istep)

use Numbers, only : N
use CommonBlocks, only : QMaster, QAveTh, QCorrectCutoff, cThermostatMethod, &
&   Qstdout
use CommonPI
use UnitExParam, only : rprs, cvol
use BathParam
use EwaldParam, only : Ene_Eslf
use TailCorrect, only : Ene_LJ_co, Virial_co
use CellParam, only : Volume
use TimeParam, only : itgn, lk, Timeps
use ThermoData, only : Temp, Ene_kin, Pkinp, Virial

implicit none

integer :: i, j, k, istep
real(8) :: Hamiltonian
real(8), dimension(3,3) :: Pressure
real(8) :: Kinetic, Potential, EneSystem
real(8), dimension(NHchain) :: Kin_ss
real(8), dimension(NHchain) :: Pot_ss
real(8) :: ThermoBath, Qkinetic, qdummy, ThermoBath2
! ## average
real(8), save :: Av_EnePI, Av_Qkinetic, Av_ThermoBath
real(8), save :: Av_Temp, Av_Ene_kin, fprint
real(8), dimension(3,3), save :: Av_Pkinp, Av_Vir

   if((.not.QAveTh).and.(mod(istep,itgn)/=0)) Return

!-------------------------
   call CalcTempPI
   call SumVir( Virial )
!-------------------------

!     /*  quantum kinetic energy (harmonic interaction):  *
!      *  primitive estimator                             */

   Qkinetic = 0.d0
   ThermoBath2 = 0.d0

   if(QMasterPI) then

!     /*  quantum kinetic energy (harmonic interaction):  *
!      *  primitive estimator                             */
     do j = IniBead, FinBead

       if(j==1) cycle
       do i = 1, N
         Qkinetic = Qkinetic + NmMass(i,j) * dot_product( Rnm(:,i,j),Rnm(:,i,j) )
       end do

     end do

     Qkinetic = 0.5d0 * OmegaP2 * Qkinetic * cvol

     do i = IniBead, FinBead

       if(i==1) cycle
       qdummy = Qmass(i) * 0.5d0

       do j = 1, NHchain

         do k = 1, N
           ThermoBath2 = ThermoBath2 &
           &           + qdummy * dot_product(Vbath(:,k,j,i),Vbath(:,k,j,i))     &
           &           + kT * ( Rbath(1,k,j,i) + Rbath(2,k,j,i) + Rbath(3,k,j,i) )
         end do

       end do

     end do

   end if

!--------------------------------------------
   call SumEnePI(EnePI,Qkinetic,ThermoBath2)
!--------------------------------------------

   if(QMaster) then

     EnePI = EnePI + Ene_Eslf * Pbead

     if(QCorrectCutoff) then
       Virial = Virial - Virial_co * Pbead
       EnePI  = EnePI  + Ene_LJ_co * Pbead
     end if

     if((cThermostatMethod == 'NH').or.(cThermostatMethod == 'NHC')) then

       Kin_ss = 0.5d0 * Mts * Vss * Vss * cvol
       Pot_ss = kT * Rss                * cvol
       Pot_ss(1) = gkT * Rss(1) * cvol

     else if( cThermostatMethod == 'MNHC' ) then

       Kin_ss = 0.d0
       Pot_ss = 0.d0
       do i = 1, NumMNHC
         Kin_ss(:) = Kin_ss(:) + MMNHC(:,i) * VMNHC(:,i) * VMNHC(:,i)
         Pot_ss(:) = Pot_ss(:) + RMNHC(:,i)
       end do
       Kin_ss = Kin_ss * 0.5d0 * cvol
       Pot_ss = Pot_ss * kT    * cvol

     end if

     ThermoBath = sum( Pot_ss ) + sum( Kin_ss ) + ThermoBath2 * cvol

     if(QAveTh) then

       if(mod(istep-lk,itgn)==0) then

         Av_Ene_kin  = 0.d0
         Av_EnePI    = 0.d0
         Av_Qkinetic = 0.d0

         Av_Temp  = 0.d0
         Av_Pkinp = 0.d0
         Av_Vir   = 0.d0

         Av_ThermoBath = 0.d0

       end if

       Av_Ene_kin  = Av_Ene_kin  + Ene_kin
       Av_EnePI    = Av_EnePI    + EnePI
       Av_Qkinetic = Av_Qkinetic + Qkinetic

       Av_Temp  = Av_Temp  + Temp
       Av_Pkinp = Av_Pkinp + Pkinp
       Av_Vir   = Av_Vir   + Virial

       Av_ThermoBath = Av_ThermoBath + ThermoBath

     end if

! ---------------------------------------------
     if(mod(istep,itgn)==0) then

       Kinetic    = Ene_kin * cvol * 0.5d0

       Potential  = EnePI * cvol * InvP

       EneSystem = Kinetic + Qkinetic + Potential

       Pressure = ( Pkinp + Virial ) * InvP / Volume / rprs * 1.d-6 ! [MPa]

       Hamiltonian = EneSystem + ThermoBath

#ifdef BMONI
       write(11) Timeps, Temp
#else
       write(11,'(f12.4,f10.2)') Timeps, Temp
#endif
       if(Qstdout) then
         write( 6,'(f12.4,f10.2)') Timeps, Temp
       end if

#ifdef BMONI
       write(11) Potential, Kinetic, Qkinetic, EneSystem, &
       &         ThermoBath, Hamiltonian, Pressure
#else
       write(11,'(5d13.5/d16.8/3(3f12.4/))') &
       &          Potential, Kinetic, Qkinetic, EneSystem, &
       &          ThermoBath, Hamiltonian,       &
       &          ( Pressure(1,i) , i = 1 , 3 ), &
       &          ( Pressure(2,i) , i = 1 , 3 ), &
       &          ( Pressure(3,i) , i = 1 , 3 )
#endif
       if(Qstdout) then
         write( 6,'(5d13.5/d16.8/3(3f12.4/))') &
         &          Potential, Kinetic, Qkinetic, EneSystem, &
         &          ThermoBath, Hamiltonian,       &
         &          ( Pressure(1,i) , i = 1 , 3 ), &
         &          ( Pressure(2,i) , i = 1 , 3 ), &
         &          ( Pressure(3,i) , i = 1 , 3 )
       end if

       if(QAveTh) then

         fprint = dble(lk) / dble(itgn)

         Av_Ene_kin  = Av_Ene_kin  * fprint
         Av_EnePI    = Av_EnePI    * fprint
         Av_Qkinetic = Av_Qkinetic * fprint

         Av_Temp  = Av_Temp  * fprint
         Av_Pkinp = Av_Pkinp * fprint
         Av_Vir   = Av_Vir   * fprint

         Av_ThermoBath = Av_ThermoBath * fprint

         Kinetic    = Av_Ene_kin * cvol * 0.5d0

         Potential  = Av_EnePI * cvol * InvP

         EneSystem = Kinetic + Av_Qkinetic + Potential

         Pressure = ( Av_Pkinp + Av_Vir ) * InvP / Volume / rprs * 1.d-6 ! [MPa]

         Hamiltonian = EneSystem + Av_ThermoBath

#ifdef BMONI
         write(12) Timeps, Av_Temp
         write(12) Potential, Kinetic, Av_Qkinetic, EneSystem, &
         &         Av_ThermoBath, Hamiltonian,Pressure
#else
         write(12,'(f12.4,f10.2)') Timeps, Av_Temp
         write(12,'(5d13.5/d16.8/3(3f12.4/))') &
         &          Potential, Kinetic, Av_Qkinetic, EneSystem, &
         &          Av_ThermoBath, Hamiltonian,   &
         &          ( Pressure(1,i) , i = 1, 3 ), &
         &          ( Pressure(2,i) , i = 1, 3 ), &
         &          ( Pressure(3,i) , i = 1, 3 )
#endif

       end if

     end if

   end if


end subroutine Print_Energy_NV_PI


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_NP_PI(istep)

use Numbers, only : N
use CommonBlocks, only : QMaster, Qstdout, cBarostatMethod, ForceField, &
&   QAveTh, QCorrectCutoff, cThermostatMethod
use CommonPI
use UnitExParam, only : rprs, pi, cvol
use BathParam
use EwaldParam, only : Ene_Eslf
use CellParam, only : H, Volume
use TailCorrect, only : Ene_LJ_co, Virial_co
use AtomParam, only : Mass
use TimeParam, only : itgn, lk, Timeps
use ThermoData, only : Temp, Ene_kin, Pkinp, Virial

implicit none

integer :: i, j, k, istep
real(8) :: Hamiltonian
real(8), dimension(3,3) :: Pressure
real(8) :: Kinetic, Potential, EneSystem
real(8), dimension(NHchain) :: Kin_ss
real(8), dimension(NHchain) :: Pot_ss
real(8) :: ThermoBath, ThermoBath2
real(8) :: Pot_p, Kin_p
real(8) :: density, Area, TotalM
real(8), dimension(3) :: VecA, VecB, VecC
real(8) :: LenA, LenB, LenC, csAB, csBC, csCA
real(8) :: AngAB, AngBC, AngCA
! ## Stress >>
real(8), dimension(3,3) :: ElasM, Htrans, G
real(8) :: Pot_elastic
! ## << Stress
real(8) :: Qkinetic, qdummy
real(8), dimension(3,3) :: Sckin
! ## average
real(8), save :: Av_EnePI, Av_Qkinetic, Av_ThermoBath
real(8), save :: Av_Temp, Av_Ene_kin, fprint
real(8), save :: Av_Volume, Av_Baro_kin, Av_Pot_elastic
real(8), dimension(3,3), save :: Av_Sckin, Av_Vir, Av_H

   if((.not.QAveTh).and.(mod(istep,itgn)/=0)) Return

!-------------------------
   call CalcTempPI
!-------------------------

!     /*  quantum kinetic energy (harmonic interaction):  *
!      *  primitive estimator                             */

   Qkinetic = 0.d0
   ThermoBath2 = 0.d0

   if(QMasterPI) then

     do j = IniBead, FinBead

       if(j==1) cycle
       do i = 1, N
         Qkinetic = Qkinetic + NmMass(i,j) * dot_product( Rnm(:,i,j),Rnm(:,i,j) )
       end do

     end do

     Qkinetic = 0.5d0 * OmegaP2 * Qkinetic * cvol

     do i = IniBead, FinBead

       if(i==1) cycle
       qdummy = 0.5d0 * Qmass(i)

       do j = 1, NHchain
         do k = 1, N
           ThermoBath2 = ThermoBath2 &
           &           + qdummy * dot_product(Vbath(:,k,j,i),Vbath(:,k,j,i)) &
           &           + kT * ( Rbath(1,k,j,i) + Rbath(2,k,j,i) + Rbath(3,k,j,i) )
         end do
       end do

     end do

   end if

!--------------------------------------------
   call SumEnePI(EnePI,Qkinetic,ThermoBath2)
!--------------------------------------------

   if(QMaster) then

!-------------------------
     call BathTemp
!-------------------------


     EnePI = EnePI + Ene_Eslf * Pbead

     if(QCorrectCutoff) then
       EnePI = EnePI + Ene_LJ_co * Pbead
     end if

! ## Stress >>
     if( cBarostatMethod == 'ST' ) then
       Htrans = transpose( H )
       G = matmul(Htrans,H)
       ElasM = matmul(SigmaS,G)
       Pot_elastic = 0.5 * ( ElasM(1,1) + ElasM(2,2) + ElasM(3,3) ) * cvol
     end if
! ## << Stress

     Sckin = 0.d0

     do i = 1, N
       do j = 1, 3
         do k = j, 3
           Sckin(j,k) = Sckin(j,k) + FictMass(i,1) * Vnm(j,i,1) * Vnm(k,i,1)
         end do
       end do
     end do

     Sckin(2,1) = Sckin(1,2)
     Sckin(3,1) = Sckin(1,3)
     Sckin(3,2) = Sckin(2,3)

     if((cThermostatMethod == 'NH').or.(cThermostatMethod == 'NHC')) then

       Kin_ss = 0.5d0 * Mts * Vss * Vss * cvol
       Pot_ss = kT * Rss                * cvol
       Pot_ss(1) = gkT * Rss(1) * cvol

     else if( cThermostatMethod == 'MNHC' ) then

       Kin_ss = 0.d0
       Pot_ss = 0.d0
       do i = 1, NumMNHC
         Kin_ss(:) = Kin_ss(:) + MMNHC(:,i) * VMNHC(:,i) * VMNHC(:,i)
         Pot_ss(:) = Pot_ss(:) + RMNHC(:,i)
       end do
       Kin_ss = Kin_ss * 0.5d0 * cvol
       Pot_ss = Pot_ss * kT    * cvol

     else if( cThermostatMethod == 'VSCALE' ) then

       Kin_ss = 0.d0
       Pot_ss = 0.d0

     end if

     ThermoBath = sum( Pot_ss ) + sum( Kin_ss ) + ThermoBath2 * cvol

     if(QAveTh) then

       if(mod(istep-lk,itgn)==0) then

         Av_Ene_kin  = 0.d0
         Av_EnePI    = 0.d0
         Av_Qkinetic = 0.d0

         Av_Temp  = 0.d0
         Av_Sckin = 0.d0
         Av_Vir   = 0.d0

         Av_ThermoBath = 0.d0

         Av_H        = 0.d0
         Av_Volume   = 0.d0
         Av_Baro_kin = 0.d0

         Av_Pot_elastic = 0.d0

       end if

       Av_Ene_kin  = Av_Ene_kin  + Ene_kin
       Av_EnePI    = Av_EnePI    + EnePI
       Av_Qkinetic = Av_Qkinetic + Qkinetic

       Av_Temp  = Av_Temp  + Temp
       Av_Sckin = Av_Sckin + Sckin
       Av_Vir   = Av_Vir   + Virial

       Av_ThermoBath = Av_ThermoBath + ThermoBath

       Av_H        = Av_H        + H       
       Av_Volume   = Av_Volume   + Volume  
       Av_Baro_kin = Av_Baro_kin + Baro_kin

       if( cBarostatMethod == 'ST' ) then
         Av_Pot_elastic = Av_Pot_elastic + Pot_elastic
       end if

     end if

! ---------------------------------------------
     if(mod(istep,itgn)==0) then

       Kinetic    = Ene_kin   * cvol * 0.5d0

       Potential  = EnePI * cvol * InvP

       EneSystem = Kinetic + Qkinetic + Potential

       Kin_p = 0.5d0 * Baro_kin    * cvol
       Pot_p = Pressure_o * Volume * cvol

! ## Stress >>
       if( cBarostatMethod == 'ST' ) then
         Pot_p = Pot_p + Pot_elastic
       end if
! ## << Stress

       Pressure = ( Sckin + Virial ) / Volume / rprs * 1.d-6 ! [MPa]

       TotalM  = Sum(Mass)
       density = ( TotalM * 1.d3 ) / ( Volume * 1.d-24 ) ! [g/cm^3]

       VecA = H(:,1)
       VecB = H(:,2)
       VecC = H(:,3)

       LenA = sqrt( dot_product(VecA,VecA) )
       LenB = sqrt( dot_product(VecB,VecB) )
       LenC = sqrt( dot_product(VecC,VecC) )

       csAB = dot_product(VecA,VecB)/(LenA*LenB)
       csBC = dot_product(VecB,VecC)/(LenB*LenC)
       csCA = dot_product(VecC,VecA)/(LenC*LenA)
       AngAB = acos( csAB ) / pi * 180.d0
       AngBC = acos( csBC ) / pi * 180.d0
       AngCA = acos( csCA ) / pi * 180.d0

       Area = LenA * LenB * sin( acos( csAB ) )

       if(ForceField(1:3) == 'EAM') then
         write(13) sngl(Timeps),H
       end if

#ifdef BMONI
       write(11) Timeps, Temp
#else
       write(11,'(f12.4,f10.2)') Timeps, Temp
#endif
       if(Qstdout) then
         write( 6,'(f12.4,f10.2)') Timeps, Temp
       end if

       Hamiltonian = EneSystem + ThermoBath + Kin_p + Pot_p

#ifdef BMONI
       write(11) Potential, Kinetic, Qkinetic, EneSystem, ThermoBath
       if( cBarostatMethod == 'ST' ) then
         write(11) Pot_p
       end if
       write(11)  Hamiltonian, Area, Volume, density,      &
       &          LenA, LenB, LenC, AngBC, AngCA, AngAB
       write(11)  H, Pressure
#else
       write(11,'(5d13.5)') Potential, Kinetic, Qkinetic, EneSystem, ThermoBath
       if( cBarostatMethod == 'ST' ) then
         write(11,'(d13.5)') Pot_p
       end if
       write(11,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')      &
       &          Hamiltonian, Area, Volume, density,      &
       &          LenA, LenB, LenC, AngBC, AngCA, AngAB
       write(11,'(3(3f8.3,5x,3f12.4/))')                   &
       &          (H(1,i),i=1,3) , (Pressure(1,i),i=1,3),  &
       &          (H(2,i),i=1,3) , (Pressure(2,i),i=1,3),  &
       &          (H(3,i),i=1,3) , (Pressure(3,i),i=1,3)
#endif
       if(Qstdout) then
         write( 6,'(5d13.5)') Potential, Kinetic, Qkinetic, EneSystem, ThermoBath
         if( cBarostatMethod == 'ST' ) then
           write( 6,'(d13.5)') Pot_p
         end if
         write( 6,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')      &
         &          Hamiltonian, Area, Volume, density,      &
         &          LenA, LenB, LenC, AngBC, AngCA, AngAB
         write( 6,'(3(3f8.3,5x,3f12.4/))')                   &
         &          (H(1,i),i=1,3) , (Pressure(1,i),i=1,3),  &
         &          (H(2,i),i=1,3) , (Pressure(2,i),i=1,3),  &
         &          (H(3,i),i=1,3) , (Pressure(3,i),i=1,3)
       end if

       if(QAveTh) then

         fprint = dble(lk) / dble(itgn)

         Av_Ene_kin  = Av_Ene_kin  * fprint
         Av_EnePI    = Av_EnePI    * fprint
         Av_Qkinetic = Av_Qkinetic * fprint

         Av_Temp  = Av_Temp  * fprint
         Av_Sckin = Av_Sckin * fprint
         Av_Vir   = Av_Vir   * fprint

         Av_ThermoBath = Av_ThermoBath * fprint

         Av_H         = Av_H         * fprint
         Av_Volume    = Av_Volume    * fprint
         Av_Baro_kin  = Av_Baro_kin  * fprint

         Av_Pot_elastic = Av_Pot_elastic * fprint


         Kinetic    = Av_Ene_kin   * cvol * 0.5d0

         Potential  = Av_EnePI * cvol * InvP

         EneSystem = Kinetic + Av_Qkinetic + Potential

         Kin_p = 0.5d0 * Av_Baro_kin    * cvol
         Pot_p = Pressure_o * Av_Volume * cvol

! ## Stress >>
         if( cBarostatMethod == 'ST' ) then
           Pot_p = Pot_p + Av_Pot_elastic
         end if
! ## << Stress

         Pressure = ( Av_Sckin + Av_Vir ) / Av_Volume / rprs * 1.d-6 ! [MPa]

         density = ( TotalM * 1.d3 ) / ( Av_Volume * 1.d-24 ) ! [g/cm^3]

         VecA = Av_H(:,1)
         VecB = Av_H(:,2)
         VecC = Av_H(:,3)

         LenA = sqrt( dot_product(VecA,VecA) )
         LenB = sqrt( dot_product(VecB,VecB) )
         LenC = sqrt( dot_product(VecC,VecC) )

         csAB = dot_product(VecA,VecB)/(LenA*LenB)
         csBC = dot_product(VecB,VecC)/(LenB*LenC)
         csCA = dot_product(VecC,VecA)/(LenC*LenA)
         AngAB = acos( csAB ) / pi * 180.d0
         AngBC = acos( csBC ) / pi * 180.d0
         AngCA = acos( csCA ) / pi * 180.d0

         Area = LenA * LenB * sin( acos( csAB ) )

#ifdef BMONI
         write(12) Timeps, Av_Temp

         Hamiltonian = EneSystem + Av_ThermoBath + Kin_p + Pot_p

         write(12) Potential, Kinetic, Av_Qkinetic, EneSystem, Av_ThermoBath
         if( cBarostatMethod == 'ST' ) then
           write(12) Pot_p
         end if
         write(12)  Hamiltonian, Area, Av_Volume, density,   &
         &          LenA, LenB, LenC, AngBC, AngCA, AngAB
         write(12) Av_H, Pressure
#else
         write(12,'(f12.4,f10.2)') Timeps, Av_Temp

         Hamiltonian = EneSystem + Av_ThermoBath + Kin_p + Pot_p

         write(12,'(5d13.5)') Potential, Kinetic, Av_Qkinetic, EneSystem, Av_ThermoBath
         if( cBarostatMethod == 'ST' ) then
           write(12,'(d13.5)') Pot_p
         end if
         write(12,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')      &
         &          Hamiltonian, Area, Av_Volume, density,   &
         &          LenA, LenB, LenC, AngBC, AngCA, AngAB
         write(12,'(3(3f8.3,5x,3f12.4/))')                   &
         &          (Av_H(1,i),i=1,3) , (Pressure(1,i),i=1,3),  &
         &          (Av_H(2,i),i=1,3) , (Pressure(2,i),i=1,3),  &
         &          (Av_H(3,i),i=1,3) , (Pressure(3,i),i=1,3)
#endif
       end if

     end if

   end if

end subroutine Print_Energy_NP_PI


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_iso_PI(istep)

use Numbers, only : N
use CommonBlocks, only : QMaster, Qstdout, QAveTh, cThermostatMethod
use CommonPI
use UnitExParam, only : cvol
use BathParam
use TimeParam, only : itgn, lk, Timeps
use ThermoData, only : Temp, Ene_kin

implicit none

integer :: i, j, k, istep
real(8) :: Hamiltonian
real(8) :: Kinetic, Potential, EneSystem
real(8), dimension(NHchain) :: Kin_ss
real(8), dimension(NHchain) :: Pot_ss
real(8) :: ThermoBath, Qkinetic, qdummy, ThermoBath2
! ## average
real(8), save :: Av_EnePI, Av_Qkinetic, Av_ThermoBath
real(8), save :: Av_Temp, Av_Ene_kin, fprint

!-------------------------
   call CalcTempPI
!-------------------------

   Qkinetic = 0.d0
   ThermoBath2 = 0.d0

   if(QMasterPI) then

!     /*  quantum kinetic energy (harmonic interaction):  *
!      *  primitive estimator                             */
     do j = IniBead, FinBead

       if(j==1) cycle
       do i = 1, N
         Qkinetic = Qkinetic + NmMass(i,j) * dot_product( Rnm(:,i,j),Rnm(:,i,j) )
       end do

     end do

     Qkinetic = 0.5d0 * OmegaP2 * Qkinetic * cvol

     do i = IniBead, FinBead

       if(i==1) cycle
       qdummy = Qmass(i) * 0.5d0

       do j = 1, NHchain

         do k = 1, N
           ThermoBath2 = ThermoBath2 &
           &           + qdummy * dot_product(Vbath(:,k,j,i),Vbath(:,k,j,i))     &
           &           + kT * ( Rbath(1,k,j,i) + Rbath(2,k,j,i) + Rbath(3,k,j,i) )
         end do

       end do

     end do

   end if

!--------------------------------------------
   call SumEnePI(EnePI,Qkinetic,ThermoBath2)
!--------------------------------------------

   if(QMaster) then

     if(cThermostatMethod == 'NHC') then

       Kin_ss = 0.5d0 * Mts * Vss * Vss * cvol
       Pot_ss = kT * Rss                * cvol
       Pot_ss(1) = gkT * Rss(1) * cvol

     else if( cThermostatMethod == 'MNHC' ) then

       Kin_ss = 0.d0
       Pot_ss = 0.d0
       do i = 1, NumMNHC
         Kin_ss(:) = Kin_ss(:) + MMNHC(:,i) * VMNHC(:,i) * VMNHC(:,i)
         Pot_ss(:) = Pot_ss(:) + RMNHC(:,i)
       end do
       Kin_ss = Kin_ss * 0.5d0 * cvol
       Pot_ss = Pot_ss * kT    * cvol

     end if

     ThermoBath = sum( Pot_ss ) + sum( Kin_ss )

     ThermoBath = ThermoBath + ThermoBath2 * cvol

     if(QAveTh) then

       if(mod(istep-lk,itgn)==0) then

         Av_Ene_kin  = 0.d0
         Av_EnePI    = 0.d0
         Av_Qkinetic = 0.d0

         Av_Temp  = 0.d0

         Av_ThermoBath = 0.d0

       end if

       Av_Ene_kin  = Av_Ene_kin  + Ene_kin
       Av_EnePI    = Av_EnePI    + EnePI
       Av_Qkinetic = Av_Qkinetic + Qkinetic

       Av_Temp  = Av_Temp  + Temp

       Av_ThermoBath = Av_ThermoBath + ThermoBath

     end if

! ---------------------------------------------
     if(mod(istep,itgn)==0) then

       Kinetic    = Ene_kin * cvol * 0.5d0

       Potential  = EnePI * cvol * InvP

       EneSystem = Kinetic + Qkinetic + Potential

#ifdef BMONI
       write(11) Timeps, Temp
#else
       write(11,'(f12.4,f10.2)') Timeps, Temp
#endif
       if(Qstdout) then
         write( 6,'(f12.4,f10.2)') Timeps, Temp
       end if

       Hamiltonian = EneSystem + ThermoBath

#ifdef BMONI
       write(11) Potential,Kinetic,Qkinetic, &
       &         ThermoBath,Hamiltonian
#else
       write(11,'(4d13.5,d16.8/)')            &
       &          Potential,Kinetic,Qkinetic, &
       &          ThermoBath,Hamiltonian
#endif
       if(Qstdout) then
         write( 6,'(4d13.5,d16.8/)')            &
         &          Potential,Kinetic,Qkinetic, &
         &          ThermoBath,Hamiltonian
       end if

       if(QAveTh) then

         fprint = dble(lk) / dble(itgn)

         Av_Ene_kin  = Av_Ene_kin  * fprint
         Av_EnePI    = Av_EnePI    * fprint
         Av_Qkinetic = Av_Qkinetic * fprint

         Av_Temp  = Av_Temp  * fprint

         Av_ThermoBath = Av_ThermoBath * fprint

         Kinetic    = Av_Ene_kin * cvol * 0.5d0

         Potential  = Av_EnePI * cvol * InvP

         EneSystem = Kinetic + Av_Qkinetic + Potential

#ifdef BMONI
         write(12) Timeps, Av_Temp

         Hamiltonian = EneSystem + Av_ThermoBath

         write(12)  Potential,Kinetic,Av_Qkinetic, &
         &          Av_ThermoBath,Hamiltonian
#else
         write(12,'(f12.4,f10.2)') Timeps, Av_Temp

         Hamiltonian = EneSystem + Av_ThermoBath

         write(12,'(4d13.5,d16.8/)')               &
         &          Potential,Kinetic,Av_Qkinetic, &
         &          Av_ThermoBath,Hamiltonian
#endif

       end if

     end if

   end if

end subroutine Print_Energy_iso_PI


! >> F monitor ##


!######################################################################
!######################################################################


! ***********************************************
! *   Monitoring Force exerted on Gramicidin A  *
! ***********************************************
!
subroutine Monitor_Force

use IOparam, only : DirectoryName, velocity_file, NtrjF
use Configuration, only : R
use CommonBlocks, only : Job_name
use F_monitor, only : NgA, NiniF, Fint, Fext, NfrcF, IfrcF, StoreFrc, &
&   force_file, confg_file

implicit none

integer :: i, j
real, dimension(3,NgA) :: RgA
real, dimension(3,NgA) :: Fexternal, Finternal

   IfrcF = IfrcF + 1

   if(IfrcF==1) then

     open(41,file=trim(DirectoryName)//trim(force_file),form='unformatted',status='new')
     open(51,file=trim(DirectoryName)//trim(confg_file),form='unformatted',status='new')

   end if

   do i = 1 , NgA
     j = NiniF + i
     Fexternal(:,i) = sngl(Fext(:,j))
     Finternal(:,i) = sngl(Fint(:,j))
   end do

   write(41) Fexternal , Finternal

   do i = 1 , NgA
     j = NiniF + i
     RgA(:,i) = sngl(R(:,j))
   end do

   write(51) RgA

   if(IfrcF==StoreFrc) then

     NfrcF = NfrcF + 1

     if(NfrcF<10) then

       write(force_file,'(a,a,i1)') trim(adjustl(Job_name)),'FRC000',NfrcF
       write(confg_file,'(a,a,i1)') trim(adjustl(Job_name)),'CRD000',NfrcF

     else if(NfrcF<100) then

       write(force_file,'(a,a,i2)') trim(adjustl(Job_name)),'FRC00',NfrcF
       write(confg_file,'(a,a,i2)') trim(adjustl(Job_name)),'CRD00',NfrcF

     else if(NfrcF<1000) then

       write(force_file,'(a,a,i3)') trim(adjustl(Job_name)),'FRC0',NfrcF
       write(confg_file,'(a,a,i3)') trim(adjustl(Job_name)),'CRD0',NfrcF

     else

       write(force_file,'(a,a,i4)') trim(adjustl(Job_name)),'FRC',NfrcF
       write(confg_file,'(a,a,i4)') trim(adjustl(Job_name)),'CRD',NfrcF

     end if

     close(41)
     close(51)

     IfrcF = 0

   end if

end subroutine Monitor_Force
! << F monitor ##


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_MM(FF)

use Numbers, only : N
use CommonBlocks, only : QMaster, ForceField, QCorrectCutoff
use EAM_param, only : Ene_EAM
use UnitExParam, only : rprs, cvol
use EwaldParam, only : Ene_Eksp, Ene_Eslf
use OptConstraintParam, only : Ene_OptC
use NonbondParam, only : Ene_Elec, Ene_LJ, Ene_Ersp, Ene_ELshrt, Ene_ELlong, &
&   Ene_NBshrt, Ene_NBlong
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use TailCorrect, only : Ene_LJ_co
use CellParam, only : Volume
use ThermoData, only : Virial

implicit none

integer :: i
real(8), dimension(3,3) :: Pressure
real(8) :: Potential
real(8) :: E_BondMol, E_AngleMol, E_UBMol, E_DihedMol, E_ImproMol
real(8) :: E_LJMol, E_ErMol, E_EkMol, E_EsMol, E_OptC
real(8) :: E_ElecMol, E_Int_Mol, E_Ext_Mol, E_LJcoMol
real(8), dimension(3,N) :: FF

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

! ---------------------------------------------
     E_BondMol  = Ene_Bond  * cvol
     E_AngleMol = Ene_Angle * cvol
     E_UBMol    = Ene_UB    * cvol
     E_DihedMol = Ene_Dihed * cvol
     E_ImproMol = Ene_Impro * cvol
     E_OptC     = Ene_OptC  * cvol

     E_LJMol    = Ene_LJ    * cvol

     E_ErMol    = Ene_Ersp  * cvol
     E_EkMol    = Ene_Eksp  * cvol
     E_EsMol    = Ene_Eslf  * cvol
     E_ElecMol  = E_ErMol + E_EkMol + E_EsMol

     E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
     &          + E_DihedMol + E_ImproMol + E_OptC
     E_Ext_Mol  = E_LJMol + E_ElecMol

     if(QCorrectCutoff) then

       E_LJcoMol = Ene_LJ_co * cvol
       E_Ext_Mol = E_Ext_Mol + E_LJcoMol

     end if

     Potential  = E_Int_Mol + E_Ext_Mol

     Pressure = Virial / Volume / rprs * 1.d-6 ! [MPa]

     write(11,'(3x,a/)') ' Energy [ kcal / mol ] '
     write(11,'(3x,a/,6(8x,a,d23.16/))') &
     &  '-----------------Intramolecular Energy------------------', &
     &  'U (Bond stretching)  = ',E_BondMol,   &
     &  'U (Angle bending)    = ',E_AngleMol,  &
     &  'U (UB type bending)  = ',E_UBMol,     &
     &  'U (Torsion)          = ',E_DihedMol,  &
     &  'U (Improper torsion) = ',E_ImproMol,  &
     &  'U (optional bond)    = ',E_OptC
     write(11,'(3x,a/,8x,a,d23.16/,8x,a/,4(8x,a,d23.16/))') &
     &  '-----------------Intermolecular Energy------------------', &
     &  'U (Lennard-Jones)    = ',E_LJMol,     &
     &  'U (Electrostatic)      ',             &
     &  '    (1.real space)   = ',E_ErMol,     &
     &  '    (2.k-space)      = ',E_EkMol,     &
     &  '    (3.self term)    = ',E_EsMol,     &
     &  '    (4.total)        = ',E_ElecMol
     write(11,'(3x,a/,3(8x,a,d23.16/)/8x,a,3d16.8/,2(24x,3d16.8/))') &
     &  '---------------------Total Energy-----------------------',  &
     &  'U (sum intramolecule)  = ',E_Int_Mol, &
     &  'U (sum nonbond )       = ',E_Ext_Mol, &
     &  'Total potential energy = ',Potential, &
     &  'Virial [MPa]  = ',  ( Pressure(1,i) , i = 1 , 3 ),      &
     &                       ( Pressure(2,i) , i = 1 , 3 ),      &
     &                       ( Pressure(3,i) , i = 1 , 3 )

     write( 6,'(3x,a/)') ' Energy [ kcal / mol ] '
     write( 6,'(3x,a/,6(8x,a,d23.16/))') &
     &  '-----------------Intramolecular Energy------------------', &
     &  'U (Bond stretching)  = ',E_BondMol,   &
     &  'U (Angle bending)    = ',E_AngleMol,  &
     &  'U (UB type bending)  = ',E_UBMol,     &
     &  'U (Torsion)          = ',E_DihedMol,  &
     &  'U (Improper torsion) = ',E_ImproMol,  &
     &  'U (optional bond)    = ',E_OptC
     write( 6,'(3x,a/,8x,a,d23.16/,8x,a/,4(8x,a,d23.16/))') &
     &  '-----------------Intermolecular Energy------------------', &
     &  'U (Lennard-Jones)    = ',E_LJMol,     &
     &  'U (Electrostatic)      ',             &
     &  '    (1.real space)   = ',E_ErMol,     &
     &  '    (2.k-space)      = ',E_EkMol,     &
     &  '    (3.self term)    = ',E_EsMol,     &
     &  '    (4.total)        = ',E_ElecMol
     write( 6,'(3x,a/,3(8x,a,d23.16/)/8x,a,3d16.8/,2(24x,3d16.8/))')&
     &  '---------------------Total Energy-----------------------', &
     &  'U (sum intramolecule)  = ',E_Int_Mol, &
     &  'U (sum nonbond )       = ',E_Ext_Mol, &
     &  'Total potential energy = ',Potential, &
     &  'Virial [MPa]  = ',  ( Pressure(1,i) , i = 1 , 3 ),      &
     &                       ( Pressure(2,i) , i = 1 , 3 ),      &
     &                       ( Pressure(3,i) , i = 1 , 3 )

     write(11,'(3x,a)') '---------------Force / kcal/mol/A--------------------'
     write(11,'(3x,a,3d20.12)') 'F_  1 = ',FF(:, 1)*cvol
     if(N>=100) write(11,'(3x,a,3d20.12)') 'F_100 = ',FF(:,100)*cvol
     if(N>=200) write(11,'(3x,a,3d20.12)') 'F_200 = ',FF(:,200)*cvol
     write(11,'(3x,a,3d20.12)') 'F_  N = ',FF(:,N)*cvol
     write(11,*)

     write( 6,'(3x,a)') '---------------Force / kcal/mol/A--------------------'
     write( 6,'(3x,a,3d20.12)') 'F_  1 = ',FF(:, 1)*cvol
     if(N>=100) write( 6,'(3x,a,3d20.12)') 'F_100 = ',FF(:,100)*cvol
     if(N>=200) write( 6,'(3x,a,3d20.12)') 'F_200 = ',FF(:,200)*cvol
     write( 6,'(3x,a,3d20.12)') 'F_  N = ',FF(:,N)*cvol
     write( 6,*)

   end if

end subroutine Print_Energy_MM


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_MM_iso(FF)

use Numbers, only : N
use CommonBlocks, only : QMaster, ForceField
use UnitExParam, only : cvol
use EwaldParam, only : Ene_Eksp
use OptConstraintParam, only : Ene_OptC
use NonbondParam, only : Ene_LJ, Ene_Elec, Ene_Ersp, Ene_ELshrt, Ene_ELlong, &
&   Ene_NBshrt, Ene_NBlong
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro

implicit none

real(8) :: Potential
real(8) :: E_BondMol, E_AngleMol, E_UBMol, E_DihedMol, E_ImproMol
real(8) :: E_LJMol, E_OptC
real(8) :: E_ElecMol, E_Int_Mol, E_Ext_Mol
real(8), dimension(3,N) :: FF

   if(ForceField(1:2) == 'CG') then
     Ene_Elec = Ene_ELshrt + Ene_ELlong
     Ene_LJ   = Ene_NBshrt + Ene_NBlong
   end if

!----------------------------------------------------------------------
   call SumEnergy(Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro, &
   &              Ene_LJ, Ene_Elec, Ene_Ersp, Ene_Eksp, Ene_OptC)
!----------------------------------------------------------------------

   if(QMaster) then

! ---------------------------------------------
     E_BondMol  = Ene_Bond  * cvol
     E_AngleMol = Ene_Angle * cvol
     E_UBMol    = Ene_UB    * cvol
     E_DihedMol = Ene_Dihed * cvol
     E_ImproMol = Ene_Impro * cvol
     E_OptC     = Ene_OptC  * cvol

     E_LJMol    = Ene_LJ    * cvol

     E_ElecMol  = Ene_Elec  * cvol

     E_Int_Mol  = E_BondMol  + E_AngleMol + E_UBMol &
     &          + E_DihedMol + E_ImproMol + E_OptC
     E_Ext_Mol  = E_LJMol    + E_ElecMol

     Potential  = E_Int_Mol + E_Ext_Mol

     write(11,'(3x,a/)') ' Energy [ kcal / mol ] '
     write(11,'(/3x,a/,6(8x,a,d23.16/)/3x,a/,2(8x,a,d23.16/)/3x,a/,3(8x,a,d23.16/))') &
     &  '-----------------Intramolecular Energy------------------', &
     &  'U (Bond stretching)  = ',E_BondMol,   &
     &  'U (Angle bending)    = ',E_AngleMol,  &
     &  'U (UB type bending)  = ',E_UBMol,     &
     &  'U (Torsion)          = ',E_DihedMol,  &
     &  'U (Improper torsion) = ',E_ImproMol,  &
     &  'U (optional bond)    = ',E_OptC,      &
     &  '-----------------Intermolecular Energy------------------', &
     &  'U (Lennard-Jones)    = ',E_LJMol,     &
     &  'U (Electrostatic)    = ',E_ElecMol,   &
     &  '---------------------Total Energy-----------------------', &
     &  'U (sum intramolecule)  = ',E_Int_Mol, &
     &  'U (sum nonbond )       = ',E_Ext_Mol, &
     &  'Total potential energy = ',Potential

     write( 6,'(3x,a/)') ' Energy [ kcal / mol ] '
     write( 6,'(/3x,a/,6(8x,a,d23.16/)/3x,a/,2(8x,a,d23.16/)/3x,a/,3(8x,a,d23.16/))') &
     &  '-----------------Intramolecular Energy------------------', &
     &  'U (Bond stretching)  = ',E_BondMol,   &
     &  'U (Angle bending)    = ',E_AngleMol,  &
     &  'U (UB type bending)  = ',E_UBMol,     &
     &  'U (Torsion)          = ',E_DihedMol,  &
     &  'U (Improper torsion) = ',E_ImproMol,  &
     &  'U (optional bond)    = ',E_OptC,      &
     &  '-----------------Intermolecular Energy------------------', &
     &  'U (Lennard-Jones)    = ',E_LJMol,     &
     &  'U (Electrostatic)    = ',E_ElecMol,   &
     &  '---------------------Total Energy-----------------------', &
     &  'U (sum intramolecule)  = ',E_Int_Mol, &
     &  'U (sum nonbond )       = ',E_Ext_Mol, &
     &  'Total potential energy = ',Potential

     write(11,'(3x,a)') '---------------Force / kcal/mol/A--------------------'
     write(11,'(3x,a,3d23.12)') 'F_  1 = ',FF(:, 1)*cvol
     if(N>=100) write(11,'(3x,a,3d20.12)') 'F_100 = ',FF(:,100)*cvol
     if(N>=200) write(11,'(3x,a,3d20.12)') 'F_200 = ',FF(:,200)*cvol
     write(11,'(3x,a,3d23.12)') 'F_  N = ',FF(:, N)*cvol
     write(11,*)

     write( 6,'(/3x,a)') 'Force / kcal/mol/A '
     write( 6,'(3x,a,3d20.12)') 'F_  1 = ',FF(:, 1)*cvol
     if(N>=100) write( 6,'(3x,a,3d20.12)') 'F_100 = ',FF(:,100)*cvol
     if(N>=200) write( 6,'(3x,a,3d20.12)') 'F_200 = ',FF(:,200)*cvol
     write( 6,'(3x,a,3d23.12)') 'F_  N = ',FF(:, N)*cvol
     write( 6,*)

   end if

end subroutine Print_Energy_MM_iso


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_DynaLib_NV

use CommonBlocks, only : QThermostat, cThermostatMethod
use UnitExParam, only : rprs, cvol
use BathParam
use CellParam, only : Volume
use TimeParam, only : Timeps
use ThermoData, only : Temp, Ene_kin, Pkinp

implicit none

integer :: i
real(8), dimension(3,3) :: Pressure
real(8) :: Kinetic
real(8), dimension(NHchain) :: Kin_ss
real(8), dimension(NHchain) :: Pot_ss
real(8) :: ThermoBath

!-------------------------
     call CalcTemp
!-------------------------

! ---------------------------------------------

     Kinetic    = Ene_kin   * cvol * 0.5d0

     Pressure = ( Pkinp ) / Volume / rprs * 1.d-6 ! [MPa]

     if(QThermostat) then

       if((cThermostatMethod == 'NH').or.(cThermostatMethod == 'NHC')) then

         Kin_ss = 0.5d0 * Mts * Vss * Vss * cvol
         Pot_ss = kT * Rss                * cvol
         Pot_ss(1) = gkT * Rss(1) * cvol

       else if( cThermostatMethod == 'MNHC' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0
         do i = 1, NumMNHC
           Kin_ss(:) = Kin_ss(:) + MMNHC(:,i) * VMNHC(:,i) * VMNHC(:,i)
           Pot_ss(:) = Pot_ss(:) + RMNHC(:,i)
         end do
         Kin_ss = Kin_ss * 0.5d0 * cvol
         Pot_ss = Pot_ss * kT    * cvol

       else if( cThermostatMethod == 'VSCALE' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0

       end if

       ThermoBath = sum( Pot_ss ) + sum( Kin_ss )

#ifdef BMONI
       write(11) Timeps, Temp
       write(11) Kinetic, ThermoBath, Pressure

     else

       write(11) Timeps, Temp
       write(11) Kinetic, Pressure
#else
       write(11,'(f12.4,f10.2)') Timeps, Temp
       write(11,'(2d16.8/3(3f12.4/))')                &
       &          Kinetic, ThermoBath,                &
       &          ( Pressure(1,i) , i = 1 , 3 ),      &
       &          ( Pressure(2,i) , i = 1 , 3 ),      &
       &          ( Pressure(3,i) , i = 1 , 3 )

     else

       write(11,'(f12.4,f10.2)') Timeps, Temp
       write(11,'(d16.8/3(3f12.4/))')                 &
       &          Kinetic,                            &
       &          ( Pressure(1,i) , i = 1 , 3 ),      &
       &          ( Pressure(2,i) , i = 1 , 3 ),      &
       &          ( Pressure(3,i) , i = 1 , 3 )
#endif

     end if


end subroutine Print_Energy_DynaLib_NV


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_DynaLib_NP

use CommonBlocks, only : QMaster, QThermostat, cBarostatMethod, cThermostatMethod
use UnitExParam, only : rprs, pi, cvol
use BathParam
use CellParam, only : H, Volume
use AtomParam, only : Mass
use TimeParam, only : Timeps
use ThermoData, only : Temp, Ene_kin, Pkinp

implicit none

integer :: i
real(8), dimension(3,3) :: Pressure
real(8) :: Kinetic
real(8), dimension(NHchain) :: Kin_ss
real(8), dimension(NHchain) :: Pot_ss
real(8) :: ThermoBath
real(8) :: Pot_p, Kin_p
real(8) :: density, Area, TotalM
real(8), dimension(3) :: VecA, VecB, VecC
real(8) :: LenA, LenB, LenC, csAB, csBC, csCA
real(8) :: AngAB, AngBC, AngCA
! ## Stress >>
real(8), dimension(3,3) :: ElasM, Htrans, G
real(8) :: Pot_elastic
! ## << Stress

   if(QMaster) then

!-------------------------
     call CalcTemp
     call BathTemp
!-------------------------

! ---------------------------------------------

     Kinetic    = Ene_kin   * cvol * 0.5d0

     Kin_p = 0.5d0 * Baro_kin    * cvol
     Pot_p = Pressure_o * Volume * cvol

! ## Stress >>
     if( cBarostatMethod == 'ST' ) then
       Htrans = transpose( H )
       G = matmul(Htrans,H)
       ElasM = matmul(SigmaS,G)
       Pot_elastic = 0.5 * ( ElasM(1,1) + ElasM(2,2) + ElasM(3,3) ) * cvol
       Pot_p = Pot_p + Pot_elastic
     end if
! ## << Stress

     Pressure = ( Pkinp ) / Volume / rprs * 1.d-6 ! [MPa]

     TotalM  = Sum(Mass)
     density = ( TotalM * 1.d3 ) / ( Volume * 1.d-24 ) ! [g/cm^3]

     VecA = H(:,1)
     VecB = H(:,2)
     VecC = H(:,3)

     LenA = sqrt( dot_product(VecA,VecA) )
     LenB = sqrt( dot_product(VecB,VecB) )
     LenC = sqrt( dot_product(VecC,VecC) )

     csAB = dot_product(VecA,VecB)/(LenA*LenB)
     csBC = dot_product(VecB,VecC)/(LenB*LenC)
     csCA = dot_product(VecC,VecA)/(LenC*LenA)
     AngAB = acos( csAB ) / pi * 180.d0
     AngBC = acos( csBC ) / pi * 180.d0
     AngCA = acos( csCA ) / pi * 180.d0

     Area = LenA * LenB * sin( acos( csAB ) )

#ifdef BMONI
     write(11) Timeps, Temp
#else
     write(11,'(f12.4,f10.2)') Timeps, Temp
#endif
     if(QThermostat) then

       if((cThermostatMethod == 'NH').or.(cThermostatMethod == 'NHC')) then

         Kin_ss = 0.5d0 * Mts * Vss * Vss * cvol
         Pot_ss = kT * Rss                * cvol
         Pot_ss(1) = gkT * Rss(1) * cvol

       else if( cThermostatMethod == 'MNHC' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0
         do i = 1, NumMNHC
           Kin_ss(:) = Kin_ss(:) + MMNHC(:,i) * VMNHC(:,i) * VMNHC(:,i)
           Pot_ss(:) = Pot_ss(:) + RMNHC(:,i)
         end do
         Kin_ss = Kin_ss * 0.5d0 * cvol
         Pot_ss = Pot_ss * kT    * cvol

       else if( cThermostatMethod == 'VSCALE' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0

       end if

       ThermoBath = sum( Pot_ss ) + sum( Kin_ss )

#ifdef BMONI
       write(11) Kinetic, ThermoBath, Kin_p, Pot_p
       write(11) Pressure

     else

       write(11) Kinetic, Kin_p, Pot_p
       write(11) Pressure
#else
       write(11,'(4d16.8)') Kinetic, ThermoBath, Kin_p, Pot_p

       write(11,'(3(3f12.4/))')            &
       &           (Pressure(1,i),i=1,3),  &
       &           (Pressure(2,i),i=1,3),  &
       &           (Pressure(3,i),i=1,3)

     else

       write(11,'(3d16.8)') Kinetic, Kin_p, Pot_p

       write(11,'(3(3f12.4/))')            &
       &           (Pressure(1,i),i=1,3),  &
       &           (Pressure(2,i),i=1,3),  &
       &           (Pressure(3,i),i=1,3)
#endif

     end if

   end if

end subroutine Print_Energy_DynaLib_NP


!######################################################################
!######################################################################


! **************************************
! ** output thermodynamic quantities  **
! **************************************

subroutine Print_Energy_DynaLib_iso

use CommonBlocks, only : QMaster, QThermostat, cThermostatMethod
use UnitExParam, only : cvol
use BathParam, only : NHchain, Mts, Rss, Vss, &
&   NumMNHC, MMNHC, RMNHC, VMNHC, kT, gkT
use TimeParam, only : Timeps
use ThermoData, only : Temp, Ene_kin

implicit none

integer :: i
real(8) :: Kinetic
real(8), dimension(NHchain) :: Kin_ss
real(8), dimension(NHchain) :: Pot_ss
real(8) :: ThermoBath

   if(QMaster) then

!-------------------------
     call CalcTemp
!-------------------------

! ---------------------------------------------

     Kinetic    = Ene_kin   * cvol * 0.5d0

#ifdef BMONI
     write(11) Timeps, Temp
#else
     write(11,'(f12.4,f10.2)') Timeps, Temp
#endif
     if(QThermostat) then

       if((cThermostatMethod == 'NH').or.(cThermostatMethod == 'NHC')) then

         Kin_ss = 0.5d0 * Mts * Vss * Vss * cvol
         Pot_ss = kT * Rss                * cvol
         Pot_ss(1) = gkT * Rss(1) * cvol

       else if( cThermostatMethod == 'MNHC' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0
         do i = 1, NumMNHC
           Kin_ss(:) = Kin_ss(:) + MMNHC(:,i) * VMNHC(:,i) * VMNHC(:,i)
           Pot_ss(:) = Pot_ss(:) + RMNHC(:,i)
         end do
         Kin_ss = Kin_ss * 0.5d0 * cvol
         Pot_ss = Pot_ss * kT    * cvol

       else if( cThermostatMethod == 'VSCALE' ) then

         Kin_ss = 0.d0
         Pot_ss = 0.d0

       end if

       ThermoBath = sum( Pot_ss ) + sum( Kin_ss )

#ifdef BMONI
       write(11) Kinetic,ThermoBath
#else
       write(11,'(2d16.8/)') Kinetic,ThermoBath
#endif
     else

#ifdef BMONI
       write(11) Kinetic
#else
       write(11,'(d16.8/)') Kinetic
#endif
     end if

   end if

end subroutine Print_Energy_DynaLib_iso


!######################################################################
!######################################################################


subroutine Print_Eflux(istep)

use TimeParam, only : lk, Timeps
use AtomParam, only : Mass
use Configuration, only : Vel
use Numbers, only : NumSpec, NumMol, NumAtm
use CommonBlocks, only : Job_name
use IOparam, only : DirectoryName
use Conduct

implicit none

integer :: istep, i, j, k, l, ii
integer :: Ncat, Nani
real(8), dimension(:,:), allocatable :: Vcat, Vani
real(8), dimension(3) :: Eflx
character(len=80) :: efile

   if(istep==lk) then
     write(efile,'(a,a)') trim(adjustl(Job_name)),'_Eflux.dat'
     open(59,file=trim(DirectoryName)//trim(efile),status='unknown')

     ii = 0
     do i = 1, NumSpec
       if(i==id_cation) then
         invMasscat = 0.d0
         k = ii
         do j = 1, NumAtm(i)
           k = k + 1
           invMasscat = invMasscat + Mass(k)
         end do
         invMasscat = 1.d0 / invMasscat
       else if(i==id_anion) then
         invMassani = 0.d0
         k = ii
         do j = 1, NumAtm(i)
           k = k + 1
           invMassani = invMassani + Mass(k)
         end do
         invMassani = 1.d0 / invMassani
       end if
       ii = ii + NumMol(i)*NumAtm(i)
     end do
   end if

   ii = 0
   do i = 1, NumSpec
     if(i==id_cation) then
       Ncat = NumMol(i)
       allocate(Vcat(3,Ncat))
       l = ii
       do j = 1, Ncat
         Vcat(:,j) = 0.d0
         do k = 1, NumAtm(i)
           l = l + 1
           Vcat(:,j) = Vcat(:,j) + Mass(l)*Vel(:,l)
         end do
         Vcat(:,j) = Vcat(:,j) * invMasscat
       end do
     else if(i==id_anion) then
       Nani = NumMol(i)
       allocate(Vani(3,Nani))
       l = ii
       do j = 1, Nani
         Vani(:,j) = 0.d0
         do k = 1, NumAtm(i)
           l = l + 1
           Vani(:,j) = Vani(:,j) + Mass(l)*Vel(:,l)
         end do
         Vani(:,j) = Vani(:,j) * invMassani
       end do
     end if
     ii = ii + NumMol(i)*NumAtm(i)
   end do

   Eflx = 0.d0
   do i = 1, Ncat
     Eflx = Eflx + Vcat(:,i)
   end do
   do i = 1, Nani
     Eflx = Eflx - Vani(:,i)
   end do

   write(59,'(d16.8,3d24.16)') Timeps,Eflx(:)

   deallocate(Vcat,Vani)

end subroutine Print_Eflux
