module RTR
  logical :: QThermostatNEW, QBarostatNEW
  logical :: QPathIntNEW, QPINPTNEW, QRigidBodyNEW
  character(len=6) :: cThermostatMethodNEW
  character(len=3) :: cBarostatMethodNEW
  character(len=10) :: SimMethodNEW
  real(8) :: Tau_NEW, OmegaP2NEW, GammaC2NEW
  real(8) :: Temp_NEW, Time_NEW
  integer :: NbeadNEW, NHchainNEW

  logical :: QChange_Method, QBathMomentRefresh, QChange_ThermostatMethod
  logical :: QChange_Constraint, QChange_Thermostat, QChange_Barostat
  logical :: QChange_Temperature, QChange_BarostatMethod

  real(8), dimension(:), allocatable :: Rnh,Vnh,Mnh
  real(8), dimension(:,:), allocatable :: Rmnh,Vmnh,Mmnh
  real(8), dimension(:,:,:,:), allocatable :: Rbb,Vbb
  real(8), dimension(:,:,:), allocatable :: RnmNEW, VnmNEW
  real(8), dimension(:,:), allocatable :: MnmNEW, FMnmNEW

  logical :: QExtension
  integer :: ix_ext, iy_ext, iz_ext

  logical :: QMove
  real(8) :: deltaX, deltaY, deltaZ

end module RTR


! ####################################################################
! ####################################################################


program ResTORes

use CommonBlocks, only : QMaster
use RTR
use CommonPI

implicit none

real(8) :: Omeg, sv, QmassNEW, vv, Gauss
integer :: i, j, k, l, nMNHC
character(len=1) :: Ch
external Gauss

! ## Serial job is assumed 
   QMaster = .True.

! ## reading setup routine : input.data 
   call Setup

   call Checkup

   call ReadNewCondition

   if(QExtension) call extend_sys
   if(QMove) call move_sys

! ## MD  
   if(SimMethodNEW=='MD') then
     call MDrestart
! ## PIMD 
   else if(SimMethodNEW=='PIMD') then
     call PIMDrestart
   end if

end program ResTORes


! ####################################################################
! ####################################################################


subroutine ReadNewCondition

use CommonBlocks, only : QMaster, QRigidBody, QBarostat, QThermostat, &
&   cThermostatMethod, cBarostatMethod, SimMethod, QResForm
use CommonPI
use BathParam
use RTR
use UnitExParam, only : Plank2pi
use TimeParam, only : Timeps

implicit none

integer, parameter :: MaxOption = 200
integer, parameter :: MaxScript = 2000
character(len=20), dimension(MaxOption), save :: OPTIONS
character(len=80), dimension(MaxScript), save :: Script
integer, dimension(MaxOption), save :: Ist
logical, dimension(MaxOption), save :: ReadFlag
integer , save :: NumOption, NumScript

   call Read_Script

   call Read_Simulation_Method

   call Read_Bath_Refresh

   call Read_Time_Refresh

   call Read_Data_Form

   call Read_Extension

   call Read_Move

   call Read_Constraint_Cond

   call Read_Thermostat_Cond

   call Read_Barostat_Cond

! #####################################################################

Contains

   subroutine Read_Script

   implicit none

   character(len=80) :: String, String1
   character(len=1) :: Ch

   integer :: NumCh, i
   integer :: Icount
   integer :: eofile

     open(2,file='input_convert.data',status='old')

     Icount = 0
     NumOption = 0

     do

       read(2,'(a)',iostat=eofile) String1

       String = adjustl(String1)

       if(eofile == -1) exit
       if(String(1:5)=='<end>') exit
       if(String(1:5)=='<END>') exit

       if((String(1:1) == '!').or. &
       &  (String(1:1) == '#').or. &
       &  (String(1:1) == ' ')) cycle

       Icount = Icount + 1

       if(Icount > MaxScript) then
         if(QMaster) write(*,*) 'ERROR : too many lines in "input.data"'
         stop
       end if

       NumCh = 0

id0:   do i = 1, 80

         Ch = String(i:i)

         if(Ch == '!') exit id0

         NumCh = NumCh + 1

       end do id0

       write(Script(Icount),'(a)') adjustl(String(1:NumCh))

       if(String(1:2)=='>>') then

         NumOption = NumOption + 1

         if(NumOption > MaxOption) then
           if(QMaster) write(*,*) 'ERROR : too many options in "input.data"'
           stop
         end if

         Ist(NumOption) = Icount

         String1 = adjustl(String(3:NumCh))

         write(OPTIONS(NumOption),'(a)') trim(String1)

       end if

     end do

     close(2)

     NumScript = Icount

     ReadFlag = .False.

   end subroutine Read_Script


! #####################################################################


   subroutine Read_Simulation_Method

   implicit none

   integer :: IJOB, i

     SimMethodNEW = SimMethod
     QChange_Method = .False.


     do i = 1, NumOption

       if(OPTIONS(i)=='METHOD') then

         IJOB = Ist(i) + 1
         read(Script(IJOB),*) SimMethodNEW

         ReadFlag(i) = .True.

         exit

       end if

     end do

     if((SimMethodNEW == 'PIMD').or.(SimMethodNEW == 'CMD')) then
       QPathIntNEW = .True.
     else
       QPathIntNEW = .False.
     end if

     if(SimMethodNEW /= SimMethod) QChange_Method = .True.

     if((SimMethodNEW == 'MD').and.((SimMethod == 'PIMD').or.(SimMethod == 'CMD')&
     &  .or.(SimMethod == 'HMC'))) then
       write(*,'(/a)') 'ERROR : Xconvert cannot generate the initial configuration'
       write(*,'( a)') '        in case of PIMD --> MD or CMD --> MD or HMC --> MD'
       write(*,'(/a)') '        You may use final.pdb from PIMD simulation to make'
       write(*,'( a)') '        initial.crd for a new MD simulation using         '
       write(*,'(a/)') '        MPDyn_Prep instead!                               '
       stop
     end if

   end subroutine Read_Simulation_Method


! #####################################################################


   subroutine Read_Bath_Refresh

   implicit none

   integer :: i

     QBathMomentRefresh = .False.

     do i = 1, NumOption

       if(OPTIONS(i)=='BATHREFRESH') then
         QBathMomentRefresh = .True.
       end if

     end do

   end subroutine Read_Bath_Refresh


! #####################################################################


   subroutine Read_Time_Refresh

   implicit none

   integer :: i

     Time_NEW = Timeps

     do i = 1, NumOption

       if(OPTIONS(i)=='TIMEREFRESH') then
         Time_NEW = 0.d0
       end if

     end do

   end subroutine Read_Time_Refresh


! #####################################################################


   subroutine Read_Data_Form

   implicit none

   integer :: i
   logical :: QResFormNEW

     QResFormNEW = QResForm

     do i = 1, NumOption

       if(OPTIONS(i)=='CHANGE_FORM') then
         if(QResForm) then
           QResFormNEW = .False.
         else
           QResFormNEW = .True.
         end if
       end if

     end do

     QResForm = QResFormNEW

   end subroutine Read_Data_Form


! #####################################################################


   subroutine Read_Extension

   implicit none

   integer :: i, ii
   logical :: QResFormNEW

     QExtension = .False.

     do i = 1, NumOption

       if(OPTIONS(i)=='EXTENSION') then

         QExtension = .True.

         ix_ext = 1
         iy_ext = 1
         iz_ext = 1

         ii = Ist(i)

lia10:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia10

           if(Script(ii)(1:2) == 'X=') then
             read(Script(ii)(3:),*) ix_ext
           else if(Script(ii)(1:2) == 'Y=') then
             read(Script(ii)(3:),*) iy_ext
           else if(Script(ii)(1:2) == 'Z=') then
             read(Script(ii)(3:),*) iz_ext
           end if

         end do lia10

         exit

       end if

     end do

   end subroutine Read_Extension


! #####################################################################


   subroutine Read_Move

   implicit none

   integer :: i, ii

     QMove = .False.

     do i = 1, NumOption

       if(OPTIONS(i)=='MOVE') then

         QMove = .True.

         deltaX = 0.d0
         deltaY = 0.d0
         deltaZ = 0.d0
         ii = Ist(i)

lia20:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia20

           if(Script(ii)(1:2) == 'X=') then
             read(Script(ii)(3:),*) deltaX
           else if(Script(ii)(1:2) == 'Y=') then
             read(Script(ii)(3:),*) deltaY
           else if(Script(ii)(1:2) == 'Z=') then
             read(Script(ii)(3:),*) deltaZ
           end if

         end do lia20

         exit

       end if

     end do

   end subroutine Read_Move


! #####################################################################


   subroutine Read_Constraint_Cond

   implicit none

   integer :: i, ii
   character(len=10) :: cConstraint

       QRigidBodyNEW = QRigidBody
       QChange_Constraint = .False.

       do i = 1, NumOption

         if(OPTIONS(i)=='CONSTRAINT') then

           ii = Ist(i)

lia9:      do

             ii = ii + 1

             if(Script(ii)(1:2) == '<<') exit lia9

             if(Script(ii)(1:7) == 'METHOD=') then

               read(Script(ii)(8:),*) cConstraint

               if(cConstraint=='RigidBody') then  !<<new>>

                 QRigidBodyNEW = .True.

               else

                 if(QMaster) write(*,*) 'ERROR : Constraint condition'

               end if

             end if

           end do lia9

           exit

         end if

       end do

       if((QRigidBodyNEW.and.(.not.QRigidBody)).or. &
       &  ((.not.QRigidBodyNEW).and.QRigidBody)) then
         QChange_Constraint = .True.
       end if

   end subroutine Read_Constraint_Cond


! #####################################################################


   subroutine Read_Thermostat_Cond

   implicit none

   integer :: i, ii

       QChange_Thermostat = .False.
       QChange_ThermostatMethod = .False.
       QThermostatNEW=QThermostat
       if(QThermostat) then
         cThermostatMethodNEW = cThermostatMethod
         Temp_NEW = Temp_o
         NHchainNEW= NHchain
       end if

       do i = 1, NumOption

         if(OPTIONS(i)=='TEMPERATURE') then

           Temp_NEW = 298.0
           NHchainNEW = 1
           Tau_NEW = 0.5

           ii = Ist(i)

lia8:      do

             ii = ii + 1

             if(Script(ii)(1:2) == '<<') exit lia8

             if(Script(ii)(1:7) == 'METHOD=') then

               read(Script(ii)(8:),*) cThermostatMethodNEW

               if(cThermostatMethodNEW /= cThermostatMethod) then
                 QThermostatNEW=.True.
                 QChange_ThermostatMethod = .True.
               end if
               if(cThermostatMethodNEW == 'VSCALE') then
                 QThermostatNEW=.False.
               end if

             else if(Script(ii)(1:5) == 'TEXT=') then

               read(Script(ii)(6:),*) Temp_NEW

             else if(Script(ii)(1:7) == 'LENGTH=') then

               read(Script(ii)(8:),*) NHchainNEW

             else if(Script(ii)(1:4) == 'TAU=') then

               read(Script(ii)(5:),*) Tau_NEW

             end if

           end do lia8

           ReadFlag(i) = .True.

           exit

         end if

       end do

       if(QThermostatNEW) then

         if(QPathIntNEW) then
           if((cThermostatMethodNEW /= 'MNHC').and.(cThermostatMethodNEW /= 'NHC')) then
             if(QMaster) write(*,*) &
             & 'WARNING : (massive) Nose-Hoover chain method should be', &
             & 'chosen when PIMD or CMD simulation is undertaken.'
             stop
           end if
         end if

       end if

       if((QThermostatNEW.and.(.not.QThermostat)).or. &
       &  ((.not.QThermostatNEW).and.QThermostat)) then
         QChange_Thermostat = .True.
       end if

       if(cThermostatMethodNEW/=cThermostatMethod) then
         QChange_ThermostatMethod = .True.
       end if

       if(Temp_NEW/=Temp_o) then
         QChange_Temperature = .True.
       end if

   end subroutine Read_Thermostat_Cond


! #####################################################################

   subroutine Read_Barostat_Cond

   implicit none

   integer :: i, ii

       QChange_Barostat = .False.
       QChange_BarostatMethod = .False.
       QBarostatNEW=QBarostat
       if(QBarostat) then
         cBarostatMethodNEW = cBarostatMethod
       end if

       do i = 1, NumOption

         if(OPTIONS(i)=='PRESSURE') then

           QBarostatNEW=.True.

           ii = Ist(i)

lia5:      do

             ii = ii + 1

             if(Script(ii)(1:2) == '<<') exit lia5

             if(Script(ii)(1:7) == 'METHOD=') then

               read(Script(ii)(8:),*) cBarostatMethodNEW

             end if

           end do lia5

           exit

         end if

       end do

       if((QBarostatNEW.and.(.not.QBarostat)).or. &
       &  ((.not.QBarostatNEW).and.QBarostat)) then
         QChange_Barostat = .True.
       end if

       if(cBarostatMethodNEW/=cBarostatMethod) then
         QChange_BarostatMethod = .True.
       end if

       if(cBarostatMethodNEW(1:3)=='OFF') then
         QBarostatNEW = .False.
       end if

   end subroutine Read_Barostat_Cond


! #####################################################################

   subroutine Read_Bead_Length

   implicit none

   integer :: i, ii
   real(8) :: xxx

     if(QPathIntNEW) then

       NbeadNEW = Nbead

       do i = 1, NumOption

         if(OPTIONS(i)=='PI_PARAM') then

           ii = Ist(i) + 1

           if(Script(ii)(1:12) == 'BEAD_LENGTH=') then

               read(Script(ii)(13:),*) NbeadNEW

           else

             if(QMaster) write(*,*) 'ERROR : PI_PARAM option!'
             stop

           end if

           ReadFlag(i) = .True.

           exit

         end if

       end do

       xxx  = sqrt(dble(NbeadNEW)) * kT / Plank2pi
       OmegaP2new = xxx * xxx

     end if

   end subroutine Read_Bead_Length

! #####################################################################

   subroutine Read_CMD_param

   implicit none

   integer :: i, ii
   real(8) :: TimeScalingFactor, xx

     if(QPathIntNEW) then

       TimeScalingFactor = 1.d0

       do i = 1, NumOption

         if(OPTIONS(i)=='CMD_PARAM') then

           ii = Ist(i) + 1

           if(Script(ii)(1:11) == 'TIME_SCALE=') then

               read(Script(ii)(12:),*) TimeScalingFactor

           else

             if(QMaster) write(*,*) 'ERROR : CMD_PARAM option!'
             stop

           end if

           exit

         end if

       end do

       xx = 1.d0 / TimeScalingFactor
       GammaC2NEW = xx * xx

     end if

   end subroutine Read_CMD_param

end subroutine ReadNewCondition


! ####################################################################
! ####################################################################


subroutine Checkup

use CommonBlocks, only : QInitial
use TimeParam, only : Timeps

implicit none

character(len=1) :: Ch

   write(*,'(/a)')        '--------------------------'
   write(*,'(a,f10.3,a)') '##  TIME = ',Timeps,' ##'
   write(*,'(a/)')        '--------------------------'

   if(QInitial) then

     write(*,'(/a)') '-------------------------------------------------------------'
     write(*,'(a)')  '%%%%%%%%%%%%%%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%'
     write(*,'(a)')  '  You may have to choose the flag "STATUS" in "input.data"   '
     write(*,'(a)')  '  to "Restart", though your choice was "Initial".            '
     write(*,'(a)')  '  NOW the program read "initial.crd" instead of "restart.dat"'
     write(*,'(a/)') '-------------------------------------------------------------'
     write(*,'(a/)') '-->  Will you continue this job ? [ Y / N ] '

     do
       read(*,*) Ch
       if((Ch=='Y').or.(Ch=='y')) then
         exit
       else if((Ch=='N').or.(Ch=='n')) then
         stop
       else
         write(*,*) 'try again; answer with "Y" or "N"'
         cycle
       end if
     end do

   end if

   write(*,'(/a/)')  ' ... Reading a NEW condition (input_convert.data) ... '

end subroutine Checkup


! ####################################################################
! ####################################################################


subroutine extend_sys

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration
use AtomParam, only : DefaultAtomName, DefaultResidName, ResidNum
use CellParam, only : H
use BathParam, only : InvMts
use RTR

implicit none

integer :: Ne, i, j, ii, jj, ix, iy, iz, ini
real(8), dimension(3) :: dr
character(len=4), dimension(:), allocatable :: AtomName_e, ResidName_e
integer, dimension(:), allocatable :: ResidNum_e
real(8), dimension(:,:), allocatable :: R_e, Vel_e

   Ne = N*ix_ext*iy_ext*iz_ext

   allocate(AtomName_e(N))
   allocate(ResidName_e(N))
   allocate(ResidNum_e(N))
   allocate(R_e(3,N))
   allocate(Vel_e(3,N))

   AtomName_e  = DefaultAtomName
   ResidName_e = DefaultResidName
   ResidNum_e  = ResidNum
   R_e         = R
   Vel_e       = Vel

   deallocate(DefaultAtomName, DefaultResidName)
   deallocate(ResidNum, R, Vel)

   allocate(DefaultAtomName(Ne))
   allocate(DefaultResidName(Ne))
   allocate(ResidNum(Ne))
   allocate(R(3,Ne))
   allocate(Vel(3,Ne))

   if((H(1,2)/=0.).or.(H(1,3)/=0.).or.(H(2,3)/=0.)) then
     print *, "ERROR: EXTENSION is not useful for non ortholonbic cell"
     stop
   end if

   N = Ne

   if(QThermostatNEW) then
     if((cThermostatMethodNEW == 'NHC').or. &
     &  (cThermostatMethodNEW == 'NH')) then
       InvMts(1) = InvMts(1)/(ix_ext*iy_ext*iz_ext)
     end if
   end if

   ii = 0
   jj = 0
   do i = 1, NumSpec
     if(i==1) then
       ini = 0
     else
       ini = ini + NumMol(i-1)*NumAtm(i-1)
     end if
     do ix = 1, ix_ext
       do iy = 1, iy_ext
         do iz = 1, iz_ext
           jj = jj + NumMol(i)
           dr(1) = (ix-1)*H(1,1) - H(1,1)*((ix_ext-1.d0)*0.5d0)
           dr(2) = (iy-1)*H(2,2) - H(2,2)*((iy_ext-1.d0)*0.5d0)
           dr(3) = (iz-1)*H(3,3) - H(3,3)*((iz_ext-1.d0)*0.5d0)
           do j = ini+1, ini+NumMol(i)*NumAtm(i)
             ii = ii + 1
             R(:,ii) = R_e(:,j) + dr(:)
             ResidNum(ii) = ResidNum_e(j)+(jj-NumMol(i))
             DefaultResidName(ii) = ResidName_e(j)
             DefaultAtomName(ii) = AtomName_e(j)
             Vel(:,ii) = Vel_e(:,j)
           end do
         end do
       end do
     end do
   end do

   QBathMomentRefresh = .True.

   H(1,1) = H(1,1) * ix_ext
   H(2,2) = H(2,2) * iy_ext
   H(3,3) = H(3,3) * iz_ext

   print '(/a,i9,a)', "The new system has ",N," atoms."
   print '(a)', "    Simulation Box Size is now "
   print '(5x,3f12.3)', H(1,:)
   print '(5x,3f12.3)', H(2,:)
   print '(5x,3f12.3/)', H(3,:)

   deallocate(AtomName_e, ResidName_e)
   deallocate(ResidNum_e, R_e, Vel_e)

end subroutine extend_sys


! ####################################################################
! ####################################################################


subroutine move_sys

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration
use CellParam, only : H
use RTR

implicit none

integer :: i

   do i = 1, N
     R(1,i) = R(1,i) + deltaX
     R(2,i) = R(2,i) + deltaY
     R(3,i) = R(3,i) + deltaZ
   end do

   call PBC

end subroutine move_sys


! ####################################################################
! ####################################################################


subroutine MDrestart

use Numbers, only : N, Nf
use CommonBlocks, only : QPBC, QResForm, cThermostatMethod
use Configuration, only : Vel
use RTR
use UnitExParam, only : kb, pi
use BathParam
use CellParam, only : H
use TimeParam, only : Timeps

implicit none

integer :: i, j, nMNHC
real(8) :: Scale, Gauss, Omeg, sv
external Gauss

     if(QChange_Constraint) then
       write(*,'(/a)') '-----------------------------------------------------------'
       write(*,'(a)')  ' Whenever you change the constraint method, this software  '
       write(*,'(a)')  ' does not work. You may use "MPDyn_Prep" instead!          '
       write(*,'(a/)') '-----------------------------------------------------------'
       stop
     end if

     if(QResForm) then
       open(1,file='restart.converted',form='formatted',status='unknown')
     else
       open(1,file='restart.converted',form='unformatted',status='unknown')
     end if

     call Write_Config(2)

     if(QChange_Temperature) then
       print *, ' Temperature will be changed by scaling the velocity '
       Scale = sqrt(Temp_NEW/Temp_o)
       do i = 1, N
         Vel(:,i) = Vel(:,i) * Scale
       end do
       kT = kb * Temp_NEW
       gkT  = Nf * kT
     end if

     call Write_Velocity

     if(QChange_Thermostat.or.QChange_ThermostatMethod) then

       if(QThermostatNEW) then

       if((cThermostatMethodNEW == 'NH').or.&
       &  (cThermostatMethodNEW == 'NHC')) then

         allocate( Rnh(NHchainNEW) )
         allocate( Vnh(NHchainNEW) )
         allocate( Mnh(NHchainNEW) )

         Rnh = 0.d0

         Omeg = 2.d0 * pi / Tau_NEW
         Mnh(1) = gkT / ( Omeg * Omeg )
         do i = 2, NHchainNEW
           Mnh(i) = kT / ( Omeg * Omeg )
         end do
         do i = 1, NHchainNEW
           sv = sqrt( kT / Mnh(i) )
           Vnh(i) = sv * Gauss()
         end do

         if(QResForm) then
           do i = 1 , NHchainNEW
             write(1,'(2d23.16)') Rnh(i), Vnh(i)
           end do
         else
           do i = 1 , NHchainNEW
             write(1) Rnh(i), Vnh(i)
           end do
         end if

       else if(cThermostatMethodNEW == 'MNHC') then

         Nf = N * 3 ! recover the exact number of degrees of freedom

         if(QBarostatNEW) then

           if( cBarostatMethodNEW == 'AN' ) then
             nMNHC = Nf + 1
           else if( cBarostatMethodNEW == 'A2' ) then
             nMNHC = Nf + 2
           else if( cBarostatMethodNEW == 'A3' ) then
             nMNHC = Nf + 3
           else if( (cBarostatMethodNEW == 'PR').or.( cBarostatMethodNEW == 'ST' ) ) then
             nMNHC = Nf + 6
           end if

         else

           nMNHC = Nf

         end if

         allocate( Rmnh(NHchainNEW,nMNHC) )
         allocate( Vmnh(NHchainNEW,nMNHC) )
         allocate( Mmnh(NHchainNEW,nMNHC) )

         Omeg = 2.d0 * pi / Tau_NEW

         do i = 1 , NHchainNEW
         do j = 1, nMNHC
           Mmnh(i,j) = kT / ( Omeg * Omeg )
         end do
         end do

         Rmnh = 0.d0

         do i = 1, nMNHC
           do j = 1, NHchainNEW
             sv = sqrt( kT / Mmnh(j,i) )
             Vmnh(j,i) = sv * Gauss()
           end do
         end do

         if(QResForm) then
           do i = 1 , nMNHC
             do j = 1 , NHchainNEW
               write(1,'(2d23.16)') Rmnh(j,i), Vmnh(j,i)
             end do
           end do
         else
           do i = 1 , nMNHC
             do j = 1 , NHchainNEW
               write(1) Rmnh(j,i), Vmnh(j,i)
             end do
           end do
         end if

       end if

       if(QPBC) then
         if(QChange_Barostat) Vg = 0.d0

         if(QResForm) then
           write(1,'(3d23.16)') (H(1,i), i = 1, 3)
           write(1,'(3d23.16)') (H(2,i), i = 1, 3)
           write(1,'(3d23.16)') (H(3,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(1,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(2,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(3,i), i = 1, 3)
         else
           write(1) H
           write(1) Vg
         end if
       end if

       end if

     else if(QBathMomentRefresh) then

       if((cThermostatMethod == 'NHC').or. &
       &  (cThermostatMethod == 'NH')) then
         Rss = 0.d0

         do i = 1, NHchain
           sv = sqrt( kT * InvMts(i) )
           Vss(i) = sv * Gauss()
         end do

         if(QResForm) then

           do i = 1 , NHchain
             write(1,'(2d23.16)') Rss(i), Vss(i)
           end do

         else

           do i = 1 , NHchain
             write(1) Rss(i), Vss(i)
           end do

         end if

       else if(cThermostatMethod == 'MNHC') then
         RMNHC = 0.d0

         do i = 1, NumMNHC
           do j = 1, NHchain
             sv = sqrt( kT * InvMMNHC(j,i) )
             VMNHC(j,i) = sv * Gauss()
           end do
         end do

         if(QResForm) then
           do i = 1 , NumMNHC
             do j = 1 , NHchain
               write(1,'(2d23.16)') RMNHC(j,i), VMNHC(j,i)
             end do
           end do
         else
           do i = 1 , NumMNHC
             do j = 1 , NHchain
               write(1) RMNHC(j,i), VMNHC(j,i)
             end do
           end do
         end if

       end if

       if(QPBC) then
         Vg = 0.d0

         if(QResForm) then
           write(1,'(3d23.16)') (H(1,i), i = 1, 3)
           write(1,'(3d23.16)') (H(2,i), i = 1, 3)
           write(1,'(3d23.16)') (H(3,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(1,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(2,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(3,i), i = 1, 3)
         else
           write(1) H
           write(1) Vg
         end if
       end if

     else

       if(QChange_Barostat) Vg = 0.d0

       call Write_Bath

     end if

     Timeps = Time_NEW

     call Write_Time

end subroutine MDrestart


! ####################################################################
! ####################################################################


subroutine PIMDrestart

use CommonBlocks, only : QResForm
use TimeParam, only : Timeps
use RTR

implicit none

     if(QResForm) then
       open(1,file='restart.converted',form='formatted',status='unknown')
     else
       open(1,file='restart.converted',form='unformatted',status='unknown')
     end if

     if(.not.QChange_Method) then
       call PIMD_PIMD
     else
       call MD_PIMD
     end if

     Timeps = Time_NEW

     call Write_Time

end subroutine PIMDrestart


! ####################################################################
! ####################################################################


subroutine PIMD_PIMD

use Numbers, only : N, Nf
use CommonBlocks, only : QPBC, QResForm, cThermostatMethod
use Configuration, only : R
use RTR
use CommonPI
use UnitExParam, only : Plank2pi, kb, pi2, pi
use BathParam
use CellParam, only : H
use AtomParam, only : MolName, ResidNum, DefaultAtomName, DefaultResidName, &
&   Mass

implicit none

integer :: i, j, nMNHC, k, l
character(len=10) :: ModelName
real(8) :: sv, vv, Omeg, xxx, QmassNEW, Gauss
external Gauss

     if((NbeadNEW == Nbead).and.(.not.QChange_Temperature)) then  ! no change in bead length 

       call Write_Config(0)
       call Write_Velocity

     else  ! bead length was changed 

       do i = 1 , N
         R(:,i) = Rnm(:,i,1)
       end do

       if(QChange_Temperature) then
         print *, ' Temperature will be changed '
         kT = kb * Temp_NEW
         gkT  = Nf * kT
         xxx  = sqrt(dble(NbeadNEW)) * kT / Plank2pi
         OmegaP2new = xxx * xxx
       end if

       if(QResForm) then

         ModelName = MolName(1)

         write(1,'(a/)') ModelName
         write(1,*) N

         if(N>100000) then
           do i = 1 , N
             write(1,'(2i9,2(x,a4),/3d24.16)')  &
             &     i , ResidNum(i) , DefaultResidName(i) ,  &
             &     DefaultAtomName(i) , R(:,i)
           end do
         else
           do i = 1 , N
             write(1,'(2i5,2(x,a4),/3d24.16)')  &
             &     i , ResidNum(i) , DefaultResidName(i) ,  &
             &     DefaultAtomName(i) , R(:,i)
           end do
         end if

       else

         write(1) N
         write(1) ResidNum
         write(1) DefaultResidName
         write(1) DefaultAtomName
         write(1) R

       end if

       allocate( RnmNEW(3,N,NbeadNEW) )
       allocate( VnmNEW(3,N,NbeadNEW) )
       allocate( MnmNEW(N,NbeadNEW) )
       allocate( FMnmNEW(N,NbeadNEW) )

       do j = 1, N

         MnmNEW(j,1)     = 0.d0
         MnmNEW(j,NbeadNEW) = 4.d0 * NbeadNEW * Mass(j)

         do i = 1, (NbeadNEW-2)/2

           MnmNEW(j,2*i)   = 2.d0 * (1.d0 - cos( pi2 * i / NbeadNEW )) * NbeadNEW * Mass(j)
           MnmNEW(j,2*i+1) = MnmNEW(j,2*i)

         end do

       end do

       do j = 2, NbeadNEW

         do i = 1, N

           sv = sqrt( kT / OmegaP2NEW / MnmNEW(i,j) )
           RnmNEW(1,i,j) = sv * Gauss()
           RnmNEW(2,i,j) = sv * Gauss()
           RnmNEW(3,i,j) = sv * Gauss()

         end do

       end do

! ## write 
       if(QResForm) then
         do j = 2, NbeadNEW
           do i = 1, N
             write(1,'(3d24.16)') RnmNEW(:,i,j)
           end do
         end do
       else
         do j = 2, NbeadNEW
           do i = 1, N
             write(1) RnmNEW(:,i,j)
           end do
         end do
       end if

! ## momentum for centers of mass  will be preserved 
       VnmNEW(:,:,1) = Vnm(:,:,1)

       do j = 1, N
         FMnmNEW(j,1) = Mass(j)
         do i = 2, NbeadNEW
           FMnmNEW(j,i) = GammaC2NEW * MnmNEW(j,i)
         end do
       end do

       do j = 2, NbeadNEW
!     /*  vsigma: standard devation of Maxwell distribution  */
         do i = 1, N
           vv = sqrt( kT / FMnmNEW(i,j) )
           VnmNEW(1,i,j) = vv * Gauss()
           VnmNEW(2,i,j) = vv * Gauss()
           VnmNEW(3,i,j) = vv * Gauss()
         end do
       end do

       if(QResForm) then

         do j = 1, Nbead
           do i = 1, N
             write(1,'(3d24.16)') VnmNEW(:,i,j)
           end do
         end do

       else

         write(1) VnmNEW

       end if

     end if

! ## Bath 

     if(QChange_Thermostat.or.QChange_ThermostatMethod) then

       if(QThermostatNEW) then

         if((cThermostatMethodNEW == 'NHC')) then
           allocate( Rnh(NHchainNEW) )
           allocate( Vnh(NHchainNEW) )
           allocate( Mnh(NHchainNEW) )

           Rnh = 0.d0
           Omeg = 2.d0 * pi / Tau_NEW
           Mnh(1) = gkT / ( Omeg * Omeg )

           do i = 2, NHchainNEW
             Mnh(i) = kT / ( Omeg * Omeg )
           end do
           do i = 1, NHchainNEW
             sv = sqrt( kT / Mnh(i) )
             Vnh(i) = sv * Gauss()
           end do
! ## Write 
           if(QResForm) then
             do i = 1 , NHchainNEW
               write(1,'(2d23.16)') Rnh(i), Vnh(i)
             end do
           else
             do i = 1 , NHchainNEW
               write(1) Rnh(i), Vnh(i)
             end do
           end if

         else if(cThermostatMethodNEW == 'MNHC') then

           if(QBarostatNEW) then

             if( cBarostatMethodNEW == 'AN' ) then
               nMNHC = Nf + 1
             else if( cBarostatMethodNEW == 'A2' ) then
               nMNHC = Nf + 2
             else if( cBarostatMethodNEW == 'A3' ) then
               nMNHC = Nf + 3
             else if( (cBarostatMethodNEW == 'PR').or.( cBarostatMethodNEW == 'ST' ) ) then
               nMNHC = Nf + 6
             end if

           else

             NumMNHC = Nf

           end if

           allocate( Rmnh(NHchainNEW,nMNHC) )
           allocate( Vmnh(NHchainNEW,nMNHC) )
           allocate( Mmnh(NHchainNEW,nMNHC) )

           Omeg = 2.d0 * pi / Tau_NEW

           do i = 1 , NHchainNEW
           do j = 1, nMNHC
             Mmnh(1,j) = kT / ( Omeg * Omeg )
           end do
           end do

           Rmnh = 0.d0

           do i = 1, nMNHC
             do j = 1, NHchainNEW
               sv = sqrt( kT / Mmnh(j,i) )
               Vmnh(j,i) = sv * Gauss()
             end do
           end do
! ## Write
           if(QResForm) then
             do i = 1 , nMNHC
               do j = 1 , NHchainNEW
                 write(1,'(2d23.16)') Rmnh(j,i), Vmnh(j,i)
               end do
             end do
           else
             do i = 1 , nMNHC
               do j = 1 , NHchainNEW
                 write(1) Rmnh(j,i), Vmnh(j,i)
               end do
             end do
           end if

         end if

       end if

! ##

       allocate( Vbb(3,N,NHchainNEW,2:NbeadNEW) )
       allocate( Rbb(3,N,NHchainNEW,2:NbeadNEW) )

       Rbb = 0.d0
       QmassNEW =  kT / OmegaP2NEW * GammaC2NEW

       do k = 2, Nbead
       do j = 1, NHchain
         do i = 1, N
           sv = sqrt( kT / QmassNEW )
           Vbb(1,i,j,k) = sv * Gauss()
           Vbb(2,i,j,k) = sv * Gauss()
           Vbb(3,i,j,k) = sv * Gauss()
         end do
       end do
       end do

! ## Write
       if(QResForm) then
         do k = 2, NbeadNEW
           do j = 1, NHchainNEW
             do i = 1, N
               do l = 1, 3
                 write(1,'(2d23.16)') Rbb(l,i,j,k), Vbb(l,i,j,k)
               end do
             end do
           end do
         end do
       else
         do k = 2, NbeadNEW
           do j = 1, NHchainNEW
             do i = 1, N
               do l = 1, 3
                 write(1) Rbb(l,i,j,k), Vbb(l,i,j,k)
               end do
             end do
           end do
         end do
       end if

       if(QPBC) then
         if(QChange_Barostat) Vg = 0.d0

         if(QResForm) then

           write(1,'(3d23.16)') (H(1,i), i = 1, 3)
           write(1,'(3d23.16)') (H(2,i), i = 1, 3)
           write(1,'(3d23.16)') (H(3,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(1,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(2,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(3,i), i = 1, 3)

         else

           write(1) H
           write(1) Vg

         end if

       end if

     else if(QBathMomentRefresh) then

       if((cThermostatMethod == 'NHC').or. &
       &  (cThermostatMethod == 'NH')) then

         Rss = 0.d0

         do i = 1, NHchain
           sv = sqrt( kT * InvMts(i) )
           Vss(i) = sv * Gauss()
         end do

         if(QResForm) then
           do i = 1 , NHchain
             write(1,'(2d23.16)') Rss(i), Vss(i)
           end do
         else
           do i = 1 , NHchain
             write(1) Rss(i), Vss(i)
           end do
         end if

       else if(cThermostatMethod == 'MNHC') then

         RMNHC = 0.d0

         do i = 1, NumMNHC
           do j = 1, NHchain
             sv = sqrt( kT * InvMMNHC(j,i) )
             VMNHC(j,i) = sv * Gauss()
           end do
         end do

         if(QResForm) then
           do i = 1 , NumMNHC
             do j = 1 , NHchain
               write(1,'(2d23.16)') RMNHC(j,i), VMNHC(j,i)
             end do
           end do
         else
           do i = 1 , NumMNHC
             do j = 1 , NHchain
               write(1) RMNHC(j,i), VMNHC(j,i)
             end do
           end do
         end if

       end if

       Rbath = 0.d0

       do k = 2, Nbead
         do j = 1, NHchain
           do i = 1, N
             sv = sqrt( kT * InvQmass(k) )
             Vbath(1,i,j,k) = sv * Gauss()
             Vbath(2,i,j,k) = sv * Gauss()
             Vbath(3,i,j,k) = sv * Gauss()
           end do
         end do
       end do

       if(QResForm) then
         do k = 2, Nbead
           do j = 1, NHchain
             do i = 1, N
               do l = 1, 3
                 write(1,'(2d23.16)') Rbath(l,i,j,k), Vbath(l,i,j,k)
               end do
             end do
           end do
         end do
       else
         do k = 2, Nbead
           do j = 1, NHchain
             do i = 1, N
               do l = 1, 3
                 write(1) Rbath(l,i,j,k), Vbath(l,i,j,k)
               end do
             end do
           end do
         end do
       end if

       if(QPBC) then
         Vg = 0.d0

         if(QResForm) then
           write(1,'(3d23.16)') (H(1,i), i = 1, 3)
           write(1,'(3d23.16)') (H(2,i), i = 1, 3)
           write(1,'(3d23.16)') (H(3,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(1,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(2,i), i = 1, 3)
           write(1,'(3d23.16)') (Vg(3,i), i = 1, 3)
         else
           write(1) H
           write(1) Vg
         end if
       end if

     else

       if(QChange_Barostat) Vg = 0.d0

       call Write_Bath

     end if

end subroutine PIMD_PIMD


! ####################################################################
! ####################################################################


subroutine MD_PIMD

use Numbers, only : N, Nf
use CommonBlocks, only : QPBC, QResForm
use Configuration
use RTR
use CommonPI
use UnitExParam, only : Plank2pi, kb, pi2, pi
use BathParam
use CellParam, only : H
use AtomParam, only : MolName, ResidNum, DefaultAtomName, DefaultResidName, &
&   Mass

implicit none

integer :: i, j, k, l, nMNHC
character(len=10) :: ModelName
real(8) :: sv, vv, Omeg, xxx, QmassNEW, Gauss
external Gauss

     if(QResForm) then

       ModelName = MolName(1)

       write(1,'(a/)') ModelName
       write(1,*) N

       if(N>100000) then
         do i = 1 , N
           write(1,'(2i9,2(x,a4),/3d24.16)')  &
           &     i , ResidNum(i) , DefaultResidName(i) ,  &
           &     DefaultAtomName(i) , R(:,i)
         end do
       else
         do i = 1 , N
           write(1,'(2i5,2(x,a4),/3d24.16)')  &
           &     i , ResidNum(i) , DefaultResidName(i) ,  &
           &     DefaultAtomName(i) , R(:,i)
         end do
       end if

     else

       write(1) N
       write(1) ResidNum
       write(1) DefaultResidName
       write(1) DefaultAtomName
       write(1) R

     end if

     if(QChange_Temperature) then
       print *, ' Temperature will be changed '
       kT = kb * Temp_NEW
       gkT  = Nf * kT
       xxx  = sqrt(dble(NbeadNEW)) * kT / Plank2pi
       OmegaP2new = xxx * xxx
     end if

     allocate( RnmNEW(3,N,NbeadNEW) )
     allocate( VnmNEW(3,N,NbeadNEW) )
     allocate( MnmNEW(N,NbeadNEW) )
     allocate( FMnmNEW(N,NbeadNEW) )

     do j = 1, N

       MnmNEW(j,1)     = 0.d0
       MnmNEW(j,NbeadNEW) = 4.d0 * NbeadNEW * Mass(j)

       do i = 1, (NbeadNEW-2)/2

         MnmNEW(j,2*i)   = 2.d0 * (1.d0 - cos( pi2 * i / NbeadNEW )) * NbeadNEW * Mass(j)
         MnmNEW(j,2*i+1) = MnmNEW(j,2*i)

       end do

     end do

     do j = 2, NbeadNEW
       do i = 1, N

         sv = sqrt( kT / OmegaP2NEW / MnmNEW(i,j) )
         RnmNEW(1,i,j) = sv * Gauss()
         RnmNEW(2,i,j) = sv * Gauss()
         RnmNEW(3,i,j) = sv * Gauss()

       end do
     end do

     if(QResForm) then
       do j = 2, NbeadNEW
         do i = 1, N
           write(1,'(3d24.16)') RnmNEW(:,i,j)
         end do
       end do
     else
       do j = 2, NbeadNEW
         do i = 1, N
           write(1) RnmNEW(:,i,j)
         end do
       end do
     end if

     VnmNEW(:,:,1) = Vel(:,:)

     do j = 1, N
       FMnmNEW(j,1) = Mass(j)
       do i = 2, NbeadNEW
         FMnmNEW(j,i) = GammaC2NEW * MnmNEW(j,i)
       end do
     end do

     do j = 2, NbeadNEW
!     /*  vsigma: standard devation of Maxwell distribution  */
       do i = 1, N
         vv = sqrt( kT / FMnmNEW(i,j) )
         VnmNEW(1,i,j) = vv * Gauss()
         VnmNEW(2,i,j) = vv * Gauss()
         VnmNEW(3,i,j) = vv * Gauss()
       end do
     end do

     if(QResForm) then
       do j = 1, Nbead
         do i = 1, N
           write(1,'(3d24.16)') VnmNEW(:,i,j)
         end do
       end do
     else
       write(1) VnmNEW
     end if


     if((cThermostatMethodNEW == 'NHC')) then
       allocate( Rnh(NHchainNEW) )
       allocate( Vnh(NHchainNEW) )
       allocate( Mnh(NHchainNEW) )
       Rnh = 0.d0
       Omeg = 2.d0 * pi / Tau_NEW
       Mnh(1) = gkT / ( Omeg * Omeg )
       do i = 2, NHchainNEW
         Mnh(i) = kT / ( Omeg * Omeg )
       end do
       do i = 1, NHchainNEW
         sv = sqrt( kT / Mnh(i) )
         Vnh(i) = sv * Gauss()
       end do
       if(QResForm) then
         do i = 1 , NHchainNEW
           write(1,'(2d23.16)') Rnh(i), Vnh(i)
         end do
       else
         do i = 1 , NHchainNEW
           write(1) Rnh(i), Vnh(i)
         end do
       end if

     else if(cThermostatMethodNEW == 'MNHC') then

       if(QBarostatNEW) then
         if( cBarostatMethodNEW == 'AN' ) then
           nMNHC = Nf + 1
         else if( cBarostatMethodNEW == 'A2' ) then
           nMNHC = Nf + 2
         else if( cBarostatMethodNEW == 'A3' ) then
           nMNHC = Nf + 3
         else if( (cBarostatMethodNEW == 'PR').or.( cBarostatMethodNEW == 'ST' ) ) then
           nMNHC = Nf + 6
         end if

       else

         NumMNHC = Nf

       end if

       allocate( Rmnh(NHchainNEW,nMNHC) )
       allocate( Vmnh(NHchainNEW,nMNHC) )
       allocate( Mmnh(NHchainNEW,nMNHC) )

       Omeg = 2.d0 * pi / Tau_NEW

       do i = 1 , NHchainNEW
       do j = 1, nMNHC
         Mmnh(1,j) = kT / ( Omeg * Omeg )
       end do
       end do

       Rmnh = 0.d0

       do i = 1, nMNHC
         do j = 1, NHchainNEW
           sv = sqrt( kT / Mmnh(j,i) )
           Vmnh(j,i) = sv * Gauss()
         end do
       end do

       if(QResForm) then
         do i = 1 , nMNHC
           do j = 1 , NHchainNEW
             write(1,'(2d23.16)') Rmnh(j,i), Vmnh(j,i)
           end do
         end do
       else
         do i = 1 , nMNHC
           do j = 1 , NHchainNEW
             write(1) Rmnh(j,i), Vmnh(j,i)
           end do
         end do
       end if

     end if

     allocate( Vbb(3,N,NHchainNEW,2:NbeadNEW) )
     allocate( Rbb(3,N,NHchainNEW,2:NbeadNEW) )

     Rbb = 0.d0
     QmassNEW =  kT / OmegaP2NEW * GammaC2NEW

     do k = 2, NbeadNEW
       do j = 1, NHchainNEW
         do i = 1, N
           sv = sqrt( kT / QmassNEW )
           Vbb(1,i,j,k) = sv * Gauss()
           Vbb(2,i,j,k) = sv * Gauss()
           Vbb(3,i,j,k) = sv * Gauss()
         end do
       end do
     end do

     if(QResForm) then

       do k = 2, NbeadNEW
         do j = 1, NHchainNEW
           do i = 1, N
             do l = 1, 3
               write(1,'(2d23.16)') Rbb(l,i,j,k), Vbb(l,i,j,k)
             end do
           end do
         end do
       end do

     else

       do k = 2, NbeadNEW
         do j = 1, NHchainNEW
           do i = 1, N
             do l = 1, 3
               write(1) Rbb(l,i,j,k), Vbb(l,i,j,k)
             end do
           end do
         end do
       end do

     end if

     if(QPBC) then
       Vg = 0.d0

       if(QResForm) then
         write(1,'(3d23.16)') (H(1,i), i = 1, 3)
         write(1,'(3d23.16)') (H(2,i), i = 1, 3)
         write(1,'(3d23.16)') (H(3,i), i = 1, 3)
         write(1,'(3d23.16)') (Vg(1,i), i = 1, 3)
         write(1,'(3d23.16)') (Vg(2,i), i = 1, 3)
         write(1,'(3d23.16)') (Vg(3,i), i = 1, 3)
       else
         write(1) H
         write(1) Vg
       end if
     end if

end subroutine MD_PIMD
