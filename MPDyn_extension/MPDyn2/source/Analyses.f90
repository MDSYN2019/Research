! ############################
! ## SUBROUTINE LIST 
! ## -- Select_Analyses 
! ## -- Read_RTraj 
! ## -- Read_VTraj 
! ## -- OpenTraj 
! ############################


! ***************************************************************
! **                                                           **
! **  MPDyn2                               For Analyses        **
! **                                        Ver. A.1           **
! **  Author --- Wataru Shinoda                                **
! **                                                           **
! **  Department of Applied Chemistry, Nagoya University       **
! **  Furo-cho, Chikusa-ku, Nagoya, Aichi 464-8603 Japan       **
! **  e-mail : w.shinoda@apchem.nagoya-u.ac.jp                 **
! **                                         July 20, 2014     **
! **                                                           **
! ***************************************************************


!######################################################################
!######################################################################


subroutine Select_Analyses

use ParamAnalyze
use Numbers, only : N

implicit NONE

real(8), dimension(3,N) :: Rorig

#ifdef MOLFILE
   call f77_molfile_init_
#endif

   if(cWhat=='CBPI') then

     call Cav_ParticleInsertion

   else if(cWhat=='TraceG') then

     call TraceG

   else if(cWhat=='ElecZ') then

     call ElecPot

   else if(cWhat=='DiffMSD') then

     call MSD_Diffusion

   else if(cWhat=='DiffMSDcom') then

     call MSD_DiffusionCom

   else if(cWhat=='DiffW_Z') then

     call Water_Diffusion_Zaxis

   else if(cWhat=='CorrRotW') then

     call CorrRotationW

   else if(cWhat=='RDF') then

     call GR

   else if(cWhat=='RDFG') then

     call GRG

   else if(cWhat=='LipidChAng') then

     call CoAngleChain

   else if(cWhat=='GR_LipidGG') then

     call GrGG_Lipid

   else if(cWhat=='Edens') then

     call ElecDensProfile

   else if(cWhat=='OriChain') then

     call OrientChain

   else if(cWhat=='OriLipid') then

     call OrientLipidSegments

   else if(cWhat=='RotLipid') then

     call RotLipid

   else if(cWhat=='LenChain') then

     call LengthChain

   else if(cWhat=='CoTGLipidC') then

     call CorrTG_LipidChain

   else if(cWhat=='dPP') then

     call dPP

   else if(cWhat=='SCD') then

     call SCD

   else if(cWhat=='FracGauche') then

     call FracGauche_Lipid

   else if( (cWhat=='PDB'  ).or.(cWhat=='PDBsp').or.(cWhat=='ARC').or.&
   &        (cWhat=='ARC_P').or.(cWhat=='XYZsp').or.(cWhat=='PDBex').or.&
   &        (cWhat=='DCD') ) then

     call MakePDB

   else if(cWhat=='PDBminR') then

     call MakePDBmin

   else if(cWhat=='CRD') then

     call MakeCRD

   else if(cWhat=='Cavity') then

     call CavityDist

   else if(cWhat=='RZ') then

     call RZ

   else if(cWhat=='RZG') then

     call RZG

   else if(cWhat=='RD') then

     call RD

   else if(cWhat=='RDdisk') then

     call RDdisk

   else if(cWhat=='RDdiskCOM') then

     call RDdiskG

   else if(cWhat=='RZCyl') then

     call RZCyl

   else if(cWhat=='RZRes') then

     call RZRes

   else if(cWhat=='RZResG') then

     call RZRes

   else if((cWhat=='WatOrient' ).or.&
   &       (cWhat=='WatOrientG').or.&
   &       (cWhat=='WatOriCyl' ).or.&
   &       (cWhat=='WatOriCylG')) then

     call OrientationOfWater

   else if(cWhat=='PoreWatTrj') then

     call PoreWaterTrajectroy

   else if(cWhat=='PoreWater') then

     call PoreWaterDensity

   else if(cWhat=='ResidTrj') then

     call ResidueTrajectroy

   else if(cWhat=='MSD_PmD_MS') then

     call AvConformProtein(Rorig)
     call MSD_MinimP_MS(Rorig)

   else if(cWhat=='MSD_P') then

     call MSD_protein

   else if(cWhat=='MSD_PD') then

     call MSD_proteinD

   else if(cWhat=='MSD_Pr') then

     call MSD_proteinR

   else if(cWhat=='MSD_PrD') then

     call MSD_proteinR_detail

   else if(cWhat=='MSD_Pm') then

     call MSD_proteinMinim

   else if(cWhat=='MSD_PmD') then

     call MSD_proteinMinim_detail

   else if(cWhat=='RofG_Prot') then

     call RadiusOfGyrationProtein

   else if(cWhat=='AreaOccL') then

     call AreaOccupL

   else if(cWhat=='VoronoiL') then

     call VoronoiLipid

   else if(cWhat=='VoronoiC') then

     call VoronoiLipid

   else if(cWhat=='LifeTimeW') then

     call Water_Exchange_Correlation

   else if(cWhat=='LS_anim') then

     call SingleLipidAnim

   else if(cWhat=='CGdata') then

     call CG_distributions

   else if(cWhat=='CGalkane') then

     call CoarseGrainC12

   else if(cWhat=='TGchain') then

     call TGchain

   else if(cWhat=='Neichain') then

     call NeighborChainCorr

   else if(cWhat=='PNcorr') then

     call PNcorr

   else if(cWhat=='CGW') then

     call CoarseGrainW

   else if(cWhat=='CGchainW') then

     call CG_Chain_W

   else if(cWhat=='DISTinMOL') then

     call IntraMolDist

   else if(cWhat=='ID_MICELLE') then

     call IDmicelle

   else if((cWhat=='Rot5Ring').or.&
   &       (cWhat=='RotVec')) then

     call RotRing5

   else if(cWhat=='SZZ') then

     call Szz_CG

   else if(cWhat=='INERTMOM') then

     call Anal_InerM

   else if(cWhat=='STRESSPROF') then

     call Ana_Stress

   else if(cWhat=='COHESIVE') then

     call CohesiveEne

   else if(cWhat=='COORDIS') then

     call CoorDis

#ifdef EnergyRep
   else if( cWhat == 'ENERGYHIST' ) then

     call Ana_ER

#endif
   else

     write(*,*) 'ERROR : no program was found for [',trim(cWhat),']'
     call Finalize

   end if

#ifdef MOLFILE
   call f77_molfile_finish_
#endif

end subroutine Select_Analyses


!######################################################################
!######################################################################


! *******************************
! * Instantaneous Configuration *
! *******************************
!
#ifdef MOLFILE
subroutine Read_RTraj(Nfile)
#else
subroutine Read_RTraj
#endif

use Configuration, only : R
use CommonBlocks, only : iTrjForm, QPBC
use ParamAnalyze
use IOparam, only : trajectory_file, NtrjF, ItrjF, StoreTrj
use TimeParam, only : Timeps
use Numbers, only : N
use CellParam, only : H
#ifdef MOLFILE
use MOLFILEparam
#endif

implicit none

integer :: i
character(len=1) :: Aname1
real, dimension(3,N) :: Rsn
real, dimension(3,3) :: Hsn
real :: Timesn
#ifdef MOLFILE
integer :: Nfile, j
integer :: Natom, status
character(len=10) :: Intype
real, dimension(3*N) :: XYZ
real, dimension(6) :: Box
real(8) :: hhx, hhy, hhz

   Intype = 'auto'

   ItrjF = ItrjF + 1

   if(ItrjF==1) then
     Natom = -1
     call f77_molfile_open_read_(Handle(Nfile), Natom, trajectory_file, Intype)
     if(Handle(Nfile)<0) then
        print*,'file type unknown or not registered'
!    else
!       print*,'file successfully opened:'
!       print*,'handle:',Handle(Nfile)
!       print*,'natom: ',Natom
     end if
     if(Natom /= N) then
       print *, 'ERROR : Natom in the trajectory file is not consistent'
       print *, 'Natom=',Natom
       print *, 'N=',N
       call Finalize
     end if
   end if

   status = 1
   call f77_molfile_read_next_(Handle(Nfile),Natom,XYZ(1),Box(1),status)

   if(ItrjF==NTrjStep(Nfile)) then
     call f77_molfile_close_read_(Handle(Nfile),status)
   end if

   do i = 1, N
     j = (i-1)*3
     R(1,i) = dble( XYZ(j+1) )
     R(2,i) = dble( XYZ(j+2) )
     R(3,i) = dble( XYZ(j+3) )
   end do

   H = 0.d0
   H(1,1) = dble(Box(1))
   H(2,2) = dble(Box(2))
   H(3,3) = dble(Box(3))

   if((Box(4)/=90.).or.(Box(5)/=90.).or.(Box(6)/=90.)) then
     print *, 'ERROR: Check box size'
     print *, Box(:)
     call Finalize
   end if

! ## move the whole system to have the origin at the center of box
! ## just for GROMACS trajectory
   if((iTrjForm==4).or.(iTrjForm==5)) then
     call PBC
   end if

#else

   ItrjF = ItrjF + 1

   if(iTrjForm==1) then

     if(ItrjF==1) open(21,file=trajectory_file,form='formatted',status='unknown')

     read(21,'(i9)') N
     if(QPBC) then
       read(21,'(10f9.4)') Timeps,H
     else
       read(21,'(f9.4)') Timeps
     end if

     do i = 1 , N
       read(21,'(a,3f9.3)') Aname1,Rsn(:,i)
     end do

     if(ItrjF==StoreTrj) then
       NtrjF = NtrjF + 1
       if(NtrjF<10) then
         write(trajectory_file,'(a,a,i1,a)') trim(adjustl(ReadJobName)),'.r00',NtrjF,'.xyz'
       else if(NtrjF<100) then
         write(trajectory_file,'(a,a,i2,a)') trim(adjustl(ReadJobName)),'.r0',NtrjF,'.xyz'
       else
         write(trajectory_file,'(a,a,i3,a)') trim(adjustl(ReadJobName)),'.r',NtrjF,'.xyz'
       end if
       close(21)
       ItrjF = 0
     end if

   else if(iTrjForm==2) then

     if(ItrjF==1) open(21,file=trajectory_file,form='unformatted',status='unknown')

     read(21) Timesn, N
     if(QPBC) read(21) Hsn
     read(21) Rsn

     if(ItrjF==StoreTrj) then
       NtrjF = NtrjF + 1
       if(NtrjF<10) then
         write(trajectory_file,'(a,a,i1)') trim(adjustl(ReadJobName)),'.r00',NtrjF
       else if(NtrjF<100) then
         write(trajectory_file,'(a,a,i2)') trim(adjustl(ReadJobName)),'.r0',NtrjF
       else
         write(trajectory_file,'(a,a,i3)') trim(adjustl(ReadJobName)),'.r',NtrjF
       end if
       close(21)
       ItrjF = 0
     end if

   else if(iTrjForm==3) then

     write(*,*) 'error : use MPDynS_MOLFILE for this analysis'
     call Finalize

   end if

   R = dble( Rsn )
   if(iTrjForm==2) then
     Timeps = dble(Timesn)
     H = dble( Hsn )
   end if
#endif

end subroutine Read_RTraj


!######################################################################
!######################################################################


! *******************************
! *   Instantaneous Velocity    *
! *******************************
!
subroutine Read_VTraj

use Configuration, only : Vel
use CommonBlocks, only : Job_name, iTrjForm
use IOparam, only : velocity_file, NvelF, IvelF, StoreVel
use CellParam, only : H
use TimeParam, only : Timeps
use Numbers, only : N

implicit none

real, dimension(3,N) :: Vsn
real, dimension(3,3) :: Hsn
real :: Timesn

   IvelF = IvelF + 1

   if(IvelF==1) then

     open(31,file=velocity_file,form='unformatted',status='unknown')

   end if

   read(31) Timesn, N
   read(31) Hsn
   read(31) Vsn

   if(IvelF==StoreVel) then

     NvelF = NvelF + 1

     if(NvelF<10) then

       write(velocity_file,'(a,a,i1)') trim(adjustl(Job_name)),'.v00',NvelF

     else if(NvelF<100) then

       write(velocity_file,'(a,a,i2)') trim(adjustl(Job_name)),'.v0',NvelF

     else

       write(velocity_file,'(a,a,i3)') trim(adjustl(Job_name)),'.v',NvelF

     end if

     close(31)

     IvelF = 0

   end if

   Vel = dble( Vsn )
   if(iTrjForm==2) then
     Timeps = dble(Timesn)
     H = dble( Hsn )
   end if

end subroutine Read_VTraj


!######################################################################
!######################################################################


subroutine OpenTraj(i)

use CommonBlocks, only : iTrjForm
use ParamAnalyze
use IOparam, only : trajectory_file, NtrjF, ItrjF

implicit none

integer :: i

   ItrjF = 0
   NtrjF = 1
   ReadJobName = FileHead(i)

#ifdef MOLFILE
   if(iTrjForm==3) then
     write(trajectory_file,'(a,a)') trim(adjustl(ReadJobName)),'.dcd'
   else if(iTrjForm==4) then
     write(trajectory_file,'(a,a)') trim(adjustl(ReadJobName)),'.trr'
   else if(iTrjForm==5) then
     write(trajectory_file,'(a,a)') trim(adjustl(ReadJobName)),'.xtc'
   else if(iTrjForm==1.or.iTrjForm==2) then
     write(*,*) 'error : use MPDynS for the analysis of this type of trajectory file'
     call Finalize
   end if
#else
   if(iTrjForm==1) then
     write(trajectory_file,'(a,a)') trim(adjustl(ReadJobName)),'.r001.xyz'
   else if(iTrjForm==2) then
     write(trajectory_file,'(a,a)') trim(adjustl(ReadJobName)),'.r001'
   else if(iTrjForm==3) then
     write(*,*) 'error : use MPDynS_MOLFILE for the analysis of this type of trajectory file'
     call Finalize
   end if
#endif

end subroutine OpenTraj
