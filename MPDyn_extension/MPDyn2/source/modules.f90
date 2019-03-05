

!---------------
! ## file names
!---------------
module IOparam
   character(len=80) :: DirectoryName
   character(len=80) :: PSF_file          ! charmm PSF file name
   character(len=80) :: trajectory_file   ! trajectory
   character(len=80) :: velocity_file     ! velocity
   character(len=80) :: Topology_file     ! topology file
   character(len=80) :: Parameter_file    ! parameter file
   character(len=80) :: Extmem_file       ! parameter file for external membrane potential
   character(len=80) :: Sphere_file       ! file for macrosphere vs atom interaction
   integer :: NtrjF, ItrjF, StoreTrj
   integer :: NvelF, IvelF, StoreVel
end module IOparam


module Configuration
   real(8), dimension(:,:), allocatable :: R         ! Coordinates
   real(8), dimension(:,:), allocatable :: Vel       ! Velocities
end module Configuration


module CommonBlocks
   character(len=50) :: Job_name
! --------------------
! ## operation flags
! --------------------
   logical :: QMaster
   logical :: QMinimization, QEps_r, QOption, QPLC, QREPWALL, QSTWALL, QCGWALL
   logical :: QPBC, QSHAKE, QRigidBody, QOpFix, QJarzynski, QOpJarz, QClust
   logical :: QInitial, QATaup, QHfunc, Qmixord, QLarge, QFrameMove, QFixCOM
   logical :: QBarostat, QThermostat, QPDB, QMacro, Qwcformol, QCyl, QFSCyl
   logical :: QPathInt, QFmonitor, QPINPT, QDelCellMove, QInert
   logical :: QSPCF, QHeatCap, QSwitch, QCorrectCutoff, QCorrCone
   logical :: QTICR, QMP, QCoulomb, QCellList, QGeneconf, QEflux, QEfield
   logical :: QResForm, QVolScale, Qstdout, Qdebug, QAveTh
   character(len=6) :: cThermostatMethod
   character(len=3) :: cBarostatMethod
   character(len=6) :: cCOULOMB
   character(len=10) :: ForceField
   character(len=10) :: SimMethod
   logical, dimension(:), allocatable :: PolyFlag
   integer :: iTrjForm
end module CommonBlocks


! ------------
! ## Numbers
! ------------
module Numbers
! N, Number: the number of atoms
! Nf: the number of degrees of freedom
! NfT : the number of translational degrees of freedom
! NfR : the number of rotational degrees of freedom
! NumSpec: the number of molecular species
! NumMol: the number of molecules ( for each species )
! NumMol: the number of the atoms in a molecule ( for each species)
! NumMer: the number of unit (for polymer)
   integer :: N, Nf, NfT, NfR
   integer :: Number
   integer :: NumSpec
   integer, dimension(:), allocatable :: NumMol, NumAtm
   integer, dimension(:), allocatable :: NumMer
end module Numbers


module ThermoData
! ## Temperature & Pressure
   real(8) :: Temp             ! calculated temperature
   real(8) :: Temp_Translation ! translational temperature
   real(8) :: Temp_Rotation    ! rotational temperature

   real(8) :: Ene_kin    ! Energy (kinetic)
   real(8) :: Ene_kinT   ! Energy (kinetic:translational)
   real(8) :: Ene_kinR   ! Energy (kinetic:rotational)

   real(8), dimension(3,3) :: Pkinp     ! kinetic part of pressure tensor
   real(8), dimension(3,3) :: Virial    ! Virial

   real(8), dimension(3) :: CellMomentum
end module ThermoData


module TimeParam
! ## Step Number & printout frequencies
   integer :: Nstep, itgn
   integer :: ixc, ixv, irs
! ## time step & multiple time step parameter
   integer :: lp, lk
   integer :: isampleHC
   integer :: BookFreq
   real(8) :: Timeps  ! simulation time
   real(8) :: deltat, dt2
end module TimeParam


module AtomParam
!  ## assigned potential parameters
   character(len=8), dimension(:), allocatable :: MolName  ! Name of molecules
   character(len=4), dimension(:), allocatable :: AtomName
   character(len=6), dimension(:), allocatable :: TypeName
   character(len=4), dimension(:), allocatable :: DefaultAtomName
   character(len=4), dimension(:), allocatable :: ResidName
   character(len=4), dimension(:), allocatable :: DefaultResidName
   integer, dimension(:), allocatable :: ResidNum
   real(8), dimension(:), allocatable :: Mass       ! mass
   real(8), dimension(:), allocatable :: InvMass    ! 1/mass
   integer :: NAAType
   character(len=6), dimension(:), allocatable :: AATypeList
   integer, dimension(:), allocatable :: NBAAType
end module AtomParam


module CutoffParam
! ## cut-off length for interaction
   real(8) :: Rcutoff2
   real(8) :: Ron2, swf1
   real(8) :: Rbook2
   real(8) :: Rskin
end module CutoffParam


module CellParam
! ## Cell Length, Volume
   real(8), dimension(3,3) :: H, InvH
   real(8), dimension(3) :: CellL, InvCL
   real(8) :: Volume
   integer :: CellShape
! ## Cell List for Virial calculation
   real(8), dimension(3,-13:13) :: CellShft
! ## Volume Scale
   real(8), dimension(3) :: Vsc_Rate
end module CellParam


! ## Cut-off Correction Parameter
module TailCorrect
   real(8) :: CorrectE, CorrectV
   real(8) :: Ene_LJ_co  ! Energy (cut-off correction for LJ)
   real(8), dimension(3,3) :: Virial_co ! Virial (Correction for Cutoff)
end module TailCorrect


module NonbondParam
! ## nonbonded parameters
   real(8), dimension(:), allocatable :: Charge     ! partial charge
   real(8), dimension(:), allocatable :: ScrnCR     ! screening shift parameter
   real(8), dimension(:), allocatable :: Rminh      ! Rmin/2  (LJ)
   real(8), dimension(:), allocatable :: EpsLJ      ! epsilon (LJ)
   real(8), dimension(:), allocatable :: SgmLJ      ! sigma (LJ)
   real(8), dimension(:), allocatable :: Eps14      ! epsilon (LJ 1-4)
   real(8), dimension(:), allocatable :: Rminh14    ! Rmin/2  (LJ 1-4)
! BKS parameters
   integer :: IJPair
   real(8), dimension(:), allocatable :: CoeA, CoeB, CoeC
   integer, dimension(:), allocatable :: BKStype
   integer, dimension(:,:), allocatable :: PairBack
! Force
   real(8), dimension(:,:), allocatable :: Frc_Elec  ! Forces (Electrostatic)
   real(8), dimension(:,:), allocatable :: Frc_Ersp  ! Forces (Nonbond Force R-Space)
   real(8), dimension(:,:), allocatable :: Frc_NBlong  ! Forces for CG models
   real(8), dimension(:,:), allocatable :: Frc_NBshrt  ! Forces for CG models

   real(8) :: Ene_LJ     ! Energy (Lennard-Jones or van der Waals)
   real(8) :: Ene_Elec   ! Energy (Electrostatic)
   real(8) :: Ene_Ersp   ! Energy (Ewald r-space)

   real(8) :: Ene_NBshrt  ! Energy 
   real(8) :: Ene_NBlong  ! Energy 
   real(8) :: Ene_ELshrt  ! Energy 
   real(8) :: Ene_ELlong  ! Energy 

   real(8), dimension(3,3) :: Vir_Ersp  ! Virial (NonBond R-Space)
   real(8), dimension(3,3) :: Vir_NBshrt ! Virial
   real(8), dimension(3,3) :: Vir_NBlong ! Virial
end module NonbondParam

! ## Intramolecular bond parameters
module BondedParam
   integer :: NumBond
   integer :: NumAngle
   integer :: NumUB
   integer :: NumDihedral
   integer :: NumImproper

   integer, dimension(:), allocatable :: BondI,BondJ,FTypeBond
   real(8), dimension(:), allocatable :: kBond,rBond
   real(8), dimension(:,:), allocatable :: Frc_Bond  ! Forces (Stretching)
   real(8), dimension(3,3) :: Vir_Bond  ! Virial (Bond Stretching)
   real(8) :: Ene_Bond   ! Energy (bond stretching)

   integer, dimension(:), allocatable :: AngleI,AngleJ,AngleK,FTypeAngle
   real(8), dimension(:), allocatable :: kTheta,Theta0
   logical, dimension(:), allocatable :: vdWSubtAng
   real(8), dimension(:,:), allocatable :: Frc_Angle ! Forces (Bending)
   real(8), dimension(3,3) :: Vir_Angle ! Virial (Angle Bending)
   real(8) :: Ene_Angle  ! Energy (angle bending)

   integer, dimension(:), allocatable :: UB_I,UB_J
   real(8), dimension(:), allocatable :: Kub,S0
   real(8), dimension(:,:), allocatable :: Frc_UB    ! Forces (Urey-Bradley)
   real(8), dimension(3,3) :: Vir_UB    ! Virial (Urey-Bradley)
   real(8) :: Ene_UB     ! Energy (Urey-Bradley)

   integer, dimension(:), allocatable :: DihedI,DihedJ,DihedK,DihedL,FTypeDihed
   real(8), dimension(:), allocatable :: kChi,DeltaDih,CsDelDih
   integer, dimension(:), allocatable :: NDih,DupFlag
   real(8), dimension(:,:), allocatable :: DupkChi,DupDeltaDih,DupCsDelDih
   integer, dimension(:,:), allocatable :: DupNDih
   logical, dimension(:), allocatable :: vdWSubtDih
   real(8), dimension(:), allocatable :: Ts      ! Tosion Angle List
   real(8), dimension(:,:), allocatable :: Frc_Dihed ! Forces (Torsion)
   real(8), dimension(3,3) :: Vir_Dihed ! Virial (Dihedral)
   real(8) :: Ene_Dihed  ! Energy (torsion)

   integer, dimension(:), allocatable :: ImproI,ImproJ,ImproK,ImproL,FTypeImpro
   real(8), dimension(:), allocatable :: kPsi,PsiImp
   real(8), dimension(:), allocatable :: kImp,DeltaImp,CsDelImp
   real(8), dimension(:,:), allocatable :: Frc_Impro ! Forces (Improper)
   real(8), dimension(3,3) :: Vir_Impro ! Virial (Improper)
   real(8) :: Ene_Impro  ! Energy (improper torsion)
   integer, dimension(:), allocatable :: NImp
end module BondedParam


! ----------------------------------
! ## Optional constraint parameters
! ----------------------------------
module OptConstraintParam
   integer :: NumOptC
   integer, dimension(:), allocatable :: OptCI,OptCJ
   real(8), dimension(:), allocatable :: kOptC,rOptC
   real(8), dimension(:,:), allocatable :: Frc_OptC  ! Forces
   real(8), dimension(3,3) :: Vir_OptC  ! Virial (Optional Constraint)
   real(8) :: Ene_OptC   ! Energy (optional constraint)

   integer :: NumPLC
   integer, dimension(:), allocatable :: PLCI
   real(8), dimension(:), allocatable :: kPLC, rPLC
   integer :: IaxisPLC

   integer :: NHam
   real(8) :: kHam
   real(8), dimension(:,:), allocatable :: RIni      ! Coordinates ( for Position Const. )
   real(8), dimension(:,:), allocatable :: Rrot      ! Coordinates ( for Position Const. )
   logical, dimension(:)  , allocatable :: HamR      ! Flag for the fixed particle

   integer :: Nccom
   real(8) :: kccom
   integer, dimension(:), allocatable :: Idcomf, Idcomt
   real(8), dimension(:), allocatable :: InvMscom
   real(8), dimension(:,:), allocatable :: Xcon
   real(8), dimension(:,:), allocatable :: Rcon
   real(8), dimension(:,:), allocatable :: fcom
   real(8), dimension(:,:), allocatable :: fscom

   real(8) :: Rad_Lipo, Thick_Lipo, Vol_Lipo
   real(8) :: facta, factb, factc, Invfacta
   real(8) :: ShiftCone

   integer :: NumClust, Ipartner
   integer, dimension(:), allocatable :: NAtomClus
   real(8) :: Rsh_clust, k_clust
end module OptConstraintParam

! -------------------
! ## Ewald Parameters
! -------------------
module EwaldParam
   integer :: Nh, Nel
   integer, dimension(:), allocatable :: Nelist
   real(8) :: Alpha
   real(8) :: alp2, ar2
   integer :: msh
   real(8), dimension(:), allocatable   :: EFList  ! error function list
   integer, dimension(:,:), allocatable :: ih
   real(8), dimension(:), allocatable   :: PCh
   real(8), dimension(:,:), allocatable :: Frc_Eksp  ! Forces (Ewald K-Space)
   real(8), dimension(3,3) :: Vir_Eksp  ! Virial (Ewald K-Space)
   real(8) :: Ene_Eksp   ! Energy (Ewald k-space)
   real(8) :: Ene_Eslf   ! Energy (Ewald self term)

   integer :: ih2mx, TmpNh, kmaxx, kmaxy, kmaxz
end module EwaldParam


module PMEparam

   integer :: Bsp_order
   integer :: SizeBtheta,SizeGridQ, MaxGrid
   integer :: NffTable, NffWork
   integer, dimension(3) :: Nfft, Nfftdim
   real(8), dimension(:,:), allocatable :: BthetaX, BthetaY, BthetaZ
   real(8), dimension(:,:), allocatable :: dBthetaX, dBthetaY, dBthetaZ
   real(8), dimension(:,:), allocatable :: ScRs
   real(8), dimension(:), allocatable :: BsplineModuleX, BsplineModuleY, BsplineModuleZ
   real(8), dimension(:), allocatable :: FFTtable
   real(8), dimension(:), allocatable :: FFTwork
! ## for parallel
   integer, dimension(:), allocatable :: Nscnt, Ndisp, Nrenum

#ifdef FFTW
   integer(8) :: planFx, planFy, planFz
   integer(8) :: planBx, planBy, planBz
   complex(8), dimension(:), allocatable :: Forward_in_x, Forward_out_x
   complex(8), dimension(:), allocatable :: Forward_in_y, Forward_out_y
   complex(8), dimension(:), allocatable :: Forward_in_z, Forward_out_z
   complex(8), dimension(:), allocatable :: Backward_in_x, Backward_out_x
   complex(8), dimension(:), allocatable :: Backward_in_y, Backward_out_y
   complex(8), dimension(:), allocatable :: Backward_in_z, Backward_out_z
#endif

end module PMEparam


module BathParam
! ## set points
   real(8) :: Temp_o, kT, gkT  ! set point
   real(8) :: Beta             ! = 1/kT 
   real(8) :: Pressure_o       ! Pressure
! -------------------
! ## Bath Parameters
! -------------------
   integer :: NHchain, NumMNHC
   integer :: NYoshid
   integer :: Nsc
   real(8) :: Tau_p, Tau_s0, Tau_s1
   real(8), dimension(3,3) :: prefTaup
   real(8), dimension(:), allocatable :: Mts, Rss, Vss
   real(8), dimension(:,:), allocatable :: RMNHC, VMNHC
   real(8), dimension(:,:), allocatable :: MMNHC, InvMMNHC
   real(8), dimension(:), allocatable :: InvMts
   real(8) :: Baro_kin
   real(8), dimension(3,3) :: Mp
   real(8), dimension(3,3) :: Vg
   real(8), dimension(:), allocatable :: wdti2, wdti4, wdti8
   real(8), dimension(3,3) :: InvMp
   real(8) :: InvNf
   integer, dimension(2) :: CoupleEdge
! ---------------------------------
! ## for constant stress ensemble
! ---------------------------------
   real(8), dimension(3,3) :: Stress_ext
   real(8), dimension(3,3) :: SigmaS
end module BathParam


module UnitExParam
! ------------------------
! ## Unit exchange factor
! ------------------------
   real(8), parameter :: reng = 1.d-04
   real(8), parameter :: rprs = 1.d-34
! ---------------------
! ## constants
! ---------------------
   real(8), parameter :: pi = 3.14159265358979d0
   real(8), parameter :: kb = 1.380662d-23 * reng
   real(8), parameter :: Avogadro = 6.022045d23
   real(8), parameter :: Plank2pi = 1.0545716d-34 * reng * 1.d12
! ## Unit exchange
   real(8), parameter :: ExParam = 4.184d+3 * reng / Avogadro
   real(8), parameter :: cvol = 1.d0 / ExParam
! ## Coulomb
   real(8), parameter :: e     = 1.602189d-19
   real(8), parameter :: eps_o = 8.8541878d-12
   real(8), parameter :: eps_r = 1.d0
   real(8), parameter :: ec    = (1.d0 / (4.*pi*eps_o*eps_r)) * e * e * (reng*1.d+10)

   real(8), parameter :: InvPi = 1.d0 / pi
   real(8), parameter :: sqpi = pi * pi
   real(8), parameter :: pi2 = 2.d0 * pi
end module UnitExParam

! -------------------
! ## List-Up
! -------------------
module BookParam
   integer :: MaxPair
   integer :: Npair
   integer, dimension(:,:), allocatable :: ListIJ
   integer :: Npair_short
   integer, dimension(:,:), allocatable :: List_shortIJ
end module BookParam

! -------------------------------------------
! ## constraint for Free energy calculations
! -------------------------------------------
module SMDparam
   integer :: NumFreqConst

   integer :: NumOpFixS
   integer, dimension(:), allocatable :: NumFixF, NumFixL
   integer, dimension(:), allocatable :: FixDir
   real(8), dimension(:), allocatable :: RGConst
   real(8), dimension(:), allocatable :: MassRG
   real(8), dimension(:), allocatable :: FFsm

   integer :: NumOpFixP
   integer, dimension(:), allocatable :: NumFixFi, NumFixLi
   integer, dimension(:), allocatable :: NumFixFj, NumFixLj
   real(8), dimension(:), allocatable :: ConstDis, ConstDis2
   real(8), dimension(:), allocatable :: MassRGi, MassRGj
   real(8), dimension(:,:), allocatable :: Rcomij

   real(8), dimension(3,3) :: Vir_ConstR ! Virial (Optional Constraint)
   real(8), dimension(3,3) :: Vir_ConstV ! Virial (Optional Constraint)

   integer :: NSelCyl, NSelCylI, NSelCone, icaxis, icx1, icx2
   integer, dimension(:), allocatable :: IdSelC
   real(8) :: FFcyl, FFcone

! ------------------------
! ## Jarzynski parameters
! ------------------------
   integer :: Nsample_steering
   real(8) :: k_steering, R0_steering, R0_initial
   real(8) :: Vshift_steering
   real(8), dimension(3) :: Uvec_steering
   real(8) :: Workforce
end module SMDparam


! --------------------------------------------
! ## wall constraint for a selected molecule
! --------------------------------------------
module wcparam
   real(8) :: kc_wall, InvMconst_atom
   real(8) :: dw_low, dw_upp, dw_mid
   real(8) :: rc_low, rc_upp
   integer :: atom_init, atom_fin, Iaxis_const, Nsample_wc, corder, Nmolwc
   integer, dimension(:), allocatable :: idf, ide
   real(8), dimension(:), allocatable :: invMwc
   real(8) :: Ronc
   logical :: Qwc_edge
end module wcparam


! --------------------
! ## SHAKE parameters 
! --------------------
module SHAKEparam

   real(8), dimension(3,3) :: Vir_Const ! Virial (Constraint)
   real(8), dimension(3,3) :: Vir_SHAKE ! Virial (Constraint SHAKE cubic cell)
   real(8), dimension(3,3) :: Vir_RATTL ! Virial (Constraint RATTLE cubic cell)

   real(8), dimension(:,:), allocatable :: R_o       ! Coordinates ( for SHAKE )
   real(8), dimension(:,:), allocatable :: Frc_Const ! Constraint Forces
   real(8), dimension(:,:), allocatable :: Lagmultip ! Lagrange multipliers

   integer :: NSHAKEGroup
   integer :: MaxHcBond
   integer, dimension(:),     allocatable :: NCoupleBond
   integer, dimension(:),     allocatable :: NCoupleAtom
   integer, dimension(:,:,:), allocatable :: CouplePair
   integer, dimension(:,:),   allocatable :: CoupleAtom
   real(8), dimension(:,:),   allocatable :: rSHAKE

   integer :: MaxIteration
   real(8) :: TolA
   real(8) :: TolB
end module SHAKEparam


! --------------------------
! ## Rigid Body parameters
! --------------------------
module RBparam
   integer :: NumRB
   integer, dimension(:), allocatable :: NumRBinMol
   integer, dimension(:), allocatable :: NumRBAtom
   integer :: NumRBType
   real(8), dimension(:,:), allocatable :: R_RB
   real(8), dimension(:,:), allocatable :: ScG
   real(8), dimension(:,:), allocatable :: V_RB
   real(8), dimension(:,:), allocatable :: F_RB
   real(8), dimension(:,:), allocatable :: Torque
   real(8), dimension(:,:), allocatable :: Quaternion
   real(8), dimension(:,:), allocatable :: Lmoment
   real(8), dimension(:,:), allocatable :: Omega

   integer :: MaxNumAtomRB
   real(8), dimension(:), allocatable :: MassRB
   real(8), dimension(:), allocatable :: InvMassRB
   real(8), dimension(:,:), allocatable :: InertiaRB
   real(8), dimension(:,:), allocatable :: InvInertiaRB
   integer, dimension(:), allocatable :: RBType
   logical, dimension(:), allocatable :: QLinear
   logical, dimension(:), allocatable :: QSingle

   real(8), dimension(:,:,:), allocatable :: R_onMol
   real(8), dimension(:,:,:), allocatable :: Rmolec
   real(8), dimension(:,:,:), allocatable :: Rotation

   integer, dimension(:), allocatable :: InitAtomRB
   integer, dimension(:), allocatable :: AtomUnitNum
   character(len=6), dimension(:), allocatable :: AtomUnitName
end module RBparam


module NoLJparam
   integer, dimension(:), allocatable :: NumNoLJ
   integer, dimension(:,:), allocatable :: NoLJ
end module NoLJparam


! ## Energy Minimization
module Eminparam
   integer :: MinTry
   real(8) :: dRmax
   real(8) :: dev_relative
end module Eminparam


! ## Wall parameter
module WallParam
   character(len=6) :: Ccap
   integer :: Lcap, nwall
   real(8), dimension(:), allocatable :: Rpcap
   integer, dimension(:), allocatable :: Ifcdir
   real(8) :: kcap,Rcap,Rpore
   real(8), dimension(3) :: Ocap
   integer :: NSelR
   integer, dimension(:), allocatable :: IdSelR

   integer :: nstwall, STdir, nlayer
   real(8) :: sg_wall, ep_wall, dlattice_wall, Delt_z
   real(8), dimension(:), allocatable :: Rstwall
   real(8), dimension(:), allocatable :: AijST, BijST
   real(8), dimension(:), allocatable :: FAijST, FBijST
   integer, dimension(:), allocatable :: FdirST

   integer :: ngrid_wall
   character(len=80) :: WallFile
   real(8), dimension(:), allocatable :: Rwmin, Rwmax
   real(8), dimension(:), allocatable :: invdz_size
   real(8), dimension(:,:,:), allocatable :: TabWall
end module WallParam


module CylParam
   integer :: icaxis, ngrid_cyl, ngrid_cylZ
   character(len=80) :: CylFile
   integer :: icx1, icx2
   real(8), dimension(:,:,:), allocatable :: TabCyl
   real(8), dimension(:,:,:,:,:), allocatable :: TabCylZ
   real(8), dimension(:,:), allocatable :: FscylR, FscylZ
   real(8), dimension(:), allocatable :: Rcylmin2, Rcylmax2
   real(8), dimension(:), allocatable :: Zcylmax
   real(8), dimension(:), allocatable :: invdr_size
   real(8) :: FCylRs, FCylH, PMFcyl, PMFcylH
end module CylParam


module ParamInertia
   integer :: NumInert
   integer, dimension(:), allocatable :: ListInert
   real(8), dimension(:), allocatable :: Msel
   real(8), dimension(:,:), allocatable :: Rsel
   real(8) :: TargetD, k_D, Dop, DPMF
end module ParamInertia


module EAM_param
   real(8), dimension(:,:), allocatable :: Frc_EAM
   real(8), dimension(3,3) :: Vir_EAM
   real(8) :: Ene_EAM
   integer, parameter :: Nref = 1000
   integer :: NumSpecPair
   real(8), parameter :: dh = 1.d-5
   integer, dimension(:,:), allocatable :: ResidPair
   real(8), dimension(:,:), allocatable :: Rrefa, Rhorefa
   real(8), dimension(:,:), allocatable :: Rhoy2, dRhodRref, dRhodRy2
   real(8), dimension(:),   allocatable :: Rmaxa, Rmina
   real(8), dimension(:,:), allocatable :: Rhorefb, FRhoref
   real(8), dimension(:,:), allocatable :: FRhoy2, dFdRhoref, dFdRhoy2
   real(8), dimension(:),   allocatable :: Rhomax, Rhomin
   real(8), dimension(:,:), allocatable :: Rrefb, PhiRef
   real(8), dimension(:,:), allocatable :: Phiy2, dPhidRref, dPhidRy2
   real(8), dimension(:),   allocatable :: Rmaxb, Rminb
#ifdef MEAM
   integer :: NumEAM, NumEAMSpec
#endif
end module EAM_param


! ## for Analyses 
module ParamAnalyze
   character(len=20) :: cWhat
   logical :: QConfig, QVeloc
   integer :: NJobs
   real(8) :: dtime

   character(len=50), dimension(:), allocatable :: FileHead   ! JobNames
   character(len=50) :: ReadJobName
   integer, dimension(:), allocatable :: NTrjStep, Interval
   character(len=50), dimension(:), allocatable :: ResultFile

   integer :: Nlipid
   integer :: NComp
   integer, dimension(:), allocatable :: NumComp
   logical :: Qregion
   real(8) :: Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
   integer :: Kcomp
   integer :: Nafiles, Nsnap

! ## MSD
   integer :: MSDdim
   integer, dimension(:), allocatable :: MSDdir

! ## radial distribution function
   logical :: QIntraSubt
   integer :: NumGR, IRcut
   character(len=4), dimension(:), allocatable :: GRAtomI, GRResiI
   character(len=4), dimension(:), allocatable :: GRAtomJ, GRResiJ
   character(len=40), dimension(:), allocatable :: GRfile
   integer, dimension(:,:), allocatable :: igr, ListGR_I, ListGR_J
   integer, dimension(:), allocatable :: NumGR_I, NumGR_J
   integer, dimension(:), allocatable :: NumGR_Subt, MoleculeID
   real(8) :: DRgrid

! ## distribution function of atoms along the bilayer normal
   integer :: NumRZ, RZdir
   character(len=4), dimension(:), allocatable :: RZAtom, RZResi
   character(len=30), dimension(:), allocatable :: RZfile
   character(len=30), dimension(:), allocatable :: RZfileI, RZfileO

! ## mean square displacement of atoms
   integer :: NumMSD
   character(len=4), dimension(:), allocatable :: MSDAtom, MSDResi
   character(len=8), dimension(:), allocatable :: MSDMol
   character(len=30), dimension(:), allocatable :: MSDfile

! ## cavity distribution analysis
   real(8) :: R_cav
   real(8) :: Rxmax, Rymax, Rzmax

! ## cavity biased particle insertion
   character(len=80) :: CBPI_FILE
   character(len=8), dimension(:), allocatable  :: MolNameCBPI
   integer :: NumInsertion, NumSpecInsert, NumOriTry
   integer, dimension(:), allocatable :: NumAtmInsert
   real(8), dimension(:), allocatable :: EneIns, EneInsert, EneInsertPT
!   real(8), dimension(:), allocatable :: EneInsLJ, EneInsertLJ
!   real(8), dimension(:), allocatable :: EneInsEl, EneInsertEl
   real(8), dimension(:,:), allocatable :: EpsLJPI, RminhPI, chargePI
   integer, dimension(:), allocatable :: ListI

! ## within a Cylinder 
   real(8) :: RadCyl, RadCyl2
   integer, dimension(1000) :: ResidFlag
   real(8) :: Zsh

   integer, dimension(:), allocatable :: ResRZ
   integer :: NumResRZ

   real(8), dimension(:,:), allocatable :: Rav

   integer, dimension(-40:40) :: NInCyl, NOutCyl
   integer :: NTotalInCyl, NTotalOutCyl
   real(8), dimension(-40:40,-20:20) :: WcsInCyl, WcsOutCyl

   integer, dimension(:), allocatable :: PoreWater

! ## electron density profile
   real(8), dimension(:), allocatable :: elecnum

   real(8) :: AreaL

! ## water diffusion coefficient in slab z+-dz 
   real(8) :: deltaZ
   integer :: ilength, interv, BlockStep
   integer :: maxdz, mindz

! ## 
   integer :: NBdAna, NAnAna, NDhAna
   real(8) :: DEGgrid
   character(len=4), dimension(:), allocatable :: BdAtomI, BdAtomJ
   character(len=4), dimension(:), allocatable :: AnAtomI, AnAtomJ, AnAtomK
   character(len=4), dimension(:), allocatable :: DhAtomI, DhAtomJ, DhAtomK, DhAtomL
   character(len=8), dimension(:), allocatable :: BdMol, AnMol, DhMol
   real(8) :: Rrangemin, Rrangemax

! ## stress profile
   logical :: QIK, QEW, Qslab, QLC
   integer :: Nbin, Nblock
   real(8) :: dZbin, RcutC2, dRbin, InvdRbin
   real(8), dimension(:), allocatable :: VirProXY,VirProZZ
   real(8), dimension(:), allocatable :: VirProN
   real(8), dimension(:), allocatable :: VirTProN
   real(8), dimension(:), allocatable :: VirTProT, VirProT
   real(8), dimension(:), allocatable :: Zbin
   integer, dimension(:), allocatable :: Ibin

! ## Cohesive Enery
   integer, dimension(:), allocatable :: MolID

! ## COORDIS
    character(len=6) :: MOLNAME1,MOLNAME2
    character(len=6), dimension(:), allocatable :: NatomI, NatomJ
    integer :: iatom, jatom, imol, jmol
    integer :: mex, nex, ngrid, ncgrid
    real(8) :: dAB, DNgrid

end module ParamAnalyze

! ------------------------------------------------------
! ## parameters in the CHARMM / OPLS / BKS force field
! ------------------------------------------------------
module FFParameters
   integer :: NumAtomTypeParam
   integer :: NumBondParam
   integer :: NumAngleParam
   integer :: NumDihedralParam
   integer :: NumImproperParam
   integer :: NumLJParam
   integer :: NumBKSParam

   character(len=6), dimension(2000) :: AtomTypeParam
   real(8), dimension(2000) :: MassParam

   character(len=6), dimension(2,3000) :: BondPairAtoms
   real(8), dimension(3000) :: kBondParam, rBondParam

   character(len=6), dimension(3,5000) :: AnglePairAtoms
   real(8), dimension(5000) :: kThetaParam, Theta0Param
! ## CHARMM
   real(8), dimension(5000) :: KubParam, S0Param
! ##

   character(len=6), dimension(4,7000) :: DihedralPairAtoms
   real(8), dimension(7000) :: kChiParam, DeltaDihParam
   integer, dimension(7000) :: NdihParam

! ## CHARMM
   character(len=6), dimension(4,700) :: ImproperPairAtoms
   real(8), dimension(700) :: kPsiParam, PsiImpParam
! ## OPLS
   real(8), dimension(700) :: kImpParam, DeltaImpParam
   integer, dimension(700) :: NImpParam
! ##

   character(len=6), dimension(3000) :: LJAtoms
   real(8), dimension(3000) :: EpsLJParam
! ## CHARMM
   real(8), dimension(3000) :: RminhParam,Eps14Param,Rminh14Param
   real, dimension(3000) :: ignoredParam,ignored14Param
! ## OPLS
   real(8), dimension(3000) :: SigmaLJParam
! ##

! ## BKS
   character(len=6), dimension(50) :: BKSAtomI, BKSAtomJ
   real(8), dimension(50) :: BKS_Aij, BKS_Bij, BKS_Cij

! ----------------------
! ## Links : Atom Name
! ----------------------
   integer :: NumResidueParam
   integer, dimension(200) :: NumAtom_inResi
   character(len=4), dimension(200) :: ResiNameParam
   character(len=4), dimension(400,200) :: AtomNameParam
   character(len=6), dimension(400,200) :: AtomType

! ## for RB UNIT
   integer, dimension(200) :: NumUNIT_inResi
   character(len=6), dimension(200,200) :: UNITNameParam
   integer, dimension(200,200) :: NumAtom_inUNITParam

! ## To make a topology data
   integer, dimension(200) :: NumBond_inResi
   integer, dimension(200) :: NumDoub_inResi
   integer, dimension(200) :: NumIMPR_inResi

   real(8), dimension(400,200) :: ChargeParam
   character(len=4), dimension(2,400,200) :: BondPair_inResi
   character(len=4), dimension(2,200,200) :: DoubPair_inResi
   character(len=4), dimension(4,200,200) :: ImprPair_inResi
! ##

end module FFParameters


module CGParameters

! ## parameters

   integer :: NumBondParam, NumAngleParam, NumDihedParam, NumImproParam
   integer :: NumNonBondParam

   character(len=6), dimension(:,:), allocatable :: BondPairAtoms
   integer, dimension(:), allocatable :: BondType
   real(8), dimension(:), allocatable :: kBondParam, rBondParam

   character(len=6), dimension(:,:), allocatable :: AnglePairAtoms
   integer, dimension(:), allocatable :: AngleType
   real(8), dimension(:), allocatable :: kThetaParam, Theta0Param

   character(len=6), dimension(:,:), allocatable :: DihedPairAtoms
   integer, dimension(:), allocatable :: DihedType
   real(8), dimension(:), allocatable :: kChiParam, DeltaDihParam
   integer, dimension(:), allocatable :: NdihParam

   character(len=6), dimension(:,:), allocatable :: ImproPairAtoms
   integer, dimension(:), allocatable :: ImproType
   real(8), dimension(:), allocatable :: kPsiParam, Psi0Param
   integer, dimension(:), allocatable :: NimpParam

   character(len=6), dimension(:,:), allocatable :: NonBondPairAtoms
   integer, dimension(:), allocatable :: NonBondType
   real(8), dimension(:), allocatable :: AijParam, BijParam, CijParam
   real(8), dimension(:), allocatable :: RminParam, RmaxParam
   real(8), dimension(:), allocatable :: R13cutParam, EpsParam
   character(len=13), dimension(:), allocatable :: NonBondFileHeader

! ## topology 

   integer :: NumMolParam

   integer, dimension(:), allocatable :: NumAtom_inMol
   character(len=8), dimension(:), allocatable :: MolNameParam
   character(len=4), dimension(:,:), allocatable :: AtomNameParam
   character(len=6), dimension(:,:), allocatable :: AtomType
   real(8), dimension(:,:), allocatable :: AtomMass, ChargeParam, ScrnParam

   integer, dimension(:), allocatable :: NumBond_inMol
   integer, dimension(:), allocatable :: NumIMPR_inMol
   integer, dimension(:,:,:), allocatable :: BondPair_inMol
   integer, dimension(:,:,:), allocatable :: ImprPair_inMol

   integer, dimension(:), allocatable :: NumUnit_inMol
   character(len=6), dimension(:,:), allocatable :: UNITNameParam
   integer, dimension(:,:), allocatable :: NumAtom_inUNITParam

   integer, dimension(:), allocatable :: NumResid_inMol
   character(len=4), dimension(:,:), allocatable :: ResidNameParam
   integer, dimension(:,:), allocatable :: NumAtom_inResidParam

end module CGParameters


module CommonMPI
   integer :: ierror
   integer :: MyRank, NProcs
   real(8), dimension(:), allocatable :: Buff1, Buff2
end module CommonMPI


module CommonDPD
   character(len=10) :: IntegrMethod
   integer :: Iequil, Iterate

   logical, dimension(:), allocatable :: ColloidFlag
   integer :: NumColloid, NumCollAtm
   logical :: QSheared, QColloid

   real(8), dimension(:,:), allocatable :: Velt
   real(8), dimension(:,:), allocatable :: V_RBt
   real(8), dimension(:,:), allocatable :: FrcDP, FrcDPt, FrcDPd
   real(8), dimension(:,:), allocatable :: FrcCo, FrcCot, FrcCod
   real(8), dimension(:,:), allocatable :: Torqued
   real(8), dimension(:,:), allocatable :: Lmomentt
   real(8), dimension(:,:), allocatable :: Omegat

   real(8), dimension(:,:), allocatable :: VelP
   real(8) :: ErrIter

! ---------------------
! ## model parameters
! ---------------------
   real(8), dimension(:,:), allocatable :: a
!   integer :: mfact, nfact
   real(8), dimension(:,:), allocatable :: e_morse
   real(8) :: a_morse
   real(8) :: sigma, gamma, sigmt, gammt
   real(8), dimension(:,:), allocatable :: GammDP ! for Peters only

! --------------
! ## time step
! --------------
   real(8) :: Lambda

! --------------
! ## shear rate
! --------------
   real(8) :: ShearRate, SlideGap

   integer, dimension(:), allocatable :: TypeNum

! --------------
! ## Numbers
! --------------
   integer, dimension(:), allocatable :: Ncal

! -------------------
! ## density
! -------------------
   real(8) :: dens

! -------------------
! ## Energy
! -------------------
   real(8) :: PotDP

! -----------
! ## average
! -----------
   real(8) :: Pressure_av, Vir_av, Pxy_av, Pr_xy_av
   real(8), dimension(3) :: Temperature_av
   real(8), dimension(3) :: Energy_av
   real(8), dimension(:), allocatable :: AveVslab
   integer :: Nslab
   real(8) :: ThSlab

   real(8), dimension(:,:), allocatable :: dR, SdR
   real(8), dimension(:,:), allocatable :: SdRc, dRc

   integer, dimension(:), allocatable :: irdf
   integer :: irc

! -----------
! ## Bond 
! -----------
   integer, dimension(:),     allocatable :: NumBondSpec
   integer, dimension(:,:,:), allocatable :: BondPairS
   real(8), dimension(:,:),   allocatable :: kBondSpec
   real(8), dimension(:,:),   allocatable :: rBondSpec

   integer, dimension(:,:), allocatable :: NpBond

! --------------------------------
! ## Self-Consistent calculation
! --------------------------------
   integer :: Ndm
   real(8), dimension(:)  , allocatable :: R1List, pfList, SLList
   real(8), dimension(:,:), allocatable :: dRList

! -----------------------------
! ## Check : temporary
! -----------------------------
   real(8), dimension(3,3) :: VirialC, VirialR, VirialD
   real(8) :: VirC_av, VirR_av, VirD_av

! ## H-function
   integer, dimension(:,:), allocatable :: Vlist

! ## mix order parameter
   real(8) :: sf1

end module CommonDPD


! ## Cell List Method
module CellListMethod
   integer, dimension(3) :: Ndiv
   integer :: Ncell
   integer :: Maps

   integer, dimension(:), allocatable :: Head
   integer, dimension(:), allocatable :: NextP
   integer, dimension(:), allocatable :: Map
   integer, dimension(:), allocatable :: MapShft
   integer, dimension(:), allocatable :: SelfShft
end module CellListMethod


module CommonHMC
  integer :: MDstep
  integer :: TimeMC
  integer :: NumAccept
  logical :: QPartial
  real(8) :: ThetaPartial
end module CommonHMC


! ## path integral MD
module CommonPI
   integer :: Nbead
   real(8) :: Pbead, InvP

   integer :: Ncref
   real(8) :: deltat_ref

   real(8), dimension(:,:,:), allocatable :: Rpi, Rnm, Vnm
   real(8), dimension(:,:,:), allocatable :: Fpi, Fnm, Fnm_ref
   real(8), dimension(:,:), allocatable :: Rcentroid
   real(8), dimension(:), allocatable :: QMass, InvQMass
   real(8), dimension(:,:), allocatable :: NmMass, FictMass, InvFictMass
   real(8), dimension(:,:), allocatable :: Tnm, InvTnm
   real(8) :: EnePI
   real(8) :: OmegaP2

   real(8), dimension(:,:,:,:), allocatable :: Rbath, Vbath, Fbath

! ## for CMD
   real(8) :: GammaC2

! ## for parallel calculation
   logical :: QBead, QMasterPI
   integer :: IniBead, FinBead
   integer :: NumProcess, BeadNum, MyRankPI
   integer, dimension(:), allocatable :: Nassi
   integer, dimension(:), allocatable :: BeadOrder
end module CommonPI


! >> F monitor ##
module F_monitor
   integer :: NgA
   real(8), dimension(:,:), allocatable :: Fint, Fext
   integer :: NumBondgA, NumAnglegA, NumDihedralgA, NumImpropergA, NumSHK

   integer, dimension(:), allocatable :: BondIgA, BondJgA
   integer, dimension(:), allocatable :: AngleIgA, AngleJgA, AngleKgA
   integer, dimension(:), allocatable :: DihedIgA, DihedJgA, DihedKgA, DihedLgA
   integer, dimension(:), allocatable :: SKpairI, SKpairJ
   logical, dimension(:), allocatable :: vdWSubtDihgA
   integer :: NlistgA
   integer, dimension(160000) :: ListIgA, ListJgA

   character(len=50) :: force_file   ! force
   character(len=50) :: confg_file   ! configuration
   integer :: NfrcF, IfrcF
   integer :: StoreFrc
   integer :: isampleF, NiniF, NfinF

end module F_monitor
! << F monitor ##


module QMDynamics
   integer :: Icurrent
   real(8), dimension(:,:), allocatable :: FrcQM
end module QMDynamics


module FEparam
   logical :: QCreation, QLJ, QCharge
   logical, dimension(:), allocatable :: QTICRcalc
   integer :: isampleTI, iTIstep
   integer, dimension(:), allocatable :: iequil_TICR, isample_TICR
   real(8), dimension(:), allocatable :: TIlambda
   real(8) :: Shift_Param
   integer :: NumTIpart, NiniTI, NfinTI
   integer :: NTIspec,NumTIMol
   real(8) :: E_TICR, AcAvE_TICR
   real(8), dimension(:), allocatable :: ChargeOrg
end module FEparam


module HeatCap
   real(8), dimension(:,:), allocatable :: Hessian

   real(8) :: acc_K_T, acc_K_V, acc_K_CV, acc_Vbar
   real(8) :: acc_E_T, acc_E_V, acc_E_CV, acc_E_T_E_T
   real(8) :: acc_E_T_E_V, acc_E_T_E_CV, acc_E_V_E_V
   real(8) :: acc_E_CV_E_CV, acc_dE_T_dbeta
   real(8) :: acc_B_V_E_V, acc_B_CV_E_CV
   real(8) :: acc_ROH(1000), acc_RHH(1000)
end module HeatCap


module MembrParam
   real(8), dimension(:), allocatable :: dEz, amem, bmem
   integer, dimension(:), allocatable :: Memptype
end module MembrParam


! ForCG 
module TableFuncs
   integer :: NumTables
   real(8), dimension(:), allocatable :: Invgs
   real(8), dimension(:,:,:), allocatable :: TabFunc
end module TableFuncs


! ForCG 
module CGdata
   integer :: NumAtype
   integer, dimension(:), allocatable :: NBAtomType
   character(len=6), dimension(:), allocatable :: AtomTypeList
   integer, dimension(:,:), allocatable :: NBFuncType
   integer, dimension(:,:), allocatable :: IDtable
   real(8), dimension(:,:), allocatable :: CoefAtype
   real(8), dimension(:,:), allocatable :: CoefBtype
   real(8), dimension(:,:), allocatable :: CoefCtype
   real(8), dimension(:,:), allocatable :: Rmin2, Rcut2, Rbk2, Rsw2, Swch, Rc13, Eps13
   real(8) :: Kappa ! Screening factor for electrostatic interactions
   real(8) :: Eps_relative ! Epsilon_r is not 1 as in AA simulation
   character(len=13), dimension(:), allocatable :: TabFileHead
   real(8) :: Rcut_MAX, InvRcut_MAX, Rbk_short2, Rcut_short2, Rrespa, Rheal, Rheal2
   real(8) :: Fsw1
   real(8) :: Ush0, Ush3, Ush4, Fsh2, Fsh3
end module CGdata


!------------------------
! ## for macroparticle
!------------------------
module CGball
   integer :: NumSphere ! the number of macroparticle
   integer :: NumTypeSphere ! the number of macroparticle types (how many different radius)
   real(8), dimension(:), allocatable :: RadiusSphere ! radius of macroparticle
   real(8) :: RhoSphere ! number density of atoms in the macroparticle
   real(8), dimension(:,:), allocatable :: Rcutmacro, Rcutmacro2, Rlstmacro2
   ! cutoff length between macroparticle and atoms
   real(8), dimension(:,:), allocatable :: Rminmacro, Rminmacro2
   ! minimum distance between macroparticle and atoms
   integer, dimension(:), allocatable :: IDsphere ! .true. for macroparticles
   integer :: MaxPairMac
   integer, dimension(:), allocatable :: NmacList
   integer, dimension(:,:,:), allocatable :: Listmac
   integer :: NumSphTab
   real(8), dimension(:,:), allocatable :: InvddSph
   real(8), dimension(:,:,:,:), allocatable :: SphFunc
   integer :: NMacSample
   real(8), dimension(:), allocatable :: PMFball, FSphRs
   integer, dimension(:), allocatable :: IdSph, ItypeSph
   integer :: Nox
   integer, dimension(:), allocatable :: IdOx
   real(8) :: eps_CO ! epsllon for SPCE oxygen vs. Colloid 
   integer :: MaxPairMM, NumPairMM
   integer, dimension(:,:), allocatable :: ListMM
   real(8), dimension(:,:), allocatable :: Rminmm2 ! square of minimum distance between macroparticles
   real(8), dimension(:,:), allocatable :: Rcutmm2 ! square of cutoff distance between macroparticles
   real(8), dimension(:,:), allocatable :: Rbkmm2  ! square of book-keeping distance between macroparticles
   real(8) :: sigMM, epsMM
   real(8), dimension(:,:), allocatable :: InvddMM
   real(8), dimension(:,:,:,:), allocatable :: UFtableMM
end module CGball

! -----------------------
! ## Simulated Annealing
! -----------------------
module SimAnneal
   logical :: QSimAnneal
   character(len=6) :: cAnnealMethod
   real(8) :: factA, factB, factC
   real(8) :: switch_time, T_target
   integer :: Nswitch
   real(8) :: kT_0, gkT_0
   real(8) :: TempSet
end module SimAnneal


!----------------------
! conductivity
!----------------------
module Conduct
   integer :: id_cation
   integer :: id_anion
   real(8) :: invMasscat, invMassani
   real(8), dimension(3) :: Eforce
   real(8) :: Efield
end module Conduct

!------------------------
! ## for script reading
!------------------------
module ScriptData
  integer :: NumOption, NumScript
  integer :: mfile
  logical :: Qcheck, QdescSW
  character(len=20), dimension(:), allocatable :: OPTIONS
  character(len=80), dimension(:), allocatable :: Script
  integer, dimension(:), allocatable :: Ist
  logical, dimension(:), allocatable :: ReadFlag
end module ScriptData


#ifdef FASTF
module invsqrt_table

integer(4), parameter :: fractmask = 8388607
integer(4), parameter :: expmask = 2139095040
integer(4), parameter :: explsb = 8388608
integer(4), parameter :: fractshift=12
integer(4), parameter :: expshift=23

integer(4), parameter, dimension(256) :: exptab =      &
&   (/z'5f000000',z'5e800000',z'5e800000',z'5e000000', &
&     z'5e000000',z'5d800000',z'5d800000',z'5d000000', &
&     z'5d000000',z'5c800000',z'5c800000',z'5c000000', &
&     z'5c000000',z'5b800000',z'5b800000',z'5b000000', &
&     z'5b000000',z'5a800000',z'5a800000',z'5a000000', &
&     z'5a000000',z'59800000',z'59800000',z'59000000', &
&     z'59000000',z'58800000',z'58800000',z'58000000', &
&     z'58000000',z'57800000',z'57800000',z'57000000', &
&     z'57000000',z'56800000',z'56800000',z'56000000', &
&     z'56000000',z'55800000',z'55800000',z'55000000', &
&     z'55000000',z'54800000',z'54800000',z'54000000', &
&     z'54000000',z'53800000',z'53800000',z'53000000', &
&     z'53000000',z'52800000',z'52800000',z'52000000', &
&     z'52000000',z'51800000',z'51800000',z'51000000', &
&     z'51000000',z'50800000',z'50800000',z'50000000', &
&     z'50000000',z'4f800000',z'4f800000',z'4f000000', &
&     z'4f000000',z'4e800000',z'4e800000',z'4e000000', &
&     z'4e000000',z'4d800000',z'4d800000',z'4d000000', &
&     z'4d000000',z'4c800000',z'4c800000',z'4c000000', &
&     z'4c000000',z'4b800000',z'4b800000',z'4b000000', &
&     z'4b000000',z'4a800000',z'4a800000',z'4a000000', &
&     z'4a000000',z'49800000',z'49800000',z'49000000', &
&     z'49000000',z'48800000',z'48800000',z'48000000', &
&     z'48000000',z'47800000',z'47800000',z'47000000', &
&     z'47000000',z'46800000',z'46800000',z'46000000', &
&     z'46000000',z'45800000',z'45800000',z'45000000', &
&     z'45000000',z'44800000',z'44800000',z'44000000', &
&     z'44000000',z'43800000',z'43800000',z'43000000', &
&     z'43000000',z'42800000',z'42800000',z'42000000', &
&     z'42000000',z'41800000',z'41800000',z'41000000', &
&     z'41000000',z'40800000',z'40800000',z'40000000', &
&     z'40000000',z'3f800000',z'3f800000',z'3f000000', &
&     z'3f000000',z'3e800000',z'3e800000',z'3e000000', &
&     z'3e000000',z'3d800000',z'3d800000',z'3d000000', &
&     z'3d000000',z'3c800000',z'3c800000',z'3c000000', &
&     z'3c000000',z'3b800000',z'3b800000',z'3b000000', &
&     z'3b000000',z'3a800000',z'3a800000',z'3a000000', &
&     z'3a000000',z'39800000',z'39800000',z'39000000', &
&     z'39000000',z'38800000',z'38800000',z'38000000', &
&     z'38000000',z'37800000',z'37800000',z'37000000', &
&     z'37000000',z'36800000',z'36800000',z'36000000', &
&     z'36000000',z'35800000',z'35800000',z'35000000', &
&     z'35000000',z'34800000',z'34800000',z'34000000', &
&     z'34000000',z'33800000',z'33800000',z'33000000', &
&     z'33000000',z'32800000',z'32800000',z'32000000', &
&     z'32000000',z'31800000',z'31800000',z'31000000', &
&     z'31000000',z'30800000',z'30800000',z'30000000', &
&     z'30000000',z'2f800000',z'2f800000',z'2f000000', &
&     z'2f000000',z'2e800000',z'2e800000',z'2e000000', &
&     z'2e000000',z'2d800000',z'2d800000',z'2d000000', &
&     z'2d000000',z'2c800000',z'2c800000',z'2c000000', &
&     z'2c000000',z'2b800000',z'2b800000',z'2b000000', &
&     z'2b000000',z'2a800000',z'2a800000',z'2a000000', &
&     z'2a000000',z'29800000',z'29800000',z'29000000', &
&     z'29000000',z'28800000',z'28800000',z'28000000', &
&     z'28000000',z'27800000',z'27800000',z'27000000', &
&     z'27000000',z'26800000',z'26800000',z'26000000', &
&     z'26000000',z'25800000',z'25800000',z'25000000', &
&     z'25000000',z'24800000',z'24800000',z'24000000', &
&     z'24000000',z'23800000',z'23800000',z'23000000', &
&     z'23000000',z'22800000',z'22800000',z'22000000', &
&     z'22000000',z'21800000',z'21800000',z'21000000', &
&     z'21000000',z'20800000',z'20800000',z'20000000', &
&     z'20000000',z'1f800000',z'1f800000',z'1f000000' /)


integer(4), parameter, dimension(4096) :: fractab =    &
& (/z'3504f3',z'34f9a4',z'34ee57',z'34e30c',z'34d7c3',z'34cc7c',z'34c137',z'34b5f5', &
&   z'34aab4',z'349f76',z'34943a',z'348900',z'347dc7',z'347291',z'34675e',z'345c2c', &
&   z'3450fc',z'3445ce',z'343aa3',z'342f79',z'342452',z'34192c',z'340e09',z'3402e8', &
&   z'33f7c9',z'33ecac',z'33e191',z'33d678',z'33cb61',z'33c04c',z'33b539',z'33aa28', &
&   z'339f19',z'33940d',z'338902',z'337df9',z'3372f3',z'3367ee',z'335cec',z'3351eb', &
&   z'3346ed',z'333bf0',z'3330f6',z'3325fd',z'331b07',z'331013',z'330520',z'32fa30', &
&   z'32ef41',z'32e455',z'32d96b',z'32ce82',z'32c39c',z'32b8b7',z'32add5',z'32a2f5', &
&   z'329816',z'328d3a',z'32825f',z'327787',z'326cb0',z'3261dc',z'325709',z'324c38', &
&   z'32416a',z'32369d',z'322bd2',z'32210a',z'321643',z'320b7e',z'3200bb',z'31f5fa', &
&   z'31eb3b',z'31e07e',z'31d5c3',z'31cb0a',z'31c053',z'31b59d',z'31aaea',z'31a038', &
&   z'319589',z'318adb',z'318030',z'317586',z'316ade',z'316038',z'315594',z'314af2', &
&   z'314052',z'3135b4',z'312b18',z'31207d',z'3115e5',z'310b4e',z'3100b9',z'30f627', &
&   z'30eb96',z'30e107',z'30d67a',z'30cbee',z'30c165',z'30b6dd',z'30ac58',z'30a1d4', &
&   z'309752',z'308cd2',z'308254',z'3077d8',z'306d5e',z'3062e5',z'30586e',z'304dfa', &
&   z'304387',z'303916',z'302ea7',z'302439',z'3019ce',z'300f64',z'3004fc',z'2ffa96', &
&   z'2ff032',z'2fe5d0',z'2fdb6f',z'2fd111',z'2fc6b4',z'2fbc59',z'2fb200',z'2fa7a9', &
&   z'2f9d53',z'2f9300',z'2f88ae',z'2f7e5e',z'2f7410',z'2f69c3',z'2f5f79',z'2f5530', &
&   z'2f4ae9',z'2f40a4',z'2f3661',z'2f2c1f',z'2f21df',z'2f17a1',z'2f0d65',z'2f032b', &
&   z'2ef8f2',z'2eeebc',z'2ee487',z'2eda53',z'2ed022',z'2ec5f2',z'2ebbc5',z'2eb199', &
&   z'2ea76e',z'2e9d46',z'2e931f',z'2e88fa',z'2e7ed7',z'2e74b5',z'2e6a96',z'2e6078', &
&   z'2e565c',z'2e4c41',z'2e4229',z'2e3812',z'2e2dfd',z'2e23e9',z'2e19d8',z'2e0fc8', &
&   z'2e05ba',z'2dfbad',z'2df1a3',z'2de79a',z'2ddd93',z'2dd38d',z'2dc989',z'2dbf87', &
&   z'2db587',z'2dab89',z'2da18c',z'2d9791',z'2d8d97',z'2d83a0',z'2d79aa',z'2d6fb6', &
&   z'2d65c3',z'2d5bd2',z'2d51e3',z'2d47f6',z'2d3e0a',z'2d3420',z'2d2a38',z'2d2051', &
&   z'2d166c',z'2d0c89',z'2d02a8',z'2cf8c8',z'2ceeea',z'2ce50d',z'2cdb33',z'2cd15a', &
&   z'2cc782',z'2cbdad',z'2cb3d9',z'2caa06',z'2ca036',z'2c9667',z'2c8c99',z'2c82ce', &
&   z'2c7904',z'2c6f3b',z'2c6575',z'2c5bb0',z'2c51ed',z'2c482b',z'2c3e6b',z'2c34ad', &
&   z'2c2af0',z'2c2135',z'2c177b',z'2c0dc4',z'2c040e',z'2bfa59',z'2bf0a6',z'2be6f5', &
&   z'2bdd46',z'2bd398',z'2bc9eb',z'2bc041',z'2bb698',z'2bacf0',z'2ba34b',z'2b99a6', &
&   z'2b9004',z'2b8663',z'2b7cc4',z'2b7326',z'2b698a',z'2b5ff0',z'2b5657',z'2b4cc0', &
&   z'2b432a',z'2b3996',z'2b3004',z'2b2673',z'2b1ce4',z'2b1357',z'2b09cb',z'2b0040', &
&   z'2af6b7',z'2aed30',z'2ae3ab',z'2ada27',z'2ad0a4',z'2ac724',z'2abda4',z'2ab427', &
&   z'2aaaab',z'2aa130',z'2a97b7',z'2a8e40',z'2a84ca',z'2a7b56',z'2a71e3',z'2a6872', &
&   z'2a5f03',z'2a5595',z'2a4c29',z'2a42be',z'2a3955',z'2a2fed',z'2a2687',z'2a1d23', &
&   z'2a13c0',z'2a0a5e',z'2a00fe',z'29f7a0',z'29ee43',z'29e4e8',z'29db8e',z'29d236', &
&   z'29c8e0',z'29bf8b',z'29b637',z'29ace5',z'29a395',z'299a46',z'2990f8',z'2987ad', &
&   z'297e62',z'297519',z'296bd2',z'29628c',z'295948',z'295005',z'2946c4',z'293d85', &
&   z'293446',z'292b0a',z'2921cf',z'291895',z'290f5d',z'290626',z'28fcf1',z'28f3be', &
&   z'28ea8c',z'28e15b',z'28d82c',z'28cefe',z'28c5d2',z'28bca8',z'28b37f',z'28aa57', &
&   z'28a131',z'28980c',z'288ee9',z'2885c7',z'287ca7',z'287389',z'286a6b',z'286150', &
&   z'285835',z'284f1c',z'284605',z'283cef',z'2833db',z'282ac8',z'2821b7',z'2818a7', &
&   z'280f98',z'28068b',z'27fd80',z'27f475',z'27eb6d',z'27e266',z'27d960',z'27d05c', &
&   z'27c759',z'27be57',z'27b557',z'27ac59',z'27a35c',z'279a60',z'279166',z'27886d', &
&   z'277f76',z'277680',z'276d8c',z'276499',z'275ba7',z'2752b7',z'2749c9',z'2740db', &
&   z'2737f0',z'272f05',z'27261c',z'271d35',z'27144f',z'270b6a',z'270287',z'26f9a5', &
&   z'26f0c4',z'26e7e5',z'26df08',z'26d62c',z'26cd51',z'26c477',z'26bba0',z'26b2c9', &
&   z'26a9f4',z'26a120',z'26984e',z'268f7d',z'2686ad',z'267ddf',z'267512',z'266c47', &
&   z'26637d',z'265ab4',z'2651ed',z'264927',z'264063',z'2637a0',z'262ede',z'26261e', &
&   z'261d5f',z'2614a2',z'260be6',z'26032b',z'25fa72',z'25f1ba',z'25e903',z'25e04e', &
&   z'25d79a',z'25cee7',z'25c636',z'25bd87',z'25b4d8',z'25ac2b',z'25a37f',z'259ad5', &
&   z'25922c',z'258985',z'2580de',z'257839',z'256f96',z'2566f4',z'255e53',z'2555b3', &
&   z'254d15',z'254479',z'253bdd',z'253343',z'252aaa',z'252213',z'25197d',z'2510e8', &
&   z'250855',z'24ffc3',z'24f732',z'24eea3',z'24e615',z'24dd88',z'24d4fc',z'24cc72', &
&   z'24c3ea',z'24bb62',z'24b2dc',z'24aa57',z'24a1d4',z'249952',z'2490d1',z'248852', &
&   z'247fd3',z'247756',z'246edb',z'246661',z'245de8',z'245570',z'244cfa',z'244485', &
&   z'243c11',z'24339f',z'242b2e',z'2422be',z'241a4f',z'2411e2',z'240976',z'24010c', &
&   z'23f8a2',z'23f03a',z'23e7d4',z'23df6e',z'23d70a',z'23cea7',z'23c646',z'23bde6', &
&   z'23b587',z'23ad29',z'23a4cc',z'239c71',z'239417',z'238bbf',z'238368',z'237b12', &
&   z'2372bd',z'236a69',z'236217',z'2359c6',z'235177',z'234928',z'2340db',z'23388f', &
&   z'233045',z'2327fb',z'231fb3',z'23176c',z'230f27',z'2306e2',z'22fe9f',z'22f65e', &
&   z'22ee1d',z'22e5de',z'22dda0',z'22d563',z'22cd28',z'22c4ed',z'22bcb4',z'22b47c', &
&   z'22ac46',z'22a411',z'229bdd',z'2293aa',z'228b78',z'228348',z'227b19',z'2272eb', &
&   z'226abe',z'226293',z'225a69',z'225240',z'224a18',z'2241f2',z'2239cc',z'2231a8', &
&   z'222985',z'222164',z'221944',z'221124',z'220907',z'2200ea',z'21f8ce',z'21f0b4', &
&   z'21e89b',z'21e083',z'21d86d',z'21d057',z'21c843',z'21c030',z'21b81e',z'21b00e', &
&   z'21a7fe',z'219ff0',z'2197e3',z'218fd8',z'2187cd',z'217fc4',z'2177bc',z'216fb5', &
&   z'2167af',z'215faa',z'2157a7',z'214fa5',z'2147a4',z'213fa4',z'2137a5',z'212fa8', &
&   z'2127ac',z'211fb1',z'2117b7',z'210fbe',z'2107c7',z'20ffd0',z'20f7db',z'20efe7', &
&   z'20e7f5',z'20e003',z'20d813',z'20d023',z'20c835',z'20c048',z'20b85d',z'20b072', &
&   z'20a889',z'20a0a1',z'2098ba',z'2090d4',z'2088ef',z'20810b',z'207929',z'207148', &
&   z'206968',z'206189',z'2059ab',z'2051cf',z'2049f3',z'204219',z'203a40',z'203268', &
&   z'202a91',z'2022bb',z'201ae7',z'201313',z'200b41',z'200370',z'1ffba0',z'1ff3d1', &
&   z'1fec04',z'1fe437',z'1fdc6c',z'1fd4a2',z'1fccd9',z'1fc511',z'1fbd4a',z'1fb584', &
&   z'1fadc0',z'1fa5fc',z'1f9e3a',z'1f9679',z'1f8eb9',z'1f86fa',z'1f7f3c',z'1f777f', &
&   z'1f6fc4',z'1f680a',z'1f6050',z'1f5898',z'1f50e1',z'1f492b',z'1f4176',z'1f39c3', &
&   z'1f3210',z'1f2a5f',z'1f22af',z'1f1aff',z'1f1351',z'1f0ba4',z'1f03f8',z'1efc4e', &
&   z'1ef4a4',z'1eecfb',z'1ee554',z'1eddae',z'1ed608',z'1ece64',z'1ec6c1',z'1ebf1f', &
&   z'1eb77f',z'1eafdf',z'1ea840',z'1ea0a3',z'1e9906',z'1e916b',z'1e89d1',z'1e8238', &
&   z'1e7aa0',z'1e7309',z'1e6b73',z'1e63de',z'1e5c4a',z'1e54b8',z'1e4d26',z'1e4596', &
&   z'1e3e06',z'1e3678',z'1e2eeb',z'1e275f',z'1e1fd4',z'1e184a',z'1e10c1',z'1e0939', &
&   z'1e01b3',z'1dfa2d',z'1df2a8',z'1deb25',z'1de3a2',z'1ddc21',z'1dd4a1',z'1dcd22', &
&   z'1dc5a3',z'1dbe26',z'1db6aa',z'1daf2f',z'1da7b6',z'1da03d',z'1d98c5',z'1d914e', &
&   z'1d89d9',z'1d8264',z'1d7af1',z'1d737e',z'1d6c0d',z'1d649c',z'1d5d2d',z'1d55bf', &
&   z'1d4e52',z'1d46e5',z'1d3f7a',z'1d3810',z'1d30a7',z'1d293f',z'1d21d8',z'1d1a73', &
&   z'1d130e',z'1d0baa',z'1d0447',z'1cfce6',z'1cf585',z'1cee25',z'1ce6c7',z'1cdf69', &
&   z'1cd80d',z'1cd0b1',z'1cc957',z'1cc1fe',z'1cbaa5',z'1cb34e',z'1cabf8',z'1ca4a2', &
&   z'1c9d4e',z'1c95fb',z'1c8ea9',z'1c8758',z'1c8008',z'1c78b8',z'1c716a',z'1c6a1d', &
&   z'1c62d1',z'1c5b86',z'1c543c',z'1c4cf3',z'1c45ab',z'1c3e65',z'1c371f',z'1c2fda', &
&   z'1c2896',z'1c2153',z'1c1a11',z'1c12d0',z'1c0b90',z'1c0452',z'1bfd14',z'1bf5d7', &
&   z'1bee9b',z'1be760',z'1be027',z'1bd8ee',z'1bd1b6',z'1bca7f',z'1bc349',z'1bbc15', &
&   z'1bb4e1',z'1badae',z'1ba67c',z'1b9f4c',z'1b981c',z'1b90ed',z'1b89bf',z'1b8292', &
&   z'1b7b67',z'1b743c',z'1b6d12',z'1b65e9',z'1b5ec1',z'1b579a',z'1b5074',z'1b4950', &
&   z'1b422c',z'1b3b09',z'1b33e7',z'1b2cc6',z'1b25a6',z'1b1e87',z'1b1769',z'1b104c', &
&   z'1b0930',z'1b0215',z'1afafb',z'1af3e2',z'1aecc9',z'1ae5b2',z'1ade9c',z'1ad787', &
&   z'1ad073',z'1ac95f',z'1ac24d',z'1abb3c',z'1ab42b',z'1aad1c',z'1aa60d',z'1a9f00', &
&   z'1a97f3',z'1a90e8',z'1a89dd',z'1a82d4',z'1a7bcb',z'1a74c3',z'1a6dbd',z'1a66b7', &
&   z'1a5fb2',z'1a58ae',z'1a51ab',z'1a4aa9',z'1a43a8',z'1a3ca8',z'1a35a9',z'1a2eab', &
&   z'1a27ae',z'1a20b1',z'1a19b6',z'1a12bc',z'1a0bc2',z'1a04ca',z'19fdd2',z'19f6dc', &
&   z'19efe6',z'19e8f2',z'19e1fe',z'19db0b',z'19d419',z'19cd28',z'19c638',z'19bf49', &
&   z'19b85b',z'19b16e',z'19aa82',z'19a396',z'199cac',z'1995c3',z'198eda',z'1987f3', &
&   z'19810c',z'197a26',z'197342',z'196c5e',z'19657b',z'195e99',z'1957b8',z'1950d8', &
&   z'1949f8',z'19431a',z'193c3d',z'193560',z'192e85',z'1927aa',z'1920d1',z'1919f8', &
&   z'191320',z'190c49',z'190573',z'18fe9e',z'18f7ca',z'18f0f7',z'18ea24',z'18e353', &
&   z'18dc82',z'18d5b3',z'18cee4',z'18c816',z'18c149',z'18ba7d',z'18b3b2',z'18ace8', &
&   z'18a61f',z'189f56',z'18988f',z'1891c8',z'188b03',z'18843e',z'187d7a',z'1876b7', &
&   z'186ff5',z'186934',z'186274',z'185bb4',z'1854f6',z'184e38',z'18477c',z'1840c0', &
&   z'183a05',z'18334b',z'182c92',z'1825da',z'181f23',z'18186c',z'1811b7',z'180b02', &
&   z'18044e',z'17fd9b',z'17f6e9',z'17f038',z'17e988',z'17e2d9',z'17dc2a',z'17d57d', &
&   z'17ced0',z'17c824',z'17c179',z'17bacf',z'17b426',z'17ad7e',z'17a6d6',z'17a030', &
&   z'17998a',z'1792e5',z'178c41',z'17859e',z'177efc',z'17785b',z'1771ba',z'176b1b', &
&   z'17647c',z'175dde',z'175741',z'1750a5',z'174a0a',z'17436f',z'173cd6',z'17363d', &
&   z'172fa5',z'17290f',z'172278',z'171be3',z'17154f',z'170ebb',z'170829',z'170197', &
&   z'16fb06',z'16f476',z'16ede7',z'16e759',z'16e0cb',z'16da3e',z'16d3b3',z'16cd28', &
&   z'16c69e',z'16c014',z'16b98c',z'16b305',z'16ac7e',z'16a5f8',z'169f73',z'1698ef', &
&   z'16926c',z'168be9',z'168568',z'167ee7',z'167867',z'1671e8',z'166b6a',z'1664ec', &
&   z'165e70',z'1657f4',z'165179',z'164aff',z'164486',z'163e0d',z'163796',z'16311f', &
&   z'162aa9',z'162434',z'161dc0',z'16174d',z'1610da',z'160a68',z'1603f8',z'15fd88', &
&   z'15f718',z'15f0aa',z'15ea3c',z'15e3d0',z'15dd64',z'15d6f9',z'15d08e',z'15ca25', &
&   z'15c3bc',z'15bd55',z'15b6ee',z'15b087',z'15aa22',z'15a3be',z'159d5a',z'1596f7', &
&   z'159095',z'158a34',z'1583d3',z'157d74',z'157715',z'1570b7',z'156a5a',z'1563fd', &
&   z'155da2',z'155747',z'1550ed',z'154a94',z'15443c',z'153de4',z'15378e',z'153138', &
&   z'152ae3',z'15248e',z'151e3b',z'1517e8',z'151197',z'150b45',z'1504f5',z'14fea6', &
&   z'14f857',z'14f209',z'14ebbc',z'14e570',z'14df25',z'14d8da',z'14d290',z'14cc47', &
&   z'14c5ff',z'14bfb7',z'14b971',z'14b32b',z'14ace6',z'14a6a1',z'14a05e',z'149a1b', &
&   z'1493d9',z'148d98',z'148758',z'148118',z'147ada',z'14749c',z'146e5f',z'146822', &
&   z'1461e7',z'145bac',z'145572',z'144f38',z'144900',z'1442c8',z'143c91',z'14365b', &
&   z'143026',z'1429f1',z'1423be',z'141d8b',z'141758',z'141127',z'140af6',z'1404c6', &
&   z'13fe97',z'13f869',z'13f23b',z'13ec0f',z'13e5e3',z'13dfb7',z'13d98d',z'13d363', &
&   z'13cd3a',z'13c712',z'13c0eb',z'13bac4',z'13b49e',z'13ae79',z'13a855',z'13a231', &
&   z'139c0e',z'1395ec',z'138fcb',z'1389ab',z'13838b',z'137d6c',z'13774e',z'137130', &
&   z'136b13',z'1364f8',z'135edc',z'1358c2',z'1352a8',z'134c8f',z'134677',z'134060', &
&   z'133a49',z'133433',z'132e1e',z'13280a',z'1321f6',z'131be3',z'1315d1',z'130fc0', &
&   z'1309af',z'13039f',z'12fd90',z'12f782',z'12f174',z'12eb67',z'12e55b',z'12df50', &
&   z'12d945',z'12d33b',z'12cd32',z'12c72a',z'12c122',z'12bb1b',z'12b515',z'12af10', &
&   z'12a90b',z'12a307',z'129d04',z'129702',z'129100',z'128aff',z'1284ff',z'127eff', &
&   z'127900',z'127302',z'126d05',z'126708',z'12610d',z'125b11',z'125517',z'124f1d', &
&   z'124925',z'12432c',z'123d35',z'12373e',z'123148',z'122b53',z'12255e',z'121f6b', &
&   z'121978',z'121385',z'120d94',z'1207a3',z'1201b3',z'11fbc3',z'11f5d4',z'11efe6', &
&   z'11e9f9',z'11e40d',z'11de21',z'11d836',z'11d24b',z'11cc62',z'11c679',z'11c090', &
&   z'11baa9',z'11b4c2',z'11aedc',z'11a8f7',z'11a312',z'119d2e',z'11974b',z'119168', &
&   z'118b87',z'1185a6',z'117fc5',z'1179e5',z'117407',z'116e28',z'11684b',z'11626e', &
&   z'115c92',z'1156b6',z'1150dc',z'114b02',z'114529',z'113f50',z'113978',z'1133a1', &
&   z'112dca',z'1127f5',z'112220',z'111c4b',z'111678',z'1110a5',z'110ad3',z'110501', &
&   z'10ff30',z'10f960',z'10f391',z'10edc2',z'10e7f4',z'10e226',z'10dc5a',z'10d68e', &
&   z'10d0c3',z'10caf8',z'10c52e',z'10bf65',z'10b99c',z'10b3d5',z'10ae0e',z'10a847', &
&   z'10a281',z'109cbc',z'1096f8',z'109134',z'108b72',z'1085af',z'107fee',z'107a2d', &
&   z'10746d',z'106ead',z'1068ee',z'106330',z'105d73',z'1057b6',z'1051fa',z'104c3e', &
&   z'104684',z'1040ca',z'103b10',z'103558',z'102fa0',z'1029e8',z'102432',z'101e7c', &
&   z'1018c6',z'101312',z'100d5e',z'1007ab',z'1001f8',z'ffc46',z'ff695',z'ff0e4', &
&   z'feb35',z'fe585',z'fdfd7',z'fda29',z'fd47c',z'fcecf',z'fc923',z'fc378', &
&   z'fbdce',z'fb824',z'fb27b',z'facd2',z'fa72a',z'fa183',z'f9bdd',z'f9637', &
&   z'f9092',z'f8aed',z'f854a',z'f7fa6',z'f7a04',z'f7462',z'f6ec1',z'f6920', &
&   z'f6381',z'f5de1',z'f5843',z'f52a5',z'f4d08',z'f476b',z'f41cf',z'f3c34', &
&   z'f369a',z'f3100',z'f2b66',z'f25ce',z'f2036',z'f1a9f',z'f1508',z'f0f72', &
&   z'f09dd',z'f0448',z'efeb4',z'ef921',z'ef38e',z'eedfc',z'ee86b',z'ee2da', &
&   z'edd4a',z'ed7ba',z'ed22b',z'ecc9d',z'ec710',z'ec183',z'ebbf7',z'eb66b', &
&   z'eb0e0',z'eab56',z'ea5cc',z'ea043',z'e9abb',z'e9533',z'e8fac',z'e8a26', &
&   z'e84a0',z'e7f1b',z'e7996',z'e7413',z'e6e8f',z'e690d',z'e638b',z'e5e0a', &
&   z'e5889',z'e5309',z'e4d8a',z'e480b',z'e428d',z'e3d0f',z'e3792',z'e3216', &
&   z'e2c9b',z'e2720',z'e21a5',z'e1c2c',z'e16b3',z'e113a',z'e0bc3',z'e064c', &
&   z'e00d5',z'dfb5f',z'df5ea',z'df075',z'deb01',z'de58e',z'de01b',z'ddaa9', &
&   z'dd538',z'dcfc7',z'dca57',z'dc4e7',z'dbf78',z'dba0a',z'db49c',z'daf2f', &
&   z'da9c2',z'da457',z'd9eeb',z'd9981',z'd9417',z'd8ead',z'd8945',z'd83dc', &
&   z'd7e75',z'd790e',z'd73a8',z'd6e42',z'd68dd',z'd6379',z'd5e15',z'd58b2', &
&   z'd534f',z'd4ded',z'd488c',z'd432b',z'd3dcb',z'd386c',z'd330d',z'd2dae', &
&   z'd2851',z'd22f4',z'd1d97',z'd183b',z'd12e0',z'd0d86',z'd082c',z'd02d2', &
&   z'cfd79',z'cf821',z'cf2ca',z'ced73',z'ce81c',z'ce2c7',z'cdd72',z'cd81d', &
&   z'cd2c9',z'ccd76',z'cc823',z'cc2d1',z'cbd7f',z'cb82f',z'cb2de',z'cad8f', &
&   z'ca83f',z'ca2f1',z'c9da3',z'c9856',z'c9309',z'c8dbd',z'c8871',z'c8326', &
&   z'c7ddc',z'c7892',z'c7349',z'c6e01',z'c68b9',z'c6372',z'c5e2b',z'c58e5', &
&   z'c539f',z'c4e5a',z'c4916',z'c43d2',z'c3e8f',z'c394c',z'c340a',z'c2ec9', &
&   z'c2988',z'c2448',z'c1f08',z'c19c9',z'c148b',z'c0f4d',z'c0a10',z'c04d3', &
&   z'bff97',z'bfa5b',z'bf521',z'befe6',z'beaad',z'be573',z'be03b',z'bdb03', &
&   z'bd5cb',z'bd095',z'bcb5e',z'bc629',z'bc0f4',z'bbbbf',z'bb68b',z'bb158', &
&   z'bac25',z'ba6f3',z'ba1c1',z'b9c90',z'b9760',z'b9230',z'b8d01',z'b87d2', &
&   z'b82a4',z'b7d76',z'b7849',z'b731d',z'b6df1',z'b68c6',z'b639b',z'b5e71', &
&   z'b5948',z'b541f',z'b4ef6',z'b49cf',z'b44a7',z'b3f81',z'b3a5b',z'b3535', &
&   z'b3010',z'b2aec',z'b25c8',z'b20a5',z'b1b82',z'b1660',z'b113e',z'b0c1d', &
&   z'b06fd',z'b01dd',z'afcbe',z'af79f',z'af281',z'aed64',z'ae847',z'ae32a', &
&   z'ade0e',z'ad8f3',z'ad3d8',z'acebe',z'ac9a4',z'ac48b',z'abf73',z'aba5b', &
&   z'ab544',z'ab02d',z'aab17',z'aa601',z'aa0ec',z'a9bd7',z'a96c3',z'a91b0', &
&   z'a8c9d',z'a878a',z'a8279',z'a7d67',z'a7857',z'a7347',z'a6e37',z'a6928', &
&   z'a641a',z'a5f0c',z'a59fe',z'a54f2',z'a4fe5',z'a4ada',z'a45ce',z'a40c4', &
&   z'a3bba',z'a36b0',z'a31a7',z'a2c9f',z'a2797',z'a2290',z'a1d89',z'a1883', &
&   z'a137d',z'a0e78',z'a0974',z'a0470',z'9ff6c',z'9fa69',z'9f567',z'9f065', &
&   z'9eb64',z'9e663',z'9e163',z'9dc63',z'9d764',z'9d266',z'9cd68',z'9c86a', &
&   z'9c36d',z'9be71',z'9b975',z'9b47a',z'9af7f',z'9aa85',z'9a58b',z'9a092', &
&   z'99b9a',z'996a1',z'991aa',z'98cb3',z'987bd',z'982c7',z'97dd1',z'978dc', &
&   z'973e8',z'96ef4',z'96a01',z'9650e',z'9601c',z'95b2b',z'9563a',z'95149', &
&   z'94c59',z'94769',z'9427a',z'93d8c',z'9389e',z'933b1',z'92ec4',z'929d8', &
&   z'924ec',z'92001',z'91b16',z'9162c',z'91142',z'90c59',z'90770',z'90288', &
&   z'8fda1',z'8f8ba',z'8f3d3',z'8eeed',z'8ea08',z'8e523',z'8e03e',z'8db5b', &
&   z'8d677',z'8d194',z'8ccb2',z'8c7d0',z'8c2ef',z'8be0e',z'8b92e',z'8b44e', &
&   z'8af6f',z'8aa91',z'8a5b2',z'8a0d5',z'89bf8',z'8971b',z'8923f',z'88d64', &
&   z'88889',z'883ae',z'87ed4',z'879fb',z'87522',z'87049',z'86b71',z'8669a', &
&   z'861c3',z'85ced',z'85817',z'85341',z'84e6d',z'84998',z'844c5',z'83ff1', &
&   z'83b1e',z'8364c',z'8317a',z'82ca9',z'827d8',z'82308',z'81e39',z'81969', &
&   z'8149b',z'80fcd',z'80aff',z'80632',z'80165',z'7fc99',z'7f7cd',z'7f302', &
&   z'7ee37',z'7e96d',z'7e4a4',z'7dfdb',z'7db12',z'7d64a',z'7d182',z'7ccbb', &
&   z'7c7f5',z'7c32f',z'7be69',z'7b9a4',z'7b4df',z'7b01b',z'7ab58',z'7a695', &
&   z'7a1d2',z'79d10',z'7984f',z'7938e',z'78ecd',z'78a0d',z'7854d',z'7808e', &
&   z'77bd0',z'77712',z'77254',z'76d97',z'768da',z'7641e',z'75f63',z'75aa8', &
&   z'755ed',z'75133',z'74c79',z'747c0',z'74308',z'73e50',z'73998',z'734e1', &
&   z'7302a',z'72b74',z'726be',z'72209',z'71d55',z'718a0',z'713ed',z'70f3a', &
&   z'70a87',z'705d5',z'70123',z'6fc72',z'6f7c1',z'6f311',z'6ee61',z'6e9b2', &
&   z'6e503',z'6e055',z'6dba7',z'6d6f9',z'6d24d',z'6cda0',z'6c8f4',z'6c449', &
&   z'6bf9e',z'6baf4',z'6b64a',z'6b1a0',z'6acf7',z'6a84f',z'6a3a7',z'69eff', &
&   z'69a58',z'695b2',z'6910c',z'68c66',z'687c1',z'6831d',z'67e78',z'679d5', &
&   z'67532',z'6708f',z'66bed',z'6674b',z'662aa',z'65e09',z'65969',z'654c9', &
&   z'65029',z'64b8a',z'646ec',z'6424e',z'63db1',z'63914',z'63477',z'62fdb', &
&   z'62b40',z'626a5',z'6220a',z'61d70',z'618d6',z'6143d',z'60fa4',z'60b0c', &
&   z'60674',z'601dd',z'5fd46',z'5f8b0',z'5f41a',z'5ef85',z'5eaf0',z'5e65b', &
&   z'5e1c7',z'5dd34',z'5d8a1',z'5d40e',z'5cf7c',z'5caea',z'5c659',z'5c1c9', &
&   z'5bd38',z'5b8a9',z'5b419',z'5af8a',z'5aafc',z'5a66e',z'5a1e1',z'59d54', &
&   z'598c7',z'5943b',z'58fb0',z'58b24',z'5869a',z'58210',z'57d86',z'578fd', &
&   z'57474',z'56feb',z'56b64',z'566dc',z'56255',z'55dcf',z'55949',z'554c3', &
&   z'5503e',z'54bb9',z'54735',z'542b1',z'53e2e',z'539ab',z'53529',z'530a7', &
&   z'52c25',z'527a4',z'52324',z'51ea4',z'51a24',z'515a5',z'51126',z'50ca8', &
&   z'5082a',z'503ad',z'4ff30',z'4fab4',z'4f638',z'4f1bc',z'4ed41',z'4e8c6', &
&   z'4e44c',z'4dfd3',z'4db59',z'4d6e0',z'4d268',z'4cdf0',z'4c979',z'4c502', &
&   z'4c08b',z'4bc15',z'4b79f',z'4b32a',z'4aeb5',z'4aa41',z'4a5cd',z'4a15a', &
&   z'49ce7',z'49874',z'49402',z'48f91',z'48b1f',z'486af',z'4823e',z'47dce', &
&   z'4795f',z'474f0',z'47082',z'46c14',z'467a6',z'46339',z'45ecc',z'45a60', &
&   z'455f4',z'45189',z'44d1e',z'448b3',z'44449',z'43fdf',z'43b76',z'4370d', &
&   z'432a5',z'42e3d',z'429d6',z'4256f',z'42108',z'41ca2',z'4183c',z'413d7', &
&   z'40f72',z'40b0e',z'406aa',z'40247',z'3fde4',z'3f981',z'3f51f',z'3f0bd', &
&   z'3ec5c',z'3e7fb',z'3e39b',z'3df3b',z'3dadb',z'3d67c',z'3d21d',z'3cdbf', &
&   z'3c961',z'3c504',z'3c0a7',z'3bc4a',z'3b7ee',z'3b393',z'3af37',z'3aadd', &
&   z'3a682',z'3a228',z'39dcf',z'39976',z'3951d',z'390c5',z'38c6d',z'38816', &
&   z'383bf',z'37f69',z'37b13',z'376bd',z'37268',z'36e13',z'369bf',z'3656b', &
&   z'36117',z'35cc4',z'35872',z'3541f',z'34fce',z'34b7c',z'3472b',z'342db', &
&   z'33e8b',z'33a3b',z'335ec',z'3319d',z'32d4f',z'32901',z'324b3',z'32066', &
&   z'31c1a',z'317cd',z'31381',z'30f36',z'30aeb',z'306a1',z'30256',z'2fe0d', &
&   z'2f9c3',z'2f57a',z'2f132',z'2ecea',z'2e8a2',z'2e45b',z'2e014',z'2dbce', &
&   z'2d788',z'2d343',z'2cefd',z'2cab9',z'2c675',z'2c231',z'2bded',z'2b9aa', &
&   z'2b568',z'2b125',z'2ace4',z'2a8a2',z'2a461',z'2a021',z'29be1',z'297a1', &
&   z'29362',z'28f23',z'28ae4',z'286a6',z'28269',z'27e2c',z'279ef',z'275b2', &
&   z'27176',z'26d3b',z'26900',z'264c5',z'2608b',z'25c51',z'25817',z'253de', &
&   z'24fa6',z'24b6d',z'24735',z'242fe',z'23ec7',z'23a90',z'2365a',z'23224', &
&   z'22def',z'229ba',z'22585',z'22151',z'21d1d',z'218ea',z'214b7',z'21084', &
&   z'20c52',z'20821',z'203ef',z'1ffbe',z'1fb8e',z'1f75e',z'1f32e',z'1eeff', &
&   z'1ead0',z'1e6a1',z'1e273',z'1de45',z'1da18',z'1d5eb',z'1d1bf',z'1cd93', &
&   z'1c967',z'1c53c',z'1c111',z'1bce6',z'1b8bc',z'1b493',z'1b069',z'1ac40', &
&   z'1a818',z'1a3f0',z'19fc8',z'19ba1',z'1977a',z'19354',z'18f2d',z'18b08', &
&   z'186e2',z'182be',z'17e99',z'17a75',z'17651',z'1722e',z'16e0b',z'169e9', &
&   z'165c6',z'161a5',z'15d83',z'15963',z'15542',z'15122',z'14d02',z'148e3', &
&   z'144c4',z'140a5',z'13c87',z'13869',z'1344c',z'1302f',z'12c12',z'127f6', &
&   z'123da',z'11fbf',z'11ba4',z'11789',z'1136f',z'10f55',z'10b3c',z'10723', &
&   z'1030a',z'fef2',z'fada',z'f6c2',z'f2ab',z'ee95',z'ea7e',z'e668', &
&   z'e253',z'de3e',z'da29',z'd614',z'd200',z'cded',z'c9da',z'c5c7', &
&   z'c1b4',z'bda2',z'b990',z'b57f',z'b16e',z'ad5e',z'a94e',z'a53e', &
&   z'a12e',z'9d1f',z'9911',z'9503',z'90f5',z'8ce7',z'88da',z'84ce', &
&   z'80c1',z'7cb5',z'78aa',z'749f',z'7094',z'6c89',z'687f',z'6476', &
&   z'606d',z'5c64',z'585b',z'5453',z'504b',z'4c44',z'483d',z'4436', &
&   z'4030',z'3c2a',z'3825',z'3420',z'301b',z'2c17',z'2813',z'240f', &
&   z'200c',z'1c09',z'1807',z'1405',z'1003',z'c02',z'801',z'400', &
&   z'7fffff',z'7ff001',z'7fe006',z'7fd00d',z'7fc018',z'7fb025',z'7fa036',z'7f9049', &
&   z'7f8060',z'7f7079',z'7f6095',z'7f50b5',z'7f40d7',z'7f30fc',z'7f2124',z'7f114f', &
&   z'7f017e',z'7ef1af',z'7ee1e2',z'7ed219',z'7ec253',z'7eb290',z'7ea2d0',z'7e9312', &
&   z'7e8358',z'7e73a0',z'7e63eb',z'7e543a',z'7e448b',z'7e34df',z'7e2536',z'7e1590', &
&   z'7e05ec',z'7df64c',z'7de6ae',z'7dd714',z'7dc77c',z'7db7e7',z'7da855',z'7d98c6', &
&   z'7d893a',z'7d79b0',z'7d6a2a',z'7d5aa6',z'7d4b25',z'7d3ba7',z'7d2c2c',z'7d1cb3', &
&   z'7d0d3e',z'7cfdcb',z'7cee5b',z'7cdeee',z'7ccf84',z'7cc01d',z'7cb0b8',z'7ca156', &
&   z'7c91f7',z'7c829b',z'7c7342',z'7c63eb',z'7c5497',z'7c4546',z'7c35f8',z'7c26ad', &
&   z'7c1764',z'7c081e',z'7bf8db',z'7be99b',z'7bda5d',z'7bcb23',z'7bbbeb',z'7bacb5', &
&   z'7b9d83',z'7b8e53',z'7b7f26',z'7b6ffc',z'7b60d4',z'7b51b0',z'7b428e',z'7b336e', &
&   z'7b2452',z'7b1538',z'7b0621',z'7af70c',z'7ae7fb',z'7ad8ec',z'7ac9e0',z'7abad6', &
&   z'7aabcf',z'7a9ccb',z'7a8dca',z'7a7ecb',z'7a6fcf',z'7a60d5',z'7a51df',z'7a42eb', &
&   z'7a33f9',z'7a250b',z'7a161f',z'7a0735',z'79f84f',z'79e96b',z'79da89',z'79cbab', &
&   z'79bccf',z'79adf5',z'799f1f',z'79904a',z'798179',z'7972aa',z'7963de',z'795515', &
&   z'79464e',z'793789',z'7928c8',z'791a09',z'790b4c',z'78fc92',z'78eddb',z'78df27', &
&   z'78d075',z'78c1c5',z'78b319',z'78a46e',z'7895c7',z'788722',z'78787f',z'7869e0', &
&   z'785b42',z'784ca8',z'783e10',z'782f7a',z'7820e7',z'781257',z'7803c9',z'77f53e', &
&   z'77e6b5',z'77d82f',z'77c9ab',z'77bb2a',z'77acac',z'779e30',z'778fb6',z'77813f', &
&   z'7772cb',z'776459',z'7755ea',z'77477d',z'773913',z'772aab',z'771c46',z'770de3', &
&   z'76ff83',z'76f125',z'76e2ca',z'76d472',z'76c61b',z'76b7c8',z'76a977',z'769b28', &
&   z'768cdc',z'767e92',z'76704b',z'766206',z'7653c4',z'764584',z'763747',z'76290c', &
&   z'761ad3',z'760c9d',z'75fe6a',z'75f039',z'75e20a',z'75d3de',z'75c5b5',z'75b78e', &
&   z'75a969',z'759b46',z'758d27',z'757f09',z'7570ee',z'7562d6',z'7554bf',z'7546ac', &
&   z'75389a',z'752a8c',z'751c7f',z'750e75',z'75006d',z'74f268',z'74e465',z'74d665', &
&   z'74c867',z'74ba6b',z'74ac72',z'749e7b',z'749087',z'748295',z'7474a5',z'7466b8', &
&   z'7458cd',z'744ae4',z'743cfe',z'742f1a',z'742139',z'74135a',z'74057d',z'73f7a3', &
&   z'73e9cb',z'73dbf5',z'73ce22',z'73c051',z'73b282',z'73a4b6',z'7396ec',z'738925', &
&   z'737b60',z'736d9d',z'735fdc',z'73521e',z'734462',z'7336a9',z'7328f1',z'731b3c', &
&   z'730d8a',z'72ffd9',z'72f22c',z'72e480',z'72d6d7',z'72c92f',z'72bb8b',z'72ade8', &
&   z'72a048',z'7292aa',z'72850f',z'727775',z'7269de',z'725c4a',z'724eb7',z'724127', &
&   z'723399',z'72260e',z'721884',z'720afd',z'71fd79',z'71eff6',z'71e276',z'71d4f8', &
&   z'71c77c',z'71ba02',z'71ac8b',z'719f16',z'7191a3',z'718433',z'7176c5',z'716959', &
&   z'715bef',z'714e87',z'714122',z'7133bf',z'71265e',z'711900',z'710ba3',z'70fe49', &
&   z'70f0f1',z'70e39b',z'70d648',z'70c8f6',z'70bba7',z'70ae5a',z'70a110',z'7093c7', &
&   z'708681',z'70793d',z'706bfb',z'705ebb',z'70517d',z'704442',z'703709',z'7029d2', &
&   z'701c9d',z'700f6a',z'70023a',z'6ff50c',z'6fe7e0',z'6fdab6',z'6fcd8e',z'6fc068', &
&   z'6fb345',z'6fa624',z'6f9904',z'6f8be7',z'6f7ecd',z'6f71b4',z'6f649d',z'6f5789', &
&   z'6f4a77',z'6f3d67',z'6f3059',z'6f234d',z'6f1643',z'6f093c',z'6efc36',z'6eef33', &
&   z'6ee232',z'6ed533',z'6ec836',z'6ebb3b',z'6eae42',z'6ea14c',z'6e9457',z'6e8765', &
&   z'6e7a74',z'6e6d86',z'6e609a',z'6e53b0',z'6e46c8',z'6e39e3',z'6e2cff',z'6e201d', &
&   z'6e133e',z'6e0661',z'6df985',z'6decac',z'6ddfd5',z'6dd300',z'6dc62d',z'6db95c', &
&   z'6dac8d',z'6d9fc0',z'6d92f5',z'6d862d',z'6d7966',z'6d6ca2',z'6d5fdf',z'6d531f', &
&   z'6d4660',z'6d39a4',z'6d2cea',z'6d2032',z'6d137c',z'6d06c7',z'6cfa15',z'6ced65', &
&   z'6ce0b7',z'6cd40b',z'6cc761',z'6cbab9',z'6cae14',z'6ca170',z'6c94ce',z'6c882e', &
&   z'6c7b90',z'6c6ef5',z'6c625b',z'6c55c3',z'6c492d',z'6c3c9a',z'6c3008',z'6c2378', &
&   z'6c16ea',z'6c0a5f',z'6bfdd5',z'6bf14d',z'6be4c8',z'6bd844',z'6bcbc2',z'6bbf42', &
&   z'6bb2c5',z'6ba649',z'6b99cf',z'6b8d57',z'6b80e2',z'6b746e',z'6b67fc',z'6b5b8c', &
&   z'6b4f1e',z'6b42b2',z'6b3648',z'6b29e0',z'6b1d7a',z'6b1116',z'6b04b4',z'6af854', &
&   z'6aebf5',z'6adf99',z'6ad33f',z'6ac6e6',z'6aba90',z'6aae3b',z'6aa1e9',z'6a9598', &
&   z'6a8949',z'6a7cfd',z'6a70b2',z'6a6469',z'6a5822',z'6a4bdd',z'6a3f9a',z'6a3359', &
&   z'6a271a',z'6a1adc',z'6a0ea1',z'6a0267',z'69f630',z'69e9fa',z'69ddc6',z'69d195', &
&   z'69c565',z'69b937',z'69ad0b',z'69a0e0',z'6994b8',z'698892',z'697c6d',z'69704a', &
&   z'69642a',z'69580b',z'694bee',z'693fd3',z'6933ba',z'6927a2',z'691b8d',z'690f79', &
&   z'690368',z'68f758',z'68eb4a',z'68df3e',z'68d334',z'68c72b',z'68bb25',z'68af20', &
&   z'68a31d',z'68971d',z'688b1d',z'687f20',z'687325',z'68672c',z'685b34',z'684f3e', &
&   z'68434a',z'683758',z'682b68',z'681f7a',z'68138d',z'6807a2',z'67fbb9',z'67efd2', &
&   z'67e3ed',z'67d80a',z'67cc28',z'67c048',z'67b46a',z'67a88e',z'679cb4',z'6790dc', &
&   z'678505',z'677930',z'676d5d',z'67618c',z'6755bd',z'6749ef',z'673e23',z'673259', &
&   z'672691',z'671acb',z'670f06',z'670343',z'66f782',z'66ebc3',z'66e006',z'66d44a', &
&   z'66c891',z'66bcd8',z'66b122',z'66a56e',z'6699bb',z'668e0a',z'66825b',z'6676ae', &
&   z'666b02',z'665f58',z'6653b0',z'66480a',z'663c66',z'6630c3',z'662522',z'661983', &
&   z'660de5',z'66024a',z'65f6b0',z'65eb17',z'65df81',z'65d3ec',z'65c859',z'65bcc8', &
&   z'65b139',z'65a5ab',z'659a1f',z'658e95',z'65830d',z'657786',z'656c01',z'65607e', &
&   z'6554fc',z'65497c',z'653dfe',z'653282',z'652707',z'651b8e',z'651017',z'6504a2', &
&   z'64f92e',z'64edbc',z'64e24c',z'64d6dd',z'64cb70',z'64c005',z'64b49c',z'64a934', &
&   z'649dce',z'64926a',z'648707',z'647ba6',z'647047',z'6464ea',z'64598e',z'644e34', &
&   z'6442db',z'643784',z'642c2f',z'6420dc',z'64158a',z'640a3a',z'63feec',z'63f39f', &
&   z'63e854',z'63dd0b',z'63d1c3',z'63c67d',z'63bb39',z'63aff7',z'63a4b6',z'639976', &
&   z'638e39',z'6382fd',z'6377c3',z'636c8a',z'636153',z'63561e',z'634aea',z'633fb8', &
&   z'633488',z'632959',z'631e2c',z'631301',z'6307d7',z'62fcaf',z'62f189',z'62e664', &
&   z'62db41',z'62d01f',z'62c500',z'62b9e1',z'62aec5',z'62a3aa',z'629890',z'628d79', &
&   z'628263',z'62774e',z'626c3b',z'62612a',z'62561b',z'624b0d',z'624000',z'6234f6', &
&   z'6229ed',z'621ee5',z'6213df',z'6208db',z'61fdd8',z'61f2d7',z'61e7d8',z'61dcda', &
&   z'61d1de',z'61c6e3',z'61bbea',z'61b0f3',z'61a5fd',z'619b09',z'619016',z'618525', &
&   z'617a36',z'616f48',z'61645b',z'615971',z'614e88',z'6143a0',z'6138ba',z'612dd6', &
&   z'6122f3',z'611812',z'610d32',z'610254',z'60f778',z'60ec9d',z'60e1c4',z'60d6ec', &
&   z'60cc16',z'60c141',z'60b66e',z'60ab9c',z'60a0cc',z'6095fe',z'608b31',z'608066', &
&   z'60759c',z'606ad4',z'60600e',z'605549',z'604a85',z'603fc3',z'603503',z'602a44', &
&   z'601f87',z'6014cb',z'600a11',z'5fff58',z'5ff4a1',z'5fe9eb',z'5fdf37',z'5fd485', &
&   z'5fc9d4',z'5fbf24',z'5fb476',z'5fa9ca',z'5f9f1f',z'5f9476',z'5f89ce',z'5f7f28', &
&   z'5f7483',z'5f69df',z'5f5f3e',z'5f549d',z'5f49ff',z'5f3f62',z'5f34c6',z'5f2a2c', &
&   z'5f1f93',z'5f14fc',z'5f0a66',z'5effd2',z'5ef53f',z'5eeaae',z'5ee01f',z'5ed591', &
&   z'5ecb04',z'5ec079',z'5eb5ef',z'5eab67',z'5ea0e0',z'5e965b',z'5e8bd8',z'5e8155', &
&   z'5e76d5',z'5e6c55',z'5e61d8',z'5e575c',z'5e4ce1',z'5e4268',z'5e37f0',z'5e2d79', &
&   z'5e2305',z'5e1891',z'5e0e1f',z'5e03af',z'5df940',z'5deed3',z'5de467',z'5dd9fc', &
&   z'5dcf93',z'5dc52b',z'5dbac5',z'5db061',z'5da5fd',z'5d9b9c',z'5d913b',z'5d86dc', &
&   z'5d7c7f',z'5d7223',z'5d67c9',z'5d5d70',z'5d5318',z'5d48c2',z'5d3e6d',z'5d341a', &
&   z'5d29c8',z'5d1f78',z'5d1529',z'5d0adc',z'5d0090',z'5cf645',z'5cebfc',z'5ce1b4', &
&   z'5cd76e',z'5ccd29',z'5cc2e6',z'5cb8a4',z'5cae63',z'5ca424',z'5c99e6',z'5c8faa', &
&   z'5c856f',z'5c7b36',z'5c70fe',z'5c66c7',z'5c5c92',z'5c525e',z'5c482c',z'5c3dfb', &
&   z'5c33cc',z'5c299d',z'5c1f71',z'5c1546',z'5c0b1c',z'5c00f3',z'5bf6cc',z'5beca7', &
&   z'5be282',z'5bd85f',z'5bce3e',z'5bc41e',z'5bb9ff',z'5bafe2',z'5ba5c6',z'5b9bac', &
&   z'5b9193',z'5b877b',z'5b7d65',z'5b7350',z'5b693d',z'5b5f2a',z'5b551a',z'5b4b0a', &
&   z'5b40fd',z'5b36f0',z'5b2ce5',z'5b22db',z'5b18d3',z'5b0ecc',z'5b04c6',z'5afac2', &
&   z'5af0bf',z'5ae6bd',z'5adcbd',z'5ad2be',z'5ac8c1',z'5abec5',z'5ab4ca',z'5aaad1', &
&   z'5aa0d9',z'5a96e2',z'5a8ced',z'5a82f9',z'5a7906',z'5a6f15',z'5a6525',z'5a5b37', &
&   z'5a514a',z'5a475e',z'5a3d74',z'5a338b',z'5a29a3',z'5a1fbd',z'5a15d8',z'5a0bf4', &
&   z'5a0212',z'59f831',z'59ee51',z'59e473',z'59da96',z'59d0ba',z'59c6e0',z'59bd07', &
&   z'59b330',z'59a959',z'599f84',z'5995b1',z'598bde',z'59820e',z'59783e',z'596e70', &
&   z'5964a3',z'595ad7',z'59510d',z'594744',z'593d7c',z'5933b6',z'5929f1',z'59202d', &
&   z'59166b',z'590caa',z'5902ea',z'58f92b',z'58ef6e',z'58e5b3',z'58dbf8',z'58d23f', &
&   z'58c887',z'58bed0',z'58b51b',z'58ab67',z'58a1b4',z'589803',z'588e53',z'5884a4', &
&   z'587af7',z'58714b',z'5867a0',z'585df6',z'58544e',z'584aa7',z'584101',z'58375d', &
&   z'582dba',z'582418',z'581a77',z'5810d8',z'58073a',z'57fd9d',z'57f402',z'57ea68', &
&   z'57e0cf',z'57d737',z'57cda1',z'57c40c',z'57ba78',z'57b0e6',z'57a754',z'579dc5', &
&   z'579436',z'578aa9',z'57811c',z'577792',z'576e08',z'576480',z'575af9',z'575173', &
&   z'5747ee',z'573e6b',z'5734e9',z'572b68',z'5721e9',z'57186b',z'570eee',z'570572', &
&   z'56fbf8',z'56f27e',z'56e906',z'56df90',z'56d61a',z'56cca6',z'56c333',z'56b9c1', &
&   z'56b051',z'56a6e2',z'569d74',z'569407',z'568a9b',z'568131',z'5677c8',z'566e60', &
&   z'5664fa',z'565b95',z'565231',z'5648ce',z'563f6c',z'56360c',z'562cad',z'56234f', &
&   z'5619f2',z'561097',z'56073c',z'55fde3',z'55f48c',z'55eb35',z'55e1e0',z'55d88c', &
&   z'55cf39',z'55c5e7',z'55bc97',z'55b347',z'55a9f9',z'55a0ad',z'559761',z'558e17', &
&   z'5584cd',z'557b86',z'55723f',z'5568f9',z'555fb5',z'555672',z'554d30',z'5543ef', &
&   z'553ab0',z'553171',z'552834',z'551ef8',z'5515be',z'550c84',z'55034c',z'54fa15', &
&   z'54f0df',z'54e7aa',z'54de77',z'54d544',z'54cc13',z'54c2e3',z'54b9b4',z'54b087', &
&   z'54a75a',z'549e2f',z'549505',z'548bdc',z'5482b5',z'54798e',z'547069',z'546745', &
&   z'545e22',z'545500',z'544be0',z'5442c0',z'5439a2',z'543085',z'542769',z'541e4f', &
&   z'541535',z'540c1d',z'540306',z'53f9f0',z'53f0db',z'53e7c7',z'53deb5',z'53d5a3', &
&   z'53cc93',z'53c384',z'53ba76',z'53b169',z'53a85e',z'539f54',z'53964a',z'538d42', &
&   z'53843b',z'537b36',z'537231',z'53692e',z'53602b',z'53572a',z'534e2a',z'53452b', &
&   z'533c2e',z'533331',z'532a36',z'53213b',z'531842',z'530f4a',z'530654',z'52fd5e', &
&   z'52f469',z'52eb76',z'52e284',z'52d993',z'52d0a3',z'52c7b4',z'52bec6',z'52b5d9', &
&   z'52acee',z'52a404',z'529b1b',z'529233',z'52894c',z'528066',z'527781',z'526e9e', &
&   z'5265bb',z'525cda',z'5253fa',z'524b1b',z'52423d',z'523960',z'523084',z'5227aa', &
&   z'521ed0',z'5215f8',z'520d21',z'52044b',z'51fb76',z'51f2a2',z'51e9cf',z'51e0fe', &
&   z'51d82d',z'51cf5e',z'51c68f',z'51bdc2',z'51b4f6',z'51ac2b',z'51a361',z'519a98', &
&   z'5191d1',z'51890a',z'518045',z'517780',z'516ebd',z'5165fb',z'515d3a',z'51547a', &
&   z'514bbb',z'5142fd',z'513a41',z'513185',z'5128cb',z'512011',z'511759',z'510ea2', &
&   z'5105ec',z'50fd36',z'50f483',z'50ebd0',z'50e31e',z'50da6d',z'50d1be',z'50c90f', &
&   z'50c062',z'50b7b5',z'50af0a',z'50a660',z'509db7',z'50950f',z'508c68',z'5083c2', &
&   z'507b1d',z'507279',z'5069d7',z'506135',z'505894',z'504ff5',z'504757',z'503eb9', &
&   z'50361d',z'502d82',z'5024e8',z'501c4f',z'5013b7',z'500b20',z'50028a',z'4ff9f5', &
&   z'4ff162',z'4fe8cf',z'4fe03d',z'4fd7ad',z'4fcf1d',z'4fc68f',z'4fbe01',z'4fb575', &
&   z'4facea',z'4fa460',z'4f9bd7',z'4f934e',z'4f8ac7',z'4f8241',z'4f79bc',z'4f7139', &
&   z'4f68b6',z'4f6034',z'4f57b3',z'4f4f33',z'4f46b5',z'4f3e37',z'4f35bb',z'4f2d3f', &
&   z'4f24c5',z'4f1c4b',z'4f13d3',z'4f0b5b',z'4f02e5',z'4efa70',z'4ef1fb',z'4ee988', &
&   z'4ee116',z'4ed8a5',z'4ed035',z'4ec7c6',z'4ebf58',z'4eb6ea',z'4eae7e',z'4ea613', &
&   z'4e9daa',z'4e9541',z'4e8cd9',z'4e8472',z'4e7c0c',z'4e73a7',z'4e6b43',z'4e62e1', &
&   z'4e5a7f',z'4e521e',z'4e49be',z'4e4160',z'4e3902',z'4e30a5',z'4e284a',z'4e1fef', &
&   z'4e1796',z'4e0f3d',z'4e06e5',z'4dfe8f',z'4df639',z'4dede5',z'4de591',z'4ddd3f', &
&   z'4dd4ed',z'4dcc9d',z'4dc44d',z'4dbbff',z'4db3b1',z'4dab65',z'4da319',z'4d9acf', &
&   z'4d9285',z'4d8a3d',z'4d81f5',z'4d79af',z'4d7169',z'4d6925',z'4d60e2',z'4d589f', &
&   z'4d505e',z'4d481d',z'4d3fde',z'4d379f',z'4d2f62',z'4d2725',z'4d1eea',z'4d16af', &
&   z'4d0e76',z'4d063d',z'4cfe05',z'4cf5cf',z'4ced99',z'4ce565',z'4cdd31',z'4cd4fe', &
&   z'4ccccd',z'4cc49c',z'4cbc6c',z'4cb43e',z'4cac10',z'4ca3e3',z'4c9bb8',z'4c938d', &
&   z'4c8b63',z'4c833a',z'4c7b12',z'4c72eb',z'4c6ac6',z'4c62a1',z'4c5a7d',z'4c525a', &
&   z'4c4a38',z'4c4217',z'4c39f7',z'4c31d7',z'4c29b9',z'4c219c',z'4c1980',z'4c1165', &
&   z'4c094b',z'4c0131',z'4bf919',z'4bf102',z'4be8eb',z'4be0d6',z'4bd8c1',z'4bd0ae', &
&   z'4bc89b',z'4bc089',z'4bb879',z'4bb069',z'4ba85a',z'4ba04d',z'4b9840',z'4b9034', &
&   z'4b8829',z'4b801f',z'4b7816',z'4b700e',z'4b6807',z'4b6001',z'4b57fc',z'4b4ff7', &
&   z'4b47f4',z'4b3ff2',z'4b37f0',z'4b2ff0',z'4b27f0',z'4b1ff2',z'4b17f4',z'4b0ff7', &
&   z'4b07fc',z'4b0001',z'4af807',z'4af00e',z'4ae816',z'4ae01f',z'4ad829',z'4ad034', &
&   z'4ac83f',z'4ac04c',z'4ab85a',z'4ab068',z'4aa878',z'4aa088',z'4a989a',z'4a90ac', &
&   z'4a88bf',z'4a80d3',z'4a78e8',z'4a70fe',z'4a6915',z'4a612d',z'4a5946',z'4a5160', &
&   z'4a497a',z'4a4196',z'4a39b2',z'4a31d0',z'4a29ee',z'4a220d',z'4a1a2d',z'4a124f', &
&   z'4a0a71',z'4a0294',z'49fab7',z'49f2dc',z'49eb02',z'49e328',z'49db50',z'49d378', &
&   z'49cba2',z'49c3cc',z'49bbf7',z'49b423',z'49ac50',z'49a47e',z'499cad',z'4994dd', &
&   z'498d0d',z'49853f',z'497d71',z'4975a5',z'496dd9',z'49660e',z'495e44',z'49567b', &
&   z'494eb3',z'4946ec',z'493f25',z'493760',z'492f9b',z'4927d8',z'492015',z'491853', &
&   z'491092',z'4908d2',z'490113',z'48f955',z'48f198',z'48e9db',z'48e21f',z'48da65', &
&   z'48d2ab',z'48caf2',z'48c33a',z'48bb83',z'48b3cd',z'48ac18',z'48a463',z'489cb0', &
&   z'4894fd',z'488d4b',z'48859a',z'487dea',z'48763b',z'486e8d',z'4866df',z'485f33', &
&   z'485787',z'484fdd',z'484833',z'48408a',z'4838e2',z'48313b',z'482994',z'4821ef', &
&   z'481a4a',z'4812a6',z'480b04',z'480362',z'47fbc1',z'47f420',z'47ec81',z'47e4e3', &
&   z'47dd45',z'47d5a8',z'47ce0c',z'47c672',z'47bed7',z'47b73e',z'47afa6',z'47a80e', &
&   z'47a078',z'4798e2',z'47914d',z'4789b9',z'478226',z'477a93',z'477302',z'476b71', &
&   z'4763e2',z'475c53',z'4754c5',z'474d37',z'4745ab',z'473e20',z'473695',z'472f0b', &
&   z'472783',z'471ffa',z'471873',z'4710ed',z'470968',z'4701e3',z'46fa5f',z'46f2dc', &
&   z'46eb5a',z'46e3d9',z'46dc59',z'46d4d9',z'46cd5a',z'46c5dd',z'46be60',z'46b6e4', &
&   z'46af68',z'46a7ee',z'46a074',z'4698fb',z'469184',z'468a0c',z'468296',z'467b21', &
&   z'4673ac',z'466c39',z'4664c6',z'465d54',z'4655e3',z'464e72',z'464703',z'463f94', &
&   z'463826',z'4630b9',z'46294d',z'4621e2',z'461a77',z'46130e',z'460ba5',z'46043d', &
&   z'45fcd6',z'45f56f',z'45ee0a',z'45e6a5',z'45df41',z'45d7de',z'45d07c',z'45c91a', &
&   z'45c1ba',z'45ba5a',z'45b2fb',z'45ab9d',z'45a440',z'459ce4',z'459588',z'458e2d', &
&   z'4586d3',z'457f7a',z'457822',z'4570ca',z'456974',z'45621e',z'455ac9',z'455374', &
&   z'454c21',z'4544ce',z'453d7d',z'45362c',z'452edb',z'45278c',z'45203e',z'4518f0', &
&   z'4511a3',z'450a57',z'45030c',z'44fbc1',z'44f477',z'44ed2e',z'44e5e6',z'44de9f', &
&   z'44d759',z'44d013',z'44c8ce',z'44c18a',z'44ba47',z'44b305',z'44abc3',z'44a482', &
&   z'449d42',z'449603',z'448ec5',z'448787',z'44804a',z'44790e',z'4471d3',z'446a99', &
&   z'44635f',z'445c26',z'4454ee',z'444db7',z'444681',z'443f4b',z'443816',z'4430e2', &
&   z'4429af',z'44227c',z'441b4b',z'44141a',z'440cea',z'4405ba',z'43fe8c',z'43f75e', &
&   z'43f031',z'43e905',z'43e1da',z'43daaf',z'43d385',z'43cc5c',z'43c534',z'43be0d', &
&   z'43b6e6',z'43afc0',z'43a89b',z'43a177',z'439a54',z'439331',z'438c0f',z'4384ee', &
&   z'437dcd',z'4376ae',z'436f8f',z'436871',z'436154',z'435a37',z'43531b',z'434c00', &
&   z'4344e6',z'433dcd',z'4336b4',z'432f9c',z'432885',z'43216f',z'431a5a',z'431345', &
&   z'430c31',z'43051e',z'42fe0b',z'42f6f9',z'42efe9',z'42e8d8',z'42e1c9',z'42daba', &
&   z'42d3ad',z'42cca0',z'42c593',z'42be88',z'42b77d',z'42b073',z'42a96a',z'42a261', &
&   z'429b59',z'429452',z'428d4c',z'428647',z'427f42',z'42783e',z'42713b',z'426a39', &
&   z'426337',z'425c36',z'425536',z'424e37',z'424738',z'42403a',z'42393d',z'423241', &
&   z'422b45',z'42244a',z'421d50',z'421657',z'420f5e',z'420866',z'42016f',z'41fa79', &
&   z'41f383',z'41ec8e',z'41e59a',z'41dea7',z'41d7b4',z'41d0c2',z'41c9d1',z'41c2e1', &
&   z'41bbf1',z'41b503',z'41ae14',z'41a727',z'41a03a',z'41994e',z'419263',z'418b79', &
&   z'41848f',z'417da6',z'4176be',z'416fd7',z'4168f0',z'41620a',z'415b25',z'415440', &
&   z'414d5c',z'414679',z'413f97',z'4138b6',z'4131d5',z'412af5',z'412415',z'411d37', &
&   z'411659',z'410f7c',z'41089f',z'4101c3',z'40fae9',z'40f40e',z'40ed35',z'40e65c', &
&   z'40df84',z'40d8ad',z'40d1d6',z'40cb00',z'40c42b',z'40bd57',z'40b683',z'40afb0', &
&   z'40a8de',z'40a20c',z'409b3b',z'40946b',z'408d9c',z'4086cd',z'408000',z'407932', &
&   z'407266',z'406b9a',z'4064cf',z'405e05',z'40573b',z'405072',z'4049aa',z'4042e3', &
&   z'403c1c',z'403556',z'402e91',z'4027cc',z'402109',z'401a45',z'401383',z'400cc1', &
&   z'400600',z'3fff40',z'3ff880',z'3ff1c2',z'3feb03',z'3fe446',z'3fdd89',z'3fd6cd', &
&   z'3fd012',z'3fc957',z'3fc29d',z'3fbbe4',z'3fb52c',z'3fae74',z'3fa7bd',z'3fa107', &
&   z'3f9a51',z'3f939c',z'3f8ce8',z'3f8634',z'3f7f81',z'3f78cf',z'3f721e',z'3f6b6d', &
&   z'3f64bd',z'3f5e0e',z'3f575f',z'3f50b1',z'3f4a04',z'3f4357',z'3f3cac',z'3f3601', &
&   z'3f2f56',z'3f28ac',z'3f2203',z'3f1b5b',z'3f14b3',z'3f0e0c',z'3f0766',z'3f00c1', &
&   z'3efa1c',z'3ef377',z'3eecd4',z'3ee631',z'3edf8f',z'3ed8ee',z'3ed24d',z'3ecbad', &
&   z'3ec50e',z'3ebe6f',z'3eb7d1',z'3eb134',z'3eaa97',z'3ea3fb',z'3e9d60',z'3e96c6', &
&   z'3e902c',z'3e8993',z'3e82fa',z'3e7c62',z'3e75cb',z'3e6f35',z'3e689f',z'3e620a', &
&   z'3e5b76',z'3e54e2',z'3e4e4f',z'3e47bd',z'3e412b',z'3e3a9a',z'3e340a',z'3e2d7a', &
&   z'3e26eb',z'3e205d',z'3e19cf',z'3e1342',z'3e0cb6',z'3e062b',z'3dffa0',z'3df916', &
&   z'3df28c',z'3dec03',z'3de57b',z'3ddef4',z'3dd86d',z'3dd1e7',z'3dcb61',z'3dc4dc', &
&   z'3dbe58',z'3db7d5',z'3db152',z'3daad0',z'3da44f',z'3d9dce',z'3d974e',z'3d90ce', &
&   z'3d8a4f',z'3d83d1',z'3d7d54',z'3d76d7',z'3d705b',z'3d69e0',z'3d6365',z'3d5ceb', &
&   z'3d5671',z'3d4ff9',z'3d4980',z'3d4309',z'3d3c92',z'3d361c',z'3d2fa7',z'3d2932', &
&   z'3d22be',z'3d1c4a',z'3d15d7',z'3d0f65',z'3d08f4',z'3d0283',z'3cfc13',z'3cf5a3', &
&   z'3cef34',z'3ce8c6',z'3ce259',z'3cdbec',z'3cd57f',z'3ccf14',z'3cc8a9',z'3cc23f', &
&   z'3cbbd5',z'3cb56c',z'3caf04',z'3ca89c',z'3ca235',z'3c9bcf',z'3c9569',z'3c8f04', &
&   z'3c889f',z'3c823c',z'3c7bd8',z'3c7576',z'3c6f14',z'3c68b3',z'3c6253',z'3c5bf3', &
&   z'3c5593',z'3c4f35',z'3c48d7',z'3c427a',z'3c3c1d',z'3c35c1',z'3c2f66',z'3c290b', &
&   z'3c22b1',z'3c1c57',z'3c15ff',z'3c0fa7',z'3c094f',z'3c02f8',z'3bfca2',z'3bf64c', &
&   z'3beff7',z'3be9a3',z'3be34f',z'3bdcfc',z'3bd6aa',z'3bd058',z'3bca07',z'3bc3b7', &
&   z'3bbd67',z'3bb718',z'3bb0c9',z'3baa7b',z'3ba42e',z'3b9de1',z'3b9795',z'3b914a', &
&   z'3b8aff',z'3b84b5',z'3b7e6c',z'3b7823',z'3b71db',z'3b6b93',z'3b654c',z'3b5f06', &
&   z'3b58c0',z'3b527b',z'3b4c36',z'3b45f3',z'3b3faf',z'3b396d',z'3b332b',z'3b2cea', &
&   z'3b26a9',z'3b2069',z'3b1a2a',z'3b13eb',z'3b0dad',z'3b076f',z'3b0132',z'3afaf6', &
&   z'3af4ba',z'3aee7f',z'3ae845',z'3ae20b',z'3adbd2',z'3ad599',z'3acf61',z'3ac92a', &
&   z'3ac2f3',z'3abcbd',z'3ab688',z'3ab053',z'3aaa1f',z'3aa3eb',z'3a9db8',z'3a9786', &
&   z'3a9154',z'3a8b23',z'3a84f2',z'3a7ec2',z'3a7893',z'3a7264',z'3a6c36',z'3a6609', &
&   z'3a5fdc',z'3a59b0',z'3a5384',z'3a4d59',z'3a472f',z'3a4105',z'3a3adc',z'3a34b4', &
&   z'3a2e8c',z'3a2864',z'3a223e',z'3a1c18',z'3a15f2',z'3a0fcd',z'3a09a9',z'3a0385', &
&   z'39fd62',z'39f740',z'39f11e',z'39eafd',z'39e4dc',z'39debc',z'39d89d',z'39d27e', &
&   z'39cc60',z'39c642',z'39c025',z'39ba09',z'39b3ed',z'39add2',z'39a7b7',z'39a19d', &
&   z'399b84',z'39956b',z'398f53',z'39893b',z'398324',z'397d0e',z'3976f8',z'3970e3', &
&   z'396ace',z'3964ba',z'395ea7',z'395894',z'395282',z'394c70',z'39465f',z'39404f', &
&   z'393a3f',z'393430',z'392e21',z'392813',z'392206',z'391bf9',z'3915ed',z'390fe1', &
&   z'3909d6',z'3903cb',z'38fdc1',z'38f7b8',z'38f1af',z'38eba7',z'38e5a0',z'38df99', &
&   z'38d993',z'38d38d',z'38cd88',z'38c783',z'38c17f',z'38bb7c',z'38b579',z'38af77', &
&   z'38a975',z'38a374',z'389d73',z'389774',z'389174',z'388b76',z'388577',z'387f7a', &
&   z'38797d',z'387381',z'386d85',z'38678a',z'38618f',z'385b95',z'38559b',z'384fa2', &
&   z'3849aa',z'3843b2',z'383dbb',z'3837c5',z'3831cf',z'382bd9',z'3825e4',z'381ff0', &
&   z'3819fd',z'381409',z'380e17',z'380825',z'380234',z'37fc43',z'37f653',z'37f063', &
&   z'37ea74',z'37e485',z'37de97',z'37d8aa',z'37d2bd',z'37ccd1',z'37c6e5',z'37c0fa', &
&   z'37bb10',z'37b526',z'37af3d',z'37a954',z'37a36c',z'379d84',z'37979d',z'3791b6', &
&   z'378bd0',z'3785eb',z'378006',z'377a22',z'37743e',z'376e5b',z'376879',z'376297', &
&   z'375cb5',z'3756d5',z'3750f4',z'374b15',z'374535',z'373f57',z'373979',z'37339b', &
&   z'372dbf',z'3727e2',z'372206',z'371c2b',z'371651',z'371077',z'370a9d',z'3704c4', &
&   z'36feec',z'36f914',z'36f33d',z'36ed66',z'36e790',z'36e1ba',z'36dbe5',z'36d611', &
&   z'36d03d',z'36ca69',z'36c497',z'36bec4',z'36b8f3',z'36b321',z'36ad51',z'36a781', &
&   z'36a1b1',z'369be2',z'369614',z'369046',z'368a79',z'3684ac',z'367ee0',z'367915', &
&   z'36734a',z'366d7f',z'3667b5',z'3661ec',z'365c23',z'36565b',z'365093',z'364acc', &
&   z'364505',z'363f3f',z'363979',z'3633b4',z'362df0',z'36282c',z'362269',z'361ca6', &
&   z'3616e4',z'361122',z'360b61',z'3605a0',z'35ffe0',z'35fa20',z'35f461',z'35eea3', &
&   z'35e8e5',z'35e328',z'35dd6b',z'35d7af',z'35d1f3',z'35cc38',z'35c67d',z'35c0c3', &
&   z'35bb09',z'35b550',z'35af98',z'35a9e0',z'35a429',z'359e72',z'3598bb',z'359306', &
&   z'358d50',z'35879c',z'3581e8',z'357c34',z'357681',z'3570ce',z'356b1c',z'35656b', &
&   z'355fba',z'355a09',z'355459',z'354eaa',z'3548fb',z'35434d',z'353d9f',z'3537f2', &
&   z'353245',z'352c99',z'3526ee',z'352143',z'351b98',z'3515ee',z'351045',z'350a9c'/)

end module invsqrt_table
#endif


#ifdef MOLFILE
module MOLFILEparam
integer, dimension(:), allocatable :: Handle
end module MOLFILEparam
#endif

