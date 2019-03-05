! ############################
! ## SUBROUTINE LIST 
! ## -- CG_distributions 
! ## -- RadialW 
! ## -- RadialC 
! ## -- Write_CGSites 
! ## -- SingleLipidAnim 
! ## -- Write_SingleLipid 
! ############################


!######################################################################
!######################################################################


module CGdist
   integer :: NumTypes
   integer :: NumGroup
   real(8) :: DRgrid, InvGrid
   integer :: Kw, Kcomp, NMcomp, NAcomp
   integer, dimension(:), allocatable :: GroupNumber, SecondG, NumAtm_Group
   real(8), dimension(:), allocatable :: Mass_Group
   real(8), dimension(:,:,:), allocatable :: Rcom_Group
   real(8), dimension(:,:), allocatable :: Rg_Group
   integer, dimension(:), allocatable :: NtypeGroup
   integer, dimension(:), allocatable :: NinTypeGroup
   integer, dimension(:,:), allocatable :: NonInteractPair
   integer, dimension(:,:), allocatable :: NGpairIJ
   integer, dimension(:,:), allocatable :: igr
end module CGdist


!######################################################################
!######################################################################


subroutine CG_distributions

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonBlocks, only : QPBC
use ParamAnalyze, only : NJobs, NTrjStep, Interval
use CGdist
use TimeParam, only : Timeps
use UnitExParam, only : kb, cvol, pi
use BathParam, only : Temp_o
use CellParam, only : H, Volume, InvH
use CutoffParam, only : Rcutoff2
use AtomParam, only : ResidName, Mass

implicit none

character(len=72) :: FileName
integer :: i, j, k, l, ii, Numa, NumF, jj, kk, IRcut, ll
integer :: NumCGBond, NumCGAngle, NumGR, nn, totalframe
real(8), dimension(3) :: Rij
integer, dimension(:), allocatable :: BondCGI, BondCGJ
integer, dimension(:), allocatable :: AngleCGI, AngleCGJ, AngleCGK
integer, dimension(:,:), allocatable :: CGLength
integer, dimension(:,:), allocatable :: CGAngle
real(8), dimension(:), allocatable :: AveBondL, Ave2BondL
real(8), dimension(:), allocatable :: AveAngle, Ave2Angle
real(8), dimension(800) :: DisBond
real(8), dimension(180) :: DisAngle
real(8) :: Vol, det, qk, gg, drr, xx, R2, R1, cst, RR1, RR2, deg, gn
real(8) :: cr, pp, ppuw, StdBond, StdAngle
integer :: NGpair, IR, Ics
real(8), dimension(3) :: Rkj
integer, dimension(:), allocatable :: Count
character(len=4) :: Aname
integer, dimension(20) :: icntrl
integer :: nstr, nnall
character(len=72) :: String1, String
external det

   Kw = 0
   ii = 1
   do i = 1 , NumSpec
     if(ResidName(ii)=='TIP3') Kw = i
     ii = ii + NumMol(i)*NumAtm(i)
   end do

   print *, 'Kw=', Kw

   open(21,file='CGmodel.data',status='unknown')

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     select case(nn)

     case(1)
       read(String,*) Kcomp

       NAcomp = NumAtm(Kcomp)
       NMcomp = NumMol(Kcomp)

! ## to get 'Numa'

       ii = 1

       do i = 1 , NumSpec

         if(i == Kcomp) then
           Numa = ii - 1
         end if

         ii = ii + NumMol(i) * NumAtm(i)

       end do

       allocate( GroupNumber(NAcomp) ) ! Group Number ( atom number ) within a molecule 
       allocate( SecondG(NAcomp) ) ! Group Number ( atom number ) within a molecule 

! ## Number of Group in a molecule, Number of Group Types in a molecule 

     case(2)

       read(String,*) NumGroup, NumTypes

       allocate(NumAtm_Group(NumGroup))
       allocate(Mass_Group(NumGroup))
       allocate(Rcom_Group(3,NumGroup,NMcomp))
       allocate(Rg_Group(NumGroup,NMcomp))

       allocate( NTypeGroup(NumGroup) )   ! ## Type number of a group in a molecule
       allocate( NinTypeGroup(NumTypes) ) ! ## Number of groups of a type

       exit

     end select

   end do

! ## read group number of an atom in a molecule

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     ll = 0
     do l = 2, 72
       if((String(l:l)==' ').and.(String(l-1:l-1)/=' ')) then
         ll = ll + 1
       end if
     end do
     if(ll==2) then
       read(String,*) j, GroupNumber(j)
       SecondG(j) = 0
     else if(ll==3) then
       read(String,*) j, GroupNumber(j), SecondG(j)
     else
       write(*,*) "ERROR: GroupNumbers"
       stop
     end if

     if(nn == NAcomp) then
       if(j /= NAcomp) then
         write(*,*) 'ERROR : reading GroupNumber'
         call Finalize
       end if
       exit
     end if

   end do

   write(*,*) 'atom vs. groupnumver'
   do i = 1, nn
     print *, i, GroupNumber(i), SecondG(i)
   end do

! ## read type number of a group

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     read(String,*) j, NTypeGroup(j)

     if(nn == NumGroup) then
       if(j /= NumGroup) then
         write(*,*) 'ERROR : reading NTypeGroup'
         call Finalize
       end if
       exit
     end if

   end do

   write(*,*) 'group vs. type'
   do i = 1, nn
     print *, i, NTypeGroup(i)
   end do

! ## Number of Bond in a molecule

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     read(String,*) NumCGBond

     exit

   end do

   allocate( BondCGI(NumCGBond) )
   allocate( BondCGJ(NumCGBond) )
   allocate( AveBondL(NumCGBond) )
   allocate( Ave2BondL(NumCGBond) )
   AveBondL = 0.d0
   Ave2BondL = 0.d0

! ## Bond pair

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     read(String,*) BondCGI(nn), BondCGJ(nn)

     if(nn == NumCGBond) exit

   end do

   print *, 'bond pairs'
   do i = 1, NumCGBond
     print *, BondCGI(i), BondCGJ(i)
   end do

! ## Number of bending pairs in a molecule

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     read(String,*) NumCGAngle

     exit

   end do

   if(NumCGAngle/=0) then

   allocate( AngleCGI(NumCGAngle) )
   allocate( AngleCGJ(NumCGAngle) )
   allocate( AngleCGK(NumCGAngle) )
   allocate( AveAngle(NumCGAngle) )
   allocate( Ave2Angle(NumCGAngle) )

   AveAngle = 0.d0
   Ave2Angle = 0.d0

! ## Angle pair

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     read(String,*) AngleCGI(nn), AngleCGJ(nn), AngleCGK(nn)

     if(nn == NumCGAngle) exit

   end do

   print *, 'bond pairs'
   do i = 1, NumCGAngle
     print *, AngleCGI(i), AngleCGJ(i), AngleCGK(i)
   end do

   end if

   close(21)

   allocate(CGLength(NumCGBond,800))
   allocate(CGAngle(NumCGAngle,180))

   CGLength = 0
   CGAngle = 0

   DRgrid = 0.01d0
   InvGrid = 1.d0 / DRgrid

   IRcut = int(sqrt(Rcutoff2)*InvGrid) + 1

   NinTypeGroup(:) = 0
   do i = 1, NumGroup
     j = NTypeGroup(i)
     NinTypeGroup(j) = NinTypeGroup(j) + 1
   end do

   NGpair = (NumTypes*(NumTypes-1))/2 + NumTypes

   allocate( NGpairIJ(NumTypes, NumTypes) )

   k = NumTypes
   do i = 1, NumTypes
     do j = i, NumTypes
       k = k + 1
       NGpairIJ(i,j) = k
     end do
   end do

   do j = 1 , NumTypes-1
     do i = j + 1, NumTypes
       NGpairIJ(i,j) = NGpairIJ(j,i)
     end do
  end do

   NumGR = NumTypes+NGPair ! NumTypes is for gr with Water
   allocate( igr(IRcut,NumGR) )
   igr = 0

! ## Number of atoms in a group, Mass of a group
   NumAtm_Group = 0
   Mass_Group = 0.d0

   do i = 1, NAcomp

     j = Numa + i
     k = GroupNumber(i)
     l = SecondG(i)

     NumAtm_Group(k) = NumAtm_Group(k) + 1
     if(SecondG(i) == 0) then
       Mass_Group(k)   = Mass_Group(k)   + Mass(j)
     else
       NumAtm_Group(l) = NumAtm_Group(l) + 1
       Mass_Group(k)   = Mass_Group(k)   + Mass(j) * 0.5d0
       Mass_Group(l)   = Mass_Group(l)   + Mass(j) * 0.5d0
     end if

   end do

   allocate( NonInteractPair(7,NumGroup) )

   NonInteractPair = 0

   allocate( Count(NumGroup) ) ! ## for temporary use

   Count = 0

   do i = 1, NumCGBond

     ii = BondCGI(i)
     jj = BondCGJ(i)

     Count(ii) = Count(ii) + 1
     Count(jj) = Count(jj) + 1

     NonInteractPair(Count(ii),ii) = jj
     NonInteractPair(Count(jj),jj) = ii

   end do

   do i = 1 , NumCGAngle

     ii = AngleCGI(i)
     jj = AngleCGK(i)

     Count(ii) = Count(ii) + 1
     Count(jj) = Count(jj) + 1

     NonInteractPair(Count(ii),ii) = jj
     NonInteractPair(Count(jj),jj) = ii

   end do

   open(41,file='./Analy/CGconfig.dcd',status='unknown',form='unformatted')

   totalframe = 0
   do i = 1, NJobs
     totalframe = totalframe + NTrjStep(i)/Interval(i)
   end do

   Aname = 'CORD'
   icntrl = 0
   nstr = 0
   nnall = NumGroup*NMcomp
   icntrl(1) = totalframe  ! number of frames
   icntrl(2) = 1           ! number of steps in previous run
   icntrl(3) = 1           ! frequency of saving
   icntrl(4) = totalframe  ! total number of steps
   icntrl(8) = NumGroup*NMcomp*3       ! number of degrees of freedm
   icntrl(10) = 981668463   ! coded time step
   icntrl(11) = 1           ! coded crystallographic group (or zro)
   icntrl(20) = 27        ! CHARMM version number

   write(41) Aname,icntrl
   write(41) nstr
   write(41) nnall

   NumF = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ------------------

       if( mod(j,Interval(i)) == 0 ) then

         NumF = NumF + 1

#ifdef MOLFILE
         Timeps = dble(NumF)
         call InversMatrix(H,InvH)
#else
         if(QPBC) call CellTransform
#endif

         call CalcCOM_Gyr(Numa)

         do l = 1 , NumCGBond

           ii = BondCGI(l)
           jj = BondCGJ(l)

           do k = 1 , NMcomp

             Rij = Rcom_Group(:,ii,k) - Rcom_Group(:,jj,k)
             R2  = dot_product( Rij, Rij )
             R1  = sqrt( R2 )
             IR  = int( R1 * InvGrid ) + 1
             if(IR<1) then
               write(*,*) 'bond length error,',IR
               IR=1
             else if(IR>800) then
               write(*,*) 'bond length error,',IR
               IR=800
             end if
             CGLength( l, IR ) = CGLength( l, IR ) + 1
             AveBondL(l) = AveBondL(l) + R1
             Ave2BondL(l) = Ave2BondL(l) + R2

           end do

         end do


         do l = 1 , NumCGAngle

           ii = AngleCGI(l)
           jj = AngleCGJ(l)
           kk = AngleCGK(l)

           do k = 1 , NMcomp

             Rij = Rcom_Group(:,ii,k) - Rcom_Group(:,jj,k)
             Rkj = Rcom_Group(:,kk,k) - Rcom_Group(:,jj,k)
             cst = dot_product( Rij, Rkj )
             RR1 = dot_product( Rij, Rij )
             RR2 = dot_product( Rkj, Rkj )
             cst = cst / sqrt( RR1*RR2 )

             deg = acos( cst ) * 180. / pi
             Ics = int( deg ) + 1
             if(Ics<1) then
               write(*,*) 'ANGLE error,',Ics
               Ics = 1
             else if(Ics>180) then
               write(*,*) 'ANGLE error,',Ics
               Ics = 180
             end if
             CGAngle(l,Ics) = CGAngle(l,Ics) + 1
             AveAngle(l) = AveAngle(l) + deg
             Ave2Angle(l) = Ave2Angle(l) + deg*deg

           end do

         end do

         if(QPBC) then

         if(Kw/=0) then
           call RadialW
         end if
         call RadialC

         Volume = det(H)
         Vol = Vol + Volume

         end if

         call Write_CGSites

       end if

     end do

   end do

   close(41)

   open(111,file='./Analy/IntraSummary.dat')

   do i = 1 , NumCGBond

     ii = BondCGI(i)
     jj = BondCGJ(i)

     if((ii<10).and.(jj<10)) then
       write(FileName,'(a,i1,a,i1,a)') './Analy/Bond0',ii,'--0',jj,'.dat'
     else if(jj<10) then
       write(FileName,'(a,i2,a,i1,a)') './Analy/Bond',ii,'--0',jj,'.dat'
     else if(ii<10) then
       write(FileName,'(a,i1,a,i2,a)') './Analy/Bond0',ii,'--',jj,'.dat'
     else
       write(FileName,'(a,i2,a,i2,a)') './Analy/Bond',ii,'--',jj,'.dat'
     end if

     AveBondL(i)  = AveBondL(i)  / dble(NMcomp*NumF)
     Ave2BondL(i) = Ave2BondL(i) / dble(NMcomp*NumF)
     StdBond = sqrt( Ave2BondL(i) - AveBondL(i)*AveBondL(i) )

     write(111,'(/a,i2,a,i2)') 'For Bond',ii,'--',jj
     write(111,'(a,f12.7)') '   Average  =',AveBondL(i)
     write(111,'(a,f12.7/)') '   St. dev  =',StdBond

     open(42,file=trim(FileName),status='unknown')

     write(42,'(a)') '# length / weighted dist. / potential / unweighted dist.'

     do j = 1, 800
       xx = (j-0.5)*DRgrid
       DisBond(j) = dble(CGLength(i,j)) / dble(NMcomp*NumF) / (xx*xx) ! probability
     end do
     xx = 0.d0
     do j = 1, 800
       xx = xx + DisBond(j)
     end do
     do j = 1, 800
       DisBond(j) = DisBond(j) / xx
     end do
     do j = 1 , 800
       xx = dble(CGLength(i,j)) / dble(NMcomp*NumF)
       if(DisBond(j)==0.) then
         pp = 0.
         ppuw = 0.
       else
         pp = - cvol * kb * Temp_o * log(DisBond(j))
         ppuw = - cvol * kb * Temp_o * log(xx)
       end if
       write(42,'(f7.3,2(f15.10,e15.8))') real(j-0.5)*DRgrid,DisBond(j)/DRgrid,pp,xx/DRgrid,ppuw
     end do

     close(42)

   end do

   do i = 1 , NumCGAngle

     ii = AngleCGI(i)
     jj = AngleCGJ(i)
     kk = AngleCGK(i)

     if((ii<10).and.(jj<10).and.(kk<10)) then
       write(FileName,'(a,i1,a,i1,a,i1,a)') './Analy/Angle0',ii,'--0',jj,'--0',kk,'.dat'
     else if((jj<10).and.(kk<10)) then
       write(FileName,'(a,i2,a,i1,a,i1,a)') './Analy/Angle',ii,'--0',jj,'--0',kk,'.dat'
     else if((jj<10).and.(ii<10)) then
       write(FileName,'(a,i1,a,i1,a,i2,a)') './Analy/Angle0',ii,'--0',jj,'--',kk,'.dat'
     else if((ii<10).and.(kk<10)) then
       write(FileName,'(a,i1,a,i2,a,i1,a)') './Analy/Angle0',ii,'--',jj,'--0',kk,'.dat'
     else if(ii<10) then
       write(FileName,'(a,i1,a,i2,a,i2,a)') './Analy/Angle0',ii,'--',jj,'--',kk,'.dat'
     else if(jj<10) then
       write(FileName,'(a,i2,a,i1,a,i2,a)') './Analy/Angle',ii,'--0',jj,'--',kk,'.dat'
     else if(kk<10) then
       write(FileName,'(a,i2,a,i2,a,i1,a)') './Analy/Angle',ii,'--',jj,'--0',kk,'.dat'
     else
       write(FileName,'(a,i2,a,i2,a,i2,a)') './Analy/Angle',ii,'--',jj,'--',kk,'.dat'
     end if

     AveAngle(i)  = AveAngle(i)  / dble(NMcomp*NumF)
     Ave2Angle(i) = Ave2Angle(i) / dble(NMcomp*NumF)
     StdAngle = sqrt( Ave2Angle(i) - AveAngle(i)*AveAngle(i) )

     write(111,'(/a,i2,2(a,i2))') 'For Angle',ii,'--',jj,'--',kk
     write(111,'(a,f12.7)') '   Average  =',AveAngle(i)
     write(111,'(a,f12.7/)') '   St. dev  =',StdAngle

     open(43,file=trim(FileName),status='unknown')

     write(43,'(a)') '# length / weighted dist. / potential / unweighted dist.'

     do j = 1 , 180
       xx = (j - 0.5d0)/180.d0*pi
       DisAngle(j) = dble(CGAngle(i,j)) / dble(NMcomp*NumF) / sin(xx)
     end do
     xx = 0.d0
     do j = 1, 180
       xx = xx + DisAngle(j)
     end do
     do j = 1, 180
       DisAngle(j) = DisAngle(j) / xx
     end do
     do j = 1 , 180
       xx = dble(CGAngle(i,j)) / dble(NMcomp*NumF)
       if(DisAngle(j)==0.) then
         pp = 0.d0
       else
         pp = - cvol * kb * Temp_o * log(DisAngle(j))
       end if
       write(43,'(f6.1,f15.10,e15.8,f15.10)') real(j)-0.5,DisAngle(j),pp,xx
     end do

     close(43)

   end do

   if(QPBC) then

   Vol = Vol / dble(NumF)

   if(Kw/=0) then

   do i = 1, NumTypes

     if(i<10) then
       write(FileName,'(a,i1,a)') './Analy/GR_Type0',i,'--W.dat'
     else
       write(FileName,'(a,i2,a)') './Analy/GR_Type',i,'--W.dat'
     end if

     open(44,file = trim(FileName), status='unknown')

     cr = 0.d0

     do j = 1 , IRcut

       drr = (j-0.5d0) * DRgrid
       qk = Vol / (NumMol(Kw) * 4.d0 * pi * drr * drr * DRgrid)
       gn = dble( igr(j,i) ) / dble(NumF*NinTypeGroup(i)*NMcomp)
       gg = gn * qk
       cr = cr + gn

       write(44,'(f8.3,2f12.6,e15.7)') drr, gg, gn, cr

     end do

     close(44)

   end do

   end if

   do i = 1 , NumTypes

     do j = i, NumTypes

       if((i<10).and.(j<10)) then
         write(FileName,'(a,i1,a,i1,a)') './Analy/GR_Type0',i,'--Type0',j,'.dat'
       else if(i<10) then
         write(FileName,'(a,i1,a,i2,a)') './Analy/GR_Type0',i,'--Type',j,'.dat'
       else if(j<10) then
         write(FileName,'(a,i2,a,i1,a)') './Analy/GR_Type',i,'--Type0',j,'.dat'
       else
         write(FileName,'(a,i2,a,i2,a)') './Analy/GR_Type',i,'--Type',j,'.dat'
       end if

       open(45,file = trim(FileName), status='unknown')

       k = NGpairIJ(i,j)

       cr = 0.d0
       do l = 1 , IRcut

         drr = (l-0.5d0) * DRgrid
         qk = Vol / (NinTypeGroup(j) * NMcomp * 4.d0 * pi * drr * drr * DRgrid)
         gn = dble( igr(l,k) ) / dble(NumF * NinTypeGroup(i) * NMcomp)
         if(i==j) gn = gn * 2.d0
         gg = gn * qk
         cr = cr + gn

         write(45,'(f8.3,2f12.6,e15.7)') drr, gg, gn, cr

       end do

       close(45)

     end do

   end do

   end if

end subroutine CG_distributions


!######################################################################
!######################################################################


subroutine RadialW

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use CGdist
use CutoffParam, only : Rcutoff2
use CellParam, only : H, InvH

implicit none

integer :: i , j, k, ii
integer :: ir, Numb
real(8), dimension(3) :: Rij, Sij
real(8), dimension(3,NumMol(Kw)) :: ScRw
real(8), dimension(3,NumGroup,NMcomp) :: ScRcom
real(8) :: R2, R1

   Numb = 0
   if( Kw/=1 ) then
     do i = 1, Kw-1
       Numb = Numb + NumMol(i)*NumAtm(i)
     end do
   end if

   do i = 1 , NumMol(Kw)

     do j = 1 , NumAtm(Kw)

       Numb = Numb + 1
       if(j==1) then
         ScRw(:,i) = matmul( InvH , R(:,Numb) )
       end if

     end do

   end do

   do i = 1, NumGroup

     do j = 1, NMcomp

       ScRcom(:,i,j) = matmul( InvH, Rcom_Group(:,i,j) )

     end do

   end do

   do i = 1, NumGroup

     ii = NTypeGroup(i)

     do j = 1 , NMcomp

       do k = 1 , NumMol(Kw)

         Sij = ScRcom(:,i,j) - ScRw(:,k)
         Sij = Sij - nint(Sij)
         Rij = matmul( H, Sij )
         R2 = dot_product( Rij, Rij )

         if(R2 < Rcutoff2) then

           R1 = sqrt(R2)
           ir = int(R1*InvGrid) + 1
           igr(ir,ii) = igr(ir,ii) + 1

         end if

       end do

     end do

   end do


end subroutine RadialW


!######################################################################
!######################################################################


subroutine RadialC

use Numbers, only : NumMol
use CGdist
use CutoffParam, only : Rcutoff2
use CellParam, only : H, InvH

implicit none

integer :: i , j, k, ii, jj, kk, l
integer :: ir
real(8), dimension(3) :: Rij, Sij
real(8), dimension(3,NumGroup,NMcomp) :: ScRcom
real(8) :: R2, R1

   do i = 1, NumGroup

     do j = 1, NMcomp

       ScRcom(:,i,j) = matmul( InvH, Rcom_Group(:,i,j) )

     end do

   end do

   do i = 1, NumGroup-1

     ii = NtypeGroup(i)

     do j = i+1 , NumGroup

       jj = NtypeGroup(j)

       if((j==NonInteractPair(1,i)).or. &
       &  (j==NonInteractPair(2,i)).or. &
       &  (j==NonInteractPair(3,i)).or. &
       &  (j==NonInteractPair(4,i)).or. &
       &  (j==NonInteractPair(5,i)).or. &
       &  (j==NonInteractPair(6,i)).or. &
       &  (j==NonInteractPair(7,i))) cycle

       kk = NGpairIJ(ii,jj)

       do k = 1 , NMcomp

         Sij = ScRcom(:,i,k) - ScRcom(:,j,k)
         Sij = Sij - nint(Sij)
         Rij = matmul( H, Sij )
         R2 = dot_product( Rij, Rij )

         if(R2 < Rcutoff2) then

           R1 = sqrt(R2)
           ir = int(R1*InvGrid) + 1
           igr(ir,kk) = igr(ir,kk) + 1

         end if

       end do

     end do

   end do

   do i = 1, NumGroup

     ii = NtypeGroup(i)

     do j = 1 , NumGroup

       jj = NtypeGroup(j)

       kk = NGpairIJ(ii,jj)

       do k = 1 , NMcomp-1

         do l = k+1 , NMcomp

           Sij = ScRcom(:,i,k) - ScRcom(:,j,l)
           Sij = Sij - nint(Sij)
           Rij = matmul( H, Sij )
           R2 = dot_product( Rij, Rij )

           if(R2 < Rcutoff2) then

             R1 = sqrt(R2)
             ir = int(R1*InvGrid) + 1
             igr(ir,kk) = igr(ir,kk) + 1

           end if

         end do

       end do

     end do

   end do


end subroutine RadialC


!######################################################################
!######################################################################


subroutine Write_CGSites

use CommonBlocks, only : QPBC
use CGdist
use CellParam, only : H

implicit none

integer :: ii, k, l, i, nnall
real(8), dimension(6) :: xcell
real, dimension(NMcomp*NumGroup) :: X, Y, Z

   xcell(1) = H(1,1)
   xcell(2) = 0.d0
   xcell(3) = H(2,2)
   xcell(4) = 0.d0
   xcell(5) = 0.d0
   xcell(6) = H(3,3)

   nnall = NMcomp*NumGroup

   ii = 0
   do k = 1 , NMcomp
     do l = 1 , NumGroup
       ii = ii + 1
       X(ii) = sngl(Rcom_Group(1,l,k))
       Y(ii) = sngl(Rcom_Group(2,l,k))
       Z(ii) = sngl(Rcom_Group(3,l,k))
     end do
   end do

   write(41) xcell
   write(41) (X(i),i=1,nnall)
   write(41) (Y(i),i=1,nnall)
   write(41) (Z(i),i=1,nnall)

end subroutine Write_CGSites


!######################################################################
!######################################################################


subroutine CalcCOM_Gyr(Numa)

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use CGdist
use AtomParam, only : Mass

implicit none

integer :: k, ii, ll, jj, l, Numa, nn
real(8), dimension(3) :: Rij

   do k = 1 , NMcomp

     ii = Numa + (k-1)*NAcomp

! ##  Center of Mass

     Rcom_Group(:,:,k) = 0.d0

     do l = 1 , NAcomp

       ll = GroupNumber(l)
       jj = ii + l

       nn = SecondG(l)

       if(nn==0) then
         Rcom_Group(:,ll,k) = Rcom_Group(:,ll,k) + R(:,jj) * Mass(jj)
       else
         Rcom_Group(:,ll,k) = Rcom_Group(:,ll,k) + R(:,jj) * Mass(jj) * 0.5
         Rcom_Group(:,nn,k) = Rcom_Group(:,nn,k) + R(:,jj) * Mass(jj) * 0.5
       end if

     end do

     do l = 1 , NumGroup

       Rcom_Group(:,l,k) = Rcom_Group(:,l,k) / Mass_Group(l)

     end do

! ## Radius of gyration

     Rg_Group(:,k) = 0.d0

     do l = 1 , NAcomp

       ll = GroupNumber(l)
       jj = ii + l

       nn = SecondG(l)

       if(nn==0) then
         Rij = R(:,jj) - Rcom_Group(:,ll,k)
         Rg_Group(ll,k) = Rg_Group(ll,k) + Mass(jj) * dot_product( Rij, Rij )
       else
         Rij = R(:,jj) - Rcom_Group(:,ll,k)
         Rg_Group(ll,k) = Rg_Group(ll,k) + 0.5d0 * Mass(jj) * dot_product( Rij, Rij )
         Rij = R(:,jj) - Rcom_Group(:,nn,k)
         Rg_Group(nn,k) = Rg_Group(nn,k) + 0.5d0 * Mass(jj) * dot_product( Rij, Rij )
       end if

     end do

     do l = 1 , NumGroup

       Rg_Group(l,k) = Rg_Group(l,k) / Mass_Group(l)
       Rg_Group(l,k) = sqrt( Rg_Group(l,k) )

     end do

   end do

end subroutine CalcCOM_Gyr


!######################################################################
!######################################################################


subroutine SingleLipidAnim

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze, only : NJobs, NTrjStep, Interval, Nlipid
use CGdist
use AtomParam, only : ResidName, Mass

implicit none

integer :: i, j, k, l, ii, Numa, NumF, jj, ll
real(8) :: Mass_Lipid
real(8), dimension(3) :: Rcom, Rij

   ii = 1
   do i = 1 , NumSpec
     if((ResidName(ii)=='DPPC').or.(ResidName(ii)=='PhPC')) then
       Kcomp = i
       Numa  = ii - 1
     end if
     ii = ii + NumMol(i)*NumAtm(i)
   end do

   NAcomp = NumAtm(Kcomp)
   NMcomp = NumMol(Kcomp)

   allocate( GroupNumber(NAcomp) )

   call DefineGroup(Numa)

   NumF = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ------------------

       if( mod(j,Interval(i)) == 0 ) then

         NumF = NumF + 1

         call CellTransform

! ##     do k = 1 , NMcomp
         k = Nlipid
         ii = Numa + (k-1)*NAcomp

         Rcom = 0.d0
         Mass_Lipid = sum(Mass_Group)
         do l = 1 , NAcomp
           jj = ii + l
           Rcom = Rcom + R(:,jj) * Mass(jj)
         end do
         Rcom = Rcom / Mass_Lipid

         do l = 1 , NAcomp
           jj = ii + l
           R(:,jj) = R(:,jj) - Rcom
         end do

         Rcom_Group(:,:,k) = 0.
         do l = 1 , NAcomp
           ll = GroupNumber(l)
           jj = ii + l
           Rcom_Group(:,ll,k) = Rcom_Group(:,ll,k) + R(:,jj) * Mass(jj)
         end do
         do l = 1 , NumGroup
           Rcom_Group(:,l,k) = Rcom_Group(:,l,k) / Mass_Group(l)
         end do

         Rg_Group(:,k) = 0.d0
         do l = 1 , NAcomp
           ll = GroupNumber(l)
           jj = ii + l
           Rij = R(:,jj) - Rcom_Group(:,ll,k)
           Rg_Group(ll,k) = Rg_Group(ll,k) + Mass(jj) * dot_product( Rij, Rij )
         end do
         do l = 1 , NumGroup
           Rg_Group(l,k) = Rg_Group(l,k) / Mass_Group(l)
           Rg_Group(l,k) = sqrt( Rg_Group(l,k) )
         end do

         call Write_SingleLipid(NumF,Numa)

       end if

     end do

   end do

end subroutine SingleLipidAnim


!######################################################################
!######################################################################


subroutine Write_SingleLipid(NumF,Numa)

use Configuration, only : R
use ParamAnalyze, only : Nlipid
use CGdist
use AtomParam, only : AtomName, ResidName
use Timeparam, only : Timeps
use Numbers, only : NumAtm

implicit none

integer :: i, j, ii, NumF, Numa
character(len=1) :: Ch
character(len=4) :: AName
character(len=72) :: FilePDB, FileRG
real, parameter :: zero = 0.
integer, parameter :: one = 1

   if(NumF<10) then
     write(FilePDB,'(a,i1,a)') './Analy/snap0000',NumF,'.pdb'
     write(FileRG ,'(a,i1  )') './Analy/RofG0000',NumF
   else if(NumF<100) then
     write(FilePDB,'(a,i2,a)') './Analy/snap000',NumF,'.pdb'
     write(FileRG ,'(a,i2  )') './Analy/RofG000',NumF
   else if(NumF<1000) then
     write(FilePDB,'(a,i3,a)') './Analy/snap00',NumF,'.pdb'
     write(FileRG ,'(a,i3  )') './Analy/RofG00',NumF
   else if(NumF<10000) then
     write(FilePDB,'(a,i4,a)') './Analy/snap0',NumF,'.pdb'
     write(FileRG ,'(a,i4  )') './Analy/RofG0',NumF
   else if(NumF<100000) then
     write(FilePDB,'(a,i5,a)') './Analy/snap',NumF,'.pdb'
     write(FileRG ,'(a,i5  )') './Analy/RofG',NumF
   else
     write(*,*) 'ERROR :: too many files / not supported !'
     call Finalize
   end if

   open(51,file=FilePDB,status='unknown')

   write(51,'(a,f11.3,a)') 'REMARK   1  Time = ',Timeps,' ps'

   ii = Numa + (Nlipid-1)*NAcomp

   do i = 1 , NAcomp

     j = ii + i
     Ch = AtomName(j)(4:4)
     write(AName,'(a1,a3)') Ch,AtomName(j)(1:3)

     write(51,'(a4,2x,i5,x,a4,x,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
     & 'ATOM',i,AName,ResidName(j),one,R(:,j),zero,zero,AtomName(j)(1:1)

   end do

   write(51,'(a)') 'END'

   close(51)

   open(52,file=FileRG,status='unknown')

   do i = 1 , NumGroup
     write(52,'(4f8.3)') Rcom_Group(:,i,Nlipid), Rg_Group(i,Nlipid)
   end do

   close(52)

end subroutine Write_SingleLipid


!######################################################################
!######################################################################


subroutine DefineGroup(Numa)

use Numbers, only : NumMol, NumAtm
use CGdist
use AtomParam, only : MolName, AtomName, Mass

implicit none

character(len=4) :: AName
integer :: i, j, k, Numa


   NumGroup = 13
   allocate(NumAtm_Group(NumGroup))
   allocate(Mass_Group(NumGroup))
   allocate(Rcom_Group(3,NumGroup,NMcomp))
   allocate(Rg_Group(NumGroup,NMcomp))

   if(MolName(Kcomp)=='PhPC'.or.MolName(Kcomp)=='DPPC') then

   do i = 1, NAcomp
     j = Numa + i
     AName = AtomName(j)

     if((trim(AName) == 'N'   ).or.&
     &  (trim(AName) == 'C11' ).or.&
     &  (trim(AName) == 'H11A').or.&
     &  (trim(AName) == 'H11B').or.&
     &  (trim(AName) == 'C12' ).or.&
     &  (trim(AName) == 'H12A').or.&
     &  (trim(AName) == 'H12B').or.&
     &  (trim(AName) == 'C13' ).or.&
     &  (trim(AName) == 'H13A').or.&
     &  (trim(AName) == 'H13B').or.&
     &  (trim(AName) == 'H13C').or.&
     &  (trim(AName) == 'C14' ).or.&
     &  (trim(AName) == 'H14A').or.&
     &  (trim(AName) == 'H14B').or.&
     &  (trim(AName) == 'H14C').or.&
     &  (trim(AName) == 'C15' ).or.&
     &  (trim(AName) == 'H15A').or.&
     &  (trim(AName) == 'H15B').or.&
     &  (trim(AName) == 'H15C')) then

       GroupNumber(i) = 1

     else if((trim(AName) == 'P'   ).or.&
     &       (trim(AName) == 'O11' ).or.&
     &       (trim(AName) == 'O12' ).or.&
     &       (trim(AName) == 'O13' ).or.&
     &       (trim(AName) == 'O14' )) then

       GroupNumber(i) = 2

     else if((trim(AName) == 'C1'  ).or.&
     &       (trim(AName) == 'HA'  ).or.&
     &       (trim(AName) == 'HB'  ).or.&
     &       (trim(AName) == 'C2'  ).or.&
     &       (trim(AName) == 'HS'  ).or.&
     &       (trim(AName) == 'C3'  ).or.&
     &       (trim(AName) == 'HX'  ).or.&
     &       (trim(AName) == 'HY'  )) then

       GroupNumber(i) = 3

     else if((trim(AName) == 'C21' ).or.&
     &       (trim(AName) == 'O21' ).or.&
     &       (trim(AName) == 'O22' )) then

       GroupNumber(i) = 4

     else if((trim(AName) == 'C22' ).or.&
     &       (trim(AName) == 'H2R' ).or.&
     &       (trim(AName) == 'H2S' ).or.&
     &       (trim(AName) == 'C23' ).or.&
     &       (trim(AName) == 'H3R' ).or.&
     &       (trim(AName) == 'H3S' ).or.&
     &       (trim(AName) == 'C24' ).or.&
     &       (trim(AName) == 'H4R' ).or.&
     &       (trim(AName) == 'H4S' ).or.&
     &       (trim(AName) == 'C25' ).or.&
     &       (trim(AName) == 'H5R' ).or.&
     &       (trim(AName) == 'H5S' ).or.&
     &       (trim(AName) == 'C217').or.&
     &       (trim(AName) == 'H17R').or.&
     &       (trim(AName) == 'H17S').or.&
     &       (trim(AName) == 'H17T')) then

       GroupNumber(i) = 5

     else if((trim(AName) == 'C26' ).or.&
     &       (trim(AName) == 'H6R' ).or.&
     &       (trim(AName) == 'H6S' ).or.&
     &       (trim(AName) == 'C27' ).or.&
     &       (trim(AName) == 'H7R' ).or.&
     &       (trim(AName) == 'H7S' ).or.&
     &       (trim(AName) == 'C28' ).or.&
     &       (trim(AName) == 'H8R' ).or.&
     &       (trim(AName) == 'H8S' ).or.&
     &       (trim(AName) == 'C29' ).or.&
     &       (trim(AName) == 'H9R' ).or.&
     &       (trim(AName) == 'H9S' ).or.&
     &       (trim(AName) == 'C218').or.&
     &       (trim(AName) == 'H18R').or.&
     &       (trim(AName) == 'H18S').or.&
     &       (trim(AName) == 'H18T')) then

       GroupNumber(i) = 6

     else if((trim(AName) == 'C210').or.&
     &       (trim(AName) == 'H10R').or.&
     &       (trim(AName) == 'H10S').or.&
     &       (trim(AName) == 'C211').or.&
     &       (trim(AName) == 'H11R').or.&
     &       (trim(AName) == 'H11S').or.&
     &       (trim(AName) == 'C212').or.&
     &       (trim(AName) == 'H12R').or.&
     &       (trim(AName) == 'H12S').or.&
     &       (trim(AName) == 'C213').or.&
     &       (trim(AName) == 'H13R').or.&
     &       (trim(AName) == 'H13S').or.&
     &       (trim(AName) == 'C219').or.&
     &       (trim(AName) == 'H19R').or.&
     &       (trim(AName) == 'H19S').or.&
     &       (trim(AName) == 'H19T')) then

       GroupNumber(i) = 7

     else if((trim(AName) == 'C214').or.&
     &       (trim(AName) == 'H14R').or.&
     &       (trim(AName) == 'H14S').or.&
     &       (trim(AName) == 'C215').or.&
     &       (trim(AName) == 'H15R').or.&
     &       (trim(AName) == 'H15S').or.&
     &       (trim(AName) == 'C216').or.&
     &       (trim(AName) == 'H16R').or.&
     &       (trim(AName) == 'H16S').or.&
     &       (trim(AName) == 'H16T').or.&
     &       (trim(AName) == 'C220').or.&
     &       (trim(AName) == 'H20R').or.&
     &       (trim(AName) == 'H20S').or.&
     &       (trim(AName) == 'H20T')) then

       GroupNumber(i) = 8

     else if((trim(AName) == 'C31' ).or.&
     &       (trim(AName) == 'O31' ).or.&
     &       (trim(AName) == 'O32' )) then

       GroupNumber(i) = 9

     else if((trim(AName) == 'C32' ).or.&
     &       (trim(AName) == 'H2X' ).or.&
     &       (trim(AName) == 'H2Y' ).or.&
     &       (trim(AName) == 'C33' ).or.&
     &       (trim(AName) == 'H3X' ).or.&
     &       (trim(AName) == 'H3Y' ).or.&
     &       (trim(AName) == 'C34' ).or.&
     &       (trim(AName) == 'H4X' ).or.&
     &       (trim(AName) == 'H4Y' ).or.&
     &       (trim(AName) == 'C35' ).or.&
     &       (trim(AName) == 'H5X' ).or.&
     &       (trim(AName) == 'H5Y' ).or.&
     &       (trim(AName) == 'C317').or.&
     &       (trim(AName) == 'H17X').or.&
     &       (trim(AName) == 'H17Y').or.&
     &       (trim(AName) == 'H17Z')) then

       GroupNumber(i) = 10

     else if((trim(AName) == 'C36' ).or.&
     &       (trim(AName) == 'H6X' ).or.&
     &       (trim(AName) == 'H6Y' ).or.&
     &       (trim(AName) == 'C37' ).or.&
     &       (trim(AName) == 'H7X' ).or.&
     &       (trim(AName) == 'H7Y' ).or.&
     &       (trim(AName) == 'C38' ).or.&
     &       (trim(AName) == 'H8X' ).or.&
     &       (trim(AName) == 'H8Y' ).or.&
     &       (trim(AName) == 'C39' ).or.&
     &       (trim(AName) == 'H9X' ).or.&
     &       (trim(AName) == 'H9Y' ).or.&
     &       (trim(AName) == 'C318').or.&
     &       (trim(AName) == 'H18X').or.&
     &       (trim(AName) == 'H18Y').or.&
     &       (trim(AName) == 'H18Z')) then

       GroupNumber(i) = 11

     else if((trim(AName) == 'C310').or.&
     &       (trim(AName) == 'H10X').or.&
     &       (trim(AName) == 'H10Y').or.&
     &       (trim(AName) == 'C311').or.&
     &       (trim(AName) == 'H11X').or.&
     &       (trim(AName) == 'H11Y').or.&
     &       (trim(AName) == 'C312').or.&
     &       (trim(AName) == 'H12X').or.&
     &       (trim(AName) == 'H12Y').or.&
     &       (trim(AName) == 'C313').or.&
     &       (trim(AName) == 'H13X').or.&
     &       (trim(AName) == 'H13Y').or.&
     &       (trim(AName) == 'C319').or.&
     &       (trim(AName) == 'H19X').or.&
     &       (trim(AName) == 'H19Y').or.&
     &       (trim(AName) == 'H19Z')) then

       GroupNumber(i) = 12

     else if((trim(AName) == 'C314').or.&
     &       (trim(AName) == 'H14X').or.&
     &       (trim(AName) == 'H14Y').or.&
     &       (trim(AName) == 'C315').or.&
     &       (trim(AName) == 'H15X').or.&
     &       (trim(AName) == 'H15Y').or.&
     &       (trim(AName) == 'C316').or.&
     &       (trim(AName) == 'H16X').or.&
     &       (trim(AName) == 'H16Y').or.&
     &       (trim(AName) == 'H16Z').or.&
     &       (trim(AName) == 'C320').or.&
     &       (trim(AName) == 'H20X').or.&
     &       (trim(AName) == 'H20Y').or.&
     &       (trim(AName) == 'H20Z')) then

       GroupNumber(i) = 13

     else

       write(*,*) 'ERROR : no atom'
       write(*,*) AName
       call Finalize

     end if

   end do

   else if(MolName(Kcomp)=='DMPC') then

   do i = 1, NAcomp
     j = Numa + i
     AName = AtomName(j)

     if((trim(AName) == 'N'   ).or.&
     &  (trim(AName) == 'C11' ).or.&
     &  (trim(AName) == 'H11A').or.&
     &  (trim(AName) == 'H11B').or.&
     &  (trim(AName) == 'C12' ).or.&
     &  (trim(AName) == 'H12A').or.&
     &  (trim(AName) == 'H12B').or.&
     &  (trim(AName) == 'C13' ).or.&
     &  (trim(AName) == 'H13A').or.&
     &  (trim(AName) == 'H13B').or.&
     &  (trim(AName) == 'H13C').or.&
     &  (trim(AName) == 'C14' ).or.&
     &  (trim(AName) == 'H14A').or.&
     &  (trim(AName) == 'H14B').or.&
     &  (trim(AName) == 'H14C').or.&
     &  (trim(AName) == 'C15' ).or.&
     &  (trim(AName) == 'H15A').or.&
     &  (trim(AName) == 'H15B').or.&
     &  (trim(AName) == 'H15C')) then

       GroupNumber(i) = 1

     else if((trim(AName) == 'P'   ).or.&
     &       (trim(AName) == 'O11' ).or.&
     &       (trim(AName) == 'O12' ).or.&
     &       (trim(AName) == 'O13' ).or.&
     &       (trim(AName) == 'O14' )) then

       GroupNumber(i) = 2

     else if((trim(AName) == 'C1'  ).or.&
     &       (trim(AName) == 'HA'  ).or.&
     &       (trim(AName) == 'HB'  ).or.&
     &       (trim(AName) == 'C2'  ).or.&
     &       (trim(AName) == 'HS'  ).or.&
     &       (trim(AName) == 'C3'  ).or.&
     &       (trim(AName) == 'HX'  ).or.&
     &       (trim(AName) == 'HY'  )) then

       GroupNumber(i) = 3

     else if((trim(AName) == 'C21' ).or.&
     &       (trim(AName) == 'O21' ).or.&
     &       (trim(AName) == 'O22' ).or.&
     &       (trim(AName) == 'C22' ).or.&
     &       (trim(AName) == 'H2R' ).or.&
     &       (trim(AName) == 'H2S' )) then

       GroupNumber(i) = 4

     else if((trim(AName) == 'C23' ).or.&
     &       (trim(AName) == 'H3R' ).or.&
     &       (trim(AName) == 'H3S' ).or.&
     &       (trim(AName) == 'C24' ).or.&
     &       (trim(AName) == 'H4R' ).or.&
     &       (trim(AName) == 'H4S' ).or.&
     &       (trim(AName) == 'C25' ).or.&
     &       (trim(AName) == 'H5R' ).or.&
     &       (trim(AName) == 'H5S' )) then

       GroupNumber(i) = 5

     else if((trim(AName) == 'C26' ).or.&
     &       (trim(AName) == 'H6R' ).or.&
     &       (trim(AName) == 'H6S' ).or.&
     &       (trim(AName) == 'C27' ).or.&
     &       (trim(AName) == 'H7R' ).or.&
     &       (trim(AName) == 'H7S' ).or.&
     &       (trim(AName) == 'C28' ).or.&
     &       (trim(AName) == 'H8R' ).or.&
     &       (trim(AName) == 'H8S' )) then

       GroupNumber(i) = 6

     else if((trim(AName) == 'C29' ).or.&
     &       (trim(AName) == 'H9R' ).or.&
     &       (trim(AName) == 'H9S' ).or.&
     &       (trim(AName) == 'C210').or.&
     &       (trim(AName) == 'H10R').or.&
     &       (trim(AName) == 'H10S').or.&
     &       (trim(AName) == 'C211').or.&
     &       (trim(AName) == 'H11R').or.&
     &       (trim(AName) == 'H11S')) then

       GroupNumber(i) = 7

     else if((trim(AName) == 'C212').or.&
     &       (trim(AName) == 'H12R').or.&
     &       (trim(AName) == 'H12S').or.&
     &       (trim(AName) == 'C213').or.&
     &       (trim(AName) == 'H13R').or.&
     &       (trim(AName) == 'H13S').or.&
     &       (trim(AName) == 'C214').or.&
     &       (trim(AName) == 'H14R').or.&
     &       (trim(AName) == 'H14S').or.&
     &       (trim(AName) == 'H14T')) then

       GroupNumber(i) = 8

     else if((trim(AName) == 'C31' ).or.&
     &       (trim(AName) == 'O31' ).or.&
     &       (trim(AName) == 'O32' ).or.&
     &       (trim(AName) == 'C32' ).or.&
     &       (trim(AName) == 'H2X' ).or.&
     &       (trim(AName) == 'H2Y' )) then

       GroupNumber(i) = 9

     else if((trim(AName) == 'C33' ).or.&
     &       (trim(AName) == 'H3X' ).or.&
     &       (trim(AName) == 'H3Y' ).or.&
     &       (trim(AName) == 'C34' ).or.&
     &       (trim(AName) == 'H4X' ).or.&
     &       (trim(AName) == 'H4Y' ).or.&
     &       (trim(AName) == 'C35' ).or.&
     &       (trim(AName) == 'H5X' ).or.&
     &       (trim(AName) == 'H5Y' )) then

       GroupNumber(i) = 10

     else if((trim(AName) == 'C36' ).or.&
     &       (trim(AName) == 'H6X' ).or.&
     &       (trim(AName) == 'H6Y' ).or.&
     &       (trim(AName) == 'C37' ).or.&
     &       (trim(AName) == 'H7X' ).or.&
     &       (trim(AName) == 'H7Y' ).or.&
     &       (trim(AName) == 'C38' ).or.&
     &       (trim(AName) == 'H8X' ).or.&
     &       (trim(AName) == 'H8Y' )) then

       GroupNumber(i) = 11

     else if((trim(AName) == 'C39' ).or.&
     &       (trim(AName) == 'H9X' ).or.&
     &       (trim(AName) == 'H9Y' ).or.&
     &       (trim(AName) == 'C310').or.&
     &       (trim(AName) == 'H10X').or.&
     &       (trim(AName) == 'H10Y').or.&
     &       (trim(AName) == 'C311').or.&
     &       (trim(AName) == 'H11X').or.&
     &       (trim(AName) == 'H11Y')) then

       GroupNumber(i) = 12

     else if((trim(AName) == 'C312').or.&
     &       (trim(AName) == 'H12X').or.&
     &       (trim(AName) == 'H12Y').or.&
     &       (trim(AName) == 'C313').or.&
     &       (trim(AName) == 'H13X').or.&
     &       (trim(AName) == 'H13Y').or.&
     &       (trim(AName) == 'C314').or.&
     &       (trim(AName) == 'H14X').or.&
     &       (trim(AName) == 'H14Y').or.&
     &       (trim(AName) == 'H14Z')) then

       GroupNumber(i) = 13

     else

       write(*,*) 'ERROR : no atom'
       write(*,*) AName
       call Finalize

     end if

   end do

   else

     write(*,*) 'error : no lipid'
     stop

   end if

   NumAtm_Group = 0
   Mass_Group = 0.d0

   do i = 1, NAcomp
     j = Numa + i
     k = GroupNumber(i)
     NumAtm_Group(k) = NumAtm_Group(k) + 1
     Mass_Group(k)   = Mass_Group(k)   + Mass(j)
   end do

end subroutine DefineGroup
