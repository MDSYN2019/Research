

!######################################################################
!######################################################################


! **************************
! ** BKS Parameter File   **
! **************************

subroutine Read_BKS_Parameter

use CommonBlocks, only : Qdebug
use FFParameters
use IOparam, only : Parameter_file
use UnitExParam, only : reng, e, Avogadro

implicit none

integer :: i
character(len=72) :: String, String1
character(len=1) :: Wk1
character(len=1) :: Wk2
character(len=1) :: Wk3
character(len=1) :: Wk4
real(8), parameter :: ScaleParam = e * reng ! unit exchanging rate

! parameters

open(1,file=trim(Parameter_file),status='old')

   Wk1 = '*'
   Wk2 = '!'
   Wk3 = ' '
   Wk4 = '	'

   do

     read(1,'(72a)') String

     if(String(1:7)=='NONBOND') exit

   end do

   NumBondParam=0
   NumAngleParam=0
   NumDihedralParam=0
   NumImproperParam=0
   NumBKSParam = 0

   do

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:3)=='END') exit

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumBKSParam = NumBKSParam + 1

     read(String,*) BKSAtomI(NumBKSParam),BKSAtomJ(NumBKSParam),      &
     &              BKS_Aij(NumBKSParam),BKS_Bij(NumBKSParam),BKS_Cij(NumBKSParam)

   end do


! ----
! Unit
! ----

   do i = 1, NumBKSParam

     if(Qdebug) then
     print *, BKSAtomI(i),BKSAtomJ(i),BKS_Aij(i)*e/4184.*Avogadro,BKS_Bij(i),&
     & BKS_Cij(i)*e/4184.*Avogadro
     end if

     BKS_Aij(i) = BKS_Aij(i) * ScaleParam
     BKS_Cij(i) = BKS_Cij(i) * ScaleParam


   end do

! ----------------------------------------------------------------------

close(1)

end subroutine Read_BKS_Parameter


!######################################################################
!######################################################################


subroutine AllocateBKS

use Numbers, only : N, NumSpec, NumMol, NumAtm
use CommonBlocks, only : QMaster, Qstdout, Qdebug, QCoulomb
use FFParameters
use UnitExParam, only : Avogadro, ec
use EwaldParam, only : Nel, Nelist, PCh
use NonbondParam, only : Charge, CoeA, CoeB, CoeC, BKStype, PairBack, IJPair
use BondedParam, only : NumBond, NumAngle, NumUB, NumDihedral, NumImproper
use AtomParam, only : AtomName, ResidName, Mass, InvMass

implicit none

integer :: ii, i, j, k, jj, ij
character(len=4), dimension(4) :: Name
character(len=4), dimension(4) :: TmpAName, TmpRName
integer, dimension(:), allocatable :: IspecRep
real(8) :: TotalCharge
integer :: ItotalC

   NumBond     = 0
   NumAngle    = 0
   NumUB       = 0
   NumDihedral = 0
   NumImproper = 0

   allocate( Mass   (N) )
   allocate( InvMass(N) )
   allocate( BKStype(N) )
   allocate( Charge (N) )

   allocate( IspecRep(NumSpec) )
   allocate( PairBack(NumSpec,NumSpec) )

   IJPair = NumSpec + NumSpec*(NumSpec-1)/2

   allocate( CoeA(IJPair) )
   allocate( CoeB(IJPair) )
   allocate( CoeC(IJPair) )

   ii = 0

   do i = 1, NumSpec

     jj = 0

     do j = 1, NumMol(i)

       do k = 1, NumAtm(i)

         ii = ii + 1
         jj = jj + 1
         BKStype(ii) = i

         if(jj==1) then
           IspecRep(i) = ii
         end if

       end do

     end do

   end do

   do i = 1 , N

     TmpAName (1) = AtomName (i)
     TmpRName (1) = ResidName(i)

!   ------------------
     call NameExch(1)
!   ------------------

     do j = 1, NumResidueParam

       if(ResiNameParam(j) == ResidName(i)) then

inn0:    do k = 1, NumAtom_inResi(j)

           if(AtomNameParam(k,j) == AtomName(i)) then

             Charge(i) = ChargeParam(k,j)
             exit inn0

           end if

           if(k==NumAtom_inResi(j)) then
             if(QMaster) then
               write(*,*) 'ERROR : missing parameter (charge)'
               write(*,*) 'atom = ',i
             end if
             call Finalize
           end if

         end do inn0

       end if

     end do

     do j = 1, NumAtomTypeParam

       if( Name(1) == AtomTypeParam(j) ) then

         Mass(i) = MassParam(j)

         exit

       end if

       if( j == NumAtomTypeParam ) then

         if(QMaster) write(*,*) 'No Mass parameter was found for',Name(1)
         if(QMaster) write(*,*) 'Atom Number =' , i

         call Finalize

       end if

     end do

   end do

! ## checking charge neutrality

   QCoulomb = .False.

   do i = 1, N
     if(Charge(i)/=0.) then
       QCoulomb = .True.
       exit
     end if
   end do

   if(QCoulomb) then

     TotalCharge = Sum( Charge )
     ItotalC = nint( TotalCharge * 1.d+5 )
     TotalCharge = dble(ItotalC) * 1.d-5

     Charge = Charge * sqrt(ec)

     Nel = 0
     do i = 1, N
       if(Charge(i)/=0.) Nel = Nel + 1
     end do

     if(Nel/=0) then
       allocate( Nelist(Nel) )
       allocate( PCh(Nel) )
       ii = 0
       do i = 1, N
         if(Charge(i)/=0.) then
           ii = ii + 1
           Nelist(ii) = i
           PCh(ii) = Charge(i)
         end if
       end do
     end if
   end if

   if(QMaster.and.Qstdout) then
     if(QCoulomb) then
       write(*,*) 'Checking charge neutrality: total charge =',TotalCharge
     else
       write(*,*) 'No Coulomb interaction'
     end if
   end if

! ----
! Unit
! ----
   Mass = Mass * 1.d-3 / Avogadro

   InvMass = 1.d0 / Mass

! ----------------------------------------------------------------------

   ij = 0

   do i = 1, NumSpec

     ii = IspecRep(i)

     do j = i, NumSpec

       jj = IspecRep(j)

       ij = ij + 1

       PairBack(i,j) = ij
       PairBack(j,i) = ij

       TmpAName(1) = AtomName(ii)
       TmpAName(2) = AtomName(jj)

       TmpRName(1) = ResidName(ii)
       TmpRName(2) = ResidName(jj)

!     ------------------
       call NameExch(2)
!     ------------------

       do k = 1, NumBKSParam

         if(( (BKSAtomI(k) == Name(1) )  &
      & .and. (BKSAtomJ(k) == Name(2)) ) &
      & .or.( (BKSAtomI(k) == Name(2) )  &
      & .and. (BKSAtomJ(k) == Name(1)) )) then

           CoeA(ij) = BKS_Aij(k)
           CoeB(ij) = BKS_Bij(k)
           CoeC(ij) = BKS_Cij(k)
           exit

         end if

         if( k == NumBKSParam ) then

           if(QMaster) then

#ifndef BMONI
             write(11,*) 'lost parameter : BKS'
#endif
             write( 6,*) 'lost parameter : BKS'
             write( 6,*) ' Pair=',i,j
             write( 6,*) 'TmpRName=',TmpRName(1),TmpRName(2)
             write( 6,*) 'TmpAName=',TmpAName(1),TmpAName(2)
             write( 6,*) 'Name=',Name(1),Name(2)

           end if

           call Finalize

         end if

       end do

     end do

   end do

   if(Qdebug) then
   do i = 1, NumSpec
   do j = 1, NumSpec
     ij = PairBack(i,j)
     write(*,*) 'exp6',i,j,CoeA(ij),CoeB(ij),CoeC(ij)
   end do
   end do
   end if

! --------------------------------------------------------------

Contains

   subroutine NameExch(Na)

   integer :: j0, k0, l0, Na, count1

     do j0 = 1, Na

       count1 = 0

       do k0 = 1, NumResidueParam

         if( TmpRName(j0) == ResiNameParam(k0) ) then

           count1 = count1 + 1

           do l0 = 1, NumAtom_inResi(k0)

             if(TmpAName(j0) == AtomNameParam(l0,k0)) then

               Name(j0) = AtomType(l0,k0)

             end if

           end do

         end if

       end do

       if( count1 /= 1 ) then

         if(QMaster) then

           if( Na == 1 ) then

             write(*,'(a/a,a/a,a/a,i6)')             &
             &          'ERROR: no residue (ATOM)',  &
             &          'Residue = ',TmpRName(1),    &
             &          'Atom    = ',TmpAName(1),    &
             &          'count   = ',count1

           else if( Na == 2 ) then

             write(*,'(a/a,a/a,a/a,i6)')             &
             &          'ERROR: no residue (BOND)',  &
             &          'Residue = ',TmpRName(j0),    &
             &          'Atom    = ',TmpAName(j0),    &
             &          'count   = ',count1

           else if( Na == 3 ) then

             write(*,'(a/a,a/a,a/a,i6)')              &
             &          'ERROR: no residue (ANGLE)',  &
             &          'Residue = ',TmpRName(j0),     &
             &          'Atom    = ',TmpAName(j0),     &
             &          'count   = ',count1

           else if( Na == 4 ) then

             write(*,'(a/a,a/a,a/a,i6)')                           &
             &          'ERROR: no residue (DIHEDRAL or IMPROPER)',&
             &          'Residue = ',TmpRName(j0),                  &
             &          'Atom    = ',TmpAName(j0),                  &
             &          'count   = ',count1

           end if

         end if

         call Finalize

       end if

     end do

   end subroutine NameExch

! ----------------------------------------------------------------------

end subroutine AllocateBKS


!#####################################################################
!#####################################################################


subroutine ForceR_BKS

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QPathInt
use BookParam, only : Npair, ListIJ
use EwaldParam, only : Alpha, ar2
use NonbondParam, only : Charge, CoeA, CoeB, CoeC, BKStype, PairBack, &
&   Frc_Ersp, Vir_Ersp, Ene_LJ, Ene_Ersp
use CutoffParam, only : Rcutoff2
use CellParam, only : CellShft

implicit none

integer :: l, i, j, k, Itype, Jtype, ij
real(8), dimension(3) :: Rij, FijLJ, FijEL, Fij
real(8) :: R2, R1, InvR2, A_ij, B_ij, C_ij
real(8) :: SR6, SR12, fkLJ, cf, x, ErrorFunc
real(8) :: xtm, fk1, fk2, fk
real(8), dimension(3,N) :: Forcer
real(8), dimension(3,-13:13) :: Gk
real(8) :: Error_Function
external Error_Function

   Forcer = 0.d0
   Gk = 0.d0

   do l = 1 , Npair

     i = ListIJ(1,l)
     j = ListIJ(2,l)
     k = ListIJ(3,l)

     Rij(:) = R(:,i) - R(:,j) + CellShft(:,k)
     R2 = Rij(1) * Rij(1) + Rij(2) * Rij(2) + Rij(3) * Rij(3)

     if(R2 <= Rcutoff2) then

       InvR2 = 1.d0 / R2

       Itype = BKStype(i)
       Jtype = BKStype(j)
       ij = PairBack(Itype,Jtype)

       A_ij = CoeA(ij)
       B_ij = CoeB(ij)
       C_ij = CoeC(ij)

       R1  = sqrt( R2 )

       if((A_ij/=0.).and.(C_ij/=0.)) then

         SR6  = - InvR2 * InvR2 * InvR2 * C_ij
         SR12 = A_ij * exp( - B_ij * R1 )

         fkLJ = B_ij * SR12 * R1 + 6.d0 * SR6        !4e(12()^12-6()^6)
         fkLJ = fkLJ * InvR2

         FijLJ = fkLJ * Rij

         Ene_LJ = Ene_LJ + SR12 + SR6

       else

         FijLJ = 0.d0

       end if

       cf = Charge(i) * Charge(j)

       if(cf /= 0.) then

         x   = Alpha * R1

         ErrorFunc = Error_Function(x)

         xtm = -x * x
         fk1 = cf * ErrorFunc / R1
         fk2 = cf * ar2 * exp(xtm)
         fk  = ( fk1 + fk2 ) * InvR2

         FijEL = fk * Rij

         Ene_Ersp = Ene_Ersp + fk1

       else

         FijEL = 0.d0

       end if

       Fij = FijLJ + FijEL

       Forcer(:,i) = Forcer(:,i) + Fij
       Forcer(:,j) = Forcer(:,j) - Fij
       Gk(:,k)     = Gk(:,k)     + Fij

     end if

   end do

   if(QPathInt) then
     call VirialBekkerPI(Forcer,Vir_Ersp,Gk)
   else
     call VirialBekker(Forcer,Vir_Ersp,Gk)
   end if

   Frc_Ersp = Frc_Ersp + Forcer


end subroutine ForceR_BKS
