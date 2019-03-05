

!######################################################################
!######################################################################


! >> F monitor ##

subroutine ListFmoni

use CommonBlocks, only : QSHAKE, QRigidBody, QPathInt
use CommonMPI
use CommonPI
use F_monitor
use RBparam
use SHAKEparam
use BondedParam, only : NumBond, NumAngle, NumDihedral, BondI, BondJ, &
&   AngleI, AngleJ, AngleK, DihedI, DihedJ, DihedK, DihedL, vdWSubtDih

implicit none

integer :: i, j, k, l, ii, jj, kk
integer :: Nas
logical :: flag
integer :: NProcsTemp, MyRankTemp
integer :: IRB, JRB

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp + NiniF

! ## Bond Stretching

   NumBondgA = 0
   do ii = 1 , NumBond
     i = BondI(ii)
     j = BondJ(ii)
     if(((i > NiniF).and.(i <= NfinF)).and. &
     &  ((j > NiniF).and.(j <= NfinF))) then
       NumBondgA = NumBondgA + 1
     end if
   end do

   allocate( BondIgA(NumBondgA) )
   allocate( BondJgA(NumBondgA) )

   jj = 0
   do ii = 1 , NumBond
     i = BondI(ii)
     j = BondJ(ii)
     if(((i > NiniF).and.(i <= NfinF)).and. &
     &  ((j > NiniF).and.(j <= NfinF))) then
       jj = jj + 1
       BondIgA(jj) = i
       BondJgA(jj) = j
     end if
   end do

! ## Bending

   jj = 0

   do ii = 1, NumAngle
     i = AngleI(ii)
     j = AngleJ(ii)
     k = AngleK(ii)
     if(((i > NiniF).and.(i <= NfinF)).and. &
     &  ((j > NiniF).and.(j <= NfinF)).and. &
     &  ((k > NiniF).and.(k <= NfinF))) then
       jj = jj + 1
     end if
   end do
   NumAnglegA = jj

   allocate( AngleIgA(NumAnglegA) )
   allocate( AngleJgA(NumAnglegA) )
   allocate( AngleKgA(NumAnglegA) )

   jj = 0

   do ii = 1, NumAngle
     i = AngleI(ii)
     j = AngleJ(ii)
     k = AngleK(ii)
     if(((i > NiniF).and.(i <= NfinF)).and. &
     &  ((j > NiniF).and.(j <= NfinF)).and. &
     &  ((k > NiniF).and.(k <= NfinF))) then
       jj = jj + 1
       AngleIgA(jj) = i
       AngleJgA(jj) = j
       AngleKgA(jj) = k
     end if
   end do

! ## Torsion

   jj = 0

   do ii = 1, NumDihedral
     i = DihedI(ii)
     j = DihedJ(ii)
     k = DihedK(ii)
     l = DihedL(ii)
     if(((i > NiniF).and.(i <= NfinF)).and. &
     &  ((j > NiniF).and.(j <= NfinF)).and. &
     &  ((k > NiniF).and.(k <= NfinF)).and. &
     &  ((l > NiniF).and.(l <= NfinF))) then
       jj = jj + 1
     end if
   end do
   NumDihedralgA = jj

   allocate( DihedIgA(NumDihedralgA) )
   allocate( DihedJgA(NumDihedralgA) )
   allocate( DihedKgA(NumDihedralgA) )
   allocate( DihedLgA(NumDihedralgA) )
   allocate( vdWSubtDihgA(NumDihedralgA) )

   vdWSubtDihgA = .False.

   jj = 0

   do ii = 1, NumDihedral
     i = DihedI(ii)
     j = DihedJ(ii)
     k = DihedK(ii)
     l = DihedL(ii)
     if(((i > NiniF).and.(i <= NfinF)).and. &
     &  ((j > NiniF).and.(j <= NfinF)).and. &
     &  ((k > NiniF).and.(k <= NfinF)).and. &
     &  ((l > NiniF).and.(l <= NfinF))) then
       jj = jj + 1
       DihedIgA(jj) = i
       DihedJgA(jj) = j
       DihedKgA(jj) = k
       DihedLgA(jj) = l
       if(vdWSubtDih(ii)) vdWSubtDihgA(jj) = .True.
     end if
   end do

! ## SHAKE

   if(QSHAKE) then

     l = 0

     do k = 1 , NSHAKEGroup

       do kk = 1 , NCoupleBond(k)

         ii = CouplePair(k,kk,1)
         jj = CouplePair(k,kk,2)

         i  = CoupleAtom(k,ii)
         j  = CoupleAtom(k,jj)

         if(((i > NiniF).and.(i <= NfinF)).and. &
         &  ((j > NiniF).and.(j <= NfinF))) then
           l = l + 1
         end if

       end do

     end do

     NumSHK = l

     allocate( SKpairI(NumSHK) )
     allocate( SKpairJ(NumSHK) )

     l = 0

     do k = 1 , NSHAKEGroup

       do kk = 1 , NCoupleBond(k)

         ii = CouplePair(k,kk,1)
         jj = CouplePair(k,kk,2)

         i  = CoupleAtom(k,ii)
         j  = CoupleAtom(k,jj)

         if(((i > NiniF).and.(i <= NfinF)).and. &
         &  ((j > NiniF).and.(j <= NfinF))) then

           l = l + 1
           SKpairI(l) = i
           SKpairJ(l) = j

         end if

       end do

     end do

   end if

! ## NonBond 

   NlistgA = 0

   do i = Nas, NfinF, NProcsTemp

     do j = i-2 , NiniF+1, -2

       call SubtPair(i,j,flag)

       if(flag) then
         NlistgA = NlistgA + 1
         ListIgA(NlistgA) = i
         ListJgA(NlistgA) = j
       end if

     end do

     do j = i+1 , NfinF, 2

       call SubtPair(i,j,flag)

       if(flag) then
         NlistgA = NlistgA + 1
         ListIgA(NlistgA) = i
         ListJgA(NlistgA) = j
       end if

     end do

   end do

Contains

   subroutine SubtPair(ll,mm,Qflag)

   logical :: Qflag
   integer :: ll, mm

      Qflag = .True.

      do k = 1 , NumBondgA
        ii = BondIgA(k)
        jj = BondJgA(k)
        if( ((ll==ii).and.(mm==jj)) .or. ((ll==jj).and.(mm==ii)) ) then
          Qflag = .False.
        end if
      end do

      do k = 1, NumAnglegA
        ii = AngleIgA(k)
        jj = AngleKgA(k)
        if( ((ll==ii).and.(mm==jj)) .or. ((ll==jj).and.(mm==ii)) ) then
          Qflag = .False.
        end if
      end do

      do k = 1 , NumDihedralgA
        ii = DihedIgA(k)
        jj = DihedLgA(k)
        if( ((ll==ii).and.(mm==jj)) .or. ((ll==jj).and.(mm==ii)) ) then
          Qflag = .False.
        end if
      end do

      if(QSHAKE) then
        do k = 1, NumSHK
          ii = SKpairI(k)
          jj = SKpairJ(k)
          if( ((ll==ii).and.(mm==jj)) .or. ((ll==jj).and.(mm==ii)) ) then
            Qflag = .False.
          end if
        end do
      end if

      if(QRigidBody) then
        IRB = AtomUnitNum(ll)
        JRB = AtomUnitNum(mm)
        if(IRB == JRB) Qflag = .False.
      end if

   end subroutine SubtPair

end subroutine ListFmoni


!######################################################################
!######################################################################


subroutine Force_gA

use CommonBlocks, only : QPathInt, QSwitch
use Configuration, only : R
use CommonMPI
use CommonPI
use F_monitor
use NonbondParam, only : Charge, Rminh, EpsLJ, Rminh14, Eps14
use BondedParam, only : Frc_Bond, Frc_Angle, Frc_UB, Frc_Dihed, Frc_Impro
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j, l
real(8) :: R1, R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk
real(8) :: ek
real(8) :: cf
real(8), dimension(3) :: Rij
real(8), dimension(3) :: FijLJ, FijEL
integer :: Nas
real(8) :: SC2, SC6, SC12, fk14
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Fint = 0.d0

   do l = 1 , NlistgA

     i = ListIgA(l)
     j = ListJgA(l)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product( Rij, Rij )

     InvR2 = 1.d0 / R2

     cf = Charge(i) * Charge(j)

     if(cf /= 0.) then

       R1  = sqrt( R2 )

       fk1 = cf / R1
       fk  = fk1 * InvR2

       FijEL = fk * Rij

       Fint(:,i) = Fint(:,i) + FijEL
       Fint(:,j) = Fint(:,j) - FijEL

     end if

     if(R2 <= Rcutoff2) then

       Sgm   = Rminh(i) + Rminh(j)
       Sgm2  = Sgm * Sgm
       Eps   = EpsLJ(i) * EpsLJ(j)

       SR2  = Sgm2 * InvR2                    !(sigma/r)^2
       SR6  = SR2 * SR2 * SR2                 !         ^6
       SR12 = SR6 * SR6                       !         ^12
       fk   = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
       ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

       FijLJ = fk * Rij

       Fint(:,i) = Fint(:,i) + FijLJ
       Fint(:,j) = Fint(:,j) - FijLJ

     end if

   end do

! ##

   Nas = NProcsTemp - MyRankTemp

   do l = Nas , NumDihedralgA, NProcsTemp

     if(vdWSubtDihgA(l)) cycle

     i = DihedIgA(l)
     j = DihedLgA(l)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     InvR2 = 1.d0 / R2

     cf = Charge(i) * Charge(j)

     if(cf /= 0.) then

       R1  = sqrt( R2 )

       fk1 = cf / R1
       fk  = fk1 * InvR2

       FijEL = fk * Rij

       Fint(:,i) = Fint(:,i) + FijEL
       Fint(:,j) = Fint(:,j) - FijEL

     end if

     Sgm  = Rminh14(i) + Rminh14(j)
     Sgm2 = Sgm * Sgm
     Eps  = Eps14(i) * Eps14(j)

     SC2  = Sgm2 * InvR2
     SC6  = SC2 * SC2 * SC2
     SC12 = SC6 * SC6

     fk14 = Eps * 12.d0 * ( SC12 - SC6 ) * InvR2
     FijLJ = fk14 * Rij

     Fint(:,i) = Fint(:,i) + FijLJ
     Fint(:,j) = Fint(:,j) - FijLJ

   end do

   do i = NiniF+1 , NfinF
     Fint(:,i) = Fint(:,i) + ( Frc_Bond(:,i)  + Frc_Angle(:,i) + Frc_UB(:,i) &
     &                       + Frc_Dihed(:,i) + Frc_Impro(:,i) )
   end do

end subroutine Force_gA
! << F monitor ##
