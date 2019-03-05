! ############################
! ## SUBROUTINE LIST 
! ## -- MakeSHAKEList 
! ############################


!######################################################################
!######################################################################


! ********************************************
! ** Making a list for the constraint bonds **
! ********************************************

subroutine MakeSHAKEList

use FFParameters
use SHAKEparam
use BondedParam, only : NumBond, BondI, BondJ, kBond, rBond
use AtomParam, only : AtomName, ResidName

implicit none

integer :: i, j, k, NBD
integer :: ii, jj, kk
character(len=4) :: ANameI,ANameJ,RName
character(len=4) :: ANameII,ANameJJ
logical, dimension(:), allocatable :: FlagCoupleBond
integer :: Ngroup
integer, dimension(:), allocatable :: NCouple,NCoAtom
integer, dimension(:,:), allocatable :: NCpair
integer, dimension(:,:,:), allocatable :: CPair
integer, dimension(:,:), allocatable :: CAtom

integer :: TmpNBond, TotalHcBond
integer, dimension(:), allocatable :: TmpBondI, TmpBondJ
real(8), dimension(:), allocatable :: TmpkBond, TmprBond

integer :: NonHAtom, NumWM
integer, dimension(4) :: HAtom

real(8) :: rHHW
character(len=4) :: WaterName

! integer :: bugc

   NBD=NumBond

   allocate( FlagCoupleBond(NBD) )
   allocate( NCouple(NBD) )
   allocate( NCoAtom(NBD) )
   allocate( NCpair(NBD,4) )
   allocate( CPair(NBD,4,2) )
   allocate( CAtom(NBD,5) )
   allocate( TmpBondI(NBD) )
   allocate( TmpBondJ(NBD) )
   allocate( TmpkBond(NBD) )
   allocate( TmprBond(NBD) )

   FlagCoupleBond = .False.
   Ngroup  = 0
   NCouple = 0
   NCoAtom = 0
   NumWM   = 0

   do k = 1, NumBond

     if(FlagCoupleBond(k)) cycle    ! ## it's been listed already

     i = BondI(k)
     j = BondJ(k)

     ANameI = AtomName(i)
     ANameJ = AtomName(j)

     RName = ResidName(i)

! ## for water (TIP3P or SPCE or SPC model)

     if ( ( RName == 'TIP3' ).or.( RName == 'SPC' )&
     &.or.( RName == 'SPCE' )) then

         FlagCoupleBond(k)              = .True.
         Ngroup                         = Ngroup + 1
         NCouple(Ngroup)                = NCouple(Ngroup) + 1
         NCpair(Ngroup,NCouple(Ngroup)) = k
         NumWM                          = NumWM + 1

         WaterName = RName

         if(ANameI(1:1)=='H') then
           HAtom(NCouple(Ngroup)) = i
           NonHAtom               = j
         else
           NonHAtom               = i
           HAtom(NCouple(Ngroup)) = j
         end if

  ! ## extract the coupling bonds

         do kk = k + 1 , NumBond

           if(FlagCoupleBond(kk)) cycle

           ii = BondI(kk)
           jj = BondJ(kk)

           if(ii == NonHAtom) then
             FlagCoupleBond(kk)             = .True.
             NCouple(Ngroup)                = NCouple(Ngroup) + 1
             NCpair(Ngroup,NCouple(Ngroup)) = kk
             HAtom(NCouple(Ngroup))         = jj
           else if( jj == NonHAtom ) then
             FlagCoupleBond(kk)             = .True.
             NCouple(Ngroup)                = NCouple(Ngroup) + 1
             NCpair(Ngroup,NCouple(Ngroup)) = kk
             HAtom(NCouple(Ngroup))         = ii
           end if

         end do

         if( NCouple(Ngroup) /= 2 ) then
#ifndef BMONI
           write(11,*) 'ERROR : constrant number of a water mol. is ', &
           &            NCouple(Ngroup)
#endif
           write( *,*) 'ERROR : constrant number of a water mol. is ', &
           &            NCouple(Ngroup)
           call Finalize
         end if

         NCouple(Ngroup)  = 3
         NCpair(Ngroup,3) = NCPair(Ngroup,2)
         NCPair(Ngroup,2) = 0     ! take into account the consistency with 'CPair'

  ! ## assign parameters

         CPair(Ngroup,1,1) = 1
         CPair(Ngroup,1,2) = 2

         CPair(Ngroup,2,1) = 2
         CPair(Ngroup,2,2) = 3

         CPair(Ngroup,3,1) = 3
         CPair(Ngroup,3,2) = 1

         CAtom(Ngroup,1)   = NonHAtom
         CAtom(Ngroup,2)   = HAtom(1)
         CAtom(Ngroup,3)   = HAtom(2)

         NCoAtom(Ngroup) = 3

! ## other case ( C-H, N-H, ... etc. )

     else if( ( ANameI(1:1) == 'H' ) .or. ( ANameJ(1:1) == 'H' ) ) then

         FlagCoupleBond(k) = .True.
         Ngroup            = Ngroup + 1
         NCouple(Ngroup)   = NCouple(Ngroup) + 1
         NCpair(Ngroup,NCouple(Ngroup)) = k

         if( ANameI(1:1) == 'H' ) then

           HAtom(NCouple(Ngroup)) = i
           NonHAtom               = j

         else if( ANameJ(1:1) == 'H' ) then

           NonHAtom               = i
           HAtom(NCouple(Ngroup)) = j

         end if

         if( k /= NumBond ) then

           do kk = k + 1, NumBond

             if(FlagCoupleBond(kk)) cycle

             ii      = BondI(kk)
             jj      = BondJ(kk)
             ANameII = AtomName(ii)
             ANameJJ = AtomName(jj)

             if( ( ii == NonHAtom ) .and. ( ANameJJ(1:1) == 'H' ) ) then

               FlagCoupleBond(kk)             = .True.
               NCouple(Ngroup)                = NCouple(Ngroup) + 1
               NCpair(Ngroup,NCouple(Ngroup)) = kk
               HAtom(Ncouple(Ngroup))         = jj

             else if( ( jj == NonHAtom ) .and. ( ANameII(1:1) == 'H' ) ) then

               FlagCoupleBond(kk)             = .True.
               NCouple(Ngroup)                = NCouple(Ngroup) + 1
               NCpair(Ngroup,NCouple(Ngroup)) = kk
               HAtom(Ncouple(Ngroup))         = ii

             end if

           end do

         end if

         CAtom(Ngroup,1)   = NonHAtom

         do kk = 1, NCouple(Ngroup)

           CPair(Ngroup,kk,1) = 1
           CPair(Ngroup,kk,2) = kk + 1
           CAtom(Ngroup,kk+1) = HAtom(kk)

         end do

         NCoAtom(Ngroup) = NCouple(Ngroup) + 1

     end if

   end do

! ---------------------------------------------------------------------

   if( NumWM /= 0 ) then

     if(WaterName == 'TIP3') then

       do j = 1, NumBondParam
! ## in case of TIP3P
         if( (BondPairAtoms(1,j) == 'HT' ) .and. (BondPairAtoms(2,j) == 'HT' ) ) then
           rHHW = rBondParam(j)
         end if
       end do

     else if(WaterName == 'SPCE') then
! ## in case of SPCE
       do j = 1, NumBondParam
         if( (BondPairAtoms(1,j) == 'HTE' ) .and. (BondPairAtoms(2,j) == 'HTE' ) ) then
           rHHW = rBondParam(j)
         end if
       end do

     else if(WaterName == 'SPC') then
! ## in case of SPC
       do j = 1, NumBondParam
         if( (BondPairAtoms(1,j) == 'HTC' ) .and. (BondPairAtoms(2,j) == 'HTC' ) ) then
           rHHW = rBondParam(j)
         end if
       end do

     end if

   end if

   TmpBondI = BondI
   TmpBondJ = BondJ
   TmpkBond = kBond
   TmprBond = rBond
   TmpNBond = NumBond

   deallocate( BondI,BondJ,kBond,rBond )

   NumBond = 0

   do i = 1, TmpNBond

     if(.not.FlagCoupleBond(i)) then

       NumBond = NumBond + 1

     end if

   end do

   MaxHcBond   = 0
   TotalHcBond = 0

   do i = 1, Ngroup

     if( NCouple(i) > 4 ) then

#ifndef BMONI
       write(11,*) 'ERROR :: Too many coupled bonds.'
#endif
       write( *,*) 'ERROR :: Too many coupled bonds.'
       call Finalize

     end if

     if( MaxHcBond < NCouple(i) ) MaxHcBond = NCouple(i)

     TotalHcBond = TotalHcBond + NCouple(i)

   end do

   if( (NumBond + TotalHcBond - NumWM) /= TmpNBond ) then

#ifndef BMONI
     write(11,*) 'ERROR in making list of SHAKE pairs'
#endif
     write( *,*) 'ERROR in making list of SHAKE pairs'
     call Finalize

   end if

! ---------------------------------------------------------------------

   if(NumBond /= 0) then

     allocate( BondI(NumBond) )
     allocate( BondJ(NumBond) )
     allocate( kBond(NumBond) )
     allocate( rBond(NumBond) )

     j = 0
     do i = 1, TmpNBond

       if(.not.FlagCoupleBond(i)) then

         j = j + 1
         BondI(j) = TmpBondI(i)
         BondJ(j) = TmpBondJ(i)
         kBond(j) = TmpkBond(i)
         rBond(j) = TmprBond(i)

       end if

     end do

   end if

! ---------------------------------------------------------------------

   NSHAKEGroup = Ngroup

   if(NSHAKEGroup /= 0) then

     allocate( NCoupleBond(NSHAKEGroup) )
     allocate( NCoupleAtom(NSHAKEGroup) )
     allocate( CouplePair(NSHAKEGroup,MaxHcBond  ,2) )
     allocate( CoupleAtom(NSHAKEGroup,MaxHcBond+1  ) )
     allocate( rSHAKE(NSHAKEGroup,MaxHcBond) )
     allocate( Lagmultip(NSHAKEGroup,MaxHcBond) )

     Lagmultip = 0.d0

     do i = 1, NSHAKEGroup

       NCoupleBond(i) = NCouple(i)
       NCoupleAtom(i) = NCoAtom(i)

       do j = 1, NCoupleBond(i)

         k = NCpair(i,j)
         CouplePair(i,j,1) = CPair(i,j,1)
         CouplePair(i,j,2) = CPair(i,j,2)

         if( k /= 0 ) then

           rSHAKE(i,j) = TmprBond(k) * TmprBond(k)

         else

           rSHAKE(i,j) = rHHW * rHHW

         end if

       end do

       do j = 1, NCoupleAtom(i)

         CoupleAtom(i,j) = CAtom(i,j)

       end do

     end do

   end if

end subroutine MakeSHAKEList
