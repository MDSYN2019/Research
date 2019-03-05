! ############################
! ## SUBROUTINE LIST 
! ## -- RB_data 
! ## -- RB_check 
! ############################


!######################################################################
!######################################################################


subroutine RB_data

use Numbers, only : N
use CommonBlocks, only : Qdebug
use FFParameters
use RBparam
use UnitExParam, only : Avogadro
use AtomParam, only : ResidName

implicit none

integer, parameter :: MaxNumRBType = 100
integer, parameter :: MaxNumST = 100
integer, parameter :: MaxNinRBparam = 30
integer :: i, j, k
character(len=4) :: TmpRName
character(len=6), dimension(MaxNumRBType) :: NameStore
integer, dimension(MaxNumRBType) :: NumRBAtomT
logical :: NewFlag
character(len=80) :: String, String1
character(len= 6) :: NameT
character(len= 4) :: Dummy
real(8), dimension(:), allocatable :: MassTmp
character(len=6), dimension(:), allocatable :: NameUnitType
integer :: Iunit, MyType

integer :: NumST
character(len=6), dimension(MaxNumST) :: UName
integer, dimension(MaxNumST) :: NinRBParam
real(8), dimension(MaxNinRBparam,MaxNumST) :: MassST
real(8), dimension(3,MaxNinRBparam,MaxNumST) :: R_onMolST

logical, dimension(:), allocatable :: QLintmp, QSintmp

   NumRBType = 1

   if(Qdebug) then
   do i = 1 , N
     print *, i, AtomUnitName(i), AtomUnitNum(i)
   end do

   print * , 'NumResidueParam=',NumResidueParam
   do j = 1 , NumResidueParam
     write(*,*) ResiNameParam(j)
     do k = 1 , NumUnit_inResi(j)
       write(*,*) UNITNameParam(k,j), NumAtom_inUNITParam(k,j)
     end do
   end do
   end if

   do i = 1 , N

     if( i==1 ) then

       NameStore(NumRBType) = AtomUnitName(i)
       TmpRName             = ResidName(i)

ident: do j = 1, NumResidueParam

         if( TmpRName == ResiNameParam(j) ) then

           do k = 1 , NumUnit_inResi(j)

             if(UNITNameParam(k,j)==NameStore(NumRBType)) then

               NumRBAtomT(NumRBType) = NumAtom_inUNITParam(k,j)
               exit ident

             end if

           end do

         end if

       end do ident

     else

       NewFlag = .True.

       do j = 1 , NumRBType

         if(AtomUnitName(i)==NameStore(j)) NewFlag = .False.

       end do

       if(NewFlag) then

         NumRBType = NumRBType + 1

         if(NumRBType > MaxNumRBType) then
           write(*,*) 'ERROR : NumRBType exceeds its maximum limit!'
           write(*,*) 'the value of the parameter (MaxNumRBType) should be changed !'
           write(*,*) NumRBType
           call Finalize
         end if

         NameStore(NumRBType) = AtomUnitName(i)

         TmpRName = ResidName(i)

ident2:  do j = 1, NumResidueParam

           if( TmpRName == ResiNameParam(j) ) then

             do k = 1 , NumUnit_inResi(j)

               if(UNITNameParam(k,j)==NameStore(NumRBType)) then

                 NumRBAtomT(NumRBType) = NumAtom_inUNITParam(k,j)
                 exit ident2

               end if

             end do

           end if

         end do ident2

       end if

     end if

   end do

   allocate( NameUnitType(NumRBType) )
   allocate( NumRBAtom(NumRBType) )

   do i = 1 , NumRBType

     NameUnitType(i) = NameStore(i)
     NumRBAtom(i) = NumRBAtomT(i)

   end do

   allocate( RBType(NumRB) )

   do i = 1 , N

     if( i == 1 ) then

       Iunit = 1

       do j = 1 , NumRBType

         if(AtomUnitName(i)==NameUnitType(j)) then

           RBType(Iunit) = j
           exit

         end if

       end do

     else if(AtomUnitNum(i)/=AtomUnitNum(i-1)) then

       Iunit = Iunit + 1

       do j = 1 , NumRBType

         if(AtomUnitName(i)==NameUnitType(j)) then

           RBType(Iunit) = j
           exit

         end if

       end do

     end if

   end do

   MaxNumAtomRB = NumRBAtom(1)

   do i = 2 , NumRBType

     if( NumRBAtom(i) > MaxNumAtomRB ) then

       MaxNumAtomRB = NumRBAtom(i)

     end if

   end do

   allocate( R_onMol(3,MaxNumAtomRB,NumRBType) )
   allocate( MassTmp(MaxNumAtomRB) )
   allocate( MassRB(NumRBType) )
   allocate( InvMassRB(NumRBType) )
   allocate( QLintmp(NumRBType) )
   allocate( QLinear(NumRB) )
   allocate( QSintmp(NumRBType) )
   allocate( QSingle(NumRB) )
   allocate( InertiaRB(3,NumRBType) )
   allocate( InvInertiaRB(3,NumRBType) )

   allocate( Rmolec(3,MaxNumAtomRB,NumRB) )
   allocate( Rotation(3,3,NumRB) )

! ## for safety
   R_onMol = 0.d0
   Rmolec  = 0.d0

   open(18,file = './param/RigidBodyModel.prm', status='old')

   NumST = 0

   do

     read(18,'(a80)') String1
     String = adjustl(String1)

     if((String(1:1) == '!').or.(String(1:1) == '#')) cycle
     if(String(1:5) == '<end>') exit

     if(String(1:4) == 'UNIT') then

       NumST = NumST + 1

       if(NumST > MaxNumST) then
         write(*,*) 'ERROR : UName, NinRBparam'
         write(*,*) '        NumST exceeds the default allocation'
         write(*,*) 'edit (subroutine RB_data in RB_pre.f90)'
         write(*,*) 'assign a larger dimension for UName, NinRBparam, MassST, R_onMolST'
         write(*,*) 'MaxNumST should be changed to a larger value.'
         call Finalize
       end if

       read(String,*) Dummy, UName(NumST), NinRBparam(NumST)

       if(Qdebug) print *, NumST, UName(NumST), NinRBparam(NumST)

       if(NinRBparam(NumST) > MaxNinRBparam) then
         write(*,*) 'ERROR : MassST, R_onMolST '
         write(*,*) '        NinRBparam(NumST) exceeds the default allocation'
         write(*,*) 'edit (subroutine RB_data in RB_pre.f90)'
         write(*,*) 'assign a larger dimension for MassST, R_onMolST'
         write(*,*) 'MaxNinRBparam should be changed to a larger value.'
         call Finalize
       end if

       do j = 1 , NinRBparam(NumST)

         read(18,*) k,MassST(j,NumST),R_onMolST(:,j,NumST)

       end do

     end if

   end do

   close(18)

   if(Qdebug) then
   do i = 1 , NumRBType
     print *, NameUnitType(i)
   end do
   end if

   do i = 1 , NumRBType

     NameT = NameUnitType(i)

     QLintmp(i) = .False.

     if(NameT == 'SNGL') then
       QSintmp(i)        = .True.
       MassRB(i)         = 0.d0
       InvMassRB(i)      = 0.d0
       InertiaRB(:,i)    = 0.d0
       InvInertiaRB(:,i) = 0.d0
       cycle
     end if

     QSintmp(i) = .False.

     do j = 1 , NumST

       if(NameT == UName(j)) then

         if(NumRBAtom(i) /= NinRBparam(j)) then

           write(*,*) 'ERROR : defined RB model is not consistent with the model in RB file.'
           write(*,*) NameT,NumRBAtom(i), NinRBparam(j)
           call Finalize

         end if

         do k = 1 , NinRBparam(j)

            MassTmp(k) = MassST(k,j)
            R_onMol(:,k,i) = R_onMolST(:,k,j)

         end do

         exit

       end if

       if( j == NumST ) then

         write(*,*) 'ERROR : missing parameters : rigid body model'
         write(*,*) 'unit name : ',NameT
         call Finalize

       end if

     end do

     MassRB(i) = 0.d0

     do j = 1 , NumRBAtom(i)

       MassTmp(j) = MassTmp(j) * 1.d-3 / Avogadro
       MassRB(i)  = MassRB(i) + MassTmp(j)

     end do

     InvMassRB(i) = 1.d0 / MassRB(i)

     InertiaRB(:,i) = 0.d0

     do j = 1 , NumRBAtom(i)

       InertiaRB(1,i) = InertiaRB(1,i) &
       &              + MassTmp(j) * ( R_onMol(2,j,i)*R_onMol(2,j,i) &
       &                             + R_onMol(3,j,i)*R_onMol(3,j,i) )
       InertiaRB(2,i) = InertiaRB(2,i) &
       &              + MassTmp(j) * ( R_onMol(3,j,i)*R_onMol(3,j,i) &
       &                             + R_onMol(1,j,i)*R_onMol(1,j,i) )
       InertiaRB(3,i) = InertiaRB(3,i) &
       &              + MassTmp(j) * ( R_onMol(1,j,i)*R_onMol(1,j,i) &
       &                             + R_onMol(2,j,i)*R_onMol(2,j,i) )

     end do

     if((InertiaRB(1,i)==0.).and.&
     &  (InertiaRB(2,i)==InertiaRB(3,i))) QLintmp(i) = .True.

     if((InertiaRB(2,i)==0).or.(InertiaRB(3,i)==0)) then
       write(*,*) 'ERROR   : change the defined intramolecular coordinate in '
       write(*,*) '          RigidBodyModel.prm'
       call Finalize
     end if

     if(QLintmp(i)) then
       InvInertiaRB(1,i) = 0.d0
       InvInertiaRB(2,i) = 1.d0 / InertiaRB(2,i)
       InvInertiaRB(3,i) = 1.d0 / InertiaRB(3,i)
     else
       InvInertiaRB(:,i) = 1.d0 / InertiaRB(:,i)
     end if

   end do

   do i = 1 , NumRB

     MyType = RBType(i)

     QSingle(i) = QSintmp(MyType)
     QLinear(i) = QLintmp(MyType)

   end do

   if(Qdebug) then
   do i = 1 , NumRB
     MyType = RBType(i)
     print *, 'Numb =',i,QSingle(i),QLinear(i)
     print *, 'InertiaRB    = ',InertiaRB(:,MyType)
     print *, 'InvInertiaRB = ',InvInertiaRB(:,MyType)
   end do
   end if

   allocate( InitAtomRB(NumRB) )

   j = 0

   do i = 1 , NumRB

     MyType = RBType(i)

     InitAtomRB(i) = j

     j = j + NumRBAtom(MyType)

   end do

   call RB_check(NameUnitType)

end subroutine RB_data


!######################################################################
!######################################################################


subroutine RB_check(NameUnitType)

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonBlocks, only : Qdebug
use RBparam

implicit none

character(len=6), dimension(*) :: NameUnitType
integer :: i, j1, j2, ii, Nt

   allocate( NumRBinMol(NumSpec) )

   ii = 0
   j1 = 0

   do i = 1 , NumSpec

     ii = ii + NumMol(i)*NumAtm(i)
     j2 = AtomUnitNum(ii)
     NumRBinMol(i) = ( j2 - j1 ) / NumMol(i)
     j1 = j2

   end do

   Nt = 0

   do i = 1 , NumSpec

     Nt = Nt + NumRBinMol(i) * NumMol(i)

   end do

   if(Nt /= NumRB) then

     write(*,*) 'ERROR : the number of rigid-bodies'
     call Finalize

   end if


   if(Qdebug) then
   do i = 1 , NumRB
     print *, NameUnitType(RBType(i)),QLinear(RBType(i))
   end do
   end if

end subroutine RB_check
