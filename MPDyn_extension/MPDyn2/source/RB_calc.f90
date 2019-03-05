! ############################
! ## SUBROUTINE LIST 
! ## -- IntraMolVec 
! ############################


!######################################################################
!######################################################################


! **********************************************************
! ** subroutine for calculating intramolecular coordinate **
! **********************************************************

subroutine IntraMolVec

use Configuration, only : R
use RBparam, only : NumRB, QSingle, Quaternion, Rotation, Rmolec, &
&   R_RB, RBType, NumRBAtom, R_onMol

implicit NONE

integer  :: i, j, k, l, Nc, MyType
real(8), dimension(10) :: qt2

   do i = 1 , NumRB

     if(QSingle(i)) cycle

     qt2( 1) = Quaternion(1,i) * Quaternion(1,i)
     qt2( 2) = Quaternion(2,i) * Quaternion(2,i)
     qt2( 3) = Quaternion(3,i) * Quaternion(3,i)
     qt2( 4) = Quaternion(4,i) * Quaternion(4,i)

     qt2( 5) = Quaternion(1,i) * Quaternion(2,i)
     qt2( 6) = Quaternion(1,i) * Quaternion(3,i)
     qt2( 7) = Quaternion(1,i) * Quaternion(4,i)
     qt2( 8) = Quaternion(2,i) * Quaternion(3,i)
     qt2( 9) = Quaternion(2,i) * Quaternion(4,i)
     qt2(10) = Quaternion(3,i) * Quaternion(4,i)

     Rotation(1,1,i) = qt2(1) + qt2(2) - qt2(3) - qt2(4)
     Rotation(2,2,i) = qt2(1) - qt2(2) + qt2(3) - qt2(4)
     Rotation(3,3,i) = qt2(1) - qt2(2) - qt2(3) + qt2(4)

     Rotation(1,2,i) = 2.d0 * ( qt2( 8) + qt2( 7) )
     Rotation(2,1,i) = 2.d0 * ( qt2( 8) - qt2( 7) )
     Rotation(1,3,i) = 2.d0 * ( qt2( 9) - qt2( 6) )
     Rotation(3,1,i) = 2.d0 * ( qt2( 9) + qt2( 6) )
     Rotation(2,3,i) = 2.d0 * ( qt2(10) + qt2( 5) )
     Rotation(3,2,i) = 2.d0 * ( qt2(10) - qt2( 5) )

   end do

! body fixed coordinates --->  space fixed coordinates

   k = 0

   do i = 1 , NumRB

     if(QSingle(i)) then

       k = k + 1
       Rmolec(:,1,i) = 0.d0
       R(:,k) = R_RB(:,i)

     else

       MyType = RBType(i)
       Nc = NumRBAtom(MyType)

       do j = 1 , Nc

         Rmolec(1,j,i) = Rotation(1,1,i) * R_onMol(1,j,MyType) &
         &             + Rotation(2,1,i) * R_onMol(2,j,MyType) & 
         &             + Rotation(3,1,i) * R_onMol(3,j,MyType)
         Rmolec(2,j,i) = Rotation(1,2,i) * R_onMol(1,j,MyType) &
         &             + Rotation(2,2,i) * R_onMol(2,j,MyType) & 
         &             + Rotation(3,2,i) * R_onMol(3,j,MyType)
         Rmolec(3,j,i) = Rotation(1,3,i) * R_onMol(1,j,MyType) &
         &             + Rotation(2,3,i) * R_onMol(2,j,MyType) & 
         &             + Rotation(3,3,i) * R_onMol(3,j,MyType)

         l = k + j
         R(:,l) = R_RB(:,i) + Rmolec(:,j,i)

       end do

       k = k + Nc

     end if

   end do

end subroutine IntraMolVec
