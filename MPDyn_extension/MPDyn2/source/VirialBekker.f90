! ############################
! ## SUBROUTINE LIST 
! ## -- VirialBekkerRB  
! ## -- VirialBekker  
! ## -- VirialBekkerPI  
! ############################


!######################################################################
!######################################################################


subroutine VirialBekkerRB(Frc,Vir,Gk)

use Numbers, only : N
use RBparam, only : NumRB, R_RB, QSingle, RBType, NumRBAtom
use CellParam, only : CellShft

implicit none

integer :: i, j, k, ii, MyType, Nc
real(8), dimension(3,N) :: Frc
real(8), dimension(3,-13:13) :: Gk
real(8), dimension(3,3) :: Vir
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8) :: Vxx, Vxy, Vxz, Vyx, Vyy, Vyz, Vzx, Vzy, Vzz

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyx = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzx = 0.d0
   Vzy = 0.d0
   Vzz = 0.d0

   k = 0

   do i = 1, NumRB

     Rx = R_RB(1,i)
     Ry = R_RB(2,i)
     Rz = R_RB(3,i)

     if(QSingle(i)) then

       k = k + 1

       Fx = Frc(1,k)
       Fy = Frc(2,k)
       Fz = Frc(3,k)

       Vxx = Vxx + Fx * Rx
       Vxy = Vxy + Fx * Ry
       Vxz = Vxz + Fx * Rz
       Vyx = Vyx + Fy * Rx
       Vyy = Vyy + Fy * Ry
       Vyz = Vyz + Fy * Rz
       Vzx = Vzx + Fz * Rx
       Vzy = Vzy + Fz * Ry
       Vzz = Vzz + Fz * Rz

     else

       MyType = RBType(i)
       Nc = NumRBAtom(MyType)

       Fx = 0.d0
       Fy = 0.d0
       Fz = 0.d0

       do j = 1, Nc
         ii = j + k
         Fx = Fx + Frc(1,ii)
         Fy = Fy + Frc(2,ii)
         Fz = Fz + Frc(3,ii)
       end do

       Vxx = Vxx + Fx * Rx
       Vxy = Vxy + Fx * Ry
       Vxz = Vxz + Fx * Rz
       Vyx = Vyx + Fy * Rx
       Vyy = Vyy + Fy * Ry
       Vyz = Vyz + Fy * Rz
       Vzx = Vzx + Fz * Rx
       Vzy = Vzy + Fz * Ry
       Vzz = Vzz + Fz * Rz
       k = k + Nc

     end if

   end do

   do k = -13, 13

     if(k==0) cycle

     Fx = Gk(1,k)
     Fy = Gk(2,k)
     Fz = Gk(3,k)
     Rx = CellShft(1,k)
     Ry = CellShft(2,k)
     Rz = CellShft(3,k)

     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyx = Vyx + Fy * Rx
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzx = Vzx + Fz * Rx
     Vzy = Vzy + Fz * Ry
     Vzz = Vzz + Fz * Rz

   end do

   Vir(1,1) = Vir(1,1) + Vxx
   Vir(1,2) = Vir(1,2) + Vxy
   Vir(1,3) = Vir(1,3) + Vxz
   Vir(2,1) = Vir(2,1) + Vyx
   Vir(2,2) = Vir(2,2) + Vyy
   Vir(2,3) = Vir(2,3) + Vyz
   Vir(3,1) = Vir(3,1) + Vzx
   Vir(3,2) = Vir(3,2) + Vzy
   Vir(3,3) = Vir(3,3) + Vzz

end subroutine VirialBekkerRB


!######################################################################
!######################################################################


subroutine VirialBekker(Frc,Vir,Gk)

use Numbers, only : N
use Configuration, only : R
use CellParam, only : CellShft

implicit none

integer :: i, k
real(8), dimension(3,N) :: Frc
real(8), dimension(3,-13:13) :: Gk
real(8), dimension(3,3) :: Vir
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzz = 0.d0

   do i = 1, N

     Fx = Frc(1,i)
     Fy = Frc(2,i)
     Fz = Frc(3,i)
     Rx = R(1,i)
     Ry = R(2,i)
     Rz = R(3,i)

     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzz = Vzz + Fz * Rz

   end do

   do k = -13, 13

     if(k==0) cycle
     Fx = Gk(1,k)
     Fy = Gk(2,k)
     Fz = Gk(3,k)
     Rx = CellShft(1,k)
     Ry = CellShft(2,k)
     Rz = CellShft(3,k)
     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzz = Vzz + Fz * Rz

   end do

   Vir(1,1) = Vir(1,1) + Vxx
   Vir(1,2) = Vir(1,2) + Vxy
   Vir(1,3) = Vir(1,3) + Vxz
   Vir(2,2) = Vir(2,2) + Vyy
   Vir(2,3) = Vir(2,3) + Vyz
   Vir(3,3) = Vir(3,3) + Vzz
   Vir(2,1) = Vir(1,2)
   Vir(3,1) = Vir(1,3)
   Vir(3,2) = Vir(2,3)

end subroutine VirialBekker


!######################################################################
!######################################################################


subroutine VirialBekkerPI(Frc,Vir,Gk)

use Numbers, only : N
use CommonPI, only : Rcentroid
use CellParam, only : CellShft

implicit none

integer :: i, k
real(8), dimension(3,N) :: Frc
real(8), dimension(3,-13:13) :: Gk
real(8), dimension(3,3) :: Vir
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8) :: Vxx, Vxy, Vxz, Vyx, Vyy, Vyz, Vzx, Vzy, Vzz

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyx = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzx = 0.d0
   Vzy = 0.d0
   Vzz = 0.d0

   do i = 1, N

     Fx = Frc(1,i)
     Fy = Frc(2,i)
     Fz = Frc(3,i)
     Rx = Rcentroid(1,i)
     Ry = Rcentroid(2,i)
     Rz = Rcentroid(3,i)

     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyx = Vyx + Fy * Rx
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzx = Vzx + Fz * Rx
     Vzy = Vzy + Fz * Ry
     Vzz = Vzz + Fz * Rz

   end do

   do k = -13, 13

     if(k==0) cycle

     Fx = Gk(1,k)
     Fy = Gk(2,k)
     Fz = Gk(3,k)
     Rx = CellShft(1,k)
     Ry = CellShft(2,k)
     Rz = CellShft(3,k)

     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyx = Vyx + Fy * Rx
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzx = Vzx + Fz * Rx
     Vzy = Vzy + Fz * Ry
     Vzz = Vzz + Fz * Rz

   end do

   Vir(1,1) = Vir(1,1) + Vxx
   Vir(1,2) = Vir(1,2) + Vxy
   Vir(1,3) = Vir(1,3) + Vxz
   Vir(2,1) = Vir(2,1) + Vyx
   Vir(2,2) = Vir(2,2) + Vyy
   Vir(2,3) = Vir(2,3) + Vyz
   Vir(3,1) = Vir(3,1) + Vzx
   Vir(3,2) = Vir(3,2) + Vzy
   Vir(3,3) = Vir(3,3) + Vzz

end subroutine VirialBekkerPI
