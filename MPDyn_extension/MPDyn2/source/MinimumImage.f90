! ############################
! ## SUBROUTINE LIST 
! ## -- MImageF 
! ## -- MImageO 
! ############################


!######################################################################
!######################################################################


subroutine MImageF(Sij,Rij,R2,k)

use CellParam, only : CellShft

implicit none

integer :: k
real(8), dimension(3) :: Sij, Rij
integer, dimension(3) :: Nij
real(8) :: R2

     Nij(:) = 0
     if(Sij(1) >  0.5d0) then
       Nij(1) = - 9
     else if(Sij(1) < -0.5d0) then
       Nij(1) =   9
     end if
     if(Sij(2) >  0.5d0) then
       Nij(2) = - 3
     else if(Sij(2) < -0.5d0) then
       Nij(2) =   3
     end if
     if(Sij(3) >  0.5d0) then
       Nij(3) = - 1
     else if(Sij(3) < -0.5d0) then
       Nij(3) =   1
     end if
     k = Nij(1) + Nij(2) + Nij(3)
     Rij(:) = Rij(:) + CellShft(:,k)
     R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)

end subroutine MImageF


!######################################################################
!######################################################################


subroutine MImageO(Rij,L,InvL,R2,Nij)

implicit none

real(8), dimension(3) :: Sij, Rij, L, InvL
integer, dimension(3) :: Nij
real(8) :: R2

     Sij(:) = InvL(:) * Rij(:)
#ifdef PCC
     Nij(:) = 0
     if(Sij(1) >  0.5d0) Nij(1) = - 1
     if(Sij(1) < -0.5d0) Nij(1) =   1
     if(Sij(2) >  0.5d0) Nij(2) = - 1
     if(Sij(2) < -0.5d0) Nij(2) =   1
     if(Sij(3) >  0.5d0) Nij(3) = - 1
     if(Sij(3) < -0.5d0) Nij(3) =   1
     Rij(:) = Rij(:) + L(:) * Nij(:)
     R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
     Nij = - nint( Sij )
     Rij = Rij + L * Nij
     R2  = dot_product( Rij, Rij )
#endif

end subroutine MImageO
