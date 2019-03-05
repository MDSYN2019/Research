! ############################
! ## SUBROUTINE LIST 
! ## -- CheckCellShape 
! ## -- Rot_FixedPoint 
! ############################


!######################################################################
!######################################################################


! *******************************************************
! ** Calculation of the Virial term from guiding force **
! *******************************************************

subroutine CheckCellShape

use CommonBlocks, only : QMaster, QCellList, ForceField
use CGdata, only: Rcut_MAX
use UnitExParam, only : pi
use CellParam, only : H
use CutoffParam, only : Rbook2, Rskin

implicit none

real(8), dimension(3) :: a, b, c
real(8), dimension(3) :: ab, bc, ca
real(8) :: Xab, Xbc, Xca
real(8) :: dis_a, dis_b, dis_c
real(8) :: Lna, Lnb, Lnc
real(8) :: Rcrit

   if(QCellList) call CheckCellList

   if(ForceField(1:2) == 'CG') then
     Rcrit = 2.d0 * (Rcut_MAX + Rskin)
   else
     Rcrit = 2.d0 * sqrt(Rbook2)
   end if

   a = H(:,1)
   b = H(:,2)
   c = H(:,3)

   call VecProduct( a , b , ab )
   call VecProduct( b , c , bc )
   call VecProduct( c , a , ca )

   dis_a = dot_product( bc, a )
   dis_b = dot_product( ca, b )
   dis_c = dot_product( ab, c )

   if(( dis_a < Rcrit ).or.( dis_b < Rcrit ).or.( dis_c < Rcrit )) then

! ----------------------
     call SaveParam
! ----------------------

     if(QMaster) then

       write( *,*) 'ERROR : Too large cell strain'
       write( *,'(a)') '-----------------------------------------------------'
       write( *,'(a)') '  The MD run was stopped due to large cell strain !  '
       write( *,'(a)') '-----------------------------------------------------'
#ifndef BMONI
       write(11,*) 'ERROR : Too large cell strain'
       write(11,'(a)') '-----------------------------------------------------'
       write(11,'(a)') '  The MD run was stopped due to large cell strain !  '
       write(11,'(a)') '-----------------------------------------------------'
#endif

       Lna = sqrt(dot_product(a,a))
       Lnb = sqrt(dot_product(b,b))
       Lnc = sqrt(dot_product(c,c))
       Xab = acos( dot_product(a,b) / Lna / Lnb ) * 180. / pi
       Xbc = acos( dot_product(b,c) / Lnb / Lnc ) * 180. / pi
       Xca = acos( dot_product(c,a) / Lnc / Lna ) * 180. / pi

       write( *,'(a,3f7.2)') ' Cell Length = ',Lna,Lnb,Lnc
       write( *,'(a,3f7.2)') '      Angles = ',Xab,Xbc,Xca
#ifndef BMONI
       write(11,'(a,3f7.2)') ' Cell Length = ',Lna,Lnb,Lnc
       write(11,'(a,3f7.2)') '      Angles = ',Xab,Xbc,Xca
#endif

     end if

     call Finalize

   end if

end subroutine CheckCellShape


!######################################################################
!######################################################################


!****************************************
!* Rotational operation on Fixed points *
!****************************************

subroutine Rot_FixedPoint

use Numbers, only : N
use OptConstraintParam, only : Rrot, RIni
use CellParam, only : H

implicit NONE

integer :: i
real(8) :: snt, cst, snp, csp
real(8), dimension(3) :: a, b, ab
real(8) :: dab
real(8), dimension(3,3) :: Rott, A1, A2, A3, InvA3
real(8), dimension(3,N) :: ScR
real :: cs

   a = H(:,1)
   b = H(:,2)

   ab  = VecProd( a , b )
   dab = dot_product( ab, ab )
   ab  = ab / sqrt( dab )

   cst = ab(3)
   cs  = real(cst)

   if(cs == 1.) then

     Rrot = RIni

   else

     snt = sqrt( 1.d0 - cst * cst )

     csp = ab(1) / snt
     snp = ab(2) / snt

     Rott = 0.d0
     Rott(1,1) =  csp
     Rott(1,2) =  snp
     Rott(2,1) = -snp
     Rott(2,2) =  csp
     Rott(3,3) = 1.d0

     A1 = matmul( Rott, H )

     Rott = 0.d0
     Rott(1,1) =  cst
     Rott(1,3) = -snt
     Rott(2,2) = 1.d0
     Rott(3,1) =  snt
     Rott(3,3) =  cst

     A2 = matmul( Rott, A1 )

     Rott = 0.d0
     Rott(1,1) =  csp
     Rott(1,2) = -snp
     Rott(2,1) =  snp
     Rott(2,2) =  csp
     Rott(3,3) = 1.d0

     A3 = matmul( Rott, A2 )

     call InversMatrix(A3,InvA3)

     do i = 1 , N

       ScR(:,i) = matmul( InvA3 , RIni(:,i) )

     end do

     do i = 1, N

       Rrot(:,i) = matmul( H, ScR(:,i) )

     end do

   end if

! #############

Contains

   function VecProd(x,y) Result(z)

   real(8), dimension(3) :: x, y, z

     z(1) = x(2) * y(3) - y(2) * x(3)
     z(2) = x(3) * y(1) - y(3) * x(1)
     z(3) = x(1) * y(2) - y(1) * x(2)

   end function VecProd

end subroutine Rot_FixedPoint
