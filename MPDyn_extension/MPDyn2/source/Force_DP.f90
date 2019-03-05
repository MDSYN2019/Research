! ############################
! ## SUBROUTINE LIST 
! ## -- ForceDPD_Original 
! ## -- ForceDPD_Original_SC 
! ## -- ForceDPD_Smooth22 
! ## -- ForceDPD_Smooth22_SC 
! ## -- ForceDPD_Smooth23 
! ## -- ForceDPD_Smooth23_SC 
! ## -- ForceDPD_Attractive 
! ## -- ForceDPD_Attractive_SC 
! ## -- Energy_Minimization_DPD 
! ## -- ForceDPD_Minim_Original 
! ## -- ForceDPD_Minim_Smooth12 
! ## -- ForceDPD_Minim_Smooth21 
! ## -- ForceDPD_Minim_Smooth22 
! ## -- ForceDPD_Minim_Smooth23 
! ## -- ForceDPD_Minim_Attractive 
! ## -- Force_Bond_DP 
! ## -- ForceDPD_Original_Lowe
! ## -- ForceDPD_Smooth12_Lowe
! ## -- ForceDPD_Smooth21_Lowe
! ## -- ForceDPD_Smooth22_Lowe
! ## -- ForceDPD_Smooth23_Lowe
! ## -- ForceDPD_Attractive_Lowe 
! ############################


!#####################################################################
!#####################################################################


! ****************************************************
! **  Force calculation between DP particles        **
! **  Original weighting function w(r) = (1-r/rc)   **
! ****************************************************

! ## INPUT  : R, Velt
! ## OUTPUT : FrcDPt, PotDP, Virial

subroutine ForceDPD_Original

use Configuration, only : R
use CommonDPD
use CellListMethod
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Vij, Vi, Fcij, Fdij
#ifdef DPDcheck
real(8), dimension(3) :: Frij
#endif
real(8) :: InvR, R1, R2, Wr, Wd, theta
real(8) :: rv, fc, fd, fr, ff
integer :: Slide
real(8) :: ranf, Shv
external ranf

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDPt = Frc_Bond

   VirialC = Vir_Bond
   VirialD = 0.d0
#ifdef DPDcheck
   VirialR = 0.d0
#endif
! #####################################################################
   if(QSheared) then
! #####################################################################
!
     Shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         Vi    = Velt(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij    = Ri - R(:,j)

             Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
             Rij(2) = Rij(2) - nint( Rij(2) * InvCL(2) ) * CellL(2)
             Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

             R2 = dot_product( Rij, Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.d0 / R1
               Wr   = 1.d0 - R1
               Wd   = Wr * Wr

               Sij = Rij * InvR

               theta = ranf() - 0.5d0

               Vij    = Vi - Velt(:,j)
               rv     = dot_product( Vij, Sij )

               fc =   a(is,TypeNum(j)) * Wr
               fd = - gamma * Wd * rv
               fr =   sigmt * Wr * theta

               ff =   fc + fd + fr

               PotDP = PotDP + 0.5d0 * fc * Wr

               Fcij = fc * Sij
               Fdij = fd * Sij

               VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
               VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
               VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
               VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
               VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
               VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)

               VirialD(1,1) = VirialD(1,1) + Fdij(1) * Rij(1)
               VirialD(1,2) = VirialD(1,2) + Fdij(1) * Rij(2)
               VirialD(1,3) = VirialD(1,3) + Fdij(1) * Rij(3)
               VirialD(2,2) = VirialD(2,2) + Fdij(2) * Rij(2)
               VirialD(2,3) = VirialD(2,3) + Fdij(2) * Rij(3)
               VirialD(3,3) = VirialD(3,3) + Fdij(3) * Rij(3)
#ifdef DPDcheck
               Frij = fr * Sij
               VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
               VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
               VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
               VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
               VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
               VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

               FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij(:)

               FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij(:)

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j > 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij    = Ri - R(:,j)

                 Slide  = - nint( Rij(2) * InvCL(2) )
                 Rij(1) = Rij(1) + Slide * SlideGap
                 Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                 Rij(2) = Rij(2) + Slide                     * CellL(2)
                 Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                 R2     = dot_product( Rij, Rij )

                 if( R2 < 1.d0 ) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1
                   Wr   = 1.d0 - R1
                   Wd   = Wr * Wr

                   Sij = Rij * InvR

                   theta = ranf()-0.5d0

                   Vij = Vi - Velt(:,j)
                   Vij(1) = Vij(1) + Slide * shv
                   rv  = dot_product( Vij, Sij )

                   fc =   a(is,TypeNum(j)) * Wr
                   fd = - gamma * Wd * rv
                   fr =   sigmt * Wr * theta

                   ff =   fc + fd + fr

                   PotDP = PotDP + 0.5d0 * fc * Wr

                   Fcij = fc * Sij
                   Fdij = fd * Sij

                   VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                   VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                   VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                   VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                   VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                   VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)

                   VirialD(1,1) = VirialD(1,1) + Fdij(1) * Rij(1)
                   VirialD(1,2) = VirialD(1,2) + Fdij(1) * Rij(2)
                   VirialD(1,3) = VirialD(1,3) + Fdij(1) * Rij(3)
                   VirialD(2,2) = VirialD(2,2) + Fdij(2) * Rij(2)
                   VirialD(2,3) = VirialD(2,3) + Fdij(2) * Rij(3)
                   VirialD(3,3) = VirialD(3,3) + Fdij(3) * Rij(3)
#ifdef DPDcheck
                   Frij = fr * Sij
                   VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
                   VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
                   VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
                   VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
                   VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
                   VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif
                   FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij(:)

                   FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij(:)

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     VirialD(2,1) = VirialD(1,2)
     VirialD(3,1) = VirialD(1,3)
     VirialD(3,2) = VirialD(2,3)

     VirialR(2,1) = VirialR(1,2)
     VirialR(3,1) = VirialR(1,3)
     VirialR(3,2) = VirialR(2,3)

     Virial = VirialC + VirialD

! #####################################################################
   else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Vi    = Velt(:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1
               Wr   = 1.0d0 - R1
               Wd   = Wr * Wr

               Sij = Rij * InvR

               theta = ranf() - 0.5d0

               Vij = Vi - Velt(:,j)
               rv  = dot_product( Vij, Sij )

               fc =   a(Is,TypeNum(j)) * Wr
               fd = - gamma * Wd * rv
               fr =   sigmt * Wr * theta

               ff =   fc + fd + fr

               PotDP = PotDP + 0.5d0 * fc * Wr

               Fcij = fc * Sij

               VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
               VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
               VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
               VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
               VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
               VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)

#ifdef DPDcheck
               Fdij = fd * Sij
               Frij = fr * Sij

               VirialD(1,1) = VirialD(1,1) + Fdij(1) * Rij(1)
               VirialD(1,2) = VirialD(1,2) + Fdij(1) * Rij(2)
               VirialD(1,3) = VirialD(1,3) + Fdij(1) * Rij(3)
               VirialD(2,2) = VirialD(2,2) + Fdij(2) * Rij(2)
               VirialD(2,3) = VirialD(2,3) + Fdij(2) * Rij(3)
               VirialD(3,3) = VirialD(3,3) + Fdij(3) * Rij(3)

               VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
               VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
               VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
               VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
               VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
               VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

               FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij(:)

               FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij(:)

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1
                   Wr   = 1.d0 - R1
                   Wd   = Wr * Wr

                   Sij  = Rij * InvR

                   theta = ranf() - 0.5d0

                   Vij = Vi - Velt(:,j)
                   rv  = dot_product( Vij, Sij )

                   fc =   a(is,TypeNum(j)) * Wr
                   fd = - gamma * Wd * rv
                   fr =   sigmt * Wr * theta

                   ff =   fc + fd + fr

                   PotDP = PotDP + 0.5d0 * fc * Wr

                   Fcij = fc * Sij

                   VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                   VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                   VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                   VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                   VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                   VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)

#ifdef DPDcheck
                   Fdij = fd * Sij
                   Frij = fr * Sij

                   VirialD(1,1) = VirialD(1,1) + Fdij(1) * Rij(1)
                   VirialD(1,2) = VirialD(1,2) + Fdij(1) * Rij(2)
                   VirialD(1,3) = VirialD(1,3) + Fdij(1) * Rij(3)
                   VirialD(2,2) = VirialD(2,2) + Fdij(2) * Rij(2)
                   VirialD(2,3) = VirialD(2,3) + Fdij(2) * Rij(3)
                   VirialD(3,3) = VirialD(3,3) + Fdij(3) * Rij(3)

                   VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
                   VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
                   VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
                   VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
                   VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
                   VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif
                   FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij(:)

                   FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij(:)

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     VirialD(2,1) = VirialD(1,2)
     VirialD(3,1) = VirialD(1,3)
     VirialD(3,2) = VirialD(2,3)

     VirialR(2,1) = VirialR(1,2)
     VirialR(3,1) = VirialR(1,3)
     VirialR(3,2) = VirialR(2,3)

     Virial = VirialC
!
! #####################################################################
    end if
! #####################################################################

end subroutine ForceDPD_Original


!#####################################################################
!#####################################################################


subroutine ForceDPD_Original_SC

use Configuration, only : R
use CommonDPD
use CellListMethod
use BookParam, only : Npair, ListIJ
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Fcij
#ifdef DPDcheck
real(8), dimension(3) :: Frij
#endif
real(8) :: InvR, R1, R2, Wr, Wd, theta
real(8) :: fc, fr, ff
integer :: Slide
real(8) :: ranf, shv
external ranf

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDPt = Frc_Bond

   VirialC = Vir_Bond
#ifdef DPDcheck
   VirialR = 0.d0
#endif

   Npair = 0
   SLList = 0.d0

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij    = Ri - R(:,j)

              Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
              Rij(2) = Rij(2) - nint( Rij(2) * InvCL(2) ) * CellL(2)
              Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wr   = 1.d0 - R1
                Wd   = Wr * Wr

                Sij = Rij * InvR

                theta = ranf() - 0.5d0

! >> for dissipative force

                Npair = Npair + 1

                ListIJ(1,Npair) = i
                ListIJ(2,Npair) = j

                dRList(:,Npair) = Sij
                R1List(Npair) =  R1
                pfList(Npair) = -gamma * Wd
                SLList(Npair) =  0.d0

! << for dissipative force

                fc =   a(is,TypeNum(j)) * Wr
                fr =   sigmt * Wr * theta

                ff =   fc + fr

                PotDP = PotDP + 0.5d0 * fc * Wr

                Fcij = fc * Sij

                VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)

#ifdef DPDcheck
                Frij = fr * Sij
                VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
                VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
                VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
                VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
                VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
                VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

                FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij(:)

                FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij(:)

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide                     * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1
                    Wr   = 1.d0 - R1
                    Wd   = Wr * Wr

                    Sij = Rij * InvR

                    theta = ranf()-0.5d0

! >> for dissipative force

                    Npair = Npair + 1

                    ListIJ(1,Npair) = i
                    ListIJ(2,Npair) = j

                    dRList(:,Npair) = Sij
                    R1List(Npair) =  R1
                    pfList(Npair) = -gamma * Wd
                    SLList(Npair) =  Slide * shv

! << for dissipative force

                    fc =   a(is,TypeNum(j)) * Wr
                    fr =   sigmt * Wr * theta

                    ff =   fc + fr

                    PotDP = PotDP + 0.5d0 * fc * Wr

                    Fcij = fc * Sij

                    VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                    VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                    VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                    VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                    VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                    VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)

#ifdef DPDcheck
                    Frij = fr * Sij
                    VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
                    VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
                    VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
                    VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
                    VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
                    VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif
                    FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij(:)

                    FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij(:)

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

      VirialC(2,1) = VirialC(1,2)
      VirialC(3,1) = VirialC(1,3)
      VirialC(3,2) = VirialC(2,3)

      VirialR(2,1) = VirialR(1,2)
      VirialR(3,1) = VirialR(1,3)
      VirialR(3,2) = VirialR(2,3)

      Virial = VirialC

! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1
               Wr   = 1.0d0 - R1
               Wd   = Wr * Wr

               Sij = Rij * InvR

               theta = ranf() - 0.5d0

! >> for dissipative force

               Npair = Npair + 1

               ListIJ(1,Npair) = i
               ListIJ(2,Npair) = j

               dRList(:,Npair) = Sij
               R1List(Npair) =  R1
               pfList(Npair) = -gamma * Wd

! << for dissipative force

               fc =   a(Is,TypeNum(j)) * Wr
               fr =   sigmt * Wr * theta

               ff =   fc + fr

               PotDP = PotDP + 0.5d0 * fc * Wr

               Fcij = ff * Sij

               VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
               VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
               VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
               VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
               VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
               VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)
#ifdef DPDcheck
               Frij = fr * Sij
               VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
               VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
               VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
               VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
               VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
               VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

               FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij(:)

               FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij(:)

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1
                   Wr   = 1.d0 - R1
                   Wd   = Wr * Wr

                   Sij  = Rij * InvR

                   theta = ranf() - 0.5d0

! >> for dissipative force

                   Npair = Npair + 1

                   ListIJ(1,Npair) = i
                   ListIJ(2,Npair) = j

                   dRList(:,Npair) = Sij
                   R1List(Npair) =  R1
                   pfList(Npair) = -gamma * Wd

! << for dissipative force

                   fc =   a(is,TypeNum(j)) * Wr
                   fr =   sigmt * Wr * theta

                   ff =   fc + fr

                   PotDP = PotDP + 0.5d0 * fc * Wr

                   Fcij = ff * Sij

                   VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                   VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                   VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                   VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                   VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                   VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)
#ifdef DPDcheck
                   Frij = fr * Sij
                   VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
                   VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
                   VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
                   VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
                   VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
                   VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

                   FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij(:)

                   FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij(:)

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     VirialR(2,1) = VirialR(1,2)
     VirialR(3,1) = VirialR(1,3)
     VirialR(3,2) = VirialR(2,3)

     Virial = VirialC
!
! #####################################################################
    end if
! #####################################################################

    if(Npair > Ndm) then
      write(*,*) 'ERROR : Npair exceeds its limit number'
      write(*,*) 'Npair = ', Npair, ',  Ndm = ',Ndm
      call Finalize
    end if

end subroutine ForceDPD_Original_SC


!#####################################################################
!#####################################################################


subroutine ForceDPD_Smooth22

use Configuration, only : R
use CommonDPD
use CellListMethod
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Vij, Vi, Fij
real(8) :: InvR, R1, R2, Wr, Wd, theta
real(8) :: rv, fc, fd, fr, ff, Slide
real(8) :: ranf, shv
external ranf
! >>
real(8) :: ccc, R5, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDPt = Frc_Bond

! >>
   ccc = 8.d0 / 15.d0
! <<

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         Vi    = Velt(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij    = Ri - R(:,j)

              Slide  = -nint( Rij(2) * InvCL(2) )
              Rij(1) = Rij(1) + Slide * SlideGap

              Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
              Rij(2) = Rij(2) + Slide                     * CellL(2)
              Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wr   = 1.d0 - R1
                Wd   = Wr * Wr

                Sij = Rij * InvR

! >>
                Wc   = 1.d0 - R2

                theta = ranf() - 0.5d0

                Vij    = Vi - Velt(:,j)
                Vij(1) = Vij(1) + Slide * shv
                rv     = dot_product( Vij, Sij )

! >>
                aaa =   a(is,TypeNum(j)) 
                fc  =   aaa * Wc * Wc
                fd  = - gamma * Wd * rv
                fr  =   sigmt * Wr * theta

                ff =   fc + fd + fr

! >>
                R5 = - R2 * R2 / 5.d0
                R3 =   R2 * 2.d0 / 3.d0

                PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                Fij = ff * Sij

                Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                FrcDPt(:,i) = FrcDPt(:,i) + Fij

                FrcDPt(:,j) = FrcDPt(:,j) - Fij

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide                     * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1
                    Wr   = 1.d0 - R1
                    Wd   = Wr * Wr

                    Sij = Rij * InvR

! >>
                    Wc   = 1.d0 - R2

                    theta = ranf()-0.5d0

                    Vij = Vi - Velt(:,j)
                    Vij(1) = Vij(1) + Slide * shv
                    rv  = dot_product( Vij, Sij )

! >>
                    aaa =   a(is,TypeNum(j)) 
                    fc  =   aaa * Wc * Wc
                    fd = - gamma * Wd * rv
                    fr =   sigmt * Wr * theta

                    ff =   fc + fd + fr

! >>
                    R5 = - R2 * R2 / 5.d0
                    R3 =   R2 * 2.d0 / 3.d0

                    PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                    Fij = ff * Sij

                    Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                    Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                    Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                    Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                    Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                    Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                    Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                    Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                    Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                    FrcDPt(:,i) = FrcDPt(:,i) + Fij

                    FrcDPt(:,j) = FrcDPt(:,j) - Fij

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Vi    = Velt(:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1
               Wr   = 1.0d0 - R1
               Wd   = Wr * Wr

               Sij = Rij * InvR
! >>
               Wc   = 1.d0 - R2

               theta = ranf() - 0.5d0

               Vij = Vi - Velt(:,j)
               rv  = dot_product( Vij, Sij )

! >>
               aaa =   a(is,TypeNum(j)) 
               fc  =   aaa * Wc * Wc
               fd = - gamma * Wd * rv
               fr =   sigmt * Wr * theta

               ff =   fc + fd + fr

! >>
               R5 = - R2 * R2 / 5.d0
               R3 =   R2 * 2.d0 / 3.d0

               PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

               Fij = ff * Sij

               Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
               Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
               Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
               Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
               Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
               Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

               FrcDPt(:,i) = FrcDPt(:,i) + Fij

               FrcDPt(:,j) = FrcDPt(:,j) - Fij

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1
                   Wr   = 1.d0 - R1
                   Wd   = Wr * Wr

                   Sij  = Rij * InvR

! >>
                   Wc   = 1.d0 - R2
                   theta = ranf() - 0.5d0

                   Vij = Vi - Velt(:,j)
                   rv  = dot_product( Vij, Sij )

! >>
                   aaa =   a(is,TypeNum(j)) 
                   fc  =   aaa * Wc * Wc
                   fd = - gamma * Wd * rv
                   fr =   sigmt * Wr * theta

                   ff =   fc + fd + fr

! >>
                   R5 = - R2 * R2 / 5.d0
                   R3 =   R2 * 2.d0 / 3.d0

                   PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                   Fij = ff * Sij

                   Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                   Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                   Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                   Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                   Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                   Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                   FrcDPt(:,i) = FrcDPt(:,i) + Fij

                   FrcDPt(:,j) = FrcDPt(:,j) - Fij

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     Virial(2,1) = Virial(1,2)
     Virial(3,1) = Virial(1,3)
     Virial(3,2) = Virial(2,3)
!
! #####################################################################
    end if
! #####################################################################

end subroutine ForceDPD_Smooth22


!#####################################################################
!#####################################################################


subroutine ForceDPD_Smooth22_SC

use Configuration, only : R
use CommonDPD
use CellListMethod
use BookParam, only : Npair, ListIJ
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Fcij
#ifdef DPDcheck
real(8), dimension(3) :: Frij
#endif
real(8) :: InvR, R1, R2, Wr, Wd, theta
real(8) :: fc, fr, ff, Slide
real(8) :: ranf, shv
external ranf
! >>
real(8) :: ccc, R5, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   FrcDPt = Frc_Bond

   VirialC = Vir_Bond
#ifdef DPDcheck
   VirialR = 0.d0
#endif
! >>
   ccc = 8.d0 / 15.d0
! <<

   Npair = 0
   SLList = 0.d0

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij = Ri - R(:,j)
              Rij = Rij - nint( Rij * InvCL ) * CellL

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wr   = 1.d0 - R1
                Wd   = Wr * Wr

                Sij = Rij * InvR
! >>
                Wc   = 1.d0 - R2

                theta = ranf() - 0.5d0

! >> for dissipative force

                Npair = Npair + 1

                ListIJ(1,Npair) = i
                ListIJ(2,Npair) = j

                dRList(:,Npair) = Sij
                R1List(Npair) =  R1
                pfList(Npair) = -gamma * Wd

! << for dissipative force

! >>
                aaa =   a(is,TypeNum(j)) 
                fc  =   aaa * Wc * Wc
                fr  =   sigmt * Wr * theta

                ff =   fc + fr

! >>
                R5 = - R2 * R2 / 5.d0
                R3 =   R2 * 2.d0 / 3.d0

                PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                Fcij = ff * Sij

                VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)
#ifdef DPDcheck
                Frij = fr * Sij
                VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
                VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
                VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
                VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
                VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
                VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif
                FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij

                FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1
                    Wr   = 1.d0 - R1
                    Wd   = Wr * Wr

                    Sij = Rij * InvR
! >>
                    Wc   = 1.d0 - R2

                    theta = ranf()-0.5d0

! >> for dissipative force

                    Npair = Npair + 1

                    ListIJ(1,Npair) = i
                    ListIJ(2,Npair) = j

                    dRList(:,Npair) = Sij
                    R1List(Npair) =  R1
                    pfList(Npair) = -gamma * Wd
                    SLList(Npair) =  Slide * shv

! << for dissipative force

! >>
                    aaa =   a(is,TypeNum(j)) 
                    fc  =   aaa * Wc * Wc
                    fr  =   sigmt * Wr * theta

                    ff =   fc + fr

! >>
                    R5 = - R2 * R2 / 5.d0
                    R3 =   R2 * 2.d0 / 3.d0

                    PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                    Fcij = fc * Sij

                    VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                    VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                    VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                    VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                    VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                    VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)

#ifdef DPDcheck
                   Frij = fr * Sij
                   VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
                   VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
                   VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
                   VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
                   VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
                   VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

                    FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij

                    FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     Virial = VirialC
! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1
               Wr   = 1.0d0 - R1
               Wd   = Wr * Wr

               Sij = Rij * InvR
! >>
               Wc   = 1.d0 - R2

               theta = ranf() - 0.5d0

! >> for dissipative force

               Npair = Npair + 1

               ListIJ(1,Npair) = i
               ListIJ(2,Npair) = j

               dRList(:,Npair) = Sij
               R1List(Npair) =  R1
               pfList(Npair) = -gamma * Wd

! << for dissipative force

! >>
               aaa =   a(is,TypeNum(j)) 
               fc  =   aaa * Wc * Wc
               fr  =   sigmt * Wr * theta

               ff =   fc + fr

! >>
               R5 = - R2 * R2 / 5.d0
               R3 =   R2 * 2.d0 / 3.d0

               PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

               Fcij = fc * Sij

               VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
               VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
               VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
               VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
               VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
               VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)
#ifdef DPDcheck
               Frij = fr * Sij
               VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
               VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
               VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
               VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
               VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
               VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

               FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij

               FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1
                   Wr   = 1.d0 - R1
                   Wd   = Wr * Wr

                   Sij  = Rij * InvR
! >>
                   Wc   = 1.d0 - R2

                   theta = ranf() - 0.5d0

! >> for dissipative force

                   Npair = Npair + 1

                   ListIJ(1,Npair) = i
                   ListIJ(2,Npair) = j

                   dRList(:,Npair) = Sij
                   R1List(Npair) =  R1
                   pfList(Npair) = -gamma * Wd

! << for dissipative force

! >>
                   aaa =   a(is,TypeNum(j)) 
                   fc  =   aaa * Wc * Wc
                   fr  =   sigmt * Wr * theta

                   ff =   fc + fr

! >>
                   R5 = - R2 * R2 / 5.d0
                   R3 =   R2 * 2.d0 / 3.d0

                   PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                   Fcij = fc * Sij

                   VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                   VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                   VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                   VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                   VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                   VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)
#ifdef DPDcheck
                   Frij = fr * Sij
                   VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
                   VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
                   VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
                   VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
                   VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
                   VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

                   FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij

                   FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     Virial(2,1) = Virial(1,2)
     Virial(3,1) = Virial(1,3)
     Virial(3,2) = Virial(2,3)
!
! #####################################################################
    end if
! #####################################################################

    if(Npair > Ndm) then
      write(*,*) 'ERROR : Npair exceeds its limit number'
      write(*,*) 'Npair = ', Npair, ',  Ndm = ',Ndm
      call Finalize
    end if

end subroutine ForceDPD_Smooth22_SC


!#####################################################################
!#####################################################################


subroutine ForceDPD_Smooth23

use Configuration, only : R
use CommonDPD
use CellListMethod
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Vij, Vi, Fij
real(8) :: InvR, R1, R2, Wr, Wd, theta
real(8) :: rv, fc, fd, fr, ff, Slide
real(8) :: ranf, shv
external ranf
! >>
real(8) :: ccc, R7, R5, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDPt = Frc_Bond

! >>
   ccc = 16.d0 / 35.d0
! <<

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         Vi    = Velt(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij    = Ri - R(:,j)

              Slide  = -nint( Rij(2) * InvCL(2) )
              Rij(1) = Rij(1) + Slide * SlideGap

              Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
              Rij(2) = Rij(2) + Slide * CellL(2)
              Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wr   = 1.d0 - R1
                Wd   = Wr * Wr

                Sij = Rij * InvR

! >>
                Wc   = 1.d0 - R2

                theta = ranf() - 0.5d0

                Vij    = Vi - Velt(:,j)
                Vij(1) = Vij(1) + Slide * shv
                rv     = dot_product( Vij, Sij )

! >>
                aaa =   a(is,TypeNum(j)) 
                fc  =   aaa * Wc * Wc * Wc
                fd  = - gamma * Wd * rv
                fr  =   sigmt * Wr * theta

                ff =   fc + fd + fr

! >>
                R7 =   R2 * R2 * R2 / 7.d0
                R5 = - R2 * R2 * 3.d0 / 5.d0
                R3 =   R2

                PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 - 1.d0 ) + ccc )

                Fij = ff * Sij

                Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                FrcDPt(:,i) = FrcDPt(:,i) + Fij

                FrcDPt(:,j) = FrcDPt(:,j) - Fij

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1
                    Wr   = 1.d0 - R1
                    Wd   = Wr * Wr

                    Sij = Rij * InvR

! >>
                    Wc   = 1.d0 - R2

                    theta = ranf()-0.5d0

                    Vij = Vi - Velt(:,j)
                    Vij(1) = Vij(1) + Slide * shv
                    rv  = dot_product( Vij, Sij )

! >>
                    aaa =   a(is,TypeNum(j)) 
                    fc  =   aaa * Wc * Wc * Wc
                    fd = - gamma * Wd * rv
                    fr =   sigmt * Wr * theta

                    ff =   fc + fd + fr

! >>
                    R7 =   R2 * R2 * R2 / 7.d0
                    R5 = - R2 * R2 * 3.d0 / 5.d0
                    R3 =   R2

                    PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 - 1.d0 ) + ccc )

                    Fij = ff * Sij

                    Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                    Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                    Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                    Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                    Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                    Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                    Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                    Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                    Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                    FrcDPt(:,i) = FrcDPt(:,i) + Fij

                    FrcDPt(:,j) = FrcDPt(:,j) - Fij

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Vi    = Velt(:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1
               Wr   = 1.0d0 - R1
               Wd   = Wr * Wr

               Sij = Rij * InvR
! >>
               Wc   = 1.d0 - R2

               theta = ranf() - 0.5d0

               Vij = Vi - Velt(:,j)
               rv  = dot_product( Vij, Sij )

! >>
               aaa =   a(is,TypeNum(j)) 
               fc  =   aaa * Wc * Wc * Wc
               fd = - gamma * Wd * rv
               fr =   sigmt * Wr * theta

               ff =   fc + fd + fr

! >>
               R7 =   R2 * R2 * R2 / 7.d0
               R5 = - R2 * R2 * 3.d0 / 5.d0
               R3 =   R2

               PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 -1.d0 ) + ccc )

               Fij = ff * Sij

               Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
               Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
               Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
               Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
               Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
               Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

               FrcDPt(:,i) = FrcDPt(:,i) + Fij

               FrcDPt(:,j) = FrcDPt(:,j) - Fij

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1
                   Wr   = 1.d0 - R1
                   Wd   = Wr * Wr

                   Sij  = Rij * InvR

! >>
                   Wc   = 1.d0 - R2
                   theta = ranf() - 0.5d0

                   Vij = Vi - Velt(:,j)
                   rv  = dot_product( Vij, Sij )

! >>
                   aaa =   a(is,TypeNum(j)) 
                   fc  =   aaa * Wc * Wc * Wc
                   fd = - gamma * Wd * rv
                   fr =   sigmt * Wr * theta

                   ff =   fc + fd + fr

! >>
                   R7 =   R2 * R2 * R2 / 7.d0
                   R5 = - R2 * R2 * 3.d0 / 5.d0
                   R3 =   R2

                   PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 -1.d0 ) + ccc )

                   Fij = ff * Sij

                   Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                   Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                   Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                   Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                   Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                   Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                   FrcDPt(:,i) = FrcDPt(:,i) + Fij

                   FrcDPt(:,j) = FrcDPt(:,j) - Fij

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     Virial(2,1) = Virial(1,2)
     Virial(3,1) = Virial(1,3)
     Virial(3,2) = Virial(2,3)
!
! #####################################################################
    end if
! #####################################################################

end subroutine ForceDPD_Smooth23


!#####################################################################
!#####################################################################


subroutine ForceDPD_Smooth23_SC

use Configuration, only : R
use CommonDPD
use CellListMethod
use BookParam, only : Npair, ListIJ
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Fcij
#ifdef DPDcheck
real(8), dimension(3) :: Frij
#endif
real(8) :: InvR, R1, R2, Wr, Wd, theta
real(8) :: fc, fr, ff, Slide
real(8) :: ranf, shv
external ranf
! >>
real(8) :: ccc, R7, R5, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDPt = Frc_Bond

   VirialC = Vir_Bond
#ifdef DPDcheck
   VirialR = 0.d0
#endif
! >>
   ccc = 16.d0 / 35.d0
! <<

   Npair = 0
   SLList = 0.d0

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij    = Ri - R(:,j)
              Rij = Rij - nint( Rij * InvCL ) * CellL

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wr   = 1.d0 - R1
                Wd   = Wr * Wr

                Sij = Rij * InvR
! >>
                Wc   = 1.d0 - R2

                theta = ranf() - 0.5d0

! >> for dissipative force

                Npair = Npair + 1

                ListIJ(1,Npair) = i
                ListIJ(2,Npair) = j

                dRList(:,Npair) = Sij
                R1List(Npair) =  R1
                pfList(Npair) = -gamma * Wd

! << for dissipative force

! >>
                aaa =   a(is,TypeNum(j)) 
                fc  =   aaa * Wc * Wc * Wc
                fr  =   sigmt * Wr * theta

                ff =   fc + fr

! >>
                R7 =   R2 * R2 * R2 / 7.d0
                R5 = - R2 * R2 * 3.d0 / 5.d0
                R3 =   R2

                PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 - 1.d0 ) + ccc )

                Fcij = fc * Sij

                VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)

                FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij

                FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1
                    Wr   = 1.d0 - R1
                    Wd   = Wr * Wr

                    Sij = Rij * InvR
! >>
                    Wc   = 1.d0 - R2

                    theta = ranf()-0.5d0

! >> for dissipative force

                    Npair = Npair + 1

                    ListIJ(1,Npair) = i
                    ListIJ(2,Npair) = j

                    dRList(:,Npair) = Sij
                    R1List(Npair) =  R1
                    pfList(Npair) = -gamma * Wd
                    SLList(Npair) =  Slide * shv

! << for dissipative force

! >>
                    aaa =   a(is,TypeNum(j)) 
                    fc  =   aaa * Wc * Wc * Wc
                    fr  =   sigmt * Wr * theta

                    ff =   fc + fr

! >>
                    R7 =   R2 * R2 * R2 / 7.d0
                    R5 = - R2 * R2 * 3.d0 / 5.d0
                    R3 =   R2

                    PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 - 1.d0 ) + ccc )

                    Fcij = fc * Sij

                    VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                    VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                    VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                    VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                    VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                    VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)

                    FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij

                    FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     Virial = VirialC

! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1
               Wr   = 1.0d0 - R1
               Wd   = Wr * Wr

               Sij = Rij * InvR
! >>
               Wc   = 1.d0 - R2

               theta = ranf() - 0.5d0

! >> for dissipative force

               Npair = Npair + 1

               ListIJ(1,Npair) = i
               ListIJ(2,Npair) = j

               dRList(:,Npair) = Sij
               R1List(Npair) =  R1
               pfList(Npair) = -gamma * Wd

! << for dissipative force

! >>
               aaa =   a(is,TypeNum(j)) 
               fc  =   aaa * Wc * Wc * Wc
               fr  =   sigmt * Wr * theta

               ff =   fc + fr

! >>
               R7 =   R2 * R2 * R2 / 7.d0
               R5 = - R2 * R2 * 3.d0 / 5.d0
               R3 =   R2

               PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 -1.d0 ) + ccc )

               Fcij = fc * Sij

               VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
               VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
               VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
               VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
               VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
               VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)
#ifdef DPDcheck
               Frij = fr * Sij
               VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
               VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
               VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
               VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
               VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
               VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

               FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij

               FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1
                   Wr   = 1.d0 - R1
                   Wd   = Wr * Wr

                   Sij  = Rij * InvR
! >>
                   Wc   = 1.d0 - R2

                   theta = ranf() - 0.5d0

! >> for dissipative force

                   Npair = Npair + 1

                   ListIJ(1,Npair) = i
                   ListIJ(2,Npair) = j

                   dRList(:,Npair) = Sij
                   R1List(Npair) =  R1
                   pfList(Npair) = -gamma * Wd

! << for dissipative force

! >>
                   aaa =   a(is,TypeNum(j)) 
                   fc  =   aaa * Wc * Wc * Wc
                   fr  =   sigmt * Wr * theta

                   ff =   fc + fr

! >>
                   R7 =   R2 * R2 * R2 / 7.d0
                   R5 = - R2 * R2 * 3.d0 / 5.d0
                   R3 =   R2

                   PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 -1.d0 ) + ccc )

                   Fcij = fc * Sij

                   VirialC(1,1) = VirialC(1,1) + Fcij(1) * Rij(1)
                   VirialC(1,2) = VirialC(1,2) + Fcij(1) * Rij(2)
                   VirialC(1,3) = VirialC(1,3) + Fcij(1) * Rij(3)
                   VirialC(2,2) = VirialC(2,2) + Fcij(2) * Rij(2)
                   VirialC(2,3) = VirialC(2,3) + Fcij(2) * Rij(3)
                   VirialC(3,3) = VirialC(3,3) + Fcij(3) * Rij(3)
#ifdef DPDcheck
                   Frij = fr * Sij
                   VirialR(1,1) = VirialR(1,1) + Frij(1) * Rij(1)
                   VirialR(1,2) = VirialR(1,2) + Frij(1) * Rij(2)
                   VirialR(1,3) = VirialR(1,3) + Frij(1) * Rij(3)
                   VirialR(2,2) = VirialR(2,2) + Frij(2) * Rij(2)
                   VirialR(2,3) = VirialR(2,3) + Frij(2) * Rij(3)
                   VirialR(3,3) = VirialR(3,3) + Frij(3) * Rij(3)
#endif

                   FrcDPt(:,i) = FrcDPt(:,i) + ff * Sij

                   FrcDPt(:,j) = FrcDPt(:,j) - ff * Sij

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     Virial = VirialC
!
! #####################################################################
    end if
! #####################################################################

    if(Npair > Ndm) then
      write(*,*) 'ERROR : Npair exceeds its limit number'
      write(*,*) 'Npair = ', Npair, ',  Ndm = ',Ndm
      call Finalize
    end if

end subroutine ForceDPD_Smooth23_SC


!#####################################################################
!#####################################################################


subroutine ForceDPD_Attractive

use Configuration, only : R
use CommonDPD
use CellListMethod

implicit none

end subroutine ForceDPD_Attractive


!#####################################################################
!#####################################################################


subroutine ForceDPD_Attractive_SC

use Configuration, only : R
use CommonDPD
use CellListMethod

implicit none

end subroutine ForceDPD_Attractive_SC


!#####################################################################
!#####################################################################


subroutine Energy_Minimization_DPD

use Numbers, only : N, Nf
use CommonBlocks, only : ForceField
use Configuration, only : R
use CommonDPD
use Eminparam
use CellParam, only : CellL, InvCL

implicit none

!*******************************************
!*     minimize the energy                 *
!*     Steepest Decent Algorithm           *
!*******************************************

real(8), parameter :: dRmaxDP = 0.1
real(8), dimension(3,N) :: Rtmp
real(8) :: sumF2, RmsF, RmsF_o
real(8) :: deltaE, Energy, Energy_Pre
real(8) :: deltaR
integer :: Itry, i

   print *, 'Enegy minimization has just started'

   deltaR = dRmaxDP

! ------------------------------------
   if(ForceField=='F11') then
     call ForceDPD_Minim_Original
   else if(ForceField=='F12') then
     call ForceDPD_Minim_Smooth12
   else if(ForceField=='F21') then
     call ForceDPD_Minim_Smooth21
   else if(ForceField=='F22') then
     call ForceDPD_Minim_Smooth22
   else if(ForceField=='F23') then
     call ForceDPD_Minim_Smooth23
   else if(ForceField=='Morse') then
     call ForceDPD_Minim_Attractive
   end if
! ------------------------------------

   Energy_Pre = PotDP

   write( 6,*) 'Initial Energy             : ',Energy_Pre
#ifndef BMONI
   write(11,*) 'Initial Energy             : ',Energy_Pre
#endif

   sumF2 = 0.d0

   do i = 1 , N

     sumF2 = sumF2 + dot_product( FrcDP(:,i), FrcDP(:,i) )

   end do

   RmsF = sqrt( sumF2 / Nf )

!***************************************
!* Calculate Maximum Downhill Gradient *
!***************************************

   do Itry = 1 , MinTry

     RmsF_o = RmsF
     Rtmp   = R
     FrcDPt = FrcDP

     R = R + deltaR * FrcDP / RmsF

     do i = 1 , N

       R(:,i) = R(:,i) - nint( R(:,i) * InvCL ) * CellL

     end do

! -----------------------------------------
     if(ForceField=='F11') then
       call ForceDPD_Minim_Original
     else if(ForceField=='F12') then
       call ForceDPD_Minim_Smooth12
     else if(ForceField=='F21') then
       call ForceDPD_Minim_Smooth21
     else if(ForceField=='F22') then
       call ForceDPD_Minim_Smooth22
     else if(ForceField=='F23') then
       call ForceDPD_Minim_Smooth23
     else if(ForceField=='Morse') then
       call ForceDPD_Minim_Attractive
     end if
! -----------------------------------------

     Energy = PotDP

     deltaE = Energy - Energy_Pre

     sumF2 = 0.d0

     do i = 1 , N

       sumF2 = sumF2 + dot_product( FrcDP(:,i), FrcDP(:,i) )

     end do

     RmsF = sqrt( sumF2 / Nf )

     write(*,'(a,i3,a,e12.4)') 'Step = ',Itry,' :  Energy = ', Energy

     if( ( abs(deltaE/Energy) < dev_relative ).or.(Energy == 0.) ) exit

     if( Energy < Energy_Pre ) then

       Energy_Pre = Energy
       deltaR     = deltaR * 1.2d0
       deltaR     = min(deltaR , dRmaxDP)

     else

       FrcDP  = FrcDPt
       R      = Rtmp
       RmsF   = RmsF_o
       deltaR = deltaR * 0.5d0

     end if

   end do

#ifndef BMONI
   write(11,*) 'Final Energy               : ',Energy
#endif

   print *, 'Enegy minimization has just finished'
   print *, 'Final Energy               : ',Energy

end subroutine Energy_Minimization_DPD


!#####################################################################
!#####################################################################


subroutine ForceDPD_Minim_Original

use Configuration, only : R
use CommonDPD
use CellListMethod
use BondedParam, only : Frc_Bond, Ene_Bond
use CellParam, only : CellL, InvCL

implicit none

integer :: icell, i, j, Is, NumI
real(8), dimension(3) :: Ri, Rij, Fij
real(8) :: R2, R1, InvR
real(8) :: Wr, fc
integer :: jcell0, nabor, jcell
integer, dimension(4) :: NonFc

   call Force_Bond_DP

   PotDP  = Ene_Bond
   FrcDP  = Frc_Bond

! ---------------------
   call LinkCell
! ---------------------
!
   do icell = 1 , Ncell

     i = Head(icell)
!
!       ** loop over all particles in the cell **
!
     do while( i > 0 )

       Ri    = R(:,i)
       Is    = TypeNum(i)
       NumI  = Ncal(i)
       NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

       j = NextP(i)

       do while( j > 0 )

         if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
         &                       .and. ( j /= NonFc(2) ) &
         &                       .and. ( j /= NonFc(3) ) &
         &                       .and. ( j /= NonFc(4) ) ) then

           Rij = Ri - R(:,j)
           Rij = Rij - nint( Rij * InvCL ) * CellL
           R2  = dot_product( Rij , Rij )

           if(R2 < 1.0d0) then

             R1   = sqrt(R2)
             InvR = 1.0d0 / R1
             Wr   = 1.0d0 - R1

             fc = a(Is,TypeNum(j)) * Wr

             PotDP = PotDP + 0.5d0 * fc * Wr

             Fij = fc * Rij * InvR

             FrcDP(:,i) = FrcDP(:,i) + Fij

             FrcDP(:,j) = FrcDP(:,j) - Fij

           end if

         end if

         j = NextP(j)

       end do

!** loop over neighbouring cells **

       jcell0 = 16 * ( icell - 1 )

       do nabor = 1 , 16

         jcell = Map( jcell0 + nabor )

         if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

           j = Head(jcell)

           do while( j /= 0 )

             if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
             &                       .and. ( j /= NonFc(2) ) &
             &                       .and. ( j /= NonFc(3) ) &
             &                       .and. ( j /= NonFc(4) ) ) then

               Rij = Ri - R(:,j)
               Rij = Rij - nint( Rij*InvCL ) * CellL
               R2  = dot_product( Rij, Rij )

               if(R2 < 1.d0) then

                 R1   = sqrt(R2)
                 InvR = 1.d0 / R1
                 Wr   = 1.d0 - R1

                 fc = a(Is,TypeNum(j)) * Wr

                 PotDP = PotDP + 0.5d0 * fc * Wr

                 Fij = fc * Rij * InvR

                 FrcDP(:,i) = FrcDP(:,i) + Fij

                 FrcDP(:,j) = FrcDP(:,j) - Fij

               end if

             end if

             j = NextP(j)

           end do

         end if

       end do

       i = NextP(i)

     end do

   end do

end subroutine ForceDPD_Minim_Original


!#####################################################################
!#####################################################################


subroutine ForceDPD_Minim_Smooth12

use Configuration, only : R
use CommonDPD
use CellListMethod
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Fij
real(8) :: InvR, R1, R2
real(8) :: fc
! >>
real(8) :: ccc, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDP  = Frc_Bond

! >>
   ccc = 1.d0 / 3.d0
! <<

! ---------------------
   call LinkCell
! ---------------------
!
   do icell = 1 , Ncell

     i = Head(icell)
!
!    ** loop over all particles in the cell **
!
     do while( i > 0 )

       Ri    = R (:,i)
       Is    = TypeNum(i)
       NumI  = Ncal(i)
       NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

       j = NextP(i)

       do while( j > 0 )

         if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
         &                       .and. ( j /= NonFc(2) ) &
         &                       .and. ( j /= NonFc(3) ) &
         &                       .and. ( j /= NonFc(4) ) ) then

           Rij = Ri - R(:,j)
           Rij = Rij - nint( Rij * InvCL ) * CellL
           R2  = dot_product( Rij , Rij )

           if(R2 < 1.0d0) then

             R1   = sqrt(R2)
             InvR = 1.d0 / R1
! >>
             Wc   = 1.d0 - R1

! >>
             aaa =   a(is,TypeNum(j)) 
             fc  =   aaa * Wc * Wc

! >>
             R3 =   R2 / 3.d0

             PotDP = PotDP + aaa * ( R1*( -R3 + R1 -1.d0 ) + ccc )

             Fij = fc * Rij * InvR

             FrcDP(:,i) = FrcDP(:,i) + Fij

             FrcDP(:,j) = FrcDP(:,j) - Fij

           end if

         end if

         j = NextP(j)

       end do

!** loop over neighbouring cells **

       jcell0 = 16 * ( icell - 1 )

       do nabor = 1 , 16

         jcell = Map( jcell0 + nabor )

         if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

           j = Head(jcell)

           do while( j /= 0 )

             if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
             &                       .and. ( j /= NonFc(2) ) &
             &                       .and. ( j /= NonFc(3) ) &
             &                       .and. ( j /= NonFc(4) ) ) then

               Rij = Ri - R(:,j)
               Rij = Rij - nint( Rij*InvCL ) * CellL
               R2  = dot_product( Rij, Rij )

               if(R2 < 1.d0) then

                 R1   = sqrt(R2)
                 InvR = 1.d0 / R1
! >>
                 Wc   = 1.d0 - R1
! >>
                 aaa =   a(is,TypeNum(j)) 
                 fc  =   aaa * Wc * Wc
! >>
                 R3 =   R2 / 3.d0

                 PotDP = PotDP + aaa * ( R1*( - R3 + R1 - 1.d0 ) + ccc )

                 Fij = fc * Rij * InvR

                 FrcDP(:,i) = FrcDP(:,i) + Fij

                 FrcDP(:,j) = FrcDP(:,j) - Fij

               end if

             end if

             j = NextP(j)

           end do

         end if

       end do

       i = NextP(i)

     end do

   end do

end subroutine ForceDPD_Minim_Smooth12


!#####################################################################
!#####################################################################


subroutine ForceDPD_Minim_Smooth21

use Configuration, only : R
use CommonDPD
use CellListMethod
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Fij
real(8) :: InvR, R1, R2
real(8) :: fc
! >>
real(8) :: ccc, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDP  = Frc_Bond

! >>
   ccc = 2.d0 / 3.d0
! <<

! ---------------------
   call LinkCell
! ---------------------
!
   do icell = 1 , Ncell

     i = Head(icell)
!
!    ** loop over all particles in the cell **
!
     do while( i > 0 )

       Ri    = R (:,i)
       Is    = TypeNum(i)
       NumI  = Ncal(i)
       NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

       j = NextP(i)

       do while( j > 0 )

         if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
         &                       .and. ( j /= NonFc(2) ) &
         &                       .and. ( j /= NonFc(3) ) &
         &                       .and. ( j /= NonFc(4) ) ) then

           Rij = Ri - R(:,j)
           Rij = Rij - nint( Rij * InvCL ) * CellL
           R2  = dot_product( Rij , Rij )

           if(R2 < 1.0d0) then

             R1   = sqrt(R2)
             InvR = 1.d0 / R1
! >>
             Wc   = 1.d0 - R2

! >>
             aaa =   a(is,TypeNum(j)) 
             fc  =   aaa * Wc

! >>
             R3 =   R2 / 3.d0

             PotDP = PotDP + aaa * ( R1*( R3 -1.d0 ) + ccc )

             Fij = fc * Rij * InvR

             FrcDP(:,i) = FrcDP(:,i) + Fij

             FrcDP(:,j) = FrcDP(:,j) - Fij

           end if

         end if

         j = NextP(j)

       end do

!** loop over neighbouring cells **

       jcell0 = 16 * ( icell - 1 )

       do nabor = 1 , 16

         jcell = Map( jcell0 + nabor )

         if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

           j = Head(jcell)

           do while( j /= 0 )

             if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
             &                       .and. ( j /= NonFc(2) ) &
             &                       .and. ( j /= NonFc(3) ) &
             &                       .and. ( j /= NonFc(4) ) ) then

               Rij = Ri - R(:,j)
               Rij = Rij - nint( Rij*InvCL ) * CellL
               R2  = dot_product( Rij, Rij )

               if(R2 < 1.d0) then

                 R1   = sqrt(R2)
                 InvR = 1.d0 / R1
! >>
                 Wc   = 1.d0 - R2
! >>
                 aaa =   a(is,TypeNum(j)) 
                 fc  =   aaa * Wc
! >>
                 R3 =   R2 / 3.d0

                 PotDP = PotDP + aaa * ( R1*( R3 -1.d0 ) + ccc )

                 Fij = fc * Rij * InvR

                 FrcDP(:,i) = FrcDP(:,i) + Fij

                 FrcDP(:,j) = FrcDP(:,j) - Fij

               end if

             end if

             j = NextP(j)

           end do

         end if

       end do

       i = NextP(i)

     end do

   end do

end subroutine ForceDPD_Minim_Smooth21


!#####################################################################
!#####################################################################


subroutine ForceDPD_Minim_Smooth22

use Configuration, only : R
use CommonDPD
use CellListMethod
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Fij
real(8) :: InvR, R1, R2
real(8) :: fc
! >>
real(8) :: ccc, R5, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDP  = Frc_Bond

! >>
   ccc = 8.d0 / 15.d0
! <<

! ---------------------
   call LinkCell
! ---------------------
!
   do icell = 1 , Ncell

     i = Head(icell)
!
!    ** loop over all particles in the cell **
!
     do while( i > 0 )

       Ri    = R (:,i)
       Is    = TypeNum(i)
       NumI  = Ncal(i)
       NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

       j = NextP(i)

       do while( j > 0 )

         if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
         &                       .and. ( j /= NonFc(2) ) &
         &                       .and. ( j /= NonFc(3) ) &
         &                       .and. ( j /= NonFc(4) ) ) then

           Rij = Ri - R(:,j)
           Rij = Rij - nint( Rij * InvCL ) * CellL
           R2  = dot_product( Rij , Rij )

           if(R2 < 1.0d0) then

             R1   = sqrt(R2)
             InvR = 1.d0 / R1
! >>
             Wc   = 1.d0 - R2

! >>
             aaa =   a(is,TypeNum(j)) 
             fc  =   aaa * Wc * Wc

! >>
             R5 = - R2 * R2 / 5.d0
             R3 =   R2 * 2.d0 / 3.d0

             PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

             Fij = fc * Rij * InvR

             FrcDP(:,i) = FrcDP(:,i) + Fij

             FrcDP(:,j) = FrcDP(:,j) - Fij

           end if

         end if

         j = NextP(j)

       end do

!** loop over neighbouring cells **

       jcell0 = 16 * ( icell - 1 )

       do nabor = 1 , 16

         jcell = Map( jcell0 + nabor )

         if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

           j = Head(jcell)

           do while( j /= 0 )

             if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
             &                       .and. ( j /= NonFc(2) ) &
             &                       .and. ( j /= NonFc(3) ) &
             &                       .and. ( j /= NonFc(4) ) ) then

               Rij = Ri - R(:,j)
               Rij = Rij - nint( Rij*InvCL ) * CellL
               R2  = dot_product( Rij, Rij )

               if(R2 < 1.d0) then

                 R1   = sqrt(R2)
                 InvR = 1.d0 / R1
! >>
                 Wc   = 1.d0 - R2
! >>
                 aaa =   a(is,TypeNum(j)) 
                 fc  =   aaa * Wc * Wc
! >>
                 R5 = - R2 * R2 / 5.d0
                 R3 =   R2 * 2.d0 / 3.d0

                 PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                 Fij = fc * Rij * InvR

                 FrcDP(:,i) = FrcDP(:,i) + Fij

                 FrcDP(:,j) = FrcDP(:,j) - Fij

               end if

             end if

             j = NextP(j)

           end do

         end if

       end do

       i = NextP(i)

     end do

   end do

end subroutine ForceDPD_Minim_Smooth22


!#####################################################################
!#####################################################################


subroutine ForceDPD_Minim_Smooth23

use Configuration, only : R
use CommonDPD
use CellListMethod
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none


integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Fij
real(8) :: InvR, R1, R2
real(8) :: fc
! >>
real(8) :: ccc, R7, R5, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDP  = Frc_Bond

! >>
   ccc = 16.d0 / 35.d0
! <<

! ---------------------
   call LinkCell
! ---------------------
!
   do icell = 1 , Ncell

     i = Head(icell)
!
!    ** loop over all particles in the cell **
!
     do while( i > 0 )

       Ri    = R (:,i)
       Is    = TypeNum(i)
       NumI  = Ncal(i)
       NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

       j = NextP(i)

       do while( j > 0 )

         if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
         &                       .and. ( j /= NonFc(2) ) &
         &                       .and. ( j /= NonFc(3) ) &
         &                       .and. ( j /= NonFc(4) ) ) then

           Rij = Ri - R(:,j)
           Rij = Rij - nint( Rij * InvCL ) * CellL
           R2  = dot_product( Rij , Rij )

           if(R2 < 1.0d0) then

             R1   = sqrt(R2)
             InvR = 1.d0 / R1
! >>
             Wc   = 1.d0 - R2

! >>
             aaa =   a(is,TypeNum(j)) 
             fc  =   aaa * Wc * Wc * Wc

! >>
             R7 =   R2 * R2 * R2 / 7.d0
             R5 = - R2 * R2 * 3.d0 / 5.d0
             R3 =   R2

             PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 - 1.d0 ) + ccc )

             Fij = fc * Rij * InvR

             FrcDP(:,i) = FrcDP(:,i) + Fij

             FrcDP(:,j) = FrcDP(:,j) - Fij

           end if

         end if

         j = NextP(j)

       end do

!** loop over neighbouring cells **

       jcell0 = 16 * ( icell - 1 )

       do nabor = 1 , 16

         jcell = Map( jcell0 + nabor )

         if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

           j = Head(jcell)

           do while( j /= 0 )

             if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
             &                       .and. ( j /= NonFc(2) ) &
             &                       .and. ( j /= NonFc(3) ) &
             &                       .and. ( j /= NonFc(4) ) ) then

               Rij = Ri - R(:,j)
               Rij = Rij - nint( Rij*InvCL ) * CellL
               R2  = dot_product( Rij, Rij )

               if(R2 < 1.d0) then

                 R1   = sqrt(R2)
                 InvR = 1.d0 / R1
! >>
                 Wc   = 1.d0 - R2
! >>
                 aaa =   a(is,TypeNum(j)) 
                 fc  =   aaa * Wc * Wc * Wc
! >>
                 R7 =   R2 * R2 * R2 / 7.d0
                 R5 = - R2 * R2 * 3.d0 / 5.d0
                 R3 =   R2

                 PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 - 1.d0 ) + ccc )

                 Fij = fc * Rij * InvR

                 FrcDP(:,i) = FrcDP(:,i) + Fij

                 FrcDP(:,j) = FrcDP(:,j) - Fij

               end if

             end if

             j = NextP(j)

           end do

         end if

       end do

       i = NextP(i)

     end do

   end do

end subroutine ForceDPD_Minim_Smooth23


!#####################################################################
!#####################################################################


subroutine ForceDPD_Minim_Attractive

use Configuration, only : R
use CommonDPD
use CellListMethod

implicit none

end subroutine ForceDPD_Minim_Attractive


!######################################################################
!######################################################################


! ***************************
! ** Bond Stretching Force **
! ***************************

subroutine Force_Bond_DP

use Configuration, only : R
use CommonDPD
use BondedParam, only : NumBond, BondI, BondJ, kBond, rBond, &
&   Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL

implicit NONE

integer :: i, j, k
real(8) :: R2, R1, sR, Fc
real(8), dimension(3) :: Rij, Fij

   Ene_Bond = 0.d0
   Frc_Bond = 0.d0
   Vir_Bond = 0.d0

   if( NumBond == 0 ) Return

   do k = 1 , NumBond

     i= BondI(k)
     j= BondJ(k)

     Rij = R(:,i) - R(:,j)
     Rij = Rij - nint(Rij * InvCL) * CellL
     R2  = dot_product(Rij,Rij)

     R1 = sqrt(R2)
     sR = R1 - rBond(k)
     Fc = kBond(k) * sR

     Ene_Bond = Ene_Bond + Fc * sR

     Fij = - 2.d0 * Fc * Rij / R1

     Frc_Bond(:,i) = Frc_Bond(:,i) + Fij
     Frc_Bond(:,j) = Frc_Bond(:,j) - Fij

     Vir_Bond(1,1) = Vir_Bond(1,1) + Fij(1) * Rij(1)
     Vir_Bond(1,2) = Vir_Bond(1,2) + Fij(1) * Rij(2)
     Vir_Bond(1,3) = Vir_Bond(1,3) + Fij(1) * Rij(3)
     Vir_Bond(2,2) = Vir_Bond(2,2) + Fij(2) * Rij(2)
     Vir_Bond(2,3) = Vir_Bond(2,3) + Fij(2) * Rij(3)
     Vir_Bond(3,3) = Vir_Bond(3,3) + Fij(3) * Rij(3)

   end do

   Vir_Bond(2,1) = Vir_Bond(1,2)
   Vir_Bond(3,1) = Vir_Bond(1,3)
   Vir_Bond(3,2) = Vir_Bond(2,3)

end subroutine Force_Bond_DP


!#####################################################################
!#####################################################################


subroutine ForceDPD_Original_Lowe

use Configuration, only : R
use CommonDPD
use CellListMethod
use BookParam, only : Npair, ListIJ
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Fij
real(8) :: InvR, R1, R2, Wr
real(8) :: fc
real(8) :: ranf, shv
integer :: Slide
external ranf

   call Force_Bond_DP

   PotDP   = Ene_Bond
   VirialC = Vir_Bond
   FrcDPt  = Frc_Bond

   Npair = 0
   SLList = 0.d0

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij = Ri - R(:,j)
              Rij = Rij - nint( Rij * InvCL ) * CellL

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wr   = 1.d0 - R1

                Sij = Rij * InvR

! >> for dissipative force

                Npair = Npair + 1

                ListIJ(1,Npair) = i
                ListIJ(2,Npair) = j

                R1List(Npair) = R1

                dRList(:,Npair) = Sij

! << for dissipative force

                fc =   a(is,TypeNum(j)) * Wr

                PotDP = PotDP + 0.5d0 * fc * Wr

                Fij = fc * Sij

                VirialC(1,1) = VirialC(1,1) + Fij(1) * Rij(1)
                VirialC(1,2) = VirialC(1,2) + Fij(1) * Rij(2)
                VirialC(1,3) = VirialC(1,3) + Fij(1) * Rij(3)
                VirialC(2,2) = VirialC(2,2) + Fij(2) * Rij(2)
                VirialC(2,3) = VirialC(2,3) + Fij(2) * Rij(3)
                VirialC(3,3) = VirialC(3,3) + Fij(3) * Rij(3)

                FrcDPt(:,i) = FrcDPt(:,i) + Fij

                FrcDPt(:,j) = FrcDPt(:,j) - Fij

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1
                    Wr   = 1.d0 - R1

                    Sij = Rij * InvR

! >> for dissipative force

                    Npair = Npair + 1

                    ListIJ(1,Npair) = i
                    ListIJ(2,Npair) = j

                    R1List(Npair) = R1

                    dRList(:,Npair) = Sij
                    SLList(Npair) =  Slide * shv

! << for dissipative force

                    fc =   a(is,TypeNum(j)) * Wr

                    PotDP = PotDP + 0.5d0 * fc * Wr

                    Fij = fc * Sij

                    VirialC(1,1) = VirialC(1,1) + Fij(1) * Rij(1)
                    VirialC(1,2) = VirialC(1,2) + Fij(1) * Rij(2)
                    VirialC(1,3) = VirialC(1,3) + Fij(1) * Rij(3)
                    VirialC(2,2) = VirialC(2,2) + Fij(2) * Rij(2)
                    VirialC(2,3) = VirialC(2,3) + Fij(2) * Rij(3)
                    VirialC(3,3) = VirialC(3,3) + Fij(3) * Rij(3)

                    FrcDPt(:,i) = FrcDPt(:,i) + Fij

                    FrcDPt(:,j) = FrcDPt(:,j) - Fij

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

      VirialC(2,1) = VirialC(1,2)
      VirialC(3,1) = VirialC(1,3)
      VirialC(3,2) = VirialC(2,3)

      Virial = VirialC

! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1
               Wr   = 1.0d0 - R1

               Sij = Rij * InvR

! >> for dissipative force

               Npair = Npair + 1

               ListIJ(1,Npair) = i
               ListIJ(2,Npair) = j

               R1List(Npair) = R1
!               R1List(Npair) = R1

               dRList(:,Npair) = Sij

! << for dissipative force

               fc =   a(Is,TypeNum(j)) * Wr

               PotDP = PotDP + 0.5d0 * fc * Wr

               Fij = fc * Sij

               VirialC(1,1) = VirialC(1,1) + Fij(1) * Rij(1)
               VirialC(1,2) = VirialC(1,2) + Fij(1) * Rij(2)
               VirialC(1,3) = VirialC(1,3) + Fij(1) * Rij(3)
               VirialC(2,2) = VirialC(2,2) + Fij(2) * Rij(2)
               VirialC(2,3) = VirialC(2,3) + Fij(2) * Rij(3)
               VirialC(3,3) = VirialC(3,3) + Fij(3) * Rij(3)

               FrcDPt(:,i) = FrcDPt(:,i) + Fij

               FrcDPt(:,j) = FrcDPt(:,j) - Fij

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1
                   Wr   = 1.d0 - R1

                   Sij  = Rij * InvR

! >> for dissipative force

                   Npair = Npair + 1

                   ListIJ(1,Npair) = i
                   ListIJ(2,Npair) = j

                   R1List(Npair) = R1

                   dRList(:,Npair) = Sij

! << for dissipative force

                   fc =   a(is,TypeNum(j)) * Wr

                   PotDP = PotDP + 0.5d0 * fc * Wr

                   Fij = fc * Sij

                   VirialC(1,1) = VirialC(1,1) + Fij(1) * Rij(1)
                   VirialC(1,2) = VirialC(1,2) + Fij(1) * Rij(2)
                   VirialC(1,3) = VirialC(1,3) + Fij(1) * Rij(3)
                   VirialC(2,2) = VirialC(2,2) + Fij(2) * Rij(2)
                   VirialC(2,3) = VirialC(2,3) + Fij(2) * Rij(3)
                   VirialC(3,3) = VirialC(3,3) + Fij(3) * Rij(3)

                   FrcDPt(:,i) = FrcDPt(:,i) + Fij

                   FrcDPt(:,j) = FrcDPt(:,j) - Fij

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     VirialC(2,1) = VirialC(1,2)
     VirialC(3,1) = VirialC(1,3)
     VirialC(3,2) = VirialC(2,3)

     Virial = VirialC
!
! #####################################################################
    end if
! #####################################################################

    if(Npair > Ndm) then
      write(*,*) 'ERROR : Npair exceeds its limit number'
      write(*,*) 'Npair = ', Npair, ',  Ndm = ',Ndm
      call Finalize
    end if

end subroutine ForceDPD_Original_Lowe


!#####################################################################
!#####################################################################


subroutine ForceDPD_Smooth12_Lowe

use Configuration, only : R
use CommonDPD
use CellListMethod
use BookParam, only : Npair, ListIJ
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Fij
real(8) :: InvR, R1, R2
real(8) :: fc, Slide
real(8) :: ranf, shv
external ranf
! >>
real(8) :: ccc, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDPt = Frc_Bond

   Npair = 0
   SLList = 0.d0

! >>
   ccc = 1.d0 / 3.d0
! <<

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij = Ri - R(:,j)
              Rij = Rij - nint( Rij * InvCL ) * CellL

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wc   = 1.d0 - R1

                Sij = Rij * InvR

! >> for dissipative force

                Npair = Npair + 1

                ListIJ(1,Npair) = i
                ListIJ(2,Npair) = j

                R1List(Npair) = R1

                dRList(:,Npair) = Sij

! << for dissipative force

                aaa =   a(is,TypeNum(j)) 
                fc =   aaa * Wc * Wc

                R3 =   R2 / 3.d0

                PotDP = PotDP + aaa * ( R1*( - R3 + R1 - 1.d0 ) + ccc )

                Fij = fc * Sij

                Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                FrcDPt(:,i) = FrcDPt(:,i) + Fij

                FrcDPt(:,j) = FrcDPt(:,j) - Fij

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1

                    Sij = Rij * InvR

                    Wc   = 1.d0 - R1

! >> for dissipative force

                    Npair = Npair + 1

                    ListIJ(1,Npair) = i
                    ListIJ(2,Npair) = j

                    R1List(Npair) = R1

                    dRList(:,Npair) = Sij
                    SLList(Npair) =  Slide * shv

! << for dissipative force

                    aaa =   a(is,TypeNum(j)) 
                    fc =   aaa * Wc * Wc

                    R3 =   R2 / 3.d0

                    PotDP = PotDP + aaa * ( R1*( - R3 + R1 - 1.d0 ) + ccc )

                    Fij = fc * Sij

                    Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                    Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                    Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                    Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                    Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                    Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                    Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                    Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                    Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                    FrcDPt(:,i) = FrcDPt(:,i) + Fij

                    FrcDPt(:,j) = FrcDPt(:,j) - Fij

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1

               Sij = Rij * InvR
               Wc   = 1.d0 - R1

! >> for dissipative force

               Npair = Npair + 1

               ListIJ(1,Npair) = i
               ListIJ(2,Npair) = j

               R1List(Npair) = R1

               dRList(:,Npair) = Sij

! << for dissipative force

               aaa =   a(is,TypeNum(j)) 
               fc  =   aaa * Wc * Wc

               R3 =   R2 / 3.d0

               PotDP = PotDP + aaa * ( R1*( - R3 + R1 - 1.d0 ) + ccc )

               Fij = fc * Sij

               Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
               Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
               Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
               Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
               Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
               Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

               FrcDPt(:,i) = FrcDPt(:,i) + Fij

               FrcDPt(:,j) = FrcDPt(:,j) - Fij

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1

                   Sij  = Rij * InvR
                   Wc   = 1.d0 - R1

! >> for dissipative force

                   Npair = Npair + 1

                   ListIJ(1,Npair) = i
                   ListIJ(2,Npair) = j

                   R1List(Npair) = R1

                   dRList(:,Npair) = Sij

! << for dissipative force

                   aaa =   a(is,TypeNum(j)) 
                   fc  =   aaa * Wc * Wc

                   R3 =   R2 / 3.d0

                   PotDP = PotDP + aaa * ( R1*( - R3 + R1 - 1.d0 ) + ccc )

                   Fij = fc * Sij

                   Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                   Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                   Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                   Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                   Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                   Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                   FrcDPt(:,i) = FrcDPt(:,i) + Fij

                   FrcDPt(:,j) = FrcDPt(:,j) - Fij

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     Virial(2,1) = Virial(1,2)
     Virial(3,1) = Virial(1,3)
     Virial(3,2) = Virial(2,3)
!
! #####################################################################
    end if
! #####################################################################

    if(Npair > Ndm) then
      write(*,*) 'ERROR : Npair exceeds its limit number'
      write(*,*) 'Npair = ', Npair, ',  Ndm = ',Ndm
      call Finalize
    end if

end subroutine ForceDPD_Smooth12_Lowe


!#####################################################################
!#####################################################################


subroutine ForceDPD_Smooth21_Lowe

use Configuration, only : R
use CommonDPD
use CellListMethod
use BookParam, only : Npair, ListIJ
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Fij
real(8) :: InvR, R1, R2
real(8) :: fc, Slide
real(8) :: ranf, shv
external ranf
! >>
real(8) :: ccc, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDPt = Frc_Bond

   Npair = 0
   SLList = 0.d0

! >>
   ccc = 2.d0 / 3.d0
! <<

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij    = Ri - R(:,j)

              Slide  = -nint( Rij(2) * InvCL(2) )
              Rij(1) = Rij(1) + Slide * SlideGap

              Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
              Rij(2) = Rij(2) + Slide * CellL(2)
              Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wc   = 1.d0 - R2

                Sij = Rij * InvR

! >> for dissipative force

                Npair = Npair + 1

                ListIJ(1,Npair) = i
                ListIJ(2,Npair) = j

                R1List(Npair) = R1

                dRList(:,Npair) = Sij
                SLList(Npair) =  Slide * shv

! << for dissipative force

                aaa =   a(is,TypeNum(j)) 
                fc =   aaa * Wc

                R3 =   R2 / 3.d0

                PotDP = PotDP + aaa * ( R1*( R3 -1.d0 ) + ccc )

                Fij = fc * Sij

                Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                FrcDPt(:,i) = FrcDPt(:,i) + Fij

                FrcDPt(:,j) = FrcDPt(:,j) - Fij

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1

                    Sij = Rij * InvR

                    Wc   = 1.d0 - R2

! >> for dissipative force

                    Npair = Npair + 1

                    ListIJ(1,Npair) = i
                    ListIJ(2,Npair) = j

                    R1List(Npair) = R1

                    dRList(:,Npair) = Sij
                    SLList(Npair) =  Slide * shv

! << for dissipative force

                    aaa =   a(is,TypeNum(j)) 
                    fc =   aaa * Wc

                    R3 =   R2 / 3.d0

                    PotDP = PotDP + aaa * ( R1*( R3 -1.d0 ) + ccc )

                    Fij = fc * Sij

                    Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                    Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                    Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                    Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                    Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                    Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                    Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                    Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                    Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                    FrcDPt(:,i) = FrcDPt(:,i) + Fij

                    FrcDPt(:,j) = FrcDPt(:,j) - Fij

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1

               Sij = Rij * InvR
               Wc   = 1.d0 - R2

! >> for dissipative force

               Npair = Npair + 1

               ListIJ(1,Npair) = i
               ListIJ(2,Npair) = j

               R1List(Npair) = R1
               R1List(Npair) = R1

               dRList(:,Npair) = Sij

! << for dissipative force

               aaa =   a(is,TypeNum(j)) 
               fc  =   aaa * Wc

               R3 =   R2 / 3.d0

               PotDP = PotDP + aaa * ( R1*( R3 -1.d0 ) + ccc )

               Fij = fc * Sij

               Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
               Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
               Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
               Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
               Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
               Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

               FrcDPt(:,i) = FrcDPt(:,i) + Fij

               FrcDPt(:,j) = FrcDPt(:,j) - Fij

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1

                   Sij  = Rij * InvR
                   Wc   = 1.d0 - R2

! >> for dissipative force

                   Npair = Npair + 1

                   ListIJ(1,Npair) = i
                   ListIJ(2,Npair) = j

                   R1List(Npair) = R1

                   dRList(:,Npair) = Sij

! << for dissipative force

                   aaa =   a(is,TypeNum(j)) 
                   fc  =   aaa * Wc

                   R3 =   R2 / 3.d0

                   PotDP = PotDP + aaa * ( R1*( R3 -1.d0 ) + ccc )

                   Fij = fc * Sij

                   Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                   Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                   Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                   Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                   Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                   Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                   FrcDPt(:,i) = FrcDPt(:,i) + Fij

                   FrcDPt(:,j) = FrcDPt(:,j) - Fij

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     Virial(2,1) = Virial(1,2)
     Virial(3,1) = Virial(1,3)
     Virial(3,2) = Virial(2,3)
!
! #####################################################################
    end if
! #####################################################################

    if(Npair > Ndm) then
      write(*,*) 'ERROR : Npair exceeds its limit number'
      write(*,*) 'Npair = ', Npair, ',  Ndm = ',Ndm
      call Finalize
    end if

end subroutine ForceDPD_Smooth21_Lowe


!#####################################################################
!#####################################################################


subroutine ForceDPD_Smooth22_Lowe

use Configuration, only : R
use CommonDPD
use CellListMethod
use BookParam, only : Npair, ListIJ
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Fij
real(8) :: InvR, R1, R2
real(8) :: fc, Slide
real(8) :: ranf, shv
external ranf
! >>
real(8) :: ccc, R5, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDPt = Frc_Bond

   Npair = 0
   SLList = 0.d0

! >>
   ccc = 8.d0 / 15.d0
! <<

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij    = Ri - R(:,j)

              Slide  = -nint( Rij(2) * InvCL(2) )
              Rij(1) = Rij(1) + Slide * SlideGap

              Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
              Rij(2) = Rij(2) + Slide * CellL(2)
              Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wc   = 1.d0 - R2

                Sij = Rij * InvR

! >> for dissipative force

                Npair = Npair + 1

                ListIJ(1,Npair) = i
                ListIJ(2,Npair) = j

                R1List(Npair) = R1

                dRList(:,Npair) = Sij
                SLList(Npair) =  Slide * shv

! << for dissipative force

                aaa =   a(is,TypeNum(j)) 
                fc =   aaa * Wc * Wc

                R5 = - R2 * R2 / 5.d0
                R3 =   R2 * 2.d0 / 3.d0

                PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                Fij = fc * Sij

                Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                FrcDPt(:,i) = FrcDPt(:,i) + Fij

                FrcDPt(:,j) = FrcDPt(:,j) - Fij

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1

                    Sij = Rij * InvR

                    Wc   = 1.d0 - R2

! >> for dissipative force

                    Npair = Npair + 1

                    ListIJ(1,Npair) = i
                    ListIJ(2,Npair) = j

                    R1List(Npair) = R1

                    dRList(:,Npair) = Sij
                    SLList(Npair) =  Slide * shv

! << for dissipative force

                    aaa =   a(is,TypeNum(j)) 
                    fc =   aaa * Wc * Wc

                    R5 = - R2 * R2 / 5.d0
                    R3 =   R2 * 2.d0 / 3.d0

                    PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                    Fij = fc * Sij

                    Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                    Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                    Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                    Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                    Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                    Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                    Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                    Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                    Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                    FrcDPt(:,i) = FrcDPt(:,i) + Fij

                    FrcDPt(:,j) = FrcDPt(:,j) - Fij

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1

               Sij = Rij * InvR
               Wc   = 1.d0 - R2

! >> for dissipative force

               Npair = Npair + 1

               ListIJ(1,Npair) = i
               ListIJ(2,Npair) = j

               R1List(Npair) = R1

               dRList(:,Npair) = Sij

! << for dissipative force

               aaa =   a(is,TypeNum(j)) 
               fc  =   aaa * Wc * Wc

               R5 = - R2 * R2 / 5.d0
               R3 =   R2 * 2.d0 / 3.d0

               PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

               Fij = fc * Sij

               Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
               Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
               Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
               Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
               Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
               Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

               FrcDPt(:,i) = FrcDPt(:,i) + Fij

               FrcDPt(:,j) = FrcDPt(:,j) - Fij

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1

                   Sij  = Rij * InvR
                   Wc   = 1.d0 - R2

! >> for dissipative force

                   Npair = Npair + 1

                   ListIJ(1,Npair) = i
                   ListIJ(2,Npair) = j

                   R1List(Npair) = R1

                   dRList(:,Npair) = Sij

! << for dissipative force

                   aaa =   a(is,TypeNum(j)) 
                   fc  =   aaa * Wc * Wc

                   R5 = - R2 * R2 / 5.d0
                   R3 =   R2 * 2.d0 / 3.d0

                   PotDP = PotDP + aaa * ( R1*( R5 + R3 -1.d0 ) + ccc )

                   Fij = fc * Sij

                   Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                   Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                   Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                   Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                   Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                   Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                   FrcDPt(:,i) = FrcDPt(:,i) + Fij

                   FrcDPt(:,j) = FrcDPt(:,j) - Fij

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     Virial(2,1) = Virial(1,2)
     Virial(3,1) = Virial(1,3)
     Virial(3,2) = Virial(2,3)
!
! #####################################################################
    end if
! #####################################################################

    if(Npair > Ndm) then
      write(*,*) 'ERROR : Npair exceeds its limit number'
      write(*,*) 'Npair = ', Npair, ',  Ndm = ',Ndm
      call Finalize
    end if

end subroutine ForceDPD_Smooth22_Lowe


!#####################################################################
!#####################################################################


subroutine ForceDPD_Smooth23_Lowe

use Configuration, only : R
use CommonDPD
use CellListMethod
use BookParam, only : Npair, ListIJ
use BondedParam, only : Frc_Bond, Vir_Bond, Ene_Bond
use CellParam, only : CellL, InvCL
use ThermoData, only : Virial

implicit none

integer :: icell, i, j
integer :: jcell0, jcell, nabor
integer :: NumI, Is
integer, dimension(4) :: NonFc
real(8), dimension(3) :: Ri, Rij, Sij, Fij
real(8) :: InvR, R1, R2
real(8) :: fc, Slide
real(8) :: ranf, shv
external ranf
! >>
real(8) :: ccc, R7, R5, R3, Wc, aaa
! <<

   call Force_Bond_DP

   PotDP  = Ene_Bond
   Virial = Vir_Bond
   FrcDPt = Frc_Bond

   Npair = 0
   SLList = 0.d0

! >>
   ccc = 16.d0 / 35.d0
! <<

! #####################################################################
   if(QSheared) then
! #####################################################################
!
     shv = ShearRate * CellL(2)

! -------------------------
     call Mapping_TopLine
     call LinkCell
! -------------------------
!
     do icell = 1 , Ncell

       i = Head(icell)

! ** loop over all particles in the cell **

       do while( i > 0 )

         Ri    = R(:,i)
         NumI  = Ncal(i)
         Is    = TypeNum(i)
         NonFc = NpBond(:,i)

! ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

              Rij    = Ri - R(:,j)

              Slide  = -nint( Rij(2) * InvCL(2) )
              Rij(1) = Rij(1) + Slide * SlideGap

              Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
              Rij(2) = Rij(2) + Slide * CellL(2)
              Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

              R2 = dot_product( Rij, Rij )

              if(R2 < 1.0d0) then

                R1   = sqrt(R2)
                InvR = 1.d0 / R1
                Wc   = 1.d0 - R2

                Sij = Rij * InvR

! >> for dissipative force

                Npair = Npair + 1

                ListIJ(1,Npair) = i
                ListIJ(2,Npair) = j

                R1List(Npair) = R1

                dRList(:,Npair) = Sij
                SLList(Npair) =  Slide * shv

! << for dissipative force

                aaa =   a(is,TypeNum(j)) 
                fc =   aaa * Wc * Wc * Wc

                R7 =   R2 * R2 * R2 / 7.d0
                R5 = - R2 * R2 * 3.d0 / 5.d0
                R3 =   R2

                PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 - 1.d0 ) + ccc )

                Fij = fc * Sij

                Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                FrcDPt(:,i) = FrcDPt(:,i) + Fij

                FrcDPt(:,j) = FrcDPt(:,j) - Fij

              end if

            end if

            j = NextP(j)

          end do

!** loop over neighbouring cells **

          jcell0 = 16 * ( icell - 1 )

          do nabor = 1 , 16

            jcell = Map( jcell0 + nabor )

            if(jcell > 0) then

! ** loop over all particles in neighbouring cells **

              j = Head(jcell)

              do while( j > 0 )

                if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
                &                       .and. ( j /= NonFc(2) ) &
                &                       .and. ( j /= NonFc(3) ) &
                &                       .and. ( j /= NonFc(4) ) ) then

                  Rij    = Ri - R(:,j)

                  Slide  = -nint( Rij(2) * InvCL(2) )
                  Rij(1) = Rij(1) + Slide * SlideGap
                  Rij(1) = Rij(1) - nint( Rij(1) * InvCL(1) ) * CellL(1)
                  Rij(2) = Rij(2) + Slide * CellL(2)
                  Rij(3) = Rij(3) - nint( Rij(3) * InvCL(3) ) * CellL(3)

                  R2     = dot_product( Rij, Rij )

                  if( R2 < 1.d0 ) then

                    R1   = sqrt(R2)
                    InvR = 1.d0 / R1

                    Sij = Rij * InvR

                    Wc   = 1.d0 - R2

! >> for dissipative force

                    Npair = Npair + 1

                    ListIJ(1,Npair) = i
                    ListIJ(2,Npair) = j

                    R1List(Npair) = R1

                    dRList(:,Npair) = Sij
                    SLList(Npair) =  Slide * shv

! << for dissipative force

                    aaa =   a(is,TypeNum(j)) 
                    fc  =   aaa * Wc * Wc * Wc

                    R7 =   R2 * R2 * R2 / 7.d0
                    R5 = - R2 * R2 * 3.d0 / 5.d0
                    R3 =   R2

                    PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 - 1.d0 ) + ccc )

                    Fij = fc * Sij

                    Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                    Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                    Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                    Virial(2,1) = Virial(2,1) + Fij(2) * Rij(1)
                    Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                    Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                    Virial(3,1) = Virial(3,1) + Fij(3) * Rij(1)
                    Virial(3,2) = Virial(3,2) + Fij(3) * Rij(2)
                    Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                    FrcDPt(:,i) = FrcDPt(:,i) + Fij

                    FrcDPt(:,j) = FrcDPt(:,j) - Fij

                  end if

                end if

                j = NextP(j)

              end do

            end if

          end do

          i = NextP(i)

        end do

      end do

! #####################################################################
    else
! #####################################################################
!
! ---------------------
     call LinkCell
! ---------------------
!
     do icell = 1 , Ncell

       i = Head(icell)
!
!       ** loop over all particles in the cell **
!
       do while( i > 0 )

         Ri    = R (:,i)
         Is    = TypeNum(i)
         NumI  = Ncal(i)
         NonFc = NpBond(:,i)

!  ** loop over all particles below i in the current cell **

         j = NextP(i)

         do while( j > 0 )

           if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
           &                       .and. ( j /= NonFc(2) ) &
           &                       .and. ( j /= NonFc(3) ) &
           &                       .and. ( j /= NonFc(4) ) ) then

             Rij = Ri - R(:,j)
             Rij = Rij - nint( Rij * InvCL ) * CellL
             R2  = dot_product( Rij , Rij )

             if(R2 < 1.0d0) then

               R1   = sqrt(R2)
               InvR = 1.0d0 / R1

               Sij = Rij * InvR
               Wc   = 1.d0 - R2

! >> for dissipative force

               Npair = Npair + 1

               ListIJ(1,Npair) = i
               ListIJ(2,Npair) = j

               R1List(Npair) = R1

               dRList(:,Npair) = Sij

! << for dissipative force

               aaa =   a(is,TypeNum(j)) 
               fc  =   aaa * Wc * Wc * Wc

               R7 =   R2 * R2 * R2 / 7.d0
               R5 = - R2 * R2 * 3.d0 / 5.d0
               R3 =   R2

               PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 -1.d0 ) + ccc )

               Fij = fc * Sij

               Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
               Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
               Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
               Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
               Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
               Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

               FrcDPt(:,i) = FrcDPt(:,i) + Fij

               FrcDPt(:,j) = FrcDPt(:,j) - Fij

             end if

           end if

           j = NextP(j)

         end do

!** loop over neighbouring cells **

         jcell0 = 16 * ( icell - 1 )

         do nabor = 1 , 16

           jcell = Map( jcell0 + nabor )

           if(jcell > 0) then

!  ** loop over all particles in neighbouring cells **

             j = Head(jcell)

             do while( j /= 0 )

               if( ( NumI /= Ncal(j) ) .and. ( j /= NonFc(1) ) &
               &                       .and. ( j /= NonFc(2) ) &
               &                       .and. ( j /= NonFc(3) ) &
               &                       .and. ( j /= NonFc(4) ) ) then

                 Rij = Ri - R(:,j)
                 Rij = Rij - nint( Rij*InvCL ) * CellL
                 R2  = dot_product( Rij, Rij )

                 if(R2 < 1.d0) then

                   R1   = sqrt(R2)
                   InvR = 1.d0 / R1

                   Sij  = Rij * InvR
                   Wc   = 1.d0 - R2

! >> for dissipative force

                   Npair = Npair + 1

                   ListIJ(1,Npair) = i
                   ListIJ(2,Npair) = j

                   R1List(Npair) = R1

                   dRList(:,Npair) = Sij

! << for dissipative force

                   aaa =   a(is,TypeNum(j)) 
                   fc  =   aaa * Wc * Wc * Wc

                   R7 =   R2 * R2 * R2 / 7.d0
                   R5 = - R2 * R2 * 3.d0 / 5.d0
                   R3 =   R2

                   PotDP = PotDP + aaa * ( R1*( R7 + R5 + R3 -1.d0 ) + ccc )

                   Fij = fc * Sij

                   Virial(1,1) = Virial(1,1) + Fij(1) * Rij(1)
                   Virial(1,2) = Virial(1,2) + Fij(1) * Rij(2)
                   Virial(1,3) = Virial(1,3) + Fij(1) * Rij(3)
                   Virial(2,2) = Virial(2,2) + Fij(2) * Rij(2)
                   Virial(2,3) = Virial(2,3) + Fij(2) * Rij(3)
                   Virial(3,3) = Virial(3,3) + Fij(3) * Rij(3)

                   FrcDPt(:,i) = FrcDPt(:,i) + Fij

                   FrcDPt(:,j) = FrcDPt(:,j) - Fij

                 end if

               end if

               j = NextP(j)

             end do

           end if

         end do

         i = NextP(i)

       end do

     end do

     Virial(2,1) = Virial(1,2)
     Virial(3,1) = Virial(1,3)
     Virial(3,2) = Virial(2,3)
!
! #####################################################################
    end if
! #####################################################################

    if(Npair > Ndm) then
      write(*,*) 'ERROR : Npair exceeds its limit number'
      write(*,*) 'Npair = ', Npair, ',  Ndm = ',Ndm
      call Finalize
    end if

end subroutine ForceDPD_Smooth23_Lowe


!#####################################################################
!#####################################################################


subroutine ForceDPD_Attractive_Lowe

use Configuration, only : R
use CommonDPD
use CellListMethod

implicit none

end subroutine ForceDPD_Attractive_Lowe
