

!######################################################################
!######################################################################


subroutine MD_NPT

use Forplot

implicit none

character(len=80) :: String
integer :: i, j, k
integer :: Ns, TotalStep
#ifdef BMONI
integer :: eofile
#endif

   open(51,file='E_data/Ene_Intra.dat',status='unknown',form='formatted')
   open(52,file='E_data/Ene_Inter.dat',status='unknown',form='formatted')
   open(53,file='E_data/Ene_System.dat',status='unknown',form='formatted')
   open(54,file='E_data/Cell.dat',status='unknown',form='formatted')
   open(55,file='E_data/CellABC.dat',status='unknown',form='formatted')
   open(56,file='E_data/System.dat',status='unknown',form='formatted')
   open(57,file='E_data/Pressure.dat',status='unknown',form='formatted')
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   open(58,file='E_data/Strain.dat',status='unknown',form='formatted')
   end if
#ifdef SURF
   open(59,file='E_data/SurfTension.dat',status='unknown',form='formatted')
#endif

   write(51,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(52,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ         Elec(real)   Elec(recip)  Elec(self)   Electrstatic '
   if((Ensemble=='NPT').or.(Ensemble=='NtT')) then
   write(53,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E    Bath Potentia'
   else if((Ensemble=='NPH').or.(Ensemble=='NtH')) then
   write(53,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E'
   end if
   write(54,'(a/a)') '# Cell matrix h in Angstrom ', &
   & '# Time [ps] h_xx    h_xy    h_xz    h_yx    h_yy    h_yz    h_zx    h_zy    h_zz  '
   write(55,'(a/a)') '# Cell size La Lb Lc angles', &
   & '# Time [ps] La      Lb      Lc      alpha   beta    gamma   '
   write(56,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  H[kcal/mol]  '
   write(57,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   write(58,'(a/a)') '# Strain energy [kcal/mol]', '# Time [ps]  Strain E '
   end if
#ifdef SURF
   write(59,'(a/a)') '# Surface tension in dyn/cm ', &
   & '# Time [ps]   Gamma '
#endif

   if(QAccAv) then

   open(61,file='E_data/Av_Ene_Intra.dat',status='unknown',form='formatted')
   open(62,file='E_data/Av_Ene_Inter.dat',status='unknown',form='formatted')
   open(63,file='E_data/Av_Ene_System.dat',status='unknown',form='formatted')
   open(64,file='E_data/Av_Cell.dat',status='unknown',form='formatted')
   open(65,file='E_data/Av_CellABC.dat',status='unknown',form='formatted')
   open(66,file='E_data/Av_System.dat',status='unknown',form='formatted')
   open(67,file='E_data/Av_Pressure.dat',status='unknown',form='formatted')
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   open(68,file='E_data/Av_Strain.dat',status='unknown',form='formatted')
   end if
#ifdef SURF
   open(69,file='E_data/Av_SurfTension.dat',status='unknown',form='formatted')
#endif

   write(61,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(62,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ         Elec(real)   Elec(recip)  Elec(self)   Electrstatic '
   if((Ensemble=='NPT').or.(Ensemble=='NtT')) then
   write(63,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E    Bath Potentia'
   else if((Ensemble=='NPH').or.(Ensemble=='NtH')) then
   write(63,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E'
   end if
   write(64,'(a/a)') '# Cell matrix h in Angstrom ', &
   & '# Time [ps] h_xx    h_xy    h_xz    h_yx    h_yy    h_yz    h_zx    h_zy    h_zz  '
   write(65,'(a/a)') '# Cell size La Lb Lc angles', &
   & '# Time [ps] La      Lb      Lc      alpha   beta    gamma   '
   write(66,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  H[kcal/mol]  '
   write(67,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   write(68,'(a/a)') '# Strain energy [kcal/mol]', '# Time [ps]  Strain E '
   end if
#ifdef SURF
   write(69,'(a/a)') '# Surface tension in dyn/cm ', &
   & '# Time [ps]   Gamma '
#endif

   Sum_Temp       =  0.d0
   Sum_E_BondMol  =  0.d0
   Sum_E_AngleMol =  0.d0
   Sum_E_UBMol    =  0.d0
   Sum_E_DihedMol =  0.d0
   Sum_E_ImproMol =  0.d0
   Sum_E_OptC     =  0.d0
   Sum_E_LJMol    =  0.d0
   Sum_E_ErMol    =  0.d0
   Sum_E_EkMol    =  0.d0
   Sum_E_EsMol    =  0.d0
   Sum_E_ElecMol  =  0.d0
   Sum_E_Int_Mol  =  0.d0
   Sum_E_Ext_Mol  =  0.d0
   Sum_Potential  =  0.d0
   Sum_Kinetic    =  0.d0
   Sum_EneSystem  =  0.d0
   Sum_ThermoBath =  0.d0
   Sum_Pot_p      =  0.d0
   Sum_Hamiltonian=  0.d0
   Sum_Area       =  0.d0
   Sum_Volume     =  0.d0
   Sum_density    =  0.d0
   Sum_LenA       =  0.d0
   Sum_LenB       =  0.d0
   Sum_LenC       =  0.d0
   Sum_AngBC      =  0.d0
   Sum_AngCA      =  0.d0
   Sum_AngAB      =  0.d0
   Sum_H          =  0.d0
   Sum_Pressure   =  0.d0
#ifdef SURF
   Sum_SurfTens = 0.d0
#endif

   end if

   TotalStep = 0

   do i = 1 , Nfile

#ifdef BMONI
     open(11,file=trim(Energy_File(i)),status='old',form='unformatted')
#else
     open(11,file=trim(Energy_File(i)),status='old',form='formatted')
#endif
     write(*,*) 'File ',trim(Energy_File(i)),' was opened.'

#ifdef BMONI

     Ns = 0

     do

       read(11,iostat=eofile) Timeps, Temp
       if(eofile /= 0) exit

       Ns = Ns + 1
       TotalStep = TotalStep + 1

       read(11)  E_BondMol,  E_AngleMol, E_UBMol,         &
       &         E_DihedMol, E_ImproMol, E_OptC,          &
       &         E_LJMol, E_ErMol,                        &
       &         E_EkMol, E_EsMol,   E_ElecMol
       if((Ensemble=='NPT').or.(Ensemble=='NtT')) then
       read(11)  E_Int_Mol, E_Ext_Mol, Potential,         &
       &         Kinetic, EneSystem, ThermoBath
       else if((Ensemble=='NPH').or.(Ensemble=='NtH')) then
       read(11)  E_Int_Mol, E_Ext_Mol, Potential,         &
       &         Kinetic, EneSystem
       end if
       if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
       read(11) Pot_p
       end if
       read(11) Hamiltonian, Area, Volume, density,  &
       &        LenA, LenB, LenC, AngBC, AngCA, AngAB

       read(11) H, Pressure

#else

     do

       read(11,'(a80)') String

       if(String(3:56)=='#############< Start MD time evolution >##############') then

         read(11,*)

         Ns = 0

         do

           read(11,'(a80)') String

           if(String(1:31)=='-------------------------------') exit

           Ns = Ns + 1

           TotalStep = TotalStep + 1

           read(String,'(f12.4,f10.2)') &
           &          Timeps, Temp

           read(11,'(6d13.5/5d13.5)') &
           &         E_BondMol,  E_AngleMol, E_UBMol,         &
           &         E_DihedMol, E_ImproMol, E_OptC,          &
           &         E_LJMol, E_ErMol,                        &
           &         E_EkMol, E_EsMol,   E_ElecMol
           if((Ensemble=='NPT').or.(Ensemble=='NtT')) then
           read(11,'(6d13.5)') &
           &         E_Int_Mol, E_Ext_Mol, Potential,         &
           &         Kinetic, EneSystem, ThermoBath
           else if((Ensemble=='NPH').or.(Ensemble=='NtH')) then
           read(11,'(5d13.5)') &
           &         E_Int_Mol, E_Ext_Mol, Potential,         &
           &         Kinetic, EneSystem
           end if
           if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
           read(11,'(d13.5)') Pot_p
           end if
           read(11,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')  &
           &         Hamiltonian, Area, Volume, density,  &
           &         LenA, LenB, LenC, AngBC, AngCA, AngAB

           read(11,'(3(3f8.3,5x,3f12.4/))')                   &
           &         (H(1,j),j=1,3) , (Pressure(1,j),j=1,3),  &
           &         (H(2,j),j=1,3) , (Pressure(2,j),j=1,3),  &
           &         (H(3,j),j=1,3) , (Pressure(3,j),j=1,3)

#endif

#ifdef SURF
           SurfTens = 0.5d0 * H(3,3) * (Pressure(3,3) - 0.5d0*(Pressure(1,1)+Pressure(2,2))) &
           &        * 0.1d0 ! unit [dyne/cm] or [mN/m]
#endif

           if(QAccAv.and.(TotalStep>Nskip)) call Average_MDNPT

           if(mod(Ns,Nst)==0) call Out_data_MDNPT

#ifdef BMONI

     end do

#else

         end do

       end if

       if(String(1:31)==' MD simulation ended normally! ') exit

     end do

#endif

     close(11)

   end do

   close(51)
   close(52)
   close(53)
   close(54)
   close(55)
   close(56)
   close(57)
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   close(58)
   end if
#ifdef SURF
   close(59)
#endif

   if(QAccAv) then
   close(61)
   close(62)
   close(63)
   close(64)
   close(65)
   close(66)
   close(67)
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   close(68)
   end if
#ifdef SURF
   close(69)
#endif
   end if

#ifdef COMPRESS
   if(QAccAv) call Calc_Compressibility
#endif

Contains

   subroutine Out_data_MDNPT

   implicit none

   real(8) :: xinv

     write(51,'(f11.3,6d13.5)')                  &
     &  Timeps, E_BondMol,  E_AngleMol, E_UBMol, &
     &          E_DihedMol, E_ImproMol, E_OptC

     write(52,'(f11.3,5d13.5)')                  &
     &  Timeps, E_LJMol, E_ErMol,                &
     &          E_EkMol, E_EsMol,   E_ElecMol

     if((Ensemble=='NPT').or.(Ensemble=='NtT')) then
     write(53,'(f11.3,6d13.5)')                    &
     &  Timeps, E_Int_Mol, E_Ext_Mol, Potential,   &
     &          Kinetic, EneSystem, ThermoBath
     else if((Ensemble=='NPH').or.(Ensemble=='NtH')) then
     write(53,'(f11.3,5d13.5)')                    &
     &  Timeps, E_Int_Mol, E_Ext_Mol, Potential,   &
     &          Kinetic, EneSystem
     end if

        write(56,'(f11.3,f15.5,3d15.8,f10.6,d16.8)')  &
     &   Timeps, Temp, Area, Volume/Area, Volume, density, Hamiltonian

        write(54,'(f11.3,9f9.2)') &
     &   Timeps, ((H(j,k),k=1,3),j=1,3)

     write(55,'(f11.3,6f9.2)') &
     &   Timeps, LenA, LenB, LenC, AngBC, AngCA, AngAB

     write(57,'(f11.3,9d15.7)') &
     &   Timeps, ((Pressure(j,k),k=1,3),j=1,3)

     if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
     write(58,'(f11.3,d13.5)') Timeps, Pot_p
     end if

#ifdef SURF
     write(59,'(f11.3,f11.3)') Timeps, SurfTens
#endif

     if(QAccAv.and.TotalStep>Nskip) then

       xinv = 1.d0 / dble(TotalStep-Nskip)

       Av_Temp       =  Sum_Temp       * xinv
       Av_E_BondMol  =  Sum_E_BondMol  * xinv
       Av_E_AngleMol =  Sum_E_AngleMol * xinv
       Av_E_UBMol    =  Sum_E_UBMol    * xinv
       Av_E_DihedMol =  Sum_E_DihedMol * xinv
       Av_E_ImproMol =  Sum_E_ImproMol * xinv
       Av_E_OptC     =  Sum_E_OptC     * xinv
       Av_E_LJMol    =  Sum_E_LJMol    * xinv
       Av_E_ErMol    =  Sum_E_ErMol    * xinv
       Av_E_EkMol    =  Sum_E_EkMol    * xinv
       Av_E_EsMol    =  Sum_E_EsMol    * xinv
       Av_E_ElecMol  =  Sum_E_ElecMol  * xinv
       Av_E_Int_Mol  =  Sum_E_Int_Mol  * xinv
       Av_E_Ext_Mol  =  Sum_E_Ext_Mol  * xinv
       Av_Potential  =  Sum_Potential  * xinv
       Av_Kinetic    =  Sum_Kinetic    * xinv
       Av_EneSystem  =  Sum_EneSystem  * xinv
       Av_ThermoBath =  Sum_ThermoBath * xinv
       Av_Pot_p      =  Sum_Pot_p      * xinv
       Av_Hamiltonian=  Sum_Hamiltonian* xinv
       Av_Area       =  Sum_Area       * xinv
       Av_Volume     =  Sum_Volume     * xinv
       Av_density    =  Sum_density    * xinv
       Av_LenA       =  Sum_LenA       * xinv
       Av_LenB       =  Sum_LenB       * xinv
       Av_LenC       =  Sum_LenC       * xinv
       Av_AngBC      =  Sum_AngBC      * xinv
       Av_AngCA      =  Sum_AngCA      * xinv
       Av_AngAB      =  Sum_AngAB      * xinv
       Av_H          =  Sum_H          * xinv
       Av_Pressure   =  Sum_Pressure   * xinv
#ifdef SURF
       Av_SurfTens   =  Sum_SurfTens   * xinv
#endif

       write(61,'(f11.3,6d13.5)')                  &
       &  Timeps, Av_E_BondMol,  Av_E_AngleMol, Av_E_UBMol, &
       &          Av_E_DihedMol, Av_E_ImproMol, Av_E_OptC

       write(62,'(f11.3,5d13.5)')                  &
       &  Timeps, Av_E_LJMol, Av_E_ErMol,                &
       &          Av_E_EkMol, Av_E_EsMol,   Av_E_ElecMol

       if((Ensemble=='NPT').or.(Ensemble=='NtT')) then
       write(63,'(f11.3,6d13.5)')                    &
       &  Timeps, Av_E_Int_Mol, Av_E_Ext_Mol, Av_Potential,   &
       &          Av_Kinetic, Av_EneSystem, Av_ThermoBath
       else if((Ensemble=='NPH').or.(Ensemble=='NtH')) then
       write(63,'(f11.3,5d13.5)')                    &
       &  Timeps, Av_E_Int_Mol, Av_E_Ext_Mol, Av_Potential,   &
       &          Av_Kinetic, Av_EneSystem
       end if

          write(66,'(f11.3,f15.5,3d15.7,f10.6,d16.8)')  &
       &   Timeps, Av_Temp, Av_Area, Av_Volume/Av_Area, Av_Volume, Av_density, Av_Hamiltonian

          write(64,'(f11.3,9f9.2)') &
       &   Timeps, ((Av_H(j,k),k=1,3),j=1,3)

       write(65,'(f11.3,6f9.2)') &
       &   Timeps, Av_LenA, Av_LenB, Av_LenC, Av_AngBC, Av_AngCA, Av_AngAB

       write(67,'(f11.3,9d15.7)') &
       &   Timeps, ((Av_Pressure(j,k),k=1,3),j=1,3)

       if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
       write(68,'(f11.3,d13.5)') Timeps, Av_Pot_p
       end if
#ifdef SURF
       write(69,'(f11.3,f11.3)') Timeps, Av_SurfTens
#endif

     end if

   end subroutine Out_data_MDNPT

   subroutine Average_MDNPT

   implicit none

      Sum_Temp       =  Sum_Temp       + Temp

      Sum_E_BondMol  =  Sum_E_BondMol  + E_BondMol
      Sum_E_AngleMol =  Sum_E_AngleMol + E_AngleMol
      Sum_E_UBMol    =  Sum_E_UBMol    + E_UBMol
      Sum_E_DihedMol =  Sum_E_DihedMol + E_DihedMol
      Sum_E_ImproMol =  Sum_E_ImproMol + E_ImproMol
      Sum_E_OptC     =  Sum_E_OptC     + E_OptC
      Sum_E_LJMol    =  Sum_E_LJMol    + E_LJMol
      Sum_E_ErMol    =  Sum_E_ErMol    + E_ErMol
      Sum_E_EkMol    =  Sum_E_EkMol    + E_EkMol
      Sum_E_EsMol    =  Sum_E_EsMol    + E_EsMol
      Sum_E_ElecMol  =  Sum_E_ElecMol  + E_ElecMol

      Sum_E_Int_Mol  =  Sum_E_Int_Mol  + E_Int_Mol
      Sum_E_Ext_Mol  =  Sum_E_Ext_Mol  + E_Ext_Mol
      Sum_Potential  =  Sum_Potential  + Potential
      Sum_Kinetic    =  Sum_Kinetic    + Kinetic
      Sum_EneSystem  =  Sum_EneSystem  + EneSystem
      Sum_Hamiltonian=  Sum_Hamiltonian+ Hamiltonian
      Sum_Area       =  Sum_Area       + Area
      Sum_Volume     =  Sum_Volume     + Volume
      Sum_density    =  Sum_density    + density
      Sum_LenA       =  Sum_LenA       + LenA
      Sum_LenB       =  Sum_LenB       + LenB
      Sum_LenC       =  Sum_LenC       + LenC
      Sum_AngBC      =  Sum_AngBC      + AngBC
      Sum_AngCA      =  Sum_AngCA      + AngCA
      Sum_AngAB      =  Sum_AngAB      + AngAB
      Sum_H          =  Sum_H          + H
      Sum_Pressure   =  Sum_Pressure   + Pressure

#ifdef SURF
      Sum_SurfTens   =  Sum_SurfTens   + SurfTens
#endif

      if((Ensemble=='NPT').or.(Ensemble=='NtT')) then
        Sum_ThermoBath =  Sum_ThermoBath + ThermoBath
      end if

      if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
        Sum_Pot_p      =  Sum_Pot_p      + Pot_p
      end if

   end subroutine Average_MDNPT

#ifdef COMPRESS
   subroutine Calc_Compressibility

   implicit none

   real(8), parameter :: bolt_const = 1.380662d-23
   real(8) :: d, dV, dV2, DeltaV2, sDeltaV2, RT, Cbulk

      RT = bolt_const * Av_Temp

      open(56,file='E_data/System.dat',status='old',form='formatted')
      open(60,file='E_data/Av_Compressibility.dat',status='unknown',form='formatted')

      read(56,*)
      read(56,*)

      write(60,'(a/a)') '# Bulk compressibility in GPa^-1 ','# Time [ps]   K_bulk '

      sDeltaV2 = 0.d0

      do i = Nst, TotalStep, Nst
        read(56,'(f11.3,f15.5,3f15.8,f10.6,d16.8)')  &
        &   Timeps, Temp, Area, d, Volume, density, Hamiltonian
        if(i>Nskip) then
          dV  = Volume - Av_Volume
          dV2 = dV * dV
          sDeltaV2 = sDeltaV2 + dV2
          DeltaV2 = sDeltaV2 / dble(i-Nskip)
          Cbulk = DeltaV2*1.d+09/Av_Volume/RT*1.d-30
          write(60,'(f11.3,e15.6)') Timeps, Cbulk
        end if
      end do

      close(56)
      close(60)

   end subroutine Calc_Compressibility
#endif

end subroutine MD_NPT


!######################################################################
!######################################################################


subroutine MD_NVT

use Forplot

implicit none

integer :: i,j,k
integer :: Ns, TotalStep
character(len=80) :: String
#ifdef BMONI
integer :: eofile
#endif

   open(51,file='E_data/Ene_Intra.dat',status='unknown',form='formatted')
   open(52,file='E_data/Ene_Inter.dat',status='unknown',form='formatted')
   open(53,file='E_data/Ene_System.dat',status='unknown',form='formatted')
   open(56,file='E_data/System.dat',status='unknown',form='formatted')
   open(57,file='E_data/Pressure.dat',status='unknown',form='formatted')
#ifdef SURF
   open(58,file='E_data/SurfTension.dat',status='unknown',form='formatted')
#endif

   write(51,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(52,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ         Elec(real)   Elec(recip)  Elec(self)   Electrstatic '
   if(Ensemble=='NVT') then
   write(53,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E    Bath Potentia'
   else if(Ensemble=='NVE') then
   write(53,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E '
   end if
    write(56,'(a/a)') '# Temperature & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]  H[kcal/mol]  '
   write(57,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'
#ifdef SURF
   write(58,'(a/a)') '# Surface tension in dyn/cm ', &
   & '# Time [ps]   Gamma '
#endif

   if(QAccAv) then

   open(61,file='E_data/Av_Ene_Intra.dat',status='unknown',form='formatted')
   open(62,file='E_data/Av_Ene_Inter.dat',status='unknown',form='formatted')
   open(63,file='E_data/Av_Ene_System.dat',status='unknown',form='formatted')
   open(66,file='E_data/Av_System.dat',status='unknown',form='formatted')
   open(67,file='E_data/Av_Pressure.dat',status='unknown',form='formatted')
#ifdef SURF
   open(68,file='E_data/Av_SurfTension.dat',status='unknown',form='formatted')
#endif

   write(61,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(62,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ         Elec(real)   Elec(recip)  Elec(self)   Electrstatic '
   if(Ensemble=='NVT') then
   write(63,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E    Bath Potentia'
   else if(Ensemble=='NVE') then
   write(63,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E '
   end if
    write(66,'(a/a)') '# Temperature & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]  H[kcal/mol]  '
   write(67,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'
#ifdef SURF
   write(68,'(a/a)') '# Surface tension in dyn/cm ', &
   & '# Time [ps]   Gamma '
#endif

   Sum_Temp       =  0.d0
   Sum_E_BondMol  =  0.d0
   Sum_E_AngleMol =  0.d0
   Sum_E_UBMol    =  0.d0
   Sum_E_DihedMol =  0.d0
   Sum_E_ImproMol =  0.d0
   Sum_E_OptC     =  0.d0
   Sum_E_LJMol    =  0.d0
   Sum_E_ErMol    =  0.d0
   Sum_E_EkMol    =  0.d0
   Sum_E_EsMol    =  0.d0
   Sum_E_ElecMol  =  0.d0
   Sum_E_Int_Mol  =  0.d0
   Sum_E_Ext_Mol  =  0.d0
   Sum_Potential  =  0.d0
   Sum_Kinetic    =  0.d0
   Sum_EneSystem  =  0.d0
   Sum_ThermoBath =  0.d0
   Sum_Pot_p      =  0.d0
   Sum_Hamiltonian=  0.d0
   Sum_Area       =  0.d0
   Sum_Volume     =  0.d0
   Sum_density    =  0.d0
   Sum_LenA       =  0.d0
   Sum_LenB       =  0.d0
   Sum_LenC       =  0.d0
   Sum_AngBC      =  0.d0
   Sum_AngCA      =  0.d0
   Sum_AngAB      =  0.d0
   Sum_H          =  0.d0
   Sum_Pressure   =  0.d0
#ifdef SURF
   Sum_SurfTens = 0.d0
#endif

   end if

   TotalStep = 0

   do i = 1 , Nfile

#ifdef BMONI
     open(11,file=trim(Energy_File(i)),status='old',form='unformatted')
#else
     open(11,file=trim(Energy_File(i)),status='old',form='formatted')
#endif
     write(*,*) 'File ',trim(Energy_File(i)),' was opened.'

#ifdef BMONI

     Ns = 0

     do

       read(11,iostat=eofile) Timeps, Temp
       if(eofile /= 0) exit

       Ns = Ns + 1
       TotalStep = TotalStep + 1

       read(11)   E_BondMol, E_AngleMol, E_UBMol,     &
       &          E_DihedMol,E_ImproMol, E_OptC,      &
       &          E_LJMol, E_ErMol,                   &
       &          E_EkMol, E_EsMol,   E_ElecMol
       if(Ensemble=='NVT') then
       read(11)   E_Int_Mol, E_Ext_Mol, Potential,    &
       &          Kinetic, EneSystem, ThermoBath
       else if(Ensemble=='NVE') then
       read(11)   E_Int_Mol, E_Ext_Mol, Potential,    &
       &          Kinetic, EneSystem
       end if
       read(11)   Hamiltonian, Pressure

#else

     do

       Ns = 0

       read(11,'(a80)') String

       if(String(3:56)=='#############< Start MD time evolution >##############') then

         read(11,*)

         Ns = 0

         do

           read(11,'(a80)') String

           if(String(1:31)=='-------------------------------') exit

           Ns = Ns + 1

           TotalStep = TotalStep + 1

           read(String,'(f12.4,f10.2)') &
           &          Timeps, Temp

           read(11,'(6d13.5/5d13.5)') &
           &          E_BondMol, E_AngleMol, E_UBMol,     &
           &          E_DihedMol,E_ImproMol, E_OptC,      &
           &          E_LJMol, E_ErMol,                   &
           &          E_EkMol, E_EsMol,   E_ElecMol
           if(Ensemble=='NVT') then
           read(11,'(6d13.5)') &
           &          E_Int_Mol, E_Ext_Mol, Potential,    &
           &          Kinetic, EneSystem, ThermoBath
           else if(Ensemble=='NVE') then
           read(11,'(5d13.5)') &
           &          E_Int_Mol, E_Ext_Mol, Potential,    &
           &          Kinetic, EneSystem
           end if
           read(11,'(d16.8/3(3f12.4/))') &
           &          Hamiltonian,                        &
           &          ( Pressure(1,j) , j = 1 , 3 ),      &
           &          ( Pressure(2,j) , j = 1 , 3 ),      &
           &          ( Pressure(3,j) , j = 1 , 3 )
#endif

#ifdef SURF
           SurfTens = 0.5d0 * H(3,3) * (Pressure(3,3) - 0.5d0*(Pressure(1,1)+Pressure(2,2))) &
           &        * 0.1d0 ! unit [dyne/cm] or [mN/m]
#endif

           if(QAccAv.and.TotalStep>Nskip) call Average_MDNVT

           if(mod(Ns,Nst)==0) call Out_data_MDNVT

#ifdef BMONI

     end do

#else
         end do

       end if

       if(String(1:31)==' MD simulation ended normally! ') exit

     end do

#endif

     close(11)

   end do

Contains

   subroutine Out_data_MDNVT

   implicit none

   real(8) :: xinv

     write(51,'(f11.3,6d13.5)')                  &
     &  Timeps, E_BondMol,  E_AngleMol, E_UBMol, &
     &          E_DihedMol, E_ImproMol, E_OptC

     write(52,'(f11.3,5d13.5)')                  &
     &  Timeps, E_LJMol, E_ErMol,                &
     &          E_EkMol, E_EsMol,   E_ElecMol

     if(Ensemble=='NVT') then
     write(53,'(f11.3,6d13.5)')                    &
     &  Timeps, E_Int_Mol, E_Ext_Mol, Potential,   &
     &          Kinetic, EneSystem, ThermoBath
     else if(Ensemble=='NVE') then
     write(53,'(f11.3,5d13.5)')                    &
     &  Timeps, E_Int_Mol, E_Ext_Mol, Potential,   &
     &          Kinetic, EneSystem
     end if

     write(56,'(f11.3,f10.2,d16.8)')  &
     &   Timeps, Temp, Hamiltonian

     write(57,'(f11.3,9d15.7)') &
     &   Timeps, ((Pressure(j,k),k=1,3),j=1,3)

#ifdef SURF
     write(58,'(f11.3,f11.3)') Timeps, SurfTens
#endif

     if(QAccAv.and.TotalStep>Nskip) then

       xinv = 1.d0 / dble(TotalStep-Nskip)

       Av_Temp       =  Sum_Temp       * xinv
       Av_E_BondMol  =  Sum_E_BondMol  * xinv
       Av_E_AngleMol =  Sum_E_AngleMol * xinv
       Av_E_UBMol    =  Sum_E_UBMol    * xinv
       Av_E_DihedMol =  Sum_E_DihedMol * xinv
       Av_E_ImproMol =  Sum_E_ImproMol * xinv
       Av_E_OptC     =  Sum_E_OptC     * xinv
       Av_E_LJMol    =  Sum_E_LJMol    * xinv
       Av_E_ErMol    =  Sum_E_ErMol    * xinv
       Av_E_EkMol    =  Sum_E_EkMol    * xinv
       Av_E_EsMol    =  Sum_E_EsMol    * xinv
       Av_E_ElecMol  =  Sum_E_ElecMol  * xinv
       Av_E_Int_Mol  =  Sum_E_Int_Mol  * xinv
       Av_E_Ext_Mol  =  Sum_E_Ext_Mol  * xinv
       Av_Potential  =  Sum_Potential  * xinv
       Av_Kinetic    =  Sum_Kinetic    * xinv
       Av_EneSystem  =  Sum_EneSystem  * xinv
       Av_ThermoBath =  Sum_ThermoBath * xinv
       Av_Pot_p      =  Sum_Pot_p      * xinv
       Av_Hamiltonian=  Sum_Hamiltonian* xinv
       Av_Pressure   =  Sum_Pressure   * xinv
#ifdef SURF
       Av_SurfTens   =  Sum_SurfTens   * xinv
#endif

       write(61,'(f11.3,6d13.5)')                  &
       &  Timeps, Av_E_BondMol,  Av_E_AngleMol, Av_E_UBMol, &
       &          Av_E_DihedMol, Av_E_ImproMol, Av_E_OptC

       write(62,'(f11.3,5d13.5)')                  &
       &  Timeps, Av_E_LJMol, Av_E_ErMol,                &
       &          Av_E_EkMol, Av_E_EsMol,   Av_E_ElecMol

       if(Ensemble=='NVT') then
       write(63,'(f11.3,6d13.5)')                    &
       &  Timeps, Av_E_Int_Mol, Av_E_Ext_Mol, Av_Potential,   &
       &          Av_Kinetic, Av_EneSystem, Av_ThermoBath
       else if(Ensemble=='NVE') then
       write(63,'(f11.3,5d13.5)')                    &
       &  Timeps, Av_E_Int_Mol, Av_E_Ext_Mol, Av_Potential,   &
       &          Av_Kinetic, Av_EneSystem
       end if

       write(66,'(f11.3,f10.2,d16.8)')  &
       &   Timeps, Av_Temp, Av_Hamiltonian

       write(67,'(f11.3,9d15.7)') &
       &   Timeps, ((Av_Pressure(j,k),k=1,3),j=1,3)

#ifdef SURF
       write(68,'(f11.3,f11.3)') Timeps, Av_SurfTens
#endif
     end if

   end subroutine Out_data_MDNVT

   subroutine Average_MDNVT

   implicit none

      Sum_Temp       =  Sum_Temp       + Temp

      Sum_E_BondMol  =  Sum_E_BondMol  + E_BondMol
      Sum_E_AngleMol =  Sum_E_AngleMol + E_AngleMol
      Sum_E_UBMol    =  Sum_E_UBMol    + E_UBMol
      Sum_E_DihedMol =  Sum_E_DihedMol + E_DihedMol
      Sum_E_ImproMol =  Sum_E_ImproMol + E_ImproMol
      Sum_E_OptC     =  Sum_E_OptC     + E_OptC
      Sum_E_LJMol    =  Sum_E_LJMol    + E_LJMol
      Sum_E_ErMol    =  Sum_E_ErMol    + E_ErMol
      Sum_E_EkMol    =  Sum_E_EkMol    + E_EkMol
      Sum_E_EsMol    =  Sum_E_EsMol    + E_EsMol
      Sum_E_ElecMol  =  Sum_E_ElecMol  + E_ElecMol

      Sum_E_Int_Mol  =  Sum_E_Int_Mol  + E_Int_Mol
      Sum_E_Ext_Mol  =  Sum_E_Ext_Mol  + E_Ext_Mol
      Sum_Potential  =  Sum_Potential  + Potential
      Sum_Kinetic    =  Sum_Kinetic    + Kinetic
      Sum_EneSystem  =  Sum_EneSystem  + EneSystem
      Sum_Hamiltonian=  Sum_Hamiltonian+ Hamiltonian

      Sum_Pressure   =  Sum_Pressure   + Pressure

#ifdef SURF
      Sum_SurfTens   =  Sum_SurfTens   + SurfTens
#endif

      if(Ensemble=='NVT') then
        Sum_ThermoBath =  Sum_ThermoBath + ThermoBath
      end if

   end subroutine Average_MDNVT

end subroutine MD_NVT


!######################################################################
!######################################################################


subroutine MD_iso

use Forplot

implicit none

integer :: i
integer :: Ns, TotalStep
character(len=80) :: String
#ifdef BMONI
integer :: eofile
#endif

   open(51,file='E_data/Ene_Intra.dat',status='unknown',form='formatted')
   open(52,file='E_data/Ene_Inter.dat',status='unknown',form='formatted')
   open(53,file='E_data/Ene_System.dat',status='unknown',form='formatted')
   open(56,file='E_data/System.dat',status='unknown',form='formatted')

   write(51,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(52,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ       Electrstatic '
   if(Ensemble=='NT') then
   write(53,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E    Bath Potentia'
   else if(Ensemble=='NE') then
   write(53,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E '
   end if
    write(56,'(a/a)') '# Temperature & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]  H[kcal/mol]  '

   if(QAccAv) then

   open(61,file='E_data/Av_Ene_Intra.dat',status='unknown',form='formatted')
   open(62,file='E_data/Av_Ene_Inter.dat',status='unknown',form='formatted')
   open(63,file='E_data/Av_Ene_System.dat',status='unknown',form='formatted')
   open(66,file='E_data/Av_System.dat',status='unknown',form='formatted')

   write(61,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(62,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ       Electrstatic '
   if(Ensemble=='NT') then
   write(63,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E    Bath Potentia'
   else if(Ensemble=='NE') then
   write(63,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E '
   end if
    write(66,'(a/a)') '# Temperature & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]  H[kcal/mol]  '

   Sum_Temp       =  0.d0
   Sum_E_BondMol  =  0.d0
   Sum_E_AngleMol =  0.d0
   Sum_E_UBMol    =  0.d0
   Sum_E_DihedMol =  0.d0
   Sum_E_ImproMol =  0.d0
   Sum_E_OptC     =  0.d0
   Sum_E_LJMol    =  0.d0
   Sum_E_ErMol    =  0.d0
   Sum_E_EkMol    =  0.d0
   Sum_E_EsMol    =  0.d0
   Sum_E_ElecMol  =  0.d0
   Sum_E_Int_Mol  =  0.d0
   Sum_E_Ext_Mol  =  0.d0
   Sum_Potential  =  0.d0
   Sum_Kinetic    =  0.d0
   Sum_EneSystem  =  0.d0
   Sum_ThermoBath =  0.d0
   Sum_Pot_p      =  0.d0
   Sum_Hamiltonian=  0.d0
   Sum_Area       =  0.d0
   Sum_Volume     =  0.d0
   Sum_density    =  0.d0
   Sum_LenA       =  0.d0
   Sum_LenB       =  0.d0
   Sum_LenC       =  0.d0
   Sum_AngBC      =  0.d0
   Sum_AngCA      =  0.d0
   Sum_AngAB      =  0.d0
   Sum_H          =  0.d0
   Sum_Pressure   =  0.d0

   end if

   TotalStep = 0

   do i = 1 , Nfile

#ifdef BMONI
     open(11,file=trim(Energy_File(i)),status='old',form='unformatted')
#else
     open(11,file=trim(Energy_File(i)),status='old',form='formatted')
#endif
     write(*,*) 'File ',trim(Energy_File(i)),' was opened.'

#ifdef BMONI
     Ns = 0

     do

       read(11,iostat=eofile) Timeps, Temp
       if(eofile /= 0) exit

       Ns = Ns + 1
       TotalStep = TotalStep + 1

       read(11)  E_BondMol, E_AngleMol, E_UBMol, &
       &         E_DihedMol,E_ImproMol,          &
       &         E_LJMol, E_ErMol, E_OptC,       &
       &         E_Int_Mol, E_Ext_Mol
       if(Ensemble(1:2)=='NT') then
       read(11) Potential, Kinetic, EneSystem, ThermoBath, Hamiltonian
       else if(Ensemble(1:2)=='NE') then
       read(11) Potential, Kinetic, EneSystem, Hamiltonian
       end if

#else

     do

       read(11,'(a80)') String

       if(String(3:56)=='#############< Start MD time evolution >##############') then

         read(11,*)

         Ns = 0

         do

           read(11,'(a80)') String

           if(String(1:31)=='-------------------------------') exit

           Ns = Ns + 1

           TotalStep = TotalStep + 1

           read(String,'(f12.4,f10.2)') Timeps, Temp

           read(11,'(5d13.5/5d13.5)') &
           &          E_BondMol, E_AngleMol, E_UBMol, &
           &          E_DihedMol,E_ImproMol,          &
           &          E_LJMol, E_ErMol, E_OptC,       &
           &          E_Int_Mol, E_Ext_Mol
           if(Ensemble(1:2)=='NT') then
           read(11,'(4d13.5,d16.8/)') Potential, Kinetic, EneSystem, ThermoBath, Hamiltonian
           else if(Ensemble(1:2)=='NE') then
           read(11,'(3d13.5,d16.8/)') Potential, Kinetic, EneSystem, Hamiltonian
           end if

#endif

           if(QAccAv.and.TotalStep>Nskip) call Average_MDiso

           if(mod(Ns,Nst)==0) call Out_data_MDiso

#ifndef BMONI

         end do

       end if

       if(String(1:31)==' MD simulation ended normally! ') exit

#endif

     end do

     close(11)

   end do

Contains

   subroutine Out_data_MDiso

   implicit none
   real(8) :: xinv

     write(51,'(f11.3,6d13.5)')                  &
     &  Timeps, E_BondMol,  E_AngleMol, E_UBMol, &
     &          E_DihedMol, E_ImproMol, E_OptC

     write(52,'(f11.3,2d13.5)')                  &
     &  Timeps, E_LJMol, E_ErMol

     if(Ensemble(1:2)=='NT') then
     write(53,'(f11.3,6d13.5)')                    &
     &  Timeps, E_Int_Mol, E_Ext_Mol, Potential,   &
     &          Kinetic, EneSystem, ThermoBath
     else if(Ensemble(1:2)=='NE') then
     write(53,'(f11.3,5d13.5)')                    &
     &  Timeps, E_Int_Mol, E_Ext_Mol, Potential,   &
     &          Kinetic, EneSystem
     end if

     write(56,'(f11.3,f10.2,d16.8)')  &
     &   Timeps, Temp, Hamiltonian

     if(QAccAv.and.TotalStep>Nskip) then

       xinv = 1.d0 / dble(TotalStep-Nskip)

       Av_Temp       =  Sum_Temp       * xinv
       Av_E_BondMol  =  Sum_E_BondMol  * xinv
       Av_E_AngleMol =  Sum_E_AngleMol * xinv
       Av_E_UBMol    =  Sum_E_UBMol    * xinv
       Av_E_DihedMol =  Sum_E_DihedMol * xinv
       Av_E_ImproMol =  Sum_E_ImproMol * xinv
       Av_E_OptC     =  Sum_E_OptC     * xinv

       Av_E_LJMol    =  Sum_E_LJMol    * xinv
       Av_E_ErMol    =  Sum_E_ErMol    * xinv

       Av_E_Int_Mol  =  Sum_E_Int_Mol  * xinv
       Av_E_Ext_Mol  =  Sum_E_Ext_Mol  * xinv
       Av_Potential  =  Sum_Potential  * xinv
       Av_Kinetic    =  Sum_Kinetic    * xinv
       Av_EneSystem  =  Sum_EneSystem  * xinv
       Av_ThermoBath =  Sum_ThermoBath * xinv
       Av_Pot_p      =  Sum_Pot_p      * xinv

       Av_Hamiltonian=  Sum_Hamiltonian* xinv

       write(61,'(f11.3,6d13.5)')                  &
       &  Timeps, Av_E_BondMol,  Av_E_AngleMol, Av_E_UBMol, &
       &          Av_E_DihedMol, Av_E_ImproMol, Av_E_OptC

       write(62,'(f11.3,2d13.5)')                  &
       &  Timeps, Av_E_LJMol, Av_E_ErMol

       if(Ensemble(1:2)=='NT') then
       write(63,'(f11.3,6d13.5)')                    &
       &  Timeps, Av_E_Int_Mol, Av_E_Ext_Mol, Av_Potential,   &
       &          Av_Kinetic, Av_EneSystem, Av_ThermoBath
       else if(Ensemble(1:2)=='NE') then
       write(63,'(f11.3,5d13.5)')                    &
       &  Timeps, Av_E_Int_Mol, Av_E_Ext_Mol, Av_Potential,   &
       &          Av_Kinetic, Av_EneSystem
       end if

       write(66,'(f11.3,f10.2,d16.8)')  &
       &   Timeps, Av_Temp, Av_Hamiltonian

     end if

   end subroutine Out_data_MDiso

   subroutine Average_MDiso

   implicit none

      Sum_Temp       =  Sum_Temp       + Temp

      Sum_E_BondMol  =  Sum_E_BondMol  + E_BondMol
      Sum_E_AngleMol =  Sum_E_AngleMol + E_AngleMol
      Sum_E_UBMol    =  Sum_E_UBMol    + E_UBMol
      Sum_E_DihedMol =  Sum_E_DihedMol + E_DihedMol
      Sum_E_ImproMol =  Sum_E_ImproMol + E_ImproMol
      Sum_E_OptC     =  Sum_E_OptC     + E_OptC
      Sum_E_LJMol    =  Sum_E_LJMol    + E_LJMol
      Sum_E_ErMol    =  Sum_E_ErMol    + E_ErMol

      Sum_E_Int_Mol  =  Sum_E_Int_Mol  + E_Int_Mol
      Sum_E_Ext_Mol  =  Sum_E_Ext_Mol  + E_Ext_Mol
      Sum_Potential  =  Sum_Potential  + Potential
      Sum_Kinetic    =  Sum_Kinetic    + Kinetic
      Sum_EneSystem  =  Sum_EneSystem  + EneSystem
      Sum_Hamiltonian=  Sum_Hamiltonian+ Hamiltonian

      if(Ensemble=='NT') then
        Sum_ThermoBath =  Sum_ThermoBath + ThermoBath
      end if

   end subroutine Average_MDiso

end subroutine MD_iso


!######################################################################
!######################################################################


subroutine HMC_NPT

use Forplot

implicit none

integer :: i,j,k
integer :: Ns, TotalStep
character(len=80) :: String
#ifdef BMONI
integer :: eofile
#endif

   open(51,file='E_data/Ene_Intra.dat',status='unknown',form='formatted')
   open(52,file='E_data/Ene_Inter.dat',status='unknown',form='formatted')
   open(53,file='E_data/Ene_System.dat',status='unknown',form='formatted')
   open(54,file='E_data/Cell.dat',status='unknown',form='formatted')
   open(55,file='E_data/CellABC.dat',status='unknown',form='formatted')
   open(56,file='E_data/System.dat',status='unknown',form='formatted')
   open(57,file='E_data/Pressure.dat',status='unknown',form='formatted')
   if(Ensemble=='NtT') then
   open(58,file='E_data/Strain.dat',status='unknown',form='formatted')
   end if

   write(51,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(52,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ         Elec(real)   Elec(recip)  Elec(self)   Electrstatic '
   write(53,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E'
   write(54,'(a/a)') '# Cell matrix h in Angstrom ', &
   & '# Time [ps] h_xx    h_xy    h_xz    h_yx    h_yy    h_yz    h_zx    h_zy    h_zz  '
   write(55,'(a/a)') '# Cell size La Lb Lc angles', &
   & '# Time [ps] La      Lb      Lc      alpha   beta    gamma   '
   write(56,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  accept_ratio H[kcal/mol]  '
   write(57,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'
   if(Ensemble=='NtT') then
   write(58,'(a/a)') '# Strain energy [kcal/mol]', '# Time [ps]  Strain E '
   end if

   if(QAccAv) then

   open(61,file='E_data/Av_Ene_Intra.dat',status='unknown',form='formatted')
   open(62,file='E_data/Av_Ene_Inter.dat',status='unknown',form='formatted')
   open(63,file='E_data/Av_Ene_System.dat',status='unknown',form='formatted')
   open(64,file='E_data/Av_Cell.dat',status='unknown',form='formatted')
   open(65,file='E_data/Av_CellABC.dat',status='unknown',form='formatted')
   open(66,file='E_data/Av_System.dat',status='unknown',form='formatted')
   open(67,file='E_data/Av_Pressure.dat',status='unknown',form='formatted')
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   open(68,file='E_data/Av_Strain.dat',status='unknown',form='formatted')
   end if

   write(61,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(62,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ         Elec(real)   Elec(recip)  Elec(self)   Electrstatic '
   if((Ensemble=='NPT').or.(Ensemble=='NtT')) then
   write(63,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E    Bath Potentia'
   else if((Ensemble=='NPH').or.(Ensemble=='NtH')) then
   write(63,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E'
   end if
   write(64,'(a/a)') '# Cell matrix h in Angstrom ', &
   & '# Time [ps] h_xx    h_xy    h_xz    h_yx    h_yy    h_yz    h_zx    h_zy    h_zz  '
   write(65,'(a/a)') '# Cell size La Lb Lc angles', &
   & '# Time [ps] La      Lb      Lc      alpha   beta    gamma   '
   write(66,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  accept_ratio H[kcal/mol]  '
   write(67,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   write(68,'(a/a)') '# Strain energy [kcal/mol]', '# Time [ps]  Strain E '
   end if

   Sum_Temp       =  0.d0
   Sum_E_BondMol  =  0.d0
   Sum_E_AngleMol =  0.d0
   Sum_E_UBMol    =  0.d0
   Sum_E_DihedMol =  0.d0
   Sum_E_ImproMol =  0.d0
   Sum_E_OptC     =  0.d0
   Sum_E_LJMol    =  0.d0
   Sum_E_ErMol    =  0.d0
   Sum_E_EkMol    =  0.d0
   Sum_E_EsMol    =  0.d0
   Sum_E_ElecMol  =  0.d0
   Sum_E_Int_Mol  =  0.d0
   Sum_E_Ext_Mol  =  0.d0
   Sum_Potential  =  0.d0
   Sum_Kinetic    =  0.d0
   Sum_EneSystem  =  0.d0
   Sum_Pot_p      =  0.d0
   Sum_Hamiltonian=  0.d0
   Sum_Area       =  0.d0
   Sum_Volume     =  0.d0
   Sum_density    =  0.d0
   Sum_LenA       =  0.d0
   Sum_LenB       =  0.d0
   Sum_LenC       =  0.d0
   Sum_AngBC      =  0.d0
   Sum_AngCA      =  0.d0
   Sum_AngAB      =  0.d0
   Sum_H          =  0.d0
   Sum_Pressure   =  0.d0

   Sum_Accept_Ratio = 0.d0

   end if

   TotalStep = 0

   do i = 1 , Nfile

#ifdef BMONI
     open(11,file=trim(Energy_File(i)),status='old',form='unformatted')
#else
     open(11,file=trim(Energy_File(i)),status='old',form='formatted')
#endif
     write(*,*) 'File ',trim(Energy_File(i)),' was opened.'

#ifdef BMONI

     Ns = 0

     do

       read(11,iostat=eofile) Timeps, Temp, Accept_Ratio
       if(eofile /= 0) exit

       Ns = Ns + 1
       TotalStep = TotalStep + 1

       read(11)   E_BondMol,  E_AngleMol, E_UBMol,         &
       &          E_DihedMol, E_ImproMol, E_OptC,          &
       &          E_LJMol, E_ErMol,                        &
       &          E_EkMol, E_EsMol,   E_ElecMol,           &
       &          E_Int_Mol, E_Ext_Mol, Potential,         &
       &          Kinetic, EneSystem
       if(Ensemble=='NtT') then
       read(11) Pot_p
       end if
       read(11)  Hamiltonian, Area, Volume, density,   &
       &         LenA, LenB, LenC, AngBC, AngCA, AngAB
       read(11)  H, Pressure

#else

     do

       read(11,'(a80)') String

       if(String(3:58)=='#############< Start Hybrid MC iteration >##############') then

         read(11,*)

         Ns = 0

         do

           read(11,'(a80)') String

           if(String(1:31)=='-------------------------------') exit

           Ns = Ns + 1

           TotalStep = TotalStep + 1

           read(String,'(f12.4,f10.2,f10.4)') &
           &          Timeps, Temp, Accept_Ratio

           read(11,'(6d13.5/5d13.5/5d13.5)')                  &
           &          E_BondMol,  E_AngleMol, E_UBMol,         &
           &          E_DihedMol, E_ImproMol, E_OptC,          &
           &          E_LJMol, E_ErMol,                        &
           &          E_EkMol, E_EsMol,   E_ElecMol,           &
           &          E_Int_Mol, E_Ext_Mol, Potential,         &
           &          Kinetic, EneSystem
           if(Ensemble=='NtT') then
           read(11,'(d13.5)') Pot_p
           end if
           read(11,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')   &
           &         Hamiltonian, Area, Volume, density,   &
           &         LenA, LenB, LenC, AngBC, AngCA, AngAB
           read(11,'(3(3f8.3,5x,3f12.4/))')                  &
           &         (H(1,j),j=1,3) , (Pressure(1,j),j=1,3), &
           &         (H(2,j),j=1,3) , (Pressure(2,j),j=1,3), &
           &         (H(3,j),j=1,3) , (Pressure(3,j),j=1,3)

#endif

           if(QAccAv.and.TotalStep>Nskip) call Average_HMCNPT

           if(mod(Ns,Nst)==0) call Out_data_HMCNPT

#ifndef BMONI

         end do

       end if

       if(String(1:32)==' HMC simulation ended normally! ') exit

#endif

     end do

     close(11)

   end do

Contains

   subroutine Out_data_HMCNPT

   implicit none

   real(8) :: xinv

     write(51,'(f11.3,6d13.5)')                  &
     &  Timeps, E_BondMol,  E_AngleMol, E_UBMol, &
     &          E_DihedMol, E_ImproMol, E_OptC

     write(52,'(f11.3,5d13.5)')                  &
     &  Timeps, E_LJMol, E_ErMol,                &
     &          E_EkMol, E_EsMol,   E_ElecMol

     write(53,'(f11.3,5d13.5)')                    &
     &  Timeps, E_Int_Mol, E_Ext_Mol, Potential,   &
     &          Kinetic, EneSystem

        write(56,'(f11.3,f15.5,3d15.8,f10.6,f10.4,d16.8)')  &
     &   Timeps, Temp, Area, Volume/Area, Volume, density, Accept_Ratio, Hamiltonian

        write(54,'(f11.3,9f9.2)') &
     &   Timeps, ((H(j,k),k=1,3),j=1,3)

     write(55,'(f11.3,6f9.2)') &
     &   Timeps, LenA, LenB, LenC, AngBC, AngCA, AngAB

     write(57,'(f11.3,9d15.7)') &
     &   Timeps, ((Pressure(j,k),k=1,3),j=1,3)

     if(Ensemble=='NtT') then
     write(58,'(f11.3,d13.5)') Timeps, Pot_p
     end if

     if(QAccAv.and.TotalStep>Nskip) then

       xinv = 1.d0 / dble(TotalStep-Nskip)

       Av_Temp       =  Sum_Temp       * xinv
       Av_E_BondMol  =  Sum_E_BondMol  * xinv
       Av_E_AngleMol =  Sum_E_AngleMol * xinv
       Av_E_UBMol    =  Sum_E_UBMol    * xinv
       Av_E_DihedMol =  Sum_E_DihedMol * xinv
       Av_E_ImproMol =  Sum_E_ImproMol * xinv
       Av_E_OptC     =  Sum_E_OptC     * xinv
       Av_E_LJMol    =  Sum_E_LJMol    * xinv
       Av_E_ErMol    =  Sum_E_ErMol    * xinv
       Av_E_EkMol    =  Sum_E_EkMol    * xinv
       Av_E_EsMol    =  Sum_E_EsMol    * xinv
       Av_E_ElecMol  =  Sum_E_ElecMol  * xinv
       Av_E_Int_Mol  =  Sum_E_Int_Mol  * xinv
       Av_E_Ext_Mol  =  Sum_E_Ext_Mol  * xinv
       Av_Potential  =  Sum_Potential  * xinv
       Av_Kinetic    =  Sum_Kinetic    * xinv
       Av_EneSystem  =  Sum_EneSystem  * xinv
       Av_Pot_p      =  Sum_Pot_p      * xinv
       Av_Hamiltonian=  Sum_Hamiltonian* xinv
       Av_Area       =  Sum_Area       * xinv
       Av_Volume     =  Sum_Volume     * xinv
       Av_density    =  Sum_density    * xinv
       Av_LenA       =  Sum_LenA       * xinv
       Av_LenB       =  Sum_LenB       * xinv
       Av_LenC       =  Sum_LenC       * xinv
       Av_AngBC      =  Sum_AngBC      * xinv
       Av_AngCA      =  Sum_AngCA      * xinv
       Av_AngAB      =  Sum_AngAB      * xinv
       Av_H          =  Sum_H          * xinv
       Av_Pressure   =  Sum_Pressure   * xinv

       Av_Accept_Ratio = Sum_Accept_Ratio * xinv

       write(61,'(f11.3,6d13.5)')                  &
       &  Timeps, Av_E_BondMol,  Av_E_AngleMol, Av_E_UBMol, &
       &          Av_E_DihedMol, Av_E_ImproMol, Av_E_OptC

       write(62,'(f11.3,5d13.5)')                  &
       &  Timeps, Av_E_LJMol, Av_E_ErMol,                &
       &          Av_E_EkMol, Av_E_EsMol,   Av_E_ElecMol

       write(63,'(f11.3,5d13.5)')                    &
       &  Timeps, Av_E_Int_Mol, Av_E_Ext_Mol, Av_Potential,   &
       &          Av_Kinetic, Av_EneSystem

       write(66,'(f11.3,f15.5,3d15.8,f10.6,f10.4,d16.8)')  &
       &   Timeps, Av_Temp, Av_Area, Av_Volume/Av_Area, Av_Volume, Av_density, &
       &   Av_Accept_Ratio, Av_Hamiltonian

       write(64,'(f11.3,9f9.2)') &
       &   Timeps, ((Av_H(j,k),k=1,3),j=1,3)

       write(65,'(f11.3,6f9.2)') &
       &   Timeps, Av_LenA, Av_LenB, Av_LenC, Av_AngBC, Av_AngCA, Av_AngAB

       write(67,'(f11.3,9d15.7)') &
       &   Timeps, ((Av_Pressure(j,k),k=1,3),j=1,3)

       if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
       write(68,'(f11.3,d13.5)') Timeps, Av_Pot_p
       end if

     end if

   end subroutine Out_data_HMCNPT

   subroutine Average_HMCNPT

   implicit none

      Sum_Temp       =  Sum_Temp       + Temp

      Sum_E_BondMol  =  Sum_E_BondMol  + E_BondMol
      Sum_E_AngleMol =  Sum_E_AngleMol + E_AngleMol
      Sum_E_UBMol    =  Sum_E_UBMol    + E_UBMol
      Sum_E_DihedMol =  Sum_E_DihedMol + E_DihedMol
      Sum_E_ImproMol =  Sum_E_ImproMol + E_ImproMol
      Sum_E_OptC     =  Sum_E_OptC     + E_OptC
      Sum_E_LJMol    =  Sum_E_LJMol    + E_LJMol
      Sum_E_ErMol    =  Sum_E_ErMol    + E_ErMol
      Sum_E_EkMol    =  Sum_E_EkMol    + E_EkMol
      Sum_E_EsMol    =  Sum_E_EsMol    + E_EsMol
      Sum_E_ElecMol  =  Sum_E_ElecMol  + E_ElecMol

      Sum_E_Int_Mol  =  Sum_E_Int_Mol  + E_Int_Mol
      Sum_E_Ext_Mol  =  Sum_E_Ext_Mol  + E_Ext_Mol
      Sum_Potential  =  Sum_Potential  + Potential
      Sum_Kinetic    =  Sum_Kinetic    + Kinetic
      Sum_EneSystem  =  Sum_EneSystem  + EneSystem
      Sum_Hamiltonian=  Sum_Hamiltonian+ Hamiltonian
      Sum_Area       =  Sum_Area       + Area
      Sum_Volume     =  Sum_Volume     + Volume
      Sum_density    =  Sum_density    + density
      Sum_LenA       =  Sum_LenA       + LenA
      Sum_LenB       =  Sum_LenB       + LenB
      Sum_LenC       =  Sum_LenC       + LenC
      Sum_AngBC      =  Sum_AngBC      + AngBC
      Sum_AngCA      =  Sum_AngCA      + AngCA
      Sum_AngAB      =  Sum_AngAB      + AngAB
      Sum_H          =  Sum_H          + H
      Sum_Pressure   =  Sum_Pressure   + Pressure

      if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
        Sum_Pot_p      =  Sum_Pot_p      + Pot_p
      end if

   end subroutine Average_HMCNPT

end subroutine HMC_NPT


!######################################################################
!######################################################################


subroutine HMC_NVT
use Forplot

implicit none

integer :: i,j,k
integer :: Ns, TotalStep
character(len=80) :: String
#ifdef BMONI
integer :: eofile
#endif

   open(51,file='E_data/Ene_Intra.dat',status='unknown',form='formatted')
   open(52,file='E_data/Ene_Inter.dat',status='unknown',form='formatted')
   open(53,file='E_data/Ene_System.dat',status='unknown',form='formatted')
   open(56,file='E_data/System.dat',status='unknown',form='formatted')
   open(57,file='E_data/Pressure.dat',status='unknown',form='formatted')

   write(51,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(52,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ         Elec(real)   Elec(recip)  Elec(self)   Electrstatic '
   write(53,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E '
    write(56,'(a/a)') '# Temperature & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]  accept_ratio[%]  H[kcal/mol]  '
   write(57,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'

   if(QAccAv) then

   open(61,file='E_data/Av_Ene_Intra.dat',status='unknown',form='formatted')
   open(62,file='E_data/Av_Ene_Inter.dat',status='unknown',form='formatted')
   open(63,file='E_data/Av_Ene_System.dat',status='unknown',form='formatted')
   open(66,file='E_data/Av_System.dat',status='unknown',form='formatted')
   open(67,file='E_data/Av_Pressure.dat',status='unknown',form='formatted')

   write(61,'(a/a)') '# Intramolecular Energy [kcal/mol]', &
   &  '# Time [ps]  Bond       Angle        UB           Dihedral     Improper     Optional'
   write(62,'(a/a)') '# Nonbond Energy [kcal/mol]', &
   & '# Time [ps]  LJ         Elec(real)   Elec(recip)  Elec(self)   Electrstatic '
   write(63,'(a/a)') '# Total Energy [kcal/mol]', &
  & '# Time [ps]  Intramol   Intermol     Potential    Kinetic E    Internal E '
   write(66,'(a/a)') '# Temperature & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]  accept_ratio[%]  H[kcal/mol]  '
   write(67,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'

   Sum_Temp       =  0.d0
   Sum_E_BondMol  =  0.d0
   Sum_E_AngleMol =  0.d0
   Sum_E_UBMol    =  0.d0
   Sum_E_DihedMol =  0.d0
   Sum_E_ImproMol =  0.d0
   Sum_E_OptC     =  0.d0
   Sum_E_LJMol    =  0.d0
   Sum_E_ErMol    =  0.d0
   Sum_E_EkMol    =  0.d0
   Sum_E_EsMol    =  0.d0
   Sum_E_ElecMol  =  0.d0
   Sum_E_Int_Mol  =  0.d0
   Sum_E_Ext_Mol  =  0.d0
   Sum_Potential  =  0.d0
   Sum_Kinetic    =  0.d0
   Sum_EneSystem  =  0.d0
   Sum_Pot_p      =  0.d0
   Sum_Hamiltonian=  0.d0
   Sum_Area       =  0.d0
   Sum_Volume     =  0.d0
   Sum_density    =  0.d0
   Sum_LenA       =  0.d0
   Sum_LenB       =  0.d0
   Sum_LenC       =  0.d0
   Sum_AngBC      =  0.d0
   Sum_AngCA      =  0.d0
   Sum_AngAB      =  0.d0
   Sum_H          =  0.d0
   Sum_Pressure   =  0.d0

   Sum_Accept_Ratio = 0.d0

   end if

   TotalStep = 0

   do i = 1 , Nfile

#ifdef BMONI
     open(11,file=trim(Energy_File(i)),status='old',form='unformatted')
#else
     open(11,file=trim(Energy_File(i)),status='old',form='formatted')
#endif
     write(*,*) 'File ',trim(Energy_File(i)),' was opened.'

#ifdef BMONI

     Ns = 0

     do

       read(11,iostat=eofile) Timeps, Temp, Accept_Ratio
       if(eofile /= 0) exit

       Ns = Ns + 1
       TotalStep = TotalStep + 1

       read(11)  E_BondMol, E_AngleMol, E_UBMol,     &
       &         E_DihedMol,E_ImproMol, E_OptC,      &
       &         E_LJMol, E_ErMol,                   &
       &         E_EkMol, E_EsMol,   E_ElecMol,      &
       &         E_Int_Mol, E_Ext_Mol, Potential,    &
       &         Kinetic, EneSystem, Pressure

#else

     do

       read(11,'(a80)') String

       if(String(3:58)=='#############< Start Hybrid MC iteration >##############') then

         read(11,*)

         Ns = 0

         do

           read(11,'(a80)') String

           if(String(1:31)=='-------------------------------') exit

           Ns = Ns + 1

           TotalStep = TotalStep + 1

           read(String,'(f12.4,f10.2,f10.4)') &
           &          Timeps, Temp, Accept_Ratio

           read(11,'(6d13.5/5d13.5/5d13.5/3(3f12.4/))') &
           &          E_BondMol, E_AngleMol, E_UBMol,     &
           &          E_DihedMol,E_ImproMol, E_OptC,      &
           &          E_LJMol, E_ErMol,                   &
           &          E_EkMol, E_EsMol,   E_ElecMol,      &
           &          E_Int_Mol, E_Ext_Mol, Potential,    &
           &          Kinetic, EneSystem,                 &
           &          ( Pressure(1,j) , j = 1 , 3 ),      &
           &          ( Pressure(2,j) , j = 1 , 3 ),      &
           &          ( Pressure(3,j) , j = 1 , 3 )
#endif
           if(QAccAv.and.TotalStep>Nskip) call Average_HMCNVT

           if(mod(Ns,Nst)==0) call Out_data_HMCNVT

#ifndef BMONI

         end do

       end if

       if(String(1:32)==' HMC simulation ended normally! ') exit

#endif

     end do

     close(11)

   end do

Contains

   subroutine Out_data_HMCNVT

   implicit none

   real(8) :: xinv

     write(51,'(f11.3,6d13.5)')                  &
     &  Timeps, E_BondMol,  E_AngleMol, E_UBMol, &
     &          E_DihedMol, E_ImproMol, E_OptC

     write(52,'(f11.3,5d13.5)')                  &
     &  Timeps, E_LJMol, E_ErMol,                &
     &          E_EkMol, E_EsMol,   E_ElecMol

     write(53,'(f11.3,5d13.5)')                    &
     &  Timeps, E_Int_Mol, E_Ext_Mol, Potential,   &
     &          Kinetic, EneSystem

     write(56,'(f11.3,f10.2,f10.4,d16.8)')  &
     &   Timeps, Accept_Ratio, Temp, Hamiltonian

     write(57,'(f11.3,9d15.7)') &
     &   Timeps, ((Pressure(j,k),k=1,3),j=1,3)

     if(QAccAv.and.TotalStep>Nskip) then

       xinv = 1.d0 / dble(TotalStep-Nskip)

       Av_Temp       =  Sum_Temp       * xinv
       Av_E_BondMol  =  Sum_E_BondMol  * xinv
       Av_E_AngleMol =  Sum_E_AngleMol * xinv
       Av_E_UBMol    =  Sum_E_UBMol    * xinv
       Av_E_DihedMol =  Sum_E_DihedMol * xinv
       Av_E_ImproMol =  Sum_E_ImproMol * xinv
       Av_E_OptC     =  Sum_E_OptC     * xinv
       Av_E_LJMol    =  Sum_E_LJMol    * xinv
       Av_E_ErMol    =  Sum_E_ErMol    * xinv
       Av_E_EkMol    =  Sum_E_EkMol    * xinv
       Av_E_EsMol    =  Sum_E_EsMol    * xinv
       Av_E_ElecMol  =  Sum_E_ElecMol  * xinv
       Av_E_Int_Mol  =  Sum_E_Int_Mol  * xinv
       Av_E_Ext_Mol  =  Sum_E_Ext_Mol  * xinv
       Av_Potential  =  Sum_Potential  * xinv
       Av_Kinetic    =  Sum_Kinetic    * xinv
       Av_EneSystem  =  Sum_EneSystem  * xinv
       Av_Pot_p      =  Sum_Pot_p      * xinv
       Av_Hamiltonian=  Sum_Hamiltonian* xinv
       Av_Pressure   =  Sum_Pressure   * xinv

       Av_Accept_Ratio = Sum_Accept_Ratio * xinv

       write(61,'(f11.3,6d13.5)')                  &
       &  Timeps, Av_E_BondMol,  Av_E_AngleMol, Av_E_UBMol, &
       &          Av_E_DihedMol, Av_E_ImproMol, Av_E_OptC

       write(62,'(f11.3,5d13.5)')                  &
       &  Timeps, Av_E_LJMol, Av_E_ErMol,                &
       &          Av_E_EkMol, Av_E_EsMol,   Av_E_ElecMol

       write(63,'(f11.3,5d13.5)')                    &
       &  Timeps, Av_E_Int_Mol, Av_E_Ext_Mol, Av_Potential,   &
       &          Av_Kinetic, Av_EneSystem

       write(66,'(f11.3,f10.2,f10.4,d16.8)')  &
       &   Timeps, Av_Accept_Ratio, Av_Temp, Av_Hamiltonian

       write(67,'(f11.3,9d15.7)') &
       &   Timeps, ((Av_Pressure(j,k),k=1,3),j=1,3)

     end if

   end subroutine Out_data_HMCNVT

   subroutine Average_HMCNVT

   implicit none

      Sum_Temp       =  Sum_Temp       + Temp

      Sum_E_BondMol  =  Sum_E_BondMol  + E_BondMol
      Sum_E_AngleMol =  Sum_E_AngleMol + E_AngleMol
      Sum_E_UBMol    =  Sum_E_UBMol    + E_UBMol
      Sum_E_DihedMol =  Sum_E_DihedMol + E_DihedMol
      Sum_E_ImproMol =  Sum_E_ImproMol + E_ImproMol
      Sum_E_OptC     =  Sum_E_OptC     + E_OptC
      Sum_E_LJMol    =  Sum_E_LJMol    + E_LJMol
      Sum_E_ErMol    =  Sum_E_ErMol    + E_ErMol
      Sum_E_EkMol    =  Sum_E_EkMol    + E_EkMol
      Sum_E_EsMol    =  Sum_E_EsMol    + E_EsMol
      Sum_E_ElecMol  =  Sum_E_ElecMol  + E_ElecMol

      Sum_E_Int_Mol  =  Sum_E_Int_Mol  + E_Int_Mol
      Sum_E_Ext_Mol  =  Sum_E_Ext_Mol  + E_Ext_Mol
      Sum_Potential  =  Sum_Potential  + Potential
      Sum_Kinetic    =  Sum_Kinetic    + Kinetic
      Sum_EneSystem  =  Sum_EneSystem  + EneSystem
      Sum_Hamiltonian=  Sum_Hamiltonian+ Hamiltonian

      Sum_Pressure   =  Sum_Pressure   + Pressure

      Sum_Accept_Ratio = Sum_Accept_Ratio + Accept_Ratio


   end subroutine Average_HMCNVT


end subroutine HMC_NVT


!######################################################################
!######################################################################


subroutine PIMD_NPT

use Forplot

implicit none

integer :: i,j,k
integer :: Ns, TotalStep
character(len=80) :: String
#ifdef BMONI
integer :: eofile
#endif

   open(53,file='E_data/Ene_System.dat',status='unknown',form='formatted')
   open(54,file='E_data/Cell.dat',status='unknown',form='formatted')
   open(55,file='E_data/CellABC.dat',status='unknown',form='formatted')
   open(56,file='E_data/System.dat',status='unknown',form='formatted')
   open(57,file='E_data/Pressure.dat',status='unknown',form='formatted')
   if(Ensemble=='NtT') then
   open(58,file='E_data/Strain.dat',status='unknown',form='formatted')
   end if

   write(53,'(a/a)') '# Energy [kcal/mol]', &
  & '# Time [ps]  Potential  Kinetic      Qkinetic     Internal E    Bath Potentia'
   write(54,'(a/a)') '# Cell matrix h in Angstrom ', &
   & '# Time [ps] h_xx    h_xy    h_xz    h_yx    h_yy    h_yz    h_zx    h_zy    h_zz  '
   write(55,'(a/a)') '# Cell size La Lb Lc angles', &
   & '# Time [ps] La      Lb      Lc      alpha   beta    gamma   '
   write(56,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  H[kcal/mol]  '
   write(57,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   write(58,'(a/a)') '# Strain energy [kcal/mol]', '# Time [ps]  Strain E '
   end if

   if(QAccAv) then

   open(63,file='E_data/Av_Ene_System.dat',status='unknown',form='formatted')
   open(64,file='E_data/Av_Cell.dat',status='unknown',form='formatted')
   open(65,file='E_data/Av_CellABC.dat',status='unknown',form='formatted')
   open(66,file='E_data/Av_System.dat',status='unknown',form='formatted')
   open(67,file='E_data/Av_Pressure.dat',status='unknown',form='formatted')
   if(Ensemble=='NtT') then
   open(68,file='E_data/Av_Strain.dat',status='unknown',form='formatted')
   end if

   write(63,'(a/a)') '# Energy [kcal/mol]', &
  & '# Time [ps]  Potential  Kinetic      Qkinetic     Internal E    Bath Potentia'
   write(64,'(a/a)') '# Cell matrix h in Angstrom ', &
   & '# Time [ps] h_xx    h_xy    h_xz    h_yx    h_yy    h_yz    h_zx    h_zy    h_zz  '
   write(65,'(a/a)') '# Cell size La Lb Lc angles', &
   & '# Time [ps] La      Lb      Lc      alpha   beta    gamma   '
   write(66,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  H[kcal/mol]  '
   write(67,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'
   if((Ensemble=='NtT').or.(Ensemble=='NtH')) then
   write(68,'(a/a)') '# Strain energy [kcal/mol]', '# Time [ps]  Strain E '
   end if

   Sum_Temp      =  0.d0
   Sum_Qkinetic  =  0.d0
   Sum_Potential =  0.d0
   Sum_Kinetic   =  0.d0
   Sum_EneSystem =  0.d0
   Sum_ThermoBath=  0.d0
   Sum_Pot_p     =  0.d0
   Sum_Hamiltonian=  0.d0
   Sum_Area      =  0.d0
   Sum_Volume    =  0.d0
   Sum_density   =  0.d0
   Sum_LenA      =  0.d0
   Sum_LenB      =  0.d0
   Sum_LenC      =  0.d0
   Sum_AngBC     =  0.d0
   Sum_AngCA     =  0.d0
   Sum_AngAB     =  0.d0
   Sum_H         =  0.d0
   Sum_Pressure  =  0.d0

   end if

   TotalStep = 0

   do i = 1 , Nfile

#ifdef BMONI
     open(11,file=trim(Energy_File(i)),status='old',form='unformatted')
#else
     open(11,file=trim(Energy_File(i)),status='old',form='formatted')
#endif
     write(*,*) 'File ',trim(Energy_File(i)),' was opened.'

#ifdef BMONI

     Ns = 0

     do

       read(11,iostat=eofile) Timeps, Temp
       if(eofile /= 0) exit

       Ns = Ns + 1
       TotalStep = TotalStep + 1

       read(11) Potential, Kinetic, Qkinetic, EneSystem, ThermoBath

       if(Ensemble=='NtT') then
       read(11) Pot_p
       end if
       read(11) Hamiltonian, Area, Volume, density,      &
       &        LenA, LenB, LenC, AngBC, AngCA, AngAB
       read(11) H, Pressure

#else

     do

       read(11,'(a80)') String

       if((String(3:55)=='#############< Start Path Integral MD >##############').or. &
       &  (String(3:50)=='#############< Start Centroid MD >##############')) then

         read(11,*)

         Ns = 0

         do

           read(11,'(a80)') String

           if(String(1:31)=='-------------------------------') exit

           Ns = Ns + 1
           TotalStep = TotalStep + 1

           read(String,'(f12.4,f10.2)') &
           &          Timeps, Temp

           read(11,'(5d13.5)') Potential, Kinetic, Qkinetic, EneSystem, ThermoBath

           if(Ensemble=='NtT') then
           read(11,'(d13.5)') Pot_p
           end if
           read(11,'(d16.8,2x,f10.3,f10.1,f10.6/6f8.2)')       &
           &          Hamiltonian, Area, Volume, density,      &
           &          LenA, LenB, LenC, AngBC, AngCA, AngAB
           read(11,'(3(3f8.3,5x,3f12.4/))')                    &
           &          (H(1,j),j=1,3) , (Pressure(1,j),j=1,3),  &
           &          (H(2,j),j=1,3) , (Pressure(2,j),j=1,3),  &
           &          (H(3,j),j=1,3) , (Pressure(3,j),j=1,3)
#endif

           if(QAccAv.and.TotalStep>Nskip) call Average_PIMDNPT
           if(mod(Ns,Nst)==0) call Out_data_PIMDNPT

#ifndef BMONI

         end do

       end if

       if((String(1:33)==' PIMD simulation ended normally! ').or. &
       &  (String(1:32)==' CMD simulation ended normally! ')) exit

#endif

     end do

     close(11)

   end do

Contains

   subroutine Out_data_PIMDNPT

   implicit none

   real(8) :: xinv

     write(53,'(f11.3,5d13.5)')                    &
     &  Timeps,  Potential, Kinetic, Qkinetic, EneSystem, ThermoBath

     write(56,'(f11.3,f15.5,3d15.8,f10.6,d16.8)')  &
     &   Timeps, Temp, Area, Volume/Area, Volume, density, Hamiltonian

     write(54,'(f11.3,9f9.2)') &
     &   Timeps, ((H(j,k),k=1,3),j=1,3)

     write(55,'(f11.3,6f9.2)') &
     &   Timeps, LenA, LenB, LenC, AngBC, AngCA, AngAB

     write(57,'(f11.3,9d15.7)') &
     &   Timeps, ((Pressure(j,k),k=1,3),j=1,3)

     if(Ensemble=='NtT') then
     write(58,'(f11.3,d13.5)') Timeps, Pot_p
     end if

     if(QAccAv.and.TotalStep>Nskip) then

       xinv = 1.d0 / dble(TotalStep-Nskip)

       Av_Temp       =  Sum_Temp       * xinv
       Av_Qkinetic   =  Sum_Qkinetic   * xinv
       Av_Potential  =  Sum_Potential  * xinv
       Av_Kinetic    =  Sum_Kinetic    * xinv
       Av_EneSystem  =  Sum_EneSystem  * xinv
       Av_ThermoBath =  Sum_ThermoBath * xinv
       Av_Pot_p      =  Sum_Pot_p      * xinv
       Av_Hamiltonian=  Sum_Hamiltonian* xinv
       Av_Area       =  Sum_Area       * xinv
       Av_Volume     =  Sum_Volume     * xinv
       Av_density    =  Sum_density    * xinv
       Av_LenA       =  Sum_LenA       * xinv
       Av_LenB       =  Sum_LenB       * xinv
       Av_LenC       =  Sum_LenC       * xinv
       Av_AngBC      =  Sum_AngBC      * xinv
       Av_AngCA      =  Sum_AngCA      * xinv
       Av_AngAB      =  Sum_AngAB      * xinv
       Av_H          =  Sum_H          * xinv
       Av_Pressure   =  Sum_Pressure   * xinv

       write(63,'(f11.3,5d13.5)')                    &
       &  Timeps,  Av_Potential, Av_Kinetic, Av_Qkinetic, Av_EneSystem, Av_ThermoBath

       write(66,'(f11.3,f15.5,3d15.8,f10.6,d16.8)')  &
       &   Timeps, Av_Temp, Av_Area, Av_Volume/Av_Area, Av_Volume, Av_density, Av_Hamiltonian

       write(64,'(f11.3,9f9.2)') &
       &   Timeps, ((Av_H(j,k),k=1,3),j=1,3)

       write(65,'(f11.3,6f9.2)') &
       &   Timeps, Av_LenA, Av_LenB, Av_LenC, Av_AngBC, Av_AngCA, Av_AngAB

       write(67,'(f11.3,9d15.7)') &
       &   Timeps, ((Av_Pressure(j,k),k=1,3),j=1,3)

       if(Ensemble=='NtT') then
       write(68,'(f11.3,d13.5)') Timeps, Av_Pot_p
       end if

     end if

   end subroutine Out_data_PIMDNPT

   subroutine Average_PIMDNPT

   implicit none

       Sum_Temp       =  Sum_Temp        + Temp       
       Sum_Qkinetic   =  Sum_Qkinetic    + Qkinetic   
       Sum_Potential  =  Sum_Potential   + Potential  
       Sum_Kinetic    =  Sum_Kinetic     + Kinetic    
       Sum_EneSystem  =  Sum_EneSystem   + EneSystem  
       Sum_ThermoBath =  Sum_ThermoBath  + ThermoBath 
       Sum_Pot_p      =  Sum_Pot_p       + Pot_p      
       Sum_Hamiltonian=  Sum_Hamiltonian + Hamiltonian
       Sum_Area       =  Sum_Area        + Area       
       Sum_Volume     =  Sum_Volume      + Volume     
       Sum_density    =  Sum_density     + density    
       Sum_LenA       =  Sum_LenA        + LenA       
       Sum_LenB       =  Sum_LenB        + LenB       
       Sum_LenC       =  Sum_LenC        + LenC       
       Sum_AngBC      =  Sum_AngBC       + AngBC      
       Sum_AngCA      =  Sum_AngCA       + AngCA      
       Sum_AngAB      =  Sum_AngAB       + AngAB      
       Sum_H          =  Sum_H           + H          
       Sum_Pressure   =  Sum_Pressure    + Pressure   

   end subroutine Average_PIMDNPT

end subroutine PIMD_NPT


!######################################################################
!######################################################################


subroutine PIMD_NVT

use Forplot

implicit none

integer :: i,j,k
integer :: Ns, TotalStep
character(len=80) :: String
#ifdef BMONI
integer :: eofile
#endif

   open(53,file='E_data/Ene_System.dat',status='unknown',form='formatted')
   open(56,file='E_data/System.dat',status='unknown',form='formatted')
   open(57,file='E_data/Pressure.dat',status='unknown',form='formatted')

   write(53,'(a/a)') '# Energy [kcal/mol]', &
  & '# Time [ps]  Potential  Kinetic      Qkinetic     Internal E    Bath Potentia'
   write(56,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  H[kcal/mol]  '
   write(57,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'

   if(QAccAv) then

   open(63,file='E_data/Av_Ene_System.dat',status='unknown',form='formatted')
   open(66,file='E_data/Av_System.dat',status='unknown',form='formatted')
   open(67,file='E_data/Av_Pressure.dat',status='unknown',form='formatted')

   write(63,'(a/a)') '# Energy [kcal/mol]', &
  & '# Time [ps]  Potential  Kinetic      Qkinetic     Internal E    Bath Potentia'
   write(66,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  H[kcal/mol]  '
   write(67,'(a/a)') '# Pressure tensor in MPa ', &
   & '# Time [ps] P_xx      P_xy      P_xz      P_yx      P_yy      P_yz      P_zx      P_zy      P_zz'

   Sum_Temp      =  0.d0
   Sum_Qkinetic  =  0.d0
   Sum_Potential =  0.d0
   Sum_Kinetic   =  0.d0
   Sum_EneSystem =  0.d0
   Sum_ThermoBath=  0.d0
   Sum_Hamiltonian=  0.d0
   Sum_Pressure  =  0.d0

   end if

   TotalStep = 0

   do i = 1 , Nfile

#ifdef BMONI
     open(11,file=trim(Energy_File(i)),status='old',form='unformatted')
#else
     open(11,file=trim(Energy_File(i)),status='old',form='formatted')
#endif
     write(*,*) 'File ',trim(Energy_File(i)),' was opened.'

#ifdef BMONI

     Ns = 0

     do

       read(11,iostat=eofile) Timeps, Temp
       if(eofile /= 0) exit

       Ns = Ns + 1
       TotalStep = TotalStep + 1

       read(11) Potential, Kinetic, Qkinetic, EneSystem, &
       &        ThermoBath, Hamiltonian, Pressure

#else

     do

       read(11,'(a80)') String

       if((String(3:55)=='#############< Start Path Integral MD >##############').or. &
       &  (String(3:50)=='#############< Start Centroid MD >##############')) then

         read(11,*)

         Ns = 0

         do

           read(11,'(a80)') String

           if(String(1:31)=='-------------------------------') exit

           Ns = Ns + 1
           TotalStep = TotalStep + 1

           read(String,'(f12.4,f10.2)') &
           &          Timeps, Temp

           read(11,'(5d13.5/d16.8/3(3f12.4/))') &
           &          Potential, Kinetic, Qkinetic, EneSystem, &
           &          ThermoBath, Hamiltonian,       &
           &          ( Pressure(1,j) , j = 1 , 3 ), &
           &          ( Pressure(2,j) , j = 1 , 3 ), &
           &          ( Pressure(3,j) , j = 1 , 3 )

#endif

           if(QAccAv.and.TotalStep>Nskip) call Average_PIMDNVT
           if(mod(Ns,Nst)==0) call Out_data_PIMDNVT

#ifndef BMONI

         end do

       end if

       if((String(1:33)==' PIMD simulation ended normally! ').or. &
       &  (String(1:32)==' CMD simulation ended normally! ')) exit

#endif

     end do

     close(11)

   end do

Contains

   subroutine Out_data_PIMDNVT

   implicit none
   real(8) :: xinv

     write(53,'(f11.3,5d13.5)')                    &
     &  Timeps,  Potential, Kinetic, Qkinetic, EneSystem, ThermoBath

     write(56,'(f11.3,f10.2,d16.8)')  &
     &   Timeps, Temp, Hamiltonian

     write(57,'(f11.3,9d15.7)') &
     &   Timeps, ((Pressure(j,k),k=1,3),j=1,3)

     if(QAccAv.and.TotalStep>Nskip) then

       xinv = 1.d0 / dble(TotalStep-Nskip)

       Av_Temp       =  Sum_Temp       * xinv
       Av_Qkinetic   =  Sum_Qkinetic   * xinv
       Av_Potential  =  Sum_Potential  * xinv
       Av_Kinetic    =  Sum_Kinetic    * xinv
       Av_EneSystem  =  Sum_EneSystem  * xinv
       Av_ThermoBath =  Sum_ThermoBath * xinv
       Av_Hamiltonian=  Sum_Hamiltonian* xinv
       Av_Pressure   =  Sum_Pressure   * xinv

       write(63,'(f11.3,5d13.5)')                    &
       &  Timeps,  Av_Potential, Av_Kinetic, Av_Qkinetic, Av_EneSystem, Av_ThermoBath

       write(66,'(f11.3,f15.5,d16.8)')  &
       &   Timeps, Av_Temp, Av_Hamiltonian

       write(67,'(f11.3,9d15.7)') &
       &   Timeps, ((Av_Pressure(j,k),k=1,3),j=1,3)

     end if

   end subroutine Out_data_PIMDNVT

   subroutine Average_PIMDNVT

   implicit none

       Sum_Temp       =  Sum_Temp        + Temp       
       Sum_Qkinetic   =  Sum_Qkinetic    + Qkinetic   
       Sum_Potential  =  Sum_Potential   + Potential  
       Sum_Kinetic    =  Sum_Kinetic     + Kinetic    
       Sum_EneSystem  =  Sum_EneSystem   + EneSystem  
       Sum_ThermoBath =  Sum_ThermoBath  + ThermoBath 
       Sum_Hamiltonian=  Sum_Hamiltonian + Hamiltonian
       Sum_Pressure   =  Sum_Pressure    + Pressure   

   end subroutine Average_PIMDNVT


end subroutine PIMD_NVT


!######################################################################
!######################################################################


subroutine PIMD_iso

use Forplot

implicit none

integer :: i
integer :: Ns, TotalStep
character(len=80) :: String
#ifdef BMONI
integer :: eofile
#endif

   open(53,file='E_data/Ene_System.dat',status='unknown',form='formatted')
   open(56,file='E_data/System.dat',status='unknown',form='formatted')

   write(53,'(a/a)') '# Energy [kcal/mol]', &
  & '# Time [ps]  Potential  Kinetic      Qkinetic     Bath Potentia'
   write(56,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  H[kcal/mol]  '

   if(QAccAv) then

   open(63,file='E_data/Av_Ene_System.dat',status='unknown',form='formatted')
   open(66,file='E_data/Av_System.dat',status='unknown',form='formatted')

   write(63,'(a/a)') '# Energy [kcal/mol]', &
  & '# Time [ps]  Potential  Kinetic      Qkinetic     Bath Potentia'
   write(66,'(a/a)') '# Temperature & System dimension & Hamiltonian ', &
   & '# Time [ps]  Temp.[K]     A(xy)[A^2]   Lz[A]   Vol.[A^3]  dens.[g/cc]  H[kcal/mol]  '

   Sum_Temp       =  0.d0
   Sum_Qkinetic   =  0.d0
   Sum_Potential  =  0.d0
   Sum_Kinetic    =  0.d0
   Sum_EneSystem  =  0.d0
   Sum_ThermoBath =  0.d0
   Sum_Hamiltonian=  0.d0

   end if

   TotalStep = 0

   do i = 1 , Nfile

#ifdef BMONI
     open(11,file=trim(Energy_File(i)),status='old',form='unformatted')
#else
     open(11,file=trim(Energy_File(i)),status='old',form='formatted')
#endif
     write(*,*) 'File ',trim(Energy_File(i)),' was opened.'

#ifdef BMONI

     Ns = 0

     do

       read(11,iostat=eofile) Timeps, Temp
       if(eofile /= 0) exit

       Ns = Ns + 1
       TotalStep = TotalStep + 1

       read(11) Potential, Kinetic, Qkinetic,  &
       &        ThermoBath, Hamiltonian

#else

     do

       read(11,'(a80)') String

       if((String(3:55)=='#############< Start Path Integral MD >##############').or. &
       &  (String(3:50)=='#############< Start Centroid MD >##############')) then

         read(11,*)

         Ns = 0

         do

           read(11,'(a80)') String

           if(String(1:31)=='-------------------------------') exit

           Ns = Ns + 1
           TotalStep = TotalStep + 1

           read(String,'(f12.4,f10.2)') &
           &          Timeps, Temp

           read(11,'(4d13.5/d16.8/)') &
           &          Potential, Kinetic, Qkinetic,  &
           &          ThermoBath, Hamiltonian

#endif

           if(QAccAv.and.TotalStep>Nskip) call Average_PIMDiso
           if(mod(Ns,Nst)==0) call Out_data_PIMDiso

#ifndef BMONI

         end do

       end if

       if((String(1:33)==' PIMD simulation ended normally! ').or. &
       &  (String(1:32)==' CMD simulation ended normally! ')) exit
#endif

     end do

     close(11)

   end do

Contains

   subroutine Out_data_PIMDiso

   implicit none
   real(8) :: xinv

     write(53,'(f11.3,5d13.5)')                    &
     &  Timeps,  Potential, Kinetic, Qkinetic, ThermoBath

     write(56,'(f11.3,f10.2,d16.8)')  &
     &   Timeps, Temp, Hamiltonian

     if(QAccAv.and.TotalStep>Nskip) then

       xinv = 1.d0 / dble(TotalStep-Nskip)

       Av_Temp       =  Sum_Temp       * xinv
       Av_Qkinetic   =  Sum_Qkinetic   * xinv
       Av_Potential  =  Sum_Potential  * xinv
       Av_Kinetic    =  Sum_Kinetic    * xinv
       Av_EneSystem  =  Sum_EneSystem  * xinv
       Av_ThermoBath =  Sum_ThermoBath * xinv
       Av_Hamiltonian=  Sum_Hamiltonian* xinv

       write(63,'(f11.3,5d13.5)')                    &
       &  Timeps,  Av_Potential, Av_Kinetic, Av_Qkinetic, Av_EneSystem, Av_ThermoBath

       write(66,'(f11.3,f15.5,d16.8)')  &
       &   Timeps, Av_Temp, Av_Hamiltonian

     end if

   end subroutine Out_data_PIMDiso

   subroutine Average_PIMDiso

   implicit none

       Sum_Temp       =  Sum_Temp        + Temp       
       Sum_Qkinetic   =  Sum_Qkinetic    + Qkinetic   
       Sum_Potential  =  Sum_Potential   + Potential  
       Sum_Kinetic    =  Sum_Kinetic     + Kinetic    
       Sum_EneSystem  =  Sum_EneSystem   + EneSystem  
       Sum_ThermoBath =  Sum_ThermoBath  + ThermoBath 
       Sum_Hamiltonian=  Sum_Hamiltonian + Hamiltonian

   end subroutine Average_PIMDiso

end subroutine PIMD_iso
