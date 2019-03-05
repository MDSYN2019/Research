module Forplot

integer :: Nfile, Nst, Nskip
character(len=80), dimension(:), allocatable :: Energy_File
character(len=4) :: Method
character(len=3) :: Ensemble
logical :: QAccAv

real(8) :: Hamiltonian, Timeps, Temp
real(8), dimension(3,3) :: Pressure, H

real(8) :: Kinetic, Potential, Poten_co, EneSystem
real(8) :: E_BondMol, E_AngleMol, E_UBMol, E_DihedMol, E_ImproMol
real(8) :: E_LJMol, E_LJcoMol, E_ErMol, E_EkMol, E_EsMol, E_OptC
real(8) :: E_ElecMol, E_Int_Mol, E_Ext_Mol
real(8) :: ThermoBath, Pot_p
real(8) :: density, Area, Volume
real(8) :: LenA, LenB, LenC, csAB, csBC, csCA
real(8) :: AngAB, AngBC, AngCA
real(8) :: Accept_Ratio
real(8) :: Qkinetic

real(8) :: Sum_Temp, Sum_E_BondMol, Sum_E_AngleMol, Sum_E_UBMol, Sum_E_DihedMol
real(8) :: Sum_E_ImproMol, Sum_E_OptC, Sum_E_LJMol, Sum_E_ErMol, Sum_E_EkMol
real(8) :: Sum_E_EsMol, Sum_E_ElecMol, Sum_E_Int_Mol, Sum_E_Ext_Mol, Sum_Potential
real(8) :: Sum_Kinetic, Sum_EneSystem, Sum_ThermoBath, Sum_Pot_p, Sum_Hamiltonian
real(8) :: Sum_Area, Sum_Volume, Sum_density, Sum_LenA, Sum_LenB, Sum_LenC
real(8) :: Sum_AngBC, Sum_AngCA, Sum_AngAB
real(8), dimension(3,3) :: Sum_H, Sum_Pressure

real(8) :: Sum_Accept_Ratio
real(8) :: Sum_Qkinetic

real(8) :: Av_Temp, Av_E_BondMol, Av_E_AngleMol, Av_E_UBMol, Av_E_DihedMol
real(8) :: Av_E_ImproMol, Av_E_OptC, Av_E_LJMol, Av_E_ErMol, Av_E_EkMol
real(8) :: Av_E_EsMol, Av_E_ElecMol, Av_E_Int_Mol, Av_E_Ext_Mol, Av_Potential
real(8) :: Av_Kinetic, Av_EneSystem, Av_ThermoBath, Av_Pot_p, Av_Hamiltonian
real(8) :: Av_Area, Av_Volume, Av_density, Av_LenA, Av_LenB, Av_LenC
real(8) :: Av_AngBC, Av_AngCA, Av_AngAB
real(8), dimension(3,3) :: Av_H, Av_Pressure

real(8) :: Av_Accept_Ratio
real(8) :: Av_Qkinetic

#ifdef SURF
real(8) :: SurfTens
real(8) :: Sum_SurfTens
real(8) :: Av_SurfTens
#endif

end module Forplot

