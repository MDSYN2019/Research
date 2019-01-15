      Subroutine Repulsion(Dx,Dy,Dz,U)
      Implicit None

      Include 'system.inc'

      Double Precision R,Dx,Dy,Dz,U

      R = Dsqrt(Dx*Dx + Dy*Dy + Dz*Dz)
      U = 0.0d0

      If(R.Lt.Rcut) U = Prepot*((R-Rcut)**2)

      Return
      End
