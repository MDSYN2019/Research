      Subroutine Angle(Dx1,Dy1,Dz1,Dx2,Dy2,Dz2,U)
      Implicit None

      Include 'system.inc'

      Double Precision R,U,Dx1,Dx2,Dy1,Dy2,Dz1,Dz2

      R =        (Dx1*Dx2 + Dy1*Dy2 + Dz1*Dz2)/
     &     Dsqrt((Dx1*Dx1 + Dy1*Dy1 + Dz1*Dz1)*
     &           (Dx2*Dx2 + Dy2*Dy2 + Dz2*Dz2))

      U = 0.5d0*Kb*((Dacos(R)-Thetan)**2)

      Return
      End
