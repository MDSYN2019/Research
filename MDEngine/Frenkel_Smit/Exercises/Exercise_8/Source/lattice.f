      Subroutine Lattice
C
C     Place `Npart' Particles On A Simple Cubic
C     Lattice With Density 'Rho'
C

      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Integer I, J, K, Itel, N
      Double Precision Dx, Dy, Dz, Del
 
      N = Int(Npart**(1.0d0/3.0d0)) + 1
      If (N.Eq.0) N = 1
      Del = Box/Dble(N)
      Itel = 0
      Dx = -Del
      Do I = 1, N
         Dx = Dx + Del
         Dy = -Del
         Do J = 1, N
            Dy = Dy + Del
            Dz = -Del
            Do K = 1, N
               Dz = Dz + Del
               If (Itel.Lt.Npart) Then
                  Itel = Itel + 1
                  X(Itel) = Dx
                  Y(Itel) = Dy
                  Z(Itel) = Dz
               End If
            End Do
         End Do
      End Do
      Write (6, 99001) Itel
      Return
99001 Format (' Initialisation On Lattice: ', /, I10, 
     &        ' Particles Placed On A Lattice')
      End
