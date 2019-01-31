      Subroutine Lattice
      Implicit None

C     ---Place `Npart' Particles On A Lattice With Density 'Rho'
C      --Half The Number In Box 1 And The Other Half In Box 2
      
      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Integer I, J, K, Itel, N, Ib
      Double Precision Del
 
      Del      = (Box(1)**3)**(1.D0/3.D0)
      Npbox(1) = Npart/2
      Npbox(2) = Npbox(1)

      If (Npbox(1)+Npbox(2).Ne.Npart) Then
         Stop 'Error Npart'
      End If

      Write (6, *) ' Generate Simple Cubic Lattice'
      
      N = Int(Dble(Npart)**(1.D0/3.D0)) + 1
      If (N.Eq.0) N = 1
      Del = Del/Dble(N)
      Itel = 0
      Do I = 0, N - 1
         Do J = 0, N - 1
            Do K = 0, N - 1
               Do Ib = 1, 2
                  If (Itel.Lt.Npart) Then
                     Itel = Itel + 1
                     X(Itel) = Dble(K)*Del
                     Y(Itel) = Dble(J)*Del
                     Z(Itel) = Dble(I)*Del
                     Id(Itel) = Ib
                  End If
               End Do
            End Do
         End Do
      End Do
      
      Write (6, 99001) Itel
      
99001 Format (' Initialisation On Lattice: ', /, I10, 
     &        ' Particles Placed On A Lattice')
      
      Return
      End
