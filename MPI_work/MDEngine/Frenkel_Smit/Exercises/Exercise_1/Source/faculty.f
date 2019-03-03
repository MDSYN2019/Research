      Function Faculty(N)
      Implicit None

      Double Precision Faculty
      Integer N,J

      Faculty = 0.0d0

      Do J=2,N
         Faculty = Faculty + Dlog(Dble(J))
      Enddo
      
      Return
      End
