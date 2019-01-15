      Subroutine Sample(I, En, Vir)

C     Writes Quantities To File

      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Integer I, Ib
      Double Precision En(*), Enp(2), Vir(*), Press(2), Vol, 
     &                 Rho(2)
 
      Do Ib = 1, 2
         Vol = Box(Ib)**3
         Rho(Ib) = Dble(Npbox(Ib))/Vol
         Press(Ib) = Rho(Ib)/Beta + Vir(Ib)/(3.D0*Vol)
         
         If (Npbox(Ib).Ne.0) Then
            Enp(Ib) = En(Ib)/Dble(Npbox(Ib))
         Else
            Enp(Ib) = 0.D0
         End If
      End Do
      Write (66,*) I, Dble(Enp(1)), Dble(Enp(2)), Dble(Press(1)), 
     &              Dble(Press(2)), Dble(Rho(1)), Dble(Rho(2))
      Write (44,'(2(I6,F10.2))') Npbox(1), Box(1)**3, Npbox(2), Box(2)
     &                            **3
      Write(45,'(3e20.10)') Dble(I),(Dble(Npbox(1))/Box(1)**3),
     &     (Dble(Npbox(2))/Box(2)**3)
      Return
      End
