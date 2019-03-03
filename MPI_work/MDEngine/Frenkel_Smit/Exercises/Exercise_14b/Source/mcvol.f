      Subroutine Mcvol(En, Vir, Attempt, Acc, Vmax)

C     Attempts To Change The Volume

      Implicit None
      
      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'potential.inc'
      Include 'system.inc'

      Double Precision Enn(2), En, Vir, Virn(2), Yy, Vmax, 
     &     F(2), Arg, Volo(2), Volt, Dlnv, 
     &     Voln(2), Dele1, Dele2, Dlnv1, Dlnv2, 
     &     Enold,Ran_Uniform
      Dimension En(*), Vir(*)
      Integer Attempt, Acc, I, Ib, Idi
 
      Attempt = Attempt + 1

C     ---Calulate New Volume By Making Random Walk In Ln V

      Volo(1) = Box(1)**3
      Volo(2) = Box(2)**3
      Volt = Volo(1) + Volo(2)
      Dlnv = Log(Volo(1)/Volo(2)) + (Ran_Uniform()-0.5d0)*Vmax
      Voln(1) = Exp(Dlnv)*Volt/(1.0d0+Exp(Dlnv))
      Voln(2) = Volt - Voln(1)

      Do Ib = 1, 2
         Box(Ib) = Voln(Ib)**(1.D0/3.D0)
         F(Ib) = Box(Ib)/Volo(Ib)**(1d0/3d0)
         Hbox(Ib) = Box(Ib)/2
         Rc(Ib) = F(Ib)*Rc(Ib)
         Rc2(Ib) = Rc(Ib)**2
      End Do

C     ---Determine New Coordinates

      Do I = 1, Npart
         Idi = Id(I)
         X(I) = F(Idi)*X(I)
         Y(I) = F(Idi)*Y(I)
         Z(I) = F(Idi)*Z(I)
      End Do
      
C        ---Calculate New Energy Using Scaling

      Do Ib = 1, 2
         Enold = En(Ib)
         Yy = (Volo(Ib)/Voln(Ib))**2
         Enn(Ib) = Enold*Yy*(2.0d0-Yy) - 
     &        Vir(Ib)*Yy*(1.0d0-Yy)/6.0d0
         Virn(Ib) = -12.0d0*Enold*Yy*(Yy-1.0d0) + 
     &        Vir(Ib)*Yy*(2.0d0*Yy-1.0d0)
         
      End Do
      
C     ---Acceptance:

      Dele1 = Enn(1) - En(1)
      Dele2 = Enn(2) - En(2)
      Dlnv1 = Log(Voln(1)/Volo(1))
      Dlnv2 = Log(Voln(2)/Volo(2))
      Arg = Exp(-Beta*(Dele1+Dele2-Dble((Npbox(1)+1))*
     &     Dlnv1/Beta-Dble((Npbox(2)+1))
     &      *Dlnv2/Beta))
      If (Ran_Uniform().Lt.Arg) Then

C        ---Accepted

         Acc = Acc + 1
         Do Ib = 1, 2
            En(Ib) = Enn(Ib)
            Vir(Ib) = Virn(Ib)
         End Do
      Else

C        ---Restore The Old Configuration

         Do Ib = 1, 2
            F(Ib) = 1/F(Ib)
            Box(Ib) = Box(Ib)*F(Ib)
            Hbox(Ib) = 0.5d0*Box(Ib)
            Rc(Ib) = F(Ib)*Rc(Ib)
            Rc2(Ib) = Rc(Ib)**2
         End Do
         Do I = 1, Npart
            Idi = Id(I)
            X(I) = F(Idi)*X(I)
            Y(I) = F(Idi)*Y(I)
            Z(I) = F(Idi)*Z(I)
         End Do
      End If

      Return
      End
