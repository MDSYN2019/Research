      Subroutine Init
      Implicit None
 
      Include 'system.inc'
 
C     Generates Initial Positions/Velocities
C     This Is Not To Easy; Do Not Look For Errors Here !!!!
 
      Integer I,J,K,Number,Nplace
 
      Double Precision Ran_Uniform,Rangauss,Fxo(Maxpart),
     &                 Fyo(Maxpart),Fzo(Maxpart),
     &                 Uold,Testje,Place,Size,Impx,Impy,Impz
 
C     Generate Velocities From A Gaussian; Set Impulse To Zero
 
      Impx = 0.0d0
      Impy = 0.0d0
      Impz = 0.0d0
      Ukin = 0.0d0
 
      Do I = 1,Npart
         Vxx(I) = Rangauss()
         Vyy(I) = Rangauss()
         Vzz(I) = Rangauss()
         Impx = Impx + Vxx(I)
         Impy = Impy + Vyy(I)
         Impz = Impz + Vzz(I)
      Enddo
 
      Impx = Impx/Dble(Npart)
      Impy = Impy/Dble(Npart)
      Impz = Impz/Dble(Npart)

C     Calculate The Kinetic Energy Ukin
 
      Do I = 1,Npart
         Vxx(I) = Vxx(I) - Impx
         Vyy(I) = Vyy(I) - Impy
         Vzz(I) = Vzz(I) - Impz
 
         Ukin = Ukin + Vxx(I)*Vxx(I) + 
     &        Vyy(I)*Vyy(I) + Vzz(I)*Vzz(I)
      Enddo
 
C     Scale All Velocities To The Correct Temperature
 
      Testje = Dsqrt(Temp*Dble(3*Npart-3)/Ukin)
 
      Do I = 1,Npart
         Vxx(I) = Testje*Vxx(I)
         Vyy(I) = Testje*Vyy(I)
         Vzz(I) = Testje*Vzz(I)
      Enddo
 
C     Calculate Initial Positions On A Lattice
 
      Number = Idint((Dble(Npart)**(1.0d0/3.0d0)) + 1.5d0)
      Nplace = 0
 
      Size = Box/Dble(Number + 2)
 
      Place = 0.2d0*Size
 
      Do I = 1,Number
         Do J = 1,Number
            Do K = 1,Number
               Nplace = Nplace + 1
               If (Nplace.Le.Npart) Then
                  Rxx(Nplace) = (Dble(I) + 
     &                 0.01d0*(Ran_Uniform()-0.5d0))*Size
                  Ryy(Nplace) = (Dble(J) + 
     &                 0.01d0*(Ran_Uniform()-0.5d0))*Size
                  Rzz(Nplace) = (Dble(K) + 
     &                 0.01d0*(Ran_Uniform()-0.5d0))*Size
               Endif
            Enddo
         Enddo
      Enddo
 
C     Calculate Better Positions Using A Steepest Decent Algorithm
 
      Do J = 1,50
 
         If (J.Eq.1) Then
            Call Force
            Uold = Upot
            Write (6,*) 'Initial Energy        : ',Uold
         Endif
 
         Testje = 0.0d0
 
C     Calculate Maximum Downhill Gradient
 
         Do I = 1,Npart
            Rxf(I) = Rxx(I)
            Ryf(I) = Ryy(I)
            Rzf(I) = Rzz(I)
 
            Fxo(I) = Fxx(I)
            Fyo(I) = Fyy(I)
            Fzo(I) = Fzz(I)
 
            Impx = Dabs(Fxx(I))
            Impy = Dabs(Fyy(I))
            Impz = Dabs(Fzz(I))
 
            If (Impx.Gt.Testje) Testje = Impx
            If (Impy.Gt.Testje) Testje = Impy
            If (Impz.Gt.Testje) Testje = Impz
         Enddo
 
         Testje = Place/Testje

C     Calculate Improved Positions
 
         Do I = 1,Npart
 
            Rxx(I) = Rxx(I) + Testje*Fxx(I)
            Ryy(I) = Ryy(I) + Testje*Fyy(I)
            Rzz(I) = Rzz(I) + Testje*Fzz(I)
 
C     Place Particles Back In The Box
 
            If (Rxx(I).Gt.Box) Then
               Rxx(I) = Rxx(I) - Box
            Elseif (Rxx(I).Lt.0.0d0) Then
               Rxx(I) = Rxx(I) + Box
            Endif
 
            If (Ryy(I).Gt.Box) Then
               Ryy(I) = Ryy(I) - Box
            Elseif (Ryy(I).Lt.0.0d0) Then
               Ryy(I) = Ryy(I) + Box
            Endif
 
            If (Rzz(I).Gt.Box) Then
               Rzz(I) = Rzz(I) - Box
            Elseif (Rzz(I).Lt.0.0d0) Then
               Rzz(I) = Rzz(I) + Box
            Endif
         Enddo
 
C     Calculate New Potential Energy
 
         Call Force
 
C     Check If The New Positions Are Acceptable
 
         If (Upot.Lt.Uold) Then
            Uold  = Upot
            Place = Place*1.2d0
 
            If (Place.Gt.Hbox) Place = Hbox
 
         Else
            Do I = 1,Npart
               Fxx(I) = Fxo(I)
               Fyy(I) = Fyo(I)
               Fzz(I) = Fzo(I)
               Rxx(I) = Rxf(I)
               Ryy(I) = Ryf(I)
               Rzz(I) = Rzf(I)
            Enddo
            Place = Place*0.1d0
         Endif
      Enddo
 
C     Calculate Previous Position Using The Generated Velocity
 
      Do I = 1,Npart
         Rxf(I) = Rxx(I) - Tstep*Vxx(I)
         Ryf(I) = Ryy(I) - Tstep*Vyy(I)
         Rzf(I) = Rzz(I) - Tstep*Vzz(I)
         Mxx(I) = Rxx(I)
         Myy(I) = Ryy(I)
         Mzz(I) = Rzz(I)
      Enddo
 
      Write (6,*) 'Final Energy          : ',Uold
 
      Return
      End
