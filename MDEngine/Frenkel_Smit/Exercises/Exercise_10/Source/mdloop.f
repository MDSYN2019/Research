      Subroutine Mdloop
      Implicit None
 
      Include 'system.inc'
 
      Integer I
      Double Precision Av(6),Impx,Impy,Impz,Tempz
 
      Do I = 1,6
         Av(I) = 0.0d0
      Enddo
 
C     Initialize Radial Distribution Function

      Call Sample_Gyra(1)

C     Initialize Diffusion Coefficient

      Call Sample_Diff(1)

C     Loop Over All Cycles

      Do I = 1,Nstep
 
C     Calculate The Force
 
         Call Force
 
C     Integrate The Equations Of Motion
 
         Call Integrate(I,Impx,Impy,Impz)
 
         If (I.Eq.1) Then
            Write (6,*)
            Write (6,*) 'Impulse X-Dir.        : ',Impx
            Write (6,*) 'Impulse Y-Dir.        : ',Impy
            Write (6,*) 'Impulse Z-Dir.        : ',Impz
            Write (6,*)
            Write (6,*)
            Write (6,'(3a)') ' Step         Utot         Ukin',
     &           '         Upot         Temp',
     &           '        Press'
            Write (6,*)
         Endif
 
         Utot = Ukin + Upot
 
         Tempz = 2.0d0*Ukin/Dble(3*Npart - 3)
 
         If (I.Gt.Ninit.And.Mod(I,50).Eq.0) Write (6,'(I5,8(1x,E12.5))') 
     &        I,Utot,Ukin,Upot,Tempz,Press

         If (Mod(I,100).Eq.0) Call Writepdb
 
         If (I.Eq.Nstep) Then
            Write (6,*)
            Write (6,*) 'Impulse X-Dir.        : ',Impx
            Write (6,*) 'Impulse Y-Dir.        : ',Impy
            Write (6,*) 'Impulse Z-Dir.        : ',Impz
            Write (6,*)
         Endif
 
 
C     Calculate Averages
 
         If (I.Gt.Ninit) Then
            Av(1) = Av(1) + Tempz
            Av(2) = Av(2) + Press
            Av(3) = Av(3) + Ukin
            Av(4) = Av(4) + Upot
            Av(5) = Av(5) + Utot
            Av(6) = Av(6) + 1.0d0

C     Sample Radial Distribution Function

            If(Mod(I,100).Eq.0) Call Sample_Gyra(2)

            Call Sample_Diff(2)
         Endif
      Enddo
 
C     Print Averages To Screen
 
      If (Av(6).Gt.0.5d0) Then
         Av(6) = 1.0d0/Av(6)
         Do I = 1,5
            Av(I) = Av(I)*Av(6)
         Enddo
 
         Write (6,*)
         Write (6,*) 'Average Temperature   : ',Av(1)
         Write (6,*) 'Average Pressure      : ',Av(2)
         Write (6,*) 'Average Ukin          : ',Av(3)
         Write (6,*) 'Average Upot          : ',Av(4)
         Write (6,*) 'Average Utot          : ',Av(5)
 
      Endif

      Call Sample_Gyra(3)
      Call Sample_Diff(3)

      Return
      End
