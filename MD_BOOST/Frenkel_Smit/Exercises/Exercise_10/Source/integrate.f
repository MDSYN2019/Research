      Subroutine Integrate(Step,Impx,Impy,Impz)
      Implicit None
 
      Include 'system.inc'
 
C     Integrate The Equations Of Motion And Calculate The Total Impulse
 
      Integer I,Step
      Double Precision Xxx(Maxpart),Yyy(Maxpart),Zzz(Maxpart),
     &                 Scale,Impx,Impy,Impz

Ccccccccccccccccccccccccccccccccccccccccccccccccccc 
C     Set Kinetic Energy To Zero                  C
C     Verlet Integrator                           C
C                                                 C
C     Xxx/Yyy/Zzz = New Position (T + Delta T)    C
C     Rxx/Ryy/Rzz = Position     (T          )    C
C     Rxf/Ryf/Rzf = Old Position (T - Delta T)    C
C                                                 C
C     Vxx/Vyy/Vzz = Velocity     (T          )    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccc

      Ukin = 0.0d0
 
      Do I = 1,Npart
         Xxx(I) = 2.0d0*Rxx(I) - Rxf(I) + Fxx(I)*Tstep2
         Yyy(I) = 2.0d0*Ryy(I) - Ryf(I) + Fyy(I)*Tstep2
         Zzz(I) = 2.0d0*Rzz(I) - Rzf(I) + Fzz(I)*Tstep2
 
         Vxx(I) = (Xxx(I) - Rxf(I))*I2tstep
         Vyy(I) = (Yyy(I) - Ryf(I))*I2tstep
         Vzz(I) = (Zzz(I) - Rzf(I))*I2tstep
 
         Ukin = Ukin + 0.5d0*(Vxx(I)*Vxx(I) + 
     &        Vyy(I)*Vyy(I) + Vzz(I)*Vzz(I))
      Enddo
 
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     For Nstep < Ninit; Use The Velocity Scaling To Get The Exact Temperature  C
C     Otherwise: Scale = 1. Beware That The Positions/Velocities Have To Be     C
C     Recalculated !!!!                                                         C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      If (Step.Le.Ninit) Then
         Scale = Dsqrt(Temp*Dble(3*Npart-3)/(2.0d0*Ukin))
      Else
         Scale = 1.0d0
      Endif
 
      Ukin = 0.0d0
      Impx = 0.0d0
      Impy = 0.0d0
      Impz = 0.0d0
 
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Scale Velocities And Put Particles Back In The Box      C
C     Beware: The Old Positions Are Also Put Back In The Box  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      Do I = 1,Npart
         Vxx(I) = Scale*Vxx(I)
         Vyy(I) = Scale*Vyy(I)
         Vzz(I) = Scale*Vzz(I)
 
         Impx = Impx + Vxx(I)
         Impy = Impy + Vyy(I)
         Impz = Impz + Vzz(I)
 
         Xxx(I) = Rxf(I) + 2.0d0*Vxx(I)*Tstep
         Yyy(I) = Ryf(I) + 2.0d0*Vzz(I)*Tstep
         Zzz(I) = Rzf(I) + 2.0d0*Vzz(I)*Tstep

         Mxx(I) = Mxx(I) + Xxx(I) - Rxx(I)
         Myy(I) = Myy(I) + Yyy(I) - Ryy(I)
         Mzz(I) = Mzz(I) + Zzz(I) - Rzz(I)
 
         Ukin = Ukin + Vxx(I)*Vxx(I) + Vyy(I)*Vyy(I) + Vzz(I)*Vzz(I)
 
         Rxf(I) = Rxx(I)
         Ryf(I) = Ryy(I)
         Rzf(I) = Rzz(I)
 
         Rxx(I) = Xxx(I)
         Ryy(I) = Yyy(I)
         Rzz(I) = Zzz(I)

C     Put Particles Back In The Box
C     The Previous Position Has The Same
C     Transformation (Why ??)

         If (Rxx(I).Gt.Box) Then
            Rxx(I) = Rxx(I) - Box
            Rxf(I) = Rxf(I) - Box
         Elseif (Rxx(I).Lt.0.0d0) Then
            Rxx(I) = Rxx(I) + Box
            Rxf(I) = Rxf(I) + Box
         Endif
 
         If (Ryy(I).Gt.Box) Then
            Ryy(I) = Ryy(I) - Box
            Ryf(I) = Ryf(I) - Box
         Elseif (Ryy(I).Lt.0.0d0) Then
            Ryy(I) = Ryy(I) + Box
            Ryf(I) = Ryf(I) + Box
         Endif
 
         If (Rzz(I).Gt.Box) Then
            Rzz(I) = Rzz(I) - Box
            Rzf(I) = Rzf(I) - Box
         Elseif (Rzz(I).Lt.0.0d0) Then
            Rzz(I) = Rzz(I) + Box
            Rzf(I) = Rzf(I) + Box
         Endif
      Enddo
 
Ccccccccccccccccccccccccccccccccccccccccccccc
C     Add The Kinetic Part Of The Pressure  C
Ccccccccccccccccccccccccccccccccccccccccccccc
 
      Press = Press + 2.0d0*Ukin*Dble(Npart)/
     &     (Box*Box*Box*Dble(3*Npart-3))
 
      Return
      End
