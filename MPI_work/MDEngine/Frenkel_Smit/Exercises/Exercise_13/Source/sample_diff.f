      Subroutine Sample_Diff(Switch)
      Implicit None

      Include 'system.inc'

C     Tmax  = Maximum Timesteps For The Correlation Time
C     T0max = Maximum Number Of Time Origins 

      Integer   Tmax,T0max

      Parameter (Tmax=1000)
      Parameter (T0max=200)

      Integer Ttel,Tt0(T0max),T0,Switch,I,T,Tvacf,Dt

      Double Precision R2t(Tmax),Rx0(T0max),Nvacf(Tmax),Dtime
       
      Save Nvacf,Tvacf,R2t,Rx0,Tt0,T0
       
      If (Switch.Eq.1) Then

         Tvacf = 0
         T0    = 0

         Do I = 1, Tmax
            R2t(I)   = 0.0d0
            Nvacf(I) = 0.0d0
         Enddo

      Elseif(Switch.Eq.2) Then

         Tvacf = Tvacf + 1

         If (Mod(Tvacf,50).Eq.0) Then
            T0         = T0 + 1
            Ttel       = Mod(T0-1, T0max) + 1
            Tt0(Ttel)  = Tvacf
            Rx0(Ttel)  = Xpos
         Endif

         Do T = 1, Min(T0, T0max)
            Dt = Tvacf - Tt0(T) + 1

            If (Dt.Lt.Tmax) Then
               Nvacf(Dt) = Nvacf(Dt) + 1.0d0
               R2t(Dt)   = R2t(Dt)   + (Xpos-Rx0(T))**2
            Endif
         Enddo
                
      Else

         Do I = 1, Tmax-1
                           
            Dtime = Dble(I-1)*Tstep

            If (Nvacf(I).Ge.0.5d0) Then
               R2t(I)  = R2t(I) /Nvacf(I)
            Else
               R2t(I)  = 0.0d0
            Endif

            If(I.Ne.1) Then
               Write(25,*) Dtime,R2t(I),R2t(I)/(2.0d0*Dtime)
            Else
               Write(25,*) Dtime,R2t(I),' 0.0d0'
            Endif
         Enddo
      Endif

      Return
      End
