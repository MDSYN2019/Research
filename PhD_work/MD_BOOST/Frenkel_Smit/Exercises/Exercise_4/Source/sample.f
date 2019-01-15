      Subroutine Sample(Switch,Ip)
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Sample Rms Displacement                                  C
C                                                              C
C     Tmax  = Maximum Timesteps For The Correlation Time       C
C     T0max = Maximum Number Of Time Origins                   C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer   Tmax,T0max

      Parameter (Tmax=100000)
      Parameter (T0max=1000000)

      Integer Ttel,Tt0(T0max),Switch,I,T,Tvacf,Dt,Rx0(T0max),
     &     Ip

      Double Precision R2t(Tmax),Nvacf(Tmax),Dtime
       
      Save R2t,Nvacf,Rx0,Tvacf,Tt0
       
      If (Switch.Eq.1) Then

Ccccccccccccccccccccccccccccccccccc
C     Initialize Everything       C
Ccccccccccccccccccccccccccccccccccc

         Do I = 1, Tmax
            R2t(I)   = 0.0d0
            Nvacf(I) = 0.0d0
         Enddo

      Elseif(Switch.Eq.2) Then

Cccccccccccccccccccccccccccccccccccc
C     Initialize New Random Walk   C
Cccccccccccccccccccccccccccccccccccc

         Tvacf = 0
         
      Elseif(Switch.Eq.3) Then

         Tvacf     = Tvacf + 1
         Ttel      = Mod(Tvacf-1, T0max) + 1
         Tt0(Ttel) = Tvacf
         Rx0(Ttel) = Ip
         
Cccccccccccccccccccccccccccccccccccccc
C     Loop Over All Time Origins     C
Cccccccccccccccccccccccccccccccccccccc

         Do T = 1, Min(Tvacf, T0max)
            Dt = Tvacf - Tt0(T) + 1

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Only If The Time Difference Is Shorter    C
C     Than The Maximum Correlation Time         C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

            If (Dt.Lt.Tmax) Then
               Nvacf(Dt) = Nvacf(Dt) + 1.0d0
               R2t(Dt)   = R2t(Dt)   + Dble(Ip-Rx0(T))**2  
            Endif
         Enddo
                
      Else

Ccccccccccccccccccccccccccccccccccccc
C     Write Everything To Disk      C
Ccccccccccccccccccccccccccccccccccccc

         Open(25,File="rms.dat")

         Do I = 1, Tmax-1
            Dtime = Dble(I-1)

            If (Nvacf(I).Ge.0.5d0) Then
               Write(25,'(2e20.10)') Dtime,R2t(I)/Nvacf(I)
            Endif
         Enddo

         Close(25)
      Endif

      Return
      End
