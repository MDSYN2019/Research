      Subroutine Sample_Diff(Switch)
      Implicit None

      Include 'system.inc'

C     Tmax  = Maximum Timesteps For The Correlation Time
C     T0max = Maximum Number Of Time Origins 

      Integer   Tmax,T0max

      Parameter (Tmax=1000)
      Parameter (T0max=200)

      Integer Ttel,Tt0(T0max),T0,Switch,I,T,Tvacf,Dt

      Double Precision Vxt0(Maxpart,T0max),Vyt0(Maxpart,T0max),
     &                 Vzt0(Maxpart,T0max),Vacf(Tmax),R2t(Tmax),
     &                 Rx0(Maxpart,T0max),Ry0(Maxpart,T0max),
     &                 Rz0(Maxpart,T0max),Nvacf(Tmax),Intt,Dtime
       
      Save Vacf,Nvacf,Vxt0,Vyt0,Vzt0,Tvacf,R2t,Rx0,
     &     Ry0,Rz0,Tt0,T0
       
      If (Switch.Eq.1) Then

C     Initialize Everything

         Tvacf = 0
         T0    = 0

         Do I = 1, Tmax
            R2t(I)   = 0.0d0
            Nvacf(I) = 0.0d0
            Vacf(I)  = 0.0d0
         Enddo

      Elseif(Switch.Eq.2) Then

         Tvacf = Tvacf + 1

         If (Mod(Tvacf,50).Eq.0) Then

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     New Time Origin                                             C
C     Store The Positions/Velocities; The Current Velocities      C
C     Are Vxx(I) And The Current Positions Are Mxx(I).            C
C     Question: Why Do You Have To Be Careful With Pbc ?          C
C                                                                 C
C     Make Sure To Study Algorithm 8 (Page 82) Of Frenkel/Smit    C
C     Before You Start To Make Modifications                      C
C                                                                 C
C     Note That Most Variable Names Are Different Here Then In    C
C     Frenkel/Smit. In This Way, You Will Have To Think More...   C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification

C     End   Modification

         Endif

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Loop Over All Time Origins That Have Been Stored    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification

C     End   Modification

      Else

Ccccccccccccccccccccccccccccccccccc
C     Write Everything To Disk    C
Ccccccccccccccccccccccccccccccccccc

         Intt = 0.0d0

         Do I = 1, Tmax-1
                           
            Dtime = Dble(I-1)*Tstep

            If (Nvacf(I).Ge.0.5d0) Then
               Vacf(I) = Vacf(I)/(Dble(Npart)*Nvacf(I))
               R2t(I)  = R2t(I) /(Dble(Npart)*Nvacf(I))
            Else
               Vacf(I) = 0.0d0
               R2t(I)  = 0.0d0
            Endif

            Intt = Intt + Tstep*Vacf(I)/3.0d0

            Write(24,*) Dtime,Vacf(I),Intt
            
            If(I.Ne.1) Then
               Write(25,*) Dtime,R2t(I),R2t(I)/(6.0d0*Dtime)
            Else
               Write(25,*) Dtime,R2t(I),' 0.0d0'
            Endif
         Enddo
      Endif

      Return
      End
