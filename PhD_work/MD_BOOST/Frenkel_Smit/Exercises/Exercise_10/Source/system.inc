C     System.Inc
C
C     Information About The System
C
C     Maxpart     : Maximum Number Of Particles
C
C     Fxx/Fyy/Fzz : Forces
C     Rxx/Ryy/Rzz : Positions
C     Rxf/Ryf/Rzf : Old Positions
C     Mxx/Myy/Mzz : Positions That Are Not Put Back In The Box
C     Vxx/Vyy/Vzz : Velocities
C     Box         : Boxlengths
C     Hbox        : Half Of The Box-Length
C
C     Tstep       : Timestep
C     Tstep2      : Tstep*Tstep
C     I2tstep     : 1/(2*Tstep)
C
C     Nstep       : Number Of Integration Steps
C     Ninit       : Number Of Initialization Steps
C     Npart       : Number Of Particles
C
C     Ukin        : Kinetic Energy
C     Upot        : Potential Energy
C     Utot        : Total Energy
C     Temp        : Temperature
C     Press       : Pressure
C
C     Rcutsq      : Cut-Off Radius
C     Ecut        : Cut-Off Energy

      Integer Maxpart

      Parameter (Maxpart = 100)

      Double Precision Fxx(Maxpart),Fyy(Maxpart),Fzz(Maxpart),
     &                 Rxx(Maxpart),Ryy(Maxpart),Rzz(Maxpart),
     &                 Rxf(Maxpart),Ryf(Maxpart),Rzf(Maxpart),
     &                 Vxx(Maxpart),Vyy(Maxpart),Vzz(Maxpart),
     &                 Mxx(Maxpart),Myy(Maxpart),Mzz(Maxpart),
     &                 Tstep,Tstep2,I2tstep,Box,Hbox,Rcutsq,
     &                 Ecut,Ukin,Upot,Utot,Temp,Press

      Integer          Nstep,Npart,Ninit

      Common /Syss/    Fxx,Fyy,Fzz,Rxx,Ryy,Rzz,Rxf,Ryf,Rzf,
     &                 Temp,Vxx,Vyy,Vzz,Tstep,Tstep2,I2tstep,
     &                 Box,Hbox,Ukin,Upot,Mxx,Myy,Mzz,Utot,
     &                 Rcutsq,Ecut,Press,Nstep,Npart,Ninit
