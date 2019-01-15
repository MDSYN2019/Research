C     System.Inc
C
C     Information About The System; Various Constants
C
C     Maxchain = Maximum Chain Length
C     Maxtrial = Maximum Number Of Trial Positions
C
C     Lstatic  = Do We Use A Static Scheme ?
C     Lcbmc    = Do We Use Cbmc ?
C     
C     Nuall    = Chain Length
C     Nchoi    = Number Of Trial Positions
C
C     X/Y/Zpos = Coordinates Of The Chain

      Integer          Nstep,Ninit,Maxchain,Maxtrial,Nchoi,
     &                 Nuall
      Logical          Lstatic,Lcbmc

      Parameter(Maxchain = 100)
      Parameter(Maxtrial = 20)

      Double Precision Xpos(Maxchain),Ypos(Maxchain),Zpos(Maxchain),
     &                 Rcut,Prepot,Kb,Thetan,Beta

      Common /Syss/    Xpos,Ypos,Zpos,Rcut,Prepot,Kb,Thetan,Beta,
     &                 Nstep,Nchoi,Nuall,Ninit,Lstatic,Lcbmc  
