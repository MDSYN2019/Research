C     System.Inc
C
C     Information About The System; Various Constants
C

      Integer          Nstep,Ninit,Maxhover,Nnoshover,Choise

      Parameter        (Maxhover = 100)

      Double Precision Onepi,Tstep,Nu,Temp,Xpos,Vpos,Cons,
     &                 Wdti2(5),Wdti4(5),Wdti8(5),
     &                 Xlogs(Maxhover),Vlogs(Maxhover),
     &                 Glogs(Maxhover),Qmass(Maxhover),
     &                 Iqmass(Maxhover),Freqnos,Oldf

      Common /Syss/    Onepi,Tstep,Nu,Temp,Xpos,Vpos,Oldf,
     &                 Cons,Wdti2,Wdti4,Wdti8,Xlogs,Vlogs,
     &                 Glogs,Qmass,Iqmass,Freqnos,Nstep,
     &                 Ninit,Nnoshover,Choise
