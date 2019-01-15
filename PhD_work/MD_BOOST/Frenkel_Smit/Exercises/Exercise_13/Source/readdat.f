      Subroutine Readdat
      Implicit None
 
      Include 'system.inc'
 
Cccccccccccccccccccccccccccccccccccccc
C     Read In System Information     C
Cccccccccccccccccccccccccccccccccccccc
 
      Double Precision Ran_Vel,U
      Integer I

      Read (21,*) Nstep,Temp,Tstep,Ninit
      Read (21,*) Choise,Nu,Freqnos,Nnoshover
 
      Onepi = 4.0d0*Datan(1.0d0)
      Xpos  = 0.0d0
      Vpos  = Ran_Vel()
      Oldf  = 0.0d0
      U     = 0.0d0

      Call Force(Xpos,U,Oldf)

      Wdti2(1) = 1.0d0/(4.0d0 - (4.0d0**(1.0d0/3.0d0)))
      Wdti2(2) = Wdti2(1)
      Wdti2(3) = 1.0d0 - 4.0d0*Wdti2(1)
      Wdti2(4) = Wdti2(1)
      Wdti2(5) = Wdti2(1)
 
      If(Choise.Eq.3) Vpos = Dsqrt(Temp)

      Do I = 1,5
         Wdti2(I) = Wdti2(I)*Tstep*0.5d0
         Wdti4(I) = Wdti2(I)*0.5d0
         Wdti8(I) = Wdti4(I)*0.5d0
      Enddo
 
      Do I = 1,Nnoshover
         Qmass(I)  = Temp/(Freqnos*Freqnos)
         Iqmass(I) = 1.0d0/Qmass(I)
         Xlogs(I)  = 0.0d0
         Glogs(I)  = 0.0d0
         Vlogs(I)  = Iqmass(I)
      Enddo

      Do I = 2,Nnoshover
         Glogs(I) = Iqmass(I)*(Qmass(I-1)*Vlogs(I-1)*Vlogs(I-1) - Temp)
      Enddo

      Write (6,*) 'Md Of A Particle Over A Energy Barrier'
      Write (6,*)
      Write (6,*) 'Number Of Steps       : ',Nstep
      Write (6,*) 'Number Of Init Steps  : ',Ninit
      Write (6,*) 'Timestep              : ',Tstep
      Write (6,*) 'Temperature           : ',Temp

      If(Choise.Eq.2) Write (6,*) 'Collision Frequency   : ',Nu
      If(Choise.Eq.3) Write (6,*) '# Nose-Hover Chains   : ',Nnoshover
      If(Choise.Eq.3) Write (6,*) 'Frequency Omega       : ',Freqnos

      Write (6,*)
       
      Return
      End
