      Subroutine Readdat
      Implicit None
 
      Include 'system.inc'
 
C     Read In System Information
 
      Read (21,*) Box,Npart,Nstep,Temp,Tstep,Ninit
 
      If (Npart.Gt.Maxpart) Then
         Write (6,*) 'Maximum No. Of Particles Is : ',Maxpart
         Stop
      Endif
 
C     Calculate Some Parameters
 
      Hbox = 0.5d0*Box
 
      Rcutsq = (0.49999d0*Box)**2
      Ecut = 4.0d0*((Rcutsq**(-6.0d0)) - (Rcutsq**(-3.0d0)))
 
      Tstep2 = Tstep*Tstep
      I2tstep = 0.5d0/Tstep
 
C     Print Information To The Screen
 
      Write (6,*) 'Molecular Dynamics Program'
      Write (6,*)
      Write (6,*) 'Number Of Particles   : ',Npart
      Write (6,*) 'Boxlength             : ',Box
      Write (6,*) 'Density               : ',Dble(Npart)/(Box*Box*Box)
      Write (6,*) 'Temperature           : ',Temp
      Write (6,*) 'Cut-Off Radius        : ',Dsqrt(Rcutsq)
      Write (6,*) 'Cut-Off Energy        : ',Ecut
      Write (6,*) 'Number Of Steps       : ',Nstep
      Write (6,*) 'Number Of Init Steps  : ',Ninit
      Write (6,*) 'Timestep              : ',Tstep
 
      Return
      End
