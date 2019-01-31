      Subroutine Readdat(Equil, Prod, Nsamp, Ndispl, Dr, Nvol, Vmax, 
     &                   Nswap, Succ)
      Implicit None

C     Reads Input Data And Model Parameters
C
C     ---Input Parameters: File: Fort.15
C    Ibeg  =  0 : Initilaize From A Lattice
C             1 : Read Configuration From Disk
C    Equil      : Number Of Monte Carlo Cycles During Equilibration
C    Prod       : Number Of Monte Carlo Cycles During Production
C    Nsamp      : Number Of Monte Carlo Cycles Between Two Sampling Periods
C    Dr         : Maximum Displacement
C    Vmax       : Maximum Volume Change
C    Succ       : Optimal Percentance Of Accepted Attemps
C                 The Program Adjusts Vmax Or Dr In Just A Way That
C                 On Average Succ% Of The Moves Are Accepted
C    Ndispl     : Number Of Attemps To Displace A Particle Per Mc Cycle
C    Nvol       : Number Of Attemps To Change The Volume  Per Mc Cycle
C    Nswap      : Number Of Attemps To Swap Particle Between The Two Boxes Per Mc Cycle
C    Npart      : Total Numbero Fo Particles
C    Temp       : Temperature
C    Rho        : Density
C
C     ---Input Parameters: File: Fort.25
C    Eps    = Epsilon Lennard-Jones Potential
C    Sig    = Sigma Lennard-Jones Potential
C    Mass   = Mass Of The Particle
C    Rcc    = Cut-Off Radius Of The Potential
C
C     ---Input Parameters: File: Fort.11 (Restart File
C                To Continue A Simulation From Disk)
C    Box(1)   = Length Box 1 Old Configuration
C    Hbox(1)  = Box(1)/2
C    Box(2)   = Length Box 2 Old Configuration
C    Hbox(2)  = Box(2)/2
C    Npart    = Total Number Of Particles (Over Rules Fort.15!!)
C    Npbox(1) = Number Of Particles In Box 1
C    Npbox(2) = Number Of Particles In Box 2
C    Dr     = Optimized Maximum Displacement Old Configurations
C    Vmax   = Optimized Maximum Volume Change Old Configurations
C    X(1),Y(1),Z(1)            : Position First Particle 1
C                  ,Id(1)   = 1 Particle In Box 1
C                  ,Id(1)   = 2 Particle In Box 2
C       ....
C    X(Npart),Y(Npart),Z(Npart): Position Particle Last Particle

      Include 'parameter.inc'
      Include 'system.inc'
      Include 'potential.inc'
      Include 'conf.inc'

      Integer Ibeg, Equil, Prod, I, Ndispl, Nsamp, Nvol, Ib, 
     &        Nswap
      Double Precision Eps, Sig, Rho, Dr, Vmax, Succ, 
     &                 Rcc
 
      Read (15, *)
      Read (15, *) Ibeg, Equil, Prod, Nsamp
      Read (15, *)
      Read (15, *) Dr, Vmax, Succ
      Read (15, *)
      Read (15, *) Ndispl, Nvol, Nswap
      Read (15, *)
      Read (15, *) Npart, Temp, Rho

      If (Npart.Gt.Npmax) Then
         Write (6, *) ' Error: Number Of Particles Too Large'
         Stop
      Endif

C     ---Read Model Parameters

      Read (25, *)
      Read (25, *) Eps, Sig, Mass, Rcc

C     ---Read/Generate Configuration

      Box(1)  = (Dble(Npart)/(2*Rho))**(1.D00/3.D00)
      Hbox(1) = 0.5d00*Box(1)
      Box(2)  = Box(1)
      Hbox(2) = Hbox(1)

      If (Ibeg.Eq.0) Then

C        ---Generate Configuration Form Lattice

         Call Lattice
      Else
         Write (6, *) ' Read Conf From Disk '
         Read (11, *) Box(1), Hbox(1), Box(2), Hbox(2)
         Read (11, *) Npart, Npbox(1), Npbox(2)
         Read (11, *) Dr, Vmax
         Do I = 1, Npart
            Read (11, *) X(I), Y(I), Z(I), Id(I)
         Enddo
         Rewind (11)
      Endif

C     ---Write Input Data

      Write (6, 99001) Equil, Prod, Nsamp
      Write (6, 99002) Ndispl, Dr, Nvol, Vmax, Nswap
      Write (6, 99003) Npart, Temp, 0.0d0, 
     &                 Dble(Npbox(1))/Box(1)**3, Box(1), 
     &                 Dble(Npbox(2))/Box(2)**3, Box(2)
      Write (6, 99004) Eps, Sig, Mass

C     ---Calculate Parameters:

      Beta  = 1.0d0/Temp
      Eps4  = 4.0d0*Eps
      Eps48 = 48.D0*Eps
      Sig2  = Sig*Sig

C     ---Calculate Cut-Off Radius Potential

      Do Ib = 1, 2
         Rc(Ib) = Min(Rcc, Hbox(Ib))
         Rc2(Ib) = Rc(Ib)*Rc(Ib)
      Enddo

99001 Format ('  Number Of Equilibration Cycles             :', I10, /, 
     &        '  Number Of Production Cycles                :', I10, /, 
     &        '  Sample Frequency                           :', I10, /)

99002 Format ('  Number Of Att. To Displ. A Part. Per Cycle :', I10, /, 
     &        '  Maximum Displacement                       :', F10.3, 
     &        /, '  Number Of Att. To Change Volume  Per Cycle :', I10, 
     &        /, '  Maximum Change Volume                      :', 
     &        F10.3, /, '  Number Of Att. To Exch Part.  Per Cycle    :'
     &        , I10, //)
99003 Format ('  Number Of Particles                        :', I10, /, 
     &        '  Temperature                                :', F10.3, 
     &        /, '  Pressure                                   :', 
     &        F10.3, /, '  Density Box 1                              :'
     &        , F10.3, /, 
     &        '  Box 1 Length                               :', F10.3, 
     &        /, '  Density Box 1                              :', 
     &        F10.3, /, '  Box 2 Length                               :'
     &        , F10.3, /)
99004 Format ('  Model Parameters: ', /, '     Epsilon: ', F5.3, /, 
     &        '     Sigma  : ', F5.3, /, '     Mass   : ', F5.3)
      Return
      End
