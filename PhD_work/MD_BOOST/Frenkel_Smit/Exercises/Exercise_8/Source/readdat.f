      Subroutine Readdat(Equil, Prod, Nsamp, Ndispl, Dr)
      Implicit None

C     ---Read Input Data And Model Parameters
C
C     ---Input Parameters: File: Fort.15
C    Ibeg  =  0 : Initialize From A Lattice
C             1 : Read Configuration From Disk
C    Equil      : Number Of Monte Carlo Cycles During Equilibration
C    Prod       : Number Of Monte Carlo Cycles During Production
C    Nsamp      : Number Of Monte Carlo Cycles Between Two Sampling Periods
C    Dr         : Maximum Displacement
C    Ndispl     : Number Of Attemps To Displace A Particle Per Mc Cycle
C    Npart      : Total Numbero Fo Particles
C    Temp       : Temperature
C    Rho        : Density
C
C     ---Input Parameters: File: Fort.25
C    Eps    = Epsilon Lennard-Jones Potential
C    Sig    = Sigma Lennard-Jones Potential
C    Mass   = Mass Of The Particle
C    Rc     = Cut-Off Radius Of The Potential
C
C     ---Input Parameters: File: Fort.11 (Restart File
C                To Continue A Simulation From Disk)
C    Boxf   = Box Length Old Configuration (If This One
C             Does Not Correspond To The Requested Density, The Positions
C             Of The Particles Are Rescaled!
C    Npart  = Number Of Particles (Over Rules Fort.15!!)
C    Dr     = Optimized Maximum Displacement Old Configurations
C
C
C    X(1),Y(1),Z(1)            : Position First Particle 1
C        ...
C    X(Npart),Y(Npart),Z(Npart): Position Particle Last Particle
 
      Include 'parameter.inc'
      Include 'system.inc'
      Include 'potential.inc'
      Include 'conf.inc'

      Integer Ibeg, Equil, Prod, I, Ndispl, Nsamp, Sstmm
      Double Precision Eps, Sig, Boxf, Rhof, Rho, Dr, M1
 
C     ---Read Simulation Data
      Read (15, *)
      Read (15, *) Ibeg, Equil, Prod, Nsamp
      Read (15, *)
      Read (15, *) Dr
      Read (15, *)
      Read (15, *) Ndispl
      Read (15, *)
      Read (15, *) Npart, Temp, Rho

C     ---Initialise And Random Number Generator
      
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)
 
      If (Npart.Gt.Npmax) Then
         Write (6, *) ' Error: Number Of Particles Too Large'
         Stop
      End If

      Box  = (Dble(Npart)/Rho)**(1.D0/3.D0)
      Hbox = Box/2

C     ---Read Model Parameters

      Read (25, *)
      Read (25, *) Eps, Sig, Mass, Rc

C     ---Read/Generate Configuration

      If (Ibeg.Eq.0) Then

C        ---Generate Configuration Form Lattice

         Call Lattice
      Else
         Write (6, *) ' Read Conf From Disk '
         Read (11, *) Boxf
         Read (11, *) Npart
         Read (11, *) Dr
         Rhof = Dble(Npart)/Boxf**3
         If (Abs(Boxf-Box).Gt.1d-6) Then
            Write (6, 99007) Rho, Rhof
         End If
         Do I = 1, Npart
            Read (11, *) X(I), Y(I), Z(I)
            X(I) = X(I)*Box/Boxf
            Y(I) = Y(I)*Box/Boxf
            Z(I) = Z(I)*Box/Boxf
         End Do
         Rewind (11)
      End If

C     ---Write Input Data

      Write (6, 99001) Equil, Prod, Nsamp
      Write (6, 99002) Ndispl, Dr
      Write (6, 99003) Npart, Temp, Rho, Box
      Write (6, 99004) Eps, Sig, Mass

C     ---Calculate Parameters:

      Beta = 1/Temp
      
C     ---Calculate Cut-Off Radius Potential

      Rc    = Min(Rc, Hbox)
      Rc2   = Rc*Rc
      Eps4  = 4.0d0*Eps
      Eps48 = 48.0d0*Eps
      Sig2  = Sig*Sig
      
      Return

99001 Format ('  Number Of Equilibration Cycles             :', I10, /, 
     &        '  Number Of Production Cycles                :', I10, /, 
     &        '  Sample Frequency                           :', I10, /)
99002 Format ('  Number Of Att. To Displ. A Part. Per Cycle :', I10, /, 
     &        '  Maximum Displacement                       :', F10.3, 
     &        //)
99003 Format ('  Number Of Particles                        :', I10, /, 
     &        '  Temperature                                :', F10.3, 
     &        /, '  Density                                    :', 
     &        F10.3, /, '  Box Length                                 :'
     &        , F10.3, /)
99004 Format ('  Model Parameters: ', /, '     Epsilon: ', F5.3, /, 
     &        '     Sigma  : ', F5.3, /, '     Mass   : ', F5.3)
99007 Format (' Requested Density: ', F5.2, 
     &        ' Different From Density On Disk: ', F5.2, /, 
     &        ' Rescaling Of Coordinates!')
      End
