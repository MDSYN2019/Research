      Function Ran_Vel()
      Implicit None

C     Generates A Random Velocity According To A Boltzmann
C     Distribution

      Include 'system.inc'
 
      Double Precision Ran_Vel,Ran_Gauss
 
      Ran_Vel = Dsqrt(Temp)*Ran_Gauss()
 
      Return
      End
