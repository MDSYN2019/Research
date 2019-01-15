      Subroutine Sample(I, En, Vir, Press)
      Implicit None
C
C      Write Quantities (Pressure And Energy) To File
C
C
C      Ener (Input) : Total Energy
C      Vir  (Input) : Total Virial
C
C
      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Integer I
      Double Precision En, Enp, Vir, Press,Vol
 
      If (Npart.Ne.0) Then
         Enp   = En/Dble(Npart)
         Vol   = Box**3
         Press = (Dble(Npart)/Vol)/Beta + Vir/(3.0d0*Vol)
      Else
         Enp   = 0.0d0
         Press = 0.0d0
      End If

      Write (66, *) I, Enp
      Write (67, *) I, Press

      Return
      End
