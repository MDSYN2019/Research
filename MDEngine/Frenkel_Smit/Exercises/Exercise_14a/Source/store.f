      Subroutine Store(Iout, Dr)
C
C     Writes Configuration To Disk
C
C     Iout (Input) File Number
C     Dr   (Input) Maximum Displacement
C
C
      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Integer Iout, I
      Double Precision Dr
 
      Write (Iout, *) Box, Hbox
      Write (Iout, *) Npart
      Write (Iout, *) Dr

      Do I = 1, Npart
         Write (Iout, *) X(I), Y(I), Z(I)
      End Do

      Rewind (Iout)

      Return
      End
