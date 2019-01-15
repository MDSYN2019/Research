      Subroutine Store(Iout, Dr, Vmax)

C     Writes Configuration To Disk

      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Integer Iout, I
      Double Precision Dr, Vmax
 
      Write (Iout, *) Box(1), Hbox(1), Box(2), Hbox(2)
      Write (Iout, *) Npart, Npbox(1), Npbox(2)
      Write (Iout, *) Dr, Vmax
      Do I = 1, Npart
         Write (Iout, *) X(I), Y(I), Z(I), Id(I)
      End Do
      Rewind (Iout)

      Return
      End
