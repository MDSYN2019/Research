      Subroutine Writepdb
      Implicit None

      Include 'system.inc'

      Integer I,Countmodel,Countatom

      Data Countmodel/ 0/
      Data Countatom/ 0/

      Save Countmodel,Countatom

      Countmodel = Countmodel + 1

      Write(22,'(A,I9)') 'MODEL',Countmodel

      Do I=1,Nuall

         Countatom = Countatom + 1

         Write(22,'(A,I7,A,I12,4x,3f8.3)') 'ATOM',Countatom,'  O',
     &        Countatom,Xpos(I),Ypos(I),Zpos(I)
      Enddo

      Write(22,'(A)') 'ENDMDL'

      Return
      End
