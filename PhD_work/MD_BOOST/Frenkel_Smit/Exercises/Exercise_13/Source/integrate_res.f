      Subroutine Integrate_Res
      Implicit None
 
      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Integrate The Equations Of Motion Using An Explicit Nose-Hoover Chain     C
C     See The Work Of Martyna/Tuckerman Et Al.                                  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer I,J
      Double Precision U,F,Scale,Ukin,Aaa

Ccccccccccccccccccccccccccc
C     Update Thermostat   C
Ccccccccccccccccccccccccccc

      Scale    = 1.0d0
      Ukin     = Vpos**2
      Glogs(1) = (Ukin - Temp)*Iqmass(1)
               
      Do J = 1,5
 
         Vlogs(Nnoshover) = Vlogs(Nnoshover) + Glogs(Nnoshover)*Wdti4(J)
 
         Do I = 1,Nnoshover - 1
            Aaa                  = Dexp(-Wdti8(J)*Vlogs(Nnoshover-I+1))

            Vlogs(Nnoshover - I) = Vlogs(Nnoshover - I)*Aaa*Aaa + 
     &           Wdti4(J)*Glogs(Nnoshover - I)*Aaa
         Enddo
 
         Scale    = Scale*Dexp(-Wdti2(J)*Vlogs(1))
         Glogs(1) = (Scale*Scale*Ukin - Temp)*Iqmass(1)
 
         Do I = 1,Nnoshover
            Xlogs(I) = Xlogs(I) + Vlogs(I)*Wdti2(J)
         Enddo
 
         Do I = 1,Nnoshover - 1
            Aaa          = Dexp(-Wdti8(J)*Vlogs(I+1))
            Vlogs(I)     = Vlogs(I)*Aaa*Aaa + Wdti4(J)*Glogs(I)*Aaa
            Glogs(I + 1) = 
     &           (Qmass(I)*Vlogs(I)*Vlogs(I) - Temp)*Iqmass(I + 1)
         Enddo
            
         Vlogs(Nnoshover) = Vlogs(Nnoshover) + Glogs(Nnoshover)*Wdti4(J)
      Enddo
            
      Vpos = Vpos*Scale
        
Ccccccccccccccccccccccccccc
C     Normal Integration  C
Ccccccccccccccccccccccccccc

      Vpos = Vpos + 0.5d0*Tstep*Oldf
      Xpos = Xpos + Tstep*Vpos

      Call Force(Xpos,U,F)

      Oldf = F
      Vpos = Vpos + 0.5d0*Tstep*Oldf

Ccccccccccccccccccccccccccc
C     Update Thermostat   C
Ccccccccccccccccccccccccccc

      Scale    = 1.0d0
      Ukin     = Vpos**2
      Glogs(1) = (Ukin - Temp)*Iqmass(1)
               
      Do J = 1,5
 
         Vlogs(Nnoshover) = Vlogs(Nnoshover) + Glogs(Nnoshover)*Wdti4(J)
 
         Do I = 1,Nnoshover - 1
            Aaa                  = Dexp(-Wdti8(J)*Vlogs(Nnoshover-I+1))

            Vlogs(Nnoshover - I) = Vlogs(Nnoshover - I)*Aaa*Aaa + 
     &           Wdti4(J)*Glogs(Nnoshover - I)*Aaa
         Enddo
 
         Scale    = Scale*Dexp(-Wdti2(J)*Vlogs(1))
         Glogs(1) = (Scale*Scale*Ukin - Temp)*Iqmass(1)
 
         Do I = 1,Nnoshover
            Xlogs(I) = Xlogs(I) + Vlogs(I)*Wdti2(J)
         Enddo
 
         Do I = 1,Nnoshover - 1
            Aaa          = Dexp(-Wdti8(J)*Vlogs(I+1))
            Vlogs(I)     = Vlogs(I)*Aaa*Aaa + Wdti4(J)*Glogs(I)*Aaa
            Glogs(I + 1) = 
     &           (Qmass(I)*Vlogs(I)*Vlogs(I) - Temp)*Iqmass(I + 1)
         Enddo
            
         Vlogs(Nnoshover) = Vlogs(Nnoshover) + Glogs(Nnoshover)*Wdti4(J)
      Enddo
            
      Vpos = Vpos*Scale
      Cons = U + 0.5d0*Vpos*Vpos

Cccccccccccccccccccccccccccccccccccccccc
C     Calculate Conservative Quantity  C
Cccccccccccccccccccccccccccccccccccccccc

      Do J = 1,Nnoshover
         Cons = Cons + 0.5d0*Qmass(J)*Vlogs(J)*Vlogs(J) + Temp*Xlogs(J)
      Enddo

      Return
      End
