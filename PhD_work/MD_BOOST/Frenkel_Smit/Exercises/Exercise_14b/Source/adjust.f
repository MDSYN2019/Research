      Subroutine Adjust(Attemp, Nacc, Dr, Attv, Accv, Vmax, Succ)

C     Sets Maximum Displacement And Maximum Volume Change
C     Such That 50% Of The Move Will Be Accepted

      Implicit None

      Include 'system.inc'

      Integer Attemp, Nacc, Attempp, Naccp, Accvp, Attv, Accv, Attvp
      Double Precision Dro, Frac, Vmaxo, Dr, Vmax, Succ
      Save Naccp, Attempp, Attvp, Accvp

      Data Attempp/0/
      Data Attvp/0/

C     ---Displacement:

      If (Attemp.Eq.0.Or.Attempp.Ge.Attemp) Then
         Naccp = Nacc
         Attempp = Attemp
      Else
         Frac = Dble(Nacc-Naccp)/Dble(Attemp-Attempp)
         Dro = Dr
         Dr = Dr*Abs(Frac/(Succ/100.D0))

C        ---Limit The Change:

         If (Dr/Dro.Gt.1.5d0) Dr = Dro*1.5d0
         If (Dr/Dro.Lt.0.5d0) Dr = Dro*0.5d0
         If (Dr.Gt.Hbox(1)/2.D0) Dr = Hbox(1)/2.D0
         Write (6, 99001) Dr, Dro, Frac, Attemp - Attempp, Nacc - Naccp

C        ---Store Nacc And Attemp For Next Use

         Naccp = Nacc
         Attempp = Attemp
      End If

C     ---Volume:

      If (Attv.Eq.0.Or.Attvp.Ge.Attv) Then
         Accvp = Accv
         Attvp = Attv
      Else
         Frac = Dble(Accv-Accvp)/Dble(Attv-Attvp)
         Vmaxo = Vmax
         Vmax = Vmax*Abs(Frac/(Succ/100.D0))

C        ---Limit The Change:

         If (Vmax/Vmaxo.Gt.1.5d0) Vmax = Vmaxo*1.5d0
         If (Vmax/Vmaxo.Lt.0.5d0) Vmax = Vmaxo*0.5d0
         Write (6, 99002) Vmax, Vmaxo, Frac, Attv - Attvp, Accv - Accvp

C        ---Store Nacc And Attemp For Next Use

         Accvp = Accv
         Attvp = Attv
      End If
      Return
99001 Format (' Max. Displ. Set To      : ', F6.3, ' (Old : ', F6.3, 
     &        ')', /, ' Frac. Acc.: ', F4.2, ' Attempts: ', I7, 
     &        ' Succes: ', I7)
99002 Format (' Max. Vol. Chan. Set To: ', F6.3, ' (Old : ', F6.3, ')', 
     &        /, ' Frac. Acc.: ', F4.2, ' Attempts: ', I7, ' Succes: ', 
     &        I7)
      End
