      Program Polymer
      Implicit None

      Include 'system.inc'

      Integer Sstmm
      Double Precision M1

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)
 
      Read(21,*) Nstep,Ninit,Nuall,Nchoi
      Read(21,*) Lstatic,Lcbmc
      Read(21,*) Beta,Rcut,Prepot,Kb,Thetan

      Write(6,*) 'Cbmc Program Of A Single Chain'
      Write(6,*)
      
      If(Lstatic) Then
         Write(6,*) 'Static Scheme'
      Else
         Write(6,*) 'Dynamic Scheme'
      Endif

      If(Lcbmc) Then
         Write(6,*) 'Cbmc Is Used'
      Else
         Write(6,*) 'Random Chains Are Grown'
      Endif

      Write(6,*)
      Write(6,*) 'Number Of Cycles  : ',Nstep
      Write(6,*) 'Number Of Init.   : ',Ninit
      Write(6,*) 'Temperature       : ',Beta
      Write(6,*) 'Chain Length      : ',Nuall
      Write(6,*) 'Number Of Trials  : ',Nchoi
      Write(6,*)
      Write(6,*) 'Rcut              : ',Rcut
      Write(6,*) 'A                 : ',Prepot
      Write(6,*) 'Kb                : ',Kb
      Write(6,*) 'Thetan [Radials]  : ',Thetan
      Write(6,*)

      Beta   = 1.0d0/Beta
      Prepot = Prepot/(Rcut*Rcut)

      If(Nuall.Gt.Maxchain) Then
         Write(6,*) 'Chain Is Too L0ng !!!'
         Stop
      Endif

      If(Nchoi.Gt.Maxtrial) Then
         Write(6,*) 'Too Many Trial Positions !!!'
         Stop
      Endif

      Call Mcloop
 
      Stop
      End
