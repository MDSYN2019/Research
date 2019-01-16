      Integer Maxlat,Maxpart,Npart

      Parameter (Maxlat  = 65   )
      Parameter (Maxpart = 2000 )

      Integer Lattice(Maxlat,Maxlat),Ipart(Maxpart),
     &        Iipart(Maxpart),Mxx(Maxpart),Myy(Maxpart)

      Common /Sstrr/ Lattice,Ipart,Iipart,Npart,Mxx,Myy
