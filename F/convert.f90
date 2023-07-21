      program convert
      implicit none
      real*8 r,deltm,dmt,dmi,dmc,dmr,mass,time,delt,mdot,fact,ddist
      real*8 fpi,vol
      integer*4 ijunk,ilines,m,n
      character*20 junk
      dimension r(200),deltm(200),dmt(200),dmi(200),dmc(200),dmr(200)
      dimension ddist(200),vol(200)
      open(unit=1,file='blok.txt')
      open(unit=2,file='tes.dat')
!      read(1,1)junk
!1     format(a20)
!2     format(i5,i4,3e12.5)
!3     format(i5,e16.8,e12.5,4e13.6)
!7     format(i5,e13.6)
!8     format(i5)
!      read(1,1)junk
!        read(1,2)ijunk,ilines,time,delt,mdot
!        read(1,1)junk
!        read(1,1)junk
        ilines=159
        mdot=8.576E-06
        do 4 m=1,ilines
        read(1,*)n,r(m),deltm(m),dmt(m),dmi(m),dmc(m),dmr(m)
4     continue
      if(n.ne.ilines) print *, 'there is something wrong'
      fpi=4.*3.14159
        mass=0.
      mass=dmt(1)
      ddist(1)=.0001*r(1)
      do 5 m=2,ilines
        mass=mass+dmt(m)
        ddist(m)=r(m-1)-r(m)
5     continue
      fact=mdot/3.15d7/mass*5.98d27
      write(2,*)'--- BEGIN HEADER ---'
      write(2,*)'Mass input into layers of the atmosphere by planetesimal breakup.'
      write(2,*)'Columns:'
      write(2,*)'[Height (cm)]  [mass input (g/cm^3/s)]'
      write(2,*)'--- END HEADER ---'
      write(2,*)''
      write(2,*)'--- BEGIN DATA ---'
      do 6 m=1,ilines
      vol(m)=fpi*r(m)*r(m)*ddist(m)
        dmt(m)=dmt(m)*fact/vol(m)
        write(2,*)r(m),dmt(m)
6     continue
      close(2)
      close(1)
      stop
        end