program empty

use globals
use jack

implicit none
real(8) :: dt,startMass,curMass,from,until

  call oSetup('cfg2321.dat')
  call debugMass()
  startMass=sum(debug_massInAtm)
  curMass=startMass
  from=0; dt=5e6; until=dt
  
  do while (curMass/startMass>0.1)
    call jack4(from,until)
    call debugMass()
    curMass=sum(debug_massInAtm)
    from=until
    until=until+dt
    !write(*,*)curMass/startMass
  end do
  
  open(11,file='emptytest.txt',position='append')
  write(11,'(i0,es12.2,es12.2)')params_nBins,until-dt,&
                                dust_sizeBin(params_nBins)
  close(11)


stop
end program empty
