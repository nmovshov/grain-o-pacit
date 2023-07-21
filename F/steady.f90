program steady
!--------------------------------------------------------------------
! Run the grain sedimentation routine to steady state, defined by
! M_dot_in=M_dot_out.
!--------------------------------------------------------------------

use globals
use jack
use radt
                   
implicit none
character(len=20) :: configFile='' 
character(len=50) :: out1='',out2=''
character(len=10) :: dateTime='',dateDate=''
integer :: j=0,k,safety=10
real(8) :: dt,curMass,nextMass,from,until,ct1,ct2,mark=0,Mp

  call cpu_time(ct1)
    
  configFile='cfg.dat'
  call oSetup(configFile)
  call debugMass()
  curMass=sum(debug_massInAtm)
  from=0; dt=1e7; until=dt;
  
  !Mp=sum(dust_source*atm_dV)-params_accretionRate*4*3.14*atm_R(1)**2
  !dust_source(1)=1*Mp/atm_dV(1)+params_accretionRate/atm_dR(1)
  !dust_source(2:)=0

  open(11,file='mdot.txt',action='write',status='replace')
  write(11,*)'--- BEGIN HEADER ---'
  write(11,*)'This is output from steady.f90, used to test the model for steady state'
  write(11,*)'based on the condition M_dot_in=M_dot_out.'
  write(11,*)'Columns:'
  write(11,*)'[time (s)]  [delta_M_solids/delta_t (g/s)]'
  write(11,*)'--- END HEADER ---'
  write(11,*)
  write(11,*)'--- BEGIN DATA ---'
  
  do
    call jack4(from,until)
    call debugMass()
    nextMass=sum(debug_massInAtm)
    if (equals(curMass,nextMass)) then 
      j=j+1
    else
      j=0
    end if
    if (j==safety) then
      mark=until; exit
    end if
    if (until>=1e11) then
      write(*,*)'Warning: steady state not reached.'; exit
    end if
!    write(*,*)(nextMass-curMass)/dt
    write(11,*)until,(nextMass-curMass)/dt
    curMass=nextMass
    from=until; until=until+dt
  end do
  
  close(11)
  
  if (j==safety) then
    write(*,'(a,es10.3,a)')'Steady state mark at ',mark,' s'
  end if
  
  open(11,file='steady.txt',position='append')
  write(11,'(i0,es12.2,es12.2)')params_nBins,mark,&
                                dust_sizeBin(params_nBins)
  close(11)

  
  call radt1()
  
  out1='opacOutDistSteady.dat'
  out2='opacOutOpticSteady.dat'
  
!  open(11,file=out1,action='write',status='replace')
!  write(11,*)'--- BEGIN HEADER ---'
!  write(11,*)'This is output from ''opac.f90'' run on '//&
!             dateDate(7:8)//'/'//dateDate(5:6)//'/'//dateDate(1:4)//&
!             ', at '//dateTime(1:2)//':'//dateTime(3:4)//':'//&
!             dateTime(5:6)//'.'
!  write(11,*)'Evolution time: ',(from),'--',(until),' seconds.'
!  write(11,*)'Atmosphere: ',trim(cfg_atmFile)
!  write(11,*)'Used parameters:'
!  write(11,*)'monomerSize = ',params_monomerSize
!  write(11,*)'nBins = ',params_nBins
!  write(11,*)'spacing parameter = ',params_binSpacingParameter
!  write(11,*)'srcSize = ',cfg_srcSize
!  write(11,*)'coag =',cfg_cogFlg,' sed =',cfg_flxFlg,' src =',cfg_srcFlg
!  write(11,*)
!  write(11,*)'In this file, the number in the (j,k) place is the '//&
!             'number density of k-th bin grains in the j-th layer, '//&
!             'in cm^-3.'
!  write(11,*)'--- END HEADER ---'
!  write(11,*)
!  write(11,*)'--- BEGIN DATA ---' ! it helps to have this line when we
!                                  ! read from the file to some program.
!  do k=1,atm_nZones
!    write(11,*)dust_nDensity(k,1:params_nBins)
!  end do
!  close(11)
!
!  open(11,file=out2,action='write',status='replace')
!  write(11,*)'--- BEGIN HEADER ---'
!  write(11,*)'This is output from ''opac.f90'' run on '//&
!             dateDate(7:8)//'/'//dateDate(5:6)//'/'//dateDate(1:4)//&
!             ', at '//dateTime(1:2)//':'//dateTime(3:4)//':'//&
!             dateTime(5:6)//'.'
!  write(11,*)'Evolution time: ',(from),'--',(until),' seconds.'
!  write(11,*)'Atmosphere: ',trim(cfg_atmFile)
!  write(11,*)'Used parameters:'
!  write(11,*)'monomerSize = ',params_monomerSize
!  write(11,*)'nBins = ',params_nBins
!  write(11,*)'spacing parameter = ',params_binSpacingParameter
!  write(11,*)'srcSize = ',cfg_srcSize
!  write(11,*)'coag =',cfg_cogFlg,' sed =',cfg_flxFlg,' src =',cfg_srcFlg
!  write(11,*)
!  write(11,*)'Columns in this file:'
!  write(11,*)'[height (cm)]  [opacity (cm^2/g)]  [optical depth]'
!  write(11,*)'--- END HEADER ---'
!  write(11,*)
!  write(11,*)'--- BEGIN DATA ---' ! it helps to have this line when we
!                                  ! read from the file to some program.
!  do k=1,atm_nZones
!    write(11,*)atm_Z(k),atm_opacity(k),atm_opticalDepth(k)
!  end do
!  close(11)
  
  call cpu_time(ct2)
  write(*,*)'Elapsed time is ',ct2-ct1
  
  stop
  
CONTAINS

function equals (a, b)
! reasonably safe compare of two default REALs
implicit none
logical :: equals
real(8), intent(in) :: a, b
real(8) :: eps
eps = abs(a) * epsilon(a) ! scale epsilon
if (eps == 0) then
  eps = tiny (a) ! if eps underflowed to 0 use a very small
                 ! positive value for epsilon
end if
if (abs(a-b) > eps) then
  equals = .false. ! not equal if difference>eps
  return
else
  equals = .true.  ! equal otherwise
  return
endif
end function equals

end program steady
