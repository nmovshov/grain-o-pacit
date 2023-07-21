program opac

!----------------------------------------------------------------------
! This program estimates the contribution of dust grains to the opacity
! of a model atmosphere. To do this, the program tracks the evolution
! of dust particles in the atmosphere, taking into account growth of
! particles by collisions, and the sedimentation of particles through
! the atmosphere.
!
! The code is minimally commented for readability. I will write a more
! detailed help file later.
!
! The structure of the program is as follows:
!
! opac.f90 (this file) -> this is the 'driver', or 'main'.
!
!  |- globals.f90 -> module, global variables and routines.
!  |
!  |- jack.f90 -> module, grain evolution calculations.
!  |
!  |- radt.f90 -> module, radiative transfer calculations.
!
! To build, first compile globals to an object file, then compile jack
! and radt to objects, then compile opac and link.
!
! In addition to these source files, the following input files should
! be in the same directory as the executable:
!  cfg.dat -> configuration file, model parameters, file names, control
!             options
!  atm.dat -> model atmosphere, height, density, temperature, cgs
!  ind.dat -> table of wavelength-dependent refractive indices
!  tes.dat -> rate of solids input into the atmosphere at different
!             heights, from planetesimal breakup, *optional*
! dist.dat -> initial distribution of grains of all sizes in every
!             layer, *optional*
!
! Naor Movshovitz
! naormovs@post.tau.ac.il
!
!----------------------------------------------------------------------

use globals
use jack
use radt
                   
implicit none
character(len=20) :: configFile='' 
character(len=50) :: out1='',out2=''
character(len=10) :: dateTime='',dateDate=''
integer :: j,k,N
real(8) :: ct1,ct2,m1,m2,from,until


!**********************************************************
! The following "template" shows how to run the program and produce
! output files.
!**********************************************************

! First, call the setup routine.
  configFile='cfg.dat'
  call oSetup(configFile)

! Optional: call cpu_time to time the run.
  call cpu_time(ct1)

do j=2,10               !
  if(tTimes(j)==0) exit ! there can be a number of target times
end do                  !
N=j-2                   !
do j=1,N ! Run to N target times

! Get current date and time (used in file naming)
  call date_and_time(dateDate,dateTime)

! Optional: call debugMass and remember the total mass of grains
  call debugMass()
  m1=sum(debug_massInAtm)+debug_massInCore
    
! Next, call the grain evolution routine. Note that this is the time
! consuming part of the program. Progress messages will be printed on
! the screen from within the jack4 subroutine. To edit these messages
! go to 'jack.f90' and edit lines 165-169.
  from=tTimes(j); until=tTimes(j+1) ! (seconds)
  call jack4(from,until)
  
! Optional: call debugMass again and check for conservation of total
!           mass. A significant mass imbalance may indicate a problem.
  call debugMass()
  m2=sum(debug_massInAtm)+debug_massInCore
  write(*,*)'Total start mass was ',m1,' g.'
  write(*,*)'Total end mass is ',m2,' g.'
  write(*,*)'Total mass difference is ',m2-m1,' g.'
  k=0; if (cfg_srcFlg) k=1
  write(*,*)'Total mass input was ',&
             sum(dust_source*atm_dV)*(until-from)*dble(k),' g.'
             
! Next, call the radiative transfer routine.
  call radt1()
  
! Now create two output files, one with the distribution of dust grains
! and the other with the opacity parameters of the atmosphere. A unique
! file name based on the date and time is generated to avoid running
! over files from previous runs. To write into the same files each time
! edit the following two lines.
  out1='opacOutDist'//trim(dateDate)//trim(dateTime(1:6))//'.dat'
  out2='opacOutOptic'//trim(dateDate)//trim(dateTime(1:6))//'.dat'

  open(11,file=out1,action='write')
  write(11,*)'--- BEGIN HEADER ---'
  write(11,*)'This is output from ''opac.f90'' run on '//&
             dateDate(7:8)//'/'//dateDate(5:6)//'/'//dateDate(1:4)//&
             ', at '//dateTime(1:2)//':'//dateTime(3:4)//':'//&
             dateTime(5:6)//'.'
  write(11,*)'Evolution time: ',(from),'--',(until),' seconds.'
  write(11,*)'Atmosphere: ',trim(cfg_atmFile)
  write(11,*)'Used parameters:'
  write(11,*)'monomerSize = ',params_monomerSize
  write(11,*)'nBins = ',params_nBins
  write(11,*)'spacing parameter = ',params_binSpacingParameter
  write(11,*)'sticking coefficient = ',params_stickCoef
  write(11,*)'dust-to-gas ratio = ',params_dust2gasRatio
  write(11,*)'srcSize = ',cfg_srcSize
  write(11,*)'coag =',cfg_cogFlg,' sed =',cfg_flxFlg,' src =',cfg_srcFlg
  write(11,*)
  write(11,*)'In this file, the number in the (j,k) place is the '//&
             'number density of k-th bin grains in the j-th layer, '//&
             'in cm^-3.'
  write(11,*)'--- END HEADER ---'
  write(11,*)
  write(11,*)'--- BEGIN DATA ---' ! it helps to have this line when we
                                  ! read from the file to some program.
  do k=1,atm_nZones
    write(11,*)dust_nDensity(k,1:params_nBins)
  end do
  close(11)

  open(11,file=out2,action='write')
  write(11,*)'--- BEGIN HEADER ---'
  write(11,*)'This is output from ''opac.f90'' run on '//&
             dateDate(7:8)//'/'//dateDate(5:6)//'/'//dateDate(1:4)//&
             ', at '//dateTime(1:2)//':'//dateTime(3:4)//':'//&
             dateTime(5:6)//'.'
  write(11,*)'Evolution time: ',(from),'--',(until),' seconds.'
  write(11,*)'Atmosphere: ',trim(cfg_atmFile)
  write(11,*)'Used parameters:'
  write(11,*)'monomerSize = ',params_monomerSize
  write(11,*)'nBins = ',params_nBins
  write(11,*)'spacing parameter = ',params_binSpacingParameter
  write(11,*)'sticking coefficient = ',params_stickCoef
  write(11,*)'dust-to-gas ratio = ',params_dust2gasRatio
  write(11,*)'srcSize = ',cfg_srcSize
  write(11,*)'coag =',cfg_cogFlg,' sed =',cfg_flxFlg,' src =',cfg_srcFlg
  write(11,*)
  write(11,*)'Columns in this file:'
  write(11,*)'[height (cm)]  [opacity (cm^2/g)]  [optical depth]'
  write(11,*)'--- END HEADER ---'
  write(11,*)
  write(11,*)'--- BEGIN DATA ---' ! it helps to have this line when we
                                  ! read from the file to some program.
  do k=1,atm_nZones
    write(11,*)atm_Z(k),atm_opacity(k),atm_opticalDepth(k)
  end do
  close(11)

end do ! on target times...

! Optional: call cpu_time to time the run.
  call cpu_time(ct2)
  write(*,*)'Elapsed time is ',ct2-ct1,' seconds.'

! That's it.
  stop
end program opac