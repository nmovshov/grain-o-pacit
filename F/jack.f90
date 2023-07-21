module jack
!----------------------------------------------------------------------
! This module contains the subroutine that tracks the evolution of dust
! particles in the model atmosphere. This subroutine is a rewrite of
! Morris's program jack.for with both expansions and omissions.
!----------------------------------------------------------------------

use globals
implicit none

CONTAINS

subroutine jack4(sTime,eTime)
! -------------  Evolution of dust grains in the atmosphere -----------
implicit none
real(8), intent(in) :: sTime, eTime ! start/end time
real(8) :: nd(200,100),dndt(200,100),dndtFl(200,100),dndtCg(200,100),&
           dndtSr(200,100),flxCdt(100)

real(8) :: time,dt,fluxOu,fluxIn,volume,x,y,pi,tol
integer :: j,k,m,n,nBins,aBins,nZones,iCount

  nd=dust_nDensity
  dndt=0
  dndtFl=0
  dndtCg=0
  dndtSr=0
  flxCdt=0
  nBins=params_nBins
  aBins=dust_aBins
  nZones=atm_nZones
  pi=3.141592
  
  time=sTime
  iCount=0
  do while (time.lt.eTime)
      
!**********************************************************
! Find the rate of change of the number density.
!**********************************************************
  if (cfg_flxFlg) then
! **** Calculate sedimentation.
    dndtFl=0
    
! Uppermost zone - only flux out (flux in included in source term).
    do k=1,nBins
      fluxOu=nd(1,k)*dust_vSed(1,k)*4*pi*atm_R(2)**2
      volume=atm_dV(1)
      dndtFl(1,k)=(-fluxOu)/volume
    end do
! The other zones - flux in and out.
    do j=2,nZones
      do k=1,nBins
        fluxIn=nd(j-1,k)*dust_vSed(j-1,k)*4*pi*atm_R(j)**2
        fluxOu=nd(j,k)*dust_vSed(j,k)*4*pi*atm_R(j+1)**2
        volume=atm_dV(j)
        dndtFl(j,k)=(fluxIn-fluxOu)/volume
      end do
    end do
! The "fallout zone" - only flux in.*
! * The fallout zone collects all the mass that has sedimented out of
!   the last zone into a "flux capacitor", to allow a check on the
!   conservation of total mass.
    do k=1,nBins
      fluxIn=nd(nZones,k)*dust_vSed(nZones,k)*4*pi*atm_R(nZones+1)**2
      volume=1
      flxCdt(k)=(fluxIn)/volume
    end do
            
  else
! **** Do not include sedimentation.
    continue
  end if
      
  if (cfg_cogFlg) then
! **** Calculate coagulation.

    do j=1,nZones
      do k=1,nBins
        dndtCg(j,k)=0
        do m=1,k-1
          do n=1,m-1
            dndtCg(j,k)=dndtCg(j,k)+&
                 dust_kernel(j,m,n)*nd(j,m)*nd(j,n)*dust_redist(m,n,k)
          end do
          dndtCg(j,k)=dndtCg(j,k)+&
               dust_kernel(j,m,m)*0.5*nd(j,m)**2*dust_redist(m,m,k)
        end do
        do m=1,nBins
          dndtCg(j,k)=dndtCg(j,k)-&
            dust_kernel(j,k,m)*nd(j,k)*nd(j,m)*(1-dust_redist(k,m,k))
        end do
      end do
    end do
! Return "spilled out" mass to active bins.
    do j=1,nZones
      do k=aBins+1,nBins
        dndtCg(j,aBins)=dndtCg(j,aBins)+&
                        dndtCg(j,k)*dust_massBin(k)/dust_massBin(aBins)
        dndtCg(j,k)=0
      end do
    end do

  else
! **** Do not include coagulation.
    continue
  end if
      
  if (cfg_srcFlg) then
! **** Add a source term.*
  if (cfg_srcSize.eq.'small') then
    dndtSr(:,1)=dust_source/dust_massBin(1)
  else if (cfg_srcSize.eq.'equal') then
    do k=1,nZones
      dndtSr(k,1:aBins)=(dust_source(k)/aBins)/dust_massBin(1:aBins)
    end do
  else
    write(*,*)'Error: (jack.f90)Unknown source size option'; stop
  end if
  
  else
! **** Do not include a source
    continue
  end if

! **** Sum up the rate of change from all the processes.
  do j=1,nZones
    do k=1,nBins
      dndt(j,k)=dndtFl(j,k)+dndtCg(j,k)+dndtSr(j,k)
    end do
  end do      
      
!**********************************************************
! Determine an appropriate time step.
!**********************************************************
  x=0
  do j=1,nZones
    do k=1,nBins
      tol=atm_dV(j)**(-1)
      if (nd(j,k).gt.tol) then
        y=abs(dndt(j,k)/nd(j,k))
        if (y.gt.x) x=y
      end if
    end do
  end do
  dt=0.01/x
  if (dt.gt.params_maxTimeStep) dt=params_maxTimeStep
  if ((time+dt).gt.eTime) dt=eTime-time
      
!**********************************************************
! Update the nd array.
!**********************************************************
  do j=1,nZones
    do k=1,nBins
      nd(j,k)=nd(j,k)+dndt(j,k)*dt
      if (nd(j,k).lt.0) nd(j,k)=0
    end do
  end do
  dust_nDensity=nd

!**********************************************************
! Update the "flux capacitor".
!**********************************************************
  if (cfg_flxFlg) then
    do k=1,nBins
      debug_ndCore(k)=debug_ndCore(k)+flxCdt(k)*dt
    end do
  end if
            
!**********************************************************
! Print progress and debug messages (optional).
!**********************************************************
  if (mod(iCount,1000).eq.0 .and. cfg_msgFlg) then
    write(*,*) 'Time is      ',time+dt,' seconds.'
    write(*,*) 'Time step is ',dt,' seconds.'
    write(*,*)
  end if
      
!**********************************************************
! Advance the time variable and iteration counter, and go back
! to the start line, or not.
!**********************************************************
  time=time+dt
  iCount=iCount+1
    
  end do ! while

return
end subroutine jack4

end module jack