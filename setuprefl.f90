subroutine setuprefl(lunin,mype,bwork,awork,nele,nobs,is,conv_diagsave)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    setuprefl     compute rhs of oi for radar radial winds
!
! abstract: For radar radial wind observations, this routine
!              a) reads obs assigned to given mpi task (geographic region),
!              b) simulates obs from guess,
!              c) apply some quality control to obs,
!              d) load weight and innovation arrays used in minimization
!              e) collects statistics for runtime diagnostic output
!              f) writes additional diagnostic information to output file
!
!   input argument list:
!     lunin    - unit from which to read observations
!     mype     - mpi task id
!     nele     - number of data elements per observation
!     nobs     - number of observations
!
!   output argument list:
!     bwork    - array containing information about obs-ges statistics
!     awork    - array containing information for data counts and gross checks
!
!$$$
  use mpeu_util, only: die,perr
  use kinds, only: r_kind,r_single,r_double,i_kind

  use m_obsdiags, only: reflhead
  use obsmod, only: rmiss_single,i_refl_ob_type,obsdiags,lobsdiag_forenkf,&
                    lobsdiagsave,nobskeep,lobsdiag_allocated,time_offset
  use obsmod, only: netcdf_diag, binary_diag, dirname,ianldate
  use nc_diag_write_mod, only: nc_diag_init, nc_diag_header, nc_diag_metadata, &
       nc_diag_write, nc_diag_data2d
  use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_get_dim, nc_diag_read_close
  use m_obsNode, only: obsNode
  use m_reflNode, only: reflNode
  use m_obsLList, only: obsLList_appendNode
  use obsmod, only: obs_diag,luse_obsdiag
  use gsi_4dvar, only: nobs_bins,hr_obsbin
  use oneobmod, only: magoberr,maginnov,oneobtest
  use qcmod, only: npres_print,ptop,pbot,tdrerr_inflate
  use qcmod, only: vqc
  use guess_grids, only: hrdifsig,geop_hgtl,nfldsig,&
       ges_lnprsl,sfcmod_gfs,sfcmod_mm5,comp_fact10
  use gridmod, only: nsig,get_ijk
  use constants, only: flattening,semi_major_axis,grav_ratio,zero,grav,wgtlim,&
       half,one,two,grav_equator,eccentricity,somigliana,rad2deg,deg2rad
  use constants, only: tiny_r_kind,cg_term,huge_single,r2000,three,one
  use jfunc, only: jiter,last,miter,jiterstart
  use convinfo, only: nconvtype,cermin,cermax,cgross,cvar_b,cvar_pg,ictype
  use convinfo, only: icsubtype
  use m_dtime, only: dtime_setup, dtime_check, dtime_show
  use gsi_bundlemod, only : gsi_bundlegetpointer
  use gsi_metguess_mod, only : gsi_metguess_get,gsi_metguess_bundle
  use sparsearr, only: sparr2, new, size, writearray, fullarray
  use state_vectors, only: nsdim
  implicit none

! Declare passed variables
  logical                                          ,intent(in   ) :: conv_diagsave
  integer(i_kind)                                  ,intent(in   ) :: lunin,mype,nele,nobs
  real(r_kind),dimension(100+7*nsig)               ,intent(inout) :: awork
  real(r_kind),dimension(npres_print,nconvtype,5,3),intent(inout) :: bwork
  integer(i_kind)                                  ,intent(in   ) :: is ! ndat index


! Declare local parameters
  real(r_kind),parameter:: r0_001 = 0.001_r_kind
  real(r_kind),parameter:: r8     = 8.0_r_kind
  real(r_kind),parameter:: ten    = 10.0_r_kind
  real(r_kind),parameter:: r200   = 200.0_r_kind

! Declare external calls for code analysis
  external:: tintrp2a1,tintrp2a11
  external:: tintrp31
  external:: grdcrd1
  external:: stop2

! Declare local variables
  real(r_kind) rlow,rhgh,rsig
  real(r_kind) dz,factelv,factdif
  real(r_kind) dlnp,pobl,zob
  real(r_kind) sin2,termg,termr,termrg
  real(r_kind) psges,zsges
  real(r_kind),dimension(nsig):: zges,hges
  real(r_kind) prsltmp(nsig)
  real(r_kind) sfcchk  
  real(r_kind) residual,obserrlm,obserror,ratio,scale,val2
  real(r_kind) ress,ressw
  real(r_kind) val,valqc,rwgt
  real(r_kind) cg_refl,wgross,wnotgross,wgt,arg,exp_arg,term,rat_err2
  real(r_kind) dlat,dlon,dtime,dpres,ddiff,error,slat
  real(r_kind) ratio_errors,qcgross
  real(r_kind) reflges,skint,sfcr
  real(r_kind) presrefl
  real(r_kind) errinv_input,errinv_adjst,errinv_final
  real(r_kind) err_input,err_adjst,err_final
  real(r_kind),dimension(nele,nobs):: data
  real(r_single),allocatable,dimension(:,:)::rdiagbuf

  integer(i_kind) i,nchar,nreal,k,j,k1,ii
  integer(i_kind) mm1,jj,k2,isli
  integer(i_kind) jsig,ikxx,nn,ibin,ioff,ioff0
  integer(i_kind) ier,ilat,ilon,ihgt,ireflob,ikx,itime,iuse,iqc,izz
  integer(i_kind):: id,ilone,ilate
  integer(i_kind):: ier2,idomsfc,isfcr,iskint,iff10
  integer(i_kind):: itype
  real(r_kind)   :: reflob,var_jb
  
  real(r_double) rstation_id
  character(8) station_id
  equivalence(rstation_id,station_id)

  character(8),allocatable,dimension(:):: cdiagbuf

  logical,dimension(nobs):: luse,muse
  integer(i_kind),dimension(nobs):: ioid ! initial (pre-distribution) obs ID
  logical proceed
  logical include_refl

  integer(i_kind) istat

  type(sparr2) :: dhx_dx
  real(r_single), dimension(nsdim) :: dhx_dx_array
  integer(i_kind) :: nnz, nind

  logical:: in_curbin, in_anybin, save_jacobian
  integer(i_kind),dimension(nobs_bins) :: n_alloc
  integer(i_kind),dimension(nobs_bins) :: m_alloc
  class(obsNode),pointer:: my_node
  type(reflNode),pointer:: my_head
  type(obs_diag),pointer:: my_diag
  character(len=*),parameter:: myname='setuprefl'

  real(r_kind),allocatable,dimension(:,:,:  ) :: ges_ps
  real(r_kind),allocatable,dimension(:,:,:  ) :: ges_z
  real(r_kind),allocatable,dimension(:,:,:,:) :: ges_refl

  logical :: yn_1stprint1=.true., yn_1stprint2=.true., yn_1stprint3=.true. ! wangjia: to make the subroutine only print once when necessary

! lobsdiag_forenkf - if true, save linearized H operator (jacobian) in diagnostic file on 1st outer iteration.  The Jacobian can then be used by the EnKF to compute ensemble perturbations in observation space.
  save_jacobian = conv_diagsave .and. jiter==jiterstart .and. lobsdiag_forenkf

! Check to see if required guess fields are available
  call check_vars_(proceed,include_refl)
  if((.not.proceed) .or. (.not.include_refl)) return  ! not all vars available, simply return

! If require guess vars available, extract from bundle ...
  call init_vars_

  n_alloc(:)=0
  m_alloc(:)=0
!*******************************************************************************
! Read and reformat observations in work arrays.
  read(lunin)data,luse,ioid

!    index information for data array (see reading routine)
  ier=1       ! index of obs error
  ilon=2      ! index of grid relative obs location (x)
  ilat=3      ! index of grid relative obs location (y)
  ihgt=4      ! index of obs elevation, [m]
  ireflob=5   ! index of reflectivity observation, [dBZ]
  id=6        ! index of station id
  itime=7     ! index of observation time in data array
  ikxx=8      ! index of obs type in data array
  iqc=9       ! index of quality mark
  ier2=10     ! index of original obs error
  iuse=11     ! index of use parameter
  idomsfc=12  ! index of dominant surface type
  iskint=13   ! index of surface skin temperature
  iff10=14    ! index of 10 meter wind factor
  isfcr=15    ! index of surface roughness
  ilone=16    ! index of longitude (degrees)
  ilate=17    ! index of latitude (degrees)
  izz=18      ! index of surface height

! If requested, save select data for output to diagnostic file
  if(conv_diagsave)then
     ii=0
     nchar=1
     ioff0=24
     nreal=ioff0
     if (lobsdiagsave) nreal=nreal+4*miter+1
     if (save_jacobian) then
        nnz = 0 ! number of non-zero elements in dH(x)/dx profile. It is different in "setupq.f90"
        nind = 0
        call new(dhx_dx, nnz, nind)
        nreal = nreal + size(dhx_dx)
     endif
     allocate(cdiagbuf(nobs),rdiagbuf(nreal,nobs))
     if(netcdf_diag) call init_netcdf_diag_
  end if

  mm1=mype+1
  scale=one
  rsig=float(nsig)

  do i=1,nobs
     muse(i)=nint(data(iuse,i)) <= jiter
  end do
  var_jb=zero
  
! prepare observation
  call dtime_setup()
  do i=1,nobs
     rwgt=one
     dtime=data(itime,i)
     call dtime_check(dtime, in_curbin, in_anybin)
     if(.not.in_anybin) cycle

     if(in_curbin) then
        dlat=data(ilat,i)
        dlon=data(ilon,i)
        dpres=data(ihgt,i)
        error=data(ier2,i)
        slat=data(ilate,i)*deg2rad
        ikx = nint(data(ikxx,i))
        itype=ictype(ikx)
        rstation_id  = data(id,i)
        reflob=data(ireflob,i)
        obserror = max(cermin(ikx),min(cermax(ikx),data(ier,i)))
     endif

!    Link observation to appropriate observation bin
     if (nobs_bins>1) then
        ibin = NINT( dtime/hr_obsbin ) + 1
     else
        ibin = 1
     endif
     IF (ibin<1.OR.ibin>nobs_bins) write(6,*)mype,'Error nobs_bins,ibin= ',nobs_bins,ibin

!    Link obs to diagnostics structure
     if (luse_obsdiag) then
        if (.not.lobsdiag_allocated) then
           if (.not.associated(obsdiags(i_refl_ob_type,ibin)%head)) then
              obsdiags(i_refl_ob_type,ibin)%n_alloc = 0
              allocate(obsdiags(i_refl_ob_type,ibin)%head,stat=istat)
              if (istat/=0) then
                 write(6,*)'setuprefl: failure to allocate obsdiags',istat
                 call stop2(286)
              end if
              obsdiags(i_refl_ob_type,ibin)%tail => obsdiags(i_refl_ob_type,ibin)%head
           else
              allocate(obsdiags(i_refl_ob_type,ibin)%tail%next,stat=istat)
              if (istat/=0) then
                 write(6,*)'setuprefl: failure to allocate obsdiags',istat
                 call stop2(286)
              end if
              obsdiags(i_refl_ob_type,ibin)%tail => obsdiags(i_refl_ob_type,ibin)%tail%next
           end if
           obsdiags(i_refl_ob_type,ibin)%n_alloc = obsdiags(i_refl_ob_type,ibin)%n_alloc +1
    
           allocate(obsdiags(i_refl_ob_type,ibin)%tail%muse(miter+1))
           allocate(obsdiags(i_refl_ob_type,ibin)%tail%nldepart(miter+1))
           allocate(obsdiags(i_refl_ob_type,ibin)%tail%tldepart(miter))
           allocate(obsdiags(i_refl_ob_type,ibin)%tail%obssen(miter))
           obsdiags(i_refl_ob_type,ibin)%tail%indxglb=ioid(i)
           obsdiags(i_refl_ob_type,ibin)%tail%nchnperobs=-99999
           obsdiags(i_refl_ob_type,ibin)%tail%luse=luse(i)
           obsdiags(i_refl_ob_type,ibin)%tail%muse(:)=.false.
           obsdiags(i_refl_ob_type,ibin)%tail%nldepart(:)=-huge(zero)
           obsdiags(i_refl_ob_type,ibin)%tail%tldepart(:)=zero
           obsdiags(i_refl_ob_type,ibin)%tail%wgtjo=-huge(zero)
           obsdiags(i_refl_ob_type,ibin)%tail%obssen(:)=zero
    
           n_alloc(ibin) = n_alloc(ibin) +1
           my_diag => obsdiags(i_refl_ob_type,ibin)%tail
           my_diag%idv = is
           my_diag%iob = ioid(i)
           my_diag%ich = 1
           my_diag%elat= data(ilate,i)
           my_diag%elon= data(ilone,i)
    
        else
           if (.not.associated(obsdiags(i_refl_ob_type,ibin)%tail)) then
              obsdiags(i_refl_ob_type,ibin)%tail => obsdiags(i_refl_ob_type,ibin)%head
           else
              obsdiags(i_refl_ob_type,ibin)%tail => obsdiags(i_refl_ob_type,ibin)%tail%next
           end if
           if (.not.associated(obsdiags(i_refl_ob_type,ibin)%tail)) then
              call die(myname,'.not.associated(obsdiags(i_refl_ob_type,ibin)%tail)')
           end if
           if (obsdiags(i_refl_ob_type,ibin)%tail%indxglb/=ioid(i)) then
              write(6,*)'setuprefl: index error'
              call stop2(288)
           end if
        endif
     endif

     if(.not.in_curbin) cycle

!    Interpolate log(surface pressure),  
!    log(pres) at mid-layers, and geopotenital height to 
!    observation location.

     call tintrp2a11(ges_z,zsges,dlat,dlon,dtime,hrdifsig,&
          mype,nfldsig)
     if(zsges>=dpres)then
        if(yn_1stprint1) write(6,*) trim(myname),': zsges (terrain) = ',zsges,'is greater than zob ',dpres,' [m]. Rejecting ob. (ONLY print for 1 OBS, there may be many)'
        yn_1stprint1=.false.
        cycle
     endif
     dpres=dpres-zsges ! OBS height above the surface, [m]
     call tintrp2a11(ges_ps,psges,dlat,dlon,dtime,hrdifsig,&
          mype,nfldsig)
     call tintrp2a1(ges_lnprsl,prsltmp,dlat,dlon,dtime,hrdifsig,&
          nsig,mype,nfldsig) ! ges_lnprsl: log(layer midpoint pressure)
     call tintrp2a1(geop_hgtl,hges,dlat,dlon,dtime,hrdifsig,&
          nsig,mype,nfldsig) ! geop_hgtl: guess geopotential height at mid-layers

!    Convert geopotential height at layer midpoints to geometric height using
!    equations (17, 20, 23) in MJ Mahoney's note "A discussion of various
!    measures of altitude" (2001).  Available on the web at
!    http://mtp.jpl.nasa.gov/notes/altitude/altitude.html
!
!    termg  = equation 17
!    termr  = equation 21
!    termrg = first term in the denominator of equation 23
!    zges   = equation 23
     sin2  = sin(slat)*sin(slat)
     termg = grav_equator * &
          ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
     termr = semi_major_axis /(one + flattening + grav_ratio -  &
          two*flattening*sin2)
     termrg = (termg/grav)*termr
     do k=1,nsig
        zges(k) = (termr*hges(k)) / (termrg-hges(k))  ! eq (23), geometric height
     end do

!    Given observation height, (1) adjust 10 meter wind factor if
!    necessary, (2) convert height to grid relative units, (3) compute
!    compute observation pressure (for diagnostic purposes only), and
!    (4) compute location of midpoint of first model layer above surface
!    in grid relative units

!    Convert observation height (in dpres) from meters to grid relative
!    units.  Save the observation height in zob for later use.
     zob = dpres                     ! zob:   [m], above the surface
     call grdcrd1(dpres,zges,nsig,1) ! dpres: [grid unit]

!    Set indices of model levels below (k1) and above (k2) observation.
     k=dpres ! k is an integer
     k1=max(1,k)
     k2=min(k+1,nsig)

!    Compute observation pressure (only used for diagnostics)
     if(k2>k1) then    !???????????? to fix problem where k1=k2, which should only happen if k1=k2=nsig !wangjia: another possibility is that k1=k2=1
        dz     = zges(k2)-zges(k1)
        dlnp   = prsltmp(k2)-prsltmp(k1)
        pobl   = prsltmp(k1) + (dlnp/dz)*(zob-zges(k1)) ! ln(pob)
     else if(k2<k1) then
        if(yn_1stprint2) write(6,*) trim(myname),' WARNING: k1>k2,[k1,k2]=[',k1,',',k2,' (ONLY print for 1 OBS, there may be many)' ! actually, no pts
        yn_1stprint2=.false.
        cycle
     else if(k2==k1) then
        if(yn_1stprint3) write(6,*) trim(myname),' WARNING: k1=k2=',k1,' (ONLY print for 1 OBS, there may be many)'
        yn_1stprint3=.false.
        ! write(6,*) trim(myname),' WARNING: ictype,k,k1,k2,nsig,zob,zges(k1),prsltmp(k1)=',ictype(ikx),k,k1,k2,nsig,zob,zges(k1),prsltmp(k1)
        ! pobl   = prsltmp(k1)
        cycle
     end if
        
     presrefl  = ten*exp(pobl) ! [mb]

!    Determine location in terms of grid units for midpoint of
!    first layer above surface
     sfcchk=log(psges)
     call grdcrd1(sfcchk,prsltmp,nsig,-1)

!    Check to see if observation is below midpoint of first
!    above surface layer.  If so, set rlow to that difference
     rlow=max(sfcchk-dpres,zero) ! [grid unit]

!    Check to see if observation is above midpoint of layer
!    at the top of the model.  If so, set rhgh to that difference.
     rhgh=max(dpres-r0_001-nsig,zero)

!    Increment obs counter along with low and high obs counters
     if(luse(i))then
        awork(1)=awork(1)+one
        if(rhgh/=zero) awork(2)=awork(2)+one ! # of OBS above TOP (mid)
        if(rlow/=zero) awork(3)=awork(3)+one ! # of OBS below BOT (mid)
     end if
     
!    Adjust observation error.

!    Increase error for observations over high topography
     factelv=one

!    Increase error if model and observation topography too different
     factdif=one
     
     ratio_errors = factdif*factelv*error/(abs(data(ier,i) + 1.0e6_r_kind*rhgh + r8*rlow)) ! ratio_errors=origional_error/inflated_error
     error = one/error
!    Check to see if observations is above the top of the model (regional mode)
     if(dpres < zero .or. dpres > rsig)ratio_errors = zero

!    Interpolate guess u, v, and w to observation location and time.
     call tintrp31(ges_refl,reflges,dlat,dlon,dpres,dtime,&
          hrdifsig,mype,nfldsig)
     !call tintrp2a1(ges_refl,reflgesprofile,dlat,dlon,dtime,hrdifsig,&
     !     nsig,mype,nfldsig)

     ddiff = reflob - reflges

!    If requested, setup for single obs test.
     if(oneobtest) then
        ddiff=maginnov
        error=one/magoberr
        ratio_errors=one
        !Vr=ddiff+rwwind !?????
     end if

!    Gross error checks
! cermin/cermax does NOT affect the actual error, only affect the calculation of "ratio"
     obserror = one/max(ratio_errors*error,tiny_r_kind) ! inflated error
     obserrlm = max(cermin(ikx),min(cermax(ikx),obserror))
     residual = abs(ddiff)
     ratio    = residual/obserrlm ! abs(O-B)/inflated_error
     qcgross  = cgross(ikx)

     if (ratio > qcgross .or. ratio_errors < tiny_r_kind) then
        if (luse(i)) awork(4) = awork(4)+one
        error = zero
        ratio_errors = zero
     end if
     
     if (ratio_errors*error <=tiny_r_kind) muse(i)=.false.
     if (nobskeep>0.and.luse_obsdiag) muse(i)=obsdiags(i_refl_ob_type,ibin)%tail%muse(nobskeep)

!   Oberror Tuning and Perturb Obs
! when "oberror_tune=.True." or "perturb_obs=.True.", "ddiff" will be changed
     
     val     = error*ddiff ! (O-B)/org_error

!    Compute penalty terms (linear & nonlinear qc).
     if(luse(i))then
        val2     = val*val
        exp_arg  = -half*val**2
        rat_err2 = ratio_errors**2
        ! njqc  -  When true, use Purser''s non linear QC
        ! vqc   -  when true, use ECMWF's non linear QC
        ! for now, njqc & vqc=.FALSE.
        ! In "setupt.f90" : else if (vqc .and. cvar_pg(ikx)> tiny_r_kind .and. error >tiny_r_kind)
        ! In "setuprw.f90":      if (cvar_pg(ikx) > tiny_r_kind .and. error > tiny_r_kind)
        if (vqc .and. cvar_pg(ikx) > tiny_r_kind .and. error > tiny_r_kind) then
           arg  = exp(exp_arg)
           wnotgross= one-cvar_pg(ikx)
           cg_refl=cvar_b(ikx)
           wgross = cg_term*cvar_pg(ikx)/(cg_refl*wnotgross)
           term = log((arg+wgross)/(one+wgross))
           wgt  = one-wgross/(arg+wgross)
           rwgt = wgt/wgtlim
        else
           term = exp_arg
           wgt  = one ! wangjia: wgtlim. Here, use "one" follows "setupt.f90"
           rwgt = wgt/wgtlim
        endif
        valqc = -two*rat_err2*term

! wangjia: write out data to be compared with "check_bufr.f90"
! zob: OBS height above the surface, [m]
! zsges: surface height;        zob+zsges=data(ihgt,i) 
! zges(nsig): model top height
! presrefl: obs pressure [hPa]
! obserror: inflated error
! write(6,*) "kkkkkkk", data(ilone,i), data(ilate,i), data(ihgt,i), data(ireflob,i), data(ier,i), data(ier2,i), cermin(ikx), cermax(ikx), zob, zsges, zges(nsig), presrefl, ddiff, obserror, ratio, qcgross, val, ddiff, reflob, reflges
        
!       Accumulate statistics for obs belonging to this task
        if (muse(i)) then
           if(rwgt < one) awork(21) = awork(21)+one
           jsig = dpres
           jsig=max(1,min(jsig,nsig))
           awork(6*nsig+jsig+100)=awork(6*nsig+jsig+100)+val2*rat_err2
           awork(5*nsig+jsig+100)=awork(5*nsig+jsig+100)+one
           awork(3*nsig+jsig+100)=awork(3*nsig+jsig+100)+valqc
        end if

!       Loop over pressure level groupings and obs to accumulate
!       statistics as a function of observation type.
        ress  = scale*ddiff
        ressw = ress*ress
        nn=1
        if (.not. muse(i)) then
           nn=2
           if(ratio_errors*error >=tiny_r_kind)nn=3
        end if
        do k = 1,npres_print
           if(presrefl >ptop(k) .and. presrefl<=pbot(k))then
              bwork(k,ikx,1,nn) = bwork(k,ikx,1,nn)+one            ! count
              bwork(k,ikx,2,nn) = bwork(k,ikx,2,nn)+ddiff          ! bias
              bwork(k,ikx,3,nn) = bwork(k,ikx,3,nn)+ressw          ! (o-g)**2
              bwork(k,ikx,4,nn) = bwork(k,ikx,4,nn)+val2*rat_err2  ! penalty
              bwork(k,ikx,5,nn) = bwork(k,ikx,5,nn)+valqc          ! nonlin qc penalty
              
           end if
        end do
     end if

     if (luse_obsdiag) then
        obsdiags(i_refl_ob_type,ibin)%tail%muse(jiter)=muse(i)
        obsdiags(i_refl_ob_type,ibin)%tail%nldepart(jiter)=ddiff
        obsdiags(i_refl_ob_type,ibin)%tail%wgtjo= (error*ratio_errors)**2
     endif
     
!    If obs is "acceptable", load array with obs info for use
!    in inner loop minimization (int* and stp* routines)
     if ( .not. last .and. muse(i)) then

        allocate(my_head)
        m_alloc(ibin) = m_alloc(ibin) +1
        my_node => my_head        ! this is a workaround
        call obsLList_appendNode(reflhead(ibin),my_node)
        my_node => null()

        my_head%idv = is
        my_head%iob = ioid(i)
        my_head%elat= data(ilate,i)
        my_head%elon= data(ilone,i)

!       Set (i,j,k) indices of guess gridpoint that bound obs location
        my_head%dlev = dpres
        call get_ijk(mm1,dlat,dlon,dpres,my_head%ij,my_head%wij)

        my_head%res     = ddiff
        my_head%err2    = error**2
        my_head%raterr2 = ratio_errors**2  
        my_head%time    = dtime
        my_head%luse    = luse(i)
        my_head%b       = cvar_b(ikx)
        my_head%pg      = cvar_pg(ikx)
        my_head%jb      = var_jb ! wangjia: not used

        if (luse_obsdiag) then
           my_head%diags => obsdiags(i_refl_ob_type,ibin)%tail
        
           my_diag => my_head%diags
           if(my_head%idv /= my_diag%idv .or. &
              my_head%iob /= my_diag%iob ) then
              call perr(myname,'mismatching %[head,diags]%(idv,iob,ibin) =', &
                        (/is,ioid(i),ibin/))
              call perr(myname,'my_head%(idv,iob) =',(/my_head%idv,my_head%iob/))
              call perr(myname,'my_diag%(idv,iob) =',(/my_diag%idv,my_diag%iob/))
              call die(myname)
           endif
        endif

        my_head => null()
     endif

!    Save select output for diagnostic file
     if(conv_diagsave .and. luse(i)) then
        ii=ii+1
        rstation_id = data(id,i)
        err_input   = data(ier2,i)
        err_adjst   = data(ier,i)
        if (ratio_errors*error>tiny_r_kind) then
           err_final = one/(ratio_errors*error)
        else
           err_final = huge_single
        endif

        errinv_input = huge_single
        errinv_adjst = huge_single
        errinv_final = huge_single
        if (err_input>tiny_r_kind) errinv_input = one/err_input
        if (err_adjst>tiny_r_kind) errinv_adjst = one/err_adjst
        if (err_final>tiny_r_kind) errinv_final = one/err_final

        if (binary_diag) call contents_binary_diag_
        if (netcdf_diag) call contents_netcdf_diag_

     end if
  end do

! Release memory of local guess arrays
  call final_vars_

! Write information to diagnostic file
  if(conv_diagsave)then
     if(netcdf_diag) call nc_diag_write
     if(binary_diag .and. ii>0)then
        call dtime_show(myname,'diagsave:refl',i_refl_ob_type)
        write(7)' refl',nchar,nreal,ii,mype,ioff0
        write(7)cdiagbuf(1:ii),rdiagbuf(:,1:ii)
        deallocate(cdiagbuf,rdiagbuf)
     end if
  end if

! End of routine

  return
  contains

  subroutine check_vars_ (proceed, include_refl)
  logical,intent(inout) :: proceed
  logical,intent(inout) :: include_refl
  integer(i_kind) ivar, istatus
! Check to see if required guess fields are available
  call gsi_metguess_get ('var::ps', ivar, istatus )
  proceed=ivar>0
  call gsi_metguess_get ('var::z' , ivar, istatus )
  proceed=proceed.and.ivar>0
  call gsi_metguess_get ('var::u' , ivar, istatus )
  proceed=proceed.and.ivar>0
  call gsi_metguess_get ('var::v' , ivar, istatus )
  proceed=proceed.and.ivar>0
  call gsi_metguess_get ('var::refl' , ivar, istatus )
  proceed=proceed.and.ivar>0
  if (ivar>0) then
     include_refl=.true.
  else
     include_refl=.false.
  endif
  end subroutine check_vars_ 

  subroutine init_vars_

  real(r_kind),dimension(:,:  ),pointer:: rank2=>NULL()
  real(r_kind),dimension(:,:,:),pointer:: rank3=>NULL()
  character(len=5) :: varname
  integer(i_kind) ifld, istatus

! If require guess vars available, extract from bundle ...
  if(size(gsi_metguess_bundle)==nfldsig) then
!    get ps ...
     varname='ps'
     call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank2,istatus)
     if (istatus==0) then
         if(allocated(ges_ps))then
            write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
            call stop2(999)
         endif
         allocate(ges_ps(size(rank2,1),size(rank2,2),nfldsig))
         ges_ps(:,:,1)=rank2
         do ifld=2,nfldsig
            call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank2,istatus)
            ges_ps(:,:,ifld)=rank2
         enddo
     else
         write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
         call stop2(999)
     endif
!    get z ...
     varname='z'
     call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank2,istatus)
     if (istatus==0) then
         if(allocated(ges_z))then
            write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
            call stop2(999)
         endif
         allocate(ges_z(size(rank2,1),size(rank2,2),nfldsig))
         ges_z(:,:,1)=rank2
         do ifld=2,nfldsig
            call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank2,istatus)
            ges_z(:,:,ifld)=rank2
         enddo
     else
         write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
         call stop2(999)
     endif
!    get refl ...
     if(include_refl) then
        varname='refl'
        call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank3,istatus)
        if (istatus==0) then
            if(allocated(ges_refl))then
               write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
               call stop2(999)
            endif
            allocate(ges_refl(size(rank3,1),size(rank3,2),size(rank3,3),nfldsig))
            ges_refl(:,:,:,1)=rank3
            do ifld=2,nfldsig
               call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank3,istatus)
               ges_refl(:,:,:,ifld)=rank3
            enddo
        else
            write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle,ier= ',istatus
            call stop2(999)
        endif
     end if
  else
     write(6,*) trim(myname), ': inconsistent vector sizes (nfldsig,size(metguess_bundle) ',&
                 nfldsig,size(gsi_metguess_bundle)
     call stop2(999)
  endif
  end subroutine init_vars_

  subroutine init_netcdf_diag_
  character(len=80) string
  character(len=128) diag_conv_file
  integer(i_kind) ncd_fileid,ncd_nobs
  logical append_diag
  logical,parameter::verbose=.false. 
     write(string,900) jiter
900  format('conv_refl_',i2.2,'.nc4')
     diag_conv_file=trim(dirname) // trim(string)

     inquire(file=diag_conv_file, exist=append_diag)

     if (append_diag) then
        call nc_diag_read_init(diag_conv_file,ncd_fileid)
        ncd_nobs = nc_diag_read_get_dim(ncd_fileid,'nobs')
        call nc_diag_read_close(diag_conv_file)

        if (ncd_nobs > 0) then
           if(verbose) print *,'file ' // trim(diag_conv_file) // ' exists.  Appending.  nobs,mype=',ncd_nobs,mype
        else
           if(verbose) print *,'file ' // trim(diag_conv_file) // ' exists but contains no obs.  Not appending. nobs,mype=',ncd_nobs,mype
           append_diag = .false. ! if there are no obs in existing file, then do not try to append
        endif
     end if

     call nc_diag_init(diag_conv_file, append=append_diag)

     if (.not. append_diag) then ! don't write headers on append - the module will break?
        call nc_diag_header("date_time",ianldate )
        call nc_diag_header("Number_of_state_vars", nsdim          )
     endif
  end subroutine init_netcdf_diag_
  subroutine contents_binary_diag_
        cdiagbuf(ii)    = station_id         ! station id

        rdiagbuf(1,ii)  = ictype(ikx)        ! observation type
        rdiagbuf(2,ii)  = icsubtype(ikx)     ! observation subtype
    
        rdiagbuf(3,ii)  = data(ilate,i)      ! observation latitude (degrees)
        rdiagbuf(4,ii)  = data(ilone,i)      ! observation longitude (degrees)
        rdiagbuf(5,ii)  = rmiss_single       ! station elevation (meters)
        rdiagbuf(6,ii)  = presrefl           ! observation pressure (hPa)
        rdiagbuf(7,ii)  = data(ihgt,i)       ! observation height (meters)
        rdiagbuf(8,ii)  = dtime-time_offset  ! obs time (hours relative to analysis time)

        rdiagbuf(9,ii)  = data(iqc,i)        ! input prepbufr qc or event mark
        rdiagbuf(10,ii) = rmiss_single       ! setup qc or event mark
        rdiagbuf(11,ii) = data(iuse,i)       ! read_prepbufr data usage flag
        if(muse(i)) then
           rdiagbuf(12,ii) = one             ! analysis usage flag (1=use, -1=not used)
        else
           rdiagbuf(12,ii) = -one
        endif

        rdiagbuf(13,ii) = rwgt               ! nonlinear qc relative weight
        rdiagbuf(14,ii) = errinv_input       ! prepbufr inverse obs error (m/s)**-1
        rdiagbuf(15,ii) = errinv_adjst       ! read_prepbufr inverse obs error (m/s)**-1
        rdiagbuf(16,ii) = errinv_final       ! final inverse observation error (m/s)**-1

        rdiagbuf(17,ii) = data(ireflob,i)      ! radial wind speed observation (m/s)
        rdiagbuf(18,ii) = ddiff              ! obs-ges used in analysis (m/s)
        rdiagbuf(19,ii) = data(ireflob,i)-reflges  ! obs-ges w/o bias correction (m/s) (future slot)



        rdiagbuf(23,ii) = 1.e+10_r_single    ! ges ensemble spread (filled in EnKF)
        rdiagbuf(24,ii) = 1.e+10_r_single    ! ges ensemble spread (filled in EnKF)

        ioff=ioff0
        if (lobsdiagsave) then
           do jj=1,miter 
              ioff=ioff+1
              if (obsdiags(i_refl_ob_type,ibin)%tail%muse(jj)) then
                 rdiagbuf(ioff,ii) = one
              else
                 rdiagbuf(ioff,ii) = -one
              endif
           enddo
           do jj=1,miter+1
              ioff=ioff+1
              rdiagbuf(ioff,ii) = obsdiags(i_refl_ob_type,ibin)%tail%nldepart(jj)
           enddo
           do jj=1,miter
              ioff=ioff+1
              rdiagbuf(ioff,ii) = obsdiags(i_refl_ob_type,ibin)%tail%tldepart(jj)
           enddo
           do jj=1,miter
              ioff=ioff+1
              rdiagbuf(ioff,ii) = obsdiags(i_refl_ob_type,ibin)%tail%obssen(jj)
           enddo
        endif
        if (save_jacobian) then
           call writearray(dhx_dx, rdiagbuf(ioff+1:nreal,ii))
           ioff = ioff + size(dhx_dx)
        endif

  end subroutine contents_binary_diag_
  subroutine contents_netcdf_diag_
! Observation class
  character(7),parameter     :: obsclass = '     refl'
  real(r_kind),dimension(miter) :: obsdiag_iuse
           call nc_diag_metadata("Station_ID",              station_id             )
           call nc_diag_metadata("Observation_Class",       obsclass               )
           call nc_diag_metadata("Observation_Type",        ictype(ikx)            )
           call nc_diag_metadata("Observation_Subtype",     icsubtype(ikx)         )
           call nc_diag_metadata("Latitude",                sngl(data(ilate,i))    )
           call nc_diag_metadata("Longitude",               sngl(data(ilone,i))    )
           call nc_diag_metadata("Station_Elevation",       sngl(zero)             )
           call nc_diag_metadata("Pressure",                sngl(presrefl)         )
           call nc_diag_metadata("Height",                  sngl(data(ihgt,i))     )
           call nc_diag_metadata("Time",                    sngl(dtime-time_offset))
           call nc_diag_metadata("Prep_QC_Mark",            sngl(data(iqc,i))      )
           call nc_diag_metadata("Prep_Use_Flag",           sngl(data(iuse,i))     )
!          call nc_diag_metadata("Nonlinear_QC_Var_Jb",     sngl(var_jb)           )
           call nc_diag_metadata("Nonlinear_QC_Rel_Wgt",    sngl(rwgt)             )                 
           if(muse(i)) then
              call nc_diag_metadata("Analysis_Use_Flag",    sngl(one)              )
           else
              call nc_diag_metadata("Analysis_Use_Flag",    sngl(-one)             )              
           endif

           call nc_diag_metadata("Errinv_Input",            sngl(errinv_input)     )
           call nc_diag_metadata("Errinv_Adjust",           sngl(errinv_adjst)     )
           call nc_diag_metadata("Errinv_Final",            sngl(errinv_final)     )

           call nc_diag_metadata("Observation",                   sngl(data(ireflob,i))  )
           call nc_diag_metadata("Obs_Minus_Forecast_adjusted",   sngl(ddiff)          )
           call nc_diag_metadata("Obs_Minus_Forecast_unadjusted", sngl(data(ireflob,i)-reflges) )
 
           if (lobsdiagsave) then
              do jj=1,miter
                 if (obsdiags(i_refl_ob_type,ibin)%tail%muse(jj)) then
                       obsdiag_iuse(jj) =  one
                 else
                       obsdiag_iuse(jj) = -one
                 endif
              enddo
   
              call nc_diag_data2d("ObsDiagSave_iuse",     obsdiag_iuse                             )
              call nc_diag_data2d("ObsDiagSave_nldepart", obsdiags(i_refl_ob_type,ibin)%tail%nldepart )
              call nc_diag_data2d("ObsDiagSave_tldepart", obsdiags(i_refl_ob_type,ibin)%tail%tldepart )
              call nc_diag_data2d("ObsDiagSave_obssen",   obsdiags(i_refl_ob_type,ibin)%tail%obssen   )             
           endif
           if (save_jacobian) then
              call fullarray(dhx_dx, dhx_dx_array)
              call nc_diag_data2d("Observation_Operator_Jacobian", dhx_dx_array)
           endif
   
  end subroutine contents_netcdf_diag_

  subroutine final_vars_
    if(allocated(ges_refl )) deallocate(ges_refl )
    if(allocated(ges_z    )) deallocate(ges_z )
    if(allocated(ges_ps   )) deallocate(ges_ps)
  end subroutine final_vars_

end subroutine setuprefl
