subroutine read_radarref_mosaic_gridrad(nread,ndata,nodata,infile,obstype,lunout,twind,sis,hgtl_full,nobs)
!   wangjia: model "read_prepbufr"
!   input argument list:
!     infile   - unit from which to read mosaic information file
!     obstype  - observation type to process
!     lunout   - unit to which to write data for further processing
!     twind    - input group time window (hours)
!     sis      - observation variable name
!     hgtl_full- 3d geopotential height on full domain grid
!   output argument list:
!     nread    - number of type "obstype" observations read
!     ndata    - number of type "obstype" observations retained for further processing
!     nodata   - number of type "obstype" observations retained for further processing. Normally, nodata=ndata, but for U/V observation, nodata=2*ndata 
!     nobs     - array of observations on each subdomain for each processor
!_____________________________________________________________________
!
  use kinds, only: r_single,r_kind,r_double,i_kind
  use constants, only: zero,one_tenth,one,deg2rad,fv,t0c,half,&
      three,four,rad2deg,tiny_r_kind,huge_r_kind,huge_i_kind,&
      r60inv,r10,r100,r2000
  use gridmod, only: diagnostic_reg,regional,nlon,nlat,nsig,&
      tll2xy,txy2ll,rotate_wind_ll2xy,rotate_wind_xy2ll,&
      rlats,rlons,twodvar_regional,nlon_regional,nlat_regional
  use convinfo, only: nconvtype,ctwind, &
      ncmiter,ncgroup,ncnumgrp,icuse,ictype,icsubtype,ioctype, &
      ithin_conv,rmesh_conv,pmesh_conv, &
      use_prepb_satwnd
  use convinfo, only: id_drifter
  use mod_wrfmass_to_a, only: wrfmass_obs_to_a8

  use obsmod, only: iadate,oberrflg,perturb_obs,perturb_fact,ran01dom,hilbert_curve
  use obsmod, only: blacklst,offtime_data,bmiss,ext_sonde
  use aircraftinfo, only: aircraft_t_bc,aircraft_t_bc_pof,ntail,taillist,idx_tail,npredt,predt, &
      aircraft_t_bc_ext,ntail_update,max_tail,nsort,itail_sort,idx_sort,timelist
  use converr,only: etabl
  use converr_ps,only: etabl_ps,isuble_ps,maxsub_ps
  use converr_q,only: etabl_q,isuble_q,maxsub_q
  use converr_t,only: etabl_t,isuble_t,maxsub_t
  use converr_uv,only: etabl_uv,isuble_uv,maxsub_uv
  use converr_pw,only: etabl_pw,isuble_pw,maxsub_pw
  use convb_ps,only: btabl_ps
  use convb_q,only: btabl_q
  use convb_t,only: btabl_t
  use convb_uv,only: btabl_uv
  use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen,thin4d
  use qcmod, only: errormod,errormod_aircraft,noiqc,newvad,njqc
  use qcmod, only: pvis,pcldch,scale_cv,estvisoe,estcldchoe,vis_thres,cldch_thres
  use nltransf, only: nltransf_forward
  use convthin, only: make3grids,map3grids,del3grids,use_all
  use blacklist, only : blacklist_read,blacklist_destroy
  use blacklist, only : blkstns,blkkx,ibcnt
  use sfcobsqc,only: init_rjlists,get_usagerj,get_gustqm,destroy_rjlists
  use sfcobsqc,only: init_gsd_sfcuselist,apply_gsd_sfcuselist,destroy_gsd_sfcuselist        
  use hilbertcurve,only: init_hilbertcurve, accum_hilbertcurve, &
                         apply_hilbertcurve,destroy_hilbertcurve
  use ndfdgrids,only: init_ndfdgrid,destroy_ndfdgrid,relocsfcob,adjust_error
  use jfunc, only: tsensible
  use deter_sfc_mod, only: deter_sfc_type,deter_sfc2
  use gsi_nstcouplermod, only: nst_gsi,nstinfo
  use gsi_nstcouplermod, only: gsi_nstcoupler_deter
  use aircraftobsqc, only: init_aircraft_rjlists,get_aircraft_usagerj,&
                           destroy_aircraft_rjlists
  use adjust_cloudobs_mod, only: adjust_convcldobs,adjust_goescldobs
  use mpimod, only: npe
  use rapidrefresh_cldsurf_mod, only: i_gsdsfc_uselist,i_gsdqc
  use gsi_io, only: verbose

  implicit none

! Declare passed variables
  character(len=*)                      ,intent(in   ) :: infile,obstype
  character(len=20)                     ,intent(in   ) :: sis
  integer(i_kind)                       ,intent(in   ) :: lunout !nrec_start
  integer(i_kind)                       ,intent(inout) :: nread,ndata,nodata
  integer(i_kind),dimension(npe)        ,intent(inout) :: nobs
  real(r_kind)                          ,intent(in   ) :: twind
  real(r_kind),dimension(nlat,nlon,nsig),intent(in   ) :: hgtl_full

  character(len=*),parameter :: myname='READ_RADARREF_MOSAIC_GRIDRAD: '
  integer(i_kind),parameter  :: mxlv=255
  real(r_kind),parameter:: r90  = 90.0_r_kind
  real(r_kind),parameter:: r360 = 360.0_r_kind
  real(r_kind),parameter:: r1_2 = 1.2_r_kind
  real(r_kind),parameter:: r0_01 = 0.01_r_kind
  real(r_kind),parameter:: r0_01_bmiss=r0_01*bmiss ! bmiss=10^9, on IBM bmiss=10^11
  character(len=80) :: hdstr,obstr,qcstr,oestr,levstr
  real(r_double),dimension(5):: hdr
  real(r_double),dimension(2,mxlv):: qcmark,obserr,obsdat
  real(r_double),dimension(1,mxlv):: levdat
  integer(i_kind),dimension(mxlv) :: zqm,rqm

!  data statements
  data hdstr  /'SID XOB YOB DHR TYP'/
  data obstr  /'ZOB HREF'/
  data qcstr  /'ZQM RQM '/
  data oestr  /'ZOE ROE '/
  data levstr /'ZOB'/

  integer(i_kind) ireadmg,ireadsb
  integer(i_kind)  :: nreal,nchanl,ilat,ilon
  character(len=8) :: sidchr,subset
  real(r_double) rstation_id,qcmark_huge
  equivalence(rstation_id,sidchr) ! equivalence to handle character names
  integer(i_kind) levs,i,j,k,ikx,itx,ithin,kx,iout
  integer(i_kind) maxobs,nmsgmax,mxtb,nmsg,ntb,lnbufr
  integer(i_kind) idate,iy,im,idd,ihh,iret
  real(r_kind) :: toff,t4dv,timeobs,time_correction,time,timex,zeps
  real(r_kind) :: dlat,dlon,dlat_earth,dlon_earth,dlat_earth_deg,dlon_earth_deg
  real(r_kind) :: roe,usage
  character(len=10) :: date
  integer(i_kind) minobs,minan
  integer(i_kind),dimension(5):: idate5
  integer(i_kind) :: qm,numobsa
  real(r_kind),allocatable,dimension(:,:):: cdata_all,cdata_out
  integer(i_kind),allocatable,dimension(:):: isort
  integer(i_kind) :: idomsfc
  real(r_kind) tsavg,ff10,sfcr,zz
  logical :: reflob,outside,inflate_error,print_verbose

  data lnbufr /10/
  data ithin / -9 /

  print_verbose=.false.
  if(verbose) print_verbose=.true.
  qcmark_huge = huge_i_kind ! huge_i_kind ~ 10^10

!!!!!! initialize variables
! nread-># of read; ndata-># of retained; nodata=ndata (unless for U/V)
! nchanl->not sure, but no changes in "read_prepbufr.f90"
nread=0; ndata=0; nodata=0; nchanl=0
ilon=2 ; ilat=3 ! the location of LON/LAT in "cdata_all", cdata_all(2,iout)=dlon, cdata_all(3,iout)=dlat 

!!!!!! check the "obstype", should be "refl"
reflob=.false.; ikx=0
do i=1,nconvtype
   if((trim(obstype) == trim(ioctype(i))) .and. (abs(icuse(i))== 1)) then
      reflob=.true.
      ikx=i
      nreal=18 ! # of info in cdata_all(nreal,maxobs)
      ithin=ithin_conv(i)
      exit
   end if
end do
! If no data requested to be process, exit routine
if(.not.reflob) then
   write(6,*) trim(myname),'CONVINFO DOES NOT INCLUDE ANY ',trim(obstype)
   return
end if

!!!!!! get message and subset counts
call getcount_bufr(infile,nmsgmax,mxtb)

!!!!!! 1st read 
call closbf(lnbufr)
open(lnbufr,file=trim(infile),form='unformatted')
call openbf(lnbufr,'IN',lnbufr)
call datelen(10)

! wangjia: ioctype,ictype: refl,999
! nmsg-># of message; ntb-># of subsets; maxobs-># of obs (each level counted)
nmsg=0; ntb=0; maxobs=0
msg_report: do while (ireadmg(lnbufr,subset,idate) == 0)
   if(nmsg == 0) call time_4dvar(idate,toff) ! toff: Time since start of 4D-Var window (hours)
   if(trim(subset)/="ADPUPA") then
      write(6,*) trim(myname),'ERROR: subset/=ADPUPA'; call stop2(50)
   end if
   nmsg=nmsg+1
   if (nmsg>nmsgmax) then
      write(6,*) trim(myname),'ERROR: messages exceed maximum ',nmsgmax; call stop2(50)
   endif
   loop_report: do while (ireadsb(lnbufr) == 0)
      ntb = ntb+1
      if (ntb>mxtb) then
         write(6,*) trim(myname),'ERROR: reports exceed maximum ',mxtb; call stop2(50)
      endif
      call ufbint(lnbufr,hdr,5,1,iret,hdstr)
      rstation_id=hdr(1)
      if(trim(sidchr)/="GRIDRAD".and.trim(sidchr)/="MRMS") then
         write(6,*) trim(myname),'ERROR: SID != GRIDRAD or MRMS'; call stop2(50)
      end if 
      call ufbint(lnbufr,levdat,1,mxlv,levs,levstr)
      maxobs=maxobs+max(1,levs)
   end do loop_report
end do msg_report

if (nmsg==0) then
   call closbf(lnbufr)
   close(lnbufr)
   if(print_verbose) write(6,*) trim(myname),'no messages/reports '
   return
end if
if(print_verbose) write(6,*) trim(myname),'messages/reports = ',nmsg,'/',ntb
if(print_verbose) write(6,*) trim(myname),'time offset is ',toff,' hours.'

allocate(cdata_all(nreal,maxobs),isort(maxobs)) ! isort is used to sort the observation, which is not necessary in reflectivity
cdata_all=zero

!!!!!! 2nd read
call closbf(lnbufr)
open(lnbufr,file=infile,form='unformatted')
call openbf(lnbufr,'IN',lnbufr)
call datelen(10)

nmsg=0; ntb=0
loop_msg: do while (ireadmg(lnbufr,subset,idate)== 0)
   nmsg = nmsg+1
   loop_readsb: do while(ireadsb(lnbufr) == 0)
      ntb = ntb+1
      !Extract type, date, and location information
      call ufbint(lnbufr,hdr,5,1,iret,hdstr)
      rstation_id=hdr(1)
      kx=hdr(5)
      if(abs(hdr(3))>r90 .or. abs(hdr(2))>r360) cycle loop_readsb
      if(hdr(2)== r360) hdr(2)=hdr(2)-r360
      if(hdr(2) < zero) hdr(2)=hdr(2)+r360
      dlon_earth_deg=hdr(2)
      dlat_earth_deg=hdr(3)
      dlon_earth=dlon_earth_deg*deg2rad
      dlat_earth=dlat_earth_deg*deg2rad

      if(regional)then
         call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)    !convert earth lon-lat to x-y grid coordinate, dlon/dlat->x/y grid coordinate (grid units)
         if(outside) cycle loop_readsb   ! check to see if outside regional domain
      else
         dlat = dlat_earth
         dlon = dlon_earth
         call grdcrd1(dlat,rlats,nlat,1)
         call grdcrd1(dlon,rlons,nlon,1)
      endif

     if(offtime_data) then
!             in time correction for observations to account for analysis
!                      time being different from obs file time.
        write(date,'( i10)') idate
        read (date,'(i4,3i2)') iy,im,idd,ihh
        idate5(1)=iy
        idate5(2)=im
        idate5(3)=idd
        idate5(4)=ihh
        idate5(5)=0
        call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
        idate5(1)=iadate(1)
        idate5(2)=iadate(2)
        idate5(3)=iadate(3)
        idate5(4)=iadate(4)
        idate5(5)=0
        call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date
!       Add obs reference time, then subtract analysis time to get obs time relative to analysis
        time_correction=float(minobs-minan)*r60inv ! [hour]
     else
        time_correction=zero
     end if

     timeobs=real(real(hdr(4),r_single),r_double)
     t4dv=timeobs + toff
     zeps=1.0e-8_r_kind
     if (t4dv<zero  .and.t4dv>      -zeps) t4dv=zero
     if (t4dv>winlen.and.t4dv<winlen+zeps) t4dv=winlen
     t4dv=t4dv + time_correction    ! Tobs-Tbgn. ??? could be WRONG
     time=timeobs + time_correction ! Tobs-Tan [h]

     if (l4dvar.or.l4densvar) then
        if (t4dv<zero.OR.t4dv>winlen) cycle loop_readsb ! outside time window
     else
        if((real(abs(time)) > real(ctwind(ikx))) .or. (real(abs(time)) > real(twind)) ) cycle loop_readsb ! outside time window
     end if
     timex=time ! timex is used in drifted OBS

     ! Extract data information on levels
     call ufbint(lnbufr,obsdat,2,mxlv,levs,obstr)
     call ufbint(lnbufr,qcmark,2,mxlv,levs,qcstr)
     call ufbint(lnbufr,obserr,2,mxlv,levs,oestr)

     nread=nread+levs
     ! get the QC marker 
     do k=1,levs
        do i=1,2
           qcmark(i,k) = min(qcmark(i,k),qcmark_huge)
        end do
        zqm(k)=nint(qcmark(1,k))
        rqm(k)=nint(qcmark(2,k))
     end do

     loop_k_levs: do k=1,levs
        if(reflob) then
           if(abs(obsdat(2,k))>r0_01_bmiss) cycle loop_k_levs ! check missing value
           if( (abs(obsdat(1,k))>r0_01_bmiss) .or. (abs(qcmark(1,k))>r0_01_bmiss) .or. (abs(qcmark(2,k))>r0_01_bmiss) .or. (abs(obserr(1,k))>r0_01_bmiss) .or. (abs(obserr(2,k))>r0_01_bmiss) ) then
              write(6,*) abs(obsdat(1:2,k)),abs(qcmark(1:2,k)),abs(obserr(1:2,k))
              write(6,*) trim(myname),'ERROR: something is missing'; call stop2(50)
           end if
           qm=rqm(k)
        end if

        ! Set usage variable
        usage = zero
        if(icuse(ikx)<=0) usage=100._r_kind

        ! Check qc marks to see if obs should be processed or skipped
        ! Special block for data thinning - if requested
        if (ithin > 0) then
           write(6,*) trim(myname),'ERROR: ithin > 0'; call stop2(50) 
        else
           ndata=ndata+1
           nodata=nodata+1
           iout=ndata
        endif
        if(ndata > maxobs) then
           write(6,*) trim(myname),'ERROR: ndata>maxobs for ',obstype; call stop2(50)
           !wangjia ndata = maxobs
        end if

        ! Get information from surface file necessary for conventional data here
        call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)

        ! Set inflate_error logical based on qm flag
        inflate_error=.false.
        !if (qm==3 .or. qm==7) inflate_error=.true.
        if(qm/=1) then
           write(6,*) trim(myname),'ERROR: reflectivity QM !=1, and QM=',qm; call stop2(50)
        end if 

        ! reflectivity OBS 
        roe=obserr(2,k)
        if (inflate_error) roe=roe*r1_2
        cdata_all(1,iout)=roe                     ! error
        cdata_all(2,iout)=dlon                    ! grid relative longitude
        cdata_all(3,iout)=dlat                    ! grid relative latitude
        cdata_all(4,iout)=obsdat(1,k)             ! obs absolute height (m) 
        cdata_all(5,iout)=obsdat(2,k)             ! reflectivity obs (dBZ)
        cdata_all(6,iout)=rstation_id             ! station id
        cdata_all(7,iout)=t4dv                    ! time [hour]
        cdata_all(8,iout)=ikx                     ! type
        cdata_all(9,iout)=rqm(k)                 ! quality mark
        cdata_all(10,iout)=obserr(2,k)            ! original obs error  
        cdata_all(11,iout)=usage                  ! usage parameter
        cdata_all(12,iout)=idomsfc                ! dominate surface type
        cdata_all(13,iout)=tsavg                  ! skin temperature
        cdata_all(14,iout)=ff10                   ! 10 meter wind factor
        cdata_all(15,iout)=sfcr                   ! surface roughness
        cdata_all(16,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
        cdata_all(17,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
        cdata_all(18,iout)=zz                     ! terrain height at ob location

      end do loop_k_levs 
   end do loop_readsb
end do loop_msg

call closbf(lnbufr)
close(lnbufr)
if(print_verbose) write(6,*) trim(myname),'closbf(',lnbufr,')'

! Write header record and data to output file for further processing
allocate(cdata_out(nreal,ndata))
do i=1,ndata
   itx=i
   do k=1,nreal
      cdata_out(k,i)=cdata_all(k,itx)
   end do
end do
deallocate(isort,cdata_all)

!!!!!! MAP from original grid to analysis grid
! count_obs-># of OBS on each involved PE
if(ndata > 0 ) then
   if(nlon==nlon_regional .and. nlat==nlat_regional) then
   else
      write(6,*) trim(myname),'nlon/nlat != nlon_regional/nlat_regional. Could be Problematic!!!!!!'
      call wrfmass_obs_to_a8(cdata_out,nreal,ndata,ilat,ilon,numobsa)
      nread=numobsa
      ndata=numobsa
      nodata=numobsa
   endif
   call count_obs(ndata,nreal,ilat,ilon,cdata_out,nobs)
! "ndata": I think it is NOT read/used, at least in the subroutine "disobs" of "obs_para.f90". And the "ndata" in "obs_para.f90" is from "read_obs"
   write(lunout) obstype,sis,nreal,nchanl,ilat,ilon,ndata ! compared with "read_radarref_mosaic", add "ndata"
   write(lunout) cdata_out
   deallocate(cdata_out)
endif

! End of routine
  return

end subroutine read_radarref_mosaic_gridrad
