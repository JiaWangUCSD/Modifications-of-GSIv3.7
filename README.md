These codes are modifications made into GSIv3.7 to enable direct asimilation of reflectivity following the method proposed by Wang and Wang (2017).

Wang, Yongming, and Xuguang Wang. "Direct assimilation of radar reflectivity without tangent linear and adjoint of the nonlinear observation operator in the GSI-based EnVar system: Methodology and experiment with the 8 May 2003 Oklahoma City tornadic supercell." Monthly Weather Review 145.4 (2017): 1447-1471.

1.	To augment the control variables: 
(1)	anavinfo_arw_netcdf: specify new state and control variables; 
(2)	cplr_wrf_netcdf_interface.f90: extract new state variables from the background; 
(3)	cplr_get_wrf_mass_ensperts.f90: extract new state variables from the ensemble members; 
(4)	update_guess.f90: update guess fields for new variables; 
(5)	cplr_wrwrfmassa.f90 and cplr_wrf_netcdf_interface.f90: write new variables into the final analysis; 
(6)	ensctl2state.f90, ensctl2state_ad.f90, control2state.f90, and control2state_ad.f90: add new variables.

3.	To assimilate reflectivity observations:
(1)	read_obs.F90, setuprhsall.f90, m_obsdiags.F90, m_obsNodeTypeManager.F90, m_obsHeadBundle.F90, obsmod.F90, and intjo.f90: add a block related to reflectivity observations; 
(2)	read_radarref_mosaic_gridrad.f90 (added): read three-dimensional reflectivity in the BUFR format; 
(3)	m_reflNode.F90 (added): defines “reflNode”, which stores information related reflectivity observations; 
(4)	setuprefl.f90 (added): calculate reflectivity observation innovations; 
(5)	intrefl.f90 (added): reflectivity observation operator and its adjoint; 
(6)	stprefl.f90 (added): compute the penalty term related to reflectivity observations.
