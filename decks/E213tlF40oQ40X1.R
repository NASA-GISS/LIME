E213tlF40oQ40X1.R GISS Model E  atm-ocean     M. Kelley       03/14/2019

E213F40oQ40: E200F40oQ40 + updated aerosol/ozone input files based on E25TomaF40pi

E200F40oQ40: E199F40oQ40 + 5m maximum for lake ice (LAKE_ICE_MAX=5.) 
E199F40oQ40: E198F40oQ40 + no leads (OPNOCN=0.; SEAICE_E199.f); 
             AIC=1JAN3901.rsfE198F40oQ40.nc 
E198F40oQ40: E197F40oQ40 + alternative OCNMESO_DRV 
E197F40oQ40: E196F40oQ40 + alternative TOPO_OC 
E196F40oQ40: E193F40oQ32 + 40 vertical layers in ocean; 
             adjusted straits, salinity, tides 
E193F40oQ32: E190F40oQ32 + new aer.input data based on E08TomaF40pi
E190F40oQ32: GISS coupled model, based on E137F40pi

3/27/2017 updated to latest F40 atm (aerosol/ozone inputs etc.)

M. Kelley 08/12/2016 cloned E4F40 and added 1-degree
ocean R with various new mixing options activated.
This template should be considered temporary until the
new template system has been finalized.

Preprocessor Options
#define STDHYB                  ! standard hybrid vertical coordinate
#define ATM_LAYERING L40        ! 40 layers, top at .1 mb
#define NEW_IO                  ! new I/O (netcdf) on
#define NEW_IO_SUBDD
#define CACHED_SUBDD
#define IRRIGATION_ON
#define SWFIX_20151201
#define NO_HDIURN               ! exclude hdiurn diagnostics
#define MODIS_LAI
#define CHECK_OCEAN                  ! needed to compile aux/file CMPE002
#define SIMPLE_MESODIFF
#define OCN_LAYERING L40_5008m

#define ODIFF_FIXES_2017
#define EXPEL_COASTAL_ICEXS
#define NEW_BCdalbsn
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
Atm144x90                           ! horizontal resolution is 144x90 -> 2x2.5deg
AtmLayering                         ! vertical resolution
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform
ORES_1Qx1 OFFT288E                  ! ocean horiz res 1.25x1deg

IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

    ! lat-lon grid specific source codes
AtmRes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
MODEL_COM                           ! calendar, timing variables
MODELE_DRV                          ! ModelE cap
MODELE                              ! initialization and main loop
ATM_COM                             ! main atmospheric variables
ATM_DRV                             ! driver for atmosphere-grid components
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
IRRIGMOD                            ! irrigation module
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
!SEAICE_E200 SEAICE_DRV                   ! seaice modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV     ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO ocalbedo ! radiation and albedo
DIAG_COM DIAG DEFACC                ! diagnostics
OCN_DRV                             ! driver for ocean-grid components
OLAYERS                             ! ocean layering options
ODIAG_COM  OCEAN_COM  OGEOM         ! dynamic ocean modules
OCNDYN  OCNDYN2  OTIDELL            ! dynamic ocean routines
OCN_Interp                          ! dynamic ocean routines
OSTRAITS_COM  OSTRAITS              ! straits parameterization
OCNKPP                              ! ocean vertical mixing
!OCNMESO_DRV_E198 OCNTDMIX OCNGM          ! ocean mesoscale mixing
OCNMESO_DRV OCNTDMIX OCNGM          ! ocean mesoscale mixing
OCEANR_DIM  OFLUXES
ODIAG_PRT                           ! ocean diagnostic print out
OCNFUNTAB                           ! ocean function look up table
SparseCommunicator_mod              ! sparse gather/scatter module
OCNQUS                              ! QUS advection scheme
OCNGISS_TURB                        ! GISS Turbulence vertical mixing scheme
OCNGISS_SM                          ! GISS Sub-mesoscale mixing scheme
ODIAG_ZONAL
SUBDD
LIME

OCN_Int_LATLON                      ! atm-ocn regrid routines

Components:
shared MPI_Support solvers giss_LSM 
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT 
! OPTS_dd2d = NC_IO=PNETCDF

Data input files:
    ! start from the restart file of an earlier run ...                 ISTART=9
AIC=1JAN7501.rsfE200F40oQ40.nc

OFTAB=OFTABLE_NEW                   ! ocean function table
KBASIN=KB288X180.modelE.BS1.nc      ! ocean basin designations (1 cell Bering Strait)
TOPO_OC=altocnbc288x180_20170717/OZ1QX1N.BS1.BAB.PG.HB.GIB.SIC.nc    ! ocean fraction and topography
! TOPO_OC=altocnbc288x180_20170717/OZ1QX1N.BS1.BAB.PG.HB.GIB.nc        ! ocean fraction and topography
TIDES=TIDAL_e_v2_1QX1               ! ocean bottom tidal energy and velocity squared
OSTRAITS=altocnbc288x180_20170717/OSTRAITS_288x180_zspec_BAB.T4.nml  ! parameterized straits info
TOPO=Z2HX2fromZ1QX1N.BS1.nc        ! surface fractions and topography (1 cell Bering Strait)
ICEDYN_MASKFAC=iceflowmask_144x90.nc

TDISS=altocnbc288x180_20170717/TIDAL_e_v2_1QX1.HB.nc
TDISS_N=tdiss/Jayne2009_288x180.nc
POROS=altocnbc288x180_20170717/poros.nc

RVR=RD_Fd.nc             ! river direction file
NAMERVR=RD_Fd.names.txt  ! named river outlets

CDN=CD144X90.ext.nc
VEG=V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.ext.nc
LAIMAX=V144x90_EntMM16_lai_max_trimmed_scaled_ext.nc
HITEent=V144x90_EntMM16_height_trimmed_scaled_ext.nc
LAI=V144x90_EntMM16_lai_trimmed_scaled_ext.nc
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
IRRIG=Irrig144x90_1848to2100_FixedFuture_v3.nc
SOIL=S144X900098M.ext.nc
TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc             ! topographic variation for gwdrag
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
GLMELT=GLMELT_144X90_gas.OCN.nc
RADN1=sgpgxg.table8                           ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800  ! rad.tables and history files
RADN4=LWCorrTables33k                         ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-2014_CMIP6_hdr
RADN8=cloud.epsilon4.72x46
!RADN9=solar.lean2015.ann1610-2014.nc ! need KSOLAR=2
RADN9=solar.CMIP6official.ann1850-2299.nc ! need KSOLAR=2
RADNE=topcld.trscat8

ISCCP=ISCCP.tautables
GHG=GHG.CMIP6.1-2014.txt  !  GreenHouse Gases for CMIP6 runs up to 2014
CO2profile=CO2profile.Jul16-2017.txt ! scaling of CO2 in stratosphere
dH2O=dH2O_by_CH4_monthly

! Begin NINT E2.1 input files

BCdalbsn=cmip6_nint_inputs_E25Toma/BCdalbsn/1850.nc
DUSTaer=cmip6_nint_inputs_E25Toma/DUST/1850.nc
TAero_SUL=cmip6_nint_inputs_E25Toma/SUL/1850.nc
TAero_SSA=cmip6_nint_inputs_E25Toma/SSA/1850.nc
TAero_NIT=cmip6_nint_inputs_E25Toma/NIT/1850.nc
TAero_OCA=cmip6_nint_inputs_E25Toma/OCA/1850.nc
TAero_BCA=cmip6_nint_inputs_E25Toma/BCA/1850.nc
TAero_BCB=cmip6_nint_inputs_E25Toma/BCB/1850.nc

O3file=cmip6_nint_inputs_E25Toma/O3/1850.nc
Ox_ref=o3_2010_shindell_144x90x49_April1850.nc

! End NINT E2.1 input files

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E213tlF40oQ40X1 (E200F40oQ32 + updated aerosol/ozone input files based on E25TomaF40pi)

&&PARAMETERS

! parameters set for coupled ocean runs:
KOCEAN=1            ! ocn is prognostic
OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)
OTIDE = 0           !  Ocean tides are not used
variable_lk=1
init_flake=1
ocean_use_qus=1     ! Advection uses the quadratic upstream scheme
DTO=112.5
ocean_use_tdmix=1  ! tdmix scheme for meso mixing
ocean_use_gmscz=1  ! vertically variation of meso diffusivity, option 1
ocean_kvismult=2.  ! mult. factor for meso diffusivity
ocean_enhance_shallow_kmeso=1 ! stronger meso mixing in shallow water
ocean_use_tdiss=1  ! simple tidally induced diapycnal diffusivity

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000055  ! threshold (1/s) for triggering deformation waves
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.1       ! default is 0.5
CDEF=1.6       ! tuning factor for deformation -> momentum flux
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF


! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! No tau adjustment factors for aerosols
FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.3,1.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.655 ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.    ! below 850mb and MC regions; tune this last  to get rad.balance
WMUI_multiplier = 2.
use_vmp=1
radius_multiplier=1.1

PTLISO=0.        ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=1850
!crops_yr=1850  ! if -1, crops in VEG-file is used
!s0_yr=1850
!s0_day=182
!ghg_yr=1850
!ghg_day=182
!irrig_yr=1850
volc_yr=-1
!volc_day=182
!aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
!albsn_yr=1850
dalbsnX=1.
!o3_yr=-1850
!aer_int_yr=1850    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=0
!orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols

MADVOL=2

DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)
! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time
SUBDD=' '        ! no sub-daily frequency diags
NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
KCOPY=1          ! 0: no output; 1: save .acc; 2: unused; 3: include ocean data
KRSF=1           ! 0: no output; X: save rsf at the beginning of every X month
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960
&&END_PARAMETERS

 &INPUTZ
 YEARI=7501,MONTHI=1,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=7502,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=9,IRANDI=0, YEARE=7501,MONTHE=1,DATEE=1,HOURE=1,
/
