!This module models the decomposition of carbon and nitrogen through the soil. This module, together with the paramMod module contains
!relevant parameters for the flux and balance equations used in the model.
!Pool abbreviations:
!LITm - metabolic litter
!LITs - Structural litter
!SAPb - bacteria saprotrophs
!SAPf - fungi saprotrophs
!EcM  - Ectomycorrhiza
!AM   - Arbuscular mycorrhiza
!SOMp - Physically protected soil organic matter
!SOMa - Available soil organic matter
!SOMc - Chemically protected soil organic matter


module mycmimMod
  use shr_kind_mod,   only : r8 => shr_kind_r8
  use paramMod        !TODO: use only relevant stuff?
  use initMod,        only: nlevels,calc_init_NH4
  use readMod,        only: read_maxC, read_time,read_clm_model_input,read_clay,&
                            read_WATSAT_and_profiles,read_PFTs,calc_PFT,read_mort_profiles, &
                            read_clm_mortality
  use writeMod,       only: create_netcdf,create_yearly_mean_netcdf,fill_netcdf, &
                        fill_yearly_netcdf, fluxes_netcdf,store_parameters,check,fill_MMK
  use testMod,        only: respired_mass, test_mass_conservation_C,test_mass_conservation_N, &
                        total_carbon_conservation,total_nitrogen_conservation
  use dispmodule,     only: disp !External module to pretty print matrices (mainly for testing purposes)
  use netcdf,         only: nf90_nowrite,nf90_write,nf90_close,nf90_open
  
  implicit none
    private 
    public :: decomp

    !Fluxes etc:
    real(r8)              :: C_LITmSAPb, C_LITsSAPb, C_EcMSOMp, C_EcMSOMa, C_EcMSOMc, C_AMSOMp, &
    C_LITmSAPf, C_LITsSAPf, C_AMSOMa,  C_AMSOMc,  C_SOMaSAPb,C_SOMaSAPf, &
    C_SOMpSOMa, C_SOMcSOMa, C_SAPbSOMa,C_SAPbSOMp,C_SAPbSOMc,C_SAPfSOMa, &
    C_SAPfSOMp, C_SAPfSOMc, C_PlantSOMc,C_PlantSOMp, &
    N_LITmSAPb, N_LITsSAPb, N_EcMSOMp, N_EcMSOMa, N_EcMSOMc,  N_AMSOMp, &
    N_AMSOMa,N_AMSOMc, N_SOMaSAPb,N_SOMaSAPf, N_SOMpSOMa, N_SOMcSOMa, &
    N_LITmSAPf, N_LITsSAPf, N_INPlant, N_INEcM,&
    N_INAM, N_EcMPlant, N_AMPlant,N_SAPbSOMa, N_SAPbSOMp, N_SAPbSOMc, &
    N_SAPfSOMa, N_SAPfSOMp, N_SAPfSOMc,N_SOMcEcM,N_SOMpEcM, &
    C_PlantEcM,  C_PlantAM, C_PlantLITm, C_PlantLITs, C_EcMdecompSOMp, &
    C_EcMdecompSOMc,Leaching, Deposition,nitrif_rate, &
    N_demand_SAPf,N_demand_SAPb,N_INSAPb,N_INSAPf,N_PlantSOMp, &
    C_EcMenz_prod, N_PlantLITs, N_PlantLITm,N_PlantSOMc

    real(r8) :: NH4_sol_final, NH4_sorp_final,NO3_final

    !counts for how many times mineralization/immobilization occurs. Sum is printed at the end of decomp subroutine: 
    integer                             :: c1a
    integer                             :: c1b
    integer                             :: c2
    integer                             :: c3a
    integer                             :: c3b
    integer                             :: c4a
    integer                             :: c4b
    real(r8) :: save_N,save_C
  
    !For storing profiles read from CLM input file
    real(r8),dimension(:),allocatable             :: ndep_prof
    real(r8),dimension(:),allocatable             :: leaf_prof
    real(r8),dimension(:),allocatable             :: froot_prof
    real(r8),dimension(:),allocatable             :: stem_prof
    real(r8),dimension(:),allocatable             :: croot_prof
      
contains
  
  subroutine annual_mean(output_path,yearly_sumC,yearly_sumN,yearly_sumNinorg,sum_HR, year, run_name) !Calc. annual means and store them to file
    !Input 
    real(r8), dimension(nlevels,pool_types)   , intent(in):: yearly_sumC
    real(r8), dimension(nlevels,pool_types_N) , intent(in):: yearly_sumN
    real(r8), dimension(nlevels,inorg_N_pools), intent(in):: yearly_sumNinorg
    real(r8), intent(in)  :: sum_HR    !NB: Sum over all layers
    integer,  intent(in)  :: year
    CHARACTER (len = *), intent(in) :: run_name
    character (len = *), intent(in) :: output_path

    !Local
    real(r8), dimension(nlevels,pool_types)     :: yearly_meanC
    real(r8), dimension(nlevels,pool_types_N)   :: yearly_meanN
    real(r8), dimension(nlevels,inorg_N_pools)  :: yearly_meanNinorg
    real(r8)                                    :: yearly_mean_HR
    integer, parameter                          :: hr_in_year = 24*365
    
    yearly_meanC=yearly_sumC/hr_in_year
    yearly_meanN=yearly_sumN/hr_in_year
    yearly_meanNinorg=yearly_sumNinorg/hr_in_year
    
    yearly_mean_HR = sum_HR/hr_in_year
    
    call fill_yearly_netcdf(output_path,run_name, year, yearly_meanC,yearly_meanN,yearly_meanNinorg,yearly_mean_HR,delta_z(1:nlevels))
  end subroutine annual_mean

  subroutine decomp(nsteps,   &
                    run_name, &
                    write_hour,&
                    pool_C_start,pool_N_start,inorg_N_start, &
                    pool_C_final,pool_N_final,inorg_N_final,  &
                    start_year,stop_year,clm_input_path, &
                    clm_mortality_path, &
                    clm_surf_path,out_path) !Main subroutine. Calculates changes in each pool at each timestep (dC/dt & dN/dt)
    !IN:
    integer,intent(in)                        :: nsteps                               ! number of time steps to iterate over
    character (len=*) ,intent(in)             :: run_name                             ! used for naming output files
    integer,intent(in)                        :: write_hour                           ! How often output is written to file
    real(r8),intent(in)                       :: pool_C_start(nlevels,pool_types)     ! For store and output final C pool concentrations [gC/m3] 
    real(r8),intent(in)                       :: pool_N_start(nlevels,pool_types_N)   ! For storing and output final N pool concentrations [gN/m3] 
    real(r8),intent(in)                       :: inorg_N_start(nlevels,inorg_N_pools) ! For storing and output final N pool concentrations [gN/m3]      
    integer, intent(in)                       :: start_year                           ! Forcing start year
    integer, intent(in)                       :: stop_year                            ! Forcing end year, i.e. forcing loops over interval start_year-stop_year
    character (len=*) ,intent(in)             :: clm_input_path                       ! file path for input
    character (len=*) ,intent(in)             :: clm_mortality_path                       ! file path for input
    character (len=*) ,intent(in)             :: clm_surf_path                        ! file path for surface data
    character (len=*) ,intent(in)             :: out_path                             ! file path for output file

    !OUTPUT
    real(r8),intent(out)                      :: pool_C_final(nlevels,pool_types)     ! For store and output final C pool concentrations 
    real(r8),intent(out)                      :: pool_N_final(nlevels,pool_types_N)   ! For storing and output final N pool concentrations [gN/m3] 
    real(r8),intent(out)                      :: inorg_N_final(nlevels,inorg_N_pools) ! For storing and output final inorganic N pool concentrations [gN/m3] 
    
    
    !Shape of pool_matrixC/change_matrixC
    !|       LITm LITs SAPb SAPf EcM AM SOMp SOMa SOMc |
    !|level1   1   2    3    4   5   6   7   8    9    |
    !|level2                                           |
    !| .                                               |
    !| .                                               |
    !|nlevels _________________________________________|
    
    !Shape of the pool_matrixN/change_matrixN
    !|       LITm LITs SAPb SAPf EcM AM SOMp SOMa SOMc|
    !|level1   1   2    3    4   5   6   7   8    9   |
    !|level2                                          |
    !| .                                              |
    !| .                                              |
    !|nlevels ________________________________________|
    
    !Shape of inorg_N_matrix
    !|       NH4_sol NH4_sorp   NO3|
    !|level1    1       2       3  |
    !|level2                       |
    !| .                           |
    !| .                           |
    !|nlevels _____________________|
    
    !LOCAL
    real(r8)                        :: time                                 ! [hr] t*dt 
    integer,parameter               :: t_init=1                             ! initial time
    integer                         :: date
    character (len=4)               :: year_fmt                             ! For printing/writing
    character (len=4)               :: year_char                            ! For printing/writing.
    logical                         :: isVertical                           ! True if we use vertical soil layers. NOTE: Code need updates if isVertical=True should work
    logical                         :: Spinup_run                           ! Based on starting year (not ideal..) 
    integer                         :: input_steps                          ! Number of time entries in CLM infile
    integer                         :: j,i,t                                !for iterations

    !Counters
    integer                         :: ycounter, year,write_y
    integer                         :: write_counter            !used for determining when to output results 
    integer                         :: month_counter      !for determining when to read new input values (Temp. litfall etc.)
    integer                         :: day_counter
    integer                         :: spinup_counter
    integer                         :: current_month
    integer                         :: current_day
                
    !NetCDF identifier variables
    integer :: ncid
    integer :: writencid
    integer :: spinupncid
    integer :: mort_ncid
    integer :: mort_spinupncid
    
    logical :: result_file_exist
    !Pool and respiration related variables:
    real(r8)                        :: C_Loss, C_Gain, N_Gain, N_Loss
    real(r8)                        :: pool_matrixC(nlevels,pool_types)     ! For storing C pool sizes [gC/m3]
    real(r8)                        :: change_matrixC(nlevels,pool_types)   ! For storing dC/dt for each time step [gC/(m3*hour)]
    real(r8)                        :: pool_temporaryC(nlevels,pool_types)  ! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
    real(r8)                        :: pool_matrixC_previous(nlevels,pool_types) !Used for checking mass conservation
    
    real(r8)                        :: pool_matrixN_previous(nlevels,pool_types_N) !Used for checking mass conservation
    real(r8)                        :: inorg_N_matrix_previous(nlevels,inorg_N_pools) !Used for checking mass conservation
    real(r8)                        :: pool_temporaryN(nlevels,pool_types_N)! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
    real(r8)                        :: pool_matrixN(nlevels,pool_types_N)   ! For storing N pool sizes [gN/m3] parallell to C pools and  inorganic N
    real(r8)                        :: change_matrixN(nlevels,pool_types_N) ! For storing dC/dt for each time step [gN/(m3*hour)]
    real(r8)                        :: inorg_N_matrix(nlevels, inorg_N_pools) ! Inorganic N pool sizes (gN/m3)
    real(r8)                        :: inorg_Ntemporary_matrix(nlevels, inorg_N_pools) ! Inorganic N pool sizes (gN/m3)
    
    real(r8),dimension(nlevels)     :: HR                                   ! For storing the C  that is lost to respiration [gC/m3h]
    real(r8),dimension(nlevels)     :: HRb,HRf                              !For storing respiration separately for bacteria and fungi
    real(r8),dimension(nlevels)     :: HRe,HRa                              !For storing respiration separately for mycorrhiza
    real(r8)                        :: HR_mass_accumulated, HR_mass
    real(r8)                        :: HR_mass_yearly

    !For calculating means
    real(r8)                       :: sum_consN(nlevels, pool_types_N)                !gN/m3 for calculating annual mean
    real(r8)                       :: sum_consNinorg(nlevels, inorg_N_pools)          !gN/m3 for calculating annual mean
    real(r8)                       :: sum_consC(nlevels, pool_types)                  !gC/m3 for calculating annual mean
    real(r8)                       :: monthly_sum_consN(nlevels, pool_types_N)        !gN/m3 for calculating monthly mean
    real(r8)                       :: monthly_sum_consNinorg(nlevels, inorg_N_pools)  !gN/m3 for calculating monthly mean
    real(r8)                       :: monthly_sum_consC(nlevels, pool_types)          !gC/m3 for calculating monthly mean
    
    !For mass conservation checks
    real(r8)                       :: sum_input_step          ! C, Used for checking mass conservation
    real(r8)                       :: sum_input_total         ! C, Used for checking mass conservation
    real(r8)                       :: sum_N_input_total       ! N
    real(r8)                       :: sum_N_out_total         ! N
    real(r8)                       :: sum_N_out_step          ! N
    real(r8)                       :: sum_N_input_step        ! Used for checking mass conservation
    real(r8)                       :: pool_C_start_for_mass_cons(nlevels,pool_types)           ! To use in total mass conservation subroutine
    real(r8)                       :: pool_N_start_for_mass_cons(nlevels,pool_types_N)         ! To use in total mass conservation subroutine
    real(r8)                       :: pool_Ninorg_start_for_mass_cons(nlevels,inorg_N_pools)   ! To use in total mass conservation subroutine

    !Related to vertical transport
    real(r8),allocatable           :: vertC(:,:)         !Stores the vertical change in a time step, same shape as change_matrixC
    real(r8),allocatable           :: vertN(:,:)         !Stores the vertical change in a time step, same shape as change_matrixN
    real(r8),allocatable           :: vert_inorgN(:,:)         !Stores the vertical change in a time step, same shape as change_inorgN
    real(r8)                       :: soil_depth       ![m] 

    !For storing input related variables
    real(r8)                       :: N_DEPinput
    real(r8)                       :: C_MYCinput
    real(r8)                       :: C_leaf_litter
    real(r8)                       :: C_root_litter
    real(r8)                       :: N_leaf_litter
    real(r8)                       :: N_root_litter
    real(r8)                       :: C_CWD_litter(nlevels)
    real(r8)                       :: N_CWD_litter(nlevels)
    real(r8)                       :: met_mortC(nlevels)
    real(r8)                       :: met_mortN(nlevels)
    real(r8)                       :: split_mortC(nlevels)
    real(r8)                       :: split_mortN(nlevels)

    !For reading soil temperature and moisture from CLM output file
    real(r8), dimension(nlevels)    :: TSOIL
    real(r8), dimension(nlevels)    :: SOILLIQ
    real(r8), dimension(nlevels)    :: SOILICE
    real(r8), dimension(nlevels)    :: WATSAT
    real(r8), dimension(nlevels)    :: W_SCALAR
    real(r8), dimension(nlevels)    :: T_SCALAR
    real(r8)                        :: qdrain
    real(r8)                        :: qrunoff
    real(r8), dimension(nlevels)    :: h2osoi_liq
    real(r8)                        :: H2OSOI
    real(r8),dimension(15)          :: PFT_distribution
    
    !Used in different functions/subroutines related to fluxes
    real(r8)                          :: desorp           ![1/h]From Mimics, used for the transport from physically protected SOM to available SOM pool
    real(r8),dimension(:),allocatable :: norm_froot_prof             !Normalized root profile used for turnover depth dependence.
    real(r8)                          :: fMET                        ![-] Fraction determining distribution of total litter production between LITm and LITs
    real(r8)                          :: lflitcn_avg                 !Used in function calc_met_fraction
    real(r8), dimension(nlevels)      :: ROI_EcM
    real(r8), dimension(nlevels)      :: ROI_AM
    real(r8), dimension(nlevels,no_of_myc_pools)                :: f_alloc  !fraction of incoming C from plant that is allocated to EcM(1) and AM(2)
    real(r8), dimension(no_of_som_pools,no_of_sap_pools)        :: f_sapsom !Fraction of dead saprotrophic biomass to the different SOM pools.

      call system_clock(count_rate=clock_rate) !Find the time rate
      call system_clock(count=clock_start)     !Start Timer  
      
      if (nlevels>1) then 
        soil_depth=sum(delta_z(1:nlevels))
        isVertical = .True.
      else
        soil_depth=1.52_r8
        isVertical = .False.
        allocate (vertC, mold = pool_matrixC)
        allocate (vertN, mold = pool_matrixN)
        allocate (vert_inorgN, mold = inorg_N_matrix)
                    !TODO: This can be done better, and probably needs modification to work 
      end if
          
      !Allocate and initialize
      allocate(CUE_bacteria_vr(nlevels))
      CUE_bacteria_vr=CUE_0b
      allocate(CUE_fungi_vr(nlevels))
      CUE_fungi_vr=CUE_0f
      allocate(CUE_ecm_vr(nlevels))
      CUE_ecm_vr=CUE_myc_0
      allocate(CUE_am_vr(nlevels))
      CUE_am_vr=CUE_myc_0      
      allocate(r_moist(nlevels))
      allocate(f_enzprod(nlevels))
      f_enzprod=f_enzprod_0
      allocate(NH4_sorp_eq_vr(nlevels))
      NH4_sorp_eq_vr=10.0 

      r_myc=1.0 !Initialize EcM input modifier

      !For counting mineralization/immobilization occurences:
      c1a=0;c1b=0;c2=0;c3a=0;c3b=0;c4a=0;c4b=0
      save_N=0._r8;save_C=0._r8

      !Set initial concentration values:
      pool_matrixC              = pool_C_start
      pool_matrixN              = pool_N_start
      inorg_N_matrix            = inorg_N_start
      pool_matrixC_previous     = pool_C_start
      pool_matrixN_previous     = pool_N_start
      inorg_N_matrix_previous   = inorg_N_start
      pool_C_start_for_mass_cons= pool_C_start
      pool_N_start_for_mass_cons= pool_N_start
      pool_Ninorg_start_for_mass_cons=inorg_N_start
      pool_temporaryC           = pool_C_start
      pool_temporaryN           = pool_N_start
      inorg_Ntemporary_matrix   = 0.0

      !Make sure things start from zero
      sum_consN = 0;sum_consNinorg = 0; sum_consC = 0
      monthly_sum_consN = 0; monthly_sum_consNinorg = 0; monthly_sum_consC = 0
      change_matrixC = 0.0;change_matrixN = 0.0
      HR = 0.0;  HRb = 0.0; HRf = 0.0; HRe = 0.0; HRa = 0.0
      ROI_EcM        = 0.0 ; ROI_AM      = 0.0
      f_alloc        = 0.0
      
      HR_mass_accumulated = 0
      sum_input_step   =0.0
      sum_input_total  =0.0
      sum_N_input_step =0.0
      sum_N_input_total=0.0
      sum_N_out_step   =0.0
      sum_N_out_total  =0.0

      !For counting:
      year     = start_year
      year_fmt = '(I4)'
      write (year_char,year_fmt) year
      write_y=0
      current_month = 1
      current_day   = 1      
      month_counter = 0
      day_counter   = 0
      write_counter       = 0
      ycounter      = 0

      allocate(ndep_prof(nlevels),leaf_prof(nlevels),froot_prof(nlevels), norm_froot_prof(nlevels), croot_prof(nlevels),stem_prof(nlevels))

      call read_WATSAT_and_profiles(adjustr(clm_input_path)//'all.'//"1901.nc",WATSAT,ndep_prof,froot_prof,leaf_prof) !NOTE: This subroutine needs to read a file that contains WATSAT etc. 

      call read_mort_profiles(adjustr(clm_mortality_path)//"1901.nc",croot_prof,stem_prof)

 
      !read data from CLM file
      if ( start_year == 1850 ) then
        Spinup_run = .True.
        spinup_counter =1
        call check(nf90_open(trim(adjustr(clm_input_path)//'for_spinup.1850-1869.nc'), nf90_nowrite, spinupncid))  
        call check(nf90_open(trim(adjustr(clm_mortality_path)//'for_spinup.1850-1869.nc'),nf90_nowrite,mort_spinupncid))
        call read_time(spinupncid,input_steps) !Check if inputdata is daily or monthly ("steps" is output)         
        call read_clm_model_input(spinupncid,spinup_counter,CLM_version, &
        N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
        C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
        W_SCALAR,T_SCALAR,qdrain,qrunoff,h2osoi_liq,C_CWD_litter,N_CWD_litter)    
        call read_clm_mortality(mort_spinupncid,spinup_counter,froot_prof,croot_prof,leaf_prof,stem_prof,met_mortC,met_mortN,split_mortC,split_mortN)     
      else
        Spinup_run = .False. 
        call check(nf90_open(trim(adjustr(clm_input_path)//'all.'//year_char//'.nc'), nf90_nowrite, ncid)) 
        call check(nf90_open(trim(adjustr(clm_mortality_path)//year_char//'.nc'),nf90_nowrite,mort_ncid))
        call read_time(ncid,input_steps) !Check if inputdata is daily or monthly ("steps" is output) 
        call read_clm_model_input(ncid,1, CLM_version,&
        N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
        C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
        W_SCALAR,T_SCALAR,qdrain,qrunoff,h2osoi_liq,C_CWD_litter,N_CWD_litter)

        call read_clm_mortality(mort_ncid,1,froot_prof,croot_prof,leaf_prof,stem_prof,met_mortC,met_mortN,split_mortC,split_mortN)     

      end if

      call moisture_func(SOILLIQ,WATSAT,SOILICE,r_moist)                   

      call read_clay(adjustr(clm_surf_path),fCLAY)
      call calc_PFT(adjustr(clm_input_path)//'all.'//"1901.nc",lflitcn_avg)
      
      if ( .not. use_ROI ) then !use static PFT determined fractionation between EcM and AM C input
        call read_PFTs(adjustr(clm_surf_path),PFT_distribution)
        f_alloc(:,1) = calc_EcMfrac(PFT_distribution)
        f_alloc(:,2) = 1-calc_EcMfrac(PFT_distribution)
      end if
      
      if (f_alloc(1,1)==1.0 ) then !To avoid writing errors when there is no AM (alloc to AM = 0 and alloc to EcM = 1)
        pool_matrixC(:,6)=0.0
        pool_matrixN(:,6)=0.0
        pool_matrixC_previous(:,6)=0.0
        pool_matrixN_previous(:,6)=0.0
        pool_C_start_for_mass_cons=pool_matrixC
        pool_N_start_for_mass_cons=pool_matrixN
      end if
      
      if ( Spinup_run ) then
        call create_yearly_mean_netcdf(trim(out_path),run_name)  !open and prepare files to store results. Store initial values    
        max_myc_alloc = read_maxC(spinupncid,input_steps) 
      else
        max_myc_alloc = read_maxC(ncid,input_steps)        
      end if
      
      desorp          = calc_desorp(fCLAY)
      fMET            = calc_met_fraction(C_leaf_litter,C_root_litter,C_CWD_litter,lflitcn_avg)
      norm_froot_prof = (froot_prof-minval(froot_prof))/(maxval(froot_prof)-minval(froot_prof))
      
      !Create file and write initial vallues to it:
      call create_netcdf(trim(out_path),run_name,start_year)
      call check(nf90_open(trim(out_path)//trim(run_name)//".nc", nf90_write, writencid))  
      call fill_netcdf(writencid,t_init, pool_matrixC, pool_matrixN,inorg_N_matrix, &
      date, HR_mass_accumulated,HR,HRb,HRf,HRe,HRa, change_matrixC,&
      change_matrixN,inorg_Ntemporary_matrix, write_hour,current_month,TSOIL, r_moist, &
      CUE_bacteria_vr,CUE_fungi_vr,CUE_EcM_vr,CUE_am_vr,ROI_EcM=ROI_EcM, & 
      ROI_AM=ROI_AM,enz_frac=f_enzprod,f_met = fMET,f_alloc=f_alloc, NH4_eq=inorg_N_matrix(:,2))
      
      !print initial values to terminal      
      call disp("InitC", pool_matrixC)
      call disp("InitN", pool_matrixN)
      call disp("InitN inorganic", inorg_N_matrix)
      !----------------------------------------------------------------------------------------------------------------
      do t =1,nsteps !Starting time iterations
        
        time = t*dt
        !Time counters:
        write_counter  = write_counter  + 1
        ycounter = ycounter + 1           !Counts hours in a year 
        month_counter = month_counter + 1 !Counts hours in a month
        day_counter   = day_counter   + 1     !Counts hours in a day

        !------------Update forcing and environmental data at the right timesteps-----------
        if (month_counter == days_in_month(current_month)*hr_pr_day/dt+1) then
          month_counter = 1       
          
          if (current_month == 12) then
            current_month=1 !Update to new year
          else
            current_month = current_month + 1 
          end if   

          if (input_steps==12) then  !Input is given as monthly values
            call read_clm_model_input(ncid,current_month,CLM_version, &
            N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
            C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
            W_SCALAR,T_SCALAR,qdrain,qrunoff,h2osoi_liq,C_CWD_litter,N_CWD_litter)  
            call read_clm_mortality(mort_ncid,current_month,froot_prof,croot_prof,leaf_prof,stem_prof,met_mortC,met_mortN,split_mortC,split_mortN) 
            call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)  
            max_myc_alloc = read_maxC(ncid,input_steps)
            fMET = calc_met_fraction(C_leaf_litter,C_root_litter,C_CWD_litter,lflitcn_avg)
          end if     
          
          if (input_steps==240) then !Input is monthly values in the "spinup forcing file"
            spinup_counter = spinup_counter+1
            call read_clm_model_input(spinupncid,spinup_counter, CLM_version,&
            N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
            C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
            W_SCALAR,T_SCALAR,qdrain,qrunoff,h2osoi_liq,C_CWD_litter,N_CWD_litter)  
            call read_clm_mortality(mort_spinupncid,spinup_counter,froot_prof,croot_prof,leaf_prof,stem_prof,met_mortC,met_mortN,split_mortC,split_mortN)
            call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)   
            max_myc_alloc = read_maxC(spinupncid,input_steps)           
            fMET = calc_met_fraction(C_leaf_litter,C_root_litter,C_CWD_litter,lflitcn_avg)
          end if  

        end if 
          
        if (day_counter == hr_pr_day/dt+1) then
          day_counter = 1   
          if ( current_day == 365 ) then
            current_day = 1
          else 
            current_day = current_day +1                      
          end if
          
          if (input_steps==365) then  !Input is given as daily values
            call read_clm_model_input(ncid,current_day,CLM_version, &
            N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
            C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
            W_SCALAR,T_SCALAR,qdrain,qrunoff,h2osoi_liq,C_CWD_litter,N_CWD_litter)
            call read_clm_mortality(mort_ncid,current_day,froot_prof,croot_prof,leaf_prof,stem_prof,met_mortC,met_mortN,split_mortC,split_mortN)
            call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)        
            max_myc_alloc = read_maxC(ncid,input_steps)                
            fMET = calc_met_fraction(C_leaf_litter,C_root_litter,C_CWD_litter,lflitcn_avg)
          end if        
        end if         
        !-----------------------------------------------------------------------------------

        ! if ( year == 2012 .and. ycounter == 181*24 ) then !IF test for litter bag experiments. Add litter at certain date.
        !   pool_matrixC(3,1) = pool_matrixC(3,1) + 10/delta_z(3) !Add 10gC to layer 3 
        !   pool_matrixN(3,1) = pool_matrixN(3,1) + 10*(pool_matrixN(3,1)/pool_matrixC(3,1))/delta_z(3) !Add N to layer
        !   print*, "Added litter; ", 10/delta_z(3), "gC/m3, and ",  10*(pool_matrixN(3,1)/pool_matrixC(3,1))/delta_z(3), "gN/m3 to LITm, layer 3."
        ! end if

        r_myc = myc_modifier(C_MYCinput,max_myc_alloc) !calculate factor that scales mycorrhizal activity based on C "payment" from plant
        if ( abs(r_myc) > 1.0 ) then
          print*, "r_myc, :", r_myc,"time: ", time !for checking
        end if
        
        ! Fracions of SAP that goes to different SOM pools
        f_sapsom  = calc_sap_to_som_fractions(fCLAY,fMET)
        fPHYS     = f_sapsom(1,:)
        fCHEM     = f_sapsom(2,:)
        fAVAIL    = f_sapsom(3,:)
                
        do j = 1, nlevels !For each depth level:
          
          H2OSOI=SOILLIQ(j)+SOILICE(j) !Used for N sorp/desorp calculations
          
          !Michaelis Menten parameters:
          Km      = reverse_Km_function(TSOIL(j),fCLAY)
          Vmax    = Vmax_function(TSOIL(j),r_moist(j)) !  ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics.

          ![1/h] Microbial turnover rate (SAP to SOM)
          k_sapsom  = calc_sap_turnover_rate(fMET, TSOIL(j), norm_froot_prof(j)) 
          k_mycsom  = calc_myc_mortality()  
          
          !Initialize CUEs for timestep
          CUE_bacteria_vr(j)= (CUE_slope*TSOIL(j)+CUE_0b)
          CUE_fungi_vr(j)   = (CUE_slope*TSOIL(j)+CUE_0f)
          CUE_ecm_vr(j)     = CUE_myc_0
          CUE_am_vr(j)      = CUE_myc_0
          
          !From CLM data, calculate litter input rates: 
          call input_rates(j,fMET,C_leaf_litter,C_root_litter,N_leaf_litter,&
                                      N_root_litter,N_CWD_litter,C_CWD_litter,&
                                      C_PlantLITm,C_PlantLITs, &
                                      N_PlantLITm,N_PlantLITs, &
                                      C_PlantSOMp,C_PlantSOMc, &
                                      N_PlantSOMp,N_PlantSOMc, &
                                      met_mortC, met_mortN, &
                                      split_mortC, split_mortN)

          !Calculate inorganic N rates: 
          Leaching    = calc_Leaching_Runoff(qdrain,qrunoff,h2osoi_liq,inorg_N_matrix(j,3),j) !N31 !TODO: Call this Leaching_runoff or something, as it also contains RUNOFF
          Deposition  = set_N_dep(CLMdep = N_DEPinput*ndep_prof(j)) !NOTE: either const_dep = some_value or CLMdep = N_DEPinput*ndep_prof(j) !N32

          nitrif_rate = calc_nitrification((inorg_N_matrix(j,1)+Deposition*dt),W_SCALAR(j),T_SCALAR(j),TSOIL(j)) !NOTE: Uses NH4 + Deposiiton from current timestep !N34
       
          !Calculate fluxes between pools in level j at timestep:
          call calculate_fluxes(j,TSOIL(j),H2OSOI, pool_matrixC, pool_matrixN,inorg_N_matrix,Deposition, Leaching, nitrif_rate,soil_depth,desorp) 

          !Related to NH4 sol-sorp. TODO: Can be made more intuitive
          inorg_Ntemporary_matrix(j,1)=NH4_sol_final
          inorg_Ntemporary_matrix(j,2)=NH4_sorp_final          
          inorg_Ntemporary_matrix(j,3)=NO3_final
          
          if ( use_ROI ) then !Calculate size of fraction of C payment from plant that goes to EcM and AM.
            ROI_EcM(j) = ROI_function(N_INEcM+N_SOMpEcM+N_SOMcEcM,pool_matrixC(j,5),k_mycsom(1))
            ROI_AM(j)  = ROI_function(N_INAM, pool_matrixC(j,6),k_mycsom(2))
            if ( C_MYCinput .NE. 0.0  ) then !To avoid NaN when both ROI is zero
              if ( ROI_EcM(j) + ROI_AM(j) == 0.0 ) then !Too little myc in layer to contribute
                f_alloc(j,:)=0.0 
              else
                f_alloc(j,1) = ROI_EcM(j)/(ROI_EcM(j)+ROI_AM(j)) !Eq. (4), Sulman et al. 2019 (DOI: 10.1029/2018GB005973)
                f_alloc(j,2) = ROI_AM(j)/(ROI_EcM(j)+ROI_AM(j))
              end if
            else
              f_alloc(j,:)=0.5 !Value does not really matter bc. C_MYCinput is zero
            end if
          end if    

          C_PlantEcM = f_alloc(j, 1)*C_MYCinput*froot_prof(j) !C28
          C_PlantAM  = f_alloc(j, 2)*C_MYCinput*froot_prof(j) !C29    
               
          
          call myc_to_plant(CUE_EcM_vr(j),CUE_AM_vr(j),f_enzprod(j),N_AMPlant,N_EcMPlant) !N29,N30
          
          C_EcMenz_prod=CUE_ecm_vr(j)*C_PlantEcM*f_enzprod(j) !C28

          if (write_counter == write_hour/dt .or. t==1) then !Write fluxes from calculate_fluxes to file            
            call fluxes_netcdf(writencid,int(time), write_hour, j, &
                            C_LITmSAPb, C_LITsSAPb, C_EcMSOMp, C_EcMSOMa, C_EcMSOMc, C_AMSOMp, &
                            C_LITmSAPf, C_LITsSAPf, C_AMSOMa,  C_AMSOMc,  C_SOMaSAPb,C_SOMaSAPf, &
                            C_SOMpSOMa, C_SOMcSOMa, C_SAPbSOMa,C_SAPbSOMp,C_SAPbSOMc,C_SAPfSOMa, &
                            C_SAPfSOMp, C_SAPfSOMc, C_PlantSOMc,C_PlantSOMp, &
                            N_LITmSAPb, N_LITsSAPb, N_EcMSOMp, N_EcMSOMa, N_EcMSOMc,  N_AMSOMp, &
                            N_AMSOMa,N_AMSOMc, N_SOMaSAPb,N_SOMaSAPf, N_SOMpSOMa, N_SOMcSOMa, &
                            N_LITmSAPf, N_LITsSAPf, N_INPlant, N_INEcM,&
                            N_INAM, N_EcMPlant, N_AMPlant,N_SAPbSOMa, N_SAPbSOMp, N_SAPbSOMc, &
                            N_SAPfSOMa, N_SAPfSOMp, N_SAPfSOMc,N_SOMcEcM,N_SOMpEcM, &
                            C_PlantEcM,  C_PlantAM, C_PlantLITm, C_PlantLITs, C_EcMdecompSOMp, &
                            C_EcMdecompSOMc,Leaching, Deposition,nitrif_rate, &
                            N_INSAPb,N_INSAPf,N_PlantSOMp,C_EcMenz_prod, N_PlantLITs, N_PlantLITm,N_PlantSOMc)
          end if !write fluxes

          do i = 1,pool_types !loop over all the pool types, i, in depth level j . !NOTE: This does not really need to be in a do-loop..consider using fortran types.
            !This if-loop calculates dC/dt and dN/dt for the different carbon pools.

            if (i==1) then !LITm
              C_Gain = C_PlantLITm
              C_Loss = C_LITmSAPb + C_LITmSAPf
              N_Gain = N_PlantLITm
              N_Loss = N_LITmSAPb + N_LITmSAPf
              
            elseif (i==2) then !LITs
              C_Gain = C_PlantLITs
              C_Loss = C_LITsSAPb + C_LITsSAPf
              N_Gain = N_PlantLITs
              N_Loss = N_LITsSAPb + N_LITsSAPf

            elseif (i==3) then !SAPb 
              C_Gain = CUE_bacteria_vr(j)*(C_LITmSAPb + C_LITsSAPb &
                                          + C_SOMaSAPb)
              C_Loss =  C_SAPbSOMp + C_SAPbSOMa + C_SAPbSOMc
              N_Gain = (N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb)*NUE
              N_Loss = N_SAPbSOMp + N_SAPbSOMa + N_SAPbSOMc
              if ( N_INSAPb>0 ) then
                N_Gain = N_Gain + N_INSAPb              
              else
                N_Loss=N_Loss-N_INSAPb     !two minus becomes +         
              end if

            elseif (i==4) then !SAPf
              C_Gain = CUE_fungi_vr(j)*(C_LITmSAPf + C_LITsSAPf &
                                       + C_SOMaSAPf)
              C_Loss =  C_SAPfSOMp + C_SAPfSOMa + C_SAPfSOMc
              N_Gain = (N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf)*NUE
              N_Loss = N_SAPfSOMp + N_SAPfSOMa + N_SAPfSOMc
              if ( N_INSAPf>0 ) then
                N_Gain = N_Gain + N_INSAPf 
              else
                N_Loss=N_Loss-N_INSAPf !two minus becomes +
              end if

            elseif (i==5) then !EcM
              
              C_Gain = CUE_ecm_vr(j)*C_PlantEcM
              C_Loss = C_EcMSOMp + C_EcMSOMa + C_EcMSOMc + C_EcMenz_prod
              N_Gain = N_INEcM + N_SOMpEcM + N_SOMcEcM
              N_Loss = N_EcMPlant + N_EcMSOMa + N_EcMSOMp + N_EcMSOMc

           elseif (i==6) then !AM
              C_Gain = CUE_am_vr(j)*C_PlantAM
              C_Loss = C_AMSOMp + C_AMSOMa + C_AMSOMc
              N_Gain = N_INAM 
              N_Loss = N_AMPlant + N_AMSOMa + N_AMSOMp + N_AMSOMc

            elseif (i==7) then !SOMp
              C_Gain = C_SAPbSOMp + C_SAPfSOMp + C_EcMSOMp + C_AMSOMp+ C_PlantSOMp
              C_Loss = C_SOMpSOMa+C_EcMdecompSOMp 
              N_Gain = N_SAPbSOMp + N_SAPfSOMp + N_EcMSOMp + N_AMSOMp+N_PlantSOMp
              N_Loss = N_SOMpSOMa + N_SOMpEcM

            elseif (i==8) then !SOMa
              C_Gain = C_SAPbSOMa + C_SAPfSOMa + C_EcMSOMa + C_EcMdecompSOMp + C_EcMdecompSOMc + &
              C_AMSOMa + C_SOMpSOMa + C_SOMcSOMa + C_EcMenz_prod
              C_Loss = C_SOMaSAPb + C_SOMaSAPf 
              N_Gain = N_SAPbSOMa + N_SAPfSOMa + N_EcMSOMa + &
              N_AMSOMa + N_SOMpSOMa + N_SOMcSOMa 
              N_Loss = N_SOMaSAPb + N_SOMaSAPf

            elseif (i==9) then !SOMc
              C_Gain =  C_SAPbSOMc + C_SAPfSOMc + C_EcMSOMc + C_AMSOMc+C_PlantSOMc
              C_Loss = C_SOMcSOMa+C_EcMdecompSOMc
              N_Gain =  N_SAPbSOMc + N_SAPfSOMc + N_EcMSOMc + N_AMSOMc+N_PlantSOMc
              N_Loss = N_SOMcSOMa + N_SOMcEcM

            else
              print*, 'Too many pool types expected, pool_types = ',pool_types, 'i: ', i
              stop
            end if !determine total gains and losses

            change_matrixC(j,i) = C_Gain - C_Loss !net change in timestep
            change_matrixN(j,i) = N_Gain - N_loss   

            !Store these values as temporary so that they can be used in the vertical diffusion subroutine
            pool_temporaryC(j,i)=pool_matrixC(j,i) + change_matrixC(j,i)*dt
            pool_temporaryN(j,i)=pool_matrixN(j,i) + change_matrixN(j,i)*dt
            
            if ( pool_temporaryC(j,i) < trunc_value ) then !If value is smaller than a set truncation value, round to zero and keep track of the mass discared 
              !print*, pool_temporaryC(j,i),change_matrixC(j,i)*dt, j, i, time
              !call disp(pool_temporaryC)
              !stop
              save_C=save_C+pool_temporaryC(j,i)
              pool_temporaryC(j,i)=0._r8
            end if  

            if ( pool_temporaryN(j,i) < trunc_value ) then  
              save_N=save_N + pool_temporaryN(j,i)*dt*delta_z(j)
              pool_temporaryN(j,i)=0._r8                
            end if                        

          end do !i, pool_types
          
          !Calculate the heterotrophic respiration loss from depth level j in timestep t:
          
          HRb(j) = ( C_LITmSAPb + C_LITsSAPb + C_SOMaSAPb)*(1-CUE_bacteria_vr(j))*dt
          HRf(j) = ( C_LITmSAPf + C_LITsSAPf + C_SOMaSAPf)*(1-CUE_fungi_vr(j))*dt
          HRe(j) = C_PlantEcM*(1-CUE_ecm_vr(j))*dt 
          HRa(j) = C_PlantAM*(1-CUE_am_vr(j))*dt
          
          HR(j) = HRb(j) + HRf(j) + HRe(j) + HRa(j) 

          if (HR(j) < 0 ) then
            print*, 'Negative HR: ', HR(j), t,j
            print*, "Pools C", pool_matrixC(j,1),pool_matrixC(j,2),pool_matrixC(j,3),pool_matrixC(j,4), pool_matrixC(j,8)
            print*, "pools N", inorg_N_matrix(j,1), inorg_N_matrix(j,3), N_INSAPb, N_INSAPf
            print*, "Fluxes", C_LITmSAPb,C_LITsSAPb,C_SOMaSAPb,C_LITmSAPf,C_LITsSAPf,C_SOMaSAPf
            stop !for checking
          end if
          
          if (write_counter == write_hour/dt) then          
            call fill_MMK(writencid, int(time),write_hour,j,Km,Vmax)
          end if

          !Summarize in and out print timestep to check mass balance
          sum_input_step  = sum_input_step  +(C_PlantLITm+C_PlantLITs+C_PlantEcM+C_PlantAM+C_PlantSOMc+C_PlantSOMp)*dt*delta_z(j) !g/m2  
          sum_N_input_step= sum_N_input_step+(N_PlantLITm+N_PlantLITs+N_PlantSOMc+N_PlantSOMp+Deposition)*dt*delta_z(j) !g/m2
          sum_N_out_step  = sum_N_out_step  +(N_EcMPlant+N_AMPlant+N_INPlant+Leaching)*dt*delta_z(j)

        end do !j, depth_level
        
        !Store accumulated HR mass
        call respired_mass(HR, HR_mass)
        HR_mass_accumulated = HR_mass_accumulated + HR_mass
        HR_mass_yearly = HR_mass_yearly + HR_mass

        !TODO: tot_diffC, upperC, lowerC is not used and can be removed!
        if (isVertical) then
          call vertical_diffusion(pool_temporaryC,vertC,D_carbon)
          call vertical_diffusion(pool_temporaryN,vertN,D_nitrogen)
          call vertical_diffusion(inorg_Ntemporary_matrix,vert_inorgN,D_nitrogen)

          pool_matrixC    = vertC*dt        + pool_temporaryC
          pool_matrixN    = vertN*dt        + pool_temporaryN
          inorg_N_matrix  = vert_inorgN*dt  + inorg_Ntemporary_matrix   
        else
          pool_matrixC    = pool_temporaryC
          pool_matrixN    = pool_temporaryN
          inorg_N_matrix  = inorg_Ntemporary_matrix
        end if!isVertical

        sum_consN       = sum_consN       + pool_matrixN
        sum_consNinorg  = sum_consNinorg  + inorg_N_matrix
        sum_consC       = sum_consC       + pool_matrixC

        monthly_sum_consN       = monthly_sum_consN       + pool_matrixN
        monthly_sum_consNinorg  = monthly_sum_consNinorg  + inorg_N_matrix
        monthly_sum_consC       = monthly_sum_consC       + pool_matrixC
        

        if (ycounter == 365*24/dt) then !one year has passed
          ycounter = 0
          write_y =write_y+1 !For writing to annual mean file
          
          if ( Spinup_run ) then 
            call annual_mean(trim(out_path),sum_consC,sum_consN,sum_consNinorg,HR_mass_yearly,write_y , run_name) !calculates the annual mean and write the result to file
          end if
          if (year == stop_year) then !Start new cycle of forcing years
            year = start_year         
            spinup_counter=0            
          else 
            year = year + 1             
          end if
          write (year_char,year_fmt) year   

          if ( .not. Spinup_run ) then                              
            call check(nf90_close(ncid)) !Close netcdf file containing values for the past year
            call check(nf90_open(trim(adjustr(clm_input_path)//'all.'//year_char//'.nc'), nf90_nowrite, ncid)) !open netcdf containing values for the next year
            call read_time(ncid,input_steps)     
          end if           
          sum_consN =0
          sum_consC =0
          sum_consNinorg=0
          HR_mass_yearly=0
        end if

        if (write_counter == write_hour/dt) then
          write_counter = 0        
          call fill_netcdf(writencid, int(time), pool_matrixC, pool_matrixN,inorg_N_matrix,&
                          date, HR_mass_accumulated,HR,HRb,HRf,HRe,HRa,vertC,vertN,vert_inorgN, write_hour,current_month, &
                          TSOIL, r_moist,CUE_bacteria_vr,CUE_fungi_vr,CUE_ecm_vr,CUE_am_vr,ROI_EcM=ROI_EcM,&
                          ROI_AM=ROI_AM,enz_frac=f_enzprod,f_met = fMET,f_alloc=f_alloc, NH4_eq=NH4_sorp_eq_vr)
        end if!writing

        if (t == nsteps) then 
          pool_C_final  = pool_matrixC
          pool_N_final  = pool_matrixN    
          inorg_N_final = inorg_N_matrix              
          if (Spinup_run) then !Write to file
            inquire(file="Spinup_values/"//run_name//".dat", exist=result_file_exist)
            if (result_file_exist) then
              open(12, file="Spinup_values/"//run_name//".dat", status="old", action="write")
            else
              open(12, file="Spinup_values/"//run_name//".dat", status="new", action="write")
            end if
            write(12, *) "Simulation finished at year",year_char, " Spinup status: ", Spinup_run
            write(12, *) "pool_C_final: "
            write(12, *) pool_C_final
            write(12, *) "pool_N_final: "
            write(12, *) pool_N_final
            write(12, *) "inorg_N_final: "
            write(12, *) inorg_N_final
            close(12)
          end if
          call store_parameters(writencid,soil_depth,desorp,delta_z(1:nlevels))    
        end if
        
        !Mass conservation test for timestep:
        call test_mass_conservation_C(sum_input_step,HR_mass, &
                                      pool_matrixC_previous,pool_matrixC, &
                                      nlevels)
        call test_mass_conservation_N(sum_N_input_step,sum_N_out_step, &
                                      pool_matrixN_previous,inorg_N_matrix_previous,&
                                      pool_matrixN,inorg_N_matrix,nlevels)
                                      
        pool_matrixC_previous  = pool_matrixC
        pool_matrixN_previous  = pool_matrixN
        inorg_N_matrix_previous= inorg_N_matrix
        sum_input_total        = sum_input_total+sum_input_step
        sum_input_step         = 0.0
        sum_N_input_total      = sum_N_input_total+sum_N_input_step  
        sum_N_input_step       = 0.0
        sum_N_out_total        = sum_N_out_total+sum_N_out_step    
        sum_N_out_step         = 0.0

      end do !t time-loop
      
      !Check mass conservation for total simulation time: 
      call total_carbon_conservation(sum_input_total,HR_mass_accumulated, &
                                      pool_C_start_for_mass_cons,pool_C_final,&
                                      nlevels)
      call total_nitrogen_conservation(sum_N_input_total,sum_N_out_total, &
                                      pool_N_start_for_mass_cons,pool_Ninorg_start_for_mass_cons, &
                                      pool_N_final,inorg_N_final,nlevels)
      
      call check(nf90_close(writencid))

      if ( Spinup_run ) then
        call check(nf90_close(Spinupncid))
      end if

      call print_summary(save_N, save_C,c1a,c1b,c2,c3a,c3b,c4a,c4b,pool_C_final,pool_N_final,inorg_N_final)

      !deallocation
      deallocate(ndep_prof,leaf_prof,froot_prof,norm_froot_prof,croot_prof,stem_prof,r_moist, NH4_sorp_eq_vr)
      deallocate(CUE_bacteria_vr,CUE_fungi_vr, CUE_ecm_vr,CUE_am_vr,f_enzprod)    

      !For timing
      call system_clock(count=clock_stop)      ! Stop Timer
      print*, "Total time for decomp subroutine in minutes: ", (real(clock_stop-clock_start)/real(clock_rate))/60

  end subroutine decomp

  subroutine print_summary(discarded_N, discarded_C,imm_a,imm_b,min,f_min_a,f_min_b,b_min_a,b_min_b,matrixC,matrixN,matrixInorg)
    !in: 
    real(r8), intent(in) :: discarded_N, discarded_C
    integer,  intent(in) :: imm_a,imm_b,min,f_min_a,f_min_b,b_min_a,b_min_b
    real(r8), dimension(:,:) :: matrixC,matrixN,matrixInorg


    call disp("pool_matrixC gC/m3 ",matrixC, 'F12.5')
    call disp("pool_matrixN gN/m3 ",matrixN, 'F12.5')          
    call disp("Inorganic N gN/m3 ",matrixInorg, 'F9.6')          
    !call disp("C:N : ",pool_matrixC/pool_matrixN)
    print*, "AMOUNT OF N DISCARDED: ", discarded_N
    print*, "AMOUNT OF C DISCARDED: ", discarded_C   
    print*, "Immobilization, not enough N:",imm_a
    print*, "Immobilization enough N: ",imm_b
    print*, "Mineralization: ",min
    print*, "Bacteria needs, fungi mineralize, not enough N: ",f_min_a
    print*, "Bacteria needs, fungi mineralize, enough N: ",f_min_b
    print*, "Fungi needs, bacteria mineralize, not enough N: ",b_min_a
    print*, "Fungi needs, bacteria mineralize, enough N: ",b_min_b
  end subroutine print_summary

  function calc_nitrification(nh4,t_scalar,w_scalar,soil_temp) result(f_nit)
    ! Based on /CTSM/blob/master/src/soilbiogeochem/SoilBiogeochemNitrifDenitrifMod.F90
    !IN:
    real(r8),intent(in) :: nh4 !gN/m3
    real(r8),intent(in) :: t_scalar !From CLM
    real(r8),intent(in) :: w_scalar !From CLM
    real(r8),intent(in) :: soil_temp !From CLM
    
    !Out:
    real(r8) :: f_nit

     !local
    real(r8)           :: anaerobic_frac
    real(r8),parameter :: pH  = 6.5_r8
    real(r8),parameter :: rpi = 3.14159265358979323846_R8
    real(R8),parameter :: SHR_CONST_TKFRZ   = 0.0_R8! freezing T of fresh water          ~ degC
    real(r8),parameter :: k_nitr_max = 0.1_r8/24._r8 !from paramfile ctsm51_params.c210528.nc = 0.1/day, converted to /hour
    real(r8)           :: k_nitr_t_vr,k_nitr_ph_vr,k_nitr_h2o_vr,k_nitr
    ! follows CENTURY nitrification scheme (Parton et al., (2001, 1996))

    ! assume nitrification temp function equal to the HR scalar
    k_nitr_t_vr = min(t_scalar, 1._r8)

    ! ph function from Parton et al., (2001, 1996)
    k_nitr_ph_vr = 0.56_r8 + atan(rpi * 0.45_r8 * (-5._r8+pH)/rpi)

    ! moisture function-- assume the same moisture function as limits heterotrophic respiration
    ! Parton et al. base their nitrification- soil moisture rate constants based on heterotrophic rates-- can we do the same?
    k_nitr_h2o_vr = w_scalar

    ! nitrification constant is a set scalar * temp, moisture, and ph scalars
    ! note that k_nitr_max_perday is converted from 1/day to 1/s
    k_nitr = k_nitr_max * k_nitr_t_vr * k_nitr_h2o_vr * k_nitr_ph_vr
    ! first-order decay of ammonium pool with scalar defined above
    f_nit = max(nh4 * k_nitr, 0._r8) !g/m3 h

    anaerobic_frac=0._r8 !NOTE: Assume always aerobic conditions here

    ! limit to oxic fraction of soils
    f_nit  = f_nit*(1._r8 - anaerobic_frac)

    !limit to non-frozen soil layers
    if ( soil_temp <= SHR_CONST_TKFRZ ) then
        f_nit = 0._r8
    end if 
  end function calc_nitrification

  function set_N_dep(CLMdep,const_dep) result(Dep)

    !Optional in:
    real(r8), optional :: const_dep
    real(r8), optional :: CLMdep

    !OUT
    real(r8)           :: Dep
    if (present(CLMdep) .and. .not. present(const_dep)) then
      Dep = CLMdep
    elseif (present(const_dep) .and. .not. present(CLMdep)) then
      Dep = const_dep
    elseif (.not. present(CLMdep) .and. .not. present(const_dep)) then
      print*, "N dep not set correctly, stopping"
      Dep = -999
      stop
    end if
  end function set_N_dep

  function calc_Leaching_Runoff(qdrain,qrunoff,h2osoi_liq, N_NO3,layer_nr) result(Leach_Run) 
    !Based on CTSM/src/soilbiogeochem/SoilBiogeochemNLeachingMod.F90

    !IN:
    real(r8),intent(in)           :: qdrain      !mmH2O/h = kgH2O/m2 h, From CLM, converted to 1/h from 1/s
    real(r8),intent(in)           :: qrunoff     !mmH2O/h = kgH2O/m2 h, From CLM, converted to 1/h from 1/s
    real(r8),dimension(:),intent(in):: h2osoi_liq !kgH2O/m2, FROM CLM
    real(r8),intent(in)           :: N_NO3       !gN/m3
    integer,intent(in)            :: layer_nr   ! -

    !Out
    real(r8)           :: Leach_Run       !gN/m3 h, NO3 leaching and runoff 
    
    !Local: 
    real(r8)           :: DIN       !gN/m3 (=disn_conc)
    real(r8)           :: Leach       !gN/m3 h
    real(r8)           :: Runoff       !gN/m3 h
    real(r8),parameter :: depth_runoff_Nloss = 0.05 ! (m) depth over which runoff mixes with soil water for N loss to runoff
    real(r8)           :: h2o_tot     !kgH2O/m2, 
    real(r8)           :: h2o_surface     !kgH2O/m2, 
    integer            :: i !counter

    h2o_tot=0.0
    do i = 1, nlevels ! Total liquid water in column, used for calculating leaching
      h2o_tot=h2o_tot+h2osoi_liq(i) !kg/m2 !TODO: Unneccecary to run this loop every timestep
    end do  
    h2o_surface = h2osoi_liq(1)+0.75*h2osoi_liq(2) !Liquid Water in the top 5cm

    DIN = 0._r8
    if (h2osoi_liq(layer_nr) > 0._r8) then 
      DIN = (N_NO3 * delta_z(layer_nr) )/(h2osoi_liq(layer_nr) )
    end if
    !TODO: Can be shortened to N_NO3*qdrain/H2o_tot. Is that better with regards to round off errors etc?
    Leach = min(DIN*qdrain*h2osoi_liq(layer_nr)/(h2o_tot*delta_z(layer_nr)), N_NO3/dt) !Cannot be higher than the amount of NO3 in layer
    Leach = max(Leach, 0._r8) ! limit the leaching flux to a positive value

    !Calculate runoff: 
    if (layer_nr == 1 ) then 
      Runoff = DIN*qrunoff*h2osoi_liq(layer_nr)/(h2o_surface*delta_z(layer_nr))
    elseif (layer_nr == 2 ) then
      Runoff = DIN*qrunoff*h2osoi_liq(layer_nr)*((depth_runoff_Nloss-0.02)/delta_z(layer_nr))/(h2o_surface*(depth_runoff_Nloss-0.02))
    else
      Runoff = 0.0
    end if 

    Runoff = min(Runoff,N_NO3/dt - Leach)  ! ensure that runoff rate isn't larger than soil N pool
    Runoff = max(Runoff, 0._r8)  ! limit the runoff flux to a positive value

    Leach_Run = Leach + Runoff
  end function calc_Leaching_Runoff

  function forward_MMK_flux(C_SAP,C_SUBSTRATE,MMK_nr) result(flux)
    !Compute C flux from substrate pool to saprotroph pool by using Michaelis Menten Kinetics.
    !NOTE: On the way, a fraction 1-CUE is lost as respiration. This is handeled in the "decomp" subroutine.

    !IN:
    real(r8), intent(in) :: C_SAP
    real(r8), intent(in) :: C_SUBSTRATE
    integer, intent (in) :: MMK_nr
    
    !OUT
    real(r8):: flux ![gC/(m3 hr)]
    
    !TODO: this works, but should not depend on Vmax & Km from mycmim mod
    flux = C_SAP*Vmax(MMK_nr)*C_SUBSTRATE/(Km(MMK_nr)+C_SUBSTRATE)
  end function forward_MMK_flux

  function reverse_MMK_flux(C_SAP,C_SUBSTRATE,MMK_nr) result(flux)
    !Compute C flux from substrate pool to saprotroph pool by using Michaelis Menten Kinetics.
    !NOTE: On the way, a fraction 1-CUE is lost as respiration. This is handeled in the "decomp" subroutine.
    !IN:
    real(r8), intent(in) :: C_SAP
    real(r8), intent(in) :: C_SUBSTRATE
    integer, intent (in) :: MMK_nr

    !Out
    real(r8):: flux ![gC/(m3 hr)]

    !TODO: this works, but should not depend on Vmax & Km from mycmim mod
    flux = C_SAP*Vmax(MMK_nr)*C_SUBSTRATE/(Km(MMK_nr)+C_SAP)
  end function reverse_MMK_flux

  subroutine mining_rates_Sulman(C_EcM,C_substrate,N_substrate,moisture_function,T,mining_mod, D_Cmine,D_Nmine) !Sulman et al 2019 eq 34-35 + max modifier
    !NOTE: T dependence (Arrhenius) seems a bit weird, makes flux very low...
    !NOTE: V_max(T) in article, but not sure how this temperature dependence is?
    !INPUT
    real(r8),intent(in) :: C_EcM
    real(r8),intent(in) :: C_substrate
    real(r8),intent(in) :: N_substrate
    real(r8),intent(in) :: moisture_function
    real(r8),intent(in) :: mining_mod
    real(r8), intent(in) :: T !Kelvin
    
    !OUTPUT
    real(r8),intent(out) :: D_Cmine
    real(r8),intent(out):: D_Nmine
    
    !LOCAL
    real(r8),parameter :: V_max = 0.3/hr_pr_yr !Sulman 2019 supplement page 7, Assumed SOMp,SOMc ~ slow SOM
    real(r8),parameter :: K_m = 0.015 
    real(r8),parameter :: E_a = 54000 !J/mol
    real(r8),parameter :: R   = 8.31 !J/(K mol)

    D_Cmine = V_max*exp(-E_a/(R*T))*moisture_function*C_substrate*((C_EcM/C_substrate)/(C_EcM/C_substrate+K_m))*mining_mod
    D_Nmine = V_max*exp(-E_a/(R*T))*moisture_function*N_substrate*((C_EcM/C_substrate)/(C_EcM/C_substrate+K_m))*mining_mod 

  end subroutine mining_rates_Sulman

  subroutine mining_rates_Baskaran(C_EcM,C_substrate,N_substrate,mining_mod,soil_depth,D_Cmine,D_Nmine) !Baskaran + max_myc_alloc modifier
    !INPUT
    real(r8),intent(in) :: C_EcM
    real(r8),intent(in) :: C_substrate
    real(r8),intent(in) :: N_substrate
    real(r8),intent(in) :: mining_mod
    real(r8),intent(in) :: soil_depth
    
    !OUTPUT
    real(r8),intent(out):: D_Cmine !C "released" during mining, ends up in SOMa pool 
    real(r8),intent(out):: D_Nmine !N from SOM to EcM (mined N)

    if ( .not. (C_substrate < epsilon(C_substrate)) ) then !TODO: Review this epsilon thing..
      D_Cmine = K_MO*soil_depth*C_EcM*C_substrate*mining_mod
      D_Nmine = D_Cmine*N_substrate/C_substrate
    else 
      D_Cmine = 0.0_r8
      D_Nmine = 0.0_r8
    end if
  end subroutine mining_rates_Baskaran

  subroutine myc_to_plant(CUE_EcM,CUE_AM,enzyme_prod,NAMPlant,NEcMPlant) !Calculate rates of N flow from mycorrhiza to plant
    
   !INOUT: 
   real(r8), intent(inout) :: CUE_EcM ![-]
   real(r8), intent(inout) :: CUE_AM  ![-]
   real(r8), intent(inout) :: enzyme_prod ![-]
   
   !OUTPUT
   real(r8), intent(out)   :: NAMPlant  ![gN/m3 h]
   real(r8), intent(out)   :: NEcMPlant ![gN/m3 h]
   
   !LOCAL
   real(r8) ::     AM_N_demand ![gN/m3 h]
   real(r8) ::     AM_N_uptake ![gN/m3 h]
   real(r8) ::     EcM_N_demand ![gN/m3 h]
   real(r8) ::     EcM_N_uptake ![gN/m3 h]

   !All N the Mycorrhiza dont need for its own, it gives to the plant:
   AM_N_demand = CUE_AM*C_PlantAM/CN_ratio(6) 
   AM_N_uptake = N_INAM     

   if ( AM_N_uptake >= AM_N_demand ) then   
     NAMPlant = AM_N_uptake - AM_N_demand
   else !Reduce efficiency 
     NAMPlant = (1-f_growth)*AM_N_uptake
     CUE_AM = f_growth*AM_N_uptake*CN_ratio(6)/(C_PlantAM)
   end if
   if ( abs(NAMPlant) < 1e-16 ) then !TODO: How low/high should this value be?
     save_N=save_N+NAMPlant
     NAMPlant=0.0
   end if

   !All N the Mycorrhiza dont need for its own, it gives to the plant:
   EcM_N_demand = (CUE_EcM*(1-enzyme_prod)*C_PlantEcM)/CN_ratio(5)
   EcM_N_uptake = N_INEcM + N_SOMpEcM + N_SOMcEcM 
   if ( EcM_N_uptake >= EcM_N_demand ) then   
       NEcMPlant=EcM_N_uptake-EcM_N_demand      
   else !reduce efficiency or enzyme production, determined by option in namelist file.
       NEcMPlant = (1-f_growth)*EcM_N_uptake
       if ( use_ENZ ) then
         enzyme_prod = 1 - (f_growth*EcM_N_uptake*CN_ratio(5))/(CUE_EcM*C_PlantEcM)
       else
         CUE_EcM = (f_growth*EcM_N_uptake*CN_ratio(5))/((1-enzyme_prod)*C_PlantEcM)
       end if
   end if
   if ( abs(NEcMPlant) < 1e-16 ) then
     save_N=save_N+NEcMPlant
     
     NEcMPlant=0.0
   end if

 end subroutine myc_to_plant 

  subroutine calculate_fluxes(depth,Temp_Celsius,water_content,C_pool_matrix,N_pool_matrix, &
                              N_inorg_matrix,Deposition_rate, Leaching_rate, nitrification, &
                              soil_depth,desorp) !This subroutine calculates the fluxes in and out of the SOM pools.
    integer,intent(in)        :: depth !depth level number (not depth in meters!)  
    real(r8), intent(in)      :: Temp_Celsius !deg C, from CLM
    real(r8),intent(in)       :: Deposition_rate ![gN/m3 h]
    real(r8),intent(in)       :: Leaching_rate ![gN/m3 h]
    real(r8),intent(in)       :: nitrification ![gN/m3 h]
    real(r8),intent(in)       :: water_content !!m3water/m3soil, from CLM, used for sorption/desorption
    real(r8),intent(in)       :: desorp
    real(r8),intent(in)       :: soil_depth
    
    real(r8),target :: C_pool_matrix(nlevels, pool_types) !These are "targets" so that we can make "pointers" to them below- 
    real(r8),target :: N_pool_matrix(nlevels, pool_types_N)
    real(r8),target :: N_inorg_matrix(nlevels, inorg_N_pools)
    
    !LOCAL:
    real(r8)  :: NH4_sol_tmp ![gN/m3]
    real(r8)  :: NO3_tmp ![gN/m3] 
    real(r8)  :: NH4_tot ![gN/m3]
    real(r8)  :: nh4_sol_frac ![-] Fraction of avail inorganic N that is NH4
    real(r8)  :: N_IN !NO3 + NH4_sol = available inorganic N 
    real(r8)  :: f_b ![-] used to calculate N exchange rates, gives how much inorg N is avail to SAPb and SAPf, respevtely
    real(r8)  :: minedSOMp
    real(r8)  :: minedSOMc
    real(r8)  :: U_sb, U_sf ![gC/m3 h] Uptake of C by SAP
    real(r8)  :: UN_sb,UN_sf![gN/m3 h] Uptake of N by SAP

    real(r8)  :: Temp_Kelvin ![K]

    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: C_LITm, C_LITs, C_SOMp,C_SOMa,C_SOMc,C_EcM,C_AM, C_SAPb, C_SAPf, &
                         N_LITm, N_LITs, N_SOMp,N_SOMa,N_SOMc,N_EcM,N_AM, N_SAPb, N_SAPf, &
                         N_NH4_sol,N_NH4_sorp,N_NO3
    C_LITm => C_pool_matrix(depth, 1)
    C_LITs => C_pool_matrix(depth, 2)
    C_SAPb => C_pool_matrix(depth, 3)
    C_SAPf => C_pool_matrix(depth, 4)
    C_EcM =>  C_pool_matrix(depth, 5)
    C_AM =>   C_pool_matrix(depth, 6)
    C_SOMp => C_pool_matrix(depth, 7)
    C_SOMa => C_pool_matrix(depth, 8)
    C_SOMc => C_pool_matrix(depth, 9)
    
    N_LITm => N_pool_matrix(depth, 1)
    N_LITs => N_pool_matrix(depth, 2)
    N_SAPb => N_pool_matrix(depth, 3)
    N_SAPf => N_pool_matrix(depth, 4)
    N_EcM =>  N_pool_matrix(depth, 5)
    N_AM =>   N_pool_matrix(depth, 6)
    N_SOMp => N_pool_matrix(depth, 7)
    N_SOMa => N_pool_matrix(depth, 8)
    N_SOMc => N_pool_matrix(depth, 9)
    
    N_NH4_sol  => N_inorg_matrix(depth, 1)
    N_NH4_sorp => N_inorg_matrix(depth, 2)
    N_NO3      => N_inorg_matrix(depth,3)
    
    Temp_Kelvin = Temp_Celsius+abs_zero
    
    !------------------CARBON FLUXES----------------------------:
    !Decomposition of LIT and SOMa by SAP:
    !On the way, a fraction 1-CUE is lost as respiration. This is handeled in the "decomp" subroutine.
    C_LITmSAPb=reverse_MMK_flux(C_SAPb,C_LITm,1) !C5
    C_LITsSAPb=reverse_MMK_flux(C_SAPb,C_LITs,2) !C6
    C_SOMaSAPb=reverse_MMK_flux(C_SAPb,C_SOMa,3) !C7
    C_LITmSAPf=reverse_MMK_flux(C_SAPf,C_LITm,4) !C8
    C_LITsSAPf=reverse_MMK_flux(C_SAPf,C_LITs,5) !C9
    C_SOMaSAPf=reverse_MMK_flux(C_SAPf,C_SOMa,6) !C10
    
    !Oxidation from SOMc to SOMa
    !From equations for decomposing structural litter in mimics,eq. A10
    !KO modifies Km which is used in the litter->SAP equations.
    C_SOMcSOMa =  ( C_SAPb * Vmax(2) * C_SOMc / (KO(1)*Km(2) + C_SAPb)) + &
                  ( C_SAPf * Vmax(5) * C_SOMc / (KO(2)*Km(5) + C_SAPf)) !C11
    
    !Desorbtion controls transport from physically protected to available SOM
    C_SOMpSOMa=C_SOMp*desorp !C12
    
    !Turnover from SAP to SOM. Based on the turnover equations used in mimics for flux from microbial pools to SOM pools (correspond to eq A4,A8 in Wieder 2015)
    C_SAPbSOMp=C_SAPb*k_sapsom(1)*fPHYS(1)   !gC/m3h !C13
    C_SAPbSOMc=C_SAPb*k_sapsom(1)*fCHEM(1) !C14
    C_SAPbSOMa=C_SAPb*k_sapsom(1)*fAVAIL(1)!C15
    
    C_SAPfSOMp=C_SAPf*k_sapsom(2)*fPHYS(2) !C16 
    C_SAPfSOMc=C_SAPf*k_sapsom(2)*fCHEM(2) !C17
    C_SAPfSOMa=C_SAPf*k_sapsom(2)*fAVAIL(2)!C18
    
    !Dead mycorrhizal biomass enters the SOM pools:  gC/m3h
    C_EcMSOMp=C_EcM*k_mycsom(1)*fEcMSOM(1)!somp !C19
    C_EcMSOMc=C_EcM*k_mycsom(1)*fEcMSOM(2)!somc !C20
    C_EcMSOMa=C_EcM*k_mycsom(1)*fEcMSOM(3)!soma !C21
    
    C_AMSOMp=C_AM*k_mycsom(2)*fAMSOM(1) !C22
    C_AMSOMc=C_AM*k_mycsom(2)*fAMSOM(2) !C23
    C_AMSOMa=C_AM*k_mycsom(2)*fAMSOM(3) !C24

    !Ectomycorrhizal mining options (use Baskaran) :
    if ( use_Sulman ) then
      call mining_rates_Sulman(C_EcM,C_SOMc,N_SOMc,r_moist(depth),Temp_Kelvin,r_myc, minedSOMc,N_SOMcEcM)
      call mining_rates_Sulman(C_EcM,C_SOMp,N_SOMp,r_moist(depth),Temp_Kelvin, r_myc, minedSOMp,N_SOMpEcM)
    else
      call mining_rates_Baskaran(C_EcM,C_SOMp,N_SOMp,r_myc,soil_depth,minedSOMp,N_SOMpEcM) !N25
      call mining_rates_Baskaran(C_EcM,C_SOMc,N_SOMc,r_myc,soil_depth,minedSOMc,N_SOMcEcM) !N26              
    end if
    
    C_EcMdecompSOMp = minedSOMp   ![gC/m3h] !C25 !NOTE Can drop minedSOMx and define C_EcMdecompSOMx directly
    C_EcMdecompSOMc = minedSOMc   ![gC/m3h] !C26
    
    !-----------------------------------NITROGEN FLUXES----------------------------:
    !Decomposition of LIT and SOMa by SAP
    N_LITmSAPb = calc_parallel_Nrates(C_LITmSAPb,N_LITm,C_LITm) !N5
    N_LITsSAPb = calc_parallel_Nrates(C_LITsSAPb,N_LITs,C_LITs) !N6
    N_SOMaSAPb = calc_parallel_Nrates(C_SOMaSAPb,N_SOMa,C_SOMa) !N7
    N_LITmSAPf = calc_parallel_Nrates(C_LITmSAPf,N_LITm,C_LITm) !N8
    N_LITsSAPf = calc_parallel_Nrates(C_LITsSAPf,N_LITs,C_LITs) !N9
    N_SOMaSAPf = calc_parallel_Nrates(C_SOMaSAPf,N_SOMa,C_SOMa) !N10
    !Transport from SOMc to SOMa:
    N_SOMcSOMa = calc_parallel_Nrates(C_SOMcSOMa,N_SOMc,C_SOMc) !N11
    !Desorption of SOMp to SOMa
    N_SOMpSOMa = calc_parallel_Nrates(C_SOMpSOMa,N_SOMp,C_SOMp) !N12
    !Dead saphrotroph biomass enters SOM pools
    N_SAPbSOMp = calc_parallel_Nrates(C_SAPbSOMp,N_SAPb,C_SAPb) !N13
    N_SAPbSOMc = calc_parallel_Nrates(C_SAPbSOMc,N_SAPb,C_SAPb) !N14
    N_SAPbSOMa = calc_parallel_Nrates(C_SAPbSOMa,N_SAPb,C_SAPb) !N15
    N_SAPfSOMp = calc_parallel_Nrates(C_SAPfSOMp,N_SAPf,C_SAPf) !N16
    N_SAPfSOMc = calc_parallel_Nrates(C_SAPfSOMc,N_SAPf,C_SAPf) !N17
    N_SAPfSOMa = calc_parallel_Nrates(C_SAPfSOMa,N_SAPf,C_SAPf) !N18
    !Dead mycorrhizal biomass enters SOM pools
    N_EcMSOMp = calc_parallel_Nrates(C_EcMSOMp,N_EcM,C_EcM) !N19
    N_EcMSOMa = calc_parallel_Nrates(C_EcMSOMa,N_EcM,C_EcM) !N20
    N_EcMSOMc = calc_parallel_Nrates(C_EcMSOMc,N_EcM,C_EcM) !N22
    N_AMSOMp = calc_parallel_Nrates(C_AMSOMp,N_AM,C_AM) !N22
    N_AMSOMa = calc_parallel_Nrates(C_AMSOMa,N_AM,C_AM) !N23
    N_AMSOMc = calc_parallel_Nrates(C_AMSOMc,N_AM,C_AM) !N24

    !******************************Calculating fluxes related to inorganic N: ***********************************************

    !(1) Update inorganic pools to account for Leaching, deposition,  nitrification rate and gain from decomposition (1-NUE):
    NH4_sol_tmp= N_NH4_sol + (1-NUE)*(N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf+N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb)*dt + (Deposition_rate - nitrification)*dt
    NO3_tmp    = N_NO3-Leaching_rate*dt + nitrification*dt
    call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
        
    !(2) Inorganic N taken up directly by plant roots:
    N_InPlant = calc_plant_uptake(N_IN) !N34
    
    !Update inorganic pools to account for uptake by plant roots
    NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*N_InPlant*dt
    NO3_tmp = max(NO3_tmp - (1-nh4_sol_frac)*N_InPlant*dt,0._r8)
    call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
     
    !(3) Inorganic N taken up by mycorrhiza 
    N_INEcM  = calc_myc_uptake(N_IN,C_EcM,soil_depth) !N27
    N_INAM   = calc_myc_uptake(N_IN,C_AM,soil_depth) !N28
    
    !Update inorganic pools to account for uptake by mycorrhizal fungi
    NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*(N_INEcM+N_INAM)*dt
    NO3_tmp = max(NO3_tmp - (1-nh4_sol_frac)*(N_INEcM+N_INAM)*dt,0._r8)
    call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)

    !(4) Calculate exchange of N between saprotrophic pools and available inorganic N. This 
    !This can be positive or negative depending on wether the saprotrophs immobilizes or mineralizes N 
    !in the decomposition process.

    !total C uptake (growth + respiration) of saprotrophs
    U_sb = C_LITmSAPb + C_LITsSAPb + C_SOMaSAPb  
    U_sf = C_LITmSAPf + C_LITsSAPf + C_SOMaSAPf    
    ! N uptake by saprotrophs
    UN_sb = (N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb)*NUE 
    UN_sf = (N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf)*NUE

    !SAP demand for N based on target CN ratio:
    N_demand_SAPb =  CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
    N_demand_SAPf =  CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
    
    !How much N saprotrophs need from the inorganic pool to meet target CN ratio
    N_INSAPb = N_demand_SAPb-UN_sb
    N_INSAPf = N_demand_SAPf-UN_sf
    
    !Determine exchange of N between inorganic pool and saprotrophs, N_INSAPb and N_INSAPf: !N36, !N37 is determined here.
    if ( N_INSAPb >= 0. .and. N_INSAPf >= 0. ) then !immobilization
      if ( N_IN < (N_INSAPb + N_INSAPf)*dt) then !Not enough inorganic N to meet demand
        

        f_b = N_INSAPb/(N_INSAPb + N_INSAPf) ! Bac. and fungi want the same inorganic N. This fraction determines how much N is available to each pool.
        if ( U_sb ==0._r8 ) then !To avoid division by zero 
          N_INSAPb =0._r8
        else
          CUE_bacteria_vr(depth)=((f_b*N_IN+UN_sb*dt)*CN_ratio(3))/(U_sb*dt)
          N_demand_SAPb =  CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
          N_INSAPb = f_b*N_IN/dt
        end if
        if ( U_sf ==0._r8 ) then !To avoid division by zero
          N_INSAPf = 0._r8
        else
          CUE_fungi_vr(depth) = (((1-f_b)*N_IN+UN_sf*dt)*CN_ratio(4))/(U_sf*dt)
          N_demand_SAPf =  CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
          N_INSAPf = (1-f_b)*N_IN/dt
        end if
   
        c1a=c1a+1 !Count immob, not enough N occurence
      else !Enough mineral N to meet demand
        
        c1b=c1b+1 !Count immob, enough N 
        continue
      end if    

      !Update inorganic pools
      NO3_tmp = NO3_tmp - (1-nh4_sol_frac)*(N_INSAPb + N_INSAPf)
      NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*(N_INSAPb + N_INSAPf)
      call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
        
    elseif ( N_INSAPb < 0. .and. N_INSAPf < 0. ) then !mineralization
      NH4_sol_tmp = NH4_sol_tmp - (N_INSAPb + N_INSAPf)
      call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)      
      c2=c2+1 !count mineralization occurence
      continue     
      
    elseif ( N_INSAPb >= 0. .and. N_INSAPf < 0. ) then ! bacteria can use N mineralized by fungi
      if ( (N_IN +abs(N_INSAPf)*dt) < N_INSAPb*dt ) then
        if ( U_sb ==0._r8 ) then !To avoid division by zero 
          N_INSAPb =0._r8
        else
          CUE_bacteria_vr(depth)=(((N_IN + abs(N_INSAPf)*dt)+UN_sb*dt)*CN_ratio(3))/(U_sb*dt)
          N_demand_SAPb =  CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
          N_INSAPb = N_demand_SAPb/dt-UN_sb   
        end if 
        c3a=c3a+1             
      else
        c3b=c3b+1 
      end if
      
      !Update inorganic N pools
      NO3_tmp = NO3_tmp - (1-nh4_sol_frac)*(N_INSAPb + N_INSAPf)
      NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*(N_INSAPb + N_INSAPf)
      call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
      
    elseif ( N_INSAPb < 0. .and. N_INSAPf >= 0. ) then !fungi can use N mineralized by bacteria
      if ( (N_IN+ abs(N_INSAPb)*dt) < N_INSAPf*dt ) then
        if ( U_sf == 0._r8) then
          N_INSAPf = 0._r8
        else 
          CUE_fungi_vr(depth)=(( (N_IN + abs(N_INSAPb)*dt)+UN_sf*dt)*CN_ratio(4))/(U_sf*dt)
          N_demand_SAPf =  CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
          N_INSAPf = N_demand_SAPf/dt-UN_sf
          N_INSAPf=N_IN/dt+abs(N_INSAPb)
          !nh4_sol_frac=calc_nh4_frac(NH4_sol_tmp+abs(N_INSAPb)*dt,NO3_tmp)
        end if
        
        c4a=c4a+1         
      else
        c4b=c4b+1 
      end if
      !Update inorganic N pools:
      NO3_tmp = NO3_tmp - (1-nh4_sol_frac)*(N_INSAPb + N_INSAPf)
      NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*(N_INSAPb + N_INSAPf)
      call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
                  
    else 
      print*, "No condition applies for N_INSAPb, N_INSAPf calculation (this should not happen); ", C_LITmSAPf ,C_LITsSAPf , C_SOMaSAPf,C_LITmSAPb, C_LITsSAPb, C_SOMaSAPb, depth
      stop
      
    end if

    !For determining sorption/desorption of NH4:
    NH4_tot = NH4_sol_tmp + N_NH4_sorp    
    call calc_NH4_sol_sorp(NH4_tot,water_content,N_NH4_sorp,NH4_sorp_eq_vr(depth),NH4_sorp_final)
    NH4_sol_final = max(NH4_sol_tmp - (NH4_sorp_final-N_NH4_sorp),0._r8)
    NO3_final = max(NO3_tmp,0._r8)  
    
    !reset values:
    N_IN = 0._r8 
    NH4_sol_tmp = 0._r8 
    NO3_tmp = 0._r8 
    nullify( C_LITm,C_LITs,C_SOMp,C_SOMa,C_SOMc,C_EcM,C_AM, C_SAPb,C_SAPf)
    nullify( N_LITm,N_LITs,N_SOMp,N_SOMa,N_SOMc,N_EcM,N_AM, N_SAPb,N_SAPf,N_NH4_sol,N_NH4_sorp,N_NO3)
  end subroutine calculate_fluxes

  subroutine calc_NH4_sol_sorp(NH4_tot,soil_water_frac,NH4_sorp_previous,NH4_sorp_eq,NH4_sorp)
    !IN:
    real(r8), intent(in)  :: NH4_tot   !g/m3, total NH4, both in soil solution and adsorbed
    real(r8), intent(in)  :: soil_water_frac   !m3water/m3soil (input from CLM data)
    real(r8), intent(in)  :: NH4_sorp_previous  !g/m3

    !Out:
    real(r8),intent(out)            :: NH4_sorp !g/m3, NH4 sorbed to particles
    real(r8),intent(out)            :: NH4_sorp_eq !g/m3, adsorbed NH4 at equilibrium at current concentration
    
    !lOCAL: 
    real(r8), parameter :: BD_soil      =1.6e6  !g/m3 (loam) soil from DOI: 10.3390/APP6100269 Table 1
    real(r8), parameter :: NH4_sorp_max = 0.09*BD_soil/mg_pr_g    !mg NH4 /g soil
    real(r8), parameter :: KL           = 0.4      !L/mg
    real(r8), parameter :: K_pseudo     = 0.0167*mg_pr_g*60./BD_soil !m3/(g hour) (60min/hr)
    real(r8)            :: KL_prime       !m3/g

    !1) Calculate NH4_sorp_eq 
    KL_prime = KL*mg_pr_g*m3_pr_L/soil_water_frac 
    
    NH4_sorp_eq=(1+KL_prime*NH4_tot+NH4_sorp_max*KL_prime)/(2*KL_prime) - sqrt((1+KL_prime*NH4_tot+NH4_sorp_max*KL_prime)**2-4*KL_prime**2*NH4_sorp_max*NH4_tot)/(2*KL_prime)
    
    !2) Calculate NH4_sorp after adjusting towards equilibrium for 1 timestep
    if ( NH4_sorp_eq==NH4_sorp_previous ) then !Already at equilibrium
      NH4_sorp = NH4_sorp_previous
    elseif (NH4_sorp_eq > NH4_sorp_previous ) then !Adsorption
      NH4_sorp = NH4_sorp_eq - 1_r8/(1._r8/(NH4_sorp_eq-NH4_sorp_previous) + k_pseudo*dt)
    else !Desorption
      NH4_sorp = NH4_sorp_eq + 1_r8/(1._r8/(NH4_sorp_previous-NH4_sorp_eq) + k_pseudo*dt)
    end if
  end subroutine calc_NH4_sol_sorp

  subroutine vertical_diffusion(pool_matrix,vert,Diffusivity) !This subroutine calculates the vertical transport of carbon through the soil layers.
    !IN 
    real(r8), intent(in)   :: pool_matrix(:,:)
    real(r8), intent(in)   :: Diffusivity ![m2/h] Diffusivity
    !OUT
    real(r8), allocatable, intent(out)  :: vert(:,:)
    
    !Local
    integer                :: depth, pool !For iteration
    integer,dimension(1)   :: max_pool, max_depth !For iteration
    real(r8)  :: upper_diffusion_flux, lower_diffusion_flux
    real(r8)  :: tot_diffusion_dummy ![gC/h]
    real(r8)  :: D

    allocate (vert, mold = pool_matrix)

    !Get how many depth levels and pools we will loop over.
    max_depth=shape(pool_matrix(:,1)) !TODO: Easier way to do this?
    max_pool=shape(pool_matrix(1,:))
    !In a timestep, the fluxes between pools in the same layer is calculated before the vertical diffusion. Therefore, a loop over all the entries in
    !pool_matrix is used here to calculate the diffusion.
    do depth = 1,max_depth(1)
      do pool =1, max_pool(1)
        D = Diffusivity
        if ( max_pool(1)==3 .and. pool ==2 ) then
          D = Diffusivity/3
        end if
        !eq. 6.18 and 6.20 from Soetaert & Herman, A practical guide to ecological modelling.
        if (depth == 1) then
          upper_diffusion_flux= 0._r8
          lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
        elseif (depth==max_depth(1)) then
          upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
          lower_diffusion_flux= 0._r8
        else
          upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
          lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
        end if
        tot_diffusion_dummy=(upper_diffusion_flux-lower_diffusion_flux)/delta_z(depth)
        vert(depth,pool) = tot_diffusion_dummy
      end do !pool
    end do !depth
  end subroutine vertical_diffusion

  subroutine input_rates(layer_nr, met_fraction,LEAFC_TO_LIT,FROOTC_TO_LIT, &
                        LEAFN_TO_LIT,FROOTN_TO_LIT,&
                        N_CWD,C_CWD, &
                        C_inLITm,C_inLITs,&
                        N_inLITm,N_inLITs, &
                        C_inSOMp,C_inSOMc, &
                        N_inSOMp,N_inSOMc, &
                        met_mortalityC, met_mortalityN, &
                        split_mortalityC, split_mortalityN)
                        
    !in:
    integer,  intent(in) :: layer_nr
    real(r8), intent(in) :: met_fraction
    real(r8), intent(in) :: LEAFC_TO_LIT !From CLM
    real(r8), intent(in) :: FROOTC_TO_LIT!From CLM
    real(r8), intent(in) :: LEAFN_TO_LIT !FROM CLM
    real(r8), intent(in) :: FROOTN_TO_LIT!From CLM    
    real(r8), intent(in) :: N_CWD(:) !FROM CLM
    real(r8), intent(in) :: C_CWD(:) !FROM CLM
    real(r8), intent(in) :: met_mortalityC(:)
    real(r8), intent(in) :: met_mortalityN(:)
    real(r8), intent(in) :: split_mortalityC(:)
    real(r8), intent(in) :: split_mortalityN(:)
    
    !out:
    real(r8), intent(out) :: C_inLITm
    real(r8), intent(out) :: C_inLITs
    real(r8), intent(out) :: N_inLITm
    real(r8), intent(out) :: N_inLITs
    real(r8), intent(out) :: C_inSOMp
    real(r8), intent(out) :: C_inSOMc
    real(r8), intent(out) :: N_inSOMp
    real(r8), intent(out) :: N_inSOMc
    
    !local:
    real(r8)           :: leaf_root_inputC
    real(r8)           :: leaf_root_inputN
    real(r8)           :: metabolic_litter_totC
    real(r8)           :: metabolic_litter_totN
    real(r8)           :: structural_litter_totC
    real(r8)           :: structural_litter_totN

    leaf_root_inputC = FROOTC_TO_LIT*froot_prof(layer_nr) + LEAFC_TO_LIT*leaf_prof(layer_nr) !gC/m3h !NOTE: These do not include CWD
    leaf_root_inputN = FROOTN_TO_LIT*froot_prof(layer_nr) + LEAFN_TO_LIT*leaf_prof(layer_nr)!gN/m3h

    metabolic_litter_totC = (met_fraction*(leaf_root_inputC + split_mortalityC(layer_nr)) + met_mortalityC(layer_nr))
    metabolic_litter_totN = (met_fraction*(leaf_root_inputN + split_mortalityN(layer_nr)) + met_mortalityN(layer_nr))
    
    structural_litter_totC = ((1-met_fraction)*(leaf_root_inputC + split_mortalityC(layer_nr))+  C_CWD(layer_nr))
    structural_litter_totN = ((1-met_fraction)*(leaf_root_inputN + split_mortalityN(layer_nr))+  N_CWD(layer_nr))

    C_inLITm = metabolic_litter_totC *(1-f_met_to_som) !C1
    N_inLITm = metabolic_litter_totN *(1-f_met_to_som) !N1
    
    C_inLITs = structural_litter_totC*(1-f_struct_to_som) !C2
    N_inLITs = structural_litter_totN*(1-f_struct_to_som) !N2
    
    C_inSOMp = metabolic_litter_totC *f_met_to_som !C3
    N_inSOMp = metabolic_litter_totN *f_met_to_som !N3
    
    C_inSOMc = structural_litter_totC*f_struct_to_som !C4
    N_inSOMc = structural_litter_totN*f_struct_to_som !N4
    
  end subroutine input_rates

  function calc_parallel_Nrates(C_rate,N_pool,C_pool) result(N_rate)
    !in: 
    real(r8),intent(in) :: C_rate
    real(r8),intent(in) :: N_pool
    real(r8),intent(in) :: C_pool
    
    !out
    real(r8) :: N_rate
      
    if ( C_pool <= epsilon(C_pool) ) then !Avoid division by zero
      N_rate=0.0
    else
      N_rate=C_rate*(N_pool/C_pool)
    end if
  end function calc_parallel_Nrates

  subroutine update_inorganic_N(NO3,NH4_sol,Ninorg_avail, ratio)
    real(r8), intent(in) :: NO3
    real(r8), intent(in) :: NH4_sol
    
    !out: 
    real(r8), intent(out) :: Ninorg_avail
    real(r8), intent(out) :: ratio
    
    Ninorg_avail = NO3+NH4_sol !The inorganic N available to microbes and plants (the rest, NH4_sorb_tmp is sorbed onto particles)
    if (NH4_sol+NO3 == 0._r8) Then 
      ratio = 0.5_8
    else
      ratio = NH4_sol/(NH4_sol+NO3)
    end if
  end subroutine update_inorganic_N

  function calc_plant_uptake(N_inorganic) result(N_INVeg)
    !IN 
    real(r8), intent(in) :: N_inorganic
    
    !out:
    real(r8)             :: N_INVeg

    N_INVeg = k_plant*N_inorganic
  end function calc_plant_uptake

  function calc_myc_uptake(N_inorganic,C_MYC,soil_depth) result(N_INMYC) !TODO: Fix soil_depth problem somehow..
    !IN 
    real(r8), intent(in) :: N_inorganic
    real(r8), intent(in) :: C_MYC
    real(r8), intent(in) :: soil_depth
    !OUT 
    real(r8) :: N_INMYC

    N_INMYC = V_max_myc*N_inorganic*(C_MYC/(C_MYC + Km_myc/soil_depth))*r_myc
    !NOTE: MMK parameters should maybe be specific to mycorrhizal type?
  end function calc_myc_uptake
  
end module mycmimMod
