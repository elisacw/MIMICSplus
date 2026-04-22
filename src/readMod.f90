module readMod
  use shr_kind_mod, only : r8 => shr_kind_r8
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  use netcdf
  use dispmodule , only: disp
  use paramMod,    only: sec_pr_hr,delta_z,CLM_version
  use writeMod,    only: check
  use initMod,     only: nlevels

  implicit none

  private
  public :: read_time,read_maxC,read_WATSAT_and_profiles,read_clay,read_PFTs,read_clm_model_input, &
            calc_PFT,read_namelist,read_mort_profiles, read_clm_mortality

  contains
      
    subroutine read_time(ncid,steps)
      !INPUT
      integer,intent(in):: ncid
      
      !OUTPUT
      integer, intent(out):: steps
      
      !LOCAL
      integer            :: timeid

      call check(nf90_inquire(ncid, unlimiteddimid = timeid))
      call check(nf90_inquire_dimension(ncid, timeid, len = steps))
    end subroutine read_time
      
    function read_maxC(ncid,time_steps) result(max_Cpay) !reads max C flux from plant to Myc. Used to determine myc_modifier
      !INput
      integer, intent(in)   :: ncid
      integer, intent(in)   :: time_steps
      
      !Out:
      real(r8)              :: max_Cpay
      !Local:
      integer    :: varid
      real(r8), dimension(time_steps) :: NPP_tot
      real(r8), dimension(time_steps) :: NPP_nonmyc
      real(r8), dimension(time_steps) :: NPP_myc
      
      call check(nf90_inq_varid(ncid, 'NPP_NACTIVE', varid))
      call check(nf90_get_var(ncid, varid, NPP_tot(:),start = (/1,1/),count = (/1,time_steps/)))
      call check(nf90_inq_varid(ncid, 'NPP_NNONMYC', varid))
      call check(nf90_get_var(ncid, varid, NPP_nonmyc(:),start = (/1,1/),count = (/1,time_steps/)))
      
      NPP_myc = (NPP_tot-NPP_nonmyc)*sec_pr_hr

      max_Cpay = maxval(NPP_myc)      
    end function read_maxC

    subroutine read_mort_profiles(clm_file,CROOT_PROF,STEM_PROF) !Needed bc. these vars is only given in the first outputfile of the simulation.

      !INPUT
      character (len = *),intent(in):: clm_file       
      !OUTPUT
      real(r8), intent(out),dimension(nlevels)          :: CROOT_PROF
      real(r8), intent(out),dimension(nlevels)          :: STEM_PROF
   
      !LOCAL
      integer            :: ncid, varid

      call check(nf90_open(trim(clm_file), nf90_nowrite, ncid))


      call check(nf90_inq_varid(ncid, 'CROOT_PROF', varid))
      call check(nf90_get_var(ncid, varid, CROOT_PROF,count=(/1,nlevels/))) 
      call check(nf90_inq_varid(ncid, 'STEM_PROF', varid))
      call check(nf90_get_var(ncid, varid, STEM_PROF,count=(/1,nlevels/))) 
      
      call check(nf90_close(ncid))

    end subroutine read_mort_profiles

    subroutine read_WATSAT_and_profiles(clm_history_file,WATSAT,NDEP_PROF,FROOT_PROF,LEAF_PROF) !Needed bc. these vars is only given in the first outputfile of the simulation.
      !TODO: Makes more sense to have one subroutine for profiles and one for WATSAT?
      !INPUT
      character (len = *),intent(in):: clm_history_file       
      !OUTPUT
      real(r8), intent(out),dimension(nlevels)          :: WATSAT
      real(r8), intent(out),dimension(nlevels)          :: NDEP_PROF
      real(r8), intent(out),dimension(nlevels)          :: FROOT_PROF
      real(r8), intent(out),dimension(nlevels)          :: LEAF_PROF
   
      !LOCAL
      integer            :: ncid, varid
      
      call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
      call check(nf90_inq_varid(ncid, 'watsat', varid))
      call check(nf90_get_var(ncid, varid, WATSAT,count=(/1,nlevels/)))
      
      call check(nf90_inq_varid(ncid, 'LEAF_PROF', varid))
      call check(nf90_get_var(ncid, varid, LEAF_PROF,count=(/1,nlevels/)))       
      call check(nf90_inq_varid(ncid, 'FROOT_PROF', varid))
      call check(nf90_get_var(ncid, varid, FROOT_PROF,count=(/1,nlevels/))) 
      call check(nf90_inq_varid(ncid, 'NDEP_PROF', varid))
      call check(nf90_get_var(ncid, varid, NDEP_PROF,count=(/1,nlevels/))) 
      
      call check(nf90_close(ncid))
    end subroutine read_WATSAT_and_profiles
    
    subroutine calc_PFT(clm_history_file,lflitcn_avg) !Used to calculate metabolic fraction of litter (f_met)
      !INPUT
      character (len = *),intent(in):: clm_history_file       
      !OUTPUT
      real(r8), intent(out)          :: lflitcn_avg

      !LOCAL
      integer            :: ncid, pftid
      real(r8),dimension(15)           :: PCT_NAT_PFT
      real(r8), dimension(15),parameter           :: lflitcn = (/1, 70, 80, 50, 60, 60, 50, 50, 50, 60, 50, 50, 50, 50, 50/) !from ctsm51_params.c211112.nc 
      
      call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
      call check(nf90_inq_varid(ncid, 'PCT_NAT_PFT', pftid))
      call check(nf90_get_var(ncid, pftid, PCT_NAT_PFT, count=(/1,1,15/)))
      call check(nf90_close(ncid))

      lflitcn_avg = sum(lflitcn*PCT_NAT_PFT/100.)

    end subroutine calc_PFT
    
    subroutine read_clay(clm_surface_file,mean_clay_content)
      !INPUT
      character (len = *),intent(in):: clm_surface_file
      
      !OUTPUT
      real(r8), intent(out)         :: mean_clay_content
      
      !LOCAL
      integer            :: ncid, clayid
      real(r8),dimension(10)           :: pct_clay

      call check(nf90_open(trim(clm_surface_file), nf90_nowrite, ncid))
      call check(nf90_inq_varid(ncid, 'PCT_CLAY', clayid))
      call check(nf90_get_var(ncid, clayid, pct_clay, count=(/1,1,10/)))
      call check(nf90_close(ncid))

      mean_clay_content = (sum(pct_clay)/size(pct_clay))/100.0 !output as fraction
    end subroutine read_clay
    
    subroutine read_PFTs(clm_surface_file,PFT_dist) !only used if use_ROI = False, so that myc. association is determined by PFT dist.
      !INPUT
      character (len = *),intent(in):: clm_surface_file
      
      !OUTPUT
      real(r8), dimension(15),intent(out)         :: PFT_dist
      
      !LOCAL
      integer            :: ncid, pftid

      call check(nf90_open(trim(clm_surface_file), nf90_nowrite, ncid))
      call check(nf90_inq_varid(ncid, 'PCT_NAT_PFT', pftid))
      call check(nf90_get_var(ncid, pftid, PFT_dist, count=(/1,1,15/)))
      call check(nf90_close(ncid))
    end subroutine read_PFTs

    subroutine read_clm_mortality(ncid,time_entry,froot_prof, croot_prof, leaf_prof, stem_prof, met_mortalityC,met_mortalityN, split_mortalityC, split_mortalityN)

      !INPUT
      integer,intent(in)                        :: ncid
      integer,intent(in)                        :: time_entry
      real(r8),intent(in), dimension(nlevels)   :: froot_prof
      real(r8),intent(in), dimension(nlevels)   :: croot_prof
      real(r8),intent(in), dimension(nlevels)   :: leaf_prof
      real(r8),intent(in), dimension(nlevels)   :: stem_prof

      !OUTPUT
      real(r8), dimension(nlevels) :: met_mortalityC
      real(r8), dimension(nlevels) :: split_mortalityC
      real(r8), dimension(nlevels) :: met_mortalityN
      real(r8), dimension(nlevels) :: split_mortalityN

      !LOCAL: 
      real(r8) :: met_leaf_prof_mortC
      real(r8) :: met_froot_prof_mortC
      real(r8) :: met_croot_prof_mortC
      real(r8) :: met_stem_prof_mortC
      real(r8) :: split_leaf_prof_mortC
      real(r8) :: split_froot_prof_mortC
      real(r8) :: met_leaf_prof_mortN
      real(r8) :: met_froot_prof_mortN
      real(r8) :: met_croot_prof_mortN
      real(r8) :: met_stem_prof_mortN
      real(r8) :: split_leaf_prof_mortN
      real(r8) :: split_froot_prof_mortN
      integer :: varid


      !C mortality fluxes to metabolic: 
      call check(nf90_inq_varid(ncid, 'met_leaf_prof_mortC', varid))
      call check(nf90_get_var(ncid, varid, met_leaf_prof_mortC,start=(/1, time_entry/)))
      met_leaf_prof_mortC = met_leaf_prof_mortC*sec_pr_hr  ![gC/m^2/h]
      
      call check(nf90_inq_varid(ncid, 'met_froot_prof_mortC', varid))
      call check(nf90_get_var(ncid, varid, met_froot_prof_mortC,start=(/1, time_entry/)))
      met_froot_prof_mortC = met_froot_prof_mortC*sec_pr_hr  ![gC/m^2/h]

      call check(nf90_inq_varid(ncid, 'met_croot_prof_mortC', varid))
      call check(nf90_get_var(ncid, varid, met_croot_prof_mortC,start=(/1, time_entry/)))
      met_croot_prof_mortC = met_croot_prof_mortC*sec_pr_hr  ![gC/m^2/h]

      call check(nf90_inq_varid(ncid, 'met_stem_prof_mortC', varid))
      call check(nf90_get_var(ncid, varid, met_stem_prof_mortC,start=(/1, time_entry/)))
      met_stem_prof_mortC = met_stem_prof_mortC*sec_pr_hr  ![gC/m^2/h]


      
      !N mortality fluxes to metabolic: 
      call check(nf90_inq_varid(ncid, 'met_leaf_prof_mortN', varid))
      call check(nf90_get_var(ncid, varid, met_leaf_prof_mortN,start=(/1, time_entry/)))
      met_leaf_prof_mortN = met_leaf_prof_mortN*sec_pr_hr  ![gN/m^2/h]
      
      call check(nf90_inq_varid(ncid, 'met_froot_prof_mortN', varid))
      call check(nf90_get_var(ncid, varid, met_froot_prof_mortN,start=(/1, time_entry/)))
      met_froot_prof_mortN = met_froot_prof_mortN*sec_pr_hr  ![gN/m^2/h]

      call check(nf90_inq_varid(ncid, 'met_croot_prof_mortN', varid))
      call check(nf90_get_var(ncid, varid, met_croot_prof_mortN,start=(/1, time_entry/)))
      met_croot_prof_mortN = met_croot_prof_mortN*sec_pr_hr  ![gN/m^2/h]

      call check(nf90_inq_varid(ncid, 'met_stem_prof_mortN', varid))
      call check(nf90_get_var(ncid, varid, met_stem_prof_mortN,start=(/1, time_entry/)))
      met_stem_prof_mortN = met_stem_prof_mortN*sec_pr_hr  ![gN/m^2/h]
      
      met_mortalityC = met_leaf_prof_mortC*leaf_prof+met_froot_prof_mortC*froot_prof+met_croot_prof_mortC*croot_prof+met_stem_prof_mortC*stem_prof ![gC/m^3/h]
      met_mortalityN = met_leaf_prof_mortN*leaf_prof+met_froot_prof_mortN*froot_prof+met_croot_prof_mortN*croot_prof+met_stem_prof_mortN*stem_prof ![gN/m^3/h]

      
      !Split mortality fluxes:
      call check(nf90_inq_varid(ncid, 'split_leaf_prof_mortC', varid))
      call check(nf90_get_var(ncid, varid, split_leaf_prof_mortC,start=(/1, time_entry/)))
      split_leaf_prof_mortC = split_leaf_prof_mortC*sec_pr_hr  ![gC/m^2/h]
      
      call check(nf90_inq_varid(ncid, 'split_froot_prof_mortC', varid))
      call check(nf90_get_var(ncid, varid, split_froot_prof_mortC,start=(/1, time_entry/)))
      split_froot_prof_mortC = split_froot_prof_mortC*sec_pr_hr  ![gC/m^2/h]
      
      !N mortality fluxes to splitabolic: 
      call check(nf90_inq_varid(ncid, 'split_leaf_prof_mortN', varid))
      call check(nf90_get_var(ncid, varid, split_leaf_prof_mortN,start=(/1, time_entry/)))
      split_leaf_prof_mortN = split_leaf_prof_mortN*sec_pr_hr![gN/m^2/h]
      
      call check(nf90_inq_varid(ncid, 'split_froot_prof_mortN', varid))
      call check(nf90_get_var(ncid, varid, split_froot_prof_mortN,start=(/1, time_entry/)))
      split_froot_prof_mortN = split_froot_prof_mortN*sec_pr_hr  ![gN/m^2/h]

      split_mortalityC = split_leaf_prof_mortC*leaf_prof+split_froot_prof_mortC*froot_prof ![gC/m^3/h]
      split_mortalityN = split_leaf_prof_mortN*leaf_prof+split_froot_prof_mortN*froot_prof ![gN/m^3/h]
    end subroutine read_clm_mortality
                
    subroutine read_clm_model_input(ncid,time_entry,CLM_version, &
                                    LEAFN_TO_LITTER,FROOTN_TO_LITTER, NPP_MYC,NDEP_TO_SMINN, &
                                    LEAFC_TO_LITTER,FROOTC_TO_LITTER,mcdate,TSOI,SOILLIQ,SOILICE, &
                                    W_SCALAR,T_SCALAR,QDRAI,QOVER,H2OSOI_LIQ,C_CWD,N_CWD)
      !INPUT
      integer,intent(in)                        :: ncid
      integer,intent(in)                        :: time_entry
      character (len=3),intent(in)              :: CLM_version
        
      !OUTPUT 
      real(r8),intent(out)                      :: LEAFN_TO_LITTER !litterfall N from leaves  [gN/m^2/hour] (converted from [gN/m^2/s])      
      real(r8),intent(out)                      :: FROOTN_TO_LITTER !litterfall N from fine roots [gN/m^2/hour] (converted from [gN/m^2/s])
      real(r8),intent(out)                      :: NPP_MYC                    
      real(r8),intent(out)                      :: NDEP_TO_SMINN   !atmospheric N deposition to soil mineral N [gN/m^2/hour] (converted from [gN/m^2/s]) 
      real(r8),intent(out)                      :: LEAFC_TO_LITTER
      real(r8),intent(out)                      :: FROOTC_TO_LITTER
      integer ,intent(out)                      :: mcdate  !"Current date" , end of averaging interval
      real(r8),intent(out),dimension(nlevels)   :: TSOI    !Celcius
      real(r8),intent(out),dimension(nlevels)   :: H2OSOI_LIQ !kg/m2 TODO: Name this something else, to avoid confusion with CLM var H2OSOI!
      real(r8),intent(out),dimension(nlevels)   :: SOILLIQ !m3/m3 (converted from kg/m2) !TODO: Take in H2OSOI instead?
      real(r8),intent(out),dimension(nlevels)   :: SOILICE !m3/m3 (converted from kg/m2)
      real(r8),intent(out),dimension(nlevels)   :: W_SCALAR    
      real(r8),intent(out),dimension(nlevels)   :: T_SCALAR      
      real(r8),intent(out)                      :: QDRAI   !kgH2O/(m2 h) converted from kgH2O/(m2 s)
      real(r8),intent(out)                      :: QOVER   !kgH2O/(m2 h) converted from kgH2O/(m2 s)
      
      real(r8),intent(out)                      :: C_CWD(:) !gC/(m3 h)
      real(r8),intent(out)                      :: N_CWD(:) !gN/(m3 h)

      !Local
      integer                                   :: varid
      integer                                   :: i
      real(r8)                                  :: C_CWD2(nlevels)
      real(r8)                                  :: C_CWD3(nlevels)        
      real(r8)                                  :: N_CWD2(nlevels)
      real(r8)                                  :: N_CWD3(nlevels)       
      real(r8)                                  :: NPP_NACTIVE  !Mycorrhizal N uptake used C        [gC/m^2/hour] (converted from [gC/m^2/s]) 
      real(r8)                                  :: NPP_NNONMYC  !NONMycorrhizal N uptake used C        [gC/m^2/hour] (converted from [gC/m^2/s]) 
      
      print*, 'entering read_clm_model_input, time_entry=', time_entry
      print*, 'CLM_version=', CLM_version
      select case (CLM_version)
      case("old")
        ! C in Coarse Woody Debris
        call check(nf90_inq_varid(ncid, 'CWDC_TO_LITR2C_vr', varid))
        call check(nf90_get_var(ncid, varid, C_CWD2, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))
        call check(nf90_inq_varid(ncid, 'CWDC_TO_LITR3C_vr', varid))
        call check(nf90_get_var(ncid, varid, C_CWD3, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))
          
        ! N in Coarse Woody Debris
        call check(nf90_inq_varid(ncid, 'CWDN_TO_LITR2N_vr', varid))
        call check(nf90_get_var(ncid, varid, N_CWD2, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))      
        call check(nf90_inq_varid(ncid, 'CWDN_TO_LITR3N_vr', varid))
        call check(nf90_get_var(ncid, varid, N_CWD3, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))
      
      case("new")
        ! C in Coarse Woody Debris
        call check(nf90_inq_varid(ncid, 'CWD_C_TO_LIT_CEL_C_vr', varid))
        call check(nf90_get_var(ncid, varid, C_CWD2, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))
        call check(nf90_inq_varid(ncid, 'CWD_C_TO_LIT_LIG_C_vr', varid))
        call check(nf90_get_var(ncid, varid, C_CWD3, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))
            
        ! N in Coarse Woody Debris
        call check(nf90_inq_varid(ncid, 'CWD_N_TO_LIT_CEL_N_vr', varid))
        call check(nf90_get_var(ncid, varid, N_CWD2, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))      
        call check(nf90_inq_varid(ncid, 'CWD_N_TO_LIT_LIG_N_vr', varid))
        call check(nf90_get_var(ncid, varid, N_CWD3, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))       
      end select


      C_CWD=(C_CWD2+C_CWD3)*sec_pr_hr!gC/(m3 s) to gC/(m3 h)
      N_CWD=(N_CWD2+N_CWD3)*sec_pr_hr !gN/(m3 s) to gN/(m3 h)
      
      !C and N litter from leafs and fine roots:
      print*, 'reading LEAFN_TO_LITTER'
      call check(nf90_inq_varid(ncid, 'LEAFN_TO_LITTER', varid))
      call check(nf90_get_var(ncid, varid, LEAFN_TO_LITTER,start=(/1, time_entry/)))
      LEAFN_TO_LITTER = LEAFN_TO_LITTER*sec_pr_hr  ![gC/m^2/h]
      print*, 'reading FROOTN_TO_LITTER'
      call check(nf90_inq_varid(ncid, 'FROOTN_TO_LITTER', varid))
      call check(nf90_get_var(ncid, varid, FROOTN_TO_LITTER,start=(/1, time_entry/)))
      FROOTN_TO_LITTER = FROOTN_TO_LITTER*sec_pr_hr  ![gC/m^2/h]    
      
      call check(nf90_inq_varid(ncid, 'LEAFC_TO_LITTER', varid))
      call check(nf90_get_var(ncid, varid, LEAFC_TO_LITTER,start=(/1, time_entry/)))          
      LEAFC_TO_LITTER = LEAFC_TO_LITTER*sec_pr_hr  ![gC/m^2/h]
      call check(nf90_inq_varid(ncid, 'FROOTC_TO_LITTER', varid))
      call check(nf90_get_var(ncid, varid, FROOTC_TO_LITTER,start=(/1, time_entry/)))          
      FROOTC_TO_LITTER = FROOTC_TO_LITTER*sec_pr_hr  ![gC/m^2/h]

      !C payment for total payed for N (NPP_NACTIVE) and non-mycorrhizal (NPP_NNONMYC)
      call check(nf90_inq_varid(ncid, 'NPP_NACTIVE', varid))
      call check(nf90_get_var(ncid, varid, NPP_NACTIVE,start=(/1, time_entry/)))
      NPP_NACTIVE = NPP_NACTIVE*sec_pr_hr ![gC/m^2/h]

      call check(nf90_inq_varid(ncid, 'NPP_NNONMYC', varid))
      call check(nf90_get_var(ncid, varid, NPP_NNONMYC,start=(/1, time_entry/)))
      NPP_NNONMYC = NPP_NNONMYC*sec_pr_hr ![gC/m^2/h]
      
      !Part of total NPP_NACTIVE going to mycorrhizal associations
      NPP_MYC = NPP_NACTIVE - NPP_NNONMYC
      
      !N deposition:
      call check(nf90_inq_varid(ncid, 'NDEP_TO_SMINN', varid))
      call check(nf90_get_var(ncid, varid, NDEP_TO_SMINN,start=(/1, time_entry/)))
      NDEP_TO_SMINN = NDEP_TO_SMINN*sec_pr_hr ![gN/m^2/h]
      !Environmental variables:      
      call check(nf90_inq_varid(ncid, 'TSOI', varid))
      call check(nf90_get_var(ncid, varid, TSOI, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))
      TSOI = TSOI - 273.15 !UNIT CONVERSION: Kelvin to Celcius 
      
      call check(nf90_inq_varid(ncid, 'SOILLIQ', varid))
      call check(nf90_get_var(ncid, varid, SOILLIQ, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))

      call check(nf90_inq_varid(ncid, 'SOILICE', varid))
      call check(nf90_get_var(ncid, varid, SOILICE, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))

      call check(nf90_inq_varid(ncid, 'W_SCALAR', varid))
      call check(nf90_get_var(ncid, varid, W_SCALAR, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))
      
      call check(nf90_inq_varid(ncid, 'T_SCALAR', varid))
      call check(nf90_get_var(ncid, varid, T_SCALAR, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))

      call check(nf90_inq_varid(ncid, 'QDRAI', varid)) !mmH2O/s = kg H2O/(m2 s)
      call check(nf90_get_var(ncid, varid, QDRAI,start=(/1, time_entry/)))    
      QDRAI = QDRAI*sec_pr_hr !kgH2O/(m2 h)

      call check(nf90_inq_varid(ncid, 'QOVER', varid)) !mmH2O/s = kg H2O/(m2 s)
      call check(nf90_get_var(ncid, varid, QOVER,start=(/1, time_entry/)))    
      QOVER = QOVER*sec_pr_hr !kgH2O/(m2 h)
      
      H2OSOI_LIQ = SOILLIQ  
      do i = 1, nlevels ! Unit conversion
        SOILICE(i) = SOILICE(i)/(delta_z(i)*917) !kg/m2 to m3/m3 rho_ice=917kg/m3
        SOILLIQ(i) = SOILLIQ(i)/(delta_z(i)*1000) !kg/m2 to m3/m3 rho_liq=1000kg/m3
      end do          
    
      !CLM dates
      call check(nf90_inq_varid(ncid, 'mcdate', varid))
      call check(nf90_get_var(ncid, varid, mcdate, start=(/time_entry/)))
    end subroutine read_clm_model_input

    subroutine read_namelist(file_path, use_ROI, use_Sulman, use_ENZ, timestep) !Based on https://cerfacs.fr/coop/fortran-namelist-workedex
      !! Read some parmeters,  Here we use a namelist 
      !! but if you were to change the storage format (TOML,or home-made), 
      !! this signature would not change

      character(len=*),  intent(in)  :: file_path
      logical, intent(out) :: use_ROI
      logical, intent(out) :: use_Sulman
      logical, intent(out) :: use_ENZ
      real(r8), intent(out) :: timestep
      
      !integer, intent(out) :: type_
      integer                        :: file_unit, iostat

      ! Namelist definition===============================
      namelist /OPTIONS/ &
          use_ROI , &
          use_Sulman, &
          use_ENZ, &
          timestep
      use_ROI = .False.
      use_Sulman = .False.
      use_ENZ =  .False.
      timestep = 1
      ! Namelist definition===============================

      call open_inputfile(file_path, file_unit, iostat)
      if (iostat /= 0) then
          print*, "Opening of file failed"
          !! write here what to do if opening failed"
          return
      end if

      read (nml=OPTIONS, iostat=iostat, unit=file_unit)
      call close_inputfile(file_path, file_unit, iostat)
      if (iostat /= 0) then
          print*, "Closing of file failed"          
          !! write here what to do if reading failed"
          return
      end if
    end subroutine read_namelist

    subroutine open_inputfile(file_path, file_unit, iostat) !Based on https://cerfacs.fr/coop/fortran-namelist-workedex
      !! Check whether file exists, with consitent error message
      !! return the file unit
      character(len=*),  intent(in)  :: file_path
      integer,  intent(out) :: file_unit, iostat

      inquire (file=file_path, iostat=iostat)
      if (iostat /= 0) then
          write (stderr, '(3a)') 'Error: file "', trim(file_path), '" not found!'
      end if
      open (action='read', file=file_path, iostat=iostat, newunit=file_unit)
    end subroutine open_inputfile

    subroutine close_inputfile(file_path, file_unit, iostat) !Based on https://cerfacs.fr/coop/fortran-namelist-workedex
      !! Check the reading was OK
      !! return error line IF not
      !! close the unit
      character(len=*),  intent(in)  :: file_path
      character(len=1000) :: line
      integer,  intent(in) :: file_unit, iostat

      if (iostat /= 0) then
          write (stderr, '(2a)') 'Error reading file :"', trim(file_path)
          write (stderr, '(a, i0)') 'iostat was:"', iostat
          backspace(file_unit)
          read(file_unit,fmt='(A)') line
          write(stderr,'(A)') &
              'Invalid line : '//trim(line)
      end if
      close (file_unit)   
    end subroutine close_inputfile

end module readMod
