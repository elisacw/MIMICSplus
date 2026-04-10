module writeMod
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use paramMod !, only: pool_types, pool_types_N
  use netcdf
  use dispmodule, only: disp
  use initMod, only: nlevels
  implicit none
  private
  public :: check, create_netcdf,fill_netcdf,store_parameters,fluxes_netcdf, fill_MMK,&
            create_yearly_mean_netcdf,fill_yearly_netcdf,fill_monthly_netcdf,create_monthly_mean_netcdf
    
  integer :: grid_dimid, col_dimid, t_dimid, lev_dimid,mmk_dimid,fracid
  character (len=4), dimension(pool_types)     :: variables = &
  (/  "LITm", "LITs", "SAPb","SAPf", "EcM ", "AM  ", "SOMp", "SOMa", "SOMc" /)

  character (len=*), dimension(*), parameter ::  C_name_fluxes = &
  [character(len=11) ::"LITmSAPb","LITmSAPf","LITsSAPb","LITsSAPf", &
  "SAPbSOMp","SAPfSOMp", "SAPbSOMa","SAPfSOMa", "SAPbSOMc","SAPfSOMc", &
  "EcMSOMp ", "EcMSOMa ","EcMSOMc ","AMSOMp  ","AMSOMa  ","AMSOMc  ", &
  "SOMaSAPb","SOMaSAPf","SOMpSOMa","SOMcSOMa","PlantLITm", "PlantLITs", &
  "PlantEcM","PlantAM","PlantSOMc  ","PlantSOMp  ", &
  "EcMdecoSOMp","EcMdecoSOMc","EcMenz_prod"]
  character (len=*), dimension(*), parameter ::  N_name_fluxes = &
  [character(len=11) ::"LITmSAPb","LITmSAPf","LITsSAPb","LITsSAPf", &
  "SAPbSOMp","SAPfSOMp", "SAPbSOMa","SAPfSOMa", "SAPbSOMc","SAPfSOMc", &
  "EcMSOMp ", "EcMSOMa ","EcMSOMc ","AMSOMp  ","AMSOMa  ","AMSOMc  ", &
  "SOMaSAPb","SOMaSAPf","SOMpSOMa","SOMcSOMa","PlantLITm", "PlantLITs", &
  "PlantSOMp","PlantSOMc","EcMPlant","AMPlant", "Deposition",&
  "Leaching", "INEcM","INAM","SOMpEcM","SOMcEcM","nitrif_rate",'INSAPf','INSAPb']

  contains

    subroutine check(status)
      integer, intent ( in) :: status
      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 999
      end if
    end subroutine check

    subroutine create_netcdf(output_path,run_name, start_year)
      character (len = *),intent(in) :: run_name
      character (len = *),intent(in) :: output_path
      integer,            intent(in) :: start_year

      integer :: ncid, varid
      integer, parameter :: gridcell = 1, column = 1
      integer :: v
      character(len=10) :: st_year_c
      print *,(output_path//trim(run_name)//".nc")
      call check(nf90_create(output_path//trim(run_name)//".nc",NF90_NETCDF4,ncid))


      call check(nf90_def_dim(ncid, "time", nf90_unlimited, t_dimid))
      call check(nf90_def_dim(ncid, "gridcell", gridcell, grid_dimid))
      call check(nf90_def_dim(ncid, "column", column, col_dimid))
      call check(nf90_def_dim(ncid, "levsoi", nlevels, lev_dimid))
      call check(nf90_def_dim(ncid, "NoMMKeqs", MM_eqs, mmk_dimid))
      call check(nf90_def_dim(ncid,"SAPpools",size(fPHYS),fracid))
      do v = 1, size(variables)
        call check(nf90_def_var(ncid, trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gCm^-3"))
        call check(nf90_def_var(ncid, "vert_change"//trim(variables(v)), NF90_FLOAT, (/t_dimid, lev_dimid/), varid))
        call check(nf90_put_att(ncid,varid,"unit","gCm^-3h^-1"))
        call check(nf90_def_var(ncid, "N_"//trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
        call check(nf90_def_var(ncid, "N_vert_change"//trim(variables(v)), NF90_FLOAT, (/t_dimid, lev_dimid/), varid))
        call check(nf90_put_att(ncid,varid,"unit","gNm^-3h^-1"))
        
        
      end do
      call check(nf90_def_var(ncid, "NH4_sol", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
        call check(nf90_def_var(ncid, "NH4_sorp", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
        call check(nf90_def_var(ncid, "NH4_sorp_eq", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
        call check(nf90_def_var(ncid, "NO3", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
        call check(nf90_def_var(ncid, "N_SMIN", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
        call check(nf90_def_var(ncid,"HR_sum", NF90_FLOAT, (/t_dimid /), varid ))
      call check(nf90_def_var(ncid,"HR_flux", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","gCm^-3h^-1"))
      call check(nf90_def_var(ncid,"HRb", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","gCm^-3h^-1"))
      call check(nf90_def_var(ncid,"HRf", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","gCm^-3h^-1"))
      call check(nf90_def_var(ncid,"HRe", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","gCm^-3h^-1"))
      call check(nf90_def_var(ncid,"HRa", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","gCm^-3h^-1"))
      call check(nf90_def_var(ncid, "NH4sol_vert_change",NF90_FLOAT,(/t_dimid,lev_dimid/), varid))
      call check(nf90_put_att(ncid,varid,"unit","gNm^-3h^-1"))      
      call check(nf90_def_var(ncid, "NH4sorp_vert_change",NF90_FLOAT,(/t_dimid,lev_dimid/), varid))
      call check(nf90_put_att(ncid,varid,"unit","gNm^-3h^-1"))
      call check(nf90_def_var(ncid, "NO3_vert_change",NF90_FLOAT,(/t_dimid,lev_dimid/), varid))
      call check(nf90_put_att(ncid,varid,"unit","gNm^-3h^-1"))


      call check(nf90_def_var(ncid,"r_myc", NF90_FLOAT, (/t_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))
      call check(nf90_def_var(ncid,"f_ecm", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))
      call check(nf90_def_var(ncid,"f_am", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))
      call check(nf90_def_var(ncid,"CUEb", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))
      call check(nf90_def_var(ncid,"CUEf", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))
      call check(nf90_def_var(ncid,"CUE_ecm", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))
      call check(nf90_def_var(ncid,"CUE_am", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))
      call check(nf90_def_var(ncid,"ROI_ecm", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","gN/gC"))
      call check(nf90_def_var(ncid,"ROI_am", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","gN/gC"))
      call check(nf90_def_var(ncid,"enz_frac", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))
      call check(nf90_def_var(ncid, "Temp", NF90_FLOAT, (/t_dimid, lev_dimid/),varid))
      call check(nf90_put_att(ncid,varid,"unit","deg C"))
      call check(nf90_def_var(ncid, "Moisture", NF90_FLOAT, (/t_dimid, lev_dimid/),varid))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))
      !call check(nf90_def_var(ncid, "N_changeinorganic", NF90_FLOAT,(/t_dimid, lev_dimid/), varid))
      call check(nf90_def_var(ncid, "N_InPlant", NF90_FLOAT,(/t_dimid, lev_dimid/), varid))
      call check(nf90_put_att(ncid,varid,"unit","gNm^-3h^-1"))
      
      call check(nf90_def_var(ncid, "Km", NF90_FLOAT, (/mmk_dimid,t_dimid, lev_dimid/),varid))
      call check(nf90_put_att(ncid,varid,"unit","gCm^-3"))
      call check(nf90_def_var(ncid, "Vmax", NF90_FLOAT, (/mmk_dimid,t_dimid, lev_dimid/),varid))
      call check(nf90_put_att(ncid,varid,"unit","mgC/((mgSAP)h)"))

      do v = 1, size(C_name_fluxes)
        call check(nf90_def_var(ncid, "C_"//trim(C_name_fluxes(v)), NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gCm^-3h^-1"))

      end do

      do v = 1, size(N_name_fluxes)
         call check(nf90_def_var(ncid, "N_"//trim(N_name_fluxes(v)), NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
         call check(nf90_put_att(ncid,varid,"unit","gNm^-3h^-1"))

      end do
      
      call check(nf90_def_var(ncid, "mcdate", NF90_INT,(/t_dimid/), varid))
      call check(nf90_put_att(ncid,varid,"unit","yyyymmdd"))
      call check(nf90_def_var(ncid, "time", NF90_FLOAT, (/t_dimid /), varid))
      write(st_year_c,"(I4)") start_year

      call check(nf90_put_att(ncid,varid,"units","days since " // st_year_c //"-01-01"))
      write(*,*) "days since ",st_year_c,"-01-01"
      call check(nf90_def_var(ncid, "month", NF90_INT, (/t_dimid /), varid))
      call check(nf90_def_var(ncid, "f_met", NF90_FLOAT,(/t_dimid /),varid))
      call check(nf90_put_att(ncid,varid,"unit","unitless"))

      call check(nf90_enddef(ncid))

      call check( nf90_close(ncid) )
    end subroutine create_netcdf

    subroutine fill_MMK(ncid,time,write_hour,depth,K_m,V_max)
      integer,intent(in)               :: ncid 
      integer, intent(in)              :: time
      integer, intent(in)              :: depth
      integer, intent(in)              :: write_hour
      
      real(r8),dimension(MM_eqs), intent(in)       :: K_m
      real(r8),dimension(MM_eqs), intent(in)       :: V_max
      
      integer                          :: varid, timestep 
      
      call get_timestep(time, write_hour, timestep)
      
      call check(nf90_inq_varid(ncid, "Km", varid))
      call check(nf90_put_var(ncid, varid, K_m, start = (/ 1,timestep, depth /)))
      
      call check(nf90_inq_varid(ncid, "Vmax", varid))
      call check(nf90_put_var(ncid, varid, V_max, start = (/ 1,timestep, depth /)))
    end subroutine fill_MMK
      
    subroutine fill_netcdf(ncid, time, pool_matrix, Npool_matrix,inorganic_N_matrix, &
      mcdate,HR_sum, HR_flux, HRb,HRf,HRe,HRa,vert_sum,Nvert_sum,Ninorg_vert_sum, write_hour,month, &
      TSOIL, MOIST,CUE_bacteria,CUE_fungi,CUE_ecm,CUE_am,ROI_EcM,ROI_AM,enz_frac,f_met,f_alloc, NH4_eq)
      !INPUT:
      integer,intent(in)               :: ncid 
      integer, intent(in)              :: time
      real(r8), intent(in)             :: pool_matrix(nlevels,pool_types), Npool_matrix(nlevels,pool_types_N)   ! For storing pool concentrations [gC/m3]
      real(r8), intent(in)             :: inorganic_N_matrix(nlevels,inorg_N_pools)
      real(r8), intent(in)             :: vert_sum(nlevels,pool_types)
      real(r8), intent(in)             :: Nvert_sum(nlevels,pool_types_N)
      real(r8), intent(in)             :: Ninorg_vert_sum(nlevels,inorg_N_pools)
      integer, intent(in)              :: mcdate
      integer, intent(in)              :: write_hour
      integer, intent(in)              :: month
      real(r8),intent(in)              :: HR_sum
      real(r8),dimension(nlevels,2),intent(in)              :: f_alloc
      real(r8),dimension(nlevels), intent(in)       :: HR_flux,HRb,HRf
      real(r8),dimension(nlevels), intent(in)       :: HRe,HRa
      real(r8),dimension(nlevels), intent(in)       :: ROI_EcM
      real(r8),dimension(nlevels), intent(in)       :: ROI_AM  
      real(r8),dimension(nlevels), intent(in)       :: enz_frac          
      real(r8),                    intent(in)       :: f_met
      real(r8),dimension(nlevels), intent(in)       :: TSOIL, MOIST  
      real(r8),dimension(nlevels), intent(in)       :: CUE_bacteria,CUE_fungi, CUE_ecm,CUE_am
      real(r8),dimension(nlevels), intent(in)       :: NH4_eq
      
      
          
      !OUTPUT:
      
      !LOCAL:
      integer                          :: i , j !for looping
      integer                          :: varid, timestep 
      real(r8)                         :: N_SMIN
      real(r8)                         :: time_float

      call get_timestep(time, write_hour, timestep)

      time_float = real(time)

      call check(nf90_inq_varid(ncid, "time", varid))
      call check(nf90_put_var(ncid, varid, time/24.0_r8, start = (/ timestep /)))
      call check(nf90_inq_varid(ncid, "mcdate", varid))
      call check(nf90_put_var(ncid, varid, mcdate, start = (/ timestep /)))    
      call check(nf90_inq_varid(ncid, "f_met", varid))
      call check(nf90_put_var(ncid, varid, f_met, start = (/ timestep /)))    
      call check(nf90_inq_varid(ncid, "month", varid))
      call check(nf90_put_var(ncid, varid, month, start = (/ timestep /)))
      call check(nf90_inq_varid(ncid, "HR_sum", varid))
      call check(nf90_put_var(ncid, varid, HR_sum, start = (/ timestep /)))
      call check(nf90_inq_varid(ncid, "r_myc", varid))
      call check(nf90_put_var(ncid, varid, r_myc, start = (/ timestep /)))
      
      do j=1,nlevels
        
        call check(nf90_inq_varid(ncid, "Temp", varid))
        call check(nf90_put_var(ncid, varid, TSOIL(j), start = (/ timestep, j /)))

        call check(nf90_inq_varid(ncid, "NH4sol_vert_change", varid))
        call check(nf90_put_var(ncid, varid, Ninorg_vert_sum(j,1), start = (/ timestep, j /)))
        call check(nf90_inq_varid(ncid, "NH4sorp_vert_change", varid))
        call check(nf90_put_var(ncid, varid, Ninorg_vert_sum(j,2), start = (/ timestep, j /)))
        call check(nf90_inq_varid(ncid, "NO3_vert_change", varid))
        call check(nf90_put_var(ncid, varid, Ninorg_vert_sum(j,3), start = (/ timestep, j /)))
        
        call check_and_write(ncid, "Moisture", MOIST(j),timestep,j)

        call check_and_write(ncid, "HR_flux", HR_flux(j), timestep,j)
        
        call check_and_write(ncid, "HRb", HRb(j),timestep,j)
              
        call check_and_write(ncid, "HRf", HRf(j),timestep,j)

        call check_and_write(ncid, "HRe", HRe(j),timestep,j)
        
        call check_and_write(ncid, "HRa",  HRa(j), timestep,j)
                       
        call check_and_write(ncid, "ROI_ecm", ROI_EcM(j),timestep,j)

        call check_and_write(ncid, "ROI_am",ROI_AM(j),  timestep,j)
        
        call check_and_write(ncid, "enz_frac",enz_frac(j),  timestep,j)
        
        call check_and_write(ncid, "f_ecm",f_alloc(j,1),  timestep,j)
                
        call check_and_write(ncid, "f_am", f_alloc(j,2),timestep,j)
    
        call check_and_write(ncid, "CUEb", CUE_bacteria(j), timestep,j)
        
        call check_and_write(ncid, "CUEf",CUE_fungi(j), timestep,j)
        
        call check_and_write(ncid, "CUE_ecm",CUE_ecm(j),  timestep,j)
        
        call check_and_write(ncid, "CUE_am",  CUE_am(j),timestep,j)

        call check_and_write(ncid, "NH4_sorp_eq", NH4_eq(j), timestep,j)
        
        call check_and_write(ncid, "NH4_sol", inorganic_N_matrix(j,1),timestep,j)
        call check_and_write(ncid, "NH4_sorp",inorganic_N_matrix(j,2),  timestep,j)
        call check_and_write(ncid, "NO3",  inorganic_N_matrix(j,3),timestep,j)
      
        N_SMIN = inorganic_N_matrix(j,1)+inorganic_N_matrix(j,2)+inorganic_N_matrix(j,3)
        call check_and_write(ncid, "N_SMIN",N_SMIN, timestep,j)

        do i = 1,pool_types
          !C:
          call check_and_write(ncid,trim(variables(i)),pool_matrix(j,i),timestep,j)

          !N:
          call check_and_write(ncid,"N_"//trim(variables(i)),Npool_matrix(j,i),timestep,j)

          !C change:
          call check(nf90_inq_varid(ncid,"vert_change"//trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, vert_sum(j,i), start = (/ timestep, j /)))    
          
          !N change:        
          call check(nf90_inq_varid(ncid,"N_vert_change"//trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, Nvert_sum(j,i), start = (/ timestep, j /)))      
          
        end do !pool_types
      end do ! levels
    end subroutine fill_netcdf

    subroutine store_parameters(ncid, soil_depth,desorp,dz)
      !IN:
      integer,intent(in):: ncid
      real(r8),intent(in):: soil_depth
      real(r8),intent(in):: desorp
      real(r8),dimension(:) :: dz
    
      !Local:
      integer :: clayID, desorpID, k_sapsomID, depthID, fphysID, fchemID, favailID,dzID

      call check(nf90_def_var(ncid, "f_clay", NF90_FLOAT, clayID))
      call check(nf90_put_att(ncid,clayID,"unit","unitless"))
      call check(nf90_def_var(ncid, "desorp", NF90_FLOAT, desorpID))
      call check(nf90_put_att(ncid,desorpID,"unit","h^-1"))
      call check(nf90_def_var(ncid, "k_sapsom", NF90_FLOAT,fracid,k_sapsomID))
      call check(nf90_put_att(ncid,k_sapsomID,"unit","h^-1"))
      call check(nf90_def_var(ncid, "f_phys", NF90_FLOAT,fracid,fphysID))
      call check(nf90_put_att(ncid,fphysID,"unit","unitless"))
      call check(nf90_def_var(ncid, "f_avail", NF90_FLOAT,fracid,favailID))
      call check(nf90_put_att(ncid,favailID,"unit","unitless"))
      call check(nf90_def_var(ncid, "f_chem", NF90_FLOAT,fracid,fchemID))
      call check(nf90_put_att(ncid,fchemID,"unit","unitless"))
      call check(nf90_def_var(ncid, "depth", NF90_FLOAT,depthID))
      call check(nf90_put_att(ncid,depthID,"unit","m"))
      call check(nf90_def_var(ncid, "layer_thickness", NF90_FLOAT,(/lev_dimid/),dzID))
      call check(nf90_put_att(ncid,dzID,"unit","m"))

      call check(nf90_enddef(ncid))

      call check(nf90_put_var(ncid, clayID, fCLAY))
      call check(nf90_put_var(ncid, desorpID, desorp))
      call check(nf90_put_var(ncid, fphysID, fPHYS))
      call check(nf90_put_var(ncid, fchemID, fCHEM))
      call check(nf90_put_var(ncid, favailID, fAVAIL))
      call check(nf90_put_var(ncid, depthID, soil_depth))
      call check(nf90_put_var(ncid, k_sapsomID, k_sapsom))
      call check(nf90_put_var(ncid,dzID,dz))
    end subroutine store_parameters

    subroutine get_timestep(time, write_hour, timestep)
      integer, intent(in) :: time
      integer, intent(in) :: write_hour !hours between every output-writing.
      
      integer, intent(out):: timestep

      if (write_hour == 1) then
        timestep = time/write_hour
      else
        if (time == 1) then
          timestep = 1
        else
          timestep = time/(write_hour)+1 
        end if
      end if
    end subroutine get_timestep

    subroutine fluxes_netcdf(ncid,time, write_hour, depth_level, &
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
      integer,intent(in)  :: ncid
      integer, intent(in) :: time
      integer, intent(in) :: write_hour
      integer, intent(in) :: depth_level
      real(r8),intent(in) :: C_LITmSAPb, C_LITsSAPb, C_EcMSOMp, C_EcMSOMa, C_EcMSOMc, C_AMSOMp, &
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
                            N_INSAPb,N_INSAPf,N_PlantSOMp, C_EcMenz_prod, N_PlantLITs, N_PlantLITm,  N_PlantSOMc
      !Local:
      integer :: timestep
      integer :: varid
     !---------------------------------------------
      call get_timestep(time, write_hour, timestep)

      !C fluxes
      call check_and_write(ncid, "C_PlantLITm", C_PlantLITm,timestep,depth_level)
      call check_and_write(ncid, "C_PlantLITs",C_PlantLITs ,timestep,depth_level)
      call check_and_write(ncid, "C_PlantEcM", C_PlantEcM,timestep,depth_level)
      call check_and_write(ncid, "C_PlantAM", C_PlantAM,timestep,depth_level)
      call check_and_write(ncid, "C_PlantSOMc", C_PlantSOMc,timestep,depth_level)
      call check_and_write(ncid, "C_PlantSOMp",C_PlantSOMp ,timestep,depth_level)
      call check_and_write(ncid, "C_LITmSAPb",C_LITmSAPb ,timestep,depth_level)
      call check_and_write(ncid, "C_LITsSAPb",C_LITsSAPb ,timestep,depth_level)
      call check_and_write(ncid, "C_LITmSAPf", C_LITmSAPf,timestep,depth_level)
      call check_and_write(ncid, "C_LITsSAPf", C_LITsSAPf,timestep,depth_level)
      call check_and_write(ncid, "C_SOMpSOMa", C_SOMpSOMa,timestep,depth_level)
      call check_and_write(ncid, "C_SOMcSOMa", C_SOMcSOMa,timestep,depth_level)
      call check_and_write(ncid, "C_SAPbSOMa", C_SAPbSOMa,timestep,depth_level)
      call check_and_write(ncid, "C_SAPfSOMa", C_SAPfSOMa,timestep,depth_level)
      call check_and_write(ncid, "C_EcMSOMa",C_EcMSOMa ,timestep,depth_level)
      call check_and_write(ncid, "C_EcMSOMc",C_EcMSOMc ,timestep,depth_level)
      call check_and_write(ncid, "C_EcMSOMp", C_EcMSOMp,timestep,depth_level)
      call check_and_write(ncid, "C_AMSOMp", C_AMSOMp,timestep,depth_level)
      call check_and_write(ncid, "C_AMSOMa", C_AMSOMa,timestep,depth_level)
      call check_and_write(ncid, "C_AMSOMc", C_AMSOMc,timestep,depth_level)
      call check_and_write(ncid, "C_SOMaSAPb",C_SOMaSAPb ,timestep,depth_level)
      call check_and_write(ncid, "C_SOMaSAPf", C_SOMaSAPf,timestep,depth_level)
      call check_and_write(ncid, "C_SAPbSOMp",C_SAPbSOMp,timestep,depth_level)
      call check_and_write(ncid, "C_SAPfSOMp", C_SAPfSOMp,timestep,depth_level)
      call check_and_write(ncid, "C_SAPbSOMc",C_SAPbSOMc ,timestep,depth_level)
      call check_and_write(ncid, "C_SAPfSOMc", C_SAPfSOMc,timestep,depth_level)
      call check_and_write(ncid, "C_EcMdecoSOMp",C_EcMdecompSOMp ,timestep,depth_level)
      call check_and_write(ncid, "C_EcMdecoSOMc",C_EcMdecompSOMc ,timestep,depth_level)
      call check_and_write(ncid, "C_EcMenz_prod",C_EcMenz_prod ,timestep,depth_level)
      !N fluxes:
      call check_and_write(ncid, "N_LITmSAPb", N_LITmSAPb,timestep,depth_level)
      call check_and_write(ncid, "N_LITsSAPb", N_LITsSAPb, timestep,depth_level)
      call check_and_write(ncid, "N_LITmSAPf", N_LITmSAPf, timestep,depth_level)
      call check_and_write(ncid, "N_LITsSAPf", N_LITsSAPf ,timestep,depth_level)
      call check_and_write(ncid, "N_SOMpSOMa", N_SOMpSOMa ,timestep,depth_level)
      call check_and_write(ncid, "N_SOMcSOMa", N_SOMcSOMa ,timestep,depth_level)
      call check_and_write(ncid, "N_SAPbSOMa", N_SAPbSOMa,timestep,depth_level)
      call check_and_write(ncid, "N_SAPfSOMa",N_SAPfSOMa ,timestep,depth_level)
      call check_and_write(ncid, "N_EcMSOMa", N_EcMSOMa,timestep,depth_level)
      call check_and_write(ncid, "N_EcMSOMc",N_EcMSOMc ,timestep,depth_level)
      call check_and_write(ncid, "N_EcMSOMp", N_EcMSOMp,timestep,depth_level)
      call check_and_write(ncid, "N_AMSOMa", N_AMSOMa,timestep,depth_level)
      call check_and_write(ncid, "N_AMSOMc", N_AMSOMc,timestep,depth_level)
      call check_and_write(ncid, "N_AMSOMp", N_AMSOMp,timestep,depth_level)
      call check_and_write(ncid, "N_SOMaSAPb",N_SOMaSAPb ,timestep,depth_level)
      call check_and_write(ncid, "N_SOMaSAPf", N_SOMaSAPf,timestep,depth_level)
      call check_and_write(ncid, "N_SAPbSOMp", N_SAPbSOMp,timestep,depth_level)
      call check_and_write(ncid, "N_SAPfSOMp", N_SAPfSOMp,timestep,depth_level)
      call check_and_write(ncid, "N_SAPbSOMc", N_SAPbSOMc,timestep,depth_level)
      call check_and_write(ncid, "N_SAPfSOMc",N_SAPfSOMc ,timestep,depth_level)
      call check(nf90_inq_varid(ncid, "N_INSAPf", varid))
      call check(nf90_put_var(ncid, varid, N_INSAPf, start = (/ timestep, depth_level /)))     
      call check(nf90_inq_varid(ncid, "N_INSAPb", varid))
      call check(nf90_put_var(ncid, varid, N_INSAPb, start = (/ timestep, depth_level /)))          
      call check_and_write(ncid, "N_INAM", N_INAM,timestep,depth_level)
      call check_and_write(ncid, "N_INEcM", N_INEcM,timestep,depth_level)
      call check_and_write(ncid, "N_Deposition",Deposition ,timestep,depth_level)
      call check_and_write(ncid, "N_Leaching",Leaching ,timestep,depth_level)
      call check_and_write(ncid, "N_PlantLITm",N_PlantLITm ,timestep,depth_level)
      call check_and_write(ncid, "N_PlantLITs", N_PlantLITs,timestep,depth_level)
      call check_and_write(ncid, "N_PlantSOMp",N_PlantSOMp ,timestep,depth_level)
      call check_and_write(ncid, "N_PlantSOMc",N_PlantSOMc ,timestep,depth_level)
      call check_and_write(ncid, "N_EcMPlant", N_EcMPlant,timestep,depth_level)
      call check_and_write(ncid, "N_AMPlant", N_AMPlant,timestep,depth_level)
      call check_and_write(ncid, "N_InPlant",N_InPlant ,timestep,depth_level)
      call check_and_write(ncid, "N_SOMpEcM", N_SOMpEcM,timestep,depth_level)
      call check_and_write(ncid, "N_SOMcEcM",N_SOMcEcM ,timestep,depth_level)
      call check_and_write(ncid, "N_nitrif_rate", nitrif_rate,timestep,depth_level)

    end subroutine fluxes_netcdf
    
    subroutine check_and_write(ncid,name_of_var,var,time,depth) !NB; do not use on values that can be negative!
      integer, intent(in)      :: ncid
      character(len=*),intent(in) :: name_of_var
      real(r8), intent(in)     :: var
      integer, intent(in)      :: time
      integer, intent(in)      :: depth
      
      !Local:
      integer  :: varid 
      real(r8) :: used_var      

      if ( var < 1e-30 ) then
        used_var = 0._r8 
      else
        used_var = var
      end if
      
      call check(nf90_inq_varid(ncid, name_of_var, varid))
      call check(nf90_put_var(ncid, varid, used_var, start = (/ time, depth /)))
      
    end subroutine check_and_write

    subroutine create_yearly_mean_netcdf(output_path,run_name)
      character (len = *), intent(in):: run_name
      character (len = *), intent(in) :: output_path

      integer :: ncid, varid
      integer, parameter :: gridcell = 1, column = 1
      integer :: v
      call check(nf90_create(trim(output_path)//trim(run_name)//"_yearly_mean.nc",NF90_NETCDF4,ncid))
      call check(nf90_def_dim(ncid, "time", nf90_unlimited, t_dimid))
      call check(nf90_def_dim(ncid, "gridcell", gridcell, grid_dimid))
      call check(nf90_def_dim(ncid, "column", column, col_dimid))
      call check(nf90_def_dim(ncid, "levsoi", nlevels, lev_dimid))

      do v = 1, size(variables)
        call check(nf90_def_var(ncid, trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gCm^-3"))
        call check(nf90_def_var(ncid, "N_"//trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
      end do
      
      call check(nf90_def_var(ncid, "N_NH4_sol", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
      call check(nf90_def_var(ncid, "N_NH4_sorp", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
      call check(nf90_def_var(ncid, "N_NO3", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      call check(nf90_put_att(ncid,varid,"unit","gNm^-3"))
      call check(nf90_def_var(ncid, "HR_flux", NF90_FLOAT, (/t_dimid /), varid))
      call check(nf90_put_att(ncid,varid,"unit","gCm^-3h-1"))
      
      call check(nf90_def_var(ncid, "year_since_start", NF90_FLOAT, (/t_dimid /), varid))

      call check(nf90_def_var(ncid, "layer_thickness", NF90_FLOAT,(/lev_dimid/),varid))
      call check(nf90_put_att(ncid,varid,"unit","m"))
      call check(nf90_enddef(ncid))

      call check( nf90_close(ncid) )
    end subroutine create_yearly_mean_netcdf

    subroutine create_monthly_mean_netcdf(output_path,run_name)
      character (len = *), intent(in):: run_name
      character (len = *),intent(in) :: output_path

      integer :: ncid, varid
      integer, parameter :: gridcell = 1, column = 1
      integer :: v
      call check(nf90_create(output_path//trim(run_name)//"_monthly_mean.nc",NF90_NETCDF4,ncid))

      call check(nf90_def_dim(ncid, "time", nf90_unlimited, t_dimid))
      call check(nf90_def_dim(ncid, "gridcell", gridcell, grid_dimid))
      call check(nf90_def_dim(ncid, "column", column, col_dimid))
      call check(nf90_def_dim(ncid, "levsoi", nlevels, lev_dimid))

      do v = 1, size(variables)
        call check(nf90_def_var(ncid, trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "N_"//trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
      end do

      call check(nf90_def_var(ncid, "N_NH4_sol", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      call check(nf90_def_var(ncid, "N_NH4_sorp", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      call check(nf90_def_var(ncid, "N_NO3", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      
      call check(nf90_def_var(ncid, "months_since_start", NF90_FLOAT, (/t_dimid /), varid))
      call check(nf90_enddef(ncid))
      call check(nf90_def_var(ncid, "year_since_start", NF90_FLOAT, (/t_dimid /), varid))
      call check(nf90_enddef(ncid))
      call check(nf90_def_var(ncid, "month_in_year", NF90_FLOAT, (/t_dimid /), varid))
      call check(nf90_enddef(ncid))

      call check( nf90_close(ncid) )
    end subroutine create_monthly_mean_netcdf

    subroutine fill_yearly_netcdf(output_path,run_name, year, Cpool_yearly, Npool_yearly, Ninorg_pool_yearly,HR_flux_yearly,dz)
      !INPUT
      character (len = *),intent(in):: run_name
      character (len = *),intent(in) :: output_path
      integer,intent(in)            :: year
      real(r8),intent(in)           :: HR_flux_yearly  
      real(r8), intent(in)          :: Cpool_yearly(nlevels,pool_types)  ! For storing C pool sizes [gC/m3]
      real(r8),intent(in)           :: Npool_yearly(nlevels,pool_types_N)  
      real(r8),intent(in)           :: Ninorg_pool_yearly(nlevels,inorg_N_pools)  
      real(r8),dimension(:),intent(in):: dz
      
      !LOCAL
      integer :: i,j,varid,ncid
      real(r8) :: value

      call check(nf90_open(output_path//trim(run_name)//"_yearly_mean.nc", nf90_write, ncid))
      call check(nf90_inq_varid(ncid, "year_since_start", varid))
      call check(nf90_put_var(ncid, varid, year , start = (/ year /)))

      call check(nf90_inq_varid(ncid, "layer_thickness", varid))
      call check(nf90_put_var(ncid,varid,dz))


      call check(nf90_inq_varid(ncid, "HR_flux", varid))
      call check(nf90_put_var(ncid, varid, HR_flux_yearly, start = (/year/)))

      do j=1,nlevels
        call check(nf90_inq_varid(ncid, "N_NH4_sol", varid))
        call check(nf90_put_var(ncid, varid, Ninorg_pool_yearly(j,1), start = (/year, j/)))
        
        call check(nf90_inq_varid(ncid, "N_NH4_sorp", varid))
        call check(nf90_put_var(ncid, varid, Ninorg_pool_yearly(j,2), start = (/year, j/)))
        
        call check(nf90_inq_varid(ncid, "N_NO3", varid))
        call check(nf90_put_var(ncid, varid, Ninorg_pool_yearly(j,3), start = (/year, j/)))

                
        do i = 1,pool_types
          if ( Cpool_yearly(j,i) < epsilon(Cpool_yearly) ) then
            value = 0.0_r8
          else
            value = Cpool_yearly(j,i)
          end if
          !C:
          call check(nf90_inq_varid(ncid, trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, value, start = (/ year, j /)))
          !N:
          call check(nf90_inq_varid(ncid, "N_"//trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, Npool_yearly(j,i), start = (/ year, j /)))

        end do !pool_types
      end do ! levels
      call check(nf90_close(ncid))
    end subroutine fill_yearly_netcdf

    subroutine fill_monthly_netcdf(output_path,run_name, year,month,months_since_start, Cpool_monthly, Npool_monthly, Ninorg_pool_monthly) !TODO: monthly HR and climate variables (if needed?)
      !INPUT
      character (len = *),intent(in):: run_name
      character (len = *),intent(in) :: output_path

      integer,intent(in)            :: year
      integer,intent(in)            :: month      
      integer,intent(in)            :: months_since_start            
      real(r8), intent(in)          :: Cpool_monthly(nlevels,pool_types)  ! For storing C pool sizes [gC/m3]
      real(r8),intent(in)           :: Npool_monthly(nlevels,pool_types_N)  
      real(r8),intent(in)           :: Ninorg_pool_monthly(nlevels,inorg_N_pools)  
      
      
      !LOCAL
      integer :: i,j,varid,ncid
      real(r8) :: value
      
      call check(nf90_open(output_path//trim(run_name)//"_monthly_mean.nc", nf90_write, ncid))
      call check(nf90_inq_varid(ncid, "months_since_start", varid))
      call check(nf90_put_var(ncid, varid, months_since_start , start = (/ months_since_start /)))
      call check(nf90_inq_varid(ncid, "month_in_year", varid))
      call check(nf90_put_var(ncid, varid, month , start = (/ months_since_start /)))
      call check(nf90_inq_varid(ncid, "year_since_start", varid))
      call check(nf90_put_var(ncid, varid, year , start = (/ months_since_start /)))
      do j=1,nlevels
        call check(nf90_inq_varid(ncid, "N_NH4_sol", varid))
        call check(nf90_put_var(ncid, varid, Ninorg_pool_monthly(j,1), start = (/months_since_start, j/)))
        
        call check(nf90_inq_varid(ncid, "N_NH4_sorp", varid))
        call check(nf90_put_var(ncid, varid, Ninorg_pool_monthly(j,2), start = (/months_since_start, j/)))
        
        call check(nf90_inq_varid(ncid, "N_NO3", varid))
        call check(nf90_put_var(ncid, varid, Ninorg_pool_monthly(j,3), start = (/months_since_start, j/)))
                
        do i = 1,pool_types
          if ( Cpool_monthly(j,i) < epsilon(Cpool_monthly) ) then
            value = 0.0_r8
          else
            value = Cpool_monthly(j,i)
          end if
          !C:
          call check(nf90_inq_varid(ncid, trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, value, start = (/ months_since_start, j /)))
          !N:
          call check(nf90_inq_varid(ncid, "N_"//trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, Npool_monthly(j,i), start = (/ months_since_start, j /)))

        end do !pool_types
      end do ! levels
      call check(nf90_close(ncid))
    end subroutine fill_monthly_netcdf
    
end module writeMod
