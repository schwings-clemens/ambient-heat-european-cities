get_HWMId_vCities <- function(var_name, fname_MOD, fname_out, ref_years, app_years, dir_functions) {

    #Load libraries
    .libPaths(c("/uio/kant/div-cicero-u1/clems/R/x86_64-pc-linux-gnu-library/3.6", .libPaths()))
    library(ncdf4)
    library(ncdf4.helpers)
    library(PCICt)
    library(doParallel)

    #Load HWMId function
    source(paste0(dir_functions, 'calc_HWMId.r'))

    #Open data sets
    nc_MOD <- nc_open(fname_MOD)

    #Read data
    data_MOD <- ncvar_get(nc_MOD, var_name)

    #Read variables and get information about calendar and units
    time  <- ncvar_get(nc_MOD, "time")
    cities <- ncvar_get(nc_MOD, "city")
    time_vec <- nc.get.time.series(nc_MOD, time.dim.name = "time")
    time_vec <- as.POSIXlt(time_vec)
    calendar <- nc_MOD$dim$time$calendar
    tunits  <- ncatt_get(nc_MOD, "time", "units")
    units   <- ncatt_get(nc_MOD, var_name, "units")

    #Close NetCDFs
    nc_close(nc_MOD)

    #Create time vectors for HWMId function
    time_in <- NULL
    time_in$DOY   <- time_vec$yday + 1 #Correct for start at 0
    time_in$day   <- as.numeric(format(time_vec, "%d"))
    time_in$month <- as.numeric(format(time_vec, "%m"))
    time_in$year  <- as.numeric(format(time_vec, "%Y"))
    time_in$calendar <- calendar

    #Create indicex for all grid points, on which to carry out the analysis
    N_city = length(cities)

#     #Define number of cores for parallel computing
#     registerDoParallel(cores=2)

    #Loop over all selected grid points
    HWMId_all <- foreach(n=1:N_city, .combine='cbind') %do% {
    # for (n in 1:2) {

        #Select data on grid point
        data_sel <- data_MOD[, n]

        #Calculate HWMId
        hwmid_out <- hwmid(data_sel, ref_years, app_years, time_in)

        return(hwmid_out$hwmid)

    }

    #Delete input data arrays to save memory
    remove(data_MOD)

    #Read HWMID strength, HWMID length, and start DOY
    hwmid_str <- HWMId_all[,seq(1, N_city*3, 3)]
    hwmid_len <- HWMId_all[,seq(2, N_city*3, 3)]
    hwmid_doy <- HWMId_all[,seq(3, N_city*3, 3)]

    #Select time within application years
    sel_time_out = (time_in$year>=app_years[1]) & (time_in$year<=tail(app_years, 1))

    #Create time vector in application years
    time_out <- NULL
    time_out$day   <- time_in$day[sel_time_out]
    time_out$month <- time_in$month[sel_time_out]
    time_out$year  <- time_in$year[sel_time_out]

    #Select time to put in output NetCDF
    time_app = time[sel_time_out]
    time_app <- time_app[time_out$day==30 & time_out$month==6]

    #Create empty arrays for output
    HWMID_strength <- array(NA, dim=c(N_city, length(time_app) ))
    HWMID_length   <- array(NA, dim=c(N_city, length(time_app) ))
    HWMID_DOYstart <- array(NA, dim=c(N_city, length(time_app) ))

    #Write data into output arrays
    for (n in 1:N_city) {
        HWMID_strength[n, ] <- hwmid_str[,n]
        HWMID_length[n, ]   <- hwmid_len[,n]
        HWMID_DOYstart[n, ] <- hwmid_doy[,n]
    }

    # Mask infinite values
    HWMID_strength[is.infinite(HWMID_strength)] = NA


    #### Save in NetCDF ####

    #Define dimensions
    citydim <- ncdim_def("city", "1", 1:N_city) 
    timedim <- ncdim_def("time", units=tunits$value,as.double(time_app), calendar=calendar)

    #Define variables
    fillvalue <- 1e32
    tmp_def1 <- ncvar_def("HWMID", "1", list(citydim,timedim), fillvalue, "HWMID", prec="single", compression=2)
    tmp_def2 <- ncvar_def("HWMID_length", "days", list(citydim,timedim), fillvalue, "HWMID_length", prec="single", compression=2)
    tmp_def3 <- ncvar_def("DOY_start", "DOY", list(citydim,timedim), fillvalue, "DOY_start", prec="single", compression=2)

    #Create NetCDF and put variable
    ncout <- nc_create(fname_out, list(tmp_def1, tmp_def2, tmp_def3), force_v4=TRUE)
    waldo <- ncvar_put(ncout, tmp_def1, HWMID_strength)
    waldo <- ncvar_put(ncout, tmp_def2, HWMID_length)
    waldo <- ncvar_put(ncout, tmp_def3, HWMID_DOYstart)

    #Put additional attributes
    ncatt_put(ncout, "time", "axis", "t")

    #Close file
    waldo <- nc_close(ncout)

}