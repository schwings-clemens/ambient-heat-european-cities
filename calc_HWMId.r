##------------ ---------hwmid--------------------------
## Master function hwmid to call to creat HWMId indices for a specific grid point
## Inputs:
## - data: vector containing temperature data
## - ref_years: 2-element vector containing start and end year of reference period
## - app_years: vector containing all years of the whole period (one entry per year)
## - time_in: list containing:
##            - DOY: day of year index (same size as data)
##            - day: day of month (same size as data)
##            - month: month index (same size as data)
##            - year:  year index (same size as data)
##            - calendar: type of calendar (e.g., gregorian, 365_day, 360_day,...)
## - q_thres: the percentile level for the threshold
## - win_half: Half number of days considered around it each day (actual window = 2*win_half + 1)
## - hw_length: minimum length of heat wave
## Output is a list with the following components:
## - hwmid: Matrix of size number of years times 3 containing:
##          column 1) HWMId of strongest heat wave
##          column 2) duration of strongest heat wave
##          column 3) Start day (DOY) of strongest heat wave
## - thr: vector of q_thres percentile thresholds for each day of the year
## - Tmax_25p: 25th percentile of yearly maximum values for tasmax in reference period
## - Tmax_75p: 75th percentile of yearly maximum values for tasmax in reference period
##
hwmid <- function(data, ref_years, app_years, time_in, q_thres=90, win_half=15, hw_length=3) {

    #Get number of year days and delete 29th February in leap year calendars
    if (time_in$calendar=='360_day') {

        ydays <- 360

        } else {

        #Delete 29th February (for 365_day calendars it does not affect the data)
        delete_29 <- (time_in$day==29) & (time_in$month==2)
        data          <- data[!delete_29]
        time_in$DOY   <- time_in$DOY[!delete_29]
        time_in$day   <- time_in$day[!delete_29]
        time_in$month <- time_in$month[!delete_29]
        time_in$year  <- time_in$year[!delete_29]
        
        ydays <- 365
    }

    #Select data, DOY, and years in reference period
    sel_ref  <- (time_in$year>=ref_years[1]) & (time_in$year<=ref_years[2])
    data_ref <- data[sel_ref]
    DOY_ref  <- time_in$DOY[sel_ref]
    YEAR_ref <- time_in$year[sel_ref]

    #Calculate temperature threshold for each day of the year
    thresh_out <- mvThreshold(data_ref, DOY_ref, YEAR_ref, q_thres, win_half, ydays)

    threshold <- thresh_out$threshold
    Tmax_25p  <- thresh_out$Tmax_25p
    Tmax_75p  <- thresh_out$Tmax_75p
    Tmax_IQR  <- thresh_out$Tmax_IQR
    
    #Calculate HWMId
    hwmid <- hwmidFun(data, threshold, Tmax_25p, Tmax_IQR, time_in$year, app_years, hw_length)

    return(list(hwmid=hwmid, thr=threshold, Tmax_25p=Tmax_25p, Tmax_75p=Tmax_75p))
}



##--------------------------mvThreshold------------------------##
## Function to calculate temperature thresholds of a given percentile for each day of the year
## Inputs are:
## - data_ref: vector with temperature data in reference period
## - DOY_ref: vector with days of year (DOY) in reference period
## - q_thres: the percentile level for the threshold
## - win_half: Half number of days considered around each day (actual window = 2*win_half + 1)
## - ydays: is the number of DOYs
## Output:
## - threshold: vector of q_thres percentile thresholds for each day of the year
##
mvThreshold <- function(data_ref, DOY_ref, YEAR_ref, q_thres, win_half, ydays){ 

    #Create the empty threshold vector
    threshold <- rep(0,ydays)    

    #Loop over all DOYs
    for (i in seq(1,ydays)) {
        
        #Select all DOYs in reference period
        DOY_sel <- DOY_ref
#         ind_out <- i - win_half
        
        #Correct selection index for periods that overlap years
        if (i<=win_half){
            
            ind <- ydays - win_half
            
            DOY_sel[DOY_sel>=ind] <- DOY_sel[DOY_sel>=ind] - ydays
            ind_out <- i + ydays - win_half
            
        } else if (i>=ydays-win_half) {
            
            ind <- win_half
            DOY_sel[DOY_sel<=ind] <- DOY_sel[DOY_sel<=ind] + ydays
            
        }

        #Select data within window
        sel_DOY  <- DOY_sel>=(i-win_half) & DOY_sel<=(i+win_half)
        data_sel <- data_ref[sel_DOY]
        
        #Calculate quantile
        threshold[i] <- quantile(data_sel, probs=q_thres/100., na.rm=TRUE, names=FALSE, type=8)
    }
    
    #Calculate maximum yearly temperatures in reference period
    Tmax <- NULL
    for (year in min(YEAR_ref):max(YEAR_ref)) {
        sel  <- YEAR_ref==year
        Tmax <- c(Tmax, max(data_ref[sel]))
    }

    #Get 25th and 75th percentiles and IQR of maximum yearly temperatures in reference period
    Tmax_25p <- quantile(Tmax, 0.25, na.rm=TRUE)
    Tmax_75p <- quantile(Tmax, 0.75, na.rm=TRUE)
    Tmax_IQR <- Tmax_75p - Tmax_25p
    
    return(list(threshold=threshold, Tmax_25p=Tmax_25p, Tmax_75p=Tmax_75p, Tmax_IQR=Tmax_IQR))
}


##----------------------hwmidFun------------------------##
## Function that calculates hwmid indices for a temperature vector
## Inputs are:
## - data: vector with temperature data
## - threshold: vector of q_thres percentile thresholds for each day of the year
## - Tmax_25p: 25th percentile of yearly maximum values for tasmax in reference period
## - IQR: interquartile range of yearly maximum values for tasmax in reference period
## - years: vector containing all years, same size as data
## - app_years: vector containing all years, only one entry per year
## - hw_length: minimum length of heat waves
## Output:
## - Matrix of size number of years times 3 containing:
##   column 1) HWMId of strongest heat wave first column, duration of strongest heat wave in 
##   column 2) duration of strongest heat wave
##   column 3) Start day (DOY) of strongest heat wave
##
hwmidFun <- function(data, threshold, Tmax_25p, IQR, years, app_years, hw_length) {
    
    #Create empty matrix for storing output
    hw.scale <- matrix(0, ncol=3, nrow=length(app_years))

    # Looping over each year
    ii <- 1
    for (year in app_years) {   

        #Select and mask data
        sel <- years==year
        Ts <- data[sel]
        Ts[is.na(Ts)] <- -999.9
        
        #Get heat waves
        hw <- hwyearfast(Ts, threshold, hw_length)
        
        HWD <- hw[,1]
        HWtime <- hw[,2]
        
        #Create a set of zeros with length=number of heat waves
        HWMId <- rep(0,length(HWD))

        #If first heatwave was several days:
        if (HWD[1]>0) {

            #Then looping over all the heatwaves
            for (id in 1:length(HWD)) {

                #Length of heat wave
                hwd <- HWD[id]

                #First day of heatwave (as DOY)
                hwt <- HWtime[id]

                #Temperatures during period
                thw <- Ts[hwt:(hwt + hwd - 1)]

                #Calculate heat wave magnitude, normalized by IQR
                thwN   <- thw - Tmax_25p
                fvalue <- thwN / IQR
                
                #Save heat wave magnitude sum in HWMId where thw > Tmax_25p
                HWMId[id] <- sum(fvalue[thwN>0])                

            }
        }

        #Check if there is a HWMId>0 in that year
        if (max(HWMId)>0) {

            #Find the largest HWMId
            indtime <- which(HWMId==max(HWMId))

            #If several equally strong heatwaves, select first
            if (length(indtime) > 1) indtime <- indtime[1]

            #Matrix that will be returned
            hw.scale[ii, 1] <- max(HWMId)      # Heat wave magnitude index
            hw.scale[ii, 2] <- HWD[indtime]    # Duration of strongest heat wave
            hw.scale[ii, 3] <- HWtime[indtime] # Start time of strongest heat wave (as DOY)

        } 

        #Increase index by 1
        ii <- ii + 1
        
    } # end of for 'year' loop.
    
    return(hw.scale)

}


##---------------------------hwyearfast--------------------------##
## Function to find starting dates and durations of all heatwaves in a year
## Inputs:
## - Ti: temperature vector for one year
## - threshold: vector of thresholds for each day of the year
## - nd: is the minimum number of consecutive days above "threshold"
##   to be considered as heat wave
## Outputs:
## - Matrix of size number of heatwaves times 2 containing:
##   column 1) the duration of each heatwave 
##   column 2) the starting day of each heat wave (as DOY)
##
hwyearfast <- function(Ti, threshold, nd) {

    #Finding the indices where the temperature exceeds the threshold
    index <- Ti > threshold

    #And the ones where it doesn't
    z1 <- which(index==FALSE)

    #This counts the number of heatwaves
    HWF <- 0

    #These are the values we want to calculate
    HWtime <- HWD <- NULL

    #If all the days in the year exceed threshold (if branch has special cases for 0 or one non-exceedance day a year) 
    if (sum(index)>=length(Ti)) {

        #Duration is the whole year, the start time is January 1st
        HWD <- sum(index)
        HWtime <- 1

    #The case if there is only one day below the threshold:
    } else if (length(z1)==1) {

        #If this day is after the first few days, so that the first heat period is an actual heat wave:
        if (z1[1]>nd) {

            #This is the first heatwave and starts at January 1st
            HWF <- 1
            HWtime[1] <- 1

            # Duration of the first heatwave
            HWD[1] <- sum(index[ 1:(z1[1] - 1) ])

        }

        # First checking that z1 is not the very last day of the year. Then checking
        # that there are at least nd heatwave days after the last (and only)
        # non-heatwave day. I.e. that the last heat period is an actual second heatwave
        if (max(z1) < length(Ti) && sum(index[ max(z1):length(Ti) ]) >= nd) {

            #This will be 1 or 2 depending on whether the first period was a real heat wave
            HWF <- HWF + 1

            #Duration is the above threshold days after the below threshold day:
            HWD[HWF] <- sum(index[ max(z1):length(Ti) ])

            #Start date:
            HWtime[HWF] <- max(z1) + 1
            
        }

    #All the rest with more than one non-threshold exceeding day:   
    } else {

        # Test to check if there is a full heatwave before the first non-exceedence day:
        if (z1[1]>nd) {

            #This is the first heatwave and starts at January 1st
            HWF <- 1
            HWtime[HWF] <- 1
            
            # Duration of the first heatwave
            HWD[HWF] <- sum(index[ 1:(z1[1] - 1) ])

        }

        #Summing over the non-exceedence days, to find inner year heatwaves:
        for (i in 1:(length(z1)-1)) {

            #Check first to find heatwave between this and the next non-exceedence day:
            HWD0 <- sum(index[ (z1[i] + 1):(z1[i+1] - 1) ])

            #If this time is enough for a heatwave:
            if (HWD0>=nd) {

                #New heatwave at this time and of this length
                HWF <- HWF + 1
                HWtime[HWF] <- (z1[i] + 1)
                HWD[HWF] <- HWD0

            }

        }

        #Check if there is a final heatwave at the end of the year
        if (max(z1) < length(Ti) && sum(index[ max(z1):length(Ti) ]) >= nd) {

            #Adding number, duration and date:
            HWF <- HWF + 1
            HWD[HWF] <- sum(index[ max(z1):length(Ti) ])
            HWtime[ HWF ] <- max(z1) + 1

        }

    } # end of if length of z > 1 stmts (having covered all heatwaves)

    #Making a list of the heatwave values
    HW <- cbind(HWD, HWtime)

    #Sorting out the dimensions if there are actual heatwaves in the data
    if (length(HWD)>0) {
        dim(HW) <- c(length(HWD),2)
    } else{
        #Set negative values if there is no heat wave
        HW <- c(-1111,-1111)
        dim(HW)<-c(1,2)
    }

    return(HW)

}
