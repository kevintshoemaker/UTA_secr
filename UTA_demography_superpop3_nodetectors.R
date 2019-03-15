
# TODO
   # incorporate body size element... ?? 


##########
# CLEAR WORKSPACE
##########

rm(list=ls()) 

##########
# LOAD PACKAGES
##########

library(R2jags)
library(runjags)
library(coda)
library(plyr)
library(sp)
library(raster)
library(rgeos)
#library(secr)
library(geoR)
library(spatstat)
library(rgdal)
library(Hmisc)
library(abind)
library(lubridate)

########
# LOAD FUNCTIONS
########

LongLatToUTM<-function(xy){     # zone
  xy2 <- data.frame(ID = 1:nrow(xy), X = xy[,1], Y = xy[,2])
  sub <- !is.na(xy2$X)
  out <- is.na(xy2$X)
  xy3 <- xy2[sub,]
  if(nrow(xy3)>0){
    coordinates(xy3) <- c("X", "Y")
    proj4string(xy3) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    res <- spTransform(xy3, CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
  }else{
    res=xy3
  }
  df <- xy2
  df[sub,] <- as.data.frame(res)
  return(df)
}

# ConvertToHabmat<-function(xy,oldxlim,oldylim,nrow,ncol){     # zone
#   
#   return(df)
# }

########
# READ IN DATA
########

allcaps<-read.csv("allcaptures.csv",stringsAsFactors = F,na.strings = c("NA","<NA>","n/a"))
Surveys<-read.csv("SurveyMetadata2.csv",stringsAsFactors = F,na.strings = c("NA","<NA>","n/a",""))
names(allcaps)
names(Surveys)

head(allcaps)
head(Surveys)

tail(Surveys)

keep <- !is.na(Surveys$Date)

Surveys <- Surveys[keep,]


names(Surveys)

Surveys$Date2 <- mdy(Surveys$Date)

allcaps[which(allcaps$ID==c("42")),"ID"] <- "PH42"

allcaps <- allcaps[!is.na(allcaps$Occasion),]     # remove observations with no occasion

   ### remove spatial outliers that should not really be part of the study (not really in surveyed area)

## site 7
   ##  1036 1037 1130 1197 1212

remove_site7 <- c(1036,1037,1130,1197,1212)
allcaps <- allcaps[-remove_site7,]



#######
# DATA PROCESSING
#######

allcaps$Occasion <- as.character(allcaps$Occasion)

site_data <- list()     # set up the main data storage structure for site-level data

   # ensure a common naming convention for sites
uniquesites <- unique(Surveys$Site)[!is.na(unique(Surveys$Site))]    # 9 unique sites
nuniquesites <- length(uniquesites)
newnames <- c("DH","DI","I-10","ME","MI","MV","SG","WW","PH")
sitenames_df <- data.frame(oldnames=uniquesites,newnames=newnames,stringsAsFactors = FALSE)

## number of sites and periods

global_vars <- list()
global_vars$nsites <- 9
global_vars$n_primary_occasions <- 10
global_vars$nperiods <- numeric(global_vars$nsites)

## control for data augmentation

global_vars$naug <-  c(200,200,200,200,200,200,300,200,200)   # c(200,150,200,200,200,200,300,200,200) #  

global_vars$max_naug <- max(global_vars$naug)

global_vars$n_secondary_occasions <- max(Surveys$Session,na.rm=T)
global_vars$primary_occasions <- c(1:global_vars$n_primary_occasions)
global_vars$secondary_occasions <- c(1:global_vars$n_secondary_occasions)

global_vars$occasion_names <- paste(rep(global_vars$primary_occasions,each=global_vars$n_secondary_occasions),rep(global_vars$secondary_occasions,times=global_vars$n_primary_occasions),sep="")

occasion_df <- data.frame(
  Occasion = as.character(global_vars$occasion_names),
  Period = rep(global_vars$primary_occasions,each=global_vars$n_secondary_occasions),
  Replicate = rep(global_vars$secondary_occasions,times=global_vars$n_primary_occasions),
  stringsAsFactors = F
)

occasion_df$Occasion

## first beadmarked

global_vars$first.beadmarked <- 4


## define the periods

global_vars$period_names <- c("July 2013","August 2013","March 2014","June 2014","June 2014","July 2014","August 2014","March 2015","May 2015","June 2015")
global_vars$intervals <- c(12,222,78,26,27,21,209,64,31)/30   # in months

## first sampling period


global_vars$firstdate <- min(Surveys$Date2)+days(5)


global_vars$realdates <- global_vars$firstdate + days(cumsum(c(0,12,222,78,26,27,21,209,64,31)))

i=7
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  site2 <- sitenames_df$oldnames[sitenames_df$newnames==site]
  temp <- subset(Surveys,Site==site2)
  site_data[[site]] <- list()
  site_data[[site]]$periods_surveyed <- unique(temp$Period) 
  site_data[[site]]$nperiods <- length(site_data[[site]]$periods_surveyed)
  global_vars$nperiods[i] <- site_data[[site]]$nperiods
  
  site_data[[site]]$replicate_sessions <- numeric(global_vars$n_primary_occasions)
  times_surveyed <- table(temp$Period)
  site_data[[site]]$replicate_sessions[as.numeric(names(times_surveyed))] <- as.numeric(times_surveyed)
  j=1
  site_data[[site]]$surveyeffort <- rep(NA,times=global_vars$n_primary_occasions)
  
  temp2 <- tapply(temp$SurveyTime*temp$X..surveyors,temp$Period,sum)         # survey effort
  temp2[which(is.na(temp2))] <- mean(temp2,na.rm=TRUE)
  ndx <- as.numeric(names(temp2))
  site_data[[site]]$surveyeffort[ndx] <- temp2
  
  site_data[[site]]$covsbysession <- list()
  p=1
  for(p in site_data[[site]]$periods_surveyed){
    temp2 <- subset(temp,Period==p)
    site_data[[site]]$covsbysession[[p]] <- list()
    site_data[[site]]$covsbysession[[p]]$effort <- temp2$SurveyTime*temp2$X..surveyors
    site_data[[site]]$covsbysession[[p]]$wind <- temp2$MaxWind
    site_data[[site]]$covsbysession[[p]]$is.pm <- as.numeric(temp2$Type=="PM")
    site_data[[site]]$covsbysession[[p]]$nsessions <- max(temp2$Session)
  }
  
}


#######
# GET SIZE AND BODY CONDITION 
#######




#######
# GET XY COORDINATES FOR ALL INDIVIDUALS AT EACH SITE
#######

i=1
for(i in 1:nuniquesites){
  site <- sitenames_df$newnames[i]
  temp <- subset(allcaps,Site==site)
  site_data[[site]]$capxy <- temp[,c("GPS.W","GPS.N")]
 
  
  #     # find outliers
  # keep <- !is.na(site_data[[site]]$capxy$GPS.W)
  # temp2 <- temp[keep,]
  # temp3 <- temp2[chull(site_data[[site]]$capxy[keep,]),]
  # 
  # plot(site_data[[site]]$capxy) 
  # text(temp3$GPS.W,temp3$GPS.N,labels=temp3$Lizard,col="red",cex=2)
  # 
  # temp3
  # keep <- c(1,5)
  # temp3 <- temp3[-keep,]
  # to_remove <- which((allcaps$GPS.W%in%temp3$GPS.W)&(allcaps$GPS.N%in%temp3$GPS.N))
  # 
    #   ## remove na from xy coords (not needed any more!)
    # ndx <- !is.na(site_data[[site]]$capxy[,1])
    # site_data[[site]]$capxy <- site_data[[site]]$capxy[ndx,]
  site_data[[site]]$capxy <- LongLatToUTM(site_data[[site]]$capxy)[,c("X","Y")]
    #temp <- temp[ndx,]
  site_data[[site]]$uniqueinds <- unique(temp$ID)
  site_data[[site]]$ind_locs <- list()
  ind = site_data[[site]]$uniqueinds[31]
  site_data[[site]]$init.locy <- numeric(length(site_data[[site]]$uniqueinds))
  site_data[[site]]$init.locx <- numeric(length(site_data[[site]]$uniqueinds))
  counter=1
  for(ind in site_data[[site]]$uniqueinds){
    temp2 <- subset(temp,ID==ind)
    site_data[[site]]$ind_locs[[ind]] <- LongLatToUTM(temp2[,c("GPS.W","GPS.N")])[,c("X","Y")]
    site_data[[site]]$ind_locs[[ind]]$periods <- occasion_df$Period[match(temp2$Occasion,occasion_df$Occasion)]
    site_data[[site]]$ind_locs[[ind]]$replicates <- occasion_df$Replicate[match(temp2$Occasion,occasion_df$Occasion)]
    site_data[[site]]$init.locx[counter] <- mean(site_data[[site]]$ind_locs[[ind]]$X[site_data[[site]]$ind_locs[[ind]]$periods==min(site_data[[site]]$ind_locs[[ind]]$periods)])
    site_data[[site]]$init.locy[counter] <- mean(site_data[[site]]$ind_locs[[ind]]$Y[site_data[[site]]$ind_locs[[ind]]$periods==min(site_data[[site]]$ind_locs[[ind]]$periods)])
    counter=counter+1
  }
}

site <- sitenames_df$newnames[1] 
plot(site_data[[site]]$capxy)

sapply(site_data[[site]]$capxy,max)-sapply(site_data[[site]]$capxy,min)

site_data[[1]]$init.locx
site_data[[1]]$init.locy

#######
# GET SEX AND SVL FOR ALL INDIVIDUALS
#######

## SEX AND SVL

# NOTE: SVL was generally only taken once, so we don't really get to see how it changes over time. Maybe not as important to include this in the model... 
    # look at subadults vs adults?  subadult: less than 46 mm SVL for males or 42 mm for females
    # don't look at juveniles vs adults in the model- too much to consider the recruitment of juveniles and adults? Just do post-hoc analysis

i=1
for(i in 1:nuniquesites){
  site = sitenames_df$newnames[i]
  temp = subset(allcaps,Site=site)
  site_data[[site]]$sex <- character(length(site_data[[site]]$uniqueinds))
  site_data[[site]]$is.male <- numeric(length(site_data[[site]]$uniqueinds))
  site_data[[site]]$svl <- array(NA,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions))
  site_data[[site]]$weight <- array(NA,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions))
  site_data[[site]]$age <- array(NA,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions))
  ind=site_data[[site]]$uniqueinds[1]
  counter <- 1 
  for(ind in site_data[[site]]$uniqueinds){
    temp2 <- subset(temp,ID==ind)
    sex <- names(which.max(table(temp2$Sex[temp2$Sex%in%c("m","f")])))
    if(length(sex)>0){
      allcaps$Sex[allcaps$ID==ind] <- sex
    }else{
      allcaps$Sex[allcaps$ID==ind] <- NA
    }
    site_data[[site]]$sex[counter] <- allcaps$Sex[allcaps$ID==ind][1]
    site_data[[site]]$is.male[counter] <- ifelse(site_data[[site]]$sex[counter]=="m",1,0)
    period=4
    for(period in 1:global_vars$n_primary_occasions){
      temp3 <- subset(temp2,Occasion%in%(occasion_df$Occasion[occasion_df$Period==period]))
      if(nrow(temp3)>0){
        if(any(!is.na(temp3$SVL))) site_data[[site]]$svl[counter,period] <- mean(temp3$SVL,na.rm=T)
        if(any(!is.na(temp3$SVL))) site_data[[site]]$age[counter,period] <- ifelse(site_data[[site]]$svl[counter,period] < (site_data[[site]]$is.male[counter]*46 + (1-site_data[[site]]$is.male[counter])*42),1,2)
        if(any(!is.na(temp3$Weight))) site_data[[site]]$weight[counter,period] <- mean(temp3$Weight,na.rm=T)
      }
    }
    counter <- counter + 1
  }
}


### Visualize locations of lizards within a study site

site <- sitenames_df$newnames[2] 
plot(site_data[[site]]$capxy)


### NOTE: the plots don't seem perfectly rectangular. Maybe use MCP on the captures to define study site edges? 

#######
# SET UP MCPs FOR EACH SURVEY SITE
#######

i=1
for(i in 1:nuniquesites){
  site <- sitenames_df$newnames[i] 
  keep <- !is.na(site_data[[site]]$capxy$X)
  capxy <- site_data[[site]]$capxy[keep,]
  s<- chull(capxy[,2], capxy[,1])    # convex hull points
  mcp<-capxy[c(s,s[1]),]
  site_data[[site]]$mcp <- SpatialPolygons(list(Polygons(list(Polygon(mcp)), ID=1)))   #  sampled area
  site_data[[site]]$mcp <- buffer(site_data[[site]]$mcp, width=2)   #  expand sampled area by a little...
  site_data[[site]]$mcp.b <- buffer(site_data[[site]]$mcp, width=10)           # buffer around sampled area
  site_data[[site]]$area <- site_data[[site]]$mcp.b@polygons[[1]]@area
  
  polygon <- site_data[[site]]$mcp.b
  xy <- as.data.frame(polygon@polygons[[1]]@Polygons[[1]]@coords)
  
  site_data[[site]]$xlim <- c(min(xy$x),max(xy$x))
  site_data[[site]]$ylim <- c(min(xy$y),max(xy$y))
  
  site_data[[site]]$area2 <- diff(site_data[[site]]$xlim)*diff(site_data[[site]]$ylim)
  
  
  #Create a raster of cells inside our sampled area    CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  r <- raster(crs=CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
  extent(r) <- extent(site_data[[site]]$mcp.b)
  res(r) <- c(2, 2)  #resolution of raster: 2m resolution seems good enough...
  values(r) <- 1:ncell(r)
  r <- raster::mask(r, site_data[[site]]$mcp)
  values(r) <- ifelse(is.na(values(r)), 0, 1)
  #raster::plot(r)
  
  site_data[[site]]$habraster <- r
  
  site_data[[site]]$habmat <- as.matrix(r)
  xu <- ncol(r)
  yu <- nrow(r)
  
  site_data[[site]]$xlim2 <- c(0,xu)
  site_data[[site]]$ylim2 <- c(0,yu)
  
  #Convert location coordinates into "habmat" coordinates (rows and columns in habmat)
  y <- ((site_data[[site]]$ylim[2]-site_data[[site]]$capxy$Y) /  diff(site_data[[site]]$ylim)) * yu
  x <- ((site_data[[site]]$capxy$X-site_data[[site]]$xlim[1]) /  diff(site_data[[site]]$xlim)) * xu
  site_data[[site]]$capxy2 <- data.frame(X=x,Y=y)
  
  # debug...
  # convert <- trunc(site_data[[site]]$capxy2+1)
  # site_data[[site]]$habmat[as.matrix(convert[,c(2,1)])]
  
}


######################
# VISUALIZE MCPs
######################


site=1
plot(site_data[[site]]$mcp)
plot(site_data[[site]]$mcp.b)
points(site_data[[site]]$capxy)

names(site_data)


more <- 10


graphics.off()
par(mfrow=c(2,3))
par(mai=c(.7,.7,.1,.1))
#par(ask=T)
site=1
for(site in 1:global_vars$nsites){
  sitename = sitenames_df$newnames[site]
  ext <- extent(site_data[[site]]$mcp.b)
  plot(site_data[[site]]$capxy,
       ylim=c(ext@ymin-more,ext@ymax+more),xlim=c(ext@xmin-more,ext@xmax+more),
       xlab="",ylab="",main=""
  )
  plot(site_data[[site]]$habraster,add=T,legend=F)
  points(site_data[[site]]$capxy)
  
  
  plot(ext,lwd=5,add=T)
  legend("topleft",cex=1.2,legend=sitename,bg="grey")

}



site=1
habrast <- raster(site_data[[site]]$habmat,xmn=site_data[[site]]$xlim2[1],xmx=site_data[[site]]$xlim2[2],ymn=site_data[[site]]$ylim2[1],ymx=site_data[[site]]$ylim2[2])
plot(habrast)
cap2 <- site_data[[site]]$capxy2
cap2$Y <- site_data[[site]]$ylim2[2]-cap2$Y
points(cap2)


#######
# SET UP CAPTURE HISTORY MATRICES!
#######


i=1
for(i in 1:nuniquesites){
  site = sitenames_df$newnames[i]
  site_data[[site]]$caphist.std.2d <- array(0,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions))     # for storing the standard capture history 
  site_data[[site]]$caphist.std.3d <- array(0,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions,global_vars$n_secondary_occasions))    # for storing the standardy

  site_data[[site]]$caplocationx.2d <- array(0,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions))     # for storing the standard capture history 
  site_data[[site]]$caplocationx.3d <- array(0,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions,global_vars$n_secondary_occasions)) 
 
  site_data[[site]]$caplocationy.2d <- array(0,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions))     # for storing the standard capture history 
  site_data[[site]]$caplocationy.3d <- array(0,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions,global_vars$n_secondary_occasions))
     
  site_data[[site]]$is.recap <- array(NA,dim=c(length(site_data[[site]]$uniqueinds),global_vars$n_primary_occasions))
  
  allcaps.thissite <- subset(allcaps,(Site==site))
  
  j=4
  for(j in 1:global_vars$n_primary_occasions){
    occasions <- occasion_df$Occasion[occasion_df$Period==j]
    temp <- subset(allcaps,(Site==site)&(Occasion%in%occasions))
      
    uniqueids <- unique(temp$ID)      # individuals captured at this site on this period
    nuniqueids <- length(uniqueids)
        
    if(nrow(temp)>0){
      temp2 <- LongLatToUTM(temp[,c("GPS.W","GPS.N")])[,c("X","Y")]
      temp$X <- temp2$X
      temp$Y <- temp2$Y
    }else{
      temp$X <- numeric(0)
      temp$Y <- numeric(0)
    }
    
    capsbyoccasion <- table(temp$ID,temp$Occasion)
    
    ## assemble the capture histories
    if(nrow(temp)>0){
      k=1
      for(k in 1:nuniqueids){
        thisid <- uniqueids[k]
        
        idndx <- which(site_data[[site]]$uniqueinds==thisid)
        
        site_data[[site]]$caphist.std.2d[idndx,j] <- 1
        
        caplocs <- site_data[[site]]$ind_locs[[thisid]][site_data[[site]]$ind_locs[[thisid]]$periods==j,]
        
        if(nrow(caplocs)>0){
          site_data[[site]]$caplocationx.2d[idndx,j] <- mean(caplocs$X,na.rm=T)
          site_data[[site]]$caplocationy.2d[idndx,j] <- mean(caplocs$Y,na.rm=T)
          site_data[[site]]$caplocationx.3d[idndx,j,caplocs$replicates] <- caplocs$X
          site_data[[site]]$caplocationy.3d[idndx,j,caplocs$replicates] <- caplocs$Y
        } 
        
        sessions_detected <- occasion_df$Replicate[match(as.numeric(colnames(capsbyoccasion)),occasion_df$Occasion)][capsbyoccasion[thisid,]==1]
            
        site_data[[site]]$caphist.std.3d[idndx,j,sessions_detected] <- 1
        
        occasions_captured <- occasion_df$Period[occasion_df$Occasion%in%allcaps.thissite[allcaps.thissite$ID==uniqueids[k],]$Occasion]
        
        if(length(occasions_captured)>0){
          site_data[[site]]$is.recap[idndx,j] <- ifelse(j>min(occasions_captured),1,0)
        }
        
      }
    }
  }
}


site=2
site_data[[site]]$caphist.std.2d


site_data[[site]]$sex

site_data[[site]]$is.recap

site_data[[site]]$caplocationy.2d

site_data[[3]]$caphist.std.2d

site_data[[3]]$caphist.std.3d[,,3]


###############
# CONSIDER SITE-LEVEL COVARIABLES
###############


sitecovars_df <- sitenames_df

    ### Wind power: sites DI, ME, MV, and PH
sitecovars_df$windpower <- c(0,1,0,1,0,1,0,0,1)

   ### other disturbance
sitecovars_df$disturbance_rank <- c(1,4,3,4,0,4,0,1,4)   # wind coded as 4?

sitecovars_df$disturbance <- c(24.17,34.53,67.47,26.48,22.82,49.89,2.61,40.23,59.97)

sitecovars_df$bareground <- c(82,90,72,76,63,87,84,75,85)

sitecovars_df$noiselevel <- c(48.17,49.85,54.43,47.92,45.45,60.49,41.16,47.62,66.37)

sitecovars_df$shrubdensity <- c(0.11,0.07,0.14,0.23,0.38,0.09,0.06,0.15,0.14)

sitecovars_df$shrubcanopy <- c(18.92,12.83,19.81,12.97,27.50,3.88,2.95,17.43,13.00)

sitecovars_df$shrubheight <- c(0.97,1.01,0.85,0.62,0.56,0.48,0.61,0.66,0.70)

sitecovars_df$slope <- c(4.12,1.51,8.67,8.84,5.84,1.92,17.40,2.72,13.21)

sitecovars_df$elevation <- c(500.99,318.03,390.74,711.98,640.74,366.12,822.73,468.75,464.42)

sitecovars_df$rainfall <- c(10.82,6.57,11.81,18.43,13.86,15.26,16.95,13.99,13.42)

sitecovars_df$avgwind <- c(2.58,6.84,13.13,8.40,3.64,10.15,3.80,5.65,8.25)

sitecovars_df


###############
# DATA AUGMENTATION   [change to allow different augmentation for different sites!]
###############

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i] 
  site_data[[site]]$caphist.std.3d.aug <- abind(site_data[[site]]$caphist.std.3d,array(0,dim=c((global_vars$naug[i]-length(site_data[[site]]$uniqueinds)),global_vars$n_primary_occasions,global_vars$n_secondary_occasions)),along=1)
}


###############
# SET UP DATA FOR BUGS
###############

    ## master capture history matrix

caphist4d <- array(NA, dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions,global_vars$n_secondary_occasions))
caplocsx4d <- array(NA, dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions,global_vars$n_secondary_occasions))
caplocsy4d <- array(NA, dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions,global_vars$n_secondary_occasions))

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  j=5
  for(j in 1:length(site_data[[site]]$periods_surveyed)){
    this.period <- site_data[[site]]$periods_surveyed[j]
       #ncaptures <- apply(site_data[[site]]$caphist.std.2d,2,sum)[this.period]
    caphist4d[i,1:nrow(site_data[[site]]$caphist.std.2d),j,1:site_data[[site]]$replicate_sessions[this.period]] <- site_data[[site]]$caphist.std.3d[,this.period,1:site_data[[site]]$replicate_sessions[this.period]]
    caphist4d[i,(nrow(site_data[[site]]$caphist.std.2d)+1):global_vars$naug[i],j,1:site_data[[site]]$replicate_sessions[this.period]] <- 0
    
    caplocsx4d[i,1:nrow(site_data[[site]]$caphist.std.2d),j,1:site_data[[site]]$replicate_sessions[this.period]] <- site_data[[site]]$caplocationx.3d[,this.period,1:site_data[[site]]$replicate_sessions[this.period]]
    caplocsy4d[i,1:nrow(site_data[[site]]$caphist.std.2d),j,1:site_data[[site]]$replicate_sessions[this.period]] <- site_data[[site]]$caplocationy.3d[,this.period,1:site_data[[site]]$replicate_sessions[this.period]]
  }
}

caphist4d[1,2,,1]

caplocsy4d[1,2,,1]


caplocsx4d <- apply(caplocsx4d,c(1,2,3,4),function(t) ifelse(t==0,NA,t))
caplocsy4d <- apply(caplocsy4d,c(1,2,3,4),function(t) ifelse(t==0,NA,t))

caplocsy4d[1,2,,]

caplocsy4d[1,,1,1]



## IS.MALE

is.male <- array(NA,dim=c(global_vars$nsites,global_vars$max_naug))
i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  is.male[i,1:length(site_data[[site]]$uniqueinds)] <- site_data[[site]]$is.male
}

is.male[2,]


## INIT Z

init.z <- array(0,dim=c(global_vars$nsites,global_vars$max_naug))
i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  init.z[i,1:length(site_data[[site]]$uniqueinds)] <- 1
}

init.z[1,]


## Matrix of sessions surveyed for each primary occasion

temp <- t(sapply(1:global_vars$nsites,function(t) site_data[[t]]$replicate_sessions,simplify = "array"))[,,drop=FALSE]
sessions.surveyed <- t(apply(temp,1,function(t) c(t[which(t>0)],t[which(t==0)])))


## SESSION COVARIATES

is.pm <- array(NA,dim=c(global_vars$nsites,global_vars$n_primary_occasions,global_vars$n_secondary_occasions))
effort <- array(NA,dim=c(global_vars$nsites,global_vars$n_primary_occasions,global_vars$n_secondary_occasions))
wind <- array(NA,dim=c(global_vars$nsites,global_vars$n_primary_occasions,global_vars$n_secondary_occasions))
site_data[[1]]$covsbysession[[1]]$is.pm
i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  j=4
  for(j in 1:global_vars$nperiods[i]){
    this.period <- site_data[[site]]$periods_surveyed[j]
    is.pm[i,j,1:site_data[[site]]$covsbysession[[this.period]]$nsessions] <- site_data[[site]]$covsbysession[[this.period]]$is.pm
    effort[i,j,1:site_data[[site]]$covsbysession[[this.period]]$nsessions] <- site_data[[site]]$covsbysession[[this.period]]$effort
    wind[i,j,1:site_data[[site]]$covsbysession[[this.period]]$nsessions] <- site_data[[site]]$covsbysession[[this.period]]$wind
  }
}

    
    # Standardize!

effort <- (effort-mean(effort,na.rm=T))/sd(effort,na.rm=T)

wind <- (wind-mean(effort,na.rm=T))/sd(wind,na.rm=T)

    # interpolate

effort[which(is.na(effort),arr.ind = T)] <- 0

wind[which(is.na(wind),arr.ind = T)] <- 0

wind[2,,]

effort[3,,]

## FIRST SESSION WITH CAPTURES

##firstsession[site,period,ind]

firstsession <- array(NA,dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions))

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  k=3
  for(k in 1:global_vars$max_naug){
    j=1  
    for(j in 1:global_vars$nperiods[i]){
      this.period <- site_data[[site]]$periods_surveyed[j]
      temp <- caphist4d[i,k,j,]
      if(sum(temp,na.rm=T)>0){
        val <- which(temp>0)[1]
      }else{
        val <- 3
      }
      firstsession[i,k,j] <- val
    }
  }
}

firstsession[2,,]

##### Is beadmarked

is.beadmarked <- array(0,dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions))
#is.recap <- array(0,dim=c(global_vars$nsites,global_vars$naug,global_vars$n_primary_occasions))    

beadmark_mask <- is.beadmarked
beadmark_mask[,,(global_vars$first.beadmarked+1):global_vars$n_primary_occasions] <- 1
i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  site_data[[site]]$firsts <- apply(site_data[[site]]$caphist.std.2d,1,function(t) min(which(t==1)))
  
  temp <- beadmark_mask[i,1:length(site_data[[site]]$uniqueinds),]*site_data[[site]]$caphist.std.2d
  is.beadmarked[i,1:length(site_data[[site]]$uniqueinds),] <-  t(apply(temp,1,function(t) pmin(rep(1,length(t)),cumsum(t)) ))
  
}

is.beadmarked[1,,]

## INITIAL LOCATION!

init.location.x <- array(0,dim=c(global_vars$nsites,global_vars$max_naug))
init.location.y <- array(0,dim=c(global_vars$nsites,global_vars$max_naug))


caplocsx4d.init <- array(NA,dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions,global_vars$n_secondary_occasions))  #caplocsx4d
caplocsy4d.init <- array(NA,dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions,global_vars$n_secondary_occasions))  #caplocsy4d

i=4
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  newpolygon <- gBuffer(site_data[[site]]$mcp, width=-10)
  # plot(site_data[[site]]$mcp.b)
  # plot(site_data[[site]]$mcp,add=T)
  # plot(newpolygon,add=T)
  init.location.x[i,1:length(site_data[[site]]$uniqueinds)] <- site_data[[site]]$init.locx
  init.location.x[i,(length(site_data[[site]]$uniqueinds)+1):global_vars$max_naug] <- runif(length((length(site_data[[site]]$uniqueinds)+1):global_vars$max_naug),site_data[[site]]$xlim[1],site_data[[site]]$xlim[2])
  init.location.y[i,1:length(site_data[[site]]$uniqueinds)] <- site_data[[site]]$init.locy
  init.location.y[i,(length(site_data[[site]]$uniqueinds)+1):global_vars$max_naug] <- runif(length((length(site_data[[site]]$uniqueinds)+1):global_vars$max_naug),site_data[[site]]$ylim[1],site_data[[site]]$ylim[2])
  isna <- is.na(init.location.x[i,])
  if(any(isna)){
    randpoints <- spsample(newpolygon,length(which(isna)),type="random")
    randx <- randpoints@coords[,1]
    randy <- randpoints@coords[,2]
    init.location.x[i,isna] <- randx #runif(1,site_data[[site]]$xlim[1],site_data[[site]]$xlim[2])
    init.location.y[i,isna] <- randy #runif(1,site_data[[site]]$ylim[1],site_data[[site]]$ylim[2])
  }
       # now set up the initialization for the actual capture locations... 
  
  ind=1
  for(ind in 1:length(site_data[[site]]$uniqueinds)){
    indiv <- site_data[[site]]$uniqueinds[ind]
    for(p in 1:global_vars$n_primary_occasions){
      isna <- is.na(caplocsx4d[i,ind,p,])
      if(any(!isna)){
        if(any(isna)){ 
          caplocsx4d.init[i,ind,p,isna]<- mean(caplocsx4d[i,ind,p,!isna])
          caplocsy4d.init[i,ind,p,isna]<- mean(caplocsy4d[i,ind,p,!isna])
        }
      } else{
        caplocsx4d.init[i,ind,p,] <- mean(site_data[[site]]$ind_locs[[ind]]$X,na.rm=T)
        caplocsy4d.init[i,ind,p,] <- mean(site_data[[site]]$ind_locs[[ind]]$Y,na.rm=T)
      }
    }
  }
       # fill in augmented individuals with truly random points
  isna <- which(is.na(caplocsx4d[i,,,]),arr.ind = TRUE)
  nna <- nrow(isna)
  randpoints <- spsample(newpolygon,nna,type="random")
  randx <- randpoints@coords[,1]
  randy <- randpoints@coords[,2]
  caplocsx4d.init[i,,,][isna] <- randx
  caplocsy4d.init[i,,,][isna] <- randy
  
}

init.location.x[1,]
init.location.y[1,]

caplocsx4d.init[1,,1,]
caplocsx4d[1,,1,]


site=1
habrast <- raster(site_data[[site]]$habmat[1:site_data[[site]]$ylim2[2],1:site_data[[site]]$xlim2[2]],xmn=site_data[[site]]$xlim[1],xmx=site_data[[site]]$xlim[2],ymn=site_data[[site]]$ylim[1],ymx=site_data[[site]]$ylim[2])
plot(habrast)
res(habrast)
# points(init.location.x[site,],site_data[[site]]$ylim[2]-init.location.y[site,],pch=20)
# points(init.location.x[site,1:length(site_data[[site]]$uniqueinds)],init.location.y[site,1:length(site_data[[site]]$uniqueinds)],col="red",pch=1)



#### convert to "habmat" scale

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  caplocsx4d[i,,,] <- ((caplocsx4d[i,,,]-site_data[[site]]$xlim[1]) /  diff(site_data[[site]]$xlim)) * site_data[[site]]$xlim2[2]
  caplocsy4d[i,,,] <- ((site_data[[site]]$ylim[2]-caplocsy4d[i,,,]) /  diff(site_data[[site]]$ylim)) * site_data[[site]]$ylim2[2]
  init.location.x[i,] <- ((init.location.x[i,]-site_data[[site]]$xlim[1]) /  diff(site_data[[site]]$xlim)) * site_data[[site]]$xlim2[2]
  init.location.y[i,] <- ((site_data[[site]]$ylim[2]-init.location.y[i,]) /  diff(site_data[[site]]$ylim)) * site_data[[site]]$ylim2[2]
  caplocsx4d.init[i,,,] <- ((caplocsx4d.init[i,,,]-site_data[[site]]$xlim[1]) /  diff(site_data[[site]]$xlim)) * site_data[[site]]$xlim2[2]
  caplocsy4d.init[i,,,] <- ((site_data[[site]]$ylim[2]-caplocsy4d.init[i,,,]) /  diff(site_data[[site]]$ylim)) * site_data[[site]]$ylim2[2]  
}


initx.hr <- array(1,dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions))
inity.hr <- array(1,dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions))

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  for(ind in 1:global_vars$max_naug){
    initx.hr[i,ind,1:global_vars$n_primary_occasions] <- init.location.x[i,ind]
    inity.hr[i,ind,1:global_vars$n_primary_occasions] <- init.location.y[i,ind]
  }
  
}

initx.hr[1,,]


## INIT INPOP

init.inpop <- array(0,dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions))

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  k=1
  for(k in 1:length(site_data[[site]]$uniqueinds)){
    firstndx <- which(site_data[[site]]$periods_surveyed==site_data[[site]]$firsts[k])
    lastndx <- max(site_data[[site]]$ind_locs[[site_data[[site]]$uniqueinds[k]]]$periods)
    init.inpop[i,k,firstndx:lastndx] <- 1
  }
  for(k in (length(site_data[[site]]$uniqueinds)+1):global_vars$naug[i]){
    firstndx <- sample(1:(site_data[[site]]$nperiods-1),1)
    lastndx <- sample((firstndx+1):site_data[[site]]$nperiods,1)
    init.inpop[i,k,firstndx:lastndx] <- 1
  }
}


## HABITAT MATRIX!

maxx <- max(t(sapply(1:global_vars$nsites,function(t) site_data[[t]]$xlim2[2],simplify = "array")))
maxy <- max(t(sapply(1:global_vars$nsites,function(t) site_data[[t]]$ylim2[2],simplify = "array")))

habmat <- array(0,dim=c(global_vars$nsites,maxy,maxx))

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  habmat[i,1:site_data[[site]]$ylim2[2],1:site_data[[site]]$xlim2[2]] <- site_data[[site]]$habmat
}

habmat[1,,]


### ensure that all observations are within the mask, at least to start (for initialization)

site=1
habrast <- raster(habmat[site,1:site_data[[site]]$ylim2[2],1:site_data[[site]]$xlim2[2]],xmn=site_data[[site]]$xlim2[1],xmx=site_data[[site]]$xlim2[2],ymn=site_data[[site]]$ylim2[1],ymx=site_data[[site]]$ylim2[2])
plot(habrast)
res(habrast)
points(init.location.x[site,],site_data[[site]]$ylim2[2]-init.location.y[site,],pch=20)
points(init.location.x[site,1:length(site_data[[site]]$uniqueinds)],site_data[[site]]$ylim2[2]-init.location.y[site,1:length(site_data[[site]]$uniqueinds)],col="red",pch=1)


## INTERVALS!

intervals <- array(0,dim=c(global_vars$nsites,global_vars$n_primary_occasions))
i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  j=1
  for(j in 1:(global_vars$nperiods[i]-1)){
    this.period <- site_data[[site]]$periods_surveyed[j]
    next.period <- site_data[[site]]$periods_surveyed[j+1]
    intervals[i,j] <- sum(global_vars$intervals[this.period:(next.period-1)])
  }
}


## More initialization stuff:

init.z2 <- init.z
init.inpop.all1 <- array(1,dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions))
init.didmove <- array(0,dim=c(global_vars$nsites,global_vars$max_naug,global_vars$n_primary_occasions))
init.locationx <- caplocsx4d.init
init.locationy <- caplocsy4d.init
init.initx <- init.location.x
init.inity <- init.location.y

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  k=399
  for(k in 1:global_vars$max_naug){
    if(k>global_vars$naug[i]){
      init.z2[i,k] <- NA
      #init.didmove[i,k] <- NA
      init.inity[i,k] <- NA
      init.initx[i,k] <- NA
    }
    j=1
    for(j in 1:global_vars$n_primary_occasions){
      this.period <- site_data[[site]]$periods_surveyed[j]
      if((j>site_data[[site]]$nperiods)|(k>global_vars$naug[i])){
        init.inpop.all1[i,k,j] <- NA
        init.inpop[i,k,j] <- NA
        init.didmove[i,k,j] <- NA
      }
      l=1
      for(l in 1:global_vars$n_secondary_occasions){
        if((j>site_data[[site]]$nperiods)|(k>global_vars$naug[i])|(l>site_data[[site]]$replicate_sessions[this.period])){
          init.locationy[i,k,j,l] <- NA
          init.locationx[i,k,j,l] <- NA
        }
      }
    }
  }
}

init.z2[1,]
init.didmove[1,,]
init.inpop[5,,10]
init.inpop.all1[1,,10]
init.locationy[9,,1,1]
init.initx[1,]


###############
# BUGS MODEL
###############

params <- c(p0 = 0.5, sigma = 25)
halfnorm <- function(dist,params){
  params['p0'] * exp(-1*(dist^2)/(2*params['sigma']^2))
}

curve(halfnorm(x,params),0,75)  

##START THE FILE FOR JAGS! build function into bugs model

BUGSfilename <- "Uta_BUGS2.txt"

cat("
    
model{

#############
# PRIORS
#############

  ###########
  # PROBABILITY OF CAPTURE IN SURVEYED AREA
  ###########

  p0 ~ dunif(0,1)                 # mean/intercept detection prob
  logit.p0 <- log(p0/(1-p0)) 
  p.male.eff ~ dunif(-3,3)    # effect of sex on mean capture probability
  p.pm.eff ~ dunif(-2,2)    # effect of survey time on mean capture probability
  p.effort.eff ~ dunif(-1,3)  # effect of survey effort on mean capture probability
  p.wind.eff ~ dunif(-3,3)   # logit-linear effect of wind speed on capture probability
    #p.recap.eff ~ dunif(-1,5)  # effect of being after the first capture in a primary occasion. 
  p.vismark.eff ~ dunif(-1,5)  # effect of being permanently marked
  
  #### add random effect for period to soak up any additional variance in capture probability from session to session...

  # #### determine whether or not individual has been bead tagged...
  for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through (data augmented) individuals
      is.beadmarked[site,ind,1] <- 0    # no one is beadmarked in the first four time steps
      is.beadmarked[site,ind,2] <- 0
      is.beadmarked[site,ind,3] <- 0
      is.beadmarked[site,ind,4] <- 0
      for(period in 5:nperiods[site]){
        is.beadmarked[site,ind,period] <- step(captured_this_period[site,ind,(period-1)] + is.beadmarked[site, ind, (period-1)]-1)   # was it captured any time after the fourth period?
      }
    }
  }


    #### mean capture probability
  for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through (data augmented) individuals
      for(period in 1:nperiods[site]){
        for(session in 1:nsessions[site,period]){
          is.prevcapped[site,ind,period,session] <- step(session-firstsession[site,ind,period]-1)  # after the first capture in a given session, this should be 1 for each period
          is.vismarked[site,ind,period,session] <- step(is.beadmarked[site,ind,period]+is.prevcapped[site,ind,period,session]-1)  # does it have a visual mark (so it can be re-sighted)
          logit(thisp0[site,ind,period,session]) <- logit.p0 + p.male.eff*is.male[site,ind] + p.wind.eff*wind[site,period,session] +
                                                    p.effort.eff*effort[site,period,session] + p.pm.eff*is.pm[site,period,session]  +    #p.recap.eff*is.prevcapped[site,ind,period,session] +     # p.time.eff*(period-1) + 
                                                    p.vismark.eff*is.vismarked[site,ind,period,session] #+
                                                    
        }
      }
    }
  }

  # test.vismarked1 <- is.vismarked[1,2,3,1]   # should be 0
  # test.vismarked2 <- is.vismarked[1,2,4,1]   # should be 0
  # test.vismarked3 <- is.vismarked[1,2,4,2]   # should be 1
  # test.vismarked4 <- is.vismarked[1,2,5,1]   # should be 1


    
  ###########
  # PROBABILITY OF SURVIVING
  ###########
    
  phi0 ~ dunif(0.1,1)                 # mean/intercept survival (monthly)
  phi0.logit <- log(phi0/(1-phi0)) 

  phi.site.prec ~ dgamma(0.01,0.01)               # let mean survival differ by site
  phi.site.sd <- pow((1/phi.site.prec),0.5)
  for(site in 1:nsites){
    phi.site.eff[site] ~ dnorm(0,phi.site.prec)
  }

  for(site in 1:nsites){
    logit(phi[site]) <- phi0.logit + phi.site.eff[site] 
  }

  ###########
  # PROBABILITY OF RECRUITMENT
  ###########

  # Dirichlet prior for entrance probabilities (following Royle's S-A formulation)
  for(site in 1:nsites){
      init.entranceprob[site] ~ dgamma(5,1)       # ~ dunif(0.05,0.75)    #first entrance probability is fundamentally different from the others   # probability of entering the population at time 1 is estimated separately (I don't think it should be fixed... )
      beta[site,1] <- init.entranceprob[site]
      gamma[site,1] <- beta[site,1]/sum(beta[site,1:nperiods[site]])   # init.entranceprob[site]      #
    for(period in 2:nperiods[site]){
      beta[site,period] ~ dgamma(interval[site,(period-1)],4)    # <- interval[site,(period-1)]   #close to equal probability of entrance each occasion
      gamma[site,period] <- beta[site,period]/sum(beta[site,1:nperiods[site]])    # (interval[site,(period-1)]/sum(interval[site,1:(nperiods[site]-1)]))*(1-init.entranceprob[site])    # probability of entrance...
    }
  }
  
  # convert to conditional entrance probs

  for(site in 1:nsites){  
    cprob[site,1]<-gamma[site,1]
    for(period in 2:nperiods[site]){
      cprob[site,period]<-gamma[site,period]/(1-sum(gamma[site,1:(period-1)]))
    }
  }

  ###########
  # PROBABILITY OF BEING REAL
  ###########

  # psi0 ~ dunif(0,1)                # probability of being real- related to real densities
  # psi0_logit <- log(psi0/(1-psi0))
  # psi.windfarm.eff ~ dunif(-3,3)
  # psi.disturb.eff ~ dunif(-3,3)

  for(site in 1:nsites){
    psi[site] ~ dunif(0,1)      # estimate separately for every site? Alternatively, make it a function of disturbance level?  If it is on the same terms, then augmentation should be based on size
  }


  ###########
  # MOVEMENT PROCESS
  ###########

    ###########
    # SIZE OF ACTIVITY AREA (log scale)
    ###########
    
    
  # let mean activity area size differ by site
  for(site in 1:nsites){
    hr.site.log[site] ~ dunif(0,3)   # log scale... 
  }
  
  for(sex in 1:2){                
    hr.sex.eff[sex] ~ dunif(-1,1)   # effect of sex on size of activity area 
  }  
  
  for(site in 1:nsites){
    for(sex in 1:2){                # size of activity area- related to home range size, sex dependent
      log(sigma[site,sex]) <- hr.site.log[site] + hr.sex.eff[sex] 
      tau[site,sex] <- 1/pow(sigma[site,sex],2)
    }
  }
    
    ###########
    # SHIFT IN ACTIVITY AREA CENTER AMONG PERIODS
    ###########

  probmove0 ~ dunif(0,0.25)                # probability of moving per month
  probmove0_logit <- log(probmove0/(1-probmove0))
  m.male.eff ~ dunif(-3,3)

  m.site.prec ~ dgamma(1,1)               # let mean movement rate differ by site
  m.site.sd <- pow((1/m.site.prec),0.5)
  for(site in 1:nsites){
    m.site.eff[site] ~ dnorm(0,m.site.prec)
  }

    #m.disturb.eff ~ dunif(-3,3)
    
  for(site in 1:nsites){
    for(ind in 1:ninds[site]){
      logit(prob.move[site,ind]) <- probmove0_logit + m.site.eff[site] + m.male.eff*is.male[site,ind]       #m.windfarm.eff * is.windfarm[site]   # + psi.disturb.eff * disturbance[site] 
    }
  }

  for(site in 1:nsites){
    logit(prob.move_bysite[site,1]) <- probmove0_logit + m.site.eff[site] + m.male.eff     # assume males arbitrarily...
    logit(prob.move_bysite[site,2]) <- probmove0_logit + m.site.eff[site]      # (females)
  }

  for(site in 1:nsites){
    for(ind in 1:ninds[site]){
      initx[site,ind]  ~ dunif(xlim[site,1], xlim[site,2])   # specify location of center of activity at start of study    # provide inits? 
      inity[site,ind]  ~ dunif(ylim[site,1], ylim[site,2])
      hrx[site,ind,1] <- initx[site,ind] 
      hry[site,ind,1] <- inity[site,ind]
      did.move[site,ind,1] ~ dbern(0)     # can't move in the first 
      for(period in 2:nperiods[site]){
        probmove2[site,ind,period] <- (1-pow((1-prob.move[site,ind]),interval[site,(period-1)])) * in.pop.now[site,ind,(period-1)]
        did.move[site,ind,period] ~ dbern(probmove2[site,ind,period])     
        move.dist[site,ind,period] ~ dunif(4,20)        # movement distance
        move.dir[site,ind,period] ~ dunif(0,6.283185)    # movement direction
        hrx[site,ind,period] <- (1-did.move[site,ind,period])*hrx[site,ind,(period-1)]  + did.move[site,ind,period]*max(xlim[site,1],min(xlim[site,2],hrx[site,ind,(period-1)] + move.dist[site,ind,period] * cos(move.dir[site,ind,period])))
        hry[site,ind,period] <- (1-did.move[site,ind,period])*hry[site,ind,(period-1)]  + did.move[site,ind,period]*max(ylim[site,1],min(ylim[site,2],hry[site,ind,(period-1)] + move.dist[site,ind,period] * sin(move.dir[site,ind,period])))
      }
    }
    for(ind in (ninds[site]+1):maxinds){
      initx[site,ind]  <- 0   #~ dunif(xlim[site,1], xlim[site,2])   # for initialization purposes
      inity[site,ind]  <- 0   #  ~ dunif(ylim[site,1], ylim[site,2])
    }
  }


    ###########
    # SHIFT IN LOCATION AMONG SUBOCCASIONS WITHIN A PERIOD
    ###########

  #delta <- 10  

  for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through augmented individuals
      for(period in 1:nperiods[site]){
        for(session in 1:nsessions[site,period]){
          locationx[site,ind,period,session] ~ dnorm(hrx[site,ind,period],tau[site,thissex[site,ind]]) #I(xlim[site,1], xlim[site,2])    # data node    # note that the I() construct works for WinBUGS but not JAGS here
          locationy[site,ind,period,session] ~ dnorm(hry[site,ind,period],tau[site,thissex[site,ind]]) #I(ylim[site,1], ylim[site,2]) 
          xndx[site,ind,period,session] <- max(1,min(xlim[site,2],trunc(locationx[site,ind,period,session]+1)))    # max(1,min(xlim[site,2], ))
          yndx[site,ind,period,session] <- max(1,min(ylim[site,2],trunc(locationy[site,ind,period,session]+1)))
          inplot[site,ind,period,session] <- habmat[site, yndx[site,ind,period,session], xndx[site,ind,period,session]]  #~ dbern(     # habitat check, this line was added to see if individuals are located within the sampled habitat
        }
        for(session in (nsessions[site,period]+1):3){    # for initialization purposes...
          locationx[site,ind,period,session]  <- 0   # ~ dnorm(hrx[site,ind,period],tau[site,thissex[site,ind]]) #I(xlim[site,1], xlim[site,2])    # data node    # note that the I() construct works for WinBUGS but not JAGS here
          locationy[site,ind,period,session]  <- 0   # ~ dnorm(hry[site,ind,period],tau[site,thissex[site,ind]])
        }
      }
      for(period in (nperiods[site]+1):10){   # for initialization purposes
        for(session in 1:3){
          locationx[site,ind,period,session]  <- 0   # ~ dnorm(hrx[site,ind,1],tau[site,1])    # data node    # note that the I() construct works for WinBUGS but not JAGS here
          locationy[site,ind,period,session]   <- 0   #~ dnorm(hry[site,ind,1],tau[site,1])
        }
      }
    }
    for(ind in (ninds[site]+1):maxinds){
      for(period in 1:10){
        for(session in 1:3){
          locationx[site,ind,period,session]  <- 0   # ~ dnorm(hrx[site,1,1],tau[site,1])    # for initialization
          locationy[site,ind,period,session]  <- 0   # ~ dnorm(hry[site,1,1],tau[site,1])
        }
      }
    }
  }
  
#############
# DEAL WITH MISSING DATA
#############

  for(site in 1:nsites){
    prob.male[site] ~ dunif(0,1)    # allow to vary by site?
  }
    
  prob.pm ~ dunif(0,1)
    
  for(site in 1:nsites){
    for(ind in 1:maxinds){
      is.male[site,ind] ~ dbern(prob.male[site])
    }
  }

  for(site in 1:nsites){
    for(period in 1:nperiods[site]){
      for(session in 1:nsessions[site,period]){
       is.pm[site,period,session] ~ dbern(prob.pm)
      }
    }
  }


#############
# LIKELIHOOD
#############

  for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through augmented individuals
      thissex[site,ind] <- 2-is.male[site,ind]   # male is 1, female is 2
      z[site,ind] ~ dbern(psi[site])    # is this individual real or fake? (data augmentation)
      in.pop.now[site,ind,1] ~ dbern(gamma[site,1])   # initial entrance probability

      for(session in 1:nsessions[site,1]){
        p[site,ind,1,session] <- thisp0[site,ind,1,session] * z[site,ind] * in.pop.now[site,ind,1] * inplot[site,ind,1,session]    #probability of catching an individual at a given point given it is in the population, in the study area, and real
        y[site,ind,1,session] ~ dbern(p[site,ind,1,session])
      }
      
      recruitable[site,ind,1] <- 1-in.pop.now[site,ind,1]     # 1 if the indiv is not yet in the population 
      recruited[site,ind,1] <- in.pop.now[site,ind,1]

      for(period in 2:nperiods[site]){
        survival.part[site,ind,period] <- pow(phi[site],interval[site,(period-1)]) * in.pop.now[site,ind,(period-1)]
        recruitment.part[site,ind,period] <- cprob[site,period]*recruitable[site,ind,(period-1)]
        expected.inpop[site,ind,period] <- survival.part[site,ind,period] + recruitment.part[site,ind,period]  # either it is still in the pop or it just entered
        in.pop.now[site,ind,period] ~ dbern(expected.inpop[site,ind,period])
        recruitable[site,ind,period] <- recruitable[site,ind,(period-1)] * (1-in.pop.now[site,ind,period])       # is it still not yet in the study population?
        recruited[site,ind,period] <- (1-in.pop.now[site,ind,(period-1)]) * in.pop.now[site,ind,period]   # was it recruited this year?
        
        for(session in 1:nsessions[site,period]){
          p[site,ind,period,session] <- thisp0[site,ind,period,session] * z[site,ind] * in.pop.now[site,ind,period] * inplot[site,ind,period,session]   #probability of catching an individual at a given point given it is in the population, in the study area, and real
          y[site,ind,period,session] ~ dbern(p[site,ind,period,session])                           
        }
        captured_this_period[site,ind,period] <- step(sum(y[site,ind,period,1:nsessions[site,period]])-1)  # was it captured at least once?
      }    
    }
  }


  for(site in 1:nsites){
    for(ind in 1:ninds[site]){  # loop through augmented individuals
      for(period in (nperiods[site]+1):10){
        in.pop.now[site,ind,period]  <- 0   # ~ dbern(0)
        did.move[site,ind,period]  <- 0   # ~ dbern(0)
      }
    }
    for(ind in (ninds[site]+1):maxinds){  # loop through augmented individuals
      z[site,ind]  <- 0   # ~ dbern(0)
      for(period in 1:10){
        in.pop.now[site,ind,period]  <- 0   # ~ dbern(0)
        did.move[site,ind,period]  <- 0   # ~ dbern(0)
      }
    }
  }


#############
# DERIVED TERMS
#############

  for(site in 1:nsites){
    for(period in 1:nperiods[site]){
      N[site,period] <- inprod(in.pop.now[site,1:ninds[site],period],z[site,1:ninds[site]])
      N.recruited[site,period] <- inprod(recruited[site,1:ninds[site],period],z[site,1:ninds[site]])
      Density[site,period] <- N[site,period]/A[site]       #derive density
    }
  }


}   ## end BUGS model
    
",file=BUGSfilename)



###############
# SET UP WORKSPACE FOR BUGS
###############


# get areas

first.site <- 1 # 3 # 
last.site <-  global_vars$nsites # 4 # 3 # 
first.period <- 1
last.period <- 10 # 1 # 

nsites <- 1+last.site-first.site   #global_vars$nsites
nperiods <- pmin(global_vars$nperiods[first.site:last.site,drop=FALSE],rep(1+last.period-first.period,times=nsites))
maxperiods <- max(nperiods)
maxinds <- max(global_vars$naug[first.site:last.site,drop=FALSE])

data.for.bugs <- list(
  A=sapply(first.site:last.site,function(t) site_data[[t]]$area2,simplify = "array"), 
  y=caphist4d[first.site:last.site,1:maxinds,first.period:last.period,,drop=FALSE],
  locationx=caplocsx4d[first.site:last.site,1:maxinds,,,drop=FALSE],           # first.period:last.period
  locationy=caplocsy4d[first.site:last.site,1:maxinds,,,drop=FALSE],            # first.period:last.period
  ninds=global_vars$naug[first.site:last.site,drop=FALSE],
  maxinds=maxinds,
  habmat = habmat[first.site:last.site,,,drop=FALSE],
  xlim=t(sapply(first.site:last.site,function(t) site_data[[t]]$xlim2,simplify = "array")), 
  ylim=t(sapply(first.site:last.site,function(t) site_data[[t]]$ylim2,simplify = "array")), 
  is.male = is.male[first.site:last.site,1:maxinds,drop=FALSE],
  is.pm = is.pm[first.site:last.site,first.period:last.period,,drop=FALSE],
  effort = effort[first.site:last.site,first.period:last.period,,drop=FALSE],
  wind = wind[first.site:last.site,first.period:last.period,,drop=FALSE],
  nperiods = nperiods,
 # maxperiods = maxperiods,
  nsites = nsites,
  nsessions = sessions.surveyed[first.site:last.site,first.period:last.period,drop=FALSE],
  firstsession = firstsession[first.site:last.site,,first.period:last.period,drop=FALSE],
  #is.beadmarked = is.beadmarked[first.site:last.site,,first.period:last.period,drop=FALSE],
  interval=intervals[first.site:last.site,first.period:last.period,drop=FALSE]
  #is.windfarm = site_covariates$windpower,
  #initx = init.location.x[first.site:last.site,,drop=FALSE],
  #inity = init.location.y[first.site:last.site,,drop=FALSE]
)

initz.bugs<-function(){
  list(
       z=init.z2[first.site:last.site,1:maxinds,drop=FALSE],  #array(1,dim=c(nsites,global_vars$naug)),    #
       #psi=runif(nsites,0.4,0.5),
       p0=runif(1,0.1,0.15),
       p.male.eff=runif(1,-0.01,0.01),
       p.wind.eff=runif(1,-0.01,0.01),
       p.pm.eff=runif(1,-0.01,0.01),
       p.effort.eff=runif(1,0.01,0.05),
       #p.period.prec=runif(1,100,101),
       p.vismark.eff=runif(1,1,2),
       
       phi0=runif(1,0.6,0.75),
       phi.site.prec=20,
       phi.site.eff=runif(nsites,-1,1),
       
       #beta = runif(10,0.4,1),
       
       psi = runif(nsites,0.5,0.6),
       init.entranceprob = runif(nsites,3,4),
       
       hr.site.log=runif(nsites,1.5,2),
       hr.sex.eff=runif(2,-0.1,0.1),
       # sigma=runif(2,10,15),
       
       in.pop.now=init.inpop.all1[first.site:last.site,1:maxinds,,drop=FALSE],  #init.inpop[first.site:last.site,1:maxinds,,drop=FALSE],  #   array(1,dim=c(nsites,maxinds,10)), #  initialize every indiv as being in the population?  ,  #  
       did.move=init.didmove[first.site:last.site,1:maxinds,,drop=FALSE],  #array(0,dim=c(nsites,maxinds,10)),

       probmove0 = 0.01,
       m.site.prec=20,
       m.male.eff = 0,
       m.site.eff=runif(nsites,-1,1),
       
       locationx=init.locationx[first.site:last.site,1:maxinds,,,drop=FALSE],  #caplocsx4d.init[first.site:last.site,1:maxinds,,,drop=FALSE],      #first.period:last.period
       locationy=init.locationy[first.site:last.site,1:maxinds,,,drop=FALSE],  #caplocsy4d.init[first.site:last.site,1:maxinds,,,drop=FALSE],      #first.period:last.period
       
       initx = init.initx[first.site:last.site,1:maxinds,drop=FALSE], #init.location.x[first.site:last.site,1:maxinds,drop=FALSE],
       inity = init.inity[first.site:last.site,1:maxinds,drop=FALSE],  #init.location.y[first.site:last.site,1:maxinds,drop=FALSE],
       
       prob.male = rep(.5,times=nsites)
       #psi.windfarm.eff = 0,
       
       #m.windfarm.eff = 0,
      
       # init.location.x and init.location.y?
  )
}
#initz.bugs()


system.time(
  mod<-run.jags(
    model=BUGSfilename,
    monitor=c("p0","p.male.eff","p.vismark.eff","p.wind.eff","p.pm.eff","prob.male","p.effort.eff",   #"p.recap.eff",
                         "prob.move_bysite","m.male.eff","m.site.sd","probmove0","sigma","hr.site.log","hr.sex.eff",   # "p.period.sd","p.period.eff",
                         "phi0","phi.site.sd","phi",
                         "psi","Density","N",
                         "gamma","N.recruited","init.entranceprob"),
    data=data.for.bugs,
    n.chains = 3,
    inits=initz.bugs,
    burnin = 10000,
    #"test.vismarked1","test.vismarked2","test.vismarked3","test.vismarked4"),    # "p.time.eff","psi0" ,"prob.move","m.windfarm.eff","phi.windfarm.eff","psi.windfarm.eff",
    sample=10000,
    adapt=1000,
    thin=10,
    method="parallel"
    #clearWD=FALSE
  )
)

# system.time(
#   mod<-jags(
#     data=data.for.bugs,
#     inits=initz.bugs,
#     parameters.to.save=c("p0","p.male.eff","p.vismark.eff","p.wind.eff","p.pm.eff","prob.male","p.effort.eff",   #"p.recap.eff",
#                          "prob.move_bysite","m.male.eff","m.site.sd","probmove0","sigma","hr.site.log","hr.sex.eff",   # "p.period.sd","p.period.eff",
#                          "phi0","phi.site.sd","phi",
#                          "psi","Density","N",
#                          "gamma","N.recruited","init.entranceprob"),
#                          #"test.vismarked1","test.vismarked2","test.vismarked3","test.vismarked4"),    # "p.time.eff","psi0" ,"prob.move","m.windfarm.eff","phi.windfarm.eff","psi.windfarm.eff",
#     n.iter=5000,
#     model.file=BUGSfilename, 
#     n.chains = 2,
#     n.burnin = 2000,
#     n.thin=3
#     #clearWD=FALSE
#   )
# )

##### DEBUG

# site=1
# ind=30
# perndx=9    # [1,30,9]
# period=site_data[[site]]$periods_surveyed[perndx]
# rep=1
# habrast <- raster(habmat[site,1:site_data[[site]]$ylim2[2],1:site_data[[site]]$xlim2[2]],xmn=site_data[[site]]$xlim2[1],xmx=site_data[[site]]$xlim2[2],ymn=site_data[[site]]$ylim2[1],ymx=site_data[[site]]$ylim2[2])
# plot(habrast)
# res(habrast)
# site_data[[site]]$caphist.std.2d
# site_data[[site]]$ind_locs[[site_data[[site]]$uniqueinds[ind]]]
# caphist4d[site,ind,perndx,]
# caplocsy4d[site,ind,perndx,]
# points(caplocsx4d[site,ind,perndx,rep],site_data[[site]]$ylim2[2]-caplocsy4d[site,ind,perndx,rep],col="red",pch=20)
# points(caplocsx4d.init[site,ind,period,],site_data[[site]]$ylim2[2]-caplocsy4d.init[site,ind,period,],col="blue",pch=1)
# points(init.location.x[site,ind],site_data[[site]]$ylim2[2]-init.location.y[site,ind],pch=20)
# 
# xdx <- trunc(caplocsx4d[site,ind,perndx,rep]+1)
# ydx <- trunc(caplocsy4d[site,ind,perndx,rep]+1)
# habmat[site,ydx,xdx]        # node inconsistent with parents!!
# 
# points(xdx,ydx)


# winbugs.location <- "C:\\WinBUGS\\winbugs14\\WinBUGS14"
#
# library(R2WinBUGS)
# mod<-bugs(
#   data=data.for.bugs,
#   inits=initz.bugs,
#   parameters.to.save=c("p0","p.male.eff","p.vismark.eff","p.wind.eff","p.pm.eff","prob.male","p.effort.eff",   #"p.recap.eff",
#                        "prob.move_bysite","m.male.eff","m.site.sd","probmove0","sigma","hr.site.log","hr.sex.eff",   # "p.period.sd","p.period.eff",
#                        "phi0","phi.site.sd","phi",
#                        "psi","Density","N",
#                        "gamma","N.recruited"),    # "p.time.eff","psi0" ,"prob.move","m.windfarm.eff","phi.windfarm.eff","psi.windfarm.eff",
#   n.iter=1200,
#   model.file=BUGSfilename,
#   n.chains = 1,
#   n.burnin = 500,
#   n.thin=1,
#   bugs.directory = winbugs.location,
#   debug=TRUE
#   #clearWD=FALSE
# )

### Save results to disk

filename <- "bugsresults_ALLSITES_PRELIM11.RData"

save(list = ls(all.names = TRUE), file = filename, envir = .GlobalEnv)

?par

graphics.off()

Site_mcmc <- mod$mcmc  #as.mcmc(mod)   # does this work for "runjags"?

colnames(Site_mcmc[[1]])

plot(Site_mcmc[,"psi[1]"])
plot(Site_mcmc[,"psi[2]"])
plot(Site_mcmc[,"psi[3]"])

plot(Site_mcmc[,"init.entranceprob[1]"])
plot(Site_mcmc[,"p0"])
plot(Site_mcmc[,"p.male.eff"])
plot(Site_mcmc[,"p.wind.eff"])
plot(Site_mcmc[,"p.pm.eff"])
plot(Site_mcmc[,"p.effort.eff"])
plot(Site_mcmc[,"p.vismark.eff"])

#plot(Site_mcmc[,"p.period.sd"])
#plot(Site_mcmc[,"p.period.eff[1]"])

plot(Site_mcmc[,"phi0"])
plot(Site_mcmc[,"phi.site.sd"])
plot(Site_mcmc[,"phi[2]"])

plot(Site_mcmc[,"gamma[2,3]"])

plot(Site_mcmc[,"probmove0"])
plot(Site_mcmc[,"m.male.eff"])
plot(Site_mcmc[,"prob.move_bysite[1,2]"])
plot(Site_mcmc[,"sigma[1,1]"])
plot(Site_mcmc[,"sigma[2,1]"])


plot(Site_mcmc[,"Density[1,1]"])
plot(Site_mcmc[,"N[1,1]"])
plot(Site_mcmc[,"N[2,1]"])
plot(Site_mcmc[,"prob.male[1]"])


# Site_mcmc[,"test.vismarked1"]
# Site_mcmc[,"test.vismarked2"]
# Site_mcmc[,"test.vismarked3"]
# Site_mcmc[,"test.vismarked4"]


######
  ## store site-level response variables    [missing: J:A ratio?]


siteanalysis_df <- sitecovars_df

siteanalysis_df$survival <- NA 
siteanalysis_df$recruitment <- NA     # simply per-capita?
siteanalysis_df$sexratio <- NA
siteanalysis_df$meandens <- NA
siteanalysis_df$trenddens <- NA
siteanalysis_df$movedist <- NA
siteanalysis_df$moveprob <- NA


  ## determine order for presenting sites:
newsiteorder <- siteanalysis_df[order(siteanalysis_df$disturbance_rank,decreasing=T),]$newnames

neworder <- match(newsiteorder,siteanalysis_df$newnames)


#######
# PLOT DENSITY ETC
######

filename = sprintf("densities.pdf")
pdf(filename,6,6)
par(mfrow=c(3,3))
par(mai=c(0.4,0.4,0.1,0.1))

densresults <- list()

site=1
varnames <- colnames(Site_mcmc[[1]])
for(site in neworder){    # 1:global_vars$nsites
  sitename <- sitenames_df$newnames[site]
  #filename = sprintf("density_%s.pdf",sitename)
  periods <- site_data[[site]]$periods_surveyed
  nperiods <- site_data[[site]]$nperiods
  #pdf(filename,4,4)
  
  tomatch <- sprintf("Density[%i,%i]",site,1:nperiods)
  
  thissitedata <- Site_mcmc[[1]][,tomatch]
  
  densresults[[sitename]] <- thissitedata
  
  upperlower <- apply(thissitedata,2,function(t) stats::quantile(t,c(0.025,0.975)))*10000
  mean <- apply(thissitedata,2,mean)*10000
  
  errbar(x=global_vars$realdates[periods],y=mean,yplus=upperlower[2,],yminus=upperlower[1,],
         main=sitename,ylab="",xlab="",
         ylim=c(0,140),xlim=c(global_vars$firstdate,global_vars$firstdate+days(cumsum(trunc(global_vars$intervals*30))[9]))) #,xaxt="n")
  #axis(1,at=c(1:global_vars$n_primary_occasions),labels = c(1:global_vars$n_primary_occasions))
  legend("topleft",cex=1.2,legend=sitename,bg="grey")
  #title(main=sitename)
    
  #dev.off()
  
  siteanalysis_df$meandens[site] <- mean(mean)
  dayssince <- as.numeric(global_vars$realdates[periods]-global_vars$firstdate)
  trendmod <- lm(mean~dayssince)
  siteanalysis_df$trenddens[site] <- trendmod$coefficients[2]
  temp2 <- data.frame(dayssince=seq(0,max(dayssince),5))
  temp <- predict(trendmod,newdata=temp2,interval = c("confidence"),level=0.95)
  lines(global_vars$firstdate + days(temp2[,1]),temp[,1],col="gray",lty=1,lwd=2)
  lines(global_vars$firstdate + days(temp2[,1]),temp[,2],col="gray",lty=2)
  lines(global_vars$firstdate + days(temp2[,1]),temp[,3],col="gray",lty=2)
}

dev.off()


#######
# PLOT SEX RATIO
######

graphics.off()

varnames <- colnames(Site_mcmc[[1]])
filename = sprintf("sexratio.pdf",sitename)
pdf(filename,5,3.5)

quants <- data.frame(mean=numeric(global_vars$nsites),lower=0,upper=0)
site=1
for(site in neworder){  # 1:global_vars$nsites
  sitename <- sitenames_df$newnames[site]
  tomatch <- sprintf("prob.male[%i]",site)
  thissitedata <- Site_mcmc[[1]][,tomatch]
  
  quants$mean[site] <- mean(thissitedata)
  quants[site,c(2,3)] <- quantile(thissitedata,c(0.025,0.975))
  
}

errbar(x=1:global_vars$nsites,y=quants$mean,yplus=quants$upper,yminus=quants$lower,main="",ylab="Sex ratio",xlab="Site",xaxt="n",pch=site_covariates$windpower[neworder]+1)
axis(1,at=c(1:global_vars$nsites),labels = sitenames_df$newnames[neworder])

siteanalysis_df$sexratio <- quants$mean

dev.off()





#######
# PLOT SURVIVAL
######

graphics.off()
par(mfrow=c(2,2))

varnames <- colnames(Site_mcmc[[1]])

quants <- data.frame(mean=numeric(global_vars$nsites),lower=0,upper=0)
site=1
for(site in neworder){
  sitename <- sitenames_df$newnames[site]
  tomatch <- sprintf("phi[%i]",site)
  thissitedata <- Site_mcmc[[1]][,tomatch]
  
  quants$mean[site] <- mean(thissitedata)
  quants[site,c(2,3)] <- quantile(thissitedata,c(0.025,0.975))
  
}

filename = sprintf("survival.pdf",sitename)
pdf(filename,5,3.5)
errbar(x=1:global_vars$nsites,y=quants$mean,yplus=quants$upper,yminus=quants$lower,main="monthly survival",ylab="Mean survival rate",xlab="Site",xaxt="n",pch=site_covariates$windpower[neworder]+1)
title("monthly survival")
axis(1,at=c(1:global_vars$nsites),labels = sitenames_df$newnames[neworder])

dev.off()

siteanalysis_df$survival <- quants$mean



#######
# PLOT ANNUAL RECRUITMENT
######

varnames <- colnames(Site_mcmc[[1]])


quants <- data.frame(mean=numeric(global_vars$nsites),lower=0,upper=0)
quants2 <- data.frame(mean=numeric(global_vars$nsites),lower=0,upper=0)
site=1
for(site in neworder){
  sitename <- sitenames_df$newnames[site]
  periods <- site_data[[site]]$nperiods
  tomatch <- sprintf("N.recruited[%i,%i]",site,2:periods)
    # tomatch2 <- sprintf("N[%i,%i]",site,1:periods)   # recruited per adult or per area?
  thisrecruited <- Site_mcmc[[1]][,tomatch]
  thisrecruited_area <- thisrecruited / (site_data[[site]]$area2/10000)
  thisrecruited_percap <- thisrecruited / (densresults[[sitename]][,-1]*site_data[[site]]$area2)  # per-capita recruitment
  
  thissitedata <- apply(as.matrix(thisrecruited_area),1,mean)
  quants$mean[site] <- mean(thissitedata)
  quants[site,c(2,3)] <- quantile(thissitedata,c(0.025,0.975))
  
  thissitedata2 <- apply(as.matrix(thisrecruited_percap),1,mean)
  quants2$mean[site] <- mean(thissitedata2)
  quants2[site,c(2,3)] <- quantile(thissitedata2,c(0.025,0.975))
}

filename = sprintf("Recruitment_perha.pdf",sitename)
pdf(filename,5,3.5)
errbar(x=1:global_vars$nsites,y=quants$mean,yplus=quants$upper,yminus=quants$lower,main="",ylab="Mean recruitment per ha",xlab="Site",xaxt="n",pch=site_covariates$windpower[neworder]+1)

axis(1,at=c(1:global_vars$nsites),labels = sitenames_df$newnames[neworder])

dev.off()


filename = sprintf("Recruitment_percap.pdf",sitename)
pdf(filename,5,3.5)
errbar(x=1:global_vars$nsites,y=quants2$mean,yplus=quants2$upper,yminus=quants2$lower,main="",ylab="Mean per-capita recruitment",xlab="Site",xaxt="n",pch=site_covariates$windpower[neworder]+1)
title("per-capita recruitment")
axis(1,at=c(1:global_vars$nsites),labels = sitenames_df$newnames[neworder])

dev.off()

siteanalysis_df$recruitment <- quants2$mean

#######
# PLOT MOVEMENT RATE
######

varnames <- colnames(Site_mcmc[[1]])
filename = sprintf("Movement.pdf",sitename)
pdf(filename,5,3.5)

quants <- data.frame(mean=numeric(global_vars$nsites),lower=0,upper=0)
site=1
for(site in neworder){
  sitename <- sitenames_df$newnames[site]
  tomatch <- sprintf("prob.move_bysite[%i,2]",site)
  thissitedata <- Site_mcmc[[1]][,tomatch]
  #thissitedata <- plogis(thissitedata)
  
  quants$mean[site] <- mean(thissitedata)
  quants[site,c(2,3)] <- quantile(thissitedata,c(0.025,0.975))
  
}

errbar(x=1:global_vars$nsites,y=quants$mean,yplus=quants$upper,yminus=quants$lower,main="",ylab="Female movement probability",xlab="Site",xaxt="n",pch=site_covariates$windpower[neworder]+1)
title("movement probability")
axis(1,at=c(1:global_vars$nsites),labels = sitenames_df$newnames[neworder])

dev.off()

siteanalysis_df$moveprob <- quants$mean

#######
# PLOT MOVEMENT DISTANCE
######

varnames <- colnames(Site_mcmc[[1]])
filename = sprintf("MovementDist.pdf",sitename)
pdf(filename,5,3.5)

quants <- data.frame(mean=numeric(global_vars$nsites),lower=0,upper=0)
site=1
for(site in neworder){
  sitename <- sitenames_df$newnames[site]
  tomatch <- sprintf("sigma[%i,1]",site)
  thissitedata <- Site_mcmc[[1]][,tomatch]
  
  quants$mean[site] <- mean(thissitedata)
  quants[site,c(2,3)] <- quantile(thissitedata,c(0.025,0.975))
  
}

errbar(x=1:global_vars$nsites,y=quants$mean,yplus=quants$upper,yminus=quants$lower,main="",ylab="Male movement distance \n(stdev from activity center)",xlab="Site",xaxt="n",pch=site_covariates$windpower[neworder]+1)
title("size of activity area")
axis(1,at=c(1:global_vars$nsites),labels = sitenames_df$newnames[neworder])

dev.off()

siteanalysis_df$movedist <- quants$mean

#####################
# GET ADDITIONAL RESPONSE VARIABLES
#####################

##### 
#  develop body condition model
#####

conditiondf <- data.frame(
  sex = character(0),
  svl = numeric(0),
  weight = numeric(0) 
)

site=2
for(site in 1:global_vars$nsites){
  sitename <- sitenames_df$newnames[site] 
  tempsex <- rep(site_data[[sitename]]$sex,times=ncol(site_data[[sitename]]$svl))
  tempsvl <- as.vector(site_data[[sitename]]$svl)
  tempweight <- as.vector(site_data[[sitename]]$weight)
  temp <- data.frame(
    sex = tempsex,
    svl = tempsvl,
    weight = tempweight 
  )
  conditiondf <- rbind(conditiondf,temp)
}

conditiondf <- na.omit(conditiondf)

conditiondf

conditiondf$male <- conditiondf$sex=="m"

conditiondf$is.male <- ifelse(conditiondf$sex=="m",1,0)

nrow(conditiondf)



malemodel <- loess(weight~svl, data=conditiondf[conditiondf$male,],span=1.5)
femalemodel <- loess(weight~svl, data=conditiondf[!conditiondf$male,],span=1.5)

plot(conditiondf$weight~conditiondf$svl,pch=conditiondf$is.male+1,ylab="Body mass (g)",xlab="SVL")

ndx <- order(conditiondf$svl[conditiondf$male])
lines(conditiondf$svl[conditiondf$male][ndx], predict(malemodel)[ndx], col = "blue",lwd=2)

ndx <- order(conditiondf$svl[!conditiondf$male])
lines(conditiondf$svl[!conditiondf$male][ndx], predict(femalemodel)[ndx], col = "red",lwd=2)

conditionfunc <- function(sex=thissex,svl=thissvl,weight=thisweight){
  newdata <- na.omit(data.frame(
    sex=sex,
    svl=svl,
    weight=weight,
    stringsAsFactors = F
  ))
  ifelse(newdata$sex=="m",newdata$weight-predict(malemodel,data.frame(svl=newdata$svl)),newdata$weight-predict(femalemodel,data.frame(svl=newdata$svl)))
}


lnth <- global_vars$nsites   #*global_vars$n_primary_occasions
sitedf2 <- data.frame(
  site = character(lnth),
  #year = numeric(lnth),
  meansvl = numeric(lnth),
  ageratio = numeric(lnth),
  meancondition = numeric(lnth),
  stringsAsFactors = F
)

counter <- 1
site=4
for(site in 1:global_vars$nsites){
  sitename <- sitenames_df$newnames[site] 
  thissex <- site_data[[site]]$sex
  thisage <- as.vector(site_data[[site]]$age)
  thissvl <- as.vector(site_data[[site]]$svl)
  thisweight <- as.vector(site_data[[site]]$weight)
  
  thisratio <- length(which(thisage=="1")) / length(which(!is.na(thisage)))     # percent subadult
  
  thiscondition <- conditionfunc(rep(thissex,times=global_vars$n_primary_occasions),thissvl,thisweight)

  sitedf2$site[site] <- sitename
  sitedf2$meansvl[site] <- mean(thissvl, na.rm=T)
  sitedf2$ageratio[site] <- thisratio 
  sitedf2$meancondition[site] <- mean(thiscondition)
}

sitedf2

#####################
# EVALUATE SITE COVARIATES
#####################

siteanalysis_df

siteanalysis_df <- cbind(siteanalysis_df,sitedf2[,-1])

responsevars <- c(
  "survival",
  "recruitment",
  "sexratio",
  "meandens",
  "trenddens",
  "movedist",
  "moveprob",
  "meansvl",
  "ageratio",
  "meancondition"
)    #names(siteanalysis_df)[15:21]

predictorvars <- c(
  "windpower",
   #"disturbance_rank",
  "disturbance",
  "bareground",
   #"noiselevel",
   #"shrubdensity",
   #"shrubcanopy",
   #"shrubheight",     
   #"slope",
   #"elevation",
  "rainfall"
  #"avgwind"
) #names(siteanalysis_df)[3:14]

responsevars

predictorvars

predform <- paste(predictorvars,collapse="+")


  ### correlation analysis

cor(siteanalysis_df[,predictorvars][,-1])


 ### run the analysis

bestmodels <- list()

siteanalysis_df$windpower <- as.factor(siteanalysis_df$windpower)

transform <- c("logit","log","logit","none","none","none","logit","none","none","none")   # response variable transformations
names(transform) <- responsevars 

response <- responsevars[9]
for(response in responsevars){
  if(transform[response]%in%c("log","logit")){
    respname <- sprintf("%s_%s",response,transform[response])
  }else{
    respname <- response
  }
  if(transform[response]=="logit"){
    siteanalysis_df[,respname] <- qlogis(siteanalysis_df[,response])
  }
  if(transform[response]=="log"){
    siteanalysis_df[,respname] <- log(siteanalysis_df[,response])
  }
  formula <- as.formula(sprintf("%s ~ %s",respname,predform))
  fullmodel <- lm(formula,data=siteanalysis_df)
    # summary(fullmodel)
  
  bestmodel <- step(fullmodel)
  
  bestmodels[[response]] <- bestmodel

}


summary(bestmodels[["survival"]])    # survival: openness, rainfall and wind are all important!
plot(bestmodels[["survival"]])
confint(bestmodels[["survival"]])

summary(bestmodels[["recruitment"]])    # recruitment: weak effect of rainfall?
plot(bestmodels[["recruitment"]])
confint(bestmodels[["recruitment"]])

summary(bestmodels[["sexratio"]])  # sex ratio: affected by disturbance, bareground, and rainfall? 
plot(bestmodels[["sexratio"]])
confint(bestmodels[["sexratio"]])

summary(bestmodels[["meandens"]])  # mean density: not affected by anything... 
plot(bestmodels[["meandens"]])
confint(bestmodels[["meandens"]])

summary(bestmodels[["trenddens"]])  # density trend: weakly affected by disturbance? More disturbed sites exhibited a weak decline ...
plot(bestmodels[["trenddens"]])
confint(bestmodels[["trenddens"]])

summary(bestmodels[["movedist"]])  # movement distance: stronly negatively affected by bare ground percentage
plot(bestmodels[["movedist"]])
confint(bestmodels[["movedist"]])

summary(bestmodels[["moveprob"]])  # movement probability: affected by bareground
plot(bestmodels[["moveprob"]])
confint(bestmodels[["moveprob"]])

summary(bestmodels[["meansvl"]])  # svl: affected by disturbance and rainfall
plot(bestmodels[["meansvl"]])
confint(bestmodels[["meansvl"]])

summary(bestmodels[["ageratio"]])  # age ratio: affected by disturbance
plot(bestmodels[["ageratio"]])
confint(bestmodels[["ageratio"]])

summary(bestmodels[["meancondition"]])  # body condition: affected by disturbance
plot(bestmodels[["meancondition"]])
confint(bestmodels[["meancondition"]])

######################
# Visualize relationships

response="survival";predictor="windpower"

response="sexratio";predictor="bareground"

visualizeRelation <- function(response="survival",predictor="bareground"){
  model <- bestmodels[[response]]
  pred <- siteanalysis_df[,predictor]
  otherpreds <- names(model$model)[-1][names(model$model)[-1]!=predictor]
  fac <- ifelse(length(unique(pred))<=4,T,F)
  if(fac){
    span <- unique(pred)
    newdata=data.frame(temp=numeric(length(span)))
    newdata[,predictor] <- as.factor(span)
    span <- as.factor(span)
  }else{
    toadd <- abs(max(pred) - min(pred))/10 
    span <- seq(min(pred)-toadd,max(pred)+toadd,length=100)
    newdata=data.frame(temp=numeric(length(span)))
    newdata[,predictor] <- span
  }
  other <- otherpreds[2]
  for(other in otherpreds){
    pred2 <- siteanalysis_df[,other]
    fac2 <- ifelse(length(unique(pred2))<=4,T,F)
    if(fac2){
      val <- as.numeric(names(table(pred2)))[which.max(table(pred2))]
      newdata[,other] <- factor(val,levels=names(table(pred2)))
    }else{
      newdata[,other] <- mean(siteanalysis_df[,other])
    }
  }
  temp <- predict(model,newdata=newdata,interval = c("confidence"),level=0.95)
  if(transform[response]=="logit"){
    temp <- plogis(temp)
  }
  if(transform[response]=="log"){
    temp <- exp(temp)
  }
  if(!fac){
    toaddy <- abs(max(siteanalysis_df[,response]) - min(siteanalysis_df[,response]))/10
    plot(temp[,1]~span,type="l",lty=1,lwd=3,ylab=response,xlab=predictor,
         ylim=c(min(siteanalysis_df[,response])-toaddy,max(siteanalysis_df[,response])+toaddy))
    lines(span,temp[,2],lty=2,col="gray",lwd=2)
    lines(span,temp[,3],lty=2,col="gray",lwd=2)
    text(siteanalysis_df[,predictor],siteanalysis_df[,response],siteanalysis_df$newnames)
  }else{
    errbar(x=as.numeric(span),y=temp[,1],yplus=temp[,3],yminus=temp[,2],
           ylab=response,xlab=predictor,xlim=c(0.75,max(as.numeric(span)+0.25)),xaxt="n",pch="____",cex=5)
    axis(1,at=as.numeric(span),levels(span))
  }
}


graphics.off()
pdf("survcovars.pdf",3.5,5.5)
par(mfrow=c(3,1))
par(mai=c(0.6,0.6,0.05,0.05))

visualizeRelation("survival","bareground")

visualizeRelation("survival","windpower")

visualizeRelation("survival","rainfall")

 # visualizeRelation("survival","disturbance")
dev.off()


graphics.off()
pdf("disturbance_effects.pdf",7,4)
par(mfrow=c(2,2))
par(mai=c(0.7,0.7,0.1,0.2))

#visualizeRelation("recruitment","rainfall")
 # visualizeRelation("recruitment","disturbance")

#visualizeRelation("sexratio","disturbance")

#visualizeRelation("sexratio","rainfall")

#visualizeRelation("sexratio","bareground")

visualizeRelation("trenddens","disturbance")

#visualizeRelation("movedist","bareground")

visualizeRelation("movedist","disturbance")


#visualizeRelation("moveprob","bareground")

#visualizeRelation("meansvl","disturbance")

#visualizeRelation("meansvl","rainfall")

visualizeRelation("ageratio","disturbance")

visualizeRelation("meancondition","disturbance")

dev.off()


graphics.off()
pdf("disturbance_dens.pdf",5,3.5)
visualizeRelation("trenddens","disturbance")
dev.off()


#visualizeRelation("movedist","bareground")

visualizeRelation("movedist","disturbance")


#visualizeRelation("moveprob","bareground")

#visualizeRelation("meansvl","disturbance")

#visualizeRelation("meansvl","rainfall")

visualizeRelation("ageratio","disturbance")

visualizeRelation("meancondition","disturbance")

dev.off()


####################
# Make parameter table
####################

colnames(Site_mcmc[[1]])

  # parameter table should have all key free parameters, with both prior and posterior indicated


param_table <- data.frame(parameter = character(0),prior = character(0), posterior=character(0))

  # mean detection probability

name <- "p0"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
                     parameter=name,
                     prior="Uniform(0,1)",
                     posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
                     ))

name <- "p.male.eff"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Uniform(-3,3)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))

name <- "p.pm.eff"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Uniform(-3,3)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))

name <- "p.effort.eff"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Uniform(-1,3)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))


name <- "p.wind.eff"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Uniform(-3,3)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))

name <- "p.vismark.eff"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Uniform(-3,3)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))

name <- "p.vismark.eff"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Uniform(-3,3)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))


name <- "phi0"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Uniform(0,1)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))


name <- "phi.site.prec"
name2 <- "phi.site.sd"
convert <- 1/Site_mcmc[[1]][,name2,drop=FALSE]^2
mean <- apply(convert,2,mean)
upper <- apply(convert,2,function(t) quantile(t,0.975))
lower <- apply(convert,2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Gamma(0.01,0.01)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))

name <- "psi"
names <- colnames(Site_mcmc[[1]])[grep(name,colnames(Site_mcmc[[1]]))]
mean <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=names,
  prior="Uniform(0,1)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))

name <- "hr.site.log"
names <- colnames(Site_mcmc[[1]])[grep(name,colnames(Site_mcmc[[1]]))]
mean <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=names,
  prior="Uniform(0,3)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))

name <- "hr.sex.eff"     # remember, male is 1, female is 2
names <- colnames(Site_mcmc[[1]])[grep(name,colnames(Site_mcmc[[1]]))]
mean <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=names,
  prior="Uniform(-1,1)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))

name <- "probmove0"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Uniform(0,1)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))


name <- "m.male.eff"
mean <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,name,drop=FALSE],2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Uniform(-3,3)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))


name <- "m.site.prec"
name2 <- "m.site.sd"
convert <- 1/Site_mcmc[[1]][,name2,drop=FALSE]^2
mean <- apply(convert,2,mean)
upper <- apply(convert,2,function(t) quantile(t,0.975))
lower <- apply(convert,2,function(t) quantile(t,0.025))
param_table <- rbind(param_table,data.frame(
  parameter=name,
  prior="Gamma(1,1)",
  posterior=sprintf("%1.2f (%1.2f to %1.2f)",mean,lower,upper)
))


getwd()
write.csv(param_table,file = "Parameter.table.csv")


#####################
# Write recruitment table (table of entrance probabilities)
#####################

temp <- matrix("",nrow=9,ncol=10)

sitenames <- sitenames_df$newnames
periodnames <- c(1:global_vars$n_primary_occasions)

name <- "gamma"
names <- colnames(Site_mcmc[[1]])[grep(name,colnames(Site_mcmc[[1]]))]
mean <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,mean)
upper <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.975))
lower <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.025))

 #Site_mcmc[[1]][,names,drop=FALSE] 

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  ndx <- grep(sprintf("%i,",i),names)
  names[ndx]
  temp[i,site_data[[i]]$periods_surveyed] <- sprintf("%1.2f (%1.2f to %1.2f)",mean[ndx],lower[ndx],upper[ndx])
}

colnames(temp) <- periodnames
rownames(temp) <- sitenames

getwd()
write.csv(temp,file = "Recruitment_table.csv",row.names = T)


#####################
# Write density table (table of densities)
#####################

temp <- matrix("",nrow=9,ncol=10)

sitenames <- sitenames_df$newnames
periodnames <- c(1:global_vars$n_primary_occasions)

name <- "Density"
names <- colnames(Site_mcmc[[1]])[grep(name,colnames(Site_mcmc[[1]]))]
mean <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,mean)*10000
upper <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.975))*10000
lower <- apply(Site_mcmc[[1]][,names,drop=FALSE],2,function(t) quantile(t,0.025))*10000

#Site_mcmc[[1]][,names,drop=FALSE] 

i=1
for(i in 1:global_vars$nsites){
  site = sitenames_df$newnames[i]
  ndx <- grep(sprintf("%i,",i),names)
  names[ndx]
  temp[i,site_data[[i]]$periods_surveyed] <- sprintf("%1.2f (%1.2f to %1.2f)",mean[ndx],lower[ndx],upper[ndx])
}

colnames(temp) <- periodnames
rownames(temp) <- sitenames

getwd()
write.csv(temp,file = "Density_table.csv",row.names = T)



############
# MISC
############

## compute annual survival


siteanalysis_df$survival^12




