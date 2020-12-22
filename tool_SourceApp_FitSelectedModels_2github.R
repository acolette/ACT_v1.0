rm(list = ls(all = TRUE))

library(ncdf4)
library(fields)

myargs <- commandArgs()
pol <- as.character(myargs[length(myargs)-1]) # "PM10"
day1    <- as.character(myargs[length(myargs)]) # "20161201"

# chimere output files
indir="/ccc/store/cont004/ineris/acolette/Scen_CAMS/OUTPUT/"
# surrogate fits
outdir="/ccc/scratch/cont004/ineris/acolette/Scen_CAMS/POSTPROC/tool_SourceApp/"

# function for 2D surrogate fitting
fit_a_a2_i_i2_r_t_t2_ai_ta_ti_tari=function(call_array) {
pol2d_ref <- call_array[1]
pol2d_sce <- call_array[2:12]
remi_sce <- call_array[13:29]
dz <- c(0,pol2d_ref - pol2d_sce[1],pol2d_ref - pol2d_sce[2],pol2d_ref - pol2d_sce[3],
          pol2d_ref - pol2d_sce[4],pol2d_ref - pol2d_sce[5],pol2d_ref - pol2d_sce[6],
          pol2d_ref - pol2d_sce[7],pol2d_ref - pol2d_sce[8],pol2d_ref - pol2d_sce[9],
          pol2d_ref - pol2d_sce[10],pol2d_ref - pol2d_sce[11])
da <- c(0,remi_sce[1],remi_sce[2],          0,          0,          0,          0,          0,remi_sce[8] ,remi_sce[11] ,           0,remi_sce[14])
di <- c(0,          0,          0,remi_sce[3],remi_sce[4],          0,          0,          0,remi_sce[9] ,          0  ,remi_sce[13],remi_sce[15])
dr <- c(0,          0,          0,          0,          0,remi_sce[5],          0,          0,           0,          0  ,           0,remi_sce[16])
dt <- c(0,          0,          0,          0,          0,          0,remi_sce[6],remi_sce[7],           0,remi_sce[10] ,remi_sce[12],remi_sce[17])
data = data.frame(a=da,i=di,r=dr,t=dt,z=dz)
model <- lm(dz ~ da + I(da^2) + di + I(di^2) + dr + dt + I(dt^2) + da:di + dt:da + dt:di + dt:da:dr:di - 1, data=data)
return(c(model[[1]]))}

# get grid
fic <- nc_open(paste(indir,"REF/outl.20161122_20161122_REF.nc",sep=""))
lat <- ncvar_get(fic,"lat")
lon <- ncvar_get(fic,"lon")
nc_close(fic)
nlon <- length(lon[,1])
nlat <- length(lat[1,])

# input scenarios from Chimere runs
sce_list <- c("REF","AGR60unif","AGR100unif","IND60unif","IND100unif","RH100unif","TRA60unif","TRA100unif",
              "AGR30IND60unif","TRA100AGR100unif","TRA30IND60unif","ALL100unif")

if (pol=="O3"){unit="ppb"}
if ((pol=="PM10")|(pol=="PM25")){unit="ug/m3"}

model <- list()
mod_list <- list()
dd <- 1

day2 <- day1
pol2d <- list()
#######################################################################################################
#######################################################################################################
for (sce in sce_list){
print(paste(indir,"/",sce,"/outl.",day1,"_",day2,"_",sce,".nc",sep=""))
fic <- nc_open(paste(indir,"/",sce,"/outl.",day1,"_",day2,"_",sce,".nc",sep=""))
if (pol=="O3"){
pol2d[[sce]] <- apply(ncvar_get(fic,pol,start=c(1,1,1,1),count=c(-1,-1,1,24)),1:2,max)
}
if ((pol=="PM10")|(pol=="PM25")) {
pol2d[[sce]] <- apply(ncvar_get(fic,pol,start=c(1,1,1,1),count=c(-1,-1,1,24)),1:2,mean)
}
nc_close(fic)

if (pol== "O3"){pol2d[[sce]] <- pol2d[[sce]] * 2 }

} # for sce

call_array <- numeric(nlon*nlat*29)
dim(call_array) <- c(nlon,nlat,29)
call_array[,,1] <- pol2d[["REF"]]
for (ss in c(2:(length(sce_list)))){call_array[,,ss] <- pol2d[[sce_list[ss]]]}
call_array[,,13] <- .6
call_array[,,14] <- 1.
call_array[,,15] <- .6
call_array[,,16] <- 1.
call_array[,,17] <- 1.
call_array[,,18] <- .6
call_array[,,19] <- 1.
call_array[,,20] <- .3
call_array[,,21] <- .6
call_array[,,22] <- 1.
call_array[,,23] <- 1.
call_array[,,24] <- .3
call_array[,,25] <- .6
call_array[,,26] <- 1.
call_array[,,27] <- 1.
call_array[,,28] <- 1.
call_array[,,29] <- 1.

model <- apply(call_array,1:2,fit_a_a2_i_i2_r_t_t2_ai_ta_ti_tari)

conc <- pol2d[["REF"]]
lon <- lon[,1]
lat <- lat[1,]
save(lon,lat,conc,model,file=paste(outdir,"models_FitModel2Dinter100_",pol,"_",day1,".dat",sep=""))

