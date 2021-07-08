#library("cytofCore", lib.loc="~/R/4.0.2")
#library("openxlsx", lib.loc="/opt/R/4.0.2/lib/R/library")
library("cytofCore")
library("openxlsx")
library("data.table")
library("dtwclust")
## Read imd file
setwd("~/CyTOF/PSH_data1")
#test <- cytofCore.read.imd("20200429_RPMI_FBS_0_1_01.imd", conf = "conf/Duals.conf", start_push = 0, num_pushes = NULL)
test <- cytofCore.read.imd("20200429_RPMI_FBS_0_1_01.imd", start_push = 0, num_pushes = 2^25)

# extract dual count data
test1 <- as.data.frame(test$dual)

# extract pulse data
#test2 <- as.data.frame(test$pulse)

# extract intensity data
#test3 <- as.data.frame(test$intensity)

# extract Pt and Au data
AuPt <- data.frame("Viability_Pt" = test1$`195Pt_Viability`,
                      "AuNP"      = test1$`197Au_Gold`)
# Save data for later use
# save(AuPt, file = "20200429_RPMI_FBS_0_1_01_AuPt.Rdata")

## Create a new workbook and add a worksheet
ExcelFile <- createWorkbook("CyTOF_Test_1")
addWorksheet(ExcelFile, sheetName = "Dual")
writeData(ExcelFile, sheet = 1, test1_1, rowNames = FALSE)
#addWorksheet(ExcelFile, sheetName = "Pulse")
#writeData(ExcelFile, sheet = 2, test2, rowNames = FALSE)
#addWorksheet(ExcelFile, sheetName = "Intensity")
#writeData(ExcelFile, sheet = 3, test3, rowNames = FALSE)
## Save workbook to working directory
saveWorkbook(ExcelFile, file = "20200429_RPMI_FBS_0_1_01.xlsx", overwrite = TRUE)
## End

#cytofCore.read.conf("conf/LunaTuning.conf")

## Extract Au peaks
data <- data.frame("Peak.no" = 0,
                   "Au_Signal" = AuPt$AuNP)

# build a function named "PeakAsign" to assign peak number for spICPMS Signal, the function will collect all Signal greater than t
PeakAssign <- function(x,t){
  ifelse(y <- x>t, cumsum(c(head(y, 1), tail(y, -1) - head(y, -1) == 1)), NA)  
}
# Apply function PeakAsign to collect non-zero data
data$Peak.no <- PeakAssign(data$Au_Signal,0)
data <- subset(data, Peak.no>0)
# Save data for later use
# save(data, file = "20200429_RPMI_FBS_0_1_01_Au.Rdata")


# Collect Peak Width, Peak height, Peak Area

# Count how many dwell time in one peak a.k.a. Peak.Width. "FUN=function(x){NROW(x)}" is same for every data
PW <- aggregate(data$Peak.no, by=list(Peak.no=data$Peak.no), FUN=function(x){NROW(x)})
# Name 2 columns of Peak.Width
colnames(PW) <- c("Peak.no", "Peak.Width")
DwellTime <- 13
PW$Event.Duration.us <- PW$Peak.Width * DwellTime

# Make row number of Peak.Width equal to row number of raw data and sign it as PW
#Peak.Width <- PW[rep(row.names(PW), PW$Peak.Width), 1:2]
# Numbering each dwel time in one peak
#data$Dwell.time.i <- ave(data[,1], data$Peak.no, FUN = seq_along) 
# Add Peak.Width column to raw data
#data$Peak.Width <- PW[,2]



# Calculate Peak Area (PA) for each peak group
PA <- aggregate(data$Au_Signal, by=list(Peak.no=data$Peak.no), FUN=sum)
# Name the PA table
colnames(PA) <- c("Peak.no", "Peak.Area")

# Extract peak height (PH) which is max of one peak group  
PH <- aggregate(data$Au_Signal, by=list(Peak.no=data$Peak.no), FUN=max)
colnames(PH) <- c("Peak.no", "Peak.Max")

# Combine columns Peak.no, Peak width, Peak Area and Peak Max
Table.1 <- data.frame("Peak.no" = PW$Peak.no,
                      "Peak.Width" = PW$Peak.Width,
                      "Event.duration" = PW$Event.Duration.us,
                      "Peak.Area" = PA$Peak.Area,
                      "Peak.Max" = PH$Peak.Max)
# Save data for later use
# save(Table.1, file = "20200429_RPMI_FBS_0_1_01_Table1.Rdata")

#=====================================================================================================================
## Collect peak content. Each Peak has a certain length which is number of dwell times. In each dwell time, there is signal corresponding to that dwell time denoted as DT0, DT1, DT2...  This step will collect those signals.

# Collect every peak to "PeakList", each Peak is a vector, PeakList is a list of (Peak) vectors
Table.2 <- Table.1[which(Table.1$Peak.Width > 20 & Table.1$Peak.Width < 150 & Table.1$Peak.Area>400), ]
# Save data for later use
# save(Table.2, file = "20200429_RPMI_FBS_0_1_01_Table2.Rdata")


PeakList <- lapply(Table.2$Peak.no, function(i){
  dat <- data[data$Peak.no==i, c("Au_Signal")]
})
# Save data for later use
# save(PeakList, file = "20200429_RPMI_FBS_0_1_01_PeakList.Rdata")

# Connvert each Peak (vector) to data frame
PeakList <- lapply(PeakList, as.data.frame.list)

# Put column names for each Peak with the form: DT1, DT2, ..., DTn
PeakList <- lapply(PeakList, function(i){
  names(i) <- paste("DT",1:length(i),sep="")
  return(i)
})

# Combine all Peaks, row by row (i.e. Peak 1 is row 1, Peak 2 is row 2...) to make data table "Table.3"
Table.3 <- rbindlist(PeakList, fill = TRUE)
Table.3[is.na(Table.3)] <- 0
DT_start <- data.frame("DT_start" <- 0); colnames(DT_start) <- "DT_start"
DT_end <- data.frame("DT_end" <- 0); colnames(DT_end) <- "DT_end"
Table.3 <- cbind(DT_start, Table.3, DT_end)
# Save Table.3 for later use
# save(Table.3, file = "20200429_RPMI_FBS_0_1_01_Table3.Rdata")


## Clustering peaks using Dynamic Time Wrapping method 
#=====================================================================================================================
peaksclusters <- tsclust(PeakList[1:500], type = "partitional", k = 5L, 
                         distance = "dtw_basic", centroid = "pam", 
                         seed = 3247L, trace = TRUE,
                         args = tsclust_args(dist = list(window.size = 5L)))
plot(peaksclusters, type = "centroids")
plot(peaksclusters, type = "centroids", clus=3L)
plot(peaksclusters, type = "sc")
plot(peaksclusters, type = "sc", clus=2L)
