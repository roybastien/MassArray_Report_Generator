require("RSQLite")

#Gets a list of command line arguments provided to the program,
#if set to FALSE it lists everything, if set to TRUE it only lists the
#arguments provided by the user
args=(commandArgs(FALSE))
if(length(args)==0){
  # print("No Arguments Provided")
}else{
  #Extract composite plate ID('cp'), database file name('db'), and R Script location.
  #The location is chosen by the EpiTyper GUI.
  #Dummy data provided below
  # print(args)
  plateID <- substring(args[7], 5, nchar(args[7]))
  dbname <- substring(args[8], 5, nchar(args[8]))
  # print(dbname)
}

mainDir = "I:\\MGMT_MLH1"
dir.create(file.path(mainDir, plateID))

#This is where the copied file will go
# targetFolder = paste("C:/Users/GLA656agenadb/Desktop/RoyTemp/", plateID, sep = "")
targetFolder = paste("I:\\MGMT_MLH1/", plateID, sep = "")

#Copy database file
file.copy(dbname, targetFolder)

#Move into the target folder and name the db something useful
setwd(targetFolder)
file.rename(dbname, paste(plateID, "_database.db", sep = ""))
db = substr(dbname, nchar(dbname)-10, nchar(dbname))
file.remove(paste(targetFolder, "//",db, sep = ""))

#Get the last export folder sorted by date/time
spectraFolder <- tail(sort(list.files("C:/MassARRAY/EpiTYPER 1.3/bin/EpiTYPERReports/EpiTYPERSpectraExport/")), n=1)

#Move into the reports folder
setwd(paste("C:/MassARRAY/EpiTYPER 1.3/bin/EpiTYPERReports/EpiTYPERSpectraExport/", spectraFolder, sep = ""))
csvFile = list.files(path = ".", pattern = "\\.csv$")

file.copy(csvFile, targetFolder)