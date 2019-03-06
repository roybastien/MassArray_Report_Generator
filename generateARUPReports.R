require("ggplot2")
require("ggpubr")
require("RSQLite")
require("gridExtra")

#Set current working directory
setwd("//labshare/instruments/Molecular Oncology 656/Sequenom/MGMT_MLH1")

#Include different function files to make this main file easier to follow and less redundant.
source("functions.R")

# directory = choose.dir(default = "I:\\MGMT_MLH1", caption = "Please choose the folder of the run you would like to analyze.")
directory = choose.dir(default = "//labshare/instruments/Molecular Oncology 656/Sequenom/MGMT_MLH1", caption = "Please choose the folder of the run you \n would like to analyze.")
setwd(directory)
datafile = list.files(path = directory, pattern = "\\.csv$")
dbname = list.files(path = directory, pattern = "\\.db$")

#Create a directory to place the reports
plateID = substr(dbname, 1, nchar(dbname)-12)
reportsFolder = paste(plateID, "Reports")
dir.create(file.path(directory, reportsFolder))


#Connect to SQLite database containing sample names, etc
con = dbConnect(SQLite(), dbname = dbname)

#Read peak height data
myData = read.csv(datafile, header = TRUE)

#Create a progress bar for the user
numMGMT = nrow(queryDBwithString(con, "SELECT WELL_POSITION from SPECTRAQUALITY WHERE AMPLICON_ID LIKE 'MGMT'"))*2
numMLH1 = nrow(queryDBwithString(con, "SELECT WELL_POSITION from SPECTRAQUALITY WHERE AMPLICON_ID LIKE 'MLH1'"))
total = numMGMT + numMLH1
progress = 0
pb <- winProgressBar(title = "progress bar", min = 0, max = total, width = 300)

#This is the maximum width which we fit around the CpG site for each data point
plusMinus = 5

#Increment size of the gaussian curve calculation
step = 0.25

#Check to make sure there are MGMT/SAR samples to run
if (numMGMT > 0){
  
  #Get the control sample Pass/Fail data, done only once for all samples
  controlQC = getMGMTControlPassFail()
  # print(controlQC)
  
  #Ensure all control samples pass else the sample fails
  controlQCpass = if (controlQC[1] == "Fail" || controlQC[2] == "Fail" || controlQC[3] == "Fail" || controlQC[4] == "Fail") "Fail" else "Pass"
  
  # Build the rule table to go on the bottom of the report
  controlrow1 = c("Control", "Percent Methylation", "Pass/Fail")
  controlrow2 = c("Non-methylated Control", "0-9%", controlQC[1])
  controlrow3 = c("20% Methylated Control", "10-40%", controlQC[2])
  controlrow4 = c("100% Methylated Control", "90-100%", controlQC[3])
  controlrow5 = c("NTC", "np mean <10", controlQC[4])
  ggcontroltable = ggtexttable(rbind(controlrow1, controlrow2, controlrow3, controlrow4, controlrow5), rows = NULL)
  
  #Get a list of all unique sample names in the data file
  samples = queryDBwithString(con, "SELECT WELL_POSITION from SPECTRAQUALITY WHERE AMPLICON_ID LIKE 'MGMT'")
  
  #Just generating a single report now during develpment, uncommment 
  #to generare reports for all samples in a data file
  #Uncomment when running on all samples
  for (wellID in samples$WELL_POSITION){
    
    #Comment out next line when running on all samples
    # wellID = 'E08'
    
    #The following block is used to pull accession numbers and any other necessary data from the 
    #instrument database.
    sampleID = toString(queryDBwithString(con, paste("SELECT SAMPLE_ID from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
    
    #Gather all the assay specific information
    npMean <- toString(queryDBwithString(con, paste("SELECT np_mean from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
    numSpectra <- toString(queryDBwithString(con, paste("SELECT Num_Spectra from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
    methylationValue <- queryDBwithString(con, paste("SELECT VALUE from METHYLATION WHERE SAMPLE_ID LIKE '", sampleID, "' AND AMPLICON LIKE 'MGMT'", sep = ""))
    warningFlag <- toString(queryDBwithString(con, paste("SELECT Status from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
    
    meth1 = paste(toString(methylationValue$VALUE[1]*100), "%", sep = "")
    meth2 = paste(toString(methylationValue$VALUE[2]*100), "%", sep = "")
    meth3 = paste(toString(methylationValue$VALUE[3]*100), "%", sep = "")
    meth4 = paste(toString(methylationValue$VALUE[4]*100), "%", sep = "")
    meth5 = paste(toString(methylationValue$VALUE[5]*100), "%", sep = "")
    
    #Subset data for a single sample
    sampleData = myData[ which(myData$Well==wellID), ]
    
    #Subset data based on Mass/Charge for individual regions requested
    data1 = sampleData[sampleData[3]>"4145" & sampleData[3]<"4200", ]
    data2 = sampleData[sampleData[3]>"4765" & sampleData[3]<"4835", ]
    data3 = sampleData[sampleData[3]>"2545" & sampleData[3]<"2590", ]
    data4 = sampleData[sampleData[3]>"5300" & sampleData[3]<"5375", ]
    data5 = sampleData[sampleData[3]>"3160" & sampleData[3]<"3210", ]
    
    intensity1 = toString(round(colSums(data1[4]), digits = 0))
    intensity2 = toString(round(colSums(data2[4]), digits = 0))
    intensity3 = toString(round(colSums(data3[4]), digits = 0))
    intensity4 = toString(round(colSums(data4[4]), digits = 0))
    intensity5 = toString(round(colSums(data5[4]), digits = 0))
    
    npMeanPass = if(as.integer(npMean) >= 10) "Pass" else "Fail"
    numSpectraPass = if(numSpectra >= 3) "Pass" else "Fail"
    intensityPass = if(strtoi(intensity1) >= 600 & strtoi(intensity2) >= 600 & strtoi(intensity3) >= 600 & strtoi(intensity4) >= 600 & strtoi(intensity5) >= 600) "Pass" else "Fail"
    warningPass = if(warningFlag == "ok") "Pass" else "Fail"
    methValuePass = if(methylationValue$VALUE[1]>0 && methylationValue$VALUE[2]>0 && methylationValue$VALUE[3]>0 && methylationValue$VALUE[4]>0 && methylationValue$VALUE[5]>0) "Pass" else "Review"
    totalMeth = as.integer(calcTotalMgmtMethylation(sampleID))
    methResultPass = if(totalMeth>=10) {"Pass"} else "Fail"
    
    #Detect if all MGMT metrics pass
    MGMTPass = if (npMeanPass == "Pass" & numSpectraPass == "Pass" & warningPass == "Pass" & intensityPass == "Pass") "Pass" else "Fail"
    
    #Fill the rows with identifing text and the data from above.
    row1 = c("","","","","","","Pass/Fail","QC Check")
    row2 = c("Np mean","","",npMean,"","",npMeanPass,">=10")
    row3 = c("Number of spectra","","",numSpectra,"","",numSpectraPass,">=3/5")
    row4 = c("Warning flags","","",warningFlag,"","",warningPass,"ok")
    row5 = c("CpG unit # (CpGs/unit)","1 (2)","2 (3)","3 (2)","4 (3)","5 (2)","","")
    row6 = c("Sum of peak(s) intensity per unit", intensity1, intensity2, intensity3, intensity4, intensity5, intensityPass,">=600 for each unit")
    row7 = c("Methylation value per unit (%)", meth1, meth2, meth3, meth4, meth5, methValuePass,"Calculated")
    row8 = c("Total methylation from MGMT and SAR","","",paste(totalMeth, "%", sep = ""),"","","Pass","Calculated")
    
    #Added by Roy 7/30
    if (substr(sampleID, 1, 4) == "NTC"){
      npMeanPass = if(as.integer(npMean)<10) "Pass" else "Fail"
      intensityPass = if(strtoi(intensity1) < 400 & strtoi(intensity2) < 400 & strtoi(intensity3) < 400 & strtoi(intensity4) < 400 & strtoi(intensity5) < 400) "Pass" else "Fail"
      if (npMeanPass == "Pass" & intensityPass == "Fail"){
        intensityPass = "Review"
      }
      row1 = c("","","","","","","Pass/Fail","QC Check")
      row2 = c("Np mean","","",npMean,"","",npMeanPass,"<10")
      row3 = c("Number of spectra","","",numSpectra,"","","Pass","Any Spectra")
      row4 = c("CpG unit # (CpGs/unit)","1 (2)","2 (3)","3 (2)","4 (3)","5 (2)","","")
      row5 = c("Sum of peak(s) intensity per unit", intensity1, intensity2, intensity3, intensity4, intensity5, intensityPass,"<400 for each unit")
      MGMTggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5), rows = NULL)
    } else{
      #Build the footer table with all the assay info
      MGMTggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5, row6, row7, row8), rows = NULL) 
    }
    
    #Build the footer table with all the assay info
    # MGMTggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5, row6, row7, row8), rows = NULL) 
    
    #Combine all data to feed to gauss
    allSelected = rbind(data1, data2, data3, data4, data5)
    
    #Make the single peak a nice gaussian curve
    gaussData <- makePeaksGaussian(allSelected, plusMinus, step)
    
    #Insert data points here so the lines always go to the ends of the plot
    border1a = data.frame(Mass.Charge=4145, Height=0)
    border1b = data.frame(Mass.Charge=4200, Height=0)
    border2a = data.frame(Mass.Charge=4765, Height=0)
    border2b = data.frame(Mass.Charge=4835, Height=0)
    border3a = data.frame(Mass.Charge=2545, Height=0)
    border3b = data.frame(Mass.Charge=2590, Height=0)
    border4a = data.frame(Mass.Charge=5300, Height=0)
    border4b = data.frame(Mass.Charge=5375, Height=0)
    border5a = data.frame(Mass.Charge=3160, Height=0)
    border5b = data.frame(Mass.Charge=3210, Height=0)
    gaussData = rbind(gaussData, border1a, border1b, border2a, border2b, border3a, border3b, border4a, border4b, border5a, border5b)
    
    #Sort data frame for plotting
    gaussData = gaussData[order(gaussData$Mass.Charge), ]
    
    #Partition the data to our areas of interest
    gaussData1 = gaussData[gaussData[1]>="4145" & gaussData[1]<="4200", ]
    gaussData2 = gaussData[gaussData[1]>="4765" & gaussData[1]<="4835", ]
    gaussData3 = gaussData[gaussData[1]>="2545" & gaussData[1]<="2590", ]
    gaussData4 = gaussData[gaussData[1]>="5300" & gaussData[1]<="5375", ]
    gaussData5 = gaussData[gaussData[1]>="3160" & gaussData[1]<="3210", ]
    
    #Sanity check, comment next line during production
    # plot(gaussData$Mass.Charge, gaussData$Height, type = "l", xlab = "Mass/Charge", ylab = "Intensity")
    
    #Build individual plots for specific regions.  Vertical lines indicate CpG sites
    #specific to that region.  Horizontal line indicates background.
    a = ggplot(data=gaussData1, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 4155, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 4171, color = "blue", linetype = "longdash") + 
      geom_vline(xintercept = 4187, color = "blue", linetype = "longdash") +
      ggtitle(paste("MGMT CpG Sites 72 & 73\n", meth1, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    b = ggplot(data=gaussData2, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 4774, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 4790, color = "blue", linetype = "longdash") + 
      geom_vline(xintercept = 4806, color = "blue", linetype = "longdash") + 
      geom_vline(xintercept = 4822, color = "blue", linetype = "longdash") +
      ggtitle(paste("MGMT CpG Sites  74, 75, & 76\n", meth2, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    c = ggplot(data=gaussData3, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 2549, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 2565, color = "blue", linetype = "longdash") + 
      geom_vline(xintercept = 2581, color = "blue", linetype = "longdash") +
      ggtitle(paste("MGMT CpG Sites  77 & 78\n", meth3, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    d = ggplot(data=gaussData4, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") + 
      geom_vline(xintercept = 5312, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 5328, color = "blue", linetype = "longdash") + 
      geom_vline(xintercept = 5344, color = "blue", linetype = "longdash") +
      geom_vline(xintercept = 5360, color = "blue", linetype = "longdash") +
      ggtitle(paste("MGMT CpG Sites  79, 80, & 81\n", meth4, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    e = ggplot(data=gaussData5, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") + 
      geom_vline(xintercept = 3168, color = "red", linetype = "longdash") +
      geom_vline(xintercept = 3184, color = "blue", linetype = "longdash") +
      geom_vline(xintercept = 3200, color = "blue", linetype = "longdash") +
      ggtitle(paste("MGMT CpG Sites  82 & 83\n", meth5, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    
    #Place all plots into a single object
    MGMTplots <- ggarrange(a,b,c,d,e, ncol = 3, nrow = 2)
    
    #Update the progress bar
    progress = progress + 1
    setWinProgressBar(pb, progress, title=paste(round(progress/total*100, 0), "% Complete", sep = ""))
    
    ##############################################################################################
    # Begin SAR reporting here.
    wellID <- toString(queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE '", sampleID, "' AND AMPLICON_ID LIKE 'SAR'", sep = "")))
    
    #The following block will is to pull accession numbers and any other necessary data from the 
    #instrument database.
    sampleID = toString(queryDBwithString(con, paste("SELECT SAMPLE_ID from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
    
    #Gather all the assay specific information
    npMean <- toString(queryDBwithString(con, paste("SELECT np_mean from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
    numSpectra <- toString(queryDBwithString(con, paste("SELECT Num_Spectra from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
    warningFlag <- toString(queryDBwithString(con, paste("SELECT Status from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
    methylationValue <- queryDBwithString(con, paste("SELECT VALUE from METHYLATION WHERE SAMPLE_ID LIKE '", sampleID, "' AND AMPLICON LIKE 'SAR'", sep = ""))
    
    meth1 = paste(toString(methylationValue$VALUE[2]*100), "%", sep = "")
    meth2 = paste(toString(methylationValue$VALUE[3]*100), "%", sep = "")
    meth3 = paste(toString(methylationValue$VALUE[4]*100), "%", sep = "")
    
    #Subset data for a single sample
    sampleData = myData[ which(myData$Well==wellID), ]
    
    data1 = sampleData[sampleData[3]>"2865" & sampleData[3]<"2920", ]
    data2 = sampleData[sampleData[3]>"2540" & sampleData[3]<"2575", ]
    data3 = sampleData[sampleData[3]>"3490" & sampleData[3]<"3525", ]
    
    intensity1 = toString(round(colSums(data1[4]), digits = 0))
    intensity2 = toString(round(colSums(data2[4]), digits = 0))
    intensity3 = toString(round(colSums(data3[4]), digits = 0))
    
    npMeanPass = if(as.integer(npMean) >= 10) "Pass" else "Fail"
    numSpectraPass = if(numSpectra >= 3) "Pass" else "Fail"
    intensityPass = if(strtoi(intensity1) >= 600 & strtoi(intensity2) >= 600 & strtoi(intensity3) >= 600) "Pass" else "Fail"
    warningPass = if(warningFlag == "ok") "Pass" else "Fail"
    methValuePass = if(methylationValue$VALUE[2]>0 && methylationValue$VALUE[3]>0 && methylationValue$VALUE[4]>0) "Pass" else "Review"
    
    #Determine if SAR QC metrics pass
    SARpass = if (npMeanPass == "Pass" & numSpectraPass == "Pass" & warningPass == "Pass" & intensityPass == "Pass") "Pass" else "Fail"
    
    #Fill the rows with identifing text and the data from above.
    row1 = c("","","","","Pass/Fail","QC Check")
    row2 = c("Np mean","",npMean,"",npMeanPass,">=10")
    row3 = c("Number of spectra","",numSpectra,"",numSpectraPass,">=3/5")
    row4 = c("Warning flags","",warningFlag,"",warningPass,"ok")
    row5 = c("CpG unit # (CpGs/unit)","1 (2)","2 (1)","3 (1)","","")
    row6 = c("Sum of peak(s) intensity per unit", intensity1, intensity2, intensity3, intensityPass,">=600 for each unit")
    row7 = c("Methylation value per unit (%)", meth1, meth2, meth3, methValuePass,"Calculated")
    row8 = c("Total methylation from MGMT and SAR","",paste(totalMeth, "%", sep = ""),"","Pass","Calculated")
    
    #Added by Roy 7/30
    if (substr(sampleID, 1, 4) == "NTC"){
      npMeanPass = if(as.integer(npMean)<10) "Pass" else "Fail"
      intensityPass = if(strtoi(intensity1) < 400 & strtoi(intensity2) < 400 & strtoi(intensity3) < 400) "Pass" else "Fail"
      if (npMeanPass == "Pass" & intensityPass == "Fail"){
        intensityPass = "Review"
      }
      row1 = c("","","","","Pass/Fail","QC Check")
      row2 = c("Np mean","",npMean,"",npMeanPass,"<10")
      row3 = c("Number of spectra","",numSpectra,"","Pass","Any Spectra")
      row4 = c("CpG unit # (CpGs/unit)","1 (2)","2 (1)","3 (1)","","")
      row5 = c("Sum of peak(s) intensity per unit", intensity1, intensity2, intensity3, intensityPass,"<400 for each unit")
      SARggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5), rows = NULL)
    } else{
      #Build the footer table with all the assay info
      SARggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5, row6, row7, row8), rows = NULL) 
    }
    
    #Build the footer table with all the assay info
    # SARggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5, row6, row7, row8), rows = NULL)
    
    #Combine all data to feed to gauss
    allSelected = rbind(data1, data2, data3)
    
    #Make the single peak a nice gaussian curve
    gaussData <- makePeaksGaussian(allSelected, plusMinus, step)
    
    #Insert data points here so the lines always go to the ends of the plot
    border1a = data.frame(Mass.Charge=2865, Height=0)
    border1b = data.frame(Mass.Charge=2920, Height=0)
    border2a = data.frame(Mass.Charge=2540, Height=0)
    border2b = data.frame(Mass.Charge=2575, Height=0)
    border3a = data.frame(Mass.Charge=3490, Height=0)
    border3b = data.frame(Mass.Charge=3525, Height=0)
    gaussData = rbind(gaussData, border1a, border1b, border2a, border2b, border3a, border3b)
    
    #Sort data frame for plotting
    gaussData = gaussData[order(gaussData$Mass.Charge), ]
    
    #Partition the data to our areas of interest
    gaussData1 = gaussData[gaussData[1]>="2865" & gaussData[1]<="2920", ]
    gaussData2 = gaussData[gaussData[1]>="2540" & gaussData[1]<="2575", ]
    gaussData3 = gaussData[gaussData[1]>="3490" & gaussData[1]<="3525", ]
    
    
    #Sanity check, comment next line during production
    # plot(gaussData$Mass.Charge, gaussData$Height, type = "l", xlab = "Mass/Charge", ylab = "Intensity")
    
    #Build individual plots for specific regions.  Vertical lines indicate CpG sites
    #specific to that region.  Horizontal line indicates background.
    a = ggplot(data=gaussData1, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 2879, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 2895, color = "blue", linetype = "longdash") + 
      geom_vline(xintercept = 2911, color = "blue", linetype = "longdash") +
      ggtitle(paste("SAR CpG Sites 88 & 89\n", meth1, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    b = ggplot(data=gaussData2, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 2550, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 2566, color = "blue", linetype = "longdash") + 
      ggtitle(paste("SAR CpG Site 87\n", meth2, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    c = ggplot(data=gaussData3, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 3497, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 3513, color = "blue", linetype = "longdash") +     
      ggtitle(paste("SAR CpG Site 86\n", meth3, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    
    #Place all plots into a single object
    SARplots <- ggarrange(a,b,c, ncol = 3, nrow = 1)
    
    # Build the rule table to go on the bottom of the report
    row1 = c("If...", "Then...")
    row2 = c("Total % methylation is 0 - 9%", "The sample is reported as NOT DETECTED for MGMT promoter methlyation")
    row3 = c("Total % methylation is 10 - 29%", "The sample is reported as LOW LEVEL for MGMT promoter methlyation")
    row4 = c("Total % methylation is 30 - 100%", "The sample is reported as DETECTED for MGMT promoter methlyation")
    ggruletable = ggtexttable(rbind(row1, row2, row3, row4), rows = NULL)
    
    #Determine if both MGMT and SAR pass QC
    # runPass = if (MGMTPass == "Pass" & SARpass == "Pass" & controlQCpass == "Pass") "Pass" else "Fail"
    # changed above line 9/13/18
    runPass = if (MGMTPass == "Pass" & SARpass == "Pass") "Pass" else "Fail"
    
    accession = paste("Accession Number: ", sampleID)
    assay = "Assay: MGMT/SAR Methylation"
    reportDate = paste("Report Date:", Sys.Date(), sep = " ")
    runName = paste("Run Name: ", plateID, sep = "")
    totalMeth = as.integer(calcTotalMgmtMethylation(sampleID))
    percMeth = paste("Total % Methylation: ", totalMeth, "%", sep = "")
    methResultPass = if(totalMeth>=30) {"Detected"} else if(totalMeth>=10 && totalMeth<30) {"Low-level Methylation"} else {"Not Detected"}
   
     #Added by Roy 9/12
    if (substr(sampleID, 1, 3) == "NTC"){
      methResultPass = controlQC[4]
    } else if (substr(sampleID, 1, 3) == "NEG"){
      methResultPass = controlQC[1]
    } else if (substr(sampleID, 1, 6) == "POS_20"){
      methResultPass = controlQC[2]
    } else if (substr(sampleID, 1, 7) == "POS_100"){
      methResultPass = controlQC[3]
    } else{
      methResultPass = if(runPass == "Fail") "Fail" else methResultPass
    }
    
    result = paste("Result:", methResultPass)
    
    #Build head table of all the patient information
    ggheadtable <- ggtexttable(cbind(c(accession, assay), c(reportDate, result), c(runName, percMeth)))
    
    #Build final MGMT report
    MGMTfinal <- ggarrange(ggheadtable, MGMTplots, MGMTggtailtable, ncol = 1, nrow = 3, heights = c(0.1, 0.6, 0.3))
    # ggsave(paste(paste(reportsFolder, "/", sep = ""), sampleID, "MGMT report.pdf"), MGMTfinal, width = 11, height = 8.5, dpi = 300)
  
    #Build final report
    SARfinal <- ggarrange(ggheadtable, SARplots, SARggtailtable, ggruletable, ggcontroltable, ncol = 1, nrow = 5, heights = c(0.08, 0.35, 0.25, 0.15, 0.17))
    # ggsave(paste(paste(reportsFolder, "/", sep = ""), sampleID, "SAR report.pdf"), SARfinal, width = 11, height = 8.5, dpi = 300)
    
    #Save reports in multipage pdf
    plotList = list(MGMTfinal, SARfinal)
    plotGrob <- do.call(marrangeGrob, list(grobs=plotList, nrow = 1, ncol = 1, top = NULL))
    ggsave(paste(reportsFolder, "/", sampleID, " MGMT Report.pdf", sep = ""), plotGrob, width = 11, height = 8.5, dpi = 300)
  
    #Update the progress bar
    progress = progress + 1
    setWinProgressBar(pb, progress, title=paste(round(progress/total*100, 0), "% Complete", sep = ""))
    
  #Uncomment to run on all MGMT/SAR samples in the data file
  }
}
##############################################################################################

##############################################################################################
# Begin MLH1 reporting here.

#Check to make sure there are MLH1 samples to run
if (numMLH1 > 0){
  
  #Get the control sample Pass/Fail data, done only once for all samples
  controlQC = getMLH1ControlPassFail(myData)
  # print(controlQC)
  
  #Ensure all control samples pass else the sample fails
  controlQCpass = if (controlQC[1] == "Fail" || controlQC[2] == "Fail" || controlQC[3] == "Fail" || controlQC[4] == "Fail") "Fail" else "Pass"
  
  # Build the rule table to go on the bottom of the report
  controlrow1 = c("Control", "Percent Methylation", "Pass/Fail")
  controlrow2 = c("Non-methylated Control", "0-9%", controlQC[1])
  controlrow3 = c("20% Methylated Control", "10-40%", controlQC[2])
  controlrow4 = c("100% Methylated Control", "90-100%", controlQC[3])
  controlrow5 = c("NTC", "np mean <10", controlQC[4])
  ggcontroltable = ggtexttable(rbind(controlrow1, controlrow2, controlrow3, controlrow4, controlrow5), rows = NULL)
  
  #Get a list of all unique sample names in the data file
  samples = queryDBwithString(con, "SELECT WELL_POSITION from SPECTRAQUALITY WHERE AMPLICON_ID LIKE 'MLH1'")
  
  #Just generating a single report now during develpment, uncommment 
  #to generare reports for all samples in a data file
  #Uncomment when running on all samples
  for (wellID in samples$WELL_POSITION){
    
    #Comment out next line when running on all samples
    # wellID = 'A12'
    
    #The following block will be used to pull accession numbers and any other necessary data from the 
    #instrument database.
    sampleID = toString(queryDBwithString(con, paste("SELECT SAMPLE_ID from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
    
    #Gather all the assay specific information
    npMean <- toString(queryDBwithString(con, paste("SELECT np_mean from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
    numSpectra <- toString(queryDBwithString(con, paste("SELECT Num_Spectra from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
    warningFlag <- toString(queryDBwithString(con, paste("SELECT Status from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
    methylationValue <- queryDBwithString(con, paste("SELECT VALUE from METHYLATION WHERE SAMPLE_ID LIKE '", sampleID, "' AND AMPLICON LIKE 'MLH1'", sep = ""))
    
    meth1 = paste(toString(methylationValue$VALUE[2]*100), "%", sep = "")
    meth2 = paste(toString(methylationValue$VALUE[6]*100), "%", sep = "")
    meth3 = paste(toString(methylationValue$VALUE[8]*100), "%", sep = "")
    meth4 = paste(toString(methylationValue$VALUE[9]*100), "%", sep = "")
    
    #Subset data for a single sample
    sampleData = myData[ which(myData$Well==wellID), ]
    
    # print(sarWell)
    data1 = sampleData[sampleData[3]>"2580" & sampleData[3]<"2620", ]
    data2 = sampleData[sampleData[3]>"1880" & sampleData[3]<"1930", ]
    data3 = sampleData[sampleData[3]>"2135" & sampleData[3]<"2165", ]
    data4 = sampleData[sampleData[3]>"2870" & sampleData[3]<"2905", ]
    
    intensity1 = toString(round(colSums(data1[4]), digits = 0))
    intensity2 = toString(round(colSums(data2[4]), digits = 0))
    intensity3 = toString(round(colSums(data3[4]), digits = 0))
    intensity4 = toString(round(colSums(data4[4]), digits = 0))
    
    npMeanPass = if(as.integer(npMean) >= 10) "Pass" else "Fail"
    numSpectraPass = if(numSpectra >= 3) "Pass" else "Fail"
    intensityPass = if(strtoi(intensity1) >= 600 & strtoi(intensity2) >= 600 & strtoi(intensity3) >= 600 & strtoi(intensity4) >= 600) "Pass" else "Fail"
    warningPass = if(warningFlag == "ok") "Pass" else "Fail"
    methValuePass = if(methylationValue$VALUE[2]>0 && methylationValue$VALUE[6]>0 && methylationValue$VALUE[8]>0 && methylationValue$VALUE[9]>0) "Pass" else "Review"
    
    #Detect if all MLH1 metrics pass
    # MLH1Pass = if (npMeanPass == "Pass" & numSpectraPass == "Pass" & warningPass == "Pass" & intensityPass == "Pass" & controlQCpass == "Pass") "Pass" else "Fail"
    #changed above line 9/13/18
    MLH1Pass = if (npMeanPass == "Pass" & numSpectraPass == "Pass" & warningPass == "Pass" & intensityPass == "Pass") "Pass" else "Fail"

    #Gather all information of the header table
    accession = paste("Accession Number: ", sampleID)
    assay = "Assay: MLH1 Methylation"
    reportDate = paste("Report Date:", Sys.Date(), sep = " ")
    runName = paste("Run Name: ", plateID, sep = "")
    mlh1TotalMeth = as.integer(calcTotalMlh1Methylation(sampleID))
    percMeth = paste("Total % Methylation: ", mlh1TotalMeth, "%", sep = "")
    methResultPass = if(mlh1TotalMeth>=10) "Detected" else "Not Detected"
    if (MLH1Pass == "Fail") "Fail" else methResultPass
    result = paste("Result:", methResultPass)

    # #Build head table of all the patient information
    # ggheadtable <- ggtexttable(cbind(c(accession, assay), c(reportDate, result), c(runName, percMeth)))
    
    #Fill the rows with identifing text and the data from above.
    row1 = c("","","","","","Pass/Fail","QC Check")
    row2 = c("Np mean","","",npMean,"",npMeanPass,">=10")
    row3 = c("Number of spectra","","",numSpectra,"",numSpectraPass,">=3/5")
    row4 = c("Warning flags","","",warningFlag,"",warningPass,"ok")
    row5 = c("CpG unit # (CpGs/unit)","2 (1)","6 (2)","8 (1)","9 (1)","","")
    row6 = c("Sum of peak(s) intensity per unit", intensity1, intensity2, intensity3, intensity4, intensityPass,">=600 for each unit")
    row7 = c("Methylation value per unit (%)", meth1, meth2, meth3, meth4, methValuePass,"Calculated")
    row8 = c("Total methylation from MLH1","","",paste(mlh1TotalMeth, "%", sep = ""),"","Pass","Calculated")
    
    #Added by Roy 7/30
    if (substr(sampleID, 1, 4) == "NTC"){
      intensityPass = if(strtoi(intensity1) < 400 & strtoi(intensity2) < 400 & strtoi(intensity3) < 400 & strtoi(intensity4) < 400) "Pass" else "Fail"
      npMeanPass = if(as.integer(npMean)<10) "Pass" else "Fail"
      if (npMeanPass == "Pass" & intensityPass == "Fail"){
        intensityPass = "Review"
      }
      row1 = c("","","","","","Pass/Fail","QC Check")
      row2 = c("Np mean","","",npMean,"",npMeanPass,"<10")
      row3 = c("Number of spectra","","",numSpectra,"","Pass","Any Spectra")
      row4 = c("CpG unit # (CpGs/unit)","2 (1)","6 (2)","8 (1)","9 (1)","","")
      row5 = c("Sum of peak(s) intensity per unit", intensity1, intensity2, intensity3, intensity4, intensityPass,"<400 for each unit")
      ggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5), rows = NULL)
    } else if (substr(sampleID, 1, 4) == "NEG"){
      methResultPass = controlQC[1]
      ggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5, row6, row7, row8), rows = NULL)
    } else if (substr(sampleID, 1, 7) == "POS_20"){
      methResultPass = controlQC[2]
      ggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5, row6, row7, row8), rows = NULL)
    } else if (substr(sampleID, 1, 8) == "POS_100"){
      methResultPass = controlQC[3]
      ggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5, row6, row7, row8), rows = NULL)
    } else{
      #Build the footer table with all the assay info
      ggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5, row6, row7, row8), rows = NULL) 
    }
    
    #Build head table of all the patient information
    ggheadtable <- ggtexttable(cbind(c(accession, assay), c(reportDate, result), c(runName, percMeth)))
    
    #Build the footer table with all the assay info
    # ggtailtable <- ggtexttable(rbind(row1, row2, row3, row4, row5, row6, row7, row8), rows = NULL)
    
    #Combine all data to feed to gauss
    allSelected = rbind(data1, data2, data3, data4)
    
    #Make the single peak a nice gaussian curve
    gaussData <- makePeaksGaussian(allSelected, plusMinus, step)
    
    #Insert data points here so the lines always go to the ends of the plot
    border1a = data.frame(Mass.Charge=2580, Height=0)
    border1b = data.frame(Mass.Charge=2620, Height=0)
    border2a = data.frame(Mass.Charge=1880, Height=0)
    border2b = data.frame(Mass.Charge=1930, Height=0)
    border3a = data.frame(Mass.Charge=2135, Height=0)
    border3b = data.frame(Mass.Charge=2165, Height=0)
    border4a = data.frame(Mass.Charge=2870, Height=0)
    border4b = data.frame(Mass.Charge=2905, Height=0)
    gaussData = rbind(gaussData, border1a, border1b, border2a, border2b, border3a, border3b, border4a, border4b)
    
    #Sort data frame for plotting
    gaussData = gaussData[order(gaussData$Mass.Charge), ]
    
    #Partition the data to our areas of interest
    gaussData1 = gaussData[gaussData[1]>="2580" & gaussData[1]<="2620", ]
    gaussData2 = gaussData[gaussData[1]>="1880" & gaussData[1]<="1930", ]
    gaussData3 = gaussData[gaussData[1]>="2135" & gaussData[1]<="2165", ]
    gaussData4 = gaussData[gaussData[1]>="2870" & gaussData[1]<="2905", ]
    
    
    #Sanity check, comment next line during production
    # plot(gaussData$Mass.Charge, gaussData$Height, type = "l", xlab = "Mass/Charge", ylab = "Intensity")
    
    #Build individual plots for specific regions.  Vertical lines indicate CpG sites
    #specific to that region.  Horizontal line indicates background.
    a = ggplot(data=gaussData1, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 2590, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 2606, color = "blue", linetype = "longdash") + 
      ggtitle(paste("CpGUnit_2\n", meth1, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    b = ggplot(data=gaussData2, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 1891, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 1907, color = "blue", linetype = "longdash") + 
      geom_vline(xintercept = 1923, color = "blue", linetype = "longdash") + 
      ggtitle(paste("CpGUnit_6\n", meth2, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    c = ggplot(data=gaussData3, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 2140, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 2156, color = "blue", linetype = "longdash") + 
      ggtitle(paste("CpGUnit_8\n", meth3, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    d = ggplot(data=gaussData4, aes(x=Mass.Charge, y=Height)) + geom_line() + 
      geom_hline(yintercept = 200, color = "brown4") +
      geom_vline(xintercept = 2879, color = "red", linetype = "longdash") + 
      geom_vline(xintercept = 2895, color = "blue", linetype = "longdash") + 
      ggtitle(paste("CpGUnit_9\n", meth4, " Methylated")) + 
      theme(plot.title = element_text(hjust=0.5)) + labs(x = "M/Z", y = "Intensity")
    
    #Place all plots into a single object
    plots <- ggarrange(a,b,c,d, ncol = 4, nrow = 1)
    
    # Build the rule table to go on the bottom of the report
    row1 = c("If...", "Then...")
    row2 = c("Total % methylation is 0 - 9%", "The sample is reported as NOT DETECTED for MLH1 promoter methlyation")
    row3 = c("Total % methylation is 10 - 100%", "The sample is reported as DETECTED for MLH1 promoter methlyation")
    ggruletable = ggtexttable(rbind(row1, row2, row3), rows = NULL)
    
    #Build and save final report
    final <- ggarrange(ggheadtable, plots, ggtailtable, ggruletable, ggcontroltable, ncol = 1, nrow = 5, heights = c(0.1, 0.35, 0.25, 0.12, 0.18))
    ggsave(paste(paste(reportsFolder, "/", sep = ""), sampleID, "MLH1 report.pdf"), final, width = 11, height = 8.5, dpi = 300)
  
    #Update the progress bar
    progress = progress + 1
    setWinProgressBar(pb, progress, title=paste( round(progress/total*100, 0), "% Complete"))
    
  #Uncomment to run on all MLH1 samples in the data file
  }
}

#Disconnect from the database
dbDisconnect(con)

#Close the progress bar
close(pb)


