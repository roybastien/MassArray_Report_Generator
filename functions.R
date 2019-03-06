

##############################################################################################
#Runs given query 
queryDBwithString = function(con, query){
  myQuery = dbSendQuery(con, query)
  myQueryData = dbFetch(myQuery, n = -1)
  dbClearResult(myQuery)
  value = myQueryData
}
##############################################################################################


##############################################################################################

makePeaksGaussian = function(allSelected, plusMinus, step){
  
  #Empty data frame to hold the Gaussian fitted data
  gaussData <- data.frame(mass=double(), height=double())
  
  #Build the Gaussian fitted data frame from each CpG site
  # for (row in sampleData){
  for(i in 1:nrow(allSelected)) {
    row <- allSelected[i,]
    cpgMass = row$Mass.Charge
    cpgHeight = row$Height
    cpgWidth = row$Width
    
    #This is the maximum width which we fit around the CpG site
    # plusMinus = 5
    
    #Where the Gaussian calculations are started relative to the Mass.Charge
    min = cpgMass - plusMinus
    
    #Increment size of the calculation
    # step = 0.25
    
    #This adds 0's on the outside of the curve so the line graph goes to the 0 x-axis
    lowerBorder = data.frame(Mass.Charge = min - step, Height = 0)
    upperBorder = data.frame(Mass.Charge = cpgMass + plusMinus + (step*2), Height = 0)
    gaussData = rbind(gaussData, lowerBorder, upperBorder)
    
    for (i in seq(-plusMinus, plusMinus, step)){
      newHeight = cpgHeight*exp(1)^-((8*log(2)*(i^2))/cpgWidth^2)
      min = min + step
      row = data.frame(Mass.Charge=min, Height=newHeight)
      gaussData = rbind(gaussData, row)
    }
  }
  value = gaussData
}


##############################################################################################



##############################################################################################
calcTotalMgmtMethylation = function(sampleID){
  mgmtMethylationValue <- queryDBwithString(con, paste("SELECT VALUE from METHYLATION WHERE SAMPLE_ID LIKE '", sampleID, "' AND AMPLICON LIKE 'MGMT'", sep = ""))

  mgmt1 = mgmtMethylationValue$VALUE[1]*100*2
  mgmt2 = mgmtMethylationValue$VALUE[2]*100*3
  mgmt3 = mgmtMethylationValue$VALUE[3]*100*2
  mgmt4 = mgmtMethylationValue$VALUE[4]*100*3
  mgmt5 = mgmtMethylationValue$VALUE[5]*100*2
  
  sarMethylationValue <- queryDBwithString(con, paste("SELECT VALUE from METHYLATION WHERE SAMPLE_ID LIKE '", sampleID, "' AND AMPLICON LIKE 'SAR'", sep = ""))
  
  sar1 = sarMethylationValue$VALUE[2]*100*2
  sar2 = sarMethylationValue$VALUE[3]*100
  sar3 = sarMethylationValue$VALUE[4]*100
  
  totalMethylation = toString(round((mgmt1 + mgmt2 + mgmt3 + mgmt4 + mgmt5 + sar1 + sar2 + sar3) / 16))
}


##############################################################################################



##############################################################################################
calcTotalMlh1Methylation = function(sampleID){
  mlh1MethylationValue <- queryDBwithString(con, paste("SELECT VALUE from METHYLATION WHERE SAMPLE_ID LIKE '", sampleID, "' AND AMPLICON LIKE 'MLH1'", sep = ""))
  
  mlh11 = mlh1MethylationValue$VALUE[2]*100
  mlh12 = mlh1MethylationValue$VALUE[6]*100*2
  mlh13 = mlh1MethylationValue$VALUE[8]*100
  mlh14 = mlh1MethylationValue$VALUE[9]*100

  totalMethylation = toString(round((mlh11 + mlh12 + mlh13 + mlh14) / 5))
}
##############################################################################################




##############################################################################################
getMGMTControlPassFail = function(){
  #Make sure all the controls are present on the run.
  negWell <- queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'NEG%' AND AMPLICON_ID LIKE 'MGMT'", sep = ""))
  pos20Well <- queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'POS_20%' AND AMPLICON_ID LIKE 'MGMT'", sep = ""))
  pos100Well <- queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'POS_100%' AND AMPLICON_ID LIKE 'MGMT'", sep = ""))
  ntcWell <- queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'NTC%' AND AMPLICON_ID LIKE 'MGMT'", sep = ""))
  
  if (nrow(negWell) != 1 || nrow(pos20Well) != 1 || nrow(pos100Well) != 1 || nrow(ntcWell) != 1){
    message = "One or more of the MGMT control samples is either missing, mislabeled, or there are duplicate controls.  Control data may be inaccurate."
    winDialog(type = c("ok"), message)
    c("Review", "Review", "Review", "Review")
  } else{
    # negControl = as.integer(calcTotalMgmtMethylation("NEG%"))
    # pos20Control = as.integer(calcTotalMgmtMethylation("POS_20%"))
    # pos100Control = as.integer(calcTotalMgmtMethylation("POS_100%"))
    # warningFlag <- toString(queryDBwithString(con, paste("SELECT Status from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'NTC%' AND AMPLICON_ID LIKE 'MGMT'", sep = "")))
    # npMean <- queryDBwithString(con, paste("SELECT np_mean from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'NTC%' AND AMPLICON_ID LIKE 'MGMT'", sep = ""))$np_mean
    # 
    # negPass = if (negControl<10) "Pass" else "Fail"
    # pos20Pass = if (pos20Control>=10 & pos20Control<=40) "Pass" else "Fail"
    # pos100Pass = if (pos100Control>=90) "Pass" else "Fail"
    # ntcPass = if (warningFlag == "bad recal." || as.integer(npMean)<10) "Pass" else "Fail"
    # c(negPass, pos20Pass, pos100Pass, ntcPass)
    
    
    negPass = getMGMTControlQC(myData, "NEG", 0, 10, ">=", "<")
    pos20Pass = getMGMTControlQC(myData, "POS_20", 10, 40, ">=", "<=")
    pos100Pass = getMGMTControlQC(myData, "POS_100", 90, 100, ">=", "<=")
    ntcPass = getMGMTControlQC(myData, "NTC", 0, 0, "N/A", "N/A")
    c(negPass, pos20Pass, pos100Pass, ntcPass)
  }
}
##############################################################################################
getMGMTControlQC = function(myData, sample, min, max, op1, op2){
  wellID <- toString(queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE '", sample, "%' AND AMPLICON_ID LIKE 'MGMT'", sep = "")))
  sampleID = toString(queryDBwithString(con, paste("SELECT SAMPLE_ID from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
  npMean <- toString(queryDBwithString(con, paste("SELECT np_mean from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
  numSpectra <- toString(queryDBwithString(con, paste("SELECT Num_Spectra from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
  warningFlag <- toString(queryDBwithString(con, paste("SELECT Status from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
  totalMeth = as.integer(calcTotalMgmtMethylation(sampleID))

  #Subset data for a single sample
  sampleData = myData[ which(myData$Well==wellID), ]
  
  #Subset data based on Mass/Charge for individual regions requested
  data1 = sampleData[sampleData[3]>"4145" & sampleData[3]<"4200", ]
  data2 = sampleData[sampleData[3]>"4765" & sampleData[3]<"4835", ]
  data3 = sampleData[sampleData[3]>"2545" & sampleData[3]<"2590", ]
  data4 = sampleData[sampleData[3]>"5300" & sampleData[3]<"5375", ]
  data5 = sampleData[sampleData[3]>"3160" & sampleData[3]<"3210", ]
  
  intensity1 = round(colSums(data1[4]), digits = 0)
  intensity2 = round(colSums(data2[4]), digits = 0)
  intensity3 = round(colSums(data3[4]), digits = 0)
  intensity4 = round(colSums(data4[4]), digits = 0)
  intensity5 = round(colSums(data5[4]), digits = 0)
  
  SARpass = getSARControlQC(myData, sample)
  
  if (sample == "NTC"){
    npMeanPass = if(as.integer(npMean) < 10) "Pass" else "Fail"
    intensityPass = if(intensity1 < 400 & intensity2 < 400 & intensity3 < 400) "Pass" else "Fail"
    MGMTPass = if (npMeanPass == "Pass" & intensityPass == "Pass") "Pass" else "Fail"
    if (npMeanPass == "Pass" & intensityPass == "Fail"){
      MGMTPass = "Review"
    }
    if (SARpass == "Review"){
      MGMTPass = "Review"
    }
    MGMTPass
  } else{
  npMeanPass = if(as.integer(npMean) >= 10) "Pass" else "Fail"
  numSpectraPass = if(numSpectra >= 3) "Pass" else "Fail"
  intensityPass = if(strtoi(intensity1) >= 600 & strtoi(intensity2) >= 600 & strtoi(intensity3) >= 600 & strtoi(intensity4) >= 600 & strtoi(intensity5) >= 600) "Pass" else "Fail"
  warningPass = if(warningFlag == "ok") "Pass" else "Fail"
  methPass = if (do.call(op1, list(totalMeth, min)) && do.call(op2, list(totalMeth, max))) "Pass" else "Fail"
  
  #Detect if all MGMT metrics pass
  MGMTPass = if (npMeanPass == "Pass" & numSpectraPass == "Pass" & warningPass == "Pass" & intensityPass == "Pass" & methPass == "Pass" & SARpass == "Pass") "Pass" else "Fail"
  }
}
##############################################################################################
getSARControlQC = function(myData, sample){
  wellID <- toString(queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE '", sample, "%' AND AMPLICON_ID LIKE 'SAR'", sep = "")))
  sampleID = toString(queryDBwithString(con, paste("SELECT SAMPLE_ID from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
  npMean <- toString(queryDBwithString(con, paste("SELECT np_mean from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
  numSpectra <- toString(queryDBwithString(con, paste("SELECT Num_Spectra from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
  warningFlag <- toString(queryDBwithString(con, paste("SELECT Status from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
  
  #Subset data for a single sample
  sampleData = myData[ which(myData$Well==wellID), ]
  
  data1 = sampleData[sampleData[3]>"2865" & sampleData[3]<"2920", ]
  data2 = sampleData[sampleData[3]>"2540" & sampleData[3]<"2575", ]
  data3 = sampleData[sampleData[3]>"3490" & sampleData[3]<"3525", ]
  
  intensity1 = round(colSums(data1[4]), digits = 0)
  intensity2 = round(colSums(data2[4]), digits = 0)
  intensity3 = round(colSums(data3[4]), digits = 0)
  
  
  if (sample == "NTC"){
    npMeanPass = if(as.integer(npMean) < 10) "Pass" else "Fail"
    intensityPass = if(intensity1 < 400 & intensity2 < 400 & intensity3 < 400) "Pass" else "Fail"
    SARpass = if (npMeanPass == "Pass" & intensityPass == "Pass") "Pass" else "Fail"
    if (npMeanPass == "Pass" & intensityPass == "Fail"){
      SARpass = "Review"
    }
    SARpass
  }
  else{
    npMeanPass = if(as.integer(npMean) >= 10) "Pass" else "Fail"
    numSpectraPass = if(numSpectra >= 3) "Pass" else "Fail"
    intensityPass = if(intensity1 >= 600 & intensity2 >= 600 & intensity3 >= 600) "Pass" else "Fail"
    warningPass = if(warningFlag == "ok") "Pass" else "Fail"
    
    
    #Determine if SAR QC metrics pass
    SARpass = if (npMeanPass == "Pass" & numSpectraPass == "Pass" & warningPass == "Pass" & intensityPass == "Pass") "Pass" else "Fail"
  }
}
##############################################################################################



##############################################################################################
getMLH1ControlPassFail = function(myData){
  
  #Make sure all the controls are present on the run.
  negWell <- queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'NEG%' AND AMPLICON_ID LIKE 'MLH1'", sep = ""))
  pos20Well <- queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'POS_20%' AND AMPLICON_ID LIKE 'MLH1'", sep = ""))
  pos100Well <- queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'POS_100%' AND AMPLICON_ID LIKE 'MLH1'", sep = ""))
  ntcWell <- queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE 'NTC%' AND AMPLICON_ID LIKE 'MLH1'", sep = ""))
  
  if (nrow(negWell) != 1 || nrow(pos20Well) != 1 || nrow(pos100Well) != 1 || nrow(ntcWell) != 1){
    message = "One or more of the MLH1 control samples is either missing, mislabeled, or there are duplicate controls.  Control data may be inaccurate."
    winDialog(type = c("ok"), message)
    c("Review", "Review", "Review", "Review")
  } else{
    negPass = getMLH1ControlQC(myData, "NEG", 0, 10, ">=", "<")
    pos20Pass = getMLH1ControlQC(myData, "POS_20", 10, 40, ">=", "<=")
    pos100Pass = getMLH1ControlQC(myData, "POS_100", 90, 100, ">=", "<=")
    ntcPass = getMLH1ControlQC(myData, "NTC", 0, 0, "N/A", "N/A")
    c(negPass, pos20Pass, pos100Pass, ntcPass) 
  }
}

##############################################################################################
getMLH1ControlQC = function(myData, sample, min, max, op1, op2){
  wellID <- queryDBwithString(con, paste("SELECT WELL_POSITION from SPECTRAQUALITY WHERE SAMPLE_ID LIKE '", sample, "%' AND AMPLICON_ID LIKE 'MLH1'", sep = ""))
  sampleID = toString(queryDBwithString(con, paste("SELECT SAMPLE_ID from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
  npMean <- toString(queryDBwithString(con, paste("SELECT np_mean from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
  numSpectra <- toString(queryDBwithString(con, paste("SELECT Num_Spectra from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "' LIMIT 1", sep = "")))
  warningFlag <- toString(queryDBwithString(con, paste("SELECT Status from SPECTRAQUALITY WHERE WELL_POSITION LIKE '", wellID, "'", sep = "")))
  totalMeth = as.integer(calcTotalMlh1Methylation(paste(sample, "%", sep = "")))
  
  #Subset data for a single sample
  sampleData = myData[ which(myData$Well==wellID$WELL_POSITION), ]
  
  data1 = sampleData[sampleData[3]>"2580" & sampleData[3]<"2620", ]
  data2 = sampleData[sampleData[3]>"1880" & sampleData[3]<"1930", ]
  data3 = sampleData[sampleData[3]>"2135" & sampleData[3]<"2165", ]
  data4 = sampleData[sampleData[3]>"2870" & sampleData[3]<"2905", ]
  
  intensity1 = round(colSums(data1[4]), digits = 0)
  intensity2 = round(colSums(data2[4]), digits = 0)
  intensity3 = round(colSums(data3[4]), digits = 0)
  intensity4 = round(colSums(data4[4]), digits = 0)
  
  if (sample == "NTC"){
    npMeanPass = if(as.integer(npMean) < 10) "Pass" else "Fail"
    intensityPass = if(intensity1 < 400 & intensity2 < 400 & intensity3 < 400 & intensity4 < 400) "Pass" else "Fail"
    MLH1Pass = if (npMeanPass == "Pass" & intensityPass == "Pass") "Pass" else "Fail"
    if (npMeanPass == "Pass" & intensityPass == "Fail"){
      MLH1Pass = "Review"
    }
    MLH1Pass
  }
  else{
    npMeanPass = if(as.integer(npMean) >= 10) "Pass" else "Fail"
    numSpectraPass = if(numSpectra >= 3) "Pass" else "Fail"
    intensityPass = if(strtoi(intensity1) >= 600 & strtoi(intensity2) >= 600 & strtoi(intensity3) >= 600 & strtoi(intensity4) >= 600) "Pass" else "Fail"
    warningPass = if(warningFlag == "ok") "Pass" else "Fail"
    methPass = if (do.call(op1, list(totalMeth, min)) && do.call(op2, list(totalMeth, max))) "Pass" else "Fail"
    
    #Detect if all MLH1 metrics pass
    MLH1Pass = if (npMeanPass == "Pass" & numSpectraPass == "Pass" & warningPass == "Pass" & intensityPass == "Pass" & methPass == "Pass") "Pass" else "Fail"
  }
}
##############################################################################################



