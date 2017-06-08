###################################################################################################################
# notes: modifications are numbered starting with 1, experiments starting with 0
templateFile <- readLines ("M:\\r_scripts\\r_scripts_for_different_pipelines\\fragmentation_investigation\\TMS2Independent_LUMOS.xml")


###################################################################################################################
# parameters specification


mzStandardProteins <- list (
  

  "GST"           = list ("mass1" = c(799.74, 10), # 34
                          "mass2" = c(906.20, 10), # 30
                          "mass3" = c(1045.46, 10)), # 26
  
  "myoglobin"    = list ("mass1" = c(738.01, 10),
                         "mass2" = c(808.15, 10),
                         "mass3" = c(893.16, 10)),
  
  "h2a"     =     list ("mass1" = c(609.21, 10),  # 23
                        "mass2" = c(700.45, 10), # 20
                        "mass3" = c(823.95, 10)), # 17
  
  "h2b"     =     list ("mass1" = c(600.51, 10),  # 23
                        "mass2" = c(690.39, 10), # 20
                        "mass3" = c(812.10, 10)), # 17
  
  "h3_1"     =     list ("mass1" = c(611.87, 10),# 25
                         "mass2" = c(664.98, 10),# 23
                         "mass3" = c(728.17, 10)), # 21 
  
  "h3_3"     =     list ("mass1" = c(608.87, 10),# 25
                         "mass2" = c(691.76, 10),# 22
                         "mass3" = c(800.77, 10)), # 19 
  
  "h4"     =      list ("mass1" = c(625.19, 10),# 18
                        "mass2" = c(750.03, 10),# 15
                        "mass3" = c(937.29, 10)), # 12
  
  
  "h33tail"     =     list ("mass1" = c(446.10, 10), # 12
                            "mass2" = c(486.56, 10), # 11
                            "mass3" = c(535.11, 10), # 10
                            "mass4" = c(594.46, 10), # 9
                            "mass5" = c(668.64, 10)), # 8
  
  
  "c345c_c3"     =      list ("mass1" = c(950.51, 10),# 19
                              "mass2" = c(1062.22, 10),# 17
                              "mass3" = c(1203.78, 10)) # 15
  
)


###################################################################################################################
# FUNCTIONS

modifyExperiment <- function (paramCombTable, expIndex, massList){
  
  CollisionEnergyVal <- NULL
  CollisionEnergyType <- "ETD"
  modReturn <- gsub ("substitute", expIndex, c("<Experiment ExperimentIndex=\"substitute\">"))  
  modReturn <- c (modReturn, "<TMSnScan>")
  
  for (i in 1:ncol(paramCombTable)){
    if (is.na (paramCombTable[, i])) next()
    if (colnames (paramCombTable)[i]  %in% c( "CollisionEnergyCID", "CollisionEnergyHCD")){
      CollisionEnergyVal <- paramCombTable[, i]
      CollisionEnergyType <- colnames (paramCombTable)[i]
      next()
    }
    buffer <-  gsub ("substitute", colnames (paramCombTable)[i], "<substitute>value</substitute>")
    buffer <-  gsub ("value", paramCombTable[1, i], buffer)
    modReturn  <- c (modReturn ,buffer )
  }
  
  if (CollisionEnergyType == "ETD") modReturn <- c (modReturn, "<MassList>")
  if (CollisionEnergyType == "CollisionEnergyCID") modReturn <- c (modReturn, "<MassList CollisionEnergyCID = \"true\">" )
  if (CollisionEnergyType == "CollisionEnergyHCD") modReturn <- c (modReturn, "<MassList CollisionEnergyHCD = \"true\">" )
  
    
  for (i in 1:length (massList)){
    modReturn <- c (modReturn, "<MassListRecord>")
    modReturn <- c (modReturn, gsub ("substitute", massList[[i]][1], "<MOverZ>substitute</MOverZ>"))
    modReturn <- c (modReturn, gsub ("substitute", massList[[i]][2], "<Z>substitute</Z>"))
    if (CollisionEnergyType == "CollisionEnergyCID") modReturn <- c (modReturn, gsub ("substitute", CollisionEnergyVal, "<CollisionEnergyCID>substitute</CollisionEnergyCID>"))
    if (CollisionEnergyType == "CollisionEnergyHCD") modReturn <- c (modReturn, gsub ("substitute", CollisionEnergyVal, "<CollisionEnergyHCD>substitute</CollisionEnergyHCD>"))
    modReturn <- c (modReturn, "</MassListRecord>")
  }
  
  modReturn <- c (modReturn, "</MassList>", "</TMSnScan>", "</Experiment>", "")
  return (modReturn)
}


createNewMethod <- function (parMS1, parMS2, massList, parPeriod, 
                             replPar = 1, newMethodName = "newMethod.xml", 
                             simplifyEThciD = TRUE, MS2perMS1 = 10000, labelMS2ByMass = TRUE){

  # make params table
  paramCombTable <- expand.grid (parMS2, stringsAsFactors = FALSE)
  paramCombTable <- paramCombTable[!(paramCombTable$ETDReactionTime == 0 &
                                     paramCombTable$ETDSupplementalActivationEnergy == 0),  ]
  # simplify 
  if (simplifyEThciD){
    # subsititute ETD with 0 fragmentation time with HCD and CID fragmenation
    subETD <- which (paramCombTable$ETDReactionTime == 0)
    if (length (subETD) > 0){
      valToNA <- c ("ETDReactionTime", "ETDSupplementalActivation", "ETDReagentTarget", "ETDSupplementalActivationEnergy")
      valToNA <- valToNA[valToNA %in% names (parMS2)]
      for (i in 1:length (subETD)){
        if (paramCombTable[subETD[i], "ETDSupplementalActivation"] == "ETciD"){
          ActivationType <- "CID"
          paramCombTable[subETD[i], "CollisionEnergyCID"] <- paramCombTable[subETD[i],  "ETDSupplementalActivationEnergy"]
        }  
        if (paramCombTable[subETD[i], "ETDSupplementalActivation"] == "EThcD"){
          ActivationType <- "HCD"
          paramCombTable[subETD[i], "CollisionEnergyHCD"] <- paramCombTable[subETD[i],  "ETDSupplementalActivationEnergy"]
        } 
        paramCombTable[subETD[i], "ActivationType"] <- ActivationType
        paramCombTable[subETD[i], valToNA] <- NA
      }
    }
    paramCombTable <- paramCombTable[!duplicated (paramCombTable) , ]
    subCID <- which (paramCombTable$ETDSupplementalActivationEnergy == 0)
    paramCombTable[subCID, c("ETDSupplementalActivation", "ETDSupplementalActivationEnergy")] <- NA
    paramCombTable <- paramCombTable[!duplicated (paramCombTable), ]
  }
  
  paramCombTable <- paramCombTable[rep(1:nrow(paramCombTable), times = replPar),]
  
  # assign MS1 or MS2
  MS2Groups <- split(1:nrow (paramCombTable), ceiling(1:nrow (paramCombTable)/MS2perMS1))
  experimentType <- lapply (MS2Groups, function (x) c ("MS1", rep ("MS2", length (x))))
  experimentType <- as.character (unlist (experimentType))
  
 
  ############################## 
  # start creating method
  
  # 0) declare the number of modification
  numberMod <- 1 
  
  # 1) MS1 parameters specification
  newMethod <- c( "<?xml version=\"1.0\" encoding=\"utf-8\" ?>", "<MethodModifications Version=\"1\" Model=\"OrbitrapFusion\" Family=\"Calcium\" Type=\"SL\">")
  
  newMethod <- c (newMethod, 
                  "<Modification Order=\"1\">",                                                               
                  "<Experiment ExperimentIndex=\"0\">",                                                     
                  "<FullMSScan>")
  
  for (i in 1:length (parMS1 )){
    buffer <- gsub ("substitute", names (parMS1)[i], "<substitute>value</substitute>")
    buffer <- gsub ("value", parMS1 [[i]], buffer)
    newMethod <- c (newMethod, buffer)
  }
  
  newMethod <- c (newMethod,  "</FullMSScan>" ,  "</Experiment>",  "</Modification>", "")
  
  # 2) create new experiments
  
  addPatternMS1 <- c ("<Modification Order=\"modSub\">", 
                      "<CopyAndAppendExperiment SourceExperimentIndex=\"0\"/>",
                      "</Modification>", "")
  
  addPatternMS2 <- c ("<Modification Order=\"modSub\">", 
                      "<CopyAndAppendExperiment SourceExperimentIndex=\"1\"/>",
                      "</Modification>", "")
    
  for (i in 3:length (experimentType)){
    numberMod <-  numberMod + 1
    if (experimentType[i] == "MS1") newMethod <- c (newMethod,   gsub ("modSub", numberMod, addPatternMS1))
    if (experimentType[i] == "MS2") newMethod <- c (newMethod,   gsub ("modSub", numberMod, addPatternMS2)) 
  }
  
  # 3) split time  
  StartTime <- c()
  EndTime   <- c()
  
  addPatternTime <- c ("<Modification Order=\"modSub\">",
                       "<Experiment ExperimentIndex=\"expSub\">",                                                     
                       "<StartTimeMin>valueStart</StartTimeMin>",                                                                 
                       "<EndTimeMin>valueEnd</EndTimeMin>",
                       "</Experiment>",  "</Modification>", "")
  
  
  for (i in 0:(length (experimentType) - 1)){
    numberMod <-  numberMod + 1
    # specify mod number
    buffer <- gsub ("modSub", numberMod, addPatternTime)
    # specify experiment number
    buffer <- gsub ("expSub", i, buffer)  
    # specify start and end
    buffer <- gsub ("valueStart",  i * parPeriod + 0.01, buffer)
    buffer <- gsub ("valueEnd", (i + 1) * parPeriod, buffer)
    newMethod <- c (newMethod, buffer)
    if (i == (length (experimentType) - 1)) message (paste ("the method is", (i + 1) * parPeriod, "min long"))
  }
  
  # 4) modify experiments  
  newMethod <- c (newMethod, sub ("substitute", numberMod + 1, "<Modification Order=\"substitute\">"))
  parNumber <- 0
  
  for (i in 1:length (experimentType)){
    if (experimentType[i] == "MS2"){
      parNumber <- parNumber + 1
      if (labelMS2ByMass) for (k in 1:length (massList)) massList[[k]][1] <- round (massList[[k]][1], 1) + parNumber/10000
      newMethod <- c (newMethod, modifyExperiment (paramCombTable[parNumber, ], expIndex = i-1, massList))
    } 
  }
  
  
  newMethod <- c (newMethod, "</Modification>", "</MethodModifications>")  
  
  writeLines (newMethod, newMethodName)
}
 

createNewMethods <- function (parMS1, parMS2, massList, protName = "", parMS2BetweenRuns, 
                              splitMassList = TRUE, parPeriod = 0.3, replPar = 1, simplifyEThciD = TRUE){
  
  paramCombMS2BR <- t(expand.grid (parMS2[names (parMS2) %in% parMS2BetweenRuns], stringsAsFactors = FALSE))
  parMS2 <- parMS2[!(names (parMS2) %in% parMS2BetweenRuns)]
  
  for (i in 1:ncol (paramCombMS2BR)){
    parMS2ToAdd <- setNames (list(paramCombMS2BR[i]), row.names (paramCombMS2BR))
    if (splitMassList){
      for (k in 1:length (massList)){
        newMethodName <- paste0 (c(protName, as.integer (massList[[k]][1]), paste(row.names (paramCombMS2BR), paramCombMS2BR[, i], sep = "_"), ".xml"), collapse = "_") 
        createNewMethod (parMS1 = parMS1, parMS2 = c(parMS2, parMS2ToAdd),
                         massList =  list("mass1" = massList[[k]]),
                         parPeriod = parPeriod,
                         replPar = replPar, newMethodName = newMethodName,
                         simplifyEThciD = simplifyEThciD)
      }
    } else {
      newMethodName <- paste0 (c(protName, paste(row.names (paramCombMS2BR), paramCombMS2BR[, i], sep = "_"), ".xml"), collapse = "_") 
      createNewMethod (parMS1 = parMS1 , parMS2 = c(parMS2, parMS2ToAdd),
                       massList = massList,
                       parPeriod = parPeriod,
                       replPar = replPar, newMethodName = newMethodName,
                       simplifyEThciD = simplifyEThciD)
    }
  }
  
}