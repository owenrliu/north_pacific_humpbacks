
library(mgcv)

# =================================Catch Data=========================================
ReadCatches <- function(DatFile, BreedingOpt, FeedingOpt, Nbreed, Nfeed, BreedNames, FeedNames, Yr1, Yr2, CatchSer, ByCatchFile = "Null.csv") {
  # Read in the catches
  Index <- which(DatFile[, 1] == "Total_areas_in_catch_files:")
  MaxAreas <- as.numeric(DatFile[Index, 2])
  Index <- which(DatFile[, 1] == "First_year_with_catches:")
  Y1Cat <- as.numeric(DatFile[Index, 2])
  if (CatchSer == "L") Offset <- -2
  if (CatchSer == "B") Offset <- -1
  if (CatchSer == "H") Offset <- 0

  Index <- which(DatFile[, 1] == "Columns_with_breeding_gound_catches:" & DatFile[, 2] == BreedingOpt) + 1
  CatchLinkB <- matrix(0, nrow = Nbreed, ncol = MaxAreas)
  rownames(CatchLinkB) <- BreedNames
  for (Ibreed in 1:Nbreed) {
    for (Iarea in 1:MaxAreas) CatchLinkB[Ibreed, Iarea] <- as.numeric(DatFile[Index + Ibreed, Iarea])
  }
  Index <- which(DatFile[, 1] == "Columns_with_feeding_gound_catches:" & DatFile[, 2] == FeedingOpt) + 1
  CatchLinkF <- matrix(0, nrow = Nfeed, ncol = MaxAreas)
  rownames(CatchLinkF) <- FeedNames
  for (Ifeed in 1:Nfeed) {
    for (Iarea in 1:MaxAreas) CatchLinkF[Ifeed, Iarea] <- as.numeric(DatFile[Index + Ifeed, Iarea])
  }
  Tnames <- c("Asia", "Hawaii", "Mexico", "MX_AR", "MX_ML", "Central Am", "RUS+WAL", "EAL+BER", "WGOA", "NGOA", "SEA-NBC", "SBC+WA", "OR+CA")
  colnames(CatchLinkF) <- Tnames
  colnames(CatchLinkB) <- Tnames

  # Commercial, bycatch and ship strikes
  CatchCom <- read.csv(here("data", "CatchCom.csv"), skip = 1, head = F)[, -1]
  CatchBy <- read.csv(paste0(here("data"), "/", ByCatchFile), skip = 1, head = F)[, -1]

  # Now allocate catches by breeding and feeding groups
  CatchB <- matrix(0, nrow = (Yr2 - Yr1 + 1), ncol = Nbreed)
  CatchF <- matrix(0, nrow = (Yr2 - Yr1 + 1), ncol = Nfeed)
  for (Iyr in Yr1:Yr2)
  {
    Jyr <- Iyr - Yr1 + 1
    Kyr <- Iyr - Y1Cat + 1
    for (Ibreed in 1:Nbreed) {
      for (Iarea in 1:MaxAreas) {
        if (CatchLinkB[Ibreed, Iarea] == 1) CatchB[Jyr, Ibreed] <- CatchB[Jyr, Ibreed] + CatchCom[Kyr, 3 * Iarea + Offset] + CatchBy[Kyr, 3 * Iarea + Offset]
      }
    }
    for (Ifeed in 1:Nfeed) {
      for (Iarea in 1:MaxAreas) {
        if (CatchLinkF[Ifeed, Iarea] == 1) CatchF[Jyr, Ifeed] <- CatchF[Jyr, Ifeed] + CatchCom[Kyr, 3 * Iarea + Offset] + CatchBy[Kyr, 3 * Iarea + Offset]
      }
    }
  }

  # FUDGING THE LAST YEAR BECAUSE OF MISSING DATA---REMOVE THIS LATER
  # for(i in 1:nrow(CatchB)){
  #   if(all(is.na(CatchB[i,]))) CatchB[i,] <- CatchB[i-1,]
  # }
  # for(i in 1:nrow(CatchF)){
  #   if(all(is.na(CatchF[i,]))) CatchF[i,] <- CatchF[i-1,]
  # }
  
  Outs <- NULL
  Outs$CatchF <- CatchF
  Outs$CatchB <- CatchB
  return(Outs)
}

# =================================Survey Data=========================================
# Function to read the humpback survey data
ReadSurveyData <- function(Code, BreedingOpt, FeedingOpt, BreedNames, FeedNames, Yr1, Yr2, FullDiag, SensTest = "") {
  
  Surveys <- read.csv(here("data", "survey_data_updated_04_2026_B2bF1.csv"), fill = T, comment.char = "?", header = T, row.names = NULL)[, 1:13]
  colnames <- c("Year1", "Year2", "Estimate", "CV", "Area", "Rel", "Use", "Add.cv", "Hypothesis", "Class", "SensUse", "Component", "Reference")
  colnames(Surveys) <- colnames
  # print(head(Surveys))
  if (FullDiag == T) cat("Initial number of lines", length(Surveys[, 1]), "\n")

  # Extract the surveys to use
  # Index <- which (Surveys$Use == "Yes")
  Index <- which((Surveys$Use == "Yes" | Surveys$Use == SensTest) & Surveys$SensUse != SensTest)
  Surveys <- Surveys[Index, ]
  SearchString <- c("All", paste0(BreedingOpt, FeedingOpt), paste0(BreedingOpt, " only"), paste0(FeedingOpt, " only"))
  if (FullDiag == T) print(SearchString)
  Index <- which(Surveys$Hypothesis %in% SearchString)
  Surveys <- Surveys[Index, ]
  if (FullDiag == T) cat("Final number of lines", length(Surveys[, 1]), "\n")

  # Deal with incorrect names
  if (BreedingOpt == "B1") {
    Index <- which(Surveys$Area == "MX_ML")
    Surveys$Area[Index] <- "Mexico"
  }

  Index <- (Surveys$Add.cv == "No" | Surveys$Add.cv == "Maybe")
  Surveys$Add.cv[Index] <- 0

  # check if the filepath exists; if not, create it
  if (!dir.exists(here("Diags"))) {
    dir.create(here("Diags"))
  }
  if (!dir.exists(here("Diags", "SurveyUse"))) {
    dir.create(here("Diags", "SurveyUse"))
  }

  write.csv(Surveys, file = paste0(here("Diags", "SurveyUse"), "/", Code, SensTest, ".csv"), row.names = F)

  Surveys <- Surveys |>
    mutate(across(all_of(c("Year1", "Year2", "Estimate", "CV", "Add.cv", "Component")), as.numeric))

  # Now process the data
  Surveys <- Surveys |>
    # Indicate type (total, breeding, or feeding estimate)
    mutate(Type = case_when(
      Area == "Total" ~ 1,
      Area %in% BreedNames ~ 2,
      Area %in% FeedNames ~ 3,
      TRUE ~ -1
    )) |>
    # Indicate which zone/area the estimate is for
    mutate(
      AreaBreed = match(Area, BreedNames),
      AreaFeed = match(Area, FeedNames),
      AreaTotal = ifelse(Area == "Total", -1, NA)
    ) |>
    mutate(Area = coalesce(AreaBreed, AreaFeed, AreaTotal)) |>
    # Convert years into model years
    mutate(Year1 = Year1 - Yr1 + 1, Year2 = Year2 - Yr1 + 1, Class = Class, used = 1)
  SurveyI <- Surveys |>
    dplyr::select(all_of(c("Year1", "Year2", "Type", "Area", "Class", "used", "Add.cv", "Component"))) |>
    as.matrix()
  SurveyR <- Surveys |>
    dplyr::select(all_of(c("Estimate", "CV"))) |>
    as.matrix()
  # Number of additional variances
  NextraCV1 <- max(SurveyI[, "Add.cv"])
  # Number of data points
  NsurveyData <- length(SurveyI[, 1])
  # Number of survey series
  NsurveySeries <- length(unique(SurveyI[, 5]))
  SurveySeries <- c(1, rep(2, NsurveySeries - 1))

  Outs <- NULL
  Outs$SurveyI <- SurveyI
  Outs$SurveyR <- SurveyR
  Outs$NsurveyData <- NsurveyData
  Outs$NsurveySeries <- NsurveySeries
  Outs$SurveySeries <- SurveySeries
  Outs$NextraCV1 <- NextraCV1
  return(Outs)
}

# =================================Mixing Data=========================================
# Function to read the breeding-feeding mixing data
ReadMixingData <- function(Code, DatFile, Yr1, Yr2, Nbreed, Nfeed, BreedNames, FeedNames, MixWeights, FullDiag, SensTest = "", MaxN = -1000) {
  # Maximum effective sample size
  # MaxN <- 100

  # Read in mixing data
  MixFile1 <- read.csv(here("data", "Genetics_mixing_data_allscenarios_table_Long.csv"))
  MixFile2 <- read.csv(here("data", "Mark-Recapture_mixing_data_allscenarios_table_Long_B2bF1.csv"))
  MixFile <- rbind(MixFile1, MixFile2)
  Index <- which(MixFile$Hypothesis == Code)
  MixFile <- MixFile[Index, ]

  Type <- c("Mark-Recapture", "Genetics")
  NmixData <- 0
  Nprop <- rep(0, 4)
  Index <- which(DatFile[, 2] == "Minimum_CV_for_the_mixing_data")
  MinMixCV <- as.numeric(DatFile[Index + 1, 1])
  Index <- which(DatFile[, 2] == "Minimum_SD_for_the_mixing_data")
  MinMixSD <- as.numeric(DatFile[Index + 1, 1])

  # check if the filepath exists; if not, create it
  if (!dir.exists(here("Diags"))) {
    dir.create(here("Diags"))
  }
  if (!dir.exists(here("Diags", "MixUse"))) {
    dir.create(here("Diags", "MixUse"))
  }
  FileName <- paste0(here("Diags"), "/MixUse/", Code, SensTest, ".csv")

  # Now pull out results (breeding to feeding)
  ObsMixBtoFE <- array(0, dim = c(2, Nbreed, Nfeed))
  ObsMixBtoFP <- array(0, dim = c(2, Nbreed, Nfeed))
  ObsMixBtoFO <- matrix(0, nrow = 2, ncol = Nbreed)
  write("Breeding to feeding", FileName)
  if (FullDiag == T) print("Breeding to feeding mixing")
  for (Itype in 1:2)
  {
    write(Type[Itype], FileName, append = T)
    Index <- MixFile$Method == Type[Itype] & MixFile$Direction == "BreedingtoFeeding"
    MixFile2 <- MixFile[Index, ]
    for (Ibreed in 1:length(BreedNames))
    {
      Total <- 0
      Top <- 0
      Bot <- 0
      for (Ifeed in 1:length(FeedNames))
      {
        Index <- MixFile2$Feeding %in% FeedNames[Ifeed] & MixFile2$Breeding %in% BreedNames[Ibreed]
        Est <- as.numeric(MixFile2$Estimate[Index])
        CV <- as.numeric(MixFile2$CV[Index])
        SD <- CV * Est
        if (is.na(CV)) {
          Lower95 <- as.numeric(MixFile2$lowerCI[Index])
          Upper95 <- as.numeric(MixFile2$upperCI[Index])
          CV <- (Upper95 - Lower95) / (2.0 * 1.96) / Est
          if (Est <= 0) CV <- 0
          SD <- CV * Est
        }
        Top <- Top + Est * (1 - Est)
        Bot <- Bot + SD^2
        if (CV < MinMixCV) CV <- MinMixCV
        SD <- Est * CV
        if (Est > 0 & SD < MinMixSD) SD <- MinMixSD
        if (Est > 0 || SD > 0) {
          NmixData <- NmixData + 1
          Nprop[Itype] <- Nprop[Itype] + 1
        }
        Total <- Total + Est

        # cat(Itype,Ibreed,Ifeed,Est,CV,SD,"\n")
        ObsMixBtoFE[Itype, Ibreed, Ifeed] <- Est
        ObsMixBtoFP[Itype, Ibreed, Ifeed] <- SD
      } # Feed
      ObsMixBtoFO[Itype, Ibreed] <- min(MaxN, Top / Bot * MixWeights[Itype])
      # print(Total)
    } # Breed
    write("Estimates", FileName, append = T)
    write(t(unlist(FeedNames)), FileName, append = T, sep = ",", ncol = Nfeed)
    for (Ibreed in 1:Nbreed)
    {
      xx <- paste(c(ObsMixBtoFE[Itype, Ibreed, ]))
      xx <- c(xx, ",#", unlist(BreedNames[Ibreed]))
      write(xx, FileName, append = T, ncol = Nfeed + 2, sep = ",")
    }
    write("SD", FileName, append = T)
    write(t(unlist(FeedNames)), FileName, append = T, sep = ",", ncol = Nfeed)
    for (Ibreed in 1:Nbreed)
    {
      xx <- paste(c(ObsMixBtoFP[Itype, Ibreed, ]))
      xx <- c(xx, ",#", unlist(BreedNames[Ibreed]))
      write(xx, FileName, append = T, ncol = Nfeed + 2, sep = ",")
    }
  } # Type
  write("EffN", FileName, append = T)
  write(t(unlist(BreedNames)), FileName, append = T, sep = ",", ncol = Nfeed)
  write(t(ObsMixBtoFO), FileName, append = T, ncol = Nbreed, sep = ",")

  # Now pull out results (feeding to breeding)
  ObsMixFtoBE <- array(0, dim = c(2, Nbreed, Nfeed))
  ObsMixFtoBP <- array(0, dim = c(2, Nbreed, Nfeed))
  ObsMixFtoBO <- matrix(0, nrow = 2, ncol = Nfeed)
  write("Feeding to Breeding", FileName, append = T)
  if (FullDiag == T) print("Feeding to breeding mixing")
  for (Itype in 1:2)
  {
    write(Type[Itype], FileName, append = T)
    Index <- MixFile$Method == Type[Itype] & MixFile$Direction == "FeedingtoBreeding"
    MixFile2 <- MixFile[Index, ]
    for (Ifeed in 1:length(FeedNames))
    {
      Total <- 0
      Top <- 0
      Bot <- 0
      for (Ibreed in 1:length(BreedNames))
      {
        Index <- MixFile2$Feeding %in% FeedNames[Ifeed] & MixFile2$Breeding %in% BreedNames[Ibreed]
        Est <- as.numeric(MixFile2$Estimate[Index])
        CV <- as.numeric(MixFile2$CV[Index])
        SD <- CV * Est
        if (is.na(CV)) {
          Lower95 <- as.numeric(MixFile2$lowerCI[Index])
          Upper95 <- as.numeric(MixFile2$upperCI[Index])
          CV <- (Upper95 - Lower95) / (2.0 * 1.96) / Est
          if (Est <= 0) CV <- 0
          SD <- CV * Est
        }
        Top <- Top + Est * (1 - Est)
        Bot <- Bot + SD^2
        if (CV < MinMixCV) CV <- MinMixCV
        SD <- Est * CV
        if (Est > 0 & SD < MinMixSD) SD <- MinMixSD
        if (Est > 0 || SD > 0) {
          NmixData <- NmixData + 1
          Nprop[Itype + 2] <- Nprop[Itype + 2] + 1
        }
        Total <- Total + Est

        # cat(Ibreed,Ifeed,Est,CV,SD,"\n")
        ObsMixFtoBE[Itype, Ibreed, Ifeed] <- Est
        ObsMixFtoBP[Itype, Ibreed, Ifeed] <- SD
      } # Breed
      ObsMixFtoBO[Itype, Ifeed] <- min(MaxN, Top / Bot * MixWeights[Itype])
      # print(Total)
    } # Feed
    write("Estimates", FileName, append = T)
    write(t(unlist(FeedNames)), FileName, append = T, sep = ",", ncol = Nfeed)
    for (Ibreed in 1:Nbreed)
    {
      xx <- paste(c(ObsMixFtoBE[Itype, Ibreed, ]))
      xx <- c(xx, ",#", unlist(BreedNames[Ibreed]))
      write(xx, FileName, append = T, ncol = Nfeed + 2, sep = ",")
    }
    write("SD", FileName, append = T)
    write(t(unlist(FeedNames)), FileName, append = T, sep = ",", ncol = Nfeed)
    for (Ibreed in 1:Nbreed)
    {
      xx <- paste(c(ObsMixFtoBP[Itype, Ibreed, ]))
      xx <- c(xx, ",#", unlist(BreedNames[Ibreed]))
      write(xx, FileName, append = T, ncol = Nfeed + 2, sep = ",")
    }
  } # Type
  write("EffN", FileName, append = T)
  write(t(unlist(FeedNames)), FileName, append = T, sep = ",", ncol = Nfeed)
  write(t(ObsMixFtoBO), FileName, append = T, ncol = Nfeed, sep = ",")
  # print(NmixData)
  # print(Nprop)

  Outs <- NULL
  Outs$ObsMixBtoFE <- ObsMixBtoFE
  Outs$ObsMixBtoFP <- ObsMixBtoFP
  Outs$ObsMixBtoFO <- ObsMixBtoFO
  Outs$ObsMixFtoBE <- ObsMixFtoBE
  Outs$ObsMixFtoBP <- ObsMixFtoBP
  Outs$ObsMixFtoBO <- ObsMixFtoBO
  Outs$NmixData <- NmixData

  return(Outs)
}

# =================================Environmental Data=========================================
## type= "raw" - raw covariates
## type= "anom" - zone-specific anomalies, based on a 1993-2013 mean
extract_env_var <- function(dat, var_name, FeedNames, YrStart, YrEnd, lag, vn, normalize) {
  df <- dat[[var_name]] |> mutate(zone = str_remove_all(zoneID, "_\\d*"))
  mat <- map(FeedNames, \(x) {
    df |> filter(zone == x, year >= YrStart - lag, year <= YrEnd - lag) |> pull(vn) |> as.numeric()
  }) |> do.call(rbind, args = _)
  if (normalize) mat <- (mat - mean(mat)) / sd(mat)
  mat
}

ReadGLORYS <- function(DatFile, FeedingOpt, YrStart = 2000, YrEnd = 2024, type,lag, normalize) {
  # Do you want to make raw values or anomalies?
  vn <- ""
  if (type == "raw") {
    vn <- "value"
    envdat <- read_rds(here("data", "processed covariates", "glorys_covars_by_zone.rds"))
  }
  if (type == "anom") {
    vn <- "anom"
    envdat <- read_rds(here("data", "processed covariates", "glorys_anomalies_by_zone.rds"))
  }

  # pull the right timeseries for feeding grounds
  which_feed <- which(DatFile[, 1] == "Feeding_grounds" & DatFile[, 2] == FeedingOpt)
  Nfeed <- as.numeric(DatFile[which_feed, 3])
  FeedNames <- as.character(DatFile[which_feed + 1, 1:Nfeed])
  
  sstMat <- extract_env_var(envdat,"sst",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  chlMat <- extract_env_var(envdat,"chl",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  mldMat <- extract_env_var(envdat,"mld",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  no3Mat <- extract_env_var(envdat,"no3",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  nppvMat <- extract_env_var(envdat,"nppv",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  
  Outs <- list(
    sst = sstMat,
    chl = chlMat,
    mld=mldMat,
    no3=no3Mat,
    nppv=nppvMat
  )
  Outs
}

ReadMOM6 <- function(DatFile, FeedingOpt, YrStart = 2000, YrEnd = 2024,type,lag, normalize) {
  # Do you want to make raw values or anomalies?
  vn <- ""
  if (type == "raw") {
    vn <- "value"
    envdat <- read_rds(here("data", "processed covariates", "mom6_covars_by_zone.rds"))
  }
  if (type == "anom") {
    vn <- "anom"
    envdat <- read_rds(here("data", "processed covariates", "mom6_anomalies_by_zone.rds"))
  }
  
  # pull the right timeseries for feeding grounds
  which_feed <- which(DatFile[, 1] == "Feeding_grounds" & DatFile[, 2] == FeedingOpt)
  Nfeed <- as.numeric(DatFile[which_feed, 3])
  FeedNames <- as.character(DatFile[which_feed + 1, 1:Nfeed])
  
  # out <- list(sst=sst_anoms,chl=chl_anoms,mld=mld_anoms,no3=no3_anoms,nppv=nppv_anoms,nlgz=nlgz_anoms)
  
  sstMat <- extract_env_var(envdat,"sst",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  chlMat <- extract_env_var(envdat,"chl",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  mldMat <- extract_env_var(envdat,"mld",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  no3Mat <- extract_env_var(envdat,"no3",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  nppvMat <- extract_env_var(envdat,"nppv",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  nlgzMat <- extract_env_var(envdat,"nlgz",FeedNames,YrStart,YrEnd,lag,vn,normalize=normalize)
  
  Outs <- list(
    sst = sstMat,
    chl = chlMat,
    mld=mldMat,
    no3=no3Mat,
    nppv=nppvMat,
    nlgz=nlgzMat
  )
  Outs
}

# =================================Build Final Data=========================================
# Convenience functions for building RTMB parameters and constraints
build_base_parameters <- function(Nbreed, MixPars, AddV) {
  list(
    rval = 0.09,
    logK = rep(log(20000), Nbreed),
    logBK = rep(3, Nbreed),
    InfluxP = 10,
    inert_par = 0,
    MixPars = MixPars,
    AddV = AddV
  )
}

build_base_map <- function() {
  list(
    InfluxP = factor(NA),
    inert_par = factor(NA)
  )
}

# Function to make a data list for a specific user-designated scenario
# Parameter definitions and defaults provided in the run script "run_humpback_assessment"
MakeDataScenario <- function(Code, SensCase, StrayBase, DataFileName, Yr1, Yr2,
                             UseKPrior, Kmax,
                             YrSDevs, CatchSer, envOpt, splineK,envlag,envVars,
                             Nage,ByCatchFile, AddCV, MixWeights,
                             MaxN, WithMirror,
                             SF, IAmat, SA, SC, FullDiag,
                             TimeLag, DensDepOpt, WghtTotal, Idirichlet,
                             ... # ignored
) {
  # Specific model variants
  BreedingOpt <- substr(Code, 1, 2)
  FeedingOpt <- substr(Code, 3, 4)
  
  # ==================================Basics====================================================
  # Read in the data file
  DatFile <- read.table(DataFileName, fill = T, col.names = c(1:200), comment.char = "?")

  # Set up breeding and feeding grounds
  which_breed <- which(DatFile[, 1] == "Breeding_grounds" & DatFile[, 2] == BreedingOpt)
  Nbreed <- as.numeric(DatFile[which_breed, 3])
  BreedNames <- DatFile[which_breed + 1, 1:Nbreed]
  which_feed <- which(DatFile[, 1] == "Feeding_grounds" & DatFile[, 2] == FeedingOpt)
  Nfeed <- as.numeric(DatFile[which_feed, 3])
  FeedNames <- DatFile[which_feed + 1, 1:Nfeed]

  Years <- Yr1:(Yr2 + 1)
  Nyear <- length(Years)
  which_year_feedbreed <- which(DatFile[, 2] == "Year_for_feeding_to_breeding")
  YearFeedBreed <- as.numeric(DatFile[which_year_feedbreed + 1, 1]) - Yr1 + 1 # note: gets used for calculate predicted mixing
  which_dirichlet <- which(DatFile[, 2] == "Dirichlet_(1)_or_normal_(0)")
  Idirichlet <- as.numeric(DatFile[which_dirichlet + 1, 1])

  # ==================================Catch and Abundance====================================================
  OutsCatch <- ReadCatches(DatFile, BreedingOpt, FeedingOpt, Nbreed, Nfeed, BreedNames, FeedNames, Yr1, Yr2, CatchSer, ByCatchFile)
  CatchF <- OutsCatch$CatchF
  CatchB <- OutsCatch$CatchB

  # Read the survey data (note adjustments to indices to reflect C++ to R)
  OutsSurv <- ReadSurveyData(Code, BreedingOpt, FeedingOpt, BreedNames, FeedNames, Yr1, Yr2, FullDiag = FullDiag, SensTest = SensCase)
  SurveyI <- OutsSurv$SurveyI
  SurveyR <- OutsSurv$SurveyR
  NsurveyData <- OutsSurv$NsurveyData
  NsurveySeries <- OutsSurv$NsurveySeries
  SurveySeries <- OutsSurv$SurveySeries
  NextraCV1 <- OutsSurv$NextraCV1

  # ==================================Mixing====================================================
  OutsMix <- ReadMixingData(Code, DatFile, Yr1, Yr2, Nbreed, Nfeed, BreedNames, FeedNames, MixWeights, FullDiag = FullDiag, SensTest = SensCase, MaxN = MaxN)
  ObsMixBtoFE <- OutsMix$ObsMixBtoFE
  ObsMixBtoFP <- OutsMix$ObsMixBtoFP
  ObsMixBtoFO <- OutsMix$ObsMixBtoFO
  ObsMixFtoBE <- OutsMix$ObsMixFtoBE
  ObsMixFtoBP <- OutsMix$ObsMixFtoBP
  ObsMixFtoBO <- OutsMix$ObsMixFtoBO
  NmixData <- OutsMix$NmixData
  
  # Find the correct mixing data in the data file
  mix_start <- which(DatFile[, 1] == "Mixing_matrix" & DatFile[, 2] == Code) + 1
  # make the mixing data matrix
  MixI <- map(1:Nbreed, \(x) as.numeric(DatFile[mix_start + x, 1:Nfeed]))
  MixI <- do.call(rbind, MixI)
  rownames(MixI) <- BreedNames
  colnames(MixI) <- FeedNames
  
  # Set up mixing parameters
  NmixPar <- sum(MixI > 0)
  # Mixing parameters
  MixPars <- NULL
  for (Ibreed in 1:Nbreed)
  {
    Iref <- which(MixI[Ibreed,]==-1)
    Temp <- rep(0,Nfeed)
    for (IdataS in 1:2) 
      for (Ifeed in 1:Nfeed) Temp[Ifeed] <- Temp[Ifeed] + (0.001+ObsMixBtoFE[IdataS,Ibreed,Ifeed])/2.0
    for (Ifeed in 1:Nfeed) if (MixI[Ibreed,Ifeed]!=0 && Ifeed!=Iref) MixPars <- c(MixPars,log(Temp[Ifeed]/Temp[Iref]))
  }

  # ==================================Survival and Fecundity====================================================
  
  # Survival
  nFsdevs <- sum(SF)*length(YrSDevs:Yr2) # number of random survival devs
  SFdev <- rep(0, nFsdevs)
  SFsigma <- 1
  log_SFsigma <- log(SFsigma) # variance in SFdevs

  # Carrying capacity deviates
  nKdevs <- length(YrSDevs:Yr2)*sum(SF)
  Kdev <- rep(0, nKdevs)
  Ksigma <- 1
  log_Ksigma <- log(Ksigma)
  
  # Mirror devs between feeding grounds?
  Nmirror <- 0
  Mirror <- matrix(0, nrow = 2, ncol = 2)
  if (WithMirror == 1) {
    Nmirror <- 1
    Mirror[1, 1] <- 4
    Mirror[1, 2] <- 3
    Mirror <- Mirror - 1
  }
  if (WithMirror == 2) {
    Nmirror <- 2
    Mirror[1, 1] <- 4
    Mirror[1, 2] <- 3
    Mirror[2, 1] <- 6
    Mirror[2, 2] <- 5
    Mirror <- Mirror - 1
  }
  
  # ==================================Environment====================================================
  # Specify devs for survival and fecundity
  # Make devs for survival and fecundity
  # The form and inclusion of these will differ based on option "envOpt"
  # list with sst and chl anomalies by year (YrSDevs to Yr2)
  OutsEnv <- ReadGLORYS(DatFile, FeedingOpt, YrStart = YrSDevs, YrEnd = 2023,lag=envlag,type = "anom", normalize = T)
  # OutsEnv <- ReadMOM6(DatFile, FeedingOpt, YrStart = YrSDevs, YrEnd = 2023,lag=envlag,type = "anom", normalize = T)
  
  nomegas <- sum(SF) # number of coefficients to estimate for each environmental covariate
  
  # gather data for the chosen variables
  envParams <- list()
  envData <- list()
  for(i in 1: length(envVars)){
    if(!(envVars[i]%in%names(OutsEnv))) stop("Environmental variable name mismatch. Check envVars.")
    envData[[i]] <- pluck(OutsEnv,envVars[i])
    envParams[[i]] <- rep(0,nomegas)
  }
  names(envData) <- envVars
  envParams <- simplify2array(envParams)

  # Additional random error
  NyrEnv <- length(YrSDevs:Yr2)
  nepsEnv <- sum(SF)*NyrEnv # number of params for extra environmental variation
  epsEnv <- rep(0,nepsEnv)
  sigmaEnv <- 1
  log_sigmaEnv <- log(sigmaEnv)
  
  # polynomial version
  # Function to create polynomial design matrix for multiple covariates
  # make_poly_design <- function(data, poly_specs) {
  #   # poly_specs is a named list: list(var1 = degree1, var2 = degree2, ...)
  #   # describing the polynomial degree for each variable
  #   # Example: list(temp = 2, depth = 3) means temp^0, temp^1, temp^2, depth^0, ..., depth^3
  #   X_list <- list()
  #   for (var_name in names(poly_specs)) {
  #     degree <- poly_specs[[var_name]]
  #     x <- data[[var_name]]
  #     # Create polynomial terms for this variable
  #     for (d in 0:degree) {
  #       col_name <- if (d == 0) "intercept" else paste0(var_name, "_", d)
  #       X_list[[col_name]] <- x^d
  #     }
  #   }
  #   # Combine into matrix (removing duplicate intercepts)
  #   X <- do.call(cbind, X_list)
  #   # Keep only one intercept column
  #   intercept_cols <- which(colnames(X) == "intercept")
  #   if (length(intercept_cols) > 1) {
  #     X <- X[, -intercept_cols[-1], drop = FALSE]
  #   }
  #   return(X)
  # }
  # Xpolyl <- map(1:length(FeedNames),\(x){
  #   envl <- map(envData,\(x) x[i,])
  #   specs <- rep(2,length(envVars)) # degree of polynomial for each variable
  #   names(specs) <- envVars
  #   make_poly_design(envl,specs)
  # }) |> simplify2array()
  # Xpolyl <- Xpolyl[,,which(SF==1)]
  # ncoef <- dim(Xpolyl)[2]
  # betaPoly <- matrix(0,nrow=sum(SF),ncol=ncoef) # rows are feeding grounds, columns are components of the polynomial
  # log_sigmaPoly <- 0
  
  ## Splines version (thin-plate regression spline)
  # ncoef <- splineK-1
  # X_list<- list()
  # S_list <- list()
  # for(i in 1:length(envVars)){
  #   Xtemp <- array(0,dim=c(NyrEnv,ncoef,sum(SF)))
  #   Stemp <- array(0,dim=c(ncoef,ncoef,sum(SF)))
  #   for(j in 1:sum(SF)){
  #     x <- data.frame(env=envData[[i]][j,])
  #     spline_setup_temp <- smoothCon(
  #       s(env, k=splineK,bs='cr'), # 8 knots
  #       data = x,
  #       absorb.cons = TRUE)[[1]]
  #     Xtemp[,,j] <- spline_setup_temp$X
  #     Stemp[,,j] <- spline_setup_temp$S[[1]]
  #   }
  #   X_list[[i]] <- Xtemp
  #   S_list[[i]] <- Stemp
  # }
  # Xarr <- simplify2array(X_list)
  # Sarr <- simplify2array(S_list)
  # # # Coefficients for splines (k=3 means 2 coefficients for each spline)
  # beta_splines <- array(0,dim=c(sum(SF),ncoef,length(envVars)))
  # # # SD of splines
  # sigma_splines <- matrix(1,nrow=sum(SF),ncol=length(envVars))
  # lambda_splines <- matrix(1,nrow=sum(SF),ncol=length(envVars))
  
  # ==================================Density Dependence====================================================
  # Strength of density-dependence in survival and fecundity
  alphaK <- 1
  log_alphaK <- log(alphaK)
  betaK <- 0.3
  log_betaK <- log(betaK)
  # logBKsigma <- log(4)
  
  # Set up additional variances
  NextraCV <- ifelse(NextraCV1 == 0, 1, NextraCV1)
  AddV <- rep(0, NextraCV)

  # ==================================Parameter List====================================================
  # ---- Step 1: Build common parameter base ----
  # These parameters are used by ALL model variants
  parameters <- list(
    rval = 0.09,                           # Maximum growth rate
    logK = rep(log(20000), Nbreed),        # Breeding ground K
    logBK = rep(3, Nbreed),                # Relative depletion of each breeding stock
    InfluxP = 10,                          # Transfer influx
    inert_par = 0,                         # Related to fecundity calculation
    MixPars = MixPars,                     # Mixing matrix parameters
    AddV = AddV                            # Additional variance parameters
  )
  
  # ---- Step 2: Add model-specific parameters ----
  
  if (envOpt == "rS") {
    # Recruitment-Survival model: estimate survival deviations
    parameters <- build_rS_parameters(parameters, SFdev, log_SFsigma)
    
  } else if (envOpt == "env-survival") {
    # Environment drives survival directly
    parameters <- build_env_survival_parameters(parameters, epsEnv, log_sigmaEnv, envParams)
    
  } else if (envOpt == "ddOnly") {
    # Density-dependence only (no environment)
    parameters <- build_ddOnly_parameters(parameters, SFdev, Kdev, log_Ksigma,
                                          log_alphaK, log_betaK, envParams)
    
  } else if (envOpt == "env-K") {
    # Environment drives carrying capacity
    parameters <- build_env_K_parameters(parameters, SFdev, Kdev, log_Ksigma,
                                         log_alphaK, log_betaK, envParams)
  } else if (envOpt == "envS-rvec") {
    parameters <- build_envS_rvec_parameters(parameters, epsEnv, log_sigmaEnv, envParams, Nfeed)
  } else {
    stop("Unknown model variant: ", envOpt, 
         ". Must be one of: rS, env-survival, ddOnly, env-K")
  }
  
  # ---- Step 3: Build map of RTMB constraints ----
  # Parameters that are fixed for ALL models
  mapUse <- list(
    InfluxP = factor(NA),
    inert_par = factor(NA)
  )
  
  # ---- Step 4: Add model-specific map constraints ----
  # These determine which parameters are NOT estimated
  
  if (envOpt == "rS") {
    mapUse <- build_rS_map(mapUse, AddCV = AddCV)
    
  } else if (envOpt == "env-survival") {
    mapUse <- build_env_survival_map(mapUse, AddCV = AddCV)
    
  } else if (envOpt == "ddOnly") {
    mapUse <- build_ddOnly_map(mapUse, SFdev, Kdev, envParams, AddCV = AddCV)
    
  } else if (envOpt == "env-K") {
    mapUse <- build_env_K_map(mapUse, SFdev, envParams, AddCV = AddCV)
  }
  # ==================================Output====================================================
  # Return a big list of data for the model
  datout <- list(
    Nbreed = Nbreed,
    Nfeed = Nfeed,
    BreedNames = BreedNames,
    FeedNames = FeedNames,
    UseKPrior = UseKPrior,
    WithMirror = WithMirror,
    AddCV=AddCV,
    MixWeights=MixWeights,
    Kmax = Kmax,
    Yr1 = Yr1,
    Yr2 = Yr2,
    YrSDevs = YrSDevs,
    Years = Years,
    IAmat = IAmat,
    SA = SA,
    SC = SC,
    TimeLag = TimeLag,
    SurveyI = SurveyI,
    SurveyR = SurveyR,
    NextraCV1 = NextraCV1,
    CatchB = CatchB,
    CatchF = CatchF,
    # Xpolyl = Xpolyl, # polynomial design matrix for env. vars
    # Xarr = Xarr, # design matrices, environmental splines; NyrxNcoefxFeedxvariable
    # Sarr = Sarr, # spline penalty matrices
    SF = SF,
    MixI = MixI,
    NsurveyData = NsurveyData,
    NsurveySeries = NsurveySeries,
    SurveySeries = SurveySeries,
    NmixData = NmixData,
    NmixPar = NmixPar,
    YearFeedBreed = YearFeedBreed,
    Idirichlet = Idirichlet,
    ObsMixBtoFE = ObsMixBtoFE,
    ObsMixBtoFP = ObsMixBtoFP,
    ObsMixBtoFO = ObsMixBtoFO,
    ObsMixFtoBE = ObsMixFtoBE,
    ObsMixFtoBP = ObsMixFtoBP,
    ObsMixFtoBO = ObsMixFtoBO,
    StrayBase = StrayBase,
    DensDepOpt = DensDepOpt,
    Nmirror = Nmirror,
    Mirror = Mirror,
    WghtTotal = WghtTotal,
    Nage = Nage,
    map = mapUse,
    envVars=envVars,
    envData = envData,
    parameters = parameters
  )
  
  return(datout)
}
