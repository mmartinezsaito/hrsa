preprofnc_hrsaListener_FrankeDegen16 = function(datafile, datcolnames) {

    # Check if data file exists
    if (!file.exists(datafile)) {
      stop("** Data file does not exist. Please check again. **\n")
    }
    # Load the data
    datable <- data.table::fread(file = datafile, header = TRUE, sep = "\t", data.table = TRUE,
                                 fill = TRUE, stringsAsFactors = TRUE, logical01 = FALSE)
    #separator is "\t" because fread() has trouble reading space delimited files with missing values.
    
    # Check if necessary data columns all exist (while ignoring case and underscores)
    if (!all(datcolnames %in% colnames(datable))) {
      stop("** Data file is missing one or more necessary data columns. Please check again. **\n",
           "  Required data columns are: \"", paste0(datcolnames, collapse = "\", \""), "\".\n")
    }
    
    # Remove only the rows containing NAs in required columns
    compl_rows       <- complete.cases(datable[, datcolnames, with = FALSE])
    sum_incompl_rows <- sum(!compl_rows)
    if (sum_incompl_rows > 0) {
      datable <- datable[compl_rows, ]
      cat("\n")
      cat("The following lines of the data file have NAs in required columns:\n")
      cat(paste0(head(which(!compl_rows), 100) + 1, collapse = ", "))
      if (sum_incompl_rows > 100) {
        cat(", ...")
      }
      cat(" (total", sum_incompl_rows, "lines)\n")
      cat("These rows are removed prior to modeling the data.\n")
    }
    
    # general preprocessing variables	
    sids       <- NULL # List of unique subjects (1D)
    nsid       <- NULL # Total number of subjects (0D)
    nblk_p_sid <- NULL # Number of blocks per each subject (1D)
    nblk_max   <- NULL # Maximum number of blocks across all subjects (0D)
    nt_p_sid   <- NULL # Number of trials (per block) per subject (2D or 1D)
    nt_max     <- NULL # Maximum number of trials across all blocks & subjects (0D)
    
    # To avoid NOTEs by R CMD check
    .N <- NULL
    sid <- NULL
    
    dt_ntsid <- datable[, .N, by = "sid"]
    sids     <- dt_ntsid$sid
    nsid     <- length(sids)
    nt_p_sid <- dt_ntsid$N
    nt_max   <- max(nt_p_sid)

    # Task-specific preprocessing
    reft    <- array(0, c(nsid, nt_max))
    refc    <- array(0, c(nsid, nt_max))
    refd    <- array(0, c(nsid, nt_max))
    #trigger <- array(0, c(nsid, nt_max))
    choice  <- array(0, c(nsid, nt_max))
    cond    <- array(0, c(nsid, nt_max))
 
    for (i in 1:nsid) {
      s <- sids[i]
      t <- nt_p_sid[i]
      dt <- datable[sid == s]
      
      reft[i, 1:t]    <- dt$reft
      refc[i, 1:t]    <- dt$refc
      refd[i, 1:t]    <- dt$refd
      #trigger[i, 1:t] <- dt$trigger
      choice[i, 1:t]  <- dt$choice
      cond[i, 1:t]    <- dt$cond
    }
    
    if (Sys.info()["sysname"] == "Windows") {
      if (Sys.info()["nodename"] == "AS-98") datdir <- "//tsclient/KrivokolennyServer/rsa/"
      else                                   datdir <- "C:/Users/NUS/Mario/Dropbox/NUS_YuRj/rsa/"
    } else if (Sys.info()["nodename"]=="Hydra") {
      datdir <- "~/Dropbox/NUS_YuRj/rsa/"
    } else if (Sys.info()["nodename"]=="SootyTern") {
      datdir <- "~/Dropbox/NUS_YuRj/rsa/"
    } 
    sal <- read.table(paste0(datdir, "Sb7.txt"))
    sal <- as.matrix(sal)
    
    datli <- list(
      N       = nsid,
      T       = nt_max,
      Ts      = nt_p_sid,
      reft    = reft,
      refc    = refc,
      refd    = refd,
      #trigger = trigger,
      choice  = choice,
      cond    = cond,
      sal     = sal
    )

  }
