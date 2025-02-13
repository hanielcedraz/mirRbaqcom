#!/usr/bin/env Rscript

# ########################################
# ### LOADING PACKAGES
# ########################################
if (suppressPackageStartupMessages(!require(pacman))) suppressPackageStartupMessages(install.packages("pacman"))

# if (suppressPackageStartupMessages(!require(remotes))) suppressPackageStartupMessages(install.packages("remotes"))

# if (suppressPackageStartupMessages(!require(taxize))) suppressPackageStartupMessages(remotes::install_github("ropensci/taxize"))
# suppressMessages(library(taxize))

suppressPackageStartupMessages(
    p_load(
        tools, 
        parallel, 
        optparse, 
        dplyr, 
        data.table, 
        glue, 
        stringr,
        tidyr,
        readr,
        tibble,
        # reticulate,
        mcprogress,
        progress,
        lubridate,
        # rlist,
        cli
    )
)


#########################################
### LOADING FUNCTIONS
#########################################
# Function to stop the process without message
stopQuietly <- function(Text = "Stoping the process") {
    message <- sprintf("\r%s\r", Text);
    stop(simpleError(message));
}


## Functino to check animal name given in --specie option
# checkAnimalName <- function(name = "cow") {
#     
#     animalName <- paste0(toupper(substring(name, 1, 1)), tolower(substring(name, 2)))
#     
#     animalCheck <- read.csv("https://copylists.com/downloads/animals/animals/animals.csv", header = FALSE) %>% filter(V2 == animalName) %>% pull(V2)
#     
#     if (length(animalCheck) != 1) {
#         return(NA)
#     } else {
#         return(animalCheck)
#     }
#     
# }

## Function to deal with the number of threads
prepareCore <- function(nThreads = 8){
    # if opt_procs set to 0 then expand to samples by targets
    
    if (detectCores() < nThreads) {
        write(glue("The number of cores specified ({nThreads}) is greater than the number of cores available ({detectCores()})"), stdout())
        glue("Using {detectCores()} threads")
        nThreads <- detectCores()
    }
    
    
    #if (detectCores() < opt$procs) nThreads <- detectCores()
    cat("\n\n");cli_h1("Preparing processors");cli_alert_success(" Using {nThreads} processors")
    return(nThreads)
}


# Function to get species code
# Function to get species code using taxize
# get_species_code <- function(common_name = "cow") {
#     # Get the scientific name using taxize
#     Sys.setenv("ENTREZ_KEY" = "071b0f34fa2c951d7bedeeeeeadc23a46e0a")
#     scientific_name <- taxize::comm2sci(common_name, db = "ncbi")
#     
#     spc <- str_split(scientific_name[[1]], pattern = " ")
#     
#     
#     abbName <- sapply(spc, function(x) {
#         first_letter <- substr(x[1], 1, 1)
#         first_two_letters <- substr(x[2], 1, 2)
#         tolower(paste0(first_letter, first_two_letters))
#     })
#     
#     if (length(abbName) == 0) {
#         return(NA)
#     }
#     
#     
#     return(abbName)
# }


### LOGING
# https://cran.r-project.org/web/packages/logging/index.html
# https://cran.r-project.org/web/packages/log4r/index.html

# Function to load user input
userInput <- function(question) {
    cat(question)
    con <- file("stdin")
    on.exit(close(con))
    n <- readLines(con, n = 1)
    return(n)
}


# Querying fastq files for paralelization
createFastqQuery <- function(sampleFile = NULL, inputFolder, outputSufix) {
    fastqList <- list()
    # reads <- dir(path = file.path(inputFolder), pattern = "fastq")
    
    
    if (is.null(sampleFile)) {
        sampleTable <- data.frame(
            Read_1 = list.files(inputFolder, pattern = "R1|SE")
        ) 
        
        sampleTable <- sampleTable %>% 
            separate(Read_1, into = "SAMPLE_ID", sep = "_", remove = FALSE, extra = "drop") %>% 
            relocate(SAMPLE_ID, .before = Read_1) %>% 
            mutate(basenameOutput = str_remove_all(basename(Read_1), ".fastq|.gz")) %>% 
            as.data.frame()
        
    } else {
        sampleTable <- read.table(sampleFile, header = TRUE) %>% 
            relocate(SAMPLE_ID, .before = Read_1) %>% 
            mutate(basenameOutput = str_remove_all(basename(Read_1), ".fastq|.gz")) %>% 
            as.data.frame()
    }
    
    
    # i = 1
    for (i in 1:nrow(sampleTable)) {
        map <- lapply(c("_R1_"), grep, x = sampleTable$Read_1, value = TRUE)
        names(map) <- c("R1")
        map$R1 <- sampleTable[i,"Read_1"]
        map$sampleName <-  glue("{sampleTable[i,'basenameOutput']}")
        map$output <- glue("{sampleTable[i,'basenameOutput']}_{outputSufix}")
        fastqList[[paste(map$sampleName)]] <- map
    }
    
    
    cat("\n\n");cli_h1("Preparing sample query");cli_alert_info(" Setting up {length(fastqList)} jobs");cat("\n\n")
    return(fastqList)
}

checkFastqFolder <- function(folder) {
    # Check if the directory exists
    cli_h1("Checking input folder - {folder}");
    if (!dir.exists(folder)) {
        cli_abort(c("x" = "Something went wrong! Directory {folder} does not exist!"));cat("\n\n")
    }
    
    # List files in the directory
    fastqFiles <- list.files(folder, pattern = "fastq")
    
    if (length(fastqFiles) == 0) {
        cli_abort(c("x" = "Something went wrong! Directory {folder} does not have any fastq file!"));cat("\n\n")
    } else {
        cli_alert_success(" Directory {folder} is ready");cat("\n\n")
    }
    
    
    if (!dir.exists("reports")){
        dir.create("reports", recursive = TRUE)
    }
    
    
}


########################################
### SETING PARAMETERS
########################################
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")

option_list <- list(
    make_option(
        opt_str = c("-f", "--file"),
        type = "character",
        default = NULL,
        help = "File containing sample names and paths. [default %default]",
        dest = "samplesFile"
    ),
    # make_option(
    #     opt_str = c("-m", "--multiqc"), 
    #     action = 'store_true', 
    #     type = "logical",
    #     default = FALSE,
    #     help = "Use this option if you want to run multiqc software  [default %default]",
    #     dest = "multiqc"
    # ),
    make_option(
        opt_str = c("-d", "--rawFastq"),
        #action = 'store_true',
        type = "character",
        default = "00-Fastq",
        help = "The directory raw fastq files are in  [default %default]",
        dest = "rawFastq"
    ),
    make_option(
        opt_str = c("-o", "--cleanedReads"),
        #action = 'store_true',
        type = "character",
        default = "01-CleanedReads",
        help = "The directory to store cleaned fastq files  [default %default]",
        dest = "cleanedReads"
    ),
    make_option(
        opt_str = c("-a", "--adapter"),
        #action = 'store_true',
        type = "character",
        default = "AAAAAAAAAA",
        help = "Sequence of an adapter ligated to the 3' end (paired data: of the first read). The adapter and subsequent bases are trimmed. If a '$' character is appended('anchoring'), the adapter is only found if it is a suffix of the read. [default %default]",
        dest = "adapter"
    ),
    make_option(
        opt_str = c("-p", "--processors"), 
        type = "integer", 
        default = 8,
        help = "Number of processors to use [default %default]",
        dest = "procs"
    ),
    make_option(
        opt_str = c("-q", "--sampleprocs"), 
        type = "integer", 
        default = 2,
        help = "Number of samples to process in parallel [default %default]",
        dest = "sampleToprocs"
    ),
    make_option(
        opt_str = c("-m", "--minimumLength"),
        type = "integer", 
        default = 15,
        help = "Discard reads shorter than <value> [default %default].",
        dest = "minimumLength"
    ),
    make_option(
        opt_str = c("-M", "--maximumLength"),
        type = "integer",
        default = NULL,
        help = "Discard reads longer than <value> [default %default].",
        dest = "maximumLength"
    ),
    make_option(
        opt_str = c("-u", "--cut"),
        type = "integer", 
        default = 3,
        help = "Remove LEN bases from each read (or R1 if paired; use -U option for R2). If LEN is positive, remove bases from the beginning. If LEN is negative, remove bases from the end. Can be used twice if LENs have different signs. Applied *before* adapter trimming. [ default %default]",
        dest = "cut"
    ),
    make_option(
        opt_str = c("-N", "--maxN"),
        type  = 'integer',
        default = 0,
        help = "Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length. [ default %default]",
        dest = "maxN"
    ),
    make_option(
        opt_str = c("-t", "--runTrimmomatic"),
        action  =  'store_true',
        #type  =  "logic",
        default = FALSE,
        help = "Run Trimmomatic to perform another cleaning - Check the data afterwords [default %default]",
        dest = "runTrimmomatic"
    )
    # make_option(
    #     opt_str = c("-e", "--condaEnv"),
    #     type = "character",
    #     default = "genderPurity",
    #     help = "Name of the conda env [default %default]. If new env, please provide the file containing the conda env to be created",
    #     dest = "condaEnv"
    # )
    # make_option(
    #     opt_str = c("-w", "--pmode"), 
    #     action = "store_true", 
    #     default = FALSE,
    #     help  =  "Use this option if you want to run two pass mode mapping. --mappingProgram must be hisat2 [default %default].",
    #     dest  =  "PassMode"
    # ),
    # make_option(
    #     opt_str = c("-o", "--indexFiles"), 
    #     type  = 'character', 
    #     default = 'ht2_base',
    #     help = "The basename of the index files to write. --mappingProgram must be hisat2 [%default].",
    #     dest = "indexFiles"
    #     ),
    # make_option(
    #     opt_str = c("-n", "--notification"),
    #     type = "logical",
    #     action = "store_false",
    #     default = TRUE,
    #     help = "Wheter to get notification from email. [default %default]",
    #     dest = "notification"
    # ),
    
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list, description =  paste('Authors: OLIVEIRA, H.C.; IBELLI, A.M.G.', 'Version: 0.1.0', 'E-mail: haniel.cedraz@stgen.com', sep = "\n", collapse = '\n')))


#################################################################################

cat("\n\n");cli_rule(center = "Quality Control for miRNA seq analysis");cat("\n\n")



## Check if input folder exist and it is not empty
checkFastqFolder(opt$rawFastq)


## Create folder to store the cleaned reads
cli_h1("Checking output folder - {basename(opt$cleanedReads)}")
if (!dir.exists(opt$cleanedReads)) {
    
    dir.create(opt$cleanedReads, recursive = TRUE)
    
    if (dir.exists(opt$cleanedReads)) {
        cli_alert_success(" Directory {opt$cleanedReads} created successfuly!");cat("\n\n")
        
    } else {
        cli_abort(c("x" = "Directory {opt$cleanedReads} could not be created!"));cat("\n\n")
    }
} else {
    cli_alert_success(" Directory {opt$cleanedReads} is ready");cat("\n\n")
}

# opt$samplesFile <- "samples.txt"
cutadaptQuery <- createFastqQuery(sampleFile = opt$samplesFile, inputFolder = opt$rawFastq, outputSufix = "cleaned_cutadapt")

procs <- prepareCore(opt$procs);cat("\n\n")


runCutAdapter <- mcprogress::pmclapply(cutadaptQuery, function(index) {
    try({
        system(
            glue(
                "cutadapt  --quiet -m {opt$minimumLength} -u {opt$cut}",
                " -a {opt$adapter} -j {procs}",
                "--minimum-length {opt$minimumLength}",
                ifelse(!is.null(opt$maximumLength), "--maximum-length {opt$maximumLength}", ""),
                "--max-n {opt$maxN}",
                "{opt$rawFastq}/{index$R1} -o {opt$cleanedReads}/{index$output}.fastq",
                "> reports/{index$sampleName}_cutadapt_reports.log",
                "2> reports/{index$sampleName}_cutadapt_reports.err",
                .sep = " "
            )
        )
    })
}, spinner = TRUE, eta = TRUE, title = "Running cutadapt, please wait...", mc.cores = opt$sampleToprocs)

#

if (opt$runTrimmomatic) {
    trimmomaticQuery <- createFastqQuery(sampleFile = NULL, inputFolder = opt$cleanedReads, outputSufix = "trimmomatic")
    
    runTrimmomatic <- mcprogress::pmclapply(trimmomaticQuery, function(index){
        #         write(glue("trimmomatic SE -threads {procs} ",
        #                    "{opt$cleanedReads}/{index$R1}",
        #                    "{opt$cleanedReads}/{index$output}.fastq.gz",
        #                    "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3
        # MINLEN:18",
        #                    "> reports/{index$sampleName}_trimmomatic_reports.log",
        #                    "2> reports/{index$sampleName}_trimmomatic_reports.err",
        #                    .sep = " "
        #         ), stdout())
        try({
            system(
                glue(
                    "trimmomatic SE -quiet -threads {procs} ",
                    "{opt$cleanedReads}/{index$R1}",
                    "{opt$cleanedReads}/{index$output}.fastq.gz",
                    "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:18",
                    "-summary reports/{index$sampleName}_trimmomatic_reports.log",
                    .sep = " "
                )
            )
        })
    }, spinner = TRUE, eta = TRUE, title = "Running trimmomatic, please wait...", mc.cores = opt$sampleToprocs)
    
}






