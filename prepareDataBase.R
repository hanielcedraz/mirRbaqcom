#!/usr/bin/env Rscript

# ########################################
# ### LOADING PACKAGES
# ########################################
if (suppressPackageStartupMessages(!require(pacman))) suppressPackageStartupMessages(install.packages("pacman"))

if (suppressPackageStartupMessages(!require(remotes))) suppressPackageStartupMessages(install.packages("remotes"))

if (suppressPackageStartupMessages(!require(taxize))) suppressPackageStartupMessages(remotes::install_github("ropensci/taxize"))
suppressMessages(library(taxize))
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
        reticulate,
        pbmcapply,
        progress,
        lubridate,
        rlist,
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
checkAnimalName <- function(name = "cow") {
    
    animalName <- paste0(toupper(substring(name, 1, 1)), tolower(substring(name, 2)))
    
    animalCheck <- read.csv("https://copylists.com/downloads/animals/animals/animals.csv", header = FALSE) %>% filter(V2 == animalName) %>% pull(V2)
    
    if (length(animalCheck) != 1) {
        return(NA)
    } else {
        return(animalCheck)
    }
    
}

## Function to deal with the number of threads
prepareCore <- function(nThreads = 8){
    # if opt_procs set to 0 then expand to samples by targets
    
    if (detectCores() < nThreads) {
        write(glue("The number of cores specified ({nThreads}) is greater than the number of cores available ({detectCores()})"), stdout())
        glue("Using {detectCores()} threads")
        nThreads <- detectCores()
    }
    
    
    #if (detectCores() < opt$procs) nThreads <- detectCores()
    cli_alert_success(" Using {nThreads} processors")
    return(nThreads)
}


# Function to get species code
# Function to get species code using taxize
get_species_code <- function(common_name = "cow") {
    # Get the scientific name using taxize
    Sys.setenv("ENTREZ_KEY" = "071b0f34fa2c951d7bedeeeeeadc23a46e0a")
    scientific_name <- taxize::comm2sci(common_name, db = "ncbi")
    
    spc <- str_split(scientific_name[[1]], pattern = " ")
    
    
    abbName <- sapply(spc, function(x) {
        first_letter <- substr(x[1], 1, 1)
        first_two_letters <- substr(x[2], 1, 2)
        tolower(paste0(first_letter, first_two_letters))
    })
    
    if (length(abbName) == 0) {
        return(NA)
    }
    
    
    return(abbName)
}


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



########################################
### SETING PARAMETERS
########################################
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")

option_list <- list(
    # make_option(
    #     opt_str = c("-f", "--file"), 
    #     type = "character", 
    #     default = "S1_data/gvi_skimseq_RunDate.txt",
    #     help = "The filename of the sample file. The sample file should be stored at S1_data/ as a base name gvi_skimseq_ followed by the run date(e.g. gvi_skimseq_220909) [default %default]",
    #     dest = "samplesFile"
    # ),
    # make_option(
    #     opt_str = c("-m", "--multiqc"), 
    #     action = 'store_true', 
    #     type = "logical",
    #     default = FALSE,
    #     help = "Use this option if you want to run multiqc software  [default %default]",
    #     dest = "multiqc"
    # ),
    make_option(
        opt_str = c("-d", "--database"),
        #action = 'store_true',
        type = "character",
        default = glue("database"),
        help = "The directory to store all reference files  [default %default]",
        dest = "database"
    ),
    make_option(
        opt_str = c("-s", "--specie"),
        #action = 'store_true',
        type = "character",
        default = "Cow",
        help = "Specie name - Common name (Cow) or Scientific name (Bos_taurus) [default %default]",
        dest = "specie"
    ),
    make_option(
        opt_str = c("-p", "--processors"), 
        type = "integer", 
        default = 8,
        help = "Number of processors to use [default %default]",
        dest = "procs"
    ),
    # make_option(
    #     opt_str = c("-q", "--sampleprocs"), 
    #     type = "integer", 
    #     default = 2,
    #     help = "Number of samples to process at each time [default %default]",
    #     dest = "sampleToprocs"
    # ),
    make_option(
        opt_str = c("-t", "--referenceGenome"), 
        type = "character", 
        default = "ARS-UCD1.2_Btau5.0.1Y.fa",
        help = "Path to the fasta file [reference genome] (default %default).",
        dest = "referenceGenome"
    ),
    make_option(
        opt_str = c("-M", "--useMirgeneDB"), 
        action = "store_true", 
        #type = "character", 
        default = FALSE,
        help = "Use this option for merging mature file from mirGeneDB and miRBase [%default]",
        dest = "useMirgeneDB"
    )
    # make_option(
    #     opt_str = c("-z", "--libraryType"),
    #     type  = 'character', 
    #     default = "singleEnd",
    #     help = "The library type to use. Available: 'pairEnd' or 'singleEnd'. [ default %default]",
    #     dest = "libraryType"
    # ),
    # make_option(
    #     opt_str = c("-M", "--program"),
    #     type  = 'character', 
    #     default = "bowtie",
    #     help = "Which mapping program to use. Options: 'bowtie', 'bwa'. [ default %default]",
    #     dest = "mappingProgram"
    # ),
    # make_option(
    #     opt_str = c("-x", "--useConda"),
    #     action  =  'store_true',
    #     #type  =  "logic",
    #     default = FALSE,
    #     help = "Run the pipeline inside a conda environment [default %default]",
    #     dest = "useConda"
    # ),
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

cat("\n\n");cli_rule(center = "Preparing database for miRNA seq analysis");cat("\n\n")

## Create folder to store the reference files
if (!dir.exists(opt$database)) {
    dir.create(opt$database, recursive = TRUE)
}

## Getting the spece code for filtering the mature and hairpin files
# opt$specie <- "chicken"

checkSpecieName <- checkAnimalName(opt$specie)
if (!is.na(checkSpecieName)) {
    spcCode <- get_species_code(checkSpecieName)
    if (is.na(spcCode)) {
        cli_abort(c(" x" = "'{opt$specie}' not found, please use the Specie genus name in --specie option. e.g Bos_taurus"))
    }
    
    cat("\n");cli_alert_success(" Using {spcCode}");cat("\n")
} else {
    
    cat("\n\n");cli_alert_info(" Creating code for scientific name '{opt$specie}'.")
    spc <- str_split(opt$specie, pattern = "_")
    
    if (length(spc[[1]]) == 2) {
        spcCode <- sapply(spc, function(x) {
            first_letter <- substr(x[1], 1, 1)
            first_two_letters <- substr(x[2], 1, 2)
            tolower(paste0(first_letter, first_two_letters))
        })
        
        cli_alert_success(" Using {spcCode}")
    } else {
        cli_abort(c(" x" = "'{opt$specie}' not found, please use a valid name for  Specie genus in --specie option. e.g Bos_taurus"))
    }
    
}


#################################################################################

## Prepare reference files
## Downloade reference sequences from mirBase
## Download the files mature.fa (https://www.mirbase.org/download/mature.fa) and hairpin.fa (https://www.mirbase.org/download/hairpin.fa)

cat("\n\n");cli_rule(center = "Downloading mature miRNA file from miRBase");cat("\n\n")
matureFile <- glue("{opt$database}/mature.fa")
if (!file.exists(matureFile)) {
    cli_process_start(" Downloading mature file from miRBase (https://www.mirbase.org/download/)")
    invisible(download.file("https://www.mirbase.org/download/mature.fa", destfile = matureFile, extra = "--verbose"))
    
    if (!file.exists(matureFile)) {
        cli_process_failed()
        cli_abort(c(" x" = "Something went wrong downloading the file. Please download it manually and save it in the database folder:\n wget https://www.mirbase.org/download/mature.fa -O {matureFile}"))
    } else {
        cli_process_done()
    }
}

## Filtering specie sequences from mature.fa
cat("\n\n");cli_h1("Preparing mature file");cat("\n\n")
specieMatureFile <- glue("{opt$database}/{spcCode}_mature.fa")
if (!file.exists(specieMatureFile)) {
    command <- glue("grep -A 1 '{spcCode}' {matureFile} > {specieMatureFile}")
    system(command)
}

## Sequences from mirBase have U instead of T, so it needs to be converted base U to T. Uisng mirDEEP2 suite conda install bioconda::mirdeep2
matureDNAFile <- glue("{opt$database}/{spcCode}_mature_dna.fa")
if (!file.exists(matureDNAFile)) {
    command <- glue("rna2dna.pl {specieMatureFile} > {matureDNAFile}")
    system(command)
}


## Remove characters that are not nucleotydes from the mature file
finalMatureFile <- glue("{opt$database}/{spcCode}_mature_dna_final.fa")
if (!file.exists(finalMatureFile)) {
    command <- glue("fastaparse.pl {matureDNAFile} > {finalMatureFile}")
    system(command)
}

## Check if all non nucleotydes characters were removed
fasta <- readLines(finalMatureFile) 

if (any(stringr::str_detect(fasta, "--"))) {
    newFasta <- gsub("--", "", fasta)
    if (!any(stringr::str_detect(newFasta, "--"))) {
        writeLines(newFasta, con = glue("{finalMatureFile}"))
    }
}

# opt$useMirgeneDB <- "bta_mirGeneDB_mature.fas"
if (opt$useMirgeneDB) {
    
    # Download mature file from mirGeneDB using curl
    fileLink <- glue("https://www.mirgenedb.org/fasta/{spcCode}?mat=1")
    cli_process_start(" Downloading mature file from mirGeneDB: {fileLink}")
    
    mirGeneDBFile <- glue("{opt$database}/{spcCode}_mirGeneDB_mature.fas")
    command <- glue("curl -s {fileLink}  > {mirGeneDBFile}")
    system(command)
    
    if (!file.exists(mirGeneDBFile)) {
        cli_process_failed()
        cli_abort(c(" x" = "Something wrong downloading file from {fileLink}"))
    } 
    
    cli_process_done()
    
    
    ## convert U to T first
    mirGeneFileDNA <- glue("{opt$database}/{spcCode}_mirGeneDB_mature_dna.fas")
    command <- glue("rna2dna.pl {mirGeneDBFile} > {mirGeneFileDNA}")
    system(command)
    
    
    ## merge mature fasta from mirBase with mature fasta from mirGeneDB
    finalMergedFile <- glue("{opt$database}/{spcCode}_mature_mirbase_mirgeneDB_dna.fa")
    command <- glue("cat {mirGeneFileDNA} {finalMatureFile} > {finalMergedFile}")
    system(command)

}



## Build index from mature files
if (opt$useMirgeneDB) {
    indexFolder <- glue("{opt$database}/bowtie_{spcCode}_ref_mature_mirbase_mirgeneDB")
    indexMatureFiles <- glue("{spcCode}_mature_mirbase_mirgeneDB_dna.fa")
    refMatureFile <- finalMergedFile
} else {
    indexFolder <- glue("{opt$database}/bowtie_{spcCode}_ref_mature")
    indexMatureFiles <- glue("{spcCode}_mature_dna_final")
    refMatureFile <- finalMatureFile
}

if (!dir.exists(indexFolder)) {
    dir.create(indexFolder, recursive = TRUE)
}

cli_process_start(" Indexing mature files using bowtie")
command <- glue("bowtie-build -q {refMatureFile} {indexFolder}/{indexMatureFiles}")
system(command)

bowtieFiles <- list.files(indexFolder, pattern = indexMatureFiles)

if (length(bowtieFiles) > 1){
    cli_process_done()
} else {
    cli_process_failed()
}







cli_alert_success(" Preparing mature file... done")


#################################################################################


cat("\n\n");cli_rule(center = "Downloading hairpin file from miRBase");cat("\n\n")
hairpinFile <- glue("{opt$database}/hairpin.fa")
if (!file.exists(hairpinFile)) {
    cli_process_start(" Downloading hairpin file from miRBase (https://www.mirbase.org/download/)")
    invisible(download.file("https://www.mirbase.org/download/hairpin.fa", destfile = hairpinFile, extra = "--verbose"))
    
    if (!file.exists(hairpinFile)) {
        cli_process_failed()
        cli_abort(c(" x" = "Something went wrong downloading the file. Please download it manually and save it in the database folder:\n wget https://www.mirbase.org/download/hairpin.fa -O {hairpinFile}"))
    } else {
        cli_process_done()
    }
}

## Filtering specie sequences from hairpin.fa
cat("\n\n");cli_h1("Preparing hairpin file");cat("\n\n")
specieHairpinFile <- glue("{opt$database}/{spcCode}_hairpin.fa")
if (!file.exists(specieHairpinFile)) {
    command <- glue("grep -A 1 '{spcCode}' {hairpinFile} > {specieHairpinFile}")
    system(command)
}

## Sequences from mirBase have U instead of T, so it needs to be converted base U to T. Uisng mirDEEP2 suite conda install bioconda::mirdeep2
hairpinDNAFile <- glue("{opt$database}/{spcCode}_hairpin_dna.fa")
if (!file.exists(hairpinDNAFile)) {
    command <- glue("rna2dna.pl {specieHairpinFile} > {hairpinDNAFile}")
    system(command)
}


## Remove characters that are not nucleotydes from the hairpin file
finalhairpinFile <- glue("{opt$database}/{spcCode}_hairpin_dna_final.fa")
if (!file.exists(finalhairpinFile)) {
    command <- glue("fastaparse.pl {hairpinDNAFile} > {finalhairpinFile}")
    system(command)
}


## Check if all non nucleotydes characters were removed
fasta <- readLines(finalhairpinFile) 

if (any(stringr::str_detect(fasta, "--"))) {
    newFasta <- gsub("--", "", fasta)
    if (!any(stringr::str_detect(newFasta, "--"))) {
        writeLines(newFasta, con = glue("{finalhairpinFile}"))
    }
}



## build hairpin index reference genome using bowtie
indexFolder <- glue("{opt$database}/bowtie_{spcCode}_ref_hairpin")
if(!dir.exists(indexFolder)) {
    dir.create(indexFolder, recursive = TRUE)
}

indexHairpinFiles <- glue("{spcCode}_hairpin_dna_final")

cli_process_start(" Indexing hairpin files using bowtie")
command <- glue("bowtie-build -q {finalhairpinFile} {indexFolder}/{indexHairpinFiles}")
system(command)

bowtieFiles <- list.files(indexFolder, pattern = indexHairpinFiles)

if (length(bowtieFiles) > 1){
    cli_process_done()
} else {
    cli_process_failed()
}

cli_alert_success(" Preparing hairpin file... done")

#################################################################################


## Download the bovine reference genome from ensembl: https://ftp.ensembl.org/pub/release-113/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz

## Remove spaces from reads beginning. This will help to use genomes from other sources but UCSC
cat("\n\n");cli_rule(center = "Preparing {opt$specie} reference genome");cat("\n\n")
if (file.exists(opt$referenceGenome)) {
    # opt$referenceGenome <- "Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz"
    # file_path_sans_ext(basename(opt$referenceGenome), compression = TRUE)
    cleanedRefGenome <- glue("{opt$database}/{file_path_sans_ext(basename(opt$referenceGenome), compression = TRUE)}_for_mirdeep.fa")
    
    if (!file.exists(cleanedRefGenome)) {
        cli_process_start(" Removing spaces from file {opt$referenceGenome}")
        if (file_ext(opt$referenceGenome) == "gz") {
            command <- glue("gunzip -c {opt$referenceGenome} | sed 's/ //g' > {cleanedRefGenome}")
        } else if (file_ext(opt$referenceGenome) == "fa|fasta") {
            command <- glue("sed 's/ //g' {opt$referenceGenome} > {cleanedRefGenome}")
        }
        
        system(command)
        
        if (file.exists(cleanedRefGenome)) {
            cli_process_done()
        } else {
            cli_process_failed()
        }
    }
    
    
    
    ## build index reference genome downloaded from ensembl using bowtie
    refIndexFolder <- glue("{opt$database}/{str_remove(basename(cleanedRefGenome), '_for_mirdeep.fa')}")
    
    bowtieFiles <- list.files(refIndexFolder, pattern = "ebwt")
    
    if(!any(file.exists(bowtieFiles))) {
        
        
        if(!dir.exists(refIndexFolder)) {
            dir.create(refIndexFolder, recursive = TRUE)
        }
        
        procs <- prepareCore(opt$procs)
        
        cat("\n\n");cli_h1("Running bowtie-build for {basename(cleanedRefGenome)}");cat("\n\n")
        
        command <- glue("bowtie-build --threads {procs} {opt$database}/{file_path_sans_ext(basename(opt$referenceGenome), compression = TRUE)}_for_mirdeep.fa {refIndexFolder}/{file_path_sans_ext(basename(cleanedRefGenome), compression = TRUE)}")
        # cli_progress_bar("Building index genome from {basename(cleanedRefGenome)}", total = 300)
        system(command)
        # cli_progress_update()
        
        bowtieFiles <- list.files(refIndexFolder, pattern = "ebwt")
        
        if (length(bowtieFiles) > 1){
            cli_alert_success(" bowtie-build finished successfully")
        } else {
            cli_abort(c(" x" = "Something went wrong with bowtie-build"))
        }
        
    }
    
    
    

    
    
    
} else {
    cli_abort(c(" x" = "Please use a valid reference genome file"))
}






