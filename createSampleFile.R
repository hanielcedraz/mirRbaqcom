#!/usr/bin/env Rscript


if (suppressPackageStartupMessages(!require(pacman))) suppressPackageStartupMessages(install.packages("pacman"))
suppressPackageStartupMessages(
    pacman::p_load(
        dplyr, 
        optparse,
        glue, 
        tibble, 
        stringr,
        cli,
        tidyr,
        data.table,
        tools
    )
)


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
    make_option(
        opt_str = c("-i", "--fastqFolder"), 
        type = "character", 
        default = "00-Fastq",
        help = glue::glue("The base path where the fastq files are stored. [default: %default]"),
        dest = "fastqFolder"
    ),
    make_option(
        opt_str = c("-o", "--output"), 
        type = "character", 
        default = "samples.txt",
        help = "The filename of the sample file. [default: %default]",
        dest = "output"
    ),
    make_option(
        opt_str = c("-n", "--numberOfSamples"),
        type = "numeric",
        default = NULL,
        help = "select the n first samples [default: %default]",
        dest = "numberOfSamples"
    ),
    make_option(
        opt_str = c("-r", "--ramdomSamples"),
        action = "store_true",
        default = FALSE,
        help = "select the --numberOfSamples random samples [default: %default]",
        dest = "ramdomSamples"
    ),
    make_option(
        opt_str = c("-s", "--saveFile"),
        action = "store_true",
        default = FALSE,
        help = "Save samples file [default: %default]",
        dest = "saveFile"
    )
)


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list,  description =  glue::glue('Author: OLIVEIRA, H.C.', 'Version: 1.0.0', 'E-mail: haniel.cedraz@stgen.com', 'Date: {str_sub(Sys.time(), end = 19)}', .sep = "\n"), formatter = IndentedHelpFormatter))

cat("\n\n");cli_rule(center = col_blue("Creating sample file"));cat("\n\n")

if (!dir.exists(opt$fastqFolder)) {
    cli_abort(c("x" = "Trying to read fastq files from {opt$fastqFolder} but it does not exist"))
    
} else {
    ext <- unique(file_ext(dir(file.path(opt$fastqFolder), pattern = "fastq|gz")))
    if (length(ext) == 0) {
        cli_abort(c("x" = "Cannot locate fastq file in folder {opt$fastqFolder}"))
    } else {
        cli_alert_info(" Reading fastq files from {opt$fastqFolder} ");cat("\n\n")
    }
}


if (all(str_detect(list.files(opt$fastqFolder), "PE1|PE2"))) {
    
    sampleTable <- tibble(
        Read_1 = list.files(opt$fastqFolder, pattern = "PE1"),
        Read_2 = list.files(opt$fastqFolder, pattern = "PE2")) 
    
    cli_h1("");cli_alert_success(" Found Paired-End (PE1 and PE2) fastq files");cat("\n")
    
} else if (all(str_detect(list.files(opt$fastqFolder), "R2") & str_detect(list.files(opt$fastqFolder), "R2"))) {
    
    sampleTable <- tibble(
        Read_1 = list.files(opt$fastqFolder, pattern = "R1"),
        Read_2 = list.files(opt$fastqFolder, pattern = "R2")) 
    
    cli_h1("");cli_alert_success(" Found Paired-End (R1 and R2) fastq files");cat("\n")
    
} else if (all(str_detect(list.files(opt$fastqFolder), "R1|SE"))) {
    
    
    sampleTable <- tibble(
        Read_1 = list.files(opt$fastqFolder, pattern = "R1|SE")
    )
    
    if (all(str_detect(list.files(opt$fastqFolder), "R1"))){
        cli_h1("");cli_alert_success(" Found Single-End (R1) fastq files");cat("\n")
    } else if (all(str_detect(list.files(opt$fastqFolder), "SE"))){
        cli_h1("");cli_alert_success(" Found Single-End (SE) fastq files");cat("\n")
    }
    
}

sampleTable <- sampleTable %>% 
    separate(Read_1, into = "SAMPLE_ID", sep = "_", remove = FALSE, extra = "drop") %>% 
    relocate(SAMPLE_ID, .before = Read_1) %>% 
    as.data.frame()



if (!is.null(opt$numberOfSamples)) {
    if (opt$ramdomSamples) {
        cat("\n\n");cli_h1("");cli_alert_info(" Selecting {opt$numberOfSamples} random samples");cat("\n\n")
        
        sampleTable <- sampleTable %>% 
            slice_sample(n = opt$numberOfSamples)
    } else {
        cat("\n\n");cli_h1("");cli_alert_info(" Selecting {opt$numberOfSamples} first samples");cat("\n\n")
        sampleTable <- sampleTable %>% 
            slice_head(n = opt$numberOfSamples)
    }
    
}


if (nrow(sampleTable) > 0) {
    
    
    # cat("\n\n");cli_rule(center = col_blue("Creating samples file"));
    # #cli_li(opt$output)
    
    
    
    
    
    print(as_tibble(sampleTable))
    
    
    
    if (!opt$saveFile) {
        opt$saveFile <- userInput(glue::glue("\n\nWould you like to save the file?? [Y]/N "))
        cat("\n\n")
        
        if (any(opt$saveFile %in% toupper(c("yes", "", "Y")))) {
            
            filePath <- glue("{opt$output}")
            fwrite(sampleTable, file = filePath, sep = "\t", quote = FALSE, row.names = FALSE)
            
            
            if (file.exists(filePath)) {
                cat("\n\n");cli_h1("");cli_alert_success(" {opt$output} saved successfully");cat("\n\n")
                
            } else {
                cat("\n\n");cli_abort(c("x" = "Something went wrong saving the file {opt$output}"))
                
            }
            
        } else if (any(opt$saveFile %in% toupper(c("no", "n")))) {
            break
        }
        
    } else {
        filePath <- glue("{opt$output}")
        fwrite(sampleTable, file = filePath, sep = "\t", quote = FALSE, row.names = FALSE)
        
        
        if (file.exists(filePath)) {
            cat("\n\n");cli_h1("");cli_alert_success(" {opt$output} saved successfully");cat("\n\n")
            
        } else {
            cat("\n\n");cli_abort(c("x" = "Something went wrong saving the file {opt$output}"))
            
        }
    }
    
    
    
    
    
    
    
    
}