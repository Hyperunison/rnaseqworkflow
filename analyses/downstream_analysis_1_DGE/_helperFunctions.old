## Load helper function - will be moved to UnisonRtools package
loadProjectEnvironment <- function(){
    ###############################################################################
    ##                                                                           ##
    if (!requireNamespace("remotes")){
      install.packages("remotes")
    }
    
    
    if (!requireNamespace("renv")){
      remotes::install_github("rstudio/renv")
    }
    
    if (!file.exists("renv.lock")){
        renv::init()
        print("Renv environment initiated")
    } else {
        renv::restore(prompt = FALSE)
        print("Renv environment restored from file")
    }

    
    ## Done                                                                      ##
    ###############################################################################
}


createUnifiedDesignFile <- function(
    projectList = projectList
){
    designFileList <- list()

    for (i in 1:length(projectList)){
        FNdesign <- projectList[[i]]$designFile
        print(FNdesign)
    
        if (file.exists(FNdesign)){
            dfDesign <- read.delim(
                FNdesign, 
                header = TRUE, 
                sep = "\t",
                stringsAsFactor = FALSE
            )
    
            if (!is.null(projectList[[i]]$dataset_id) & !(projectList[[i]]$dataset_id == "")){
                dfDesign[["dataset_id"]] <- projectList[[i]]$dataset_id
            } else {
                dfDesign[["dataset_id"]] <- paste0("dataset_", i)
            }
    
            if (nrow(dfDesign) > 0){
                designFileList[[names(projectList)[i]]] <-   dfDesign
            }
    
            print(
                paste0(
                    "Design file loaded for project ",
                    names(projectList)[i],
                    "."
                )
            )
        }
    
            
          
    }
    
    lapply(designFileList, dim)

    ## Carry over all DGE comparisons defined in individual datasets
    for (i in 1:length(designFileList)){
        names(designFileList[[i]]) <- gsub("comp_", paste0("comp_", LETTERS[i]), names(designFileList[[i]]))
    }
    
    ## Carry over all DGE comparisons defined in individual datasets
    for (i in 1:length(designFileList)){
        names(designFileList[[i]]) <- gsub("LRT_", paste0("LRT_", LETTERS[i], "_"), names(designFileList[[i]]))
    }
    
    ## Carry over all factors defined in individual datasets
    for (i in 1:length(designFileList)){
        names(designFileList[[i]]) <- gsub("f_", paste0("f_", LETTERS[i], "_"), names(designFileList[[i]]))
    }
    
    lapply(designFileList, names)
    
    ###############################################################################
    ## Unify design files                                                        ##
    
    # Find common columns
    colList <- lapply(designFileList, colnames)
    commonCols <- Reduce(intersect, colList)
    
    print(commonCols)
    # Now create a vector with all columns that are present in at least one dataset
    
    relevantCols <- as.vector(NULL, mode="character")
    
    for (i in 1:length(designFileList)){
        pos <- c(
            grep("comp_", names(designFileList[[i]])),
            grep("f_", names(designFileList[[i]])),
            grep("LRT_", names(designFileList[[i]]))
        )
        relevantCols <- c(
            relevantCols, 
            names(designFileList[[i]])[pos]
        )
    }

    ## Missing relevant cols to all design files in the list so they can be combined
    for (i in 1:length(designFileList)){
        missingCols <- relevantCols[!(relevantCols %in% names(designFileList[[i]]))]
        dfAdd <- data.frame(
            matrix(nrow = nrow(designFileList[[i]]), ncol=length(missingCols))
        )
        names(dfAdd) <- missingCols
    
        designFileList[[i]] <- cbind(
            designFileList[[i]],
            dfAdd
        )
    
        if (i ==1){
            mergedDesignFile <- designFileList[[i]][,c(commonCols, relevantCols)]
        } else {
            mergedDesignFile <- rbind(
                mergedDesignFile,
                designFileList[[i]][,c(commonCols, relevantCols)]
            )
        }
    }

    return(mergedDesignFile)
}


createUnifiedMatrix <- function(
    projectList = projectList,
    fileType = "geneLevelFeatureCountFile"
){
    ###############################################################################
    ## Create merged read count matrix                                           ##
    
    
    countFileList <- list()
    #fileType <- "geneLevelFeatureCountFile"
    
    for (i in 1:length(projectList)){
        FNcount <- projectList[[i]][[fileType]]
        print(FNcount)
    
        if (file.exists(FNcount)){
            dfRSEM <- read.delim(
                FNcount, 
                header = TRUE, 
                sep = "\t",
                stringsAsFactor = FALSE
            )
    
            ## Remove extra column
            dfRSEM$transcript_id.s. <- NULL
            
            if (nrow(dfRSEM) > 0){
                countFileList[[names(projectList)[i]]] <-   dfRSEM
            }
    
            print(
                paste0(
                    "Files loaded for project ",
                    names(projectList)[i],
                    "."
                )
            )
        }
    
            
          
    }
    
    lapply(countFileList, dim)
    
    ###############################################################################
    ## Create one single featureCount matrix                                     ##
    
    
    
    for (i in 1:length(countFileList)){
        dfTemp <- countFileList[[names(projectList)[i]]]
    
        if (i == 1){
            dfRes <- dfTemp
            refGeneID <- projectList[[i]]$primaryAlignmentGeneID
        } else {
            if (refGeneID != projectList[[i]]$primaryAlignmentGeneID){
                stop("Alignment gene ids in the read count files don't match.")
            }
    
            dfRes <- merge(
                dfRes, 
                dfTemp, 
                by.x = refGeneID, 
                by.y = refGeneID, 
                all = TRUE
            )
    
            dfRes[is.na(dfRes)] <- 0
        }   
    }
    
    dim(dfRes)
    print(
        paste0(
            "Merged feature count matrix for ",
            ncol(dfRes) - 1,
            " samples created."
        )
    )
    
    #dfRSEM <- dfRes
    #dim(dfRSEM)
    ## Done creating merged read count matrix  ##
    return(dfRes)

}