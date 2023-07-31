###############################################################################
## Function loadProjectEnvironment                                           ##
#' @title loadProjectEnvironment This function sets a project R-environment.
#'
#'
#' @export

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

## Function end                                                              ##
###############################################################################

###############################################################################
## Function createUnifiedDesignFile                                     ##
#' @title createUnifiedDesignFile
#'
#'
#' @param projectList A list with project parameters. 
#' @param outputfile Specifies a powerpoint output file and a location where it will be saved. 
#' @export

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

## Function end                                                              ##
###############################################################################


###############################################################################
## Function createUnifiedMatrix                                              ##
#' @title createUnifiedMatrix
#'
#'
#' @param projectList A list with paths to matrix files that need to be merged. 
#' @param fileType Must be a field in the projectList. Defines the file type that is to be merged. 
#' @export

createUnifiedMatrix <- function(
    projectList = projectList,
    fileType = "geneLevelFeatureCountFile"
){
    ###########################################################################
    ## Create merged read count matrix                                       ##
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

## Function end                                                              ##
###############################################################################

###############################################################################
## Function createPowerpointPresentation                                     ##
#' @title createPowerpointPresentationFromPlotList
#'
#'
#' @param plotList list with ggplot object to be put into the presentation. List names will appear as slide titles.
#' @param outputfile Specifies a powerpoint output file and a location where it will be saved. 
#' @import officer
#' @export

createPowerpointPresentationFromPlotList <- function(
        plotList, 
        outputfile="singleCell.powerpoint.presentation.pptx"
){
    
    doc <- officer::read_pptx()
    ## Check template presentation ##
    # layout_summary(doc)
    # layout_properties (doc)
    
    for (i in 1:length(plotList)){
        
        ##
        plot <- plotList[[i]]
        
        ## to ensure proper functionality, plot plots to temp file, then import
        # that into powerpoint. 
        tempFN <- "temp.png"
        
        png(tempFN)
            print(plot)
        dev.off()
        
        
        #doc <- add_slide(doc)
        doc <- officer::add_slide(
            doc, layout = "Title and Content", 
            master = "Office Theme"
        )
        
        #doc <- officer::ph_with(
        #    x = doc, 
        #    value = officer::external_img(
        #        src = tempFN
        #    ), 
        #    type = "sldNum", 
        #    location = officer::ph_location_type(
        #        type = "sldNum"
        #    )
        #)
        
        
        doc <- officer::ph_with(
            x = doc, 
            value =  officer::external_img(
                src = tempFN
            ), 
            location = officer::ph_location_type(
                type = "body") 
        )
        
        
        
        doc <- officer::ph_with(
            x = doc, 
            names(plotList)[i], 
            location = officer::ph_location_type(
                type = "title"
            ) 
        )
        doc <- officer::ph_with(
            doc, "Powerpoint Output", 
            location = officer::ph_location_type(type = "ftr")
        )
        
        ## delete temp file
        unlink(tempFN)
    }
    
    print(doc, target = outputfile)
    print(paste0("Presentation generated and saved as ", outputfile))
}
## Function end                                                              ##
###############################################################################


###############################################################################
## Function createHeatmapPlotList                                            ##
#' @title createHeatmapPlotList
#'
#'
#' @param HmDisplayCatsFromDb List with the gene sets to be displayed in heatmaps.
#' @param dfDesign Design file
#' @param dfData Heatmap data file

#' @import ComplexHeatmap
#' @import grid
#' @export

createHeatmapPlotList <- function(
    HmDisplayCatsFromDb,
    dfDesign,
    dfData,
    dfTPM,
    geneIDcolumn = "gene_id"
){
    HMplotList <- list()
    chnkPrefix <- "HM."
    VersionPdfExt <- VersionPdfExt <- paste0(".",chnkPrefix,"V", gsub("-", "", Sys.Date()), ".pdf")
    
    for (k in 1:length(HmDisplayCatsFromDb)){
    
    ## Select samples to display ##
    if (!is.null(HmDisplayCatsFromDb[[k]]$comparisonID)){
        dfSel <- unique(dfDesign[,c("sample.id", HmDisplayCatsFromDb[[k]]$comparisonID)])
        dfSel <- dfSel[dfSel[,HmDisplayCatsFromDb[[k]]$comparisonID] != "",]
        
        if (nrow(dfSel) > 1){
            sampleSelection <- paste0(unique(dfSel$sample.id))    
        } else {
            sampleSelection <- paste0(unique(dfDesign$sample.id))
        }
        
    } else {
        sampleSelection <- paste0(unique(dfDesign$sample.id))
    }
  
    ## Check ##
     
    sampleSelection <- unique(names(dfData)[unlist(sapply(paste0("^", sampleSelection, "$"), function(x) grep(x, names(dfData))))])
    selVec <- c(geneIDcolumn, sampleSelection )
    ## Get gene selection 
    geneSel <- HmDisplayCatsFromDb[[k]]$geneVec
    
    geneSel <- unique(geneSel)
    geneSel <- geneSel[geneSel != ""]
    
    if (length(geneSel) > 2){
        dfDataTable <- dfTPM
        dfDataTable <- unique(dfDataTable[dfDataTable[, geneIDcolumn] %in% geneSel, selVec])
        
        dfHmBase <- unique(dfDataTable[,selVec])
        
        while (sum(duplicated(dfHmBase[, geneIDcolumn])) > 0){
            dfHmBase[duplicated(dfHmBase[, geneIDcolumn]), geneIDcolumn] <- paste0(
                dfHmBase[duplicated(dfHmBase[, geneIDcolumn]), 
                geneIDcolumn], "_", i
            )
            i=i+1
        }
        
        row.names(dfHmBase) <- dfHmBase[, geneIDcolumn]
        dfHmBase[, geneIDcolumn] <- NULL
        
        ## calculate row-means ##
        rowMeans <- apply(
            dfHmBase,
            1,
            function(x) mean(x)
        )
            
        rowMeans[rowMeans ==0] <- 0.001
            
        hmMax <- 4
        for (i in 1:ncol(dfHmBase)){
            dfHmBase[,i] <- log2(dfHmBase[,i] / rowMeans)
        }
            
        dfHmBase[dfHmBase > hmMax] <- hmMax
        dfHmBase[dfHmBase < -1*hmMax] <- -1*hmMax
            
            
        names(dfHmBase) <- gsub("norm_counts_", "", names(dfHmBase))
        names(dfHmBase) <- gsub("_TPM", "", names(dfHmBase))
            
        mHmBase <- data.matrix(dfHmBase)
            
        if ( nrow(mHmBase) < 51){
            showRowNames <- TRUE
        } else {
            showRowNames <- FALSE
        }
        
        ## Create heatmap plot ##
        #library(ComplexHeatmap)
       
        f1 = circlize::colorRamp2(seq(-4, 4, length = 3), c("#3060cf", "#fffbbc","#c4463a"))    
    
        anno <- as.data.frame(colnames(mHmBase))
        colnames(anno) <- "Sample"
        anno$Group <- sapply(as.vector(anno[,1]), function(x) paste0(unlist(strsplit(x, "_"))[1], "_",unlist(strsplit(x, "_"))[2]))
        
        ## Color sample groups in line with the designated sample group color ##
        #######################################################################
        ## Add sample group colors if needed
        pos <- grep("sample.group_color", names(dfDesign))
        
        if (length(pos) == 0){
            sample.group <- unique(dfDesign$sample.group)
            sample.group_color <- sample.group
            #library(scales)
            sample.group_color = scales::hue_pal()(length(sample.group_color))
            #sample.group_color = c("#990000", "#009900")
            dfGroupColors <- unique(data.frame(sample.group, sample.group_color))
            dfDesign <- merge(dfDesign, dfGroupColors, by.x = "sample.group", "sample.group")
            if (exists("Obio")){
                Obio@dfDesign <- dfDesign
            }
            
        }
        
        
        
        #library(scales)
        #hue_pal()(2)
        names(dfDesign)
        df <- unique(data.frame(dfDesign[,c("sample.id", "sample.group", "sample.group_color", "dataset_id")]))
        head(df)
        df <- df[df$sample.id %in% colnames(mHmBase),]
                             
        df2 <- data.frame(df[,c("sample.group")])
        names(df2) <- c("Group")
        
        # Set sample group colors based on design file        
        GroupVec <- as.vector(unique(df$sample.group_color))
        names(GroupVec) <- as.vector(unique(df$sample.group))
        
        # Set dataset colors
        DatasetVec <- rainbow(length(unique(df$dataset_id)))
        names(DatasetVec) <- as.vector(unique(df$dataset_id))
        
        
        ha = ComplexHeatmap::HeatmapAnnotation(
            df = df2, 
            col = list(
                Group = GroupVec,
                Dataset = DatasetVec
            )
        )
    
        ComplexHeatmap::ht_opt(
            legend_border = "black",
            heatmap_border = TRUE,
            annotation_border = TRUE
        )
        
        hmTitle <- unlist(strsplit(names(HmDisplayCatsFromDb)[k], "_padj_"))
        if (length(hmTitle) == 2){
            hmTitle <- paste0("padj_", hmTitle[2])
        } else {
            hmTitle <- names(HmDisplayCatsFromDb)[k]
        }


                             
        HMplotList[[names(HmDisplayCatsFromDb)[k]]] = ComplexHeatmap::Heatmap(
            mHmBase,
            column_title = gsub(
                    "_", 
                    " ", 
                    hmTitle
            ),
            name = paste0("HM_", k), 
            row_km = 3,
            col = f1,
            column_split = factor(df2$Group, levels = unique(df2$Group)),
           
            show_column_names = F,
            show_row_names = showRowNames,
            border = TRUE,

         
            #Dendrogram configurations: columns
            clustering_distance_columns="euclidean",
            clustering_method_columns="complete",
            column_dend_height=grid::unit(10,"mm"),
            
            #Dendrogram configurations: rows
            clustering_distance_rows="euclidean",
            clustering_method_rows="complete",
            row_dend_width=grid::unit(10,"mm"),
            top_annotation = ha,
            show_heatmap_legend = TRUE
            #row_title = NULL,
            #show_row_dend = FALSE
        ) 
        
    ComplexHeatmap::ht_opt(RESET = TRUE)
        
    if (! is.null(HmDisplayCatsFromDb[[k]]$cat_id)){
        link <- paste0(
            'An interactive version of this heatmap with an option for further filtering can be found <a href="',
            "https://biologic.crick.ac.uk/",
            project_id,"/category-view/",
            HmDisplayCatsFromDb[[k]]$cat_id,'" target="_blank">here</a>.'
        )
        
    } else {
        link <- ""
    }
    
    ###########################################################################
    ## Save plot to file                                                     ##
    FNbase <- paste0("Heatmap.", names(HmDisplayCatsFromDb)[k],VersionPdfExt)
    FN <- paste0(figDir, FNbase)
    FNrel <- paste0("report_figures/", FNbase)
    
    pdf(FN)
        print(HMplotList[[names(HmDisplayCatsFromDb)[k]]])
    dev.off()
    ##                                                                       ##
    ###########################################################################
    
    figCap <- paste0(
    '**Figure ',
    figureCount,
    ':** Heatmap showing the gene category ', gsub('_', ' ', names(HmDisplayCatsFromDb)[k]), '. ',
        'Download a pdf of this figure <a href="',FNrel,'" target="_blank">here</a>. ',
        link
    )
    
    figureCount <- figureCount + 1 
    
    NewChnk <- paste0(
            "## HM_", names(HmDisplayCatsFromDb)[k],
            "\n```{r, results='asis', echo=F, eval=TRUE, warning=FALSE, fig.cap='",figCap,"'}\n",
            "\n",
            "\n print(HMplotList[['",names(HmDisplayCatsFromDb)[k],"']])",
            "\n cat(  '\n')",
            "\n\n\n```\n"   
    )
    
    chnkVec <- c(
        chnkVec,
        NewChnk
    )
    
    } ## End making heatmap 
}
return(HMplotList)
}
    
## Done making heatmaps                                                      ##
###############################################################################
    


