#' A function to turn non-corrected into corrected labelled peptide/protein fractions (Corr LPFs)
#'
#' This function allows users to input their non-corrected LPF matrix from "AnnotateProteins.R" and correct the fractional enrichment in individual peptides using the enrichment percentages in amino acid soluble pools from paired treatments. The function returns the corrected matrix, which can be directly used to calculate fractional synthesis rates.
#' @param qPCR_List_OR_Dir Directory where the depurated qPCR matrix is stored.
#' @param Treatment_Table_Dir Directory where the position matrix is stored.
#' @keywords qPCR deltadeltaCt
#' @export
#' @examples
#'
#' ...

DeltaDeltaCtCc <- function(qPCR_List_OR_Dir,
                           Treatment_Table_Dir) {



  ## needed functions

  ImportPlates2RawCt_DF <- function(PlatesDIR,
                                    InputPathTreatmentFile){

    qPCR_pos_table <- read.table(InputPathTreatmentFile,
                                 fill = TRUE, header = T)


    qPCR_pos_table_Look <- as.data.frame(cbind(pos_ID = paste0(qPCR_pos_table$Pos, "_", qPCR_pos_table$Plate_ID),
                                               Treat_ID = paste0(qPCR_pos_table$Gene_ID, "_", qPCR_pos_table$Sample_ID)))

    factor_pos_table <- levels(as.factor(qPCR_pos_table_Look$Treat_ID))

    out_vec_RawCt <- c()

    names_vec_RawCt <- c()

    if(class(PlatesDIR) == "list"){

      ## code for multiple plates

      for (i in 1:length(factor_pos_table)) {

        RegExpr_Look4Pattern2 <- paste0("^",
                                        factor_pos_table[i],
                                        "$", collapse = "|")


        indexOfPattern <- grep(pattern = RegExpr_Look4Pattern2,
                               x = qPCR_pos_table_Look$Treat_ID)

        positionsofSamples <- qPCR_pos_table_Look$pos_ID[indexOfPattern]


        ## need to iterate only based on the number of plates and not actually the number of samples so that the complementary information from distinct plates can then be consolidated into joined entries.

        plate_index <- list2DF(strsplit(positionsofSamples, split = "_"))[2,]

        if (class(plate_index) == "character") {

          plate_index = as.numeric(plate_index)

          #number_of_plates <- length(plate_index)

        } else if (class(plate_index) == "data.frame") {

          plate_index = unique(apply(X = plate_index, MARGIN = 2, FUN = as.numeric))

          #number_of_plates <- length(unique(apply(X = plate_index, MARGIN = 2, FUN = as.numeric)))

        }

        #if(number_of_plates != length(PlatesDIR)){

        # stop("ERROR: Number of plates in treatment file is not equal to number of plates in input list")

        #}

        for (j in 1:length(plate_index)) {

          ### import plate j

          qPCR_table_runner <- read.table(file = PlatesDIR[[plate_index[j]]], fill = TRUE)

          ### locate sample or samples in plate j

          Samples_j_DF_runner <- list2DF(strsplit(positionsofSamples,
                                                  split = "_"))

          vector_or_DF_of_samples <- Samples_j_DF_runner[1,which(Samples_j_DF_runner[2,] == plate_index[j])]

          if (class(vector_or_DF_of_samples) == "character"){

            All_Samples_plate_j <- vector_or_DF_of_samples

            RegExpr_Look4Pattern <- paste0("^",
                                           All_Samples_plate_j,
                                           "$", collapse = "|")

            indexofRawCt <- grep(pattern = RegExpr_Look4Pattern,
                                 x = qPCR_table_runner$V3)

            Raw_Ct <- as.numeric(qPCR_table_runner$V6[indexofRawCt])

            ### change output vectors

            out_vec_RawCt <- c(out_vec_RawCt, Raw_Ct)

            names_vec_RawCt <- c(names_vec_RawCt, paste0(factor_pos_table[i], ";", j))


          } else if (class(vector_or_DF_of_samples) == "data.frame"){

            All_Samples_plate_j <- apply(vector_or_DF_of_samples, 2, as.character)

            for (m in 1:length(All_Samples_plate_j)) {

              RegExpr_Look4Pattern <- paste0("^",
                                             All_Samples_plate_j[m],
                                             "$", collapse = "|")

              indexofRawCt <- grep(pattern = RegExpr_Look4Pattern,
                                   x = qPCR_table_runner$V3)

              Raw_Ct <- as.numeric(qPCR_table_runner$V6[indexofRawCt])

              ### output vectors

              out_vec_RawCt <- c(out_vec_RawCt, Raw_Ct)

              names_vec_RawCt <- c(names_vec_RawCt, paste0(factor_pos_table[i], ";", j))

            }
          }
        }

        names(out_vec_RawCt) <- names_vec_RawCt

      }

      out_df <- cbind(Raw_Ct = out_vec_RawCt,
                      gene_treat_ID = apply(list2DF(strsplit(names(out_vec_RawCt), ";"))[1,], 2, as.character),
                      plate_ID = apply(list2DF(strsplit(names(out_vec_RawCt), ";"))[2,], 2, as.character))


    } else if(class(PlatesDIR) == "character") {

      ## code for one plate

      for (i in 1:length(factor_pos_table)) {

        RegExpr_Look4Pattern2 <- paste0("^",
                                        factor_pos_table[i],
                                        "$", collapse = "|")


        indexOfPattern <- grep(pattern = RegExpr_Look4Pattern2,
                               x = qPCR_pos_table_Look$Treat_ID)

        positionsofSamples <- qPCR_pos_table_Look$pos_ID[indexOfPattern]


        ## need to iterate only based on the number of plates and not actually the number of samples so that the complementary information from distinct plates can then be consolidated into joined entries.

        plate_index <- list2DF(strsplit(positionsofSamples, split = "_"))[2,]

        if (class(plate_index) == "character") {

          plate_index = as.numeric(plate_index)

          #number_of_plates <- length(plate_index)

        } else if (class(plate_index) == "data.frame") {

          plate_index = unique(apply(X = plate_index, MARGIN = 2, FUN = as.numeric))

          #number_of_plates <- length(unique(apply(X = plate_index, MARGIN = 2, FUN = as.numeric)))

        }

        #if(number_of_plates != length(PlatesDIR)){

        # stop("ERROR: Number of plates in treatment file is not equal to number of plates in input list")

        #}

        for (j in 1:length(plate_index)) {

          ### import plate j

          qPCR_table_runner <- read.table(file = PlatesDIR, fill = TRUE)

          ### locate sample or samples in plate j

          Samples_j_DF_runner <- list2DF(strsplit(positionsofSamples,
                                                  split = "_"))

          vector_or_DF_of_samples <- Samples_j_DF_runner[1,which(Samples_j_DF_runner[2,] == plate_index[j])]

          if (class(vector_or_DF_of_samples) == "character"){

            All_Samples_plate_j <- vector_or_DF_of_samples

            RegExpr_Look4Pattern <- paste0("^",
                                           All_Samples_plate_j,
                                           "$", collapse = "|")

            indexofRawCt <- grep(pattern = RegExpr_Look4Pattern,
                                 x = qPCR_table_runner$V3)

            Raw_Ct <- as.numeric(qPCR_table_runner$V6[indexofRawCt])

            ### change output vectors

            out_vec_RawCt <- c(out_vec_RawCt, Raw_Ct)

            names_vec_RawCt <- c(names_vec_RawCt, paste0(factor_pos_table[i], ";", j))


          } else if (class(vector_or_DF_of_samples) == "data.frame"){

            All_Samples_plate_j <- apply(vector_or_DF_of_samples, 2, as.character)

            for (m in 1:length(All_Samples_plate_j)) {

              RegExpr_Look4Pattern <- paste0("^",
                                             All_Samples_plate_j[m],
                                             "$", collapse = "|")

              indexofRawCt <- grep(pattern = RegExpr_Look4Pattern,
                                   x = qPCR_table_runner$V3)

              Raw_Ct <- as.numeric(qPCR_table_runner$V6[indexofRawCt])

              ### output vectors

              out_vec_RawCt <- c(out_vec_RawCt, Raw_Ct)

              names_vec_RawCt <- c(names_vec_RawCt, paste0(factor_pos_table[i], ";", j))

            }
          }
        }

        names(out_vec_RawCt) <- names_vec_RawCt

      }

      out_df <- cbind(Raw_Ct = out_vec_RawCt,
                      gene_treat_ID = apply(list2DF(strsplit(names(out_vec_RawCt), ";"))[1,], 2, as.character),
                      plate_ID = apply(list2DF(strsplit(names(out_vec_RawCt), ";"))[2,], 2, as.character))

    }

    return(out_df)

  }


  ## main flow:

  RawCt_DF <- ImportPlates2RawCt_DF(PlatesDIR = qPCR_List_OR_Dir,
                                    InputPathTreatmentFile = Treatment_Table_Dir)


  #return(RawCt_DF)

  RawCt_DF[c(which(is.na(as.character(RawCt_DF[,1])))),1] <- 0

  ### building controls data matrices

  #### minus RT control

  MinusRevT_mat <- RawCt_DF[c(grep(pattern = MinusRevTranscriptaseC_Pattern, x = as.character(RawCt_DF[,"gene_treat_ID"]))),]

  MinusRevT_mat[c(which(as.numeric(MinusRevT_mat[,1]) > MinusRevTranscriptaseC_boundary)),1] <- 0

  #### empty wells control

  EmptyWells_mat <- RawCt_DF[c(grep(pattern = EmptyWells_Pattern, x = as.character(RawCt_DF[,"gene_treat_ID"]))),]

  EmptyWells_mat[c(which(as.numeric(EmptyWells_mat[,1]) > EmptyWells_boundary)),1] <- 0

  #### H2O control

  H2O_mat <- RawCt_DF[c(grep(pattern = H2O_Pattern, x = as.character(RawCt_DF[,"gene_treat_ID"]))),]

  H2O_mat[c(which(as.numeric(H2O_mat[,1]) > H2O_boundary)),1] <- 0


  #### remove controls from original matrix

  Control_depleted_mat <- RawCt_DF[-c(grep(pattern = MinusRevTranscriptaseC_Pattern, x = as.character(RawCt_DF[,"gene_treat_ID"])),
                                      grep(pattern = EmptyWells_Pattern, x = as.character(RawCt_DF[,"gene_treat_ID"])),
                                      grep(pattern = H2O_Pattern, x = as.character(RawCt_DF[,"gene_treat_ID"]))),]


  #### for loop to calculate Delta Cts either using one or multiple controls

  for (i in 1:length(PlatesDIRs)) {

    plate_i_subset <- Control_depleted_mat[which(as.numeric(Control_depleted_mat[,3]) == i),]

    #controlNr <- readline(prompt = "How many controls would you like to use for calculating DeltaCts (return 1 for only one control or return 2 for multiple)")

    cat(paste0(paste0(1:length(as.character(plate_i_subset[,2])), ": ", as.character(plate_i_subset[,2])), collapse = " \n "))

    control_indexes <- readline(prompt = "Could you please return the indexes of your controls ???? (please use numbers separated by spaces))")


  }

  return(control_indexes)

}


