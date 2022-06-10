library(magrittr)


parse_open_target_association_data <- function(json_fname){
  
  b <- jsonlite::fromJSON(gzfile(json_fname))
  disease_entries <- NULL
  if(!is.null(b$disease)){
    disease_entries <- data.frame('disease_efo_id' = b$disease$id, stringsAsFactors = F)
    if(!is.null(b$disease$efo_info)){
      disease_entries$disease_efo_label <- b$disease$efo_info$label
      if(!is.null(b$disease$efo_info$therapeutic_area)){
        disease_entries$disease_therapeutic_area_codes <- sapply(b$disease$efo_info$therapeutic_area$codes, paste, collapse=", ")
        disease_entries$disease_therapeutic_area_labels <- sapply(b$disease$efo_info$therapeutic_area$labels, paste, collapse=", ")
      }else{
        disease_entries$disease_therapeutic_area_codes <- NA
        disease_entries$disease_therapeutic_area_labels <- NA
      }
    }
  }else{
    disease_entries <- data.frame('disease_efo_id' = rep(NA,nrow(b)), stringsAsFactors = F)
    for(c in c('disease_efo_label','disease_therapeutic_area_code',
               'disease_therapeutic_area_label')){
      disease_entries[,c]<- NA
    }
  }
  
  target_entries <- b$target$gene_info
  target_entries$target_id <- b$target$id
  target_entries <- target_entries %>% 
    dplyr::rename(target_name = name, target_symbol = symbol)
  
  association_scores <- b$association_score$datasources
  colnames(association_scores) <- paste0('association_',colnames(association_scores))
  association_scores$association_overall <- b$association_score$overall
  for(cname in colnames(association_scores)){
    association_scores[,cname] <- round(association_scores[,cname], digits = 6)
  }
  
  association_scores$association_is_direct <- b$is_direct
  
  tractability <- data.frame('tractability_small_molecule' = b$target$tractability$smallmolecule$top_category, stringsAsFactors = F)
  tractability$tractability_antibody = b$target$tractability$antibody$top_category
  tractability <- tractability %>% 
    dplyr::mutate(tractability_antibody = dplyr::if_else(tractability_antibody == "Predicted Tractable - High confidence","Predicted_Tractable_HC",tractability_antibody)) %>%
    dplyr::mutate(tractability_antibody = dplyr::if_else(tractability_antibody == "Predicted Tractable - Medium to low confidence","Predicted_Tractable_M_LC",tractability_antibody))
  evidence_df <- dplyr::bind_cols(target_entries, disease_entries, association_scores, tractability)
  
}

parse_open_target_evidence_data <- function(json_fname){
  b <- jsonlite::fromJSON(gzfile(json_fname))
  
  pmids <- NULL
  if(!is.null(b$literature)){
    if(!is.null(b$literature$references)){
      m <- 1
      while(m <= length(b$literature$references)){
        df <- data.frame('association_pmid' = NA, stringsAsFactors = F)
        if(!is.null(b$literature$references[[m]])){
          df <- data.frame('association_pmid' = stringr::str_replace_all(paste(b$literature$references[[m]]$lit_id,collapse=","),"http://europepmc.org/abstract/MED/",""), stringsAsFactors = F)
        }
        pmids <- dplyr::bind_rows(pmids, df)
        m <- m + 1
      }
    }
  }
  
  target_entries <- b$target$gene_info
  target_entries <- dplyr::select(target_entries, -geneid)
  for(field in c('target_type','id','activity','complex_id')){
    fieldt <- field
    if(field != 'target_type'){
      fieldt <- paste0('target_',field)
      if(field == 'id'){
        fieldt <- 'target_ensembl_gene_id'
      }
    }
    if(field %in% colnames(b$target)){
      target_entries[,fieldt] <- b$target[,field]
    }else{
      target_entries[,fieldt] <- NA
    }
  }
  for(field in c('sourceID','type','id','access_level')){
    if(field %in% colnames(b)){
      fieldt <- paste0('association_',field)
      target_entries[,fieldt] <- b[,field]
    }else{
      target_entries[,fieldt] <- NA
    }
  }
  
  target2drug <- NULL
  if(!is.null(b$evidence$target2drug)){
    target2drug <- data.frame('drug_action_type' = b$evidence$target2drug$action_type, stringsAsFactors = F)
    target2drug$drug_moa <- b$evidence$target2drug$mechanism_of_action
  }
  
  drug2clinic <- NULL
  if(!is.null(b$evidence$drug2clinic)){
    drug2clinic <- data.frame('drug_clinical_status' = b$evidence$drug2clinic$status)
    drug2clinic$drug_clinical_is_associated <- b$evidence$drug2clinic$is_associated
    drug2clinic$drug_clinical_resource_score_type <- b$evidence$drug2clinic$resource_score$type
    drug2clinic$drug_clinical_resource_score_value <- b$evidence$drug2clinic$resource_score$value
    clinical_info <- data.frame()
    for(i in 1:length(b$evidence$drug2clinic$urls)){
      if(is.null(b$evidence$drug2clinic$urls[[i]])){
        df <- data.frame('drug_clinical_id' = NA, 'drug_clinical_source' = NA)
      }else{
        if("nice_name" %in% colnames(b$evidence$drug2clinic$urls[[i]])){
          if(stringr::str_detect(b$evidence$drug2clinic$urls[[i]]$nice_name, "DailyMed")){
            df <- data.frame('drug_clinical_id' = stringr::str_split(b$evidence$drug2clinic$urls[[i]]$url,"=",n=2)[[1]][2], 'drug_clinical_source' = 'DailyMed', stringsAsFactors = F)
          }
          if(stringr::str_detect(b$evidence$drug2clinic$urls[[i]]$nice_name, "ATC Information")){
            df <- data.frame('drug_clinical_id' = stringr::str_split(b$evidence$drug2clinic$urls[[i]]$url,"=",n=2)[[1]][2], 'drug_clinical_source' = 'ATC', stringsAsFactors = F)
          }
          if(stringr::str_detect(b$evidence$drug2clinic$urls[[i]]$url,"NCT[0-9]{1,}")){
            df <- data.frame('drug_clinical_id' = stringr::str_match(b$evidence$drug2clinic$urls[[i]]$url,"NCT[0-9]{1,}")[[1]], 'drug_clinical_source' = 'clinicaltrials.gov', stringsAsFactors = F)
          }
          if(stringr::str_detect(b$evidence$drug2clinic$urls[[i]]$url,"api\\.fda\\.gov")){
            df <- data.frame('drug_clinical_id' = stringr::str_split(b$evidence$drug2clinic$urls[[i]]$url,":",n=3)[[1]][3], 'drug_clinical_source' = 'FDA', stringsAsFactors = F)
          }
        }
      }
      clinical_info <- dplyr::bind_rows(clinical_info, df)
      i <- i + 1
    }
    drug2clinic <- dplyr::bind_cols(drug2clinic, clinical_info)
  }
  
  disease_entries <- NULL
  if(!is.null(b$disease)){
    disease_entries <- data.frame('disease_efo_id' = b$disease$id, stringsAsFactors = F)
    if(!is.null(b$disease$efo_info)){
      disease_entries$disease_efo_label <- b$disease$efo_info$label
      if(!is.null(b$disease$efo_info$therapeutic_area)){
        disease_entries$disease_therapeutic_area_codes <- sapply(b$disease$efo_info$therapeutic_area$codes, paste, collapse=", ")
        disease_entries$disease_therapeutic_area_labels <- sapply(b$disease$efo_info$therapeutic_area$labels, paste, collapse=", ")
      }else{
        disease_entries$disease_therapeutic_area_codes <- NA
        disease_entries$disease_therapeutic_area_labels <- NA
      }
    }
    if(!is.null(b$disease$biosample$name)){
      disease_entries$disease_biosample_name <- b$disease$biosample$name
    }else{
      disease_entries$disease_biosample_name <- NA
    }
    if(!is.null(b$disease$biosample$id)){
      disease_entries$disease_biosample_id <- b$disease$biosample$id
    }else{
      disease_entries$disease_biosample_id <- NA
    }
  }else{
    disease_entries <- data.frame('disease_efo_id' = rep(NA,nrow(b)), stringsAsFactors = F)
    for(c in c('disease_efo_label','disease_therapeutic_area_code',
               'disease_therapeutic_area_label','disease_biosample_name','disease_biosample_id')){
      disease_entries[,c]<- NA
    }
  }
  
  drug_entries <- NULL
  if(!is.null(b$drug)){
    drug_entries <- data.frame('drug_molecule_name' = b$drug$molecule_name, stringsAsFactors = F)
    drug_entries$drug_molecule_type <- b$drug$molecule_type
    drug_entries$drug_max_phase_all_diseases_numeric <- NA
    drug_entries$drug_max_phase_all_diseases_label <- NA
    if(!is.null(b$evidence$drug2clinic$clinical_trial_phase)){
      drug_entries$drug_max_phase_all_diseases_numeric <- b$evidence$drug2clinic$clinical_trial_phase$numeric_index
      drug_entries$drug_max_phase_all_diseases_label <- b$evidence$drug2clinic$clinical_trial_phase$label
    }
    chembl_compound_id <- data.frame('drug_chembl_compound_id' = stringr::str_match(b$drug$id,"CHEMBL[0-9]{1,}"))
    drug_entries <- dplyr::bind_cols(drug_entries,chembl_compound_id)
  }
  
  provenance_types <- NULL
  if(!is.null(b$evidence$provenance_type)){
    provenance_types <- data.frame('provenance_type_database_version' = b$evidence$provenance_type$database$version, 'provenance_type_database_id' = b$evidence$provenance_type$database$id)
  }
  
  
  known_mutations <- data.frame()
  k <- 1
  while(k <= length(b$evidence$known_mutations)){
    if(!is.null(b$evidence$known_mutations[[k]]) & length(b$evidence$known_mutations[[k]]) > 0){
      known_mutations <- dplyr::bind_rows(known_mutations, data.frame('known_mutations' = paste(b$evidence$known_mutations[[k]]$preferred_name,collapse=", "), stringsAsFactors = F))
    }else{
      known_mutations <- dplyr::bind_rows(known_mutations, data.frame('known_mutations' = NA))
    }
    k <- k + 1
  }
  
  association_score <- NULL
  if(!is.null(b$scores)){
    association_score <- data.frame('association_score' = b$scores$association_score, stringsAsFactors = F)
  }else{
    association_score <- data.frame('association_score' = rep(NA,nrow(b)), stringsAsFactors = F)
  }
  
    
  evidence_items <- NULL
  if(!is.null(b$evidence)){
    if(!is.null(b$evidence$resource_score)){
      evidence_items <- data.frame('evidence_resource_score_type' = b$evidence$resource_score$type, 'evidence_resource_score_value' = b$evidence$resource_score$value, 'evidence_resource_score_method' = b$evidence$resource_score$method$description, stringsAsFactors = F)
      if(!is.null(b$evidence$log2_fold_change)){
         evidence_items$evidence_log2_fold_change_percentile_rank <- b$evidence$log2_fold_change$percentile_rank
         evidence_items$evidence_log2_fold_change_value <- b$evidence$log2_fold_change$value	
      }else{
         evidence_items$evidence_log2_fold_change_percentile_rank <- NA
         evidence_items$evidence_log2_fold_change_value <- NA
      }
      for(c in c('reference_replicates_n','test_replicates_n','confidence_level','experiment_overview',
                'test_sample','comparison_name','unique_experiment_reference','is_associated')){
         col <- paste0('evidence_',c)
         if(c %in% colnames(b$evidence)){
            evidence_items[,col] <- b$evidence[,c]
         }else{
            evidence_items[,col] <- NA
         }
      }
    }
  }
  
  association_ids <- data.frame()
  if(!('species' %in% colnames(b$unique_association_fields))){
    association_ids <- data.frame('species' = rep(NA,nrow(b)), stringsAsFactors = F)
  }else{
    association_ids <- data.frame('species' = b$unique_association_fields$species, stringsAsFactors = F)
  }
  for(field in c('datasource','prediction_model','model_description','disease_phenodigm_id','model_gene_id',
                 'model_genetic_background','disease','disease_acronym','disease_url','disease_uri','original_disease_name',
                 'original_disease_name','original_disease_label','tumor_type','tumor_type_acronym',
                 'pathway_id','publicationIDs')){
    col <- paste0('unique_assoc_',field)
    if(field %in% colnames(b$unique_association_fields)){
      association_ids[,col] <- gsub("http://identifiers.org/ensembl/","",b$unique_association_fields[,field],fixed =T)
    }else{
      association_ids[,col] <- NA
    }
  }
  
  all_association_entries <- dplyr::bind_cols(target_entries, disease_entries, association_ids, association_score)
  
  if(is.null(drug_entries)){
    for(field in c('drug_molecule_name','drug_molecule_type','drug_max_phase_all_diseases_numeric',
                   'drug_max_phase_all_diseases_label','drug_chembl_compound_id')){
      all_association_entries[,field] <- NA
    }
  }else{
    all_association_entries <- dplyr::bind_cols(all_association_entries, drug_entries)
  }
  
  if(is.null(evidence_items)){
    for(field in c('evidence_resource_score_type','evidence_resource_score_method','evidence_log2_fold_change_percentile_rank',
                   'evidence_experiment_overview','evidence_is_associated','evidence_log2_fold_change_value',
                   'evidence_reference_replicates_n','evidence_test_replicates_n','evidence_comparison_name',
                   'evidence_unique_experiment_reference','evidence_test_sample','evidence_resource_score_value',
                   'evidence_confidence_level')){
      all_association_entries[,field] <- NA
    }
  }else{
    all_association_entries <- dplyr::bind_cols(all_association_entries, evidence_items)
  }
  
  if(is.null(target2drug)){
    for(field in c('drug_action_type','drug_moa')){
      all_association_entries[,field] <- NA
    }
  }else{
    all_association_entries <- dplyr::bind_cols(all_association_entries, target2drug)
  }
  
  if(is.null(pmids)){
    for(field in c('association_pmid')){
      all_association_entries[,field] <- NA
    }
  }else{
    all_association_entries <- dplyr::bind_cols(all_association_entries, pmids)
  }
  
  if(is.null(drug2clinic)){
    for(field in c('drug_clinical_status','drug_clinical_id','drug_clinical_source',
                   'drug_clinical_is_associated','drug_clinical_resource_score_type',
                   'drug_clinical_resource_score_value')){
      all_association_entries[,field] <- NA
    }
  }else{
    all_association_entries <- dplyr::bind_cols(all_association_entries, drug2clinic)
  }
  
  if(is.null(provenance_types)){
    for(field in c('provenance_type_database_version','provenance_type_database_id')){
      all_association_entries[,field] <- NA
    }
  }else{
    all_association_entries <- dplyr::bind_cols(all_association_entries, provenance_types)
  }
  
  if(nrow(known_mutations) == 0){
    for(field in c('known_mutations')){
      all_association_entries[,field] <- NA
    }
  }else{
    all_association_entries <- dplyr::bind_cols(all_association_entries, known_mutations)
  }
  
  return(all_association_entries)
}

release <- '201909'

i <- 121
json_fname <- paste0('data-raw/target_validation_evidence_chunk_',i,'.json.gz')
association_entries <- parse_open_target_evidence_data(json_fname)
#association_entries$chunk <- i
#association_entries$num_columns <- length(colnames(association_entries))
#write.table(association_entries,file = paste0("data/open_targets_associations_",release,".tsv"),sep="\t",col.names = T,row.names = F, quote = F)

#for(i in 2:290){
#  json_fname <- paste0('data-raw/target_validation_evidence_chunk_',i,'.json.gz')
#  if(file.exists(json_fname)){
#    association_entries <- parse_open_target_evidence_data(json_fname)
#    association_entries$chunk <- i
#    association_entries$num_columns <- length(colnames(association_entries))
#    if(length(colnames(association_entries)) != 66){
#      if(file.exists(paste0("data/open_targets_associations_",release,"_",length(colnames(association_entries)),".tsv"))){
#        write.table(association_entries,file=paste0("data/open_targets_associations_",release,"_",length(colnames(association_entries)),".tsv"),sep="\t",col.names = F,append = T, row.names = F,quote = F)
#      }else{
#        write.table(association_entries,file=paste0("data/open_targets_associations_",release,"_",length(colnames(association_entries)),".tsv"),sep="\t",col.names = T,row.names = F,quote = F)
#      }
#    }else{
#      write.table(association_entries,file = paste0("data/open_targets_associations_",release,".tsv"),sep="\t",append = T,col.names = F,row.names = F, quote = F)
#    }
#    cat(paste(i,length(colnames(association_entries)),sep=" - "),'\n')
#  }
#}


# i <- 1
# json_fname <- paste0('data-raw/association_data_',i,'.json.gz')
# association_entries <- parse_open_target_association_data(json_fname)
# write.table(association_entries,file = paste0("data/open_targets_disease_associations_",release,".tsv"),sep="\t",col.names = T,row.names = F, quote = F)
# 
# for(i in 2:134){
#   json_fname <- paste0('data-raw/association_data_',i,'.json.gz')
#   if(file.exists(json_fname)){
#     association_entries <- parse_open_target_association_data(json_fname)
#     write.table(association_entries,file = paste0("data/open_targets_disease_associations_",release,".tsv"),sep="\t",append = T,col.names = F,row.names = F, quote = F)
#     cat(paste(i,length(colnames(association_entries)),sep=" - "),'\n')
#   }
#   i <- i + 1
# }





