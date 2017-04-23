# feed this function "unfiltered_pype_*.sample.tab files 

apply.filter <- function(data){
  coverage <- data[data$normal_coverage_threshold != "FLAGGED",]
  somatic.filter <- coverage[coverage$BedFilter != "FLAGGED" & coverage$delly_passed_sFilter != "False",]
  tier1 <- somatic.filter[somatic.filter$delly_tier == "TIER1",]
  mapq <- tier1[tier1$MAPQ > 30,]
  mapq$distance_to_low_complexity_1 <- as.character(mapq$distance_to_low_complexity_1)
  mapq$distance_to_low_complexity_2 <- as.character(mapq$distance_to_low_complexity_2)
  
  for (i in 1:nrow(mapq)){
    if(is.na(mapq$distance_to_low_complexity_1[i]) == "TRUE"){
      mapq$distance_to_low_complexity_1[i] <- "FAR"
    } else {}
  }
  
  for (i in 1:nrow(mapq)){
    if(mapq$distance_to_low_complexity_1[i] == "0"){
      mapq$distance_to_low_complexity_1[i] <- "OVERLAP"
    } else {}
  }
  
  for (i in 1:nrow(mapq)){
    if(is.na(mapq$distance_to_low_complexity_2[i]) == "TRUE"){
      mapq$distance_to_low_complexity_2[i] <- "FAR"
    } else {}
  }
  
  for (i in 1:nrow(mapq)){
    if(mapq$distance_to_low_complexity_2[i] == "0"){
      mapq$distance_to_low_complexity_2[i] <- "OVERLAP"
    } else {}
  }
  
  mapq$key <- paste(sep=".", mapq$Chromosome1, mapq$Chromosome2, mapq$Position1, mapq$Position2, mapq$CT, mapq$Sample)
  mapq$duplicated <- duplicated(mapq$key)
  remove.dups <- mapq[mapq$duplicated == "FALSE",]
  
  #flagged by any two
  UM.HD <- remove.dups[remove.dups$unique_mapping == "FLAGGED" & remove.dups$high_depth == "FLAGGED",]
  UM.LCR <- remove.dups[remove.dups$unique_mapping =="FLAGGED" & remove.dups$distance_to_low_complexity_1 == "OVERLAP" & remove.dups$distance_to_low_complexity_2 == "OVERLAP",]
  UM.MM <- remove.dups[remove.dups$unique_mapping == "FLAGGED" & remove.dups$multi_mapping == "FLAGGED",]
  HD.MM <- remove.dups[remove.dups$high_depth == "FLAGGED" & remove.dups$multi_mapping == "FLAGGED",]
  LCR.MM <- remove.dups[remove.dups$distance_to_low_complexity_1 == "OVERLAP" & remove.dups$distance_to_low_complexity_2 == "OVERLAP" & remove.dups$multi_mapping == "FLAGGED",]
  HD.LCR <- remove.dups[remove.dups$high_depth == "FLAGGED" & remove.dups$distance_to_low_complexity_1 == "OVERLAP" & remove.dups$distance_to_low_complexity_2 == "OVERLAP",]
  
  remove.flagged.calls <- remove.dups[!remove.dups$key %in% c(UM.HD$key,UM.LCR$key,UM.MM$key, HD.MM$key, LCR.MM$key, HD.LCR$key),]
                                                
  return(remove.flagged.calls)
}
