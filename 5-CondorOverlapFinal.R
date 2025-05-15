#' Tracking solutions to a persistent threat: 
#' Spatial movement patterns reflect lead exposure in critically endangered 
#' California Condors

#' Script by: Varalika Jain

#' This script relies on 0-CondorMetricsFinal
#' The output dataset 'overlap_final.csv' is available 
#' on the Movebank data repository associated with this study

####----LIBRARIES----####
library(dplyr)
library(amt)
library(sf)
library(ggplot2)

#### (1) Prepare working space----
rm(list = ls())
dev.off()
#load("Condor_final.RData") #from 0-CondorMetricsFinal

# keep track object
rm(list=setdiff(ls(), "intrpl_track"))
#save.image("intrpl_track.RData") 

####----(2) Project and nest track----####
names(intrpl_track)
xx <- intrpl_track
coordinates(xx) <-  c("coords.x1", "coords.x2")
proj4string(xx) <- CRS("+proj=longlat +datum=WGS84")
xx <- spTransform(xx, CRS("EPSG:26712"))
#' for more info, see https://epsg.io/26712
condors_proj <- as.data.frame(xx)

xx <- split(xx, xx$id_test)

#' For each date, create track for each individual
for (i in seq_along(xx)){
  xx[[i]] <- make_track(xx[[i]],
                        .x = coords.x1, 
                        .y = coords.x2, 
                        .t = timestamps.1,
                        id = id,
                        crs = "EPSG:26712")
}

####----(3) Compute aKDE----####
kde <- list()
#homerange_sf <- list()
#area  <- list()

#' using the track data, create home ranges 
for (j in seq_along(xx)){
  kde[[j]] <- hr_akde(xx[[j]], 
                      model = fit_ctmm(xx[[j]], "auto"),
                      keep.data = TRUE,
                      trast = make_trast(xx[[j]]),
                      levels = 0.95,
                      wrap = FALSE)
  #  homerange_sf[[j]] <- hr_isopleths(kde[[j]])
  #  st_write(homerange_sf[[j]]$geometry, 
  #       paste0(names(xx)[[j]],".kml"), driver = "KML")
  # area[[j]] <- hr_area(kde[[j]])
}

####----(4) Extract values----####
area_df <- data.frame(matrix(unlist(area), nrow=length(area), byrow=TRUE))

area_df <- area_df %>%
  dplyr::select(X3,X7,X8,X9) %>%
  dplyr::rename(kde_95 = X3,
                lci_95 = X7,
                estimate_95 = X8,
                uci_95 = X9) %>%
  mutate(
    id = names(xx)
  )

#' convert to km
area_df <- area_df %>% dplyr::mutate_at(c(2:4), as.numeric) %>%
  mutate(across(.cols = c(2:4), .fns = ~.x / 1000000))

#' merge with intrpl track
cols = intrpl_track %>% group_by(id_test, id, test_date, Hatch.Date, lead_level) %>%
  tally()
names(cols)[names(cols) == "id"] <- "raw_id"

area_final <- merge(area_df, cols, by.x = "id", by.y = "id_test", all.x = TRUE)

#' exposed vs unexposed continuous KDE along weeks prior to testing
area_final %>%
  ggplot(aes(x = lead_level, y = estimate_95))+
  geom_boxplot()+
  xlab("Lead level") + ylab ("KDE area in"~km^2)+
  theme_classic()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        legend.position = "none")

write.csv(area_final, 'area_final.csv')

####----(5) Compute the overlap----####
#' add kernal density column to track list
area_final$kde <- kde

####----() KDE 95 overlaps with hr method----####
xx = area_final[c(-171),] #issues

hr_95 <- hr_overlap(xx$kde, type = "hr",
                    which = "all",
                    conditional = TRUE)

save(hr_95, file = "hr_95_overlap.RData")

names(hr_95)[names(hr_95) == "overlap"] <- "overlap_95"

hr_95$what = rep(c("lci", "estimate", "uci"), times = 32580)

overlap_final <- hr_95 %>% filter(what == "estimate")

xx$num <- seq(from = 1, to = nrow(xx))
names(xx)

overlap_final <- merge(overlap_final, xx[,c(1:9,12)], #remove hr_kde column
                       by.x = c("from"),
                       by.y = c("num"),
                       all.x = TRUE)

overlap_final <- merge(overlap_final, xx[,c(1,6,9,12)], 
                       by.x = c("to"),
                       by.y = c("num"),
                       all.x = TRUE)

names(overlap_final)
names(overlap_final)[names(overlap_final) == 'raw_id.x'] <- 'from_raw_id'
names(overlap_final)[names(overlap_final) == 'id.x'] <- 'from_id'
names(overlap_final)[names(overlap_final) == 'raw_id.y'] <- 'to_raw_id'
names(overlap_final)[names(overlap_final) == 'id.y'] <- 'to_id'
names(overlap_final)[names(overlap_final) == 'lead_level.x'] <- 'from_lead'
names(overlap_final)[names(overlap_final) == 'lead_level.y'] <- 'to_lead'


overlap_final$dyad <- apply(overlap_final[, c("from", "to")], 1, function(x) paste0(sort(x), collapse = ""))
overlap_final$lead_dyad <- paste0(overlap_final$from_lead, "-", overlap_final$to_lead)

write.csv(overlap_final, "overlap_final.csv")

#save.image("monthly_overlap.RData")
#load("monthly_overlap.RData")

xy = area_final %>% group_by(lead_level) %>% 
  dplyr::summarise (mean = mean(estimate_95), min = min(estimate_95), max = max(estimate_95))

round(xy[,2:4], 3)




