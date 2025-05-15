#' Tracking solutions to a persistent threat: 
#' Spatial movement patterns reflect lead exposure in critically endangered 
#' California Condors

#' Script by: Varalika Jain

#' The movement data on the project are not open publicly.
#' For access to this and other sensitive input datasets 
#' (i.e., blood lead testing data, demography) 
#' associated with this study, please contact me 
#' The final processed dataset 'metrics_final' is available on the Movebank data
#' repository associated with this study

####----LOAD LIBRARIES----####
library(sp)
library(terra)
library(move)
library(dplyr)
library(atlastools)
library(ggmap)
library(ggplot2)
library(ggsn)
library(amt)

####----DOWNLOAD CONDOR DATA FROM MOVEBANK----####
####----(1) Download----####
#' Permissions required to download data
#' Input movebank login details
login <- movebankLogin()

#' Source data from GPS tagged condors from movebank
condors <- getMovebankData(study="California Condor Recovery Program (AZ-UT)", 
                           login=login,
                           timestamp_start = "20141201000000000", 
                           timestamp_end = "20230130000000000",
                           removeDuplicatedTimestamps=TRUE)

#' Assigning the correct time zone
timestamps(condors) <- lubridate::with_tz(timestamps(condors), tz="US/Arizona")

#' Check the timezone
head(timestamps(condors))

#' Check individuals
levels(condors@trackId)

####----(2) Convert to dataframe and filter----####
#' convert to dataframe 
condors_df <- as(condors, "data.frame")

#' keep only GPS data for HMM 
condors_filter <- condors_df %>% subset(sensor_type == "GPS")
range(condors_filter$timestamps)

#' keep certain columns
condors_filter <- condors_filter[,c("location_lat",
                                    "location_long",
                                    "timestamps",
                                    "trackId",
                                    "tag_local_identifier",
                                    "sensor_type",
                                    "heading",
                                    "height_raw",
                                    "ground_speed")]

#' remove 'tes' rows
levels(as.factor(condors_filter$trackId))
condors_filter <- condors_filter %>% 
  filter(!grepl("Test|Test.D|Test.G|Test.I|Test.J|Test.K|Test.M| Test11|Removed.GPS|TBD.1", trackId)) %>%
  droplevels()
levels(as.factor(condors_filter$trackId))

#' Select studbook number from trackid column 
class(condors_filter$trackId)
condors_filter$trackId <- as.character(condors_filter$trackId)
condors_filter$id <- substr(condors_filter$trackId, 2, nchar(condors_filter$trackId)-3)
levels(as.factor(condors_filter$id))

#' convert to factor variable
condors_filter$id <- as.factor(condors_filter$id)
levels(condors_filter$id)
condors_filter$id <- plyr::revalue(condors_filter$id, c("9"="947",
                                                  "383F" = "383",
                                                  "447F0" = "447",
                                                  "520M" = "520")) 
unique(condors_filter$id)
#' 69 tagged birds 

#' create date column 
class(condors_filter$timestamps)
condors_filter$date <- as.Date(condors_filter$timestamps, tz = "US/Arizona")
range(condors_filter$date)
#' "2015-01-14" "2023-01-29"

head(condors_filter)
condors_filter <- subset(condors_filter, select = -c(trackId))
head(condors_filter)

####----REMOVING OUTLIERS----####
####----(1) Large-scale outliers----####
range(condors_filter$location_lat)
range(condors_filter$location_long)

hist(condors_filter$location_lat, breaks = 80, ylim = c(0,10))
hist(condors_filter$location_long, breaks = 80, ylim = c(0,10))

####----(2) Remove outliers according to spatial limits----####
#' some clear outliers which will be removed according to spatial limits
condors_atlas <- atl_filter_bounds(condors_filter,
                                   x = "location_long",
                                   y = "location_lat",
                                   x_range = c(-120, -90),
                                   y_range = c(30, 39.5),
                                   remove_inside = FALSE)

####----(3) Visualise----####
#' download map
ggmap::register_google(key = "XXX")

#' map parameters
map <- ggmap(get_googlemap(center = c(lon = -112.02, lat = 36),
                           zoom = 5, scale = 2,
                           maptype ='hybrid', color = "color"))
map

#' can change zoom parameters on map to view it better (zoom 5)
pdf("map1.pdf", height=3.5, width=4)
map1 <- map + geom_path(aes(location_long,location_lat, col = id, group=id), 
                        data = condors_atlas) + theme(legend.position = "none")+
  scale_colour_viridis_d(option = "plasma", alpha = 0.8, begin = 0.25)+
  labs(x = "Longitude", y ="Latitude")+
  ggsn::scalebar(x.min = -125, x.max = -120,
                 y.min = 25,  y.max = 25.4,
                 location = "bottomleft", dist = 200,
                 dist_unit = "km", transform = TRUE, 
                 st.bottom = FALSE, height = 1,
                 st.dist = 2, st.size = 1.75,
                 border.size = 0.5,
                 st.color = "white",
                 box.color = "white")

ggsn::north2(map1, x = 0.20, y = 0.90, scale = 0.09, symbol = 4)
dev.off()

####----(4) Finer-scale outliers----####
#' can still see that there are these 'spikes' in movement
#' to deal with this, we'll be creating additional columns:
#' 1. speed in 
#' 2. speed out 
#' 3. angle 

#' check for angle and speed
xx <- condors_atlas %>% tidyr::nest(data = -id)

condors_atlas_filtered <- xx %>% 
  mutate(data = purrr::map(data, 
                    ~ mutate(.x, 
                             speed_in = atl_get_speed(.x,
                                                      x = "location_long",
                                                      y = "location_lat",
                                                      time = "timestamps",
                                                      type = c("in")),
                             speed_out = atl_get_speed(.x,
                                                       x = "location_long",
                                                       y = "location_lat",
                                                       time = "timestamps",
                                                       type = c("out")),
                             angle = atl_turning_angle(.x,
                                                       x = "location_long",
                                                       y = "location_lat",
                                                       time = "timestamps")) 
  ))

condors_atlas_filtered <- tidyr::unnest(condors_atlas_filtered)

####----(5) Turning angle-speed plot----####
plot1 <- condors_atlas_filtered %>% ggplot() + 
  geom_point(aes(x = speed_in, y = speed_out, fill = angle),  alpha = 0.75, pch = 21)+
  scale_fill_viridis_c()+
  theme_minimal()+
  ylab("Speed in (m/s)") + xlab("Speed out (m/s)")+
  geom_hline(yintercept=0.8, linetype="dashed", color = "darkred")+
  geom_vline(xintercept=0.8, linetype="dashed", color = "darkred")+
  annotate(geom = "rect", xmin=0.8, xmax=Inf, ymin=0.8, ymax=Inf, 
           alpha=0.2, fill="darkred")

gridExtra::grid.arrange(egg::set_panel_size(p=plot1, width=unit(10, "cm"), height=unit(8, "cm")))
ggsave("plot1.pdf", dev='pdf', height=8, width=10, units="cm", dpi = 500)

####----(6) Filter by turning angle and speed----####
#' there are some points where the speed in and turning angles are high 
#' (unrealistic for the birds to make such sharp turns at high speed)
summary(condors_atlas_filtered$speed_in)
summary(condors_atlas_filtered$speed_out)
summary(condors_atlas_filtered$angle)

S <- 0.8 #' conservative threshold for high speed from plot
A <- 150 #' conservative threshold for high turning (3rd quartile)

condors_atlas_final <- atl_filter_covariates(data = condors_atlas_filtered,
                                             filters = c(
                                               "(speed_in < S & speed_out < S | angle < A)"
                                             ))

condors_atlas_final <- condors_atlas_final %>% filter(speed_in < 16 & speed_out < 16)

plot2 <- condors_atlas_final %>% ggplot() + 
  geom_point(aes(x = speed_in, y = speed_out, fill = angle),  alpha = 0.75, pch = 21)+
  scale_fill_viridis_c()+
  theme_minimal()+
  ylab("Speed in (m/s)") + xlab("Speed out (m/s)")+
  geom_hline(yintercept=0.8, linetype="dashed", color = "darkred")+
  geom_vline(xintercept=0.8, linetype="dashed", color = "darkred")+
  annotate(geom = "rect", xmin=0.8, xmax=Inf, ymin=0.8, ymax=Inf, 
           alpha=0.2, fill="darkred")

gridExtra::grid.arrange(egg::set_panel_size(p=plot1, width=unit(10, "cm"), height=unit(8, "cm")))
ggsave("plot2.pdf", dev='pdf', height=8, width=10, units="cm", dpi = 500)

####----(7) Visualise----####
pdf("map2.pdf", height=3.5, width=4)
map2 <- map + geom_path(aes(location_long,location_lat, col = id, group=id), 
                        data = condors_atlas_final) + theme(legend.position = "none")+
  scale_colour_viridis_d(option = "plasma", alpha = 0.8, begin = 0.25)+
  labs(x = "Longitude", y ="Latitude")+
  ggsn::scalebar(x.min = -125, x.max = -120,
                 y.min = 25,  y.max = 25.4,
                 location = "bottomleft", dist = 200,
                 dist_unit = "km", transform = TRUE, 
                 st.bottom = FALSE, height = 1,
                 st.dist = 2, st.size = 1.75,
                 border.size = 0.5,
                 st.color = "white",
                 box.color = "white")

ggsn::north2(map2, x = 0.20, y = 0.90, scale = 0.09, symbol = 4)
dev.off()

range(condors_atlas_final$location_lat)
range(condors_atlas_final$location_long)

#' Hours of day
condors_atlas_final$time <- format(condors_atlas_final$timestamps, format = "%H", tz = "US/Arizona")

####----DAYTIME MOVEMENT PATTERNS----####
####----(1) Set sunrise and sunset----####
#' Create a dataframe with timestamps and location for 'suncalc' package
sun <- data.frame(date= as.Date(condors_atlas_final$timestamps, tz = "US/Arizona"), 
                  lat=condors_atlas_final$location_lat, lon=condors_atlas_final$location_long,
                  id = condors_atlas_final$id)

#' Caluclate sunrise and sunset times
sunrise <-suncalc::getSunlightTimes(data=sun, keep="dawn", tz = "US/Arizona") 
sunset <- suncalc::getSunlightTimes(data=sun, keep="sunset", tz = "US/Arizona") 

setdiff(condors_atlas_final$location_lat, sunrise$lat)

condors_atlas_final$sunrise <- sunrise$dawn
condors_atlas_final$sunset <- sunset$sunset

condors_atlas_final$night_day <- ifelse(condors_atlas_final$timestamps < condors_atlas_final$sunrise | 
                                          condors_atlas_final$timestamps > condors_atlas_final$sunset,
                                         1,0)

#' check if any night time points
condors_atlas_final %>% group_by(night_day) %>% tally()

####----(2) Filter out night time points----####
#' selecting for daytime points
condors_filter_final <- condors_atlas_final %>% filter(night_day == 0)
hist(as.numeric(condors_filter_final$time))

#' remove dfs
rm(sun, sunrise, sunset)

####----READ IN BLOOD LEAD TESTING DATA ----####
####----(1) Read data up to 2020----####
bleedtag <- read.csv('./Final input datasets/Condor bleed and tag.csv')

#' remove 'Rx' rows - restests of birds in captivity (only want field data)
bleedtag_select <- bleedtag %>% filter(!grepl('Rx', X))

#' keep certain columns
col_keep <- c("Date.of.Event", "Bird..","RETEST", "ESA.1................ug.dl", "Notes")
bleedtag_select <- bleedtag_select[,col_keep]

#'rename columns
names(bleedtag_select) <- c("date", "id", "retest", "lead", "note")
head(bleedtag_select)

####----(2) Read in 2021-2022 data----####
bleedtag2 <- read.csv('./Final input datasets/Trapping Data 2021-Current.csv')
head(bleedtag2)

#' combine the datasets
bleedtag_filter <- rbind(bleedtag_select, bleedtag2)

#'change format of date column 
class(bleedtag_filter$date)
bleedtag_filter$date <- as.Date(bleedtag_filter$date, format = "%d-%b-%y")
bleedtag_filter <- bleedtag_filter[!is.na(bleedtag_filter$date),]
head(bleedtag_filter)

#' create year column
bleedtag_filter$year <- format(bleedtag_filter$date, "%Y")

####----(3) Setting lead values and categories----####
#' clean lead column 
class(bleedtag_filter$lead)
levels(as.factor(bleedtag_filter$lead))

#' get rid of character
levels(as.factor(bleedtag_filter$lead))

bleedtag_filter$lead <- plyr::revalue(bleedtag_filter$lead, c(">65"="65",
                                                              ">65.0" ="65",
                                                              ">65.00"="65",
                                                              " >65" = "65",
                                                              "<65" = "64", 
                                                              "109" = "65",
                                                              "HI"="65"))
bleedtag_filter$lead <- plyr::revalue(bleedtag_filter$lead, c("low"="14",
                                                        "Low"="14",
                                                        "LO" = "14",
                                                        "\"low\"" = "14",
                                                        "18,4" = "18.4",
                                                        "13,2" = "13.2"))

levels(as.factor(bleedtag_filter$lead))

#' make lead numeric
bleedtag_filter$lead_num <- as.numeric(bleedtag_filter$lead)
levels(as.factor(bleedtag_filter$lead_num))

#' get rid for more na rows
bleedtag_filter <- bleedtag_filter[!is.na(bleedtag_filter$lead),]
levels(as.factor(bleedtag_filter$lead_num))

#“background” lead concentrations of 0-14 µg/dl. 
#levels of 15-29 µg/dl (indicating lead exposure), 
#31-59 µg/dl, 
#60 µg/dl, the threshold at which the term “clinically affected” has been (Fry and Maurer 2003 & Cade 2007)
bleedtag_filter$lead_level <- ifelse(bleedtag_filter$lead_num < 15, "1-Background",
                                     ifelse(bleedtag_filter$lead_num >=15 & bleedtag_filter$lead_num < 30, "2-Exposure",
                                            ifelse(bleedtag_filter$lead_num >= 30 & bleedtag_filter$lead_num <60, "3-High Exposure",
                                                   "4-Clinically affected")))

table(is.na(bleedtag_filter$lead_num))
bleedtag_filter <- bleedtag_filter %>% na.omit()
bleedtag_filter %>% ggplot(aes(x = lead_level, y = (lead_num)))+geom_boxplot()

####----ASSIGNING EXPOSURES----####
####----(1) Defining retests----####
xx <- bleedtag_filter %>% group_by(id, date, year, lead_num, retest, lead_level, note)%>% 
  dplyr::summarise(n = n())

# Sort the dataset by individual, date, and condor year; compute difference in dates
xx <- xx %>%
  arrange(id, date) %>%
  group_by(id) %>% #condor year groups obs by the field season
  mutate(diff = date - lag(date))

xx$diff <- as.numeric(xx$diff)

#' four_wk column: if first diff, second diff or third diff are greater than 
#' 30, then greater than 4 weeks otherwise, less than 4 weeks

xx <- xx %>%
  mutate(
    new_diff = case_when(
      is.na(diff) ~ "na",  # if diff is NA and new_diff is not NA
      diff > 30 ~ "na",                        # if diff is greater than 60
      TRUE ~ ""
    )
  ) %>%
  ungroup()

counter <- 1

xx <- xx %>% mutate(number = ifelse(new_diff == "na", counter, NA))%>% 
  tidyr::fill(number, .direction = 'down') %>% 
  mutate(number = ifelse(new_diff == 'na', counter, number))

#' Increment the counter for each 'no' value
xx$counter <- cumsum(xx$new_diff == 'na')

#' Calculate the difference again between consecutive dates within each individual
#' as some individuals have small intervals in when they were retested
xx <- xx %>%
  group_by(id, counter) %>%
  mutate(first_diff = date - first(date), #difference to first date [grouped by condor year]
         second_diff = date - nth(date, 2), #second date
         third_diff = date - nth(date, 3), #third date 
         fourth_diff = date - nth(date, 4), #fourth date 
         #' create another count column for setting retest date differences (no. days)
         count = cumsum(as.numeric(first_diff) > 30 | lag(as.numeric(first_diff) > 30, default = FALSE))) 


xx <- xx %>%
  mutate(
    four_wk = case_when(
      count == 0 & as.numeric(first_diff) == 0 ~ "n",
      count == 1 & as.numeric(first_diff) > 30 ~ "n",
      count == 2 & as.numeric(second_diff) > 30 ~ "n",
      count == 3 & as.numeric(third_diff) > 30 ~ "n",
      count == 4 & as.numeric(fourth_diff) > 30 ~ "n",
      TRUE ~ "y"
    )
  ) %>%
  ungroup()
#' y = within 4 weeks

#' lead level difference between testing dates
xx <- xx %>% group_by(id) %>% mutate(lead_diff = lead_num - lag(lead_num))

#' Filter dates to match condor range 
range(condors_filter_final$timestamps)
range(xx$date)
retest_final <- xx %>% filter(date >= as.Date("2015-01-15") &
                                    date <= as.Date("2023-01-29"))

retest_final %>% group_by(year)%>% 
  dplyr::summarise(count = n()) %>% print(n = 24)

retest_final %>% ggplot(aes(x = lead_level, y = lead_num))+
  geom_boxplot()

####----30 DAYS PRIOR TO TESTING----####
####----(1) Selecting 30 days prior to testing, without overlaps in test date intervals---####
#' ensure that merging columns are of same type
class(condors_filter_final$id)
class(retest_final$id)
retest_final$id <- as.factor(retest_final$id)
class(condors_filter_final$date)
class(retest_final$date)

#' Merge datasets 
levels(condors_filter_final$id)
levels(retest_final$id)

condor_test <- merge(condors_filter_final, retest_final, by = c("id", "date"),
                     all.x = TRUE)

levels(as.factor(condor_test$id))

#'selecting columns
names(condor_test)

xx <- distinct(as.data.frame(condor_test[,c("id", "date", "lead_num", "lead_level",
                                                          "note",
                                            "four_wk", "year")]))
xx <- xx[!is.na(xx$lead_num),]

#' 30 days prior to lead testing
xx$startdate_four <- xx$date - 30

xx$enddate <- xx$date
xx <- droplevels(xx)

#' filtering out tests with retests that are within a 30d period
xx <- xx %>% filter(four_wk =="n") 

#' Checking the day difference between each testing period per indv
xx <- xx %>%
  group_by(id) %>%
  mutate(diff = date - lag(date, default = first(date)))

xx <- xx %>% filter(diff <= 0 | diff >= 30 )


####----(2) Filtering the 4 weeks prior from GPS data----####
#' Filter condor_merge according to start and end date in condor_test_date
#' Create an empty vector to store results in
four_wk<- list()

for (i in 1:nrow(xx)){
  four_wk[[i]] <- (condors_filter_final %>%
                     filter(id %in% xx$id[i]) %>%
                     filter(between(date, xx$startdate_four[i], 
                                    xx$enddate[i])))
  four_wk[[i]]$lead_num <- xx$lead_num[i]
  four_wk[[i]]$test_date <- xx$date[i]
  four_wk[[i]]$lead_level <- xx$lead_level[i]
  four_wk[[i]]$four_wk <- xx$four_wk[i]
}

four_wk <- do.call("rbind",four_wk)

####----(3) Check if indv were treated in any of the 30 days----####
#' treatment of birds in captivity
bleedtag_treatment <- bleedtag %>% filter(grepl('Rx', X))

#' keep certain columns
col_keep <- c("Date.of.Event", "Bird..","RETEST", "ESA.1................ug.dl", "Notes")
bleedtag_treatment <- bleedtag_treatment[,col_keep]

#'rename columns
names(bleedtag_treatment) <- c("date", "id", "retest", "lead", "note")
bleedtag_treatment$date <- as.Date(bleedtag_treatment$date, format = "%d-%b-%y")

bleedtag_treatment$id_date <- paste0(bleedtag_treatment$id, "-", bleedtag_treatment$date)

xx <- four_wk
xx$id_date <- paste0(xx$id, "-", xx$date)
xx <- xx %>% group_by(id_date, id, test_date, four_wk) %>% tally()
treatment <- merge(xx, bleedtag_treatment, by = c("id_date", "id"), all.x = TRUE)
treatment$id_test <- paste0(treatment$id, "-", treatment$test_date)

x <- treatment %>% tidyr::drop_na() %>% group_by(id_test) %>% tally()

treatment <- treatment %>% 
  filter(!(id_test %in% x$id_test)) %>% droplevels() %>%
  select("id_date", "id", "test_date", "four_wk", "n", "id_test")

levels(as.factor(treatment$id_test))

x <- which(grepl(c("barn"), bleedtag2$note))
bleedtag2$date <- as.Date(bleedtag2$date, format = "%d-%b-%y")

x <-bleedtag2 %>%
  slice(c(x, x + 1)) %>%
  group_by(id) %>%
  mutate(diff = as.numeric(date - lag(date))) %>%
  filter(diff < 30)
#' none of the testing intervals less than 30 days are in the dataset from bleedtag2

####----(4) Filter out tests----####
#' number of fixes for each indv in four week period
four_wk$id_test <- paste0(four_wk$id, "-", four_wk$test_date)
four_wk <- four_wk %>% filter(id_test %in% treatment$id_test) 

#' check times of day
four_wk %>% ggplot(aes(time))+geom_histogram(stat = "count")

#' Proportion of tests revealing lead level exposures by year
plot3 <- four_wk %>% #' can also plot with retest 
  mutate(year = format(test_date, "%Y"),
         month = format(test_date, "%m")) %>%
  group_by(month, year, lead_level, id_test)%>% 
  dplyr::summarise(count = n()) %>% 
  group_by(month, year, lead_level)%>% 
  dplyr::summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) %>%
  ggplot(aes(x = factor(month), y = perc*100, fill = factor(lead_level)))+
  geom_bar(stat="identity", position = "fill") +
  labs(x = "Month", y = "Proportion of tests", fill = "Lead level") + 
  geom_text(aes(label=paste0(count)),
            position=position_fill(vjust = 0.5),
            size = 1.7)+
  facet_grid(rows = vars(year))+
  theme_classic()+
  scale_fill_brewer(palette="Reds")+
  theme(axis.text.x = element_text(angle = 60, hjust =1))
 
gridExtra::grid.arrange(egg::set_panel_size(p=plot3, width=unit(18, "cm"), height=unit(16, "cm")))
ggsave("plot3.pdf", dev='pdf', height=16, width=18, units="cm", dpi = 500)

####----(3) Visualize----####
#plot 
map <- ggmap(get_googlemap(center = c(lon = -112.5, lat = 36.5), 
                           zoom = 8, scale = 2,
                           maptype ='hybrid'))

pdf("map3.pdf", height=6, width=6)
map3 <- map + geom_path(aes(location_long,location_lat, col = id, group=id), alpha = 1, linewidth = 0.4, 
                      data = four_wk)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = "Longitude", y = "Latitude")+
  scale_colour_viridis_d(alpha = 0.8, begin = 0.25)+
  ggsn::scalebar(x.min = -114, x.max = -113,
                 y.min = 35.3,  y.max = 35.4,
                 location = "bottomleft", dist = 20,
                 dist_unit = "km", transform = TRUE, 
                 st.bottom = FALSE, height = 1,
                 st.dist = 1, st.size = 1.75,
                 border.size = 0.5,
                 st.color = "black",
                 box.color = "white")+
  ggforce::facet_wrap_paginate(.~lead_level, ncol = 2, nrow = 2, page = 1)

ggsn::north2(map3, x = 0.14, y = 0.90, scale = 0.05, symbol = 4)
dev.off()

####----DEMOGRAPHIC DATA----####
####----(1) Read in demographic data----#### 
demog <- read.csv("./Final input datasets/Condor demog.csv")
names(demog)
levels(as.factor(demog$Studbook))
demog$Studbook <- as.factor(demog$Studbook)
levels(as.factor(four_wk$id))

#' merge with demographic data
demog_final <- merge(four_wk, demog, by.x = c("id"),
                     by.y = c("Studbook"),
                     all.x = TRUE)

####----(2) Assigning age----####
class(demog_final$date)
class(demog_final$Hatch.Date)
demog_final$Hatch.Date <- as.Date(demog_final$Hatch.Date, format = "%d/%m/%Y")

demog_final$age <- demog_final$date - demog_final$Hatch.Date
hist(as.numeric(demog_final$age))
demog_final$age_class <- ifelse(demog_final$age <= 365, "juvenile",
                                     ifelse(demog_final$age > 365 & demog_final$age <= 1460, "sub-adult", "adult"))

####----(3) Variable class----####
demog_final <- demog_final %>% 
  mutate(age = as.numeric(age),
         age_class = as.factor(age_class),
         Sex = as.factor(Sex))

####----(4) Sex unknown removed, NA removed----####
demog_final <- demog_final %>% filter(Sex != "Unknown") %>%
  filter(!is.na(Hatch.Date))

####----DATA TO PROJECTION----####
####----(1) Interpolate missing time intervals----####
#' Per individual and date
demog_final$id_date <- paste0(demog_final$id,"-",demog_final$date)

#' export demog final
write.csv(demog_final, "demog_final.csv")

nested <- demog_final %>% 
  droplevels() %>%
  tidyr::nest(data = c(-id_date)) %>%
  mutate(rows = as.numeric(map(data, nrow)))

#' remove individuals with less than 3 tracking points a day
nested <- nested %>% filter(rows > 2)

nested <- nested %>%
  mutate(data_move = map(data,  ~move::move(.x, x = .x$location_long, y= .x$location_lat, t = .x$timestamps)),
         data_intrpolt = map(data_move,  
                             ~move::interpolateTime(.x, time=as.difftime(60, units="mins"), spaceMethod='greatcircle')),
         data_final = map(data_intrpolt, as.data.frame))

#' unnest
intrpl_track <- nested %>%
  rowwise() %>%
  unnest(c(data_final))

####----(2) Convert interpolated data to projection----####
#' Duplicate the df
names(intrpl_track)
xx <- intrpl_track
coordinates(xx) <-  c("coords.x1", "coords.x2")
proj4string(xx) <- CRS("+proj=longlat +datum=WGS84")
xx <- spTransform(xx, CRS("EPSG:26712"))
#' for more info, see https://epsg.io/26712
condors_proj <- as.data.frame(xx)

#' convert to km
condors_proj$coords.x1.km <- condors_proj$coords.x1/1000
condors_proj$coords.x2.km <- condors_proj$coords.x2/1000

####----(3) Make movement track---####
condors_track <- make_track(condors_proj, .x = coords.x1.km, 
                               .y = coords.x2.km, .t = timestamps.1,
                               id = id_date,
                               crs = 26712)

#' Save the class of the track
xx <- class(condors_track)
condors_track$date <- as.Date(condors_track$t_, tz = "US/Arizona")
names(condors_track)

####----COMPUTE SL AND MSD BY DAY----####
####----(1) Calculate the mean_sl and msd per indv. for each day---####
#' nest track per id 
nested_track <- condors_track %>% 
  droplevels() %>%
  tidyr::nest(data = c(-id))

#' Attribute class 'track'
class(nested_track)<-xx

nested_track <- nested_track %>% 
  mutate(sl = map(data, step_lengths),
         msd = map(data, msd),
         nsd = map(data, nsd)) %>% 
  unnest(cols = c(data,sl, msd, nsd))

#' Note that there will be an NA at the last timestamp at the end of the day
head(nested_track)

####----DAILY MCP----####
####----(1) Track by id----####
#' Split by id 
xx <- intrpl_track
coordinates(xx) <-  c("coords.x1", "coords.x2")
proj4string(xx) <- CRS("+proj=longlat +datum=WGS84")
xx <- spTransform(xx, CRS("+proj=longlat +datum=WGS84"))

xx <- split(xx, xx$id_date)

#' For each date, create track for each individual
for (i in seq_along(xx)){
  xx[[i]] <- make_track(xx[[i]],
                                       .x = coords.x1, 
                                       .y = coords.x2, 
                                       .t = timestamps.1,
                                       id = id,
                                       crs = "+proj=longlat +datum=WGS84")
}


#' Create list for MCP
mcp_ld = rapply(xx ,length, how = "list")

sf::sf_use_s2(FALSE) # for geoms to 'connect' well (issue with duplicated vertices)

#' using the track data, create daily mcp ranges with the specific bandwidths
for (j in seq_along(xx)){
    mcp_ld[[j]] <- hr_mcp(xx[[j]], levels = c(0.95))
}

####----(2) Extract area----####
mcp_ld_area <- mcp_ld
for (m in seq_along(mcp_ld)){
    mcp_ld_area[[m]] <- hr_area(mcp_ld[[m]])
    mcp_ld_area[[m]] <- as.numeric(mcp_ld_area[[m]]$area)
  } # area is in m2

mcp_ld_area <- purrr::map_dfr(mcp_ld_area, function(x){
  unlist(x)
})

mcp_ld_area <- reshape2::melt(mcp_ld_area,
                                       na.rm = TRUE,
                                       variable.name = "id",
                                       value.name = "mcp")

####----ALL MOVEMENT METRICS----####
####----(1) Final datasets----####
metrics_final <- nested_track %>% 
  group_by(id, date, msd)  %>% 
  summarise(m.sl = mean(na.exclude(sl)),
            m.nsd = mean(nsd),
            n_points = n())

metrics_final <- merge(metrics_final, mcp_ld_area,
                  by = c("id"))

####----(2) Combining with other desired columns----####
head(demog_final)
xx <- demog_final %>% select(c("id", "lead_num", "lead_level", "test_date", 
                                "Sex", "age", "age_class", "id_date")) %>%
  unique()


metrics_final <- merge(metrics_final, xx, by.x = "id", by.y = "id_date", all.x = TRUE)
metrics_final$mcp.km <- metrics_final$mcp/1000000

####----(3) Days prior----####
#' Add a column for day prior
metrics_final$day_prior = metrics_final$test_date - metrics_final$date

#' remove day of testing 
metrics_final <- metrics_final %>% filter(day_prior != 0)

colnames(metrics_final)[colnames(metrics_final) == 'id.y'] <- 'raw_id'

metrics_final$id_test <- paste0(metrics_final$raw_id,"-", metrics_final$test_date)

#' remove duplicates
unique(as.factor(metrics_final$id_test)) # n = 156

str(metrics_final)

metrics_final <- metrics_final %>% filter(Sex != "Unknown") %>% droplevels()

####----(4) Export as csv: movement final----####
write.csv(metrics_final, "metrics_final.csv")

# extract info about testing birds
xx = metrics_final %>% group_by(id_test, test_date, lead_level) %>% tally() %>%
  mutate(id = substr(id_test, 1, nchar(id_test)-11))
retest_final = retest_final %>% mutate(test_date = date)
test_status = left_join(xx, retest_final, by = c("id", "lead_level", "test_date"))

write.csv(test_status, "test_status.csv")

#save.image('Condor_final.RData')
#load('Condor_final.RData')








