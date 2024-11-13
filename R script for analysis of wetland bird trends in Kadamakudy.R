#1. Getting started####

#> Set working directory. This is important for selecting input files and saving
#> outputs.In RStudio, you can simply click on the tab 'Session' on top and select 
#> 'Set Working Directory' and then 'Choose Directory' OR simply press Ctrl+Shift+H
#> Code for the doing the same:
#> setwd(filepath)
#> Filepath will look something like: "/home/paul/project/cnhs" in Linux or
#> "D:\\paul\\projects\\cnhs" in Windows.

setwd("/home/projects/cnhs")

#check if the intended directory/folder is correct:
getwd()

#> This coding project will makee use the documentation provided by Strimas-Mackey 
#> et al (2023) as a rough guide for some of the steps. 
#> https://ebird.github.io/ebird-best-practices/intro.html. 
#> So, download the bundled packages contained within ebird-best-practices in GitHub. 

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("ebird/ebird-best-practices")

#> The above code asks to avoid installing the "remotes" package (used for installing
#> packages from remote repositories including GitHub) and only installing if not
#> present in the R packages directory (different from the working directory in our
#> case).Use ?remotes or ? followed by a package name to learn about the package in
#> the help pane of R Studio. The installation of all the packages may time some time.
#> The ebirdbestpractices package will install the following (see details here:
#> https://github.com/ebird/ebird-best-practices/blob/main/DESCRIPTION):

# arrow,
# auk (>= 0.6.0),
# dplyr (>= 1.0.0),
# ebirdst (>= 3.2022.2),
# exactextractr,
# fields,
# ggplot2 (>= 3.4.0),
# gridExtra,
# knitr,
# landscapemetrics,
# lubridate,
# mccf1,
# precrec,
# PresenceAbsence,
# purrr,
# ranger (>= 0.11.2),
# readr,
# rlang (>= 0.4.1),
# rnaturalearth,
# scales,
# scam,
# sf (>= 1.0-0),
# smoothr,
# stringr,
# terra (>= 1.6-3),
# tibble,
# tidyr (>= 1.2.1),
# units,
# viridis

#> You won't need all of them. So, you can install only those that you require for
#> this excercise indvidually too if you want to save some space.
#> It is likely that a number of times, the installation of several packages will fail.
#> So, follow the instructions given. For example, if a certain dependenncy is not
#> available for a required package, it will probably be have to be downloaded using
#> the GitHub repository. Use the following command in that case:

install.packages(ebirdst)

#and if that doesn't work:

remotes:install_github("CornellLabofOrnithology/ebirdst")

#> This comes from the url https://github.com/CornellLabofOrnithology/ebirdst
#> Note that you most likely won't have any problem installing the "ebirdst" package
#> using install.packages(). It is just as an example.

#> Other issues can only be solved by installing packages like libudunits2-dev and
#> libgdal-dev as a root user. If you are using Linux (which you ideally should be),
#> open the terminal (ctrl+alt+t) and type the following and press enter (without
#> the hashes/#):

#>sudo apt-get install libudunits2-dev
#>sudo apt install libgdal-dev

#> To install only those packages used so far (use install.packages("packagename",
#>  dependencies = TRUE) to avoid issues regarding dependencies.
#> 
#> install.packages("auk")
#> install.packages("dplyr")
#> install.packages("osmdata")
#> install.packages("ggplot2")
#> install.packages("ggspatial")
#> install.packages("tibble")
#> install.packages("lubridate")
#> install.packages("flextable")


#2. Data preparation####

##i) Clipping EBD and SED to Polygon####

#Read the kml file for the study areas:
kk = st_read(dsn = "Kadamakudy Grama Panchayat - Boundaries.kml")
vp = st_read(dsn = "Varapuzha Grama Panchayat - Boundaries.kml")

#plot the pachayats/muncipalities (the study areas):
plot(vp)
plot(kk)

#convert to simple features in R (not to be confused with ESRI shapefiles):
kk <- st_as_sf(kk)
vp <- st_as_sf(vp)

#merge the two kmls while removing the overlapping portion:
merged <- st_union(kk, vp)
#See the structure of the file:
str(merged)
#Plot the merged file:
plot(merged)

#Save this merged vectors to the folder as a geojson. 
st_write(merged %>% st_geometry(), "merged.geojson")

#remove unwanted files from the Global Environment to save space in RAM:
remove(kk,vp)

#For clipping the EBD dataset to a bounding box ('library' loads the required package):
library(auk)
library(dplyr)

# Save the data falling within the bounding box:
f_out <- "ebd_kadamakduy_merged.txt"
auk_ebd("ebd_IN-KL-ER_smp_relSep-2024.txt") %>%
  # define filters
  auk_bbox(merged)%>%
  # compile and run filters
  auk_filter(f_out)

#> bbox is a rectangular bounding bxox which is necessary if the data is largely 
#> before the exact data is extracted using the "merged" polygon of the study area.

#> Here the ebd_IN-KL-ER_smp_relSep-2024.txt is the eBird Basic Dataset for the 
#> district Ernakulam in which Kadamakudy wetland is located. 
#> Download the EBD after making a request here:
#> https://science.ebird.org/en/use-ebird-data/download-ebird-data-products
#> The EBD used for this project contained data from January 1 to September 2024
#> for the district("county" in Ebird) Ernakulam (State = Kerala, Country = India).

# Paste the results from the Console below (use ctrl+shift+c to comment/add hashtag
# to the text)

#Input
#EBD: /home/paul/cnhs/ebd_IN-KL-ER_smp_relSep-2024.txt

#Output
#Filters not executed

#Filters
#Species: all
#Countries: all
#States: all
#Counties: all
#BCRs: all
#Bounding box: Lon 76.2 - 76.3; Lat 10 - 10.1
#Years: all
#Date: all
#Start time: all
#Last edited date: all
#Protocol: all
#Project code: all
#Duration: all
#Distance travelled: all
#Records with breeding codes only: no
#Exotic Codes: all
#Complete checklists only: no

# Importantly, it gives us the coordinates of the bounding box which will useful
# later on.

# Read the saved ebd_kadamakduy_merged.txt file from the working directory:

ebd <- read_ebd("ebd_kadamakduy_merged.txt")

#> The read_ebd does three things under the hood while reading the data:
#> Cleans up the label for use in R, taxonomic rollup (subspecies/intergrades etc
#> are rolled up to species level and individuals not identified to species levels 
#> are dropped), and only one checklist from the shared checklist
#> is kept to avoid pseudoreplicates. Read about the exact details here:
#> https://ebird.github.io/ebird-best-practices/ebird.html#sec-ebird-import

#subsetting to polygon

# convert to sf object
ebd_sf <- ebd %>%
  dplyr::select(longitude, latitude) %>%
  st_as_sf( coords = c("longitude", "latitude"), crs = 4326)

# put polygons in same crs
poly_ll <- st_transform(merged, crs = st_crs(ebd_sf))

# identify points in polygon
in_poly <- st_within(ebd_sf, poly_ll, sparse = FALSE)

# subset data frame
ebd_in_poly <- ebd[in_poly[, 1], ]

# look at the structure of ebd
str(ebd_in_poly)

# save the subsetted data:
write.table(ebd_in_poly,"ebd_in_poly.txt",sep="\t",row.names=FALSE)

#Do the same for the Sampling Event Data (SED)

f_sed <- "ebd_IN-KL-ER_smp_relSep-2024_sampling.txt"
checklists <- read_sampling(f_sed)
glimpse(checklists)

#clipping to study area polygon
f_out_sed <- "sed_kadamakduy_merged.txt"
auk_sampling("ebd_IN-KL-ER_smp_relSep-2024_sampling.txt") %>%
  # define filters
  auk_bbox(merged)%>%
  # compile and run filters
  auk_filter(f_out_sed)

# Input
# Sampling events: /home/paul/cnhs/ebd_IN-KL-ER_smp_relSep-2024_sampling.txt
#
# Output
# Sampling events: /home/paul/cnhs/sed_kadamakduy_merged.txt
#
# Filters
# Countries: all
# States: all
# Counties: all
# Bounding box: Lon 76.2 - 76.3; Lat 10 - 10.1
# Years: all
# Date: all
# Start time: all
# Last edited date: all
# Protocol: all
# Project code: all
# Duration: all
# Distance travelled: all
# Complete checklists only: no

sed <- read_sampling("sed_kadamakduy_merged.txt")

remove(kk,vp,g)

#subsetting to polygon

# convert to shapefile (sf) object
sed_sf <- sed %>%
  dplyr::select(longitude, latitude) %>%
  st_as_sf( coords = c("longitude", "latitude"), crs = 4326)

# put polygons in same crs
poly_ll <- st_transform(merged, crs = st_crs(sed_sf))

# identify points in polygon
in_poly <- st_within(sed_sf, poly_ll, sparse = FALSE)

# subset data frame
sed_in_poly <- sed[in_poly[, 1], ]

#Have a glimpse at the clipped data:
glimpse(sed_in_poly)

#Save to folder
write.table(sed_in_poly,"sed_in_poly.txt",sep="\t",row.names=FALSE)


## ii) Filtering for analysis####

# Other types of filtering needs to be done before it can be used for analysis:

# remove observations without matching checklists (due to editing of shared checklists)
observations <- semi_join(ebd_in_poly, sed_in_poly, by = "checklist_id")

# Filter out only complete observations:

# filter the checklist data
checklists <- sed_in_poly |>
  filter(all_species_reported)

# filter the observation data (do this before previous step - ideally)
observations <- observations |>
  filter(all_species_reported)

# Zero-filling
zf <- auk_zerofill(observations, checklists, collapse = TRUE)

# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

str(zf)

# clean up and add more custom variables
zf <- zf |>
  mutate(
    # convert count to integer and X to NA
    # ignore the warning "NAs introduced by coercion"
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for stationary counts
    effort_distance_km = if_else(protocol_type == "Stationary",
                                 0, effort_distance_km),
    # convert duration to hours
    effort_hours = duration_minutes / 60,
    # speed km/h
    effort_speed_kmph = effort_distance_km / effort_hours,
    # convert time to decimal hours since midnight
    hours_of_day = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

#Accounting for variation in efforts

summarise(zf$effort_speed_kmph)

# additional filtering
zf_filtered <- zf |>
  filter(protocol_type %in% c("Stationary", "Traveling"),
         effort_hours <= 6,
         effort_distance_km <= 8, #Because anything more than 8 is likely longer
         #than the longest length of the study area
         effort_speed_kmph <= 50, #being generous here
         number_observers <= 10)

write.table(zf_filtered,"filtered_merged.txt",sep="\t",row.names=FALSE)
head(zf_filtered)
str(zf_filtered)

## iii) Open Street Maps####

#> For creating the study are map, it will look better if there is there are roads,
#> buildings etc in the backdrop to get a feel for the landscape. an Open Street 
#> Map API can be used using the "osmdata" package for this purpose. 

install.packages("osmdata")
library(osmdata)

#List the features available in OSM
available_features()
#Check sub-categories of these:
available_tags(feature = "water")

#>Manually define the latitude and longitude of the bounding box for OSM data extraction
#>Use the same as the bounding box using the bbox argument of auk (used earlier):
#> A 2x2 matrix
osm_kadamakudy <- matrix(data = c(76.2, 76.3, 10, 10.1),
                         nrow = 2,
                         byrow = TRUE)
# Update column and row names
colnames(osm_kadamakudy) <- c("min", "max")
rownames(osm_kadamakudy) <- c("x", "y")
# Print the matrix to the console
osm_kadamakudy

# Check available tags for highways (if you want something more specific than all
# roads):
available_tags(feature = "highway")

osm_k_roads <- osm_kadamakudy %>%
  opq() %>%
  add_osm_feature(key = "highway") %>% #,
  #value = c("motorway", "trunk", "primary", "secondary","tertiary",
  #         "residential","pedestrian","road","footway","unclassified",
  #        "busway","path"))
  osmdata_sf()

osm_k_roads
# Object of class 'osmdata' with:
# $bbox : 10,76.2,10.1,76.3
# $overpass_call : The call submitted to the overpass API
# $meta : metadata including timestamp and version numbers
# $osm_points : 'sf' Simple Features Collection with 37413 points
# $osm_lines : 'sf' Simple Features Collection with 5794 linestrings
# $osm_polygons : 'sf' Simple Features Collection with 7 polygons
# $osm_multilines : 'sf' Simple Features Collection with 1 multilinestrings
# $osm_multipolygons : NULL

#Plot the points:
ggplot() +
  geom_sf(data = osm_k$osm_points,
          inherit.aes = FALSE,
          color = "black",
          size = 0.2)

#Plot the lines:
ggplot() +
  geom_sf(data = osm_k$osm_lines,
          inherit.aes = FALSE,
          color = "black",
          size = 0.2)

#Similarly do the same for other classes:

#For buildings:
osm_k_building <- osm_kadamakudy %>%
  opq() %>%
  add_osm_feature(key = "building") %>%
  osmdata_sf()

osm_k_building

# Object of class 'osmdata' with:
# $bbox : 10,76.2,10.1,76.3
# $overpass_call : The call submitted to the overpass API
# $meta : metadata including timestamp and version numbers
# $osm_points : 'sf' Simple Features Collection with 88683 points
# $osm_lines : NULL
# $osm_polygons : 'sf' Simple Features Collection with 20330 polygons
# $osm_multilines : NULL
# $osm_multipolygons : 'sf' Simple Features Collection with 6 multipolygons

ggplot() +
  geom_sf(data = osm_k_building$osm_polygons,
          inherit.aes = FALSE,
          color = "blue",
          size = 0.2)

#For water:
osm_k_water <- osm_kadamakudy %>%
  opq() %>%
  add_osm_feature(key = "water") %>%
  osmdata_sf()

osm_k_water

# Object of class 'osmdata' with:
# $bbox : 10,76.2,10.1,76.3
# $overpass_call : The call submitted to the overpass API
# $meta : metadata including timestamp and version numbers
# $osm_points : 'sf' Simple Features Collection with 12135 points
# $osm_lines : 'sf' Simple Features Collection with 81 linestrings
# $osm_polygons : 'sf' Simple Features Collection with 144 polygons
# $osm_multilines : NULL
# $osm_multipolygons : 'sf' Simple Features Collection with 7 multipolygons

ggplot() +
  geom_sf(data = osm_k_water$osm_polygons,
          inherit.aes = FALSE,
          color = "lightblue",
          fill = "lightblue",
          size = 0.2)

#For landuse:
osm_k_landuse <- osm_kadamakudy %>%
  opq() %>%
  add_osm_feature(key = "landuse") %>%
  osmdata_sf()

osm_k_landuse

# Object of class 'osmdata' with:
# $bbox : 10,76.2,10.1,76.3
# $overpass_call : The call submitted to the overpass API
# $meta : metadata including timestamp and version numbers
# $osm_points : 'sf' Simple Features Collection with 16164 points
# $osm_lines : NULL
# $osm_polygons : 'sf' Simple Features Collection with 504 polygons
# $osm_multilines : NULL
# $osm_multipolygons : 'sf' Simple Features Collection with 5 multipolygons

ggplot() +
  geom_sf(data = osm_k_landuse$osm_polygons,
          inherit.aes = FALSE,
          color = "brown",
          #fill = "brown",
          size = 0.2)

#For boundary:
osm_k_boundary <- osm_kadamakudy %>%
  opq() %>%
  add_osm_feature(key = "boundary") %>%
  osmdata_sf()

osm_k_boundary

# Object of class 'osmdata' with:
#   $bbox : 10,76.2,10.1,76.3
# $overpass_call : The call submitted to the overpass API
# $meta : metadata including timestamp and version numbers
# $osm_points : 'sf' Simple Features Collection with 73048 points
# $osm_lines : 'sf' Simple Features Collection with 2150 linestrings
# $osm_polygons : 'sf' Simple Features Collection with 84 polygons
# $osm_multilines : NULL
# $osm_multipolygons : 'sf' Simple Features Collection with 123 multipolygons

ggplot() +
  geom_sf(data = osm_k_boundary$osm_points,
          inherit.aes = FALSE,
          color = "brown",
          #fill = "brown",
          size = 0.2)
#Gives Kerala + Ernakulam (?) + south India (?) boundaries - useless

#For waterway
osm_k_waterway <- osm_kadamakudy %>%
  opq() %>%
  add_osm_feature(key = "waterway") %>%
  osmdata_sf()

osm_k_waterway

# Object of class 'osmdata' with:
#   $bbox : 10,76.2,10.1,76.3
# $overpass_call : The call submitted to the overpass API
# $meta : metadata including timestamp and version numbers
# $osm_points : 'sf' Simple Features Collection with 8265 points
# $osm_lines : 'sf' Simple Features Collection with 382 linestrings
# $osm_polygons : 'sf' Simple Features Collection with 0 polygons
# $osm_multilines : 'sf' Simple Features Collection with 2 multilinestrings
# $osm_multipolygons : NULL

ggplot() +
  geom_sf(data = osm_k_waterway$osm_points,
          inherit.aes = FALSE,
          color = "darkblue",
          #fill = "brown",
          size = 0.2)
#Useless

# For tourism
osm_k_tourism <- osm_kadamakudy %>%
  opq() %>%
  add_osm_feature(key = "tourism") %>%
  osmdata_sf()

osm_k_tourism

# Object of class 'osmdata' with:
#   $bbox : 10,76.2,10.1,76.3
# $overpass_call : The call submitted to the overpass API
# $meta : metadata including timestamp and version numbers
# $osm_points : 'sf' Simple Features Collection with 215 points
# $osm_lines : NULL
# $osm_polygons : 'sf' Simple Features Collection with 34 polygons
# $osm_multilines : NULL
# $osm_multipolygons : NULL

ggplot() +
  geom_sf(data = osm_k_tourism$osm_polygons,
          inherit.aes = FALSE,
          color = "darkblue",
          #fill = "brown",
          size = 0.2)
#usless

# For route
osm_k_route <- osm_kadamakudy %>%
  opq() %>%
  add_osm_feature(key = "route") %>%
  osmdata_sf()

osm_k_route

# Object of class 'osmdata' with:
#   $bbox : 10,76.2,10.1,76.3
# $overpass_call : The call submitted to the overpass API
# $meta : metadata including timestamp and version numbers
# $osm_points : 'sf' Simple Features Collection with 76532 points
# $osm_lines : 'sf' Simple Features Collection with 5136 linestrings
# $osm_polygons : 'sf' Simple Features Collection with 44 polygons
# $osm_multilines : 'sf' Simple Features Collection with 613 multilinestrings
# $osm_multipolygons : NULL

ggplot() +
  geom_sf(data = osm_k_route$osm_lines,
          inherit.aes = FALSE,
          color = "grey20",
          #fill = "brown",
          size = 0.2)

#useless

#For place
osm_k_place <- osm_kadamakudy %>%
  opq() %>%
  add_osm_feature(key = "place") %>%
  osmdata_sf()

osm_k_place

# Object of class 'osmdata' with:
#   $bbox : 10,76.2,10.1,76.3
# $overpass_call : The call submitted to the overpass API
# $meta : metadata including timestamp and version numbers
# $osm_points : 'sf' Simple Features Collection with 74002 points
# $osm_lines : 'sf' Simple Features Collection with 722 linestrings
# $osm_polygons : 'sf' Simple Features Collection with 1407 polygons
# $osm_multilines : NULL
# $osm_multipolygons : 'sf' Simple Features Collection with 5 multipolygons

ggplot() +
  geom_sf(data = osm_k_place$osm_polygons,
          inherit.aes = FALSE,
          color = "grey20",
          fill = "brown",
          size = 0.2)

# All except "highway" is useless as their coverage is low for the area and not all
# areas are represented

remove(osm_k_boundary,osm_k_landuse,osm_k_place,osm_k_route,osm_k_tourism,
       osm_k_water,osm_k_waterway)

#3. Creating the Study Area Map####

#Plot the map using "ggplot2" package

library(ggplot2)
#install.packages(ggspatial)
library(ggspatial)

ggplot()+
  geom_sf(data = merged %>% st_geometry(), aes(fill = NULL), fill = "white",
          colour= "black")+
  #Add osm data for paths/roads
  geom_sf(data = osm_k_roads$osm_lines,
          inherit.aes = FALSE,
          color = "grey70",
          size = 0.2)+
  #Add all the checklists as darkgreen points;
  geom_sf(data = ebd_sf[in_poly[, 1], ], aes(fill = NULL), col = "darkgreen",
          pch = 19, cex = 0.5)+
  #Add annotation scale:
  ggspatial::annotation_scale(
    location = "tl",
    pad_x = unit(1.25, "cm"),
    text_pad = unit(-3, "cm"))
  #Add north arrow
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20",
      text_family = "ArcherPro Book"
    )
  )+
    #Define the coordinate limits for the maps:
  coord_sf(ylim = c(10.02, 10.10), # Crop out southern part of map
           xlim = c(76.235,76.295), # Crop out western and eastern part of map
           expand = FALSE)


#4. Creating graphs from the data####
  
  str(sed_in_poly)
  class(sed_in_poly)
  attributes(sed_in_poly)
  
# To gain information about any specific variable from the SED:
#Checking out the effort distance:
  head(sed_in_poly$effort_distance_km,100)
  min(sed_in_poly$effort_distance_km)
  
  distance <- data.frame(sed_in_poly$effort_distance_km)
  distance <- na.omit(distance)
  nrow(distance)
  head(distance)
  max(distance)
  class(distance)
  write.csv(distance, "distance.csv")
  
  #count the number of values in between 0 and 1:
  table(cut(distance$sed_in_poly.effort_distance_km, breaks = c(0, 1),
            include.lowest = TRUE))
  #count the number of values in between 1 and 2:
  table(cut(distance$sed_in_poly.effort_distance_km, breaks = c(1, 2),
            include.lowest = TRUE))

#No. of checklists by distance traveled
  
ggplot(distance, aes(x = sed_in_poly.effort_distance_km)) +
    stat_bin(breaks = seq(from = 0, to = 42, by = 1)) +
    labs(x = "Distance traveled (in km)", y = "Number of checklists",
         title = "Number of checklists according to distance travelled") +
    scale_x_continuous(breaks = seq(from = 0, to = 42, by = 2)) +
    scale_y_continuous(breaks = seq(from = 0, to = 900, by = 50))+
    theme(axis.text = element_text(colour = "black", size = rel(0.8))) +
    theme(axis.title.y = element_text(size = rel(1), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1), angle = 00)) +
    theme(panel.grid.major = element_line(colour = "green", linewidth = rel(0.5)))+
    theme(panel.grid.minor = element_line(colour = "lightgreen", linewidth = rel(0.5)))
  
#No. of checklists by Protocol Type
  
  ggplot(sed_in_poly) +
    aes(x = protocol_type) +
    geom_bar(stat = "count", fill = "blue") +
    labs(x = "Protocol type", y = "Number of checklists",
         title = "Number of checklists according to protocol type") +
    scale_y_continuous(breaks = seq(from = 0, to = 1750, by = 250))+
    theme(axis.text = element_text(colour = "black", size = rel(1))) +
    theme(axis.title.y = element_text(size = rel(1.2), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
  
#count the numbers in each protocol types:
  table(sed_in_poly$protocol_type)
  
#No. of checklists by checklist duration
  
  time <- data.frame(sed_in_poly$duration_minutes)
  str(time)
  time <- na.omit(time)
  str(time)
  
#count the number of checklists with duration in between 0 and 60 minutes:
  table(cut(time$sed_in_poly.duration_minutes, breaks = c(0, 60),
            include.lowest = TRUE))
  
  #count the number of values in between 60 and 120 minutes:
  table(cut(time$sed_in_poly.duration_minutes, breaks = c(60, 120),
            include.lowest = TRUE))
  
  #count the number of values in between 120 and 180 minutes:
  table(cut(time$sed_in_poly.duration_minutes, breaks = c(120, 180),
            include.lowest = TRUE))
  
  #Full breakup
  table(time$sed_in_poly.duration_minutes)
  
#Plot

ggplot(sed_in_poly, aes(x = duration_minutes)) +
    stat_bin(breaks = seq(from = 0, to = 720, by = 30)) +
    labs(x = "Checklist duration (in minutes)", y = "Number of checklists",
         title = "Number of checklists according to duration of birding") +
    scale_x_continuous(breaks = seq(from = 0, to = 720, by = 30)) +
    scale_y_continuous(breaks = seq(from = 0, to = 1700, by = 100))+
    theme(axis.text = element_text(colour = "black", size = rel(0.8))) +
    theme(axis.title.y = element_text(size = rel(1), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1), angle = 00)) +
    theme(panel.grid.major = element_line(colour = "darkorange", linewidth = rel(0.5)))+
    theme(panel.grid.minor = element_line(colour = "orange", linewidth = rel(0.5)))
  
# The warning "Removed 88 rows containing non-finite outside the scale range (`stat_bin()`)."
# is from NAs
  
#Derive full month names from date
str(sed_in_poly)
#Date is given as "observation_date"
  month <- data.frame(as.factor(format(sed_in_poly$observation_date, "%B")))
  year <- data.frame(format(sed_in_poly$observation_date, "%Y"))
  head(month)
  head(year)
  str(year)
  #Combine both the month and year columns:
  my <- cbind(month,year)
  head(my)
  #Change the column names:
  colnames(my) <- c("month", "year")
  head(my)
  #Convert to data frame:
  my_table <- as.data.frame.matrix(table(my))
  my_table
  #For checking the maximum value in the table:
  max(table(my))
  
# Add 'Month' as label for the first column since flextable doesn't have the feature to show
# row names (using dplyr and tibble):
# install.library
  library(tibble)
  my_table <- my_table %>% rownames_to_column("Month")
  #Sort my_table by month
  library(lubridate)
  #Function for sorting dates by month (by Tyler Rinker
  #https://stackoverflow.com/questions/9769609/sorting-months-in-r
  sort.month <- function(x, dataframe = NULL){
    y <- data.frame(m1 = month.name, m2 = month.abb, n = 1:12)
    z <- if(max(nchar(x)) == 3) match(x, y[, 'm2']) else match(x, y[, 'm1'])
    x <- if(is.null(dataframe)) x else dataframe
    h <- data.frame(z, x)
    h[order(z), ][, -1]
  }
  my_table <- sort.month(my_table$Month, my_table)
  colnames(my_table) <- c("Month","2011","2012","2013","2014","2015","2016",
                          "2017","2018","2019","2020","2021","2022","2023","2024")
  str(my_table)
  
#Install flextable (if not installed) for creating a publication-ready tables:
#install.packages("flextable")
#sudo apt-get install libcairo2-dev  # Did this first in terminal (if installation
# didn't work)
  library(flextable)
  
#Set defaults (font and theme):
  get_flextable_defaults()
  set_flextable_defaults(
    font.family = "Times New Roman",
    theme_fun = theme_vanilla,
    table.layout = "autofit")
  
#Table of month and years
my_table <- flextable(my_table)
  
#To save as image (saving can be done directly from "export" option too):
save_as_image(my_table,"month and year wise data of checklists (before filtering).png")
#You can also do this in the Plots panel on the bottom right in RStudio
  
#Plot the distribution of checklists by month
  
#Arrange months chronologically and not alphabetically
  my$month = factor(my$month, levels = month.name)
  head(my$month)
  head(my)
  
ggplot(my) +
    aes(x = month) +
    geom_bar(stat = "count", fill = "purple") +
    geom_text(aes(label = after_stat(count)), stat = "count", vjust = 1.5,
              colour = "white")+
    labs(x = "Month", y = "Number of checklists",
         title = "Month-wise number of checklists (2011-2024)") +
    scale_y_continuous(breaks = seq(from = 0, to = 300, by = 25))+
    theme(axis.text = element_text(colour = "black", size = rel(1))) +
    theme(axis.title.y = element_text(size = rel(1.2), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
  
#Plot the distribution of checklists by year
  
table(my$year)
  
ggplot(my) +
    aes(x = year) +
    geom_bar(stat = "count", fill = "red") +
    geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.2,
              colour = "black")+
    labs(x = "Year", y = "Number of checklists",
         title = "Year-wise number of checklists") +
    scale_y_continuous(breaks = seq(from = 0, to = 450, by = 50))+
    theme(axis.text = element_text(colour = "black", size = rel(1))) +
    theme(axis.title.y = element_text(size = rel(1.2), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
  
#Distribution of checklists by species
str(ebd_in_poly)

#Create a table of species and number of times they were reported:
species <- data.frame(table(ebd_in_poly$common_name))
colnames(species) <- c("Common name", "No. of checklists in which the species appear")
is.data.frame(species)
str(species)
# $ Common name : Factor w/ 236 levels "Alpine Swift",..: 1 2 3 4 5 6 7 8 9 10 ...
# means that no. of species recorded in the area = 236
nrow(species)
  
#Cut it up into 10 parts (in this case) so that flextable docx handle the table:
split_species <- split(species, factor(sort(rank(row.names(species))%%10)))
str(split_species)
  
  flextable(split_species$'0')
  flextable(split_species$'1')
  flextable(split_species$'2')
  flextable(split_species$'3')
  flextable(split_species$'4')
  flextable(split_species$'5')
  flextable(split_species$'6')
  flextable(split_species$'7')
  flextable(split_species$'8')
  flextable(split_species$'9')
  
#Manually screenshotted after zooming in Viewer for now. Find a code fix for this.
  

  








