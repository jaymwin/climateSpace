
library(raster)
library(tidyverse)
library(rgdal)
library(rgeos)
library(sp)
library(ggrepel)
library(here)


# download WORLDCLIM 30-year average (~1960-1990) -------------------------

# mean temperature
annual_tmean = getData('worldclim', var='tmean', res=2.5)
avg_tmean <- mean(annual_tmean, na.rm=TRUE)

# mean precipitation
annual_prcp = getData('worldclim', var='prec', res=2.5)
avg_prcp <- mean(annual_prcp, na.rm=TRUE)


# download North America shapefile ----------------------------------------

# states and provinces
stateprov <- readOGR(here('data/10m_states_provinces'), 'ne_10m_admin_1_states_provinces')

# subset to USA, Canada
NA_only <- subset(stateprov, name %in% c('Alberta', 'Saskatchewan', 'Alaska', 'Manitoba',
                                         'Rhode Island', 'Michigan', 'Minnesota', 'Mississippi',
                                         'Connecticut', 'Arizona', 'Washington', 'Maryland', 'Oklahoma',
                                         'Idaho', 'Oregon', 'South Carolina', 'North Carolina',
                                         'North Dakota', 'South Dakota', 'British Columbia', 'Massachusetts',
                                         'Montana', 'California', 'Nevada', 'New York', 'New Jersey', 'Arkansas',
                                         'Alabama', 'Florida', 'Tennessee', 'Wisconsin', 'Texas', 'Utah',
                                         'New Mexico', 'Louisiana', 'Georgia', 'New Hampshire', 'Vermont', 
                                         'Maine', 'Ohio', 'Virginia', 'West Virginia', 'Indiana', 'Illinois',
                                         'Wyoming', 'Colorado', 'Missouri', 'Kansas', 'Nebraska', 'Iowa',
                                         'Kentucky', 'Delaware', 'Pennsylvania', 'Ontario', 'District of Columbia',
                                         'New Brunswick', 'Prince Edward Island', 'Nunavut', 'Northwest Territories',
                                         'Yukon', 'Newfoundland and Labrador', 'Québec', 'Nova Scotia'))
plot(NA_only)

# crop it
rect <- as(extent(-170, -52, 10, 74), 'SpatialPolygons')
crs(rect) <- crs(NA_only)
NA_only <- crop(NA_only, rect)
plot(NA_only)

# plot avg_tmean, avg_prcp
plot(avg_tmean)
plot(avg_prcp)

# dissolve polygon
NA_only_dissolved <- gUnaryUnion(NA_only)
plot(NA_only_dissolved)

# crop avg_tmean, avg_prcp to NA
crs(avg_tmean) <- crs(NA_only_dissolved)
avg_tmean_NA <- mask(avg_tmean, NA_only_dissolved)
plot(avg_tmean_NA)

crs(avg_prcp) <- crs(NA_only_dissolved)
avg_prcp_NA <- mask(avg_prcp, NA_only_dissolved)
plot(avg_prcp_NA)

# tmean data frame
tmean_spdf <- as(avg_tmean_NA, "SpatialPixelsDataFrame")
tmean_df <- as.data.frame(tmean_spdf)
colnames(tmean_df) <- c("tmean", "x", "y")

# tmean data frame
prcp_spdf <- as(avg_prcp_NA, "SpatialPixelsDataFrame")
prcp_df <- as.data.frame(prcp_spdf)
colnames(prcp_df) <- c("prcp", "x", "y")

# join together
climate_space <- tmean_df %>%
  left_join(., prcp_df) 

# plot it
climate_space %>%
  ggplot(., aes(prcp*10, tmean/10)) +
  geom_point(size = .05, alpha = 1, color = 'blue', stroke = 0) +
  theme_classic() +
  xlab('30-yr annual precipitation (mm)') +
  ylab('30-yr annual temperature (C)') +
  expand_limits(y = c(-30, 30)) +
  expand_limits(x = c(0, 3500)) +
  scale_y_continuous(breaks = c(-30, -20, -10, 0, 10, 20, 30)) +
  scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500))

ggsave(here('climate_space.png'), width = 4, height = 3)


dod <- read_csv(here('data/DoD_sites_latlong.csv'))

# conver to spatial points  
pointsSp <- dod %>%
  dplyr::select(lon, lat) %>%
  SpatialPoints(CRS(wgs))

# Extract covariate data to points:
dod <- bind_cols(
  dod,
  raster::extract(avg_tmean_NA, pointsSp) %>%
    as.data.frame %>%
    setNames(c('tmean'))) 
dod

dod <- bind_cols(
  dod,
  raster::extract(avg_prcp_NA, pointsSp) %>%
    as.data.frame %>%
    setNames(c('prcp'))) 
dod

# plot it
climate_space %>%
  ggplot(., aes(prcp*10, tmean/10)) +
  geom_point(size = .05, alpha = 1, color = 'blue', stroke = 0) +
  geom_point(data = dod, aes(prcp*10, tmean/10), color = 'red', size = 0.5) +
  theme_classic() +
  xlab('30-yr annual precipitation (mm)') +
  ylab('30-yr annual temperature (C)') +
  expand_limits(y = c(-30, 30)) +
  expand_limits(x = c(0, 3500)) +
  scale_y_continuous(breaks = c(-30, -20, -10, 0, 10, 20, 30)) +
  scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500))

ggsave(here('climate_space.png'), width = 4, height = 3)



range <- readOGR(here('data/kestrel_range_map'), 'kestrel_simplified_dissolved')
plot(range)
plot(NA_only_dissolved)


clip <- gIntersection(range, NA_only_dissolved) 
plot(clip)

# crop avg_tmean, avg_prcp to NA
crs(avg_tmean) <- crs(clip)
avg_tmean_range <- mask(avg_tmean, clip)
plot(avg_tmean_range)

crs(avg_prcp) <- crs(clip)
avg_prcp_range <- mask(avg_prcp, clip)
plot(avg_prcp_range)


# tmean data frame
tmean_range_spdf <- as(avg_tmean_range, "SpatialPixelsDataFrame")
tmean_range_df <- as.data.frame(tmean_range_spdf)
colnames(tmean_range_df) <- c("tmean", "x", "y")

# prcp data frame
prcp_range_spdf <- as(avg_prcp_range, "SpatialPixelsDataFrame")
prcp_range_df <- as.data.frame(prcp_range_spdf)
colnames(prcp_range_df) <- c("prcp", "x", "y")

# join together
climate_space_range <- tmean_range_df %>%
  left_join(., prcp_range_df) 

# plot it
g1 <- ggplot() +
  geom_point(data = climate_space, aes(prcp*10, tmean/10), size = .05, alpha = 1, color = 'blue', stroke = 0) +
  geom_point(data = climate_space_range, aes(prcp*10, tmean/10), size = .05, alpha = 1, color = 'aquamarine', stroke = 0) +
  geom_text_repel(data = dod, aes(prcp*10, tmean/10, label = site), size = 1.5, segment.size = 0.2) +
  geom_point(data = dod, aes(prcp*10, tmean/10), color = 'red', size = 0.5) +
  theme_classic() +
  xlab('30-yr annual precipitation (mm)') +
  ylab('30-yr annual temperature (C)') +
  expand_limits(y = c(-30, 30)) +
  expand_limits(x = c(0, 3500)) +
  scale_y_continuous(breaks = c(-30, -20, -10, 0, 10, 20, 30)) +
  scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500))

ggsave(here('climate_space2.png'), width = 4, height = 3)

  
plot(clip)
plot(NA_only_dissolved)  

clip_df <- fortify(clip)
NA_only_dissolved_df <- fortify(NA_only_dissolved)


map <- ggplotGrob(
ggplot() +
  geom_polygon(data = NA_only_dissolved_df, aes(long, lat, group = group), fill = 'blue', color = NA) +
  geom_polygon(data = clip_df, aes(long, lat, group = group), fill = 'aquamarine', color = NA) +
  geom_point(data = dod, aes(lon, lat), color = 'red', size = 0.5) +
  coord_map(ylim = c(25, 72)) +
  xlab("Longitude") + ylab("Latitude") +
  theme_classic(base_size = 5) 
)
  
g1 +
  annotation_custom(grob = map, xmin = 1700, xmax = 4000,
                    ymin = -32, ymax = -5)  

ggsave(here('climate_space3.png'), width = 4, height = 3)


  



# Get background country data:

conus <- getData('GADM', country='USA', level=1) %>%
  subset(!(NAME_1 %in% c('Alaska', 'Hawaii'))) %>%
  gUnaryUnion

mexico <- getData('GADM', country='Mexico', level=0) %>%
  gUnaryUnion

cuba <- getData('GADM', country='Cuba', level=0) %>%
  gUnaryUnion

northAmerica <- do.call(bind, list(conus, mexico, cuba))
