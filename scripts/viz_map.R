if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("dkahle/ggmap", ref = "tidyup", force=TRUE)
install.packages("ggrepel")
library(ggrepel)
library(tidyverse)
library(ggplot2)
library(ggmap)

# Register my ggmap
ggmap::register_google(key = )

dat<- readRDS("../local_data/clean_pheno.RDS")

#####
# CREATE MAP OF SAMPLES
#####
# Remove replicated genotypes
gpsDat<- subset(dat, !duplicated(dat$Genotype_ID))

# Remove entries with no river label
gpsDat<- gpsDat[-which(gpsDat$River==""),]

##TODO: Should I remove genotypes with generic gps data?

# Create df for river system labels
river<- aggregate(gpsDat[,4:5], list(gpsDat$River), median)
colnames(river)[1]<-"River"

# Create map
p <- ggmap(get_googlemap(center = c(lon = -122.1901, lat = 47.10007),
                         zoom = 7, scale = 2,
                         maptype ='terrain',
                         color = 'color'))
# Create visualization
p + geom_point(aes(Longitude, Latitude), data = gpsDat, colour= "#FF00FF", alpha=0.25, size= 2) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=18), legend.text = element_blank(),
        legend.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(title= "Source Location of Poplar Genotypes") +
  geom_point(aes(x = Longitude, y = Latitude, stroke = 2), colour="FF9900", data = river, size =5) + 
  geom_label_repel(aes(Longitude, Latitude, label = River), data=river, family = 'Times', size = 3.5, 
    box.padding = 0.2, point.padding = 0.3, segment.color = 'grey50') 




                 