# Quelites
Rscript and data for the paper entitled: Edible wild plants and traditional indigenous knowledge in Zongolica, Mexico
#########################################################
#########################################################
####R SCRIPT FOR THE DATA ANALYSIS#######################
####Edible wild plants and traditional indigenous knowledge in Zongolica, Mexico

# We set the directory 
getwd()
setwd('C:/???')

# Install packages
library(dplyr)
require(openintro)
require(knitr)
library(openintro)
library(knitr)
library(ggplot2)
library(car)
#############################################
# We add the data "Quelites"
Quelites <- read.csv("Quelites.csv", header=TRUE)
Quelites
head(Quelites)

############ We want to know if there are differencesin the number
#### of species people knows according to each ecosystem.

summary(Quelites$Index1)
summary(Quelites$Ecosystem)

tapply(Quelites$Index1, Quelites$Ecosystem, summary)
#the plot using both variables
plot(Quelites$Index1 ~ Quelites$Ecosystem)
# compare means and p-values
p.Base <- aov(Quelites$Index1 ~ Quelites$Ecosystem)
summary(p.Base)
## running some packages we need

library(Rmisc)
library(lattice)
# We must select the variables to analyze age range by ecosystem
Barra <- Quelites %>% select(Age, Index1, Ecosystem)
Barra
# The SE for both variables based on the Index1
Barras <- summarySE(Barra, measurevar="Index1", groupvars=c("Age","Ecosystem"))
Barras

# Then we plot Age range (x), the number of species each person knows
# for the two different ecosystems

figura1 <- ggplot(Barras, aes(x=Age, y=Index1,fill=Ecosystem, levels = rev(Barra2$Ecosystem))) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Index1-se, ymax=Index1 +se),
                width=.2,                    # Ancho de las barras de error
                position=position_dodge(.9)) + theme_bw()+ 
  scale_fill_manual("Ecosystem", values = c("Pine-oak forest" = "gray47", "Montane cloud forest"= "grey"))+
  labs(x = "Age range", y = " Number of quelites people knows")
  
figura1 + theme(axis.line = element_line(colour = "black"))
figura1
###############################################################
# Now we analyze the number of mentions each specie share represented
# in a heatmap and grouping the species using a cluster analysis
# First we constructed a file "que" based on the affiliation network
#Affiliation networks represents two modes, 
# the first one to the set of actors A 
# and the second to the events E, with the number of actors n 
# and the number of events m [3],. 
# So, an actor can take part in more than one event. 
# The initial set of data in this study is based on a two-mode network 
# or affiliation network where actors are the "persons"
# and events are the "species or quelites" they know. 
# We then represent an affiliation network by an affiliation 
# matrix A of size n × m, where links exists only between actors
# and species. However, actors and species can be also 
# represented independently as an "one-mode" network.

# Using the one-mode network approach, we have two set of data, 
# a network of actors and other related with species. 
# A network of species is formed by a set of species N and a set of links L. 
# We have a relationship between two species if they share at least one actor. 
# We then construct an adjacency matrix of this network XN . 
# An element of XN is noted by X_ik^N with i, k ∈ N . 
# In this case for simple data manipulation we constructed
# the data set "que" which represents the affiliation network
# persons and species and the links

# we use the igraph package to read the affiliation network"que"
## library (igraph)

quelites <- read.graph(file.choose("que"), format = "pajek")
quelites

# We convert the "two - mode" to "one-mode" network
bipartite_matrix <- as_incidence_matrix(quelites)

bipartite_matrix

# Now we can multiply bipartite_matrix by 
# its transpose: t(biparite_matrix).
# traspose the bipartite matrix

t(bipartite_matrix)

# Similar to the %in% operator we saw earlier, R gives us a 
# special operator to use for matrix multiplication: %*%.

event_matrix_prod <- t(bipartite_matrix) %*% bipartite_matrix 
diag(event_matrix_prod) <- 0
# The event matrix prod represent the number of persons each
# species "share"

event_matrix_prod

#########################################################
### Correlation heatmap for quelites
library(reshape2)
melted_cormat <- melt(event_matrix_prod)
head(melted_cormat)
library(ggplot2)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(event_matrix_prod){
  event_matrix_prod[upper.tri(event_matrix_prod)] <- NA
  return(event_matrix_prod)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(event_matrix_prod){
  event_matrix_prod[lower.tri(event_matrix_prod)]<- NA
  return(event_matrix_prod)
}

upper_tri <- get_upper_tri(event_matrix_prod)
upper_tri
# Melt the correlation matrix
library(reshape2)
melted_cor <- melt(upper_tri, na.rm = TRUE)

######################################################
# Reorder the correlation matrix
event_matrix_prod <- reorder_cor(event_matrix_prod)
upper_tri <- get_upper_tri(event_matrix_prod)
# Melt the correlation matrix
melted_cor <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cor, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low="white", high="red", mid="blue", 
                       midpoint = 170, limit = c(1,350), space = "Lab", 
                       name="Number of mentions")  +
  theme_minimal()+ # minimal th eme
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed()
# Print the heatmap 
print(ggheatmap)
## adding  numbers and legends
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

### Finally we carried out a cluster analysis
# cluster rows for quelites
hc.rows <- hclust(dist(event_matrix_prod))
plot(hc.rows)

# transpose the matrix and cluster columns
hc.cols <- hclust(dist(t(event_matrix_prod)))

# draw heatmap for first cluster
heatmap(event_matrix_prod[cutree(hc.rows,k=2)==1,], Colv=as.dendrogram(hc.cols), scale='none')

# draw heatmap for second cluster
heatmap(event_matrix_prod[cutree(hc.rows,k=2)==2,], Colv=as.dendrogram(hc.cols), scale='none')

