## Load Libraries and Data
# Load libraries
library(maptools)
library(spdep)
library(classInt)
library(RColorBrewer)
library(igraph)
library(ggplot2)
library(reshape2)

# Set the working directory and import data
fsa.shp.name = "data/Montreal_FSA_2001/Montreal_FSA_Census_2001_MTM"
fsa.visit.name = "data/montreal_ili_visits.csv"

# Read shapefile with FSA boundaries and census variables
fsa.shp = readShapePoly(fsa.shp.name, IDvar='FSA')

# Read file with data on ILI visits by FSA
fsa.visits = read.table(fsa.visit.name, header=TRUE, sep=",")

## Graph Visualization Setup
# Define nodes and edges for graph representation of Montreal regions
nodes <- c("Downtown", "Plateau", "Old Montreal", "Mile End", "Outremont")
edges <- matrix(c(
  "Downtown", "Plateau", 10,
  "Downtown", "Old Montreal", 20,
  "Plateau", "Mile End", 30,
  "Old Montreal", "Mile End", 40,
  "Mile End", "Outremont", 50
), byrow = TRUE, ncol = 3)

# Create a dataframe from edges
edges_df <- as.data.frame(edges)
names(edges_df) <- c("from", "to", "weight")
edges_df$weight <- as.numeric(edges_df$weight)

# Create a graph from the edge list
g <- graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes)

# Create adjacency matrix
adj_matrix <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)

# Convert the matrix to a format suitable for ggplot
matrix_df <- as.data.frame(adj_matrix)
matrix_df$location <- rownames(matrix_df)
matrix_long <- melt(matrix_df, id.vars = "location")

# Plotting the matrix
ggplot(matrix_long, aes(x = location, y = variable, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Connectivity Matrix of Montreal", x = "Location", y = "Connected to", fill = "Strength") +
  theme_minimal()

## 1. Plot Census Data Using Different Types of Maps
## 1.1 Plot median family income as choropleth map
nclr = 5
plotclr = brewer.pal(nclr,"Greens") 
class = classIntervals(fsa.shp@data$Median_Fam, nclr, style="quantile")
colcode = findColours(class, plotclr)

plot(fsa.shp)
plot(fsa.shp, col=colcode, add=TRUE)
title(sub="Median Family Income", main='Classed Choropleth')
legend('topleft', legend=names(attr(colcode, "table")), fill=attr(colcode,"palette"), cex=0.9, bty='n')

## 1.2 Plot median family income as a proportional symbol or 'bubble map'
max.symbol.size = 5
min.symbol.size = 1
plotvar = fsa.shp@data$Median_Fam
symbol.size = ((plotvar - min(plotvar)) / (max(plotvar) - min(plotvar)) * (max.symbol.size - min.symbol.size) + min.symbol.size)

plot(fsa.shp)
mtl.cntr = coordinates(fsa.shp)
points(mtl.cntr, pch=16, col=colcode, cex=symbol.size)
title(sub="Median Family Income", main="Bubble Plot")
legend('topleft', legend=names(attr(colcode, "table")), fill=attr(colcode,"palette"), cex=0.9, bty='n')

## 2. Plot Crude ILI Visit Rates
ili = data.frame(fsa=fsa.shp@data$FSA, pop=fsa.shp@data$POP96)
ili = merge(fsa.visits, ili, by='fsa')
ili$rate = (ili$visits / 3) / ili$pop * 1000  # annual rate per 1,000 population

groups = c(3,5,7,9)
par(mfrow=c(2,2))

for (group in groups) {
  plotclr = brewer.pal(nclr, 'Reds')
  class.crude = classIntervals(ili$rate, group, style='quantile', dataPrecision=0)
  colcode.crude = findColours(class.crude, plotclr)
  
  plot(fsa.shp)
  plot(fsa.shp, col=colcode.crude, add=T)
  title(sub="Crude Rates of ILI Visits by FSA (Annual Visits per 1,000)")
  legend('topleft', legend=names(attr(colcode.crude, "table")), fill=attr(colcode.crude,"palette"), cex=0.9, bty='n')
}
par(mfrow=c(1,1))

## 3. Connectivity Analyses
# Build a connectivity matrix and perform triangulation and nearest neighbors analysis
fsa.nb.tri = tri2nb(coordinates(fsa.shp), row.names=fsa.shp@data$FSA)
k1 = knn2nb(knearneigh(coordinates(fsa.shp)))
max.distance = max(unlist(nbdists(k1, coordinates(fsa.shp))))
fsa.nb.knn = dnearneigh(coordinates(fsa.shp), 0, max.distance, fsa.shp@data$FSA)

# Plot connectivity by triangulation
plot(fsa.shp, border='darkgrey', las=1, main='Connectivity by Triangulation')
plot(fsa.nb.tri, coordinates(fsa.shp), add=TRUE)

# Plot connectivity by nearest neighbors
plot(fsa.shp, border='darkgrey', las=1, main='Connectivity by Nearest Neighbors')
plot(fsa.nb.knn, coordinates(fsa.shp), add=TRUE)

## 4. Smooth Rates

ili.ebe = EBlocal(ili$visits/(3*52), ili$pop, fsa.nb)
ili.ebe[1:10,]

nclr = 5
plotclr = brewer.pal(nclr, 'Reds')
class.raw = classIntervals(round(ili.ebe$raw*10000,0), nclr, style='quantile')
# Note that we take the breaks from the crude rate map and pass them to the 
#  smoothed rates map so that we have the same class boundaries for both maps
class.est = classIntervals(round(ili.ebe$est*10000,0), n=nclr, style='fixed',
                           fixedBreaks=class.raw$brks)

colcode.raw = findColours(class.raw, plotclr)
colcode.est = findColours(class.est, plotclr)

par(mfrow=c(1,2))

plot(fsa.shp)
plot(fsa.shp, col=colcode.raw, add=T)
title(sub="Raw Rates of ILI Visits (Weekly Visits per 10,000)")
legend('topleft', legend=names(attr(colcode.raw, "table")), fill=attr(colcode.raw,"palette"), cex=0.9, bty='n')

plot(fsa.shp)
plot(fsa.shp, col=colcode.est, add=T)
title(sub="EBE (Local) Smoothed Rates of ILI Visits (Weekly Visits per 10,000)")
legend('topleft', legend=names(attr(colcode.est, "table")), fill=attr(colcode.est,"palette"), cex=0.9, bty='n')

par(mfrow=c(1,1))

ili.ebe$diff = ili.ebe$raw - ili.ebe$est
class.diff <- cbind(fsa.visits$fsa,ili.ebe$diff)
class.diff[,-1.683e-03]
which[class.diff==1.792e-03]

