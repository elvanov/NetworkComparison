
## Networks' comparison

This script compares several aspects of a set of networks. We start by reading the 
networks, using the `ReadNetworks` function, which gets as input the directory
in which the networks are stored (`netDir` parameter). Each network 
should be an independent file, and several formats are accepted (e.g., CSV and 
GML). Please note that we assume that the networks will have equivalent node names 
(i.e., they should all be subsets of a common set of nodes). Importantly, you 
can define whether the networks are directed or not using the `directed` parameter.

NOTE. You can also generate a predefined set of 4 toy-networks with the `GenerateToyNetworks`
function, giving as parameter the path of the folder where you want to store them.

```{r readNetworks}
# create toy networks
netDir <- "~/ToyNetworks/"
# # GenerateToyNetworks(netDir)
directed <- FALSE

networks <- ReadNetworks(netDir, directed = directed, format = "csv")

# netDir <- "~/Documents/C_elegans/LiesaData/Networks_celeganspathogenmicrobiome-subset_pos-neg_2022-01-12_1517/"
# directed <- FALSE
# 
# # read networks
# networks <- ReadNetworks(netDir, directed = directed, pattern = "neg", format = "gml")

print(paste0(length(networks), " networks read"))
print("Number of nodes per network: ")
lapply(networks, vcount)
print("")
print("Number of edges per network: ")
lapply(networks, ecount)

```

Now that we have the networks to analyze, we will proceed to calculate some basic
statistics measurements and to make a plot for each of them. In order to do so,
use the `CalculateNetworkStats` function, giving a parameter the networks you just
read with `ReadNetworks`. To generate the plots, we will use the `printStatsPlots`
function, which needs the results of the `CalculateNetworkStats` function as input
parameter.

The statistics that will be calculated and displayed are the following:

* **Density:** Radio of number of edges and the number of possible edges (i.e.,
the number of edges of the fully connected graph).
* **Diameter:** Size of the shortest path between the furthest vertices, i.e.,
the length of the longest geodesic distance.
* **Average path length:** Mean of the shortest paths between each pair of
vertices.
* **Clustering coefficient:** Number of triangles (i.e., node A connected to
B and C, and node C connected to B) in the graph. It is a measure of the degree
to which nodes in a graph tend to cluster together. In other words, the 
clustering coefficient tells you how well connected the neighborhood of 
every node is. It goes from 0 (i.e., no connections between the nodes in the 
neighborhood of each node) to 1 (i.e., the neighborhood of each node is fully 
connected).
* **Degree:** Number of edges connected to a node (i.e., its adjacent edges).
* **Closeness:** Reciprocal of the average shortest paths between a node and all
the other nodes in a graph. It is a measure of the centrality of a node. When
normalized, it allows comparing the closeness between graphs of different sizes.
NOTE. The higher the closeness of a node, the more central that node is.
* **Betweenness:** Number of shortest paths going through a node or an edge. We
calculate the betweenness of the nodes.

```{r calculateStats }
# calculate stats and generate plots
statistics <- CalculateNetworkStats(networks)
printStatsPlots(statistics)
```

We will now calculate the overlap of the nodes and edges of all of our networks,
using the `CalculateOverlap` function, which needs the networks to analyze as
input parameter. The `CalculateOverlap` function the will also generate two
upset plots, one to show the overlap of nodes, and the other for the edges. It
is to note that the overlap of edges does not take direction into account and
it will thus consider the networks as undirected.

Please note that the upset plots show an intersection of sets, but also a
difference, e.g., the upset plot of the overlap of nodes between 2 networks, will
show:
* The common nodes (i.e., intersection) of nodes in both networks
* The nodes present in the first network, but not in the second one
(difference), and
* The nodes present in the second network, but not in the first one.

```{r calculateOverlap}
# calculate overlap of nodes and edges
CalculateOverlap(networks)
```

If you want to get the list of overlapping nodes, you can use the 
`getOverlappingNodes` function, giving as input the list of networks and the 
list of index of the networks you want to take into account. 
NOTE. If no list of index is provided, all the networks are considered for the 
calculus of the overlap.

```{r getOverlappingNodes}
# get list of overlapping nodes 
getOverlappingNodes(networks)
```
