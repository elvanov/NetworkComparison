library(igraph)
library(tidyverse)
library(UpSetR)
library(grid)
library(tidyverse)
library(plyr)
library(ggplot2)
library(ggrepel)
# ################################################################################
# # rename files
# RenameFiles <- function(netDir, patternToRemove = "") {
#   netDir <- 
#     ifelse(
#       substr(netDir, nchar(netDir), nchar(netDir)) == "/", 
#       substr(netDir, 1, (nchar(netDir)-1)), 
#       netDir)
#   
#   # list all files of the given format in the corresponding directory
#   files <- 
#     list.files(
#       path = netDir, 
#       full.names = TRUE,
#       pattern =  paste0(".*", patternToRemove, "(.*).gml")
#     ) 
#   
#   # rename files
#   lapply(
#     files, 
#     function(X) {
#       newName <- str_replace(X, paste0(".*", patternToRemove, "(.*.gml)"), "\\1")
#       file.rename(X, paste0(netDir, "/", newName))
#     }
#   )
#   return()
# }


# Function to generate 4 toy networks 
# INPUT:
#   netDir - directory to store the toy networks
GenerateToyNetworks <- function(netDir) {
  n1 <- matrix(data = c("A","B", "B","C", "B","E", "C","D"), ncol = 2, byrow = TRUE)
  n2 <- matrix(data = c("A","B", "B","F", "A","F", "C","F", "C","E", "C","D"), ncol = 2, byrow = TRUE)
  n3 <- matrix(data = c("A","F", "E","F", "C","F", "C","E"), ncol = 2, byrow = TRUE)
  n4 <- matrix(data = c("A","B", "A","F", "B","F", "B","C", "A","D", "C","D", "F","D"), ncol = 2, byrow = TRUE)
  colnames(n1) <- colnames(n2) <- colnames(n3) <- colnames(n4) <- c("node1", "node2")
  write.csv(n1, file = paste0(netDir, "Network1.csv"), row.names = FALSE, quote = FALSE)
  write.csv(n2, file = paste0(netDir, "Network2.csv"), row.names = FALSE, quote = FALSE)
  write.csv(n3, file = paste0(netDir, "Network3.csv"), row.names = FALSE, quote = FALSE)
  write.csv(n4, file = paste0(netDir, "Network4.csv"), row.names = FALSE, quote = FALSE)
}


# Function to read all the networks stored in a given directory
# INPUTS:
#   netDir - directory containing all the networks. 
#            NOTE. The networks must be stored in csv format
#   directed - boolean value indicating if the networks are directed or not, 
#              the default value is FALSE (i.e., undirected networks)
#   pattern - pattern that the file names must contain to be read. If all the
#             networks in netDir are to be read, pattern can be omitted
#   format - files' format. NOTE. If the format is "csv", then the files must 
#            contain 2 columns (source - target) and each row must be an edge
# OUTPUT:
#   list of igraph objects, one per network to compare
# ReadNetworks <- function(netDir, directed = FALSE) {
  # # list all csv files in the corresponding directory
  # files <- list.files(path = netDir, pattern = "*.csv", full.names = TRUE)
  # networks <- list()
  # 
  # # loop through the files
  # for (i in files) {
  #   # read file content 
  #   n <- as.matrix(read.csv(i))
  #   
  #   # create igraph object
  #   net <- graph_from_edgelist(el = n, directed = directed)
  #   
  #   networks[[str_replace(i, ".+/(.*).csv", "\\1")]] <- net
  # }
ReadNetworks <- function(netDir, directed = FALSE, pattern = "", format = "csv") {
  
  # if exists, remove final "/" from the network directory name
  netDir <- 
    ifelse(
      substr(netDir, nchar(netDir), nchar(netDir)) == "/", 
      substr(netDir, 1, (nchar(netDir)-1)), 
      netDir)
  
  # list all files of the given format in the corresponding directory
  files <- 
    list.files(
      path = netDir, 
      pattern = paste0(".*", pattern, ".*.", format, "$"), 
      full.names = TRUE
      )
  
  networks <- list()
  
  # loop through the files
  for (i in files) {
    # verify file format
    if (format == "csv") {
      # read file content 
      n <- as.matrix(read.csv(i))
      
      # create igraph object
      net <- graph_from_edgelist(el = n, directed = directed)
      
    } else {
      net <- read_graph(i, format = format)
      
      if (directed == FALSE) {
        net <- as.undirected(net)
      }
    }
    
    # add network to named list 
    networks[[str_replace(i, paste0(".+/(.*).", format), "\\1")]] <- net
  }  
  
  return(networks)
}


# Function to calculate and create a table with and plot the following statistics:
#    - Density (no self loops are considered)
#    - Diameter (if there is more than one connected component, the largest
#        diameter will be returned, independently of the size of the connected
#        component)
#    - Average degree
#    - Average path length
#    - Clustering coefficient (the networks are considered as undirected)
# It also plots the degree distribution, and upset plots of the overlap of nodes
# and edges
# INPUT:
#   networks - list of igraph objects, as returned by ReadNetworks()
# OUTPUT: 
#   data frame with the statistics
CalculateNetworkStats <- function (networks) {
  # calculate number of networks
  noNet <- length(networks)
  
  ########################################################################################################################
  ### --- ADD CENTRALITY MEASURE - CHECK https://www.datacamp.com/community/tutorials/centrality-network-analysis-R 
  ### , WHICH SHOULD FOLLOW A POWER LAW DSITRIBUTION, I.E., SCALE-FREE NETWORK
  
  
  # initialize empty data frames
  netStats <- 
    data.frame( name = character(noNet), density = double(noNet),
      diameter = double(noNet), avDegree = double(noNet), noCC = double(noNet),
      avPathLenght = double(noNet), clustCoeff = double(noNet) )
  
  networksDegrees <- 
    data.frame(network = character(), node = character(), degree = double())
  
  networkCloseness <- 
    data.frame(network = character(), node = character(), closeness = double())

  networkBetweenness <- 
    data.frame(network = character(), node = character(), betweenness = double())
  
  # loop through all the networks
  for (i in seq_len(noNet)) {
    # calculate statistics 
    netStats$density[[i]] <- edge_density(networks[[i]])
    netStats$diameter[i] <- diameter(networks[[i]], unconnected = TRUE)
    netStats$avDegree[i] <- mean(degree(networks[[i]]))
    netStats$avPathLenght[i] <- 
      mean_distance(networks[[i]], directed = is_directed(networks[[i]]))
    netStats$clustCoeff[i] <- transitivity(networks[[i]], type = "global")
    
    # calculate the degree distribution
    degreeD <- degree.distribution(networks[[i]])

    # calculate and save degrees
    networksDegrees <- 
      rbind(
        networksDegrees,
        data.frame(
          network = names(networks)[i], 
          node = names(V(networks[[i]])), 
          degree = degree(networks[[i]])
          )
      )
    
    # get connected components of the network
    connectedComponents <- components(networks[[i]])
    
    # save total number of connected components
    netStats$noCC[[i]] <- connectedComponents$no

    # calculate the closeness in each component
    closeness <-
      unlist(
        lapply(
          seq_len(connectedComponents$no), 
          function(X) {
            # get nodes' IDs of the nodes in current connected component
            nodesIDs <- 
              V(networks[[i]])[which(names(V(networks[[i]])) %in% 
              names(connectedComponents$membership[connectedComponents$membership == X]))]
            
            # calculate closeness of current Connected Component (CC)
            # NOTE. If the closeness is normalized, the most central node (i.e., 
            # the one with the shortest paths to the rest of the nodes) is the node with the 
            # highest closeness, or the lowest value otherwise
            # NOTE 2. Each CC is considered as an independent graph
            closeness(induced_subgraph(networks[[i]], vids = nodesIDs), normalized = TRUE)
          }
        )   
      )
    
    networkCloseness <- 
      rbind(
        networkCloseness,
        data.frame(
          network = names(networks)[i],
          node = names(closeness),
          closeness = closeness
          )
        ) 
    
    # calculate and save the betweenness of the nodes
    networkBetweenness <- 
      rbind(
        networkBetweenness,
        data.frame(
          network = names(networks)[i], 
          node = names(V(networks[[i]])), 
          betweenness = betweenness(networks[[i]])
        )
      )
    
    # add network name
    netStats$name[[i]] <- names(networks)[i]
  }

  allStats <- 
    list(netStats = netStats, networksDegrees = networksDegrees, 
      networkCloseness = networkCloseness,
      networkBetweenness = networkBetweenness)
  
  return(allStats)
}


# Function to make and print several plots
# INPUT:
#   stats - list of results, as returned by the CalculateNetworkStats function
# OUTPUT: 
#   none, but prints several plots
printStatsPlots <- function(stats) {
  netStats <- stats[['netStats']]
  networksDegrees <- stats[['networksDegrees']]
  networkCloseness <- stats[['networkCloseness']]
  networkBetweenness <- stats[['networkBetweenness']]
  
  makeScatterPlot(netStats, statToPlot = "noCC", title = "Number of connected components")
  makeScatterPlot(netStats, statToPlot = "density")
  makeScatterPlot(netStats, statToPlot = "diameter")
  makeScatterPlot(netStats, statToPlot = "avPathLenght", title = "Average path length")
  makeScatterPlot(netStats, statToPlot = "clustCoeff", title = "Clustering coefficient")
  makeScatterPlot(netStats, statToPlot = "avDegree", title = "Average degree")
  makeBoxPlot(networksDegrees, statToPlot = "degree")
  makeHistogram(networksDegrees, statToPlot = "degree", binWidth = 1, minDegree = 0, maxDegree = max(networksDegrees$degree))  
  makeBoxPlot(networkCloseness, statToPlot = "closeness", bestValue = "max", showLabels = TRUE)
  makeBoxPlot(networkBetweenness, statToPlot = "betweenness", bestValue = "max", showLabels = TRUE)
  
  return()
}


# Function to make and print a scatter plot
# INPUTS:
#   netStats - data frame, as returned as part of the result of CalculateNetworkStats 
#   statToPlot - string defining the stat that wishes to be plotted. It can be 
#               "density", "diameter", "avPathLength" or "clustCoeff"
# OUTPUT: 
#   none, but prints the corresponding scatter plot
makeScatterPlot <- function(netStats, statToPlot, title = "") {
  print("")
  print(
    ggplot(
      netStats, 
      aes(
        x = name, 
        y = eval(as.name(statToPlot)), 
        color = rainbow(nrow(netStats))
        )
      ) 
    + 
    geom_point(size = 4) 
    + 
    labs(
      title = ifelse(title != "", paste0(title, " comparison\n"), paste0(str_to_title(statToPlot), " comparison\n")), 
      x = "Network", 
      y = ifelse(title != "", title, str_to_title(statToPlot))
      ) 
    +
    theme(
      axis.text.x = element_text(size = 10), 
      axis.title.x = element_text(size = 12),
      axis.text.y = element_text(size = 10), 
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
      legend.position = "none"
    )    
  )

  ### uncomment the following lines to use a legend
  # ggplot(netStats, aes(x = name, y = eval(as.name(statToPlot)), color = rainbow(nrow(netStats)))) + 
  # geom_point(size = 4) +
  # labs(title = paste0(str_to_title(statToPlot), " comparison\n"), x = "Network", y = str_to_title(statToPlot), color = "Network: ") +
  # theme(legend.position = "right") + scale_color_manual(labels = netStats$name, values = rainbow(nrow(netStats))) +
  # theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
  #         axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
  #         plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
  #         legend.text = element_text(size = 12), legend.title = element_text(size = 14))
  return()
}



# Function to make and print a boxplot
# INPUTS:
#   data - data frame, as returned as part of the result of CalculateNetworkStats 
#   statToPlot - string defining the stat that wishes to be plotted. It can be 
#               "degree" or "closeness"
#   bestValue - string to define whether the highest or lowest values are the 
#                best ones. This is only usefull if showLabels == TRUE. Possible
#                values are "max" and "min"
#   showLabels - flag to choose whether to show the label of the highest/lowest value
# OUTPUT: 
#   none, but prints the corresponding scatter plot
makeBoxPlot <- function(data, statToPlot, bestValue = "max", showLabels = FALSE) {
  print("")
  print("")
  p <- 
    ggplot(
      data, 
      aes(
        x = network, 
        y = eval(as.name(statToPlot)), 
        fill = network
      )
    ) + 
    geom_boxplot() + 
    labs(
      title = paste0(str_to_title(statToPlot), " comparison\n"), 
      x = "Network", 
      y = str_to_title(statToPlot)
    ) +
    theme(
      axis.text.x = element_text(size = 10), 
      axis.title.x = element_text(size = 12),
      axis.text.y = element_text(size = 10), 
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
      legend.position = "none"
    )
  
  # verify if labels should be added to the plot
  if(showLabels == TRUE) {
    if(bestValue == "max") {
      plotLabels <- 
        data[
          as.logical(
            ave(
              data[, eval(statToPlot)], 
              data$network, 
              FUN = function(x) x == max(x)
            )
          ),
        ]
    } else {
      plotLabels <- 
        data[
          as.logical(
            ave(
              data[, eval(statToPlot)], 
              data$network, 
              FUN = function(x) x == min(x)
            )
          ),
        ]
    }
    
    # check if there are more than 2 nodes' names to be printed for any network
    if (any(table(plotLabels$network) > 2)) {
      newPlotLabels <- data.frame()
        
      # loop through the network names
      for (net in unique(plotLabels$network)) {
        # get rows from current network
        rows <- which(plotLabels$network == net)
        
        # print some of the labels for the current network
        print(
          paste0(
            "Network ", net, " has ", length(rows), " nodes with ", bestValue, 
            " ", statToPlot, ". Here are some of them: "
            )
          )
        print(plotLabels$node[head(rows)])
        print("")
        
        # if there are more than 2 rows
        if (length(rows) > 2) {
          # pick two rows at random
          selRows <- sample(rows, 2)
          
          # add labels to data frame
          newPlotLabels <- rbind(newPlotLabels, plotLabels[selRows,])    
        } else {
          newPlotLabels <- rbind(newPlotLabels, plotLabels[rows,])    
        }
      }
      plotLabels <- newPlotLabels
    }
    
    
    # put all the nodes' names to be printed per network together
    plotLabels <-
      ddply(
        plotLabels,
        .(network),
        function(x) {
          return(c(node = paste(x$node, collapse = ", "), x = unique(x[[statToPlot]])))
        }
      )
    colnames(plotLabels)[ncol(plotLabels)] <- statToPlot
    plotLabels[[statToPlot]] <- as.numeric(plotLabels[[statToPlot]] )

    p <-
      p +
      geom_text(
        data = plotLabels,
        aes(label = node),
        position = "identity",
        vjust = ifelse(bestValue == "max", -0.3, 0.3),
        size = 2
      )
  }
  
  print(p)
  
  return()
}


# Function to calculate and print a set of histograms
# INPUTS:
#   data - data frame, as returned as part of the result of CalculateNetworkStats 
#   statToPlot - string defining the stat that wishes to be plotted. It can be 
#               "degree"
#   binWidth - width of the bins for each histogram. Default value = 1
#   minDegree - minimun degree to consider for the plot. Default value = 0
#   maxDegree - maximun degree to consider for the plot. Default value = 10
# OUTPUT: 
#   none, but prints the corresponding histograms
makeHistogram <- function(data, statToPlot, binWidth = 1, minDegree = 0, maxDegree = 10) {
  # filter data
  data <- 
    data %>% 
    filter(
      eval(as.name(statToPlot)) >= minDegree & 
      eval(as.name(statToPlot)) <= maxDegree
      )
  
  print("")
  print(
    ggplot(
      data, 
      aes(
        x = eval(as.name(statToPlot)), 
        color = network, 
        fill = network
      )
    ) + 
    geom_histogram(
      alpha = 0.5, 
      binwidth = binWidth, 
      position = "identity"
    )
  )
    
  return()
}


# Function to calculate and print plots of the overlap of nodes and edges between
# a set of networks
# INPUT:
#   networks - list of igraph objects to analyze
# OUTPUT: 
#   none, but prints the corresponding plots
CalculateOverlap <- function(networks) {
  # get the number of networks
  noNet <- length(networks)
  
  # initialize empty list to save the nodes/edges of each network
  allNodes <- list()
  allEdges <- list()
  
  # loop through the networks to get the nodes' and edges' list
  for(i in seq_len(length(networks))) {
    # save nodes' names
    allNodes[[names(networks)[i]]] <- names(V(networks[[i]]))
    
    # get edges' list
    edges <- as_edgelist(networks[[i]])
    
    # sort edges alphabetically
    edges <- t(apply(edges, 1, sort))
    
    # paste the edges
    edges <- apply(edges, 1, paste, collapse = ".")
    
    # save edges
    allEdges[[names(networks)[i]]] <- edges
  }
  
  print("")

  # make upset plot for the overlap of nodes
  makeUpsetPlot(
    allNodes, title = "Overlap of nodes",
    xLabel = "Nodes per network", yLabel = "Nodes' intersection"
    )

  print("")

  # make upset plot for the overlap of edges
  makeUpsetPlot(
    allEdges, title = "Overlap of edges",
    xLabel = "Edges per network", yLabel = "Edges' intersection"
    )

  return()  
}



# Function to get the overlapping nodes from a given set of input networks
# INPUTS:
#   networks - list of networks (igraph objects)
#   networksIndex - list of networks' index to consider for the overlap, e.g., 
#                   c(1, 3) to take the first and third networks
# OUTPUT: 
#   list of overlapping nodes
getOverlappingNodes <- function(networks, networksIndex) {
  if (exists("networksIndex")) {
    # filter network list to keep only the networks to process
    networks <- networks[networksIndex]
  }
    
  # get nodes' names
  nodes <- lapply(networks, function(X) {names(V(X))})
  
  # obtain the intersection
  overlappingNodes <- Reduce(intersect, nodes)
  
  return(overlappingNodes)
}


# Function to make and print an upset plot 
# INPUTS:
#   dataToPlot - named list containing the data to plot
#   xLabel - string defining the label for the X axis, per default it is "Set Size"
#   yLabel - string defining the label for the Y axis, per default it is "Intersection Size"
# OUTPUT: 
#   none, but prints the corresponding upset plot
makeUpsetPlot <- function(dataToPlot, title = "Overlap", xLabel = "Set Size", yLabel = "Intersection Size") {
  print("")
  
  print(
    upset(
      fromList(dataToPlot), order.by = "freq", point.size = 2, line.size = 1,
      mainbar.y.label = yLabel, sets.x.label = xLabel, empty.intersections = "on",
      # y-axis label, y-axis ruler numbers, x-axis label, x-axis ruler numbers,
      # x-axis legend content, overlap size legend
      text.scale = c(1.3, 0.9, 1, 1, 1.1, 0.9)
    )
  )

  grid.text(title, x = 0.65, y = 0.95, gp = gpar(fontsize = 16))
  
  return()
}
