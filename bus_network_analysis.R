# Bus Network Connectivity Analysis
# Compares 2019 and 2025 Konstanz bus network adjacency matrices

library(tidyverse)
library(igraph)

# =============================================================================
# FUNCTIONS
# =============================================================================

create_adjacency_matrix <- function(data) {
    # Input validation
    if (nrow(data) == 0) {
        stop("ERROR: Data is empty!")
    }

    # Get all unique stops
    all_stops <- unique(c(data$Von_Haltestelle, data$Nach_Haltestelle))
    n_stops <- length(all_stops)

    if (n_stops == 0) {
        stop("ERROR: No stops found in data!")
    }

    # Create empty adjacency matrix
    adj_matrix <- matrix(0,
        nrow = n_stops, ncol = n_stops,
        dimnames = list(all_stops, all_stops)
    )

    # For each bus line, connect ALL stops that can be reached without changing
    # (i.e., all stops on the same line are connected to each other)
    lines <- unique(data$Linie)

    for (line in lines) {
        # Get all stops on this line
        line_data <- data %>% filter(Linie == line)
        stops_on_line <- unique(c(line_data$Von_Haltestelle, line_data$Nach_Haltestelle))

        # Connect all pairs of stops on this line (both directions)
        for (stop1 in stops_on_line) {
            for (stop2 in stops_on_line) {
                if (stop1 != stop2) {
                    adj_matrix[stop1, stop2] <- 1
                }
            }
        }
    }

    return(adj_matrix)
}

calculate_network_metrics <- function(adj_matrix, name) {
    n_nodes <- nrow(adj_matrix)

    # Since we set both directions, matrix is symmetric
    # For undirected analysis, we work with the symmetric matrix directly
    # Ensure it's numeric (not logical)
    undirected_matrix <- 1 * (adj_matrix | t(adj_matrix))

    # Count edges
    n_directed_edges <- sum(adj_matrix)
    n_undirected_edges <- sum(undirected_matrix) / 2 # Each edge counted twice

    # Degree = number of unique neighbors (for undirected graph)
    degree <- rowSums(undirected_matrix)

    # Density: actual edges / possible edges (undirected)
    max_possible_edges <- n_nodes * (n_nodes - 1) / 2
    density <- n_undirected_edges / max_possible_edges

    # Convert to igraph object for advanced metrics
    g <- graph_from_adjacency_matrix(undirected_matrix, mode = "undirected")

    # Betweenness centrality (normalized)
    betweenness_vals <- betweenness(g, normalized = FALSE)
    avg_betweenness <- mean(betweenness_vals, na.rm = TRUE)
    max_betweenness <- max(betweenness_vals, na.rm = TRUE)
    if (is.na(max_betweenness) || is.infinite(max_betweenness)) max_betweenness <- 0
    if (is.na(avg_betweenness)) avg_betweenness <- 0

    # Top betweenness stops
    top_betweenness_stops <- names(sort(betweenness_vals, decreasing = TRUE)[1:min(5, length(betweenness_vals))])

    # Closeness centrality (only for largest connected component to avoid Inf)
    comp <- components(g)
    largest_cc <- which.max(comp$csize)
    nodes_in_lcc <- which(comp$membership == largest_cc)
    lcc_size <- length(nodes_in_lcc)

    if (lcc_size > 1) {
        g_lcc <- induced_subgraph(g, nodes_in_lcc)
        closeness_vals <- closeness(g_lcc, normalized = TRUE)
        avg_closeness <- mean(closeness_vals, na.rm = TRUE)
        if (is.na(avg_closeness) || is.infinite(avg_closeness)) avg_closeness <- 0

        # Network diameter and avg path length (only meaningful for connected component)
        diameter_val <- diameter(g_lcc)
        avg_path_length <- mean_distance(g_lcc, directed = FALSE)
        if (is.na(avg_path_length) || is.infinite(avg_path_length)) avg_path_length <- 0
    } else {
        avg_closeness <- 0
        diameter_val <- 0
        avg_path_length <- 0
    }

    # Clustering coefficient (global transitivity)
    clustering_coef <- transitivity(g, type = "global")
    if (is.na(clustering_coef)) clustering_coef <- 0

    # Number of connected components
    n_components <- comp$no

    metrics <- list(
        name = name,
        n_stops = n_nodes,
        n_undirected_connections = n_undirected_edges,
        density = density,

        # Degree centrality
        avg_degree = ifelse(length(degree) > 0, mean(degree), 0),
        max_degree = ifelse(length(degree) > 0, max(degree), 0),
        min_degree = ifelse(length(degree) > 0, min(degree), 0),
        isolated_stops = sum(degree == 0),
        hub_stops = names(sort(degree, decreasing = TRUE)[1:min(10, length(degree))]),

        # Betweenness centrality
        avg_betweenness = avg_betweenness,
        max_betweenness = max_betweenness,
        top_betweenness_stops = top_betweenness_stops,

        # Closeness centrality
        avg_closeness = avg_closeness,

        # Clustering & connectivity
        clustering_coefficient = clustering_coef,
        network_diameter = diameter_val,
        avg_shortest_path = avg_path_length,
        n_connected_components = n_components,
        largest_component_size = lcc_size,
        largest_component_pct = round(100 * lcc_size / n_nodes, 1)
    )

    return(metrics)
}

print_metrics <- function(metrics) {
    cat("\n")
    cat("=======================================================\n")
    cat(sprintf("  NETWORK METRICS: %s\n", metrics$name))
    cat("=======================================================\n")
    cat(sprintf("  Number of stops:              %d\n", metrics$n_stops))
    cat(sprintf("  Undirected connections:       %.0f\n", metrics$n_undirected_connections))
    cat(sprintf("  Network density:              %.4f\n", metrics$density))

    cat("\n  --- DEGREE CENTRALITY ---\n")
    cat(sprintf("  Average degree:               %.2f\n", metrics$avg_degree))
    cat(sprintf("  Max degree:                   %.0f\n", metrics$max_degree))
    cat(sprintf("  Min degree:                   %.0f\n", metrics$min_degree))
    cat(sprintf("  Isolated stops:               %.0f\n", metrics$isolated_stops))

    cat("\n  --- BETWEENNESS CENTRALITY ---\n")
    cat(sprintf("  Average betweenness:          %.2f\n", metrics$avg_betweenness))
    cat(sprintf("  Max betweenness:              %.2f\n", metrics$max_betweenness))
    cat("  Top 5 betweenness stops:\n")
    for (i in seq_along(metrics$top_betweenness_stops)) {
        cat(sprintf("    %d. %s\n", i, metrics$top_betweenness_stops[i]))
    }

    cat("\n  --- CLOSENESS CENTRALITY ---\n")
    cat(sprintf("  Average closeness:            %.4f\n", metrics$avg_closeness))

    cat("\n  --- CLUSTERING & CONNECTIVITY ---\n")
    cat(sprintf("  Clustering coefficient:       %.4f\n", metrics$clustering_coefficient))
    cat(sprintf("  Network diameter:             %.0f\n", metrics$network_diameter))
    cat(sprintf("  Avg shortest path length:     %.4f\n", metrics$avg_shortest_path))
    cat(sprintf("  Connected components:         %.0f\n", metrics$n_connected_components))
    cat(sprintf(
        "  Largest component:            %.0f stops (%.1f%%)\n",
        metrics$largest_component_size, metrics$largest_component_pct
    ))

    cat("\n  --- TOP 10 HUB STOPS (by degree) ---\n")
    for (i in seq_along(metrics$hub_stops)) {
        cat(sprintf("    %2d. %s\n", i, metrics$hub_stops[i]))
    }
    cat("=======================================================\n")
}

compare_networks <- function(metrics_2019, metrics_2025) {
    cat("\n")
    cat("=======================================================\n")
    cat("  COMPARISON: 2019 vs 2025\n")
    cat("=======================================================\n")

    # Calculate differences
    stop_diff <- metrics_2025$n_stops - metrics_2019$n_stops
    undirected_diff <- metrics_2025$n_undirected_connections - metrics_2019$n_undirected_connections
    density_diff <- metrics_2025$density - metrics_2019$density
    degree_diff <- metrics_2025$avg_degree - metrics_2019$avg_degree
    clustering_diff <- metrics_2025$clustering_coefficient - metrics_2019$clustering_coefficient
    diameter_diff <- metrics_2025$network_diameter - metrics_2019$network_diameter
    path_diff <- metrics_2025$avg_shortest_path - metrics_2019$avg_shortest_path

    cat(sprintf("  Stop difference:              %+.0f\n", stop_diff))
    cat(sprintf("  Connection difference:        %+.0f\n", undirected_diff))
    cat(sprintf("  Density difference:           %+.4f\n", density_diff))
    cat(sprintf("  Avg degree difference:        %+.2f\n", degree_diff))
    cat(sprintf("  Clustering coef difference:   %+.4f\n", clustering_diff))
    cat(sprintf("  Diameter difference:          %+.0f\n", diameter_diff))
    cat(sprintf("  Avg path length difference:   %+.4f\n", path_diff))

    cat("\n  INTERPRETATION:\n")
    if (undirected_diff > 0) {
        cat(sprintf("  - 2025 has %+.0f more reachable stop pairs\n", undirected_diff))
    } else if (undirected_diff < 0) {
        cat(sprintf("  - 2025 has %.0f fewer reachable stop pairs\n", undirected_diff))
    }

    if (density_diff > 0) {
        cat("  - 2025 is denser (better connected relative to size)\n")
    } else if (density_diff < 0) {
        cat("  - 2019 is denser (better connected relative to size)\n")
    }

    if (clustering_diff > 0) {
        cat("  - 2025 has higher clustering (more local connectivity)\n")
    } else if (clustering_diff < 0) {
        cat("  - 2019 has higher clustering (more local connectivity)\n")
    }

    if (path_diff < 0) {
        cat("  - 2025 has shorter average paths (faster to reach stops)\n")
    } else if (path_diff > 0) {
        cat("  - 2019 has shorter average paths (faster to reach stops)\n")
    }

    cat("=======================================================\n")
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

# Read data
data_2019 <- read_delim("2019fix.csv",
    delim = ";",
    show_col_types = FALSE
)

data_2025 <- read_delim("2025fix.csv",
    delim = ";",
    show_col_types = FALSE
)

# Create adjacency matrices
adj_2019 <- create_adjacency_matrix(data_2019)
adj_2025 <- create_adjacency_matrix(data_2025)

# Calculate metrics
metrics_2019 <- calculate_network_metrics(adj_2019, "Konstanz Bus Network 2019")
metrics_2025 <- calculate_network_metrics(adj_2025, "Konstanz Bus Network 2025")

# Print results
print_metrics(metrics_2019)
print_metrics(metrics_2025)
compare_networks(metrics_2019, metrics_2025)

# =============================================================================
# ADDITIONAL ANALYSIS: Stops only in one network
# =============================================================================

stops_2019 <- unique(c(data_2019$Von_Haltestelle, data_2019$Nach_Haltestelle))
stops_2025 <- unique(c(data_2025$Von_Haltestelle, data_2025$Nach_Haltestelle))

only_2019 <- setdiff(stops_2019, stops_2025)
only_2025 <- setdiff(stops_2025, stops_2019)

cat("\n")
cat("=======================================================\n")
cat("  STOPS ONLY IN 2019 (removed by 2025):\n")
cat("=======================================================\n")
if (length(only_2019) > 0) {
    for (stop in only_2019) {
        cat(sprintf("  - %s\n", stop))
    }
} else {
    cat("  (none)\n")
}

cat("\n")
cat("=======================================================\n")
cat("  STOPS ONLY IN 2025 (added since 2019):\n")
cat("=======================================================\n")
if (length(only_2025) > 0) {
    for (stop in only_2025) {
        cat(sprintf("  - %s\n", stop))
    }
} else {
    cat("  (none)\n")
}

# =============================================================================
# EXPORT ADJACENCY MATRICES (optional)
# =============================================================================

# Uncomment to save matrices as CSV:
# write.csv(adj_2019, "adjacency_matrix_2019.csv")
# write.csv(adj_2025, "adjacency_matrix_2025.csv")

# =============================================================================
# WINNERS/LOSERS ANALYSIS PER STOP
# How many connections has each stop gained/lost?
# =============================================================================

cat("\n")
cat("=======================================================\n")
cat("  WINNERS/LOSERS ANALYSIS PER STOP\n")
cat("=======================================================\n")

# Only stops that exist in both years
common_stops <- intersect(stops_2019, stops_2025)

# Calculate changes per stop
stop_changes <- data.frame(
    Stop = common_stops,
    Degree_2019 = NA_integer_,
    Degree_2025 = NA_integer_,
    Gained = NA_integer_, # New connections
    Lost = NA_integer_, # Lost connections
    Net_Change = NA_integer_
)

for (i in seq_along(common_stops)) {
    stop <- common_stops[i]

    # Neighbors in 2019 and 2025
    neighbors_2019 <- names(which(adj_2019[stop, ] == 1))
    neighbors_2025 <- names(which(adj_2025[stop, ] == 1))

    # Only common stops for fair comparison
    neighbors_2019_common <- intersect(neighbors_2019, common_stops)
    neighbors_2025_common <- intersect(neighbors_2025, common_stops)

    # Gained and lost connections
    gained <- setdiff(neighbors_2025_common, neighbors_2019_common)
    lost <- setdiff(neighbors_2019_common, neighbors_2025_common)

    stop_changes$Degree_2019[i] <- length(neighbors_2019)
    stop_changes$Degree_2025[i] <- length(neighbors_2025)
    stop_changes$Gained[i] <- length(gained)
    stop_changes$Lost[i] <- length(lost)
    stop_changes$Net_Change[i] <- length(gained) - length(lost)
}

# Sort by net change
stop_changes <- stop_changes %>% arrange(desc(Net_Change))

# Top 10 winners
cat("\n  TOP 10 WINNERS (most new connections):\n")
cat("  -------------------------------------------------------\n")
winners <- stop_changes %>%
    filter(Net_Change > 0) %>%
    head(10)
if (nrow(winners) > 0) {
    for (j in 1:nrow(winners)) {
        cat(sprintf(
            "    %2d. %s: +%d net (gained: %d, lost: %d)\n",
            j, winners$Stop[j], winners$Net_Change[j],
            winners$Gained[j], winners$Lost[j]
        ))
    }
} else {
    cat("    (no winners)\n")
}

# Top 10 losers
cat("\n  TOP 10 LOSERS (most lost connections):\n")
cat("  -------------------------------------------------------\n")
losers <- stop_changes %>%
    filter(Net_Change < 0) %>%
    arrange(Net_Change) %>%
    head(10)
if (nrow(losers) > 0) {
    for (j in 1:nrow(losers)) {
        cat(sprintf(
            "    %2d. %s: %d net (gained: %d, lost: %d)\n",
            j, losers$Stop[j], losers$Net_Change[j],
            losers$Gained[j], losers$Lost[j]
        ))
    }
} else {
    cat("    (no losers)\n")
}

# Unchanged
unchanged <- stop_changes %>% filter(Net_Change == 0)
cat(sprintf(
    "\n  Unchanged stops (net ±0): %d of %d\n",
    nrow(unchanged), length(common_stops)
))

# =============================================================================
# CONNECTION PAIR STATISTICS
# Which connections are new, which are removed?
# =============================================================================

cat("\n")
cat("=======================================================\n")
cat("  CONNECTION PAIR STATISTICS\n")
cat("=======================================================\n")

# Create sets of connection pairs (only for common stops)
get_connection_pairs <- function(adj_matrix, stops) {
    pairs <- c()
    for (i in 1:(length(stops) - 1)) {
        for (j in (i + 1):length(stops)) {
            stop1 <- stops[i]
            stop2 <- stops[j]
            if (adj_matrix[stop1, stop2] == 1) {
                # Sorted alphabetically for consistency
                pair <- paste(sort(c(stop1, stop2)), collapse = " <-> ")
                pairs <- c(pairs, pair)
            }
        }
    }
    return(pairs)
}

pairs_2019 <- get_connection_pairs(adj_2019, common_stops)
pairs_2025 <- get_connection_pairs(adj_2025, common_stops)

new_connections <- setdiff(pairs_2025, pairs_2019)
lost_connections <- setdiff(pairs_2019, pairs_2025)
unchanged_connections <- intersect(pairs_2019, pairs_2025)

cat(sprintf("  Connections 2019 (common stops): %d\n", length(pairs_2019)))
cat(sprintf("  Connections 2025 (common stops): %d\n", length(pairs_2025)))
cat(sprintf("  \n"))
cat(sprintf(
    "  New connections 2025:           %d (%+.1f%%)\n",
    length(new_connections),
    100 * length(new_connections) / max(length(pairs_2019), 1)
))
cat(sprintf(
    "  Removed connections:            %d (-%.1f%%)\n",
    length(lost_connections),
    100 * length(lost_connections) / max(length(pairs_2019), 1)
))
cat(sprintf(
    "  Unchanged connections:          %d (%.1f%%)\n",
    length(unchanged_connections),
    100 * length(unchanged_connections) / max(length(pairs_2019), 1)
))

# Show some new and lost connections
if (length(new_connections) > 0) {
    cat("\n  Examples NEW connections (max 10):\n")
    for (conn in head(new_connections, 10)) {
        cat(sprintf("    + %s\n", conn))
    }
}

if (length(lost_connections) > 0) {
    cat("\n  Examples REMOVED connections (max 10):\n")
    for (conn in head(lost_connections, 10)) {
        cat(sprintf("    - %s\n", conn))
    }
}

cat("=======================================================\n")

# Save stop changes as CSV
write.csv(stop_changes, "stop_connection_changes.csv", row.names = FALSE)
cat("\n  Stop changes saved to: stop_connection_changes.csv\n")

cat("\n\nAnalysis complete!\n")
