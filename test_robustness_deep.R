# Deep edge case tests for robustness analysis
library(tidyverse)
library(igraph)

cat("=== DEEP EDGE CASE TESTS ===\n\n")

# Recreate functions
create_adjacency_matrix <- function(data) {
    all_stops <- unique(c(data$Von_Haltestelle, data$Nach_Haltestelle))
    n_stops <- length(all_stops)
    adj_matrix <- matrix(0, nrow = n_stops, ncol = n_stops, dimnames = list(all_stops, all_stops))
    for (line in unique(data$Linie)) {
        line_data <- data %>% filter(Linie == line)
        stops_on_line <- unique(c(line_data$Von_Haltestelle, line_data$Nach_Haltestelle))
        for (stop1 in stops_on_line) {
            for (stop2 in stops_on_line) {
                if (stop1 != stop2) adj_matrix[stop1, stop2] <- 1
            }
        }
    }
    return(adj_matrix)
}

data_2019 <- read_delim("2019fix.csv", delim = ";", show_col_types = FALSE)
adj_2019 <- create_adjacency_matrix(data_2019)
g_2019 <- graph_from_adjacency_matrix(adj_2019, mode = "undirected")

# TEST 11: Check that removing top hubs actually disconnects the network eventually
cat("TEST 11: Progressive hub removal\n")
g_test <- g_2019
n_orig <- vcount(g_test)
for (i in 1:20) {
    if (vcount(g_test) == 0) {
        cat(sprintf("  Step %d: Network destroyed\n", i))
        break
    }
    degrees <- degree(g_test)
    top_node <- which.max(degrees)
    top_name <- V(g_test)$name[top_node]
    comp <- components(g_test)
    cat(sprintf(
        "  Step %d: Remove '%s' (deg=%d), LCC=%d/%d (%.1f%%), components=%d\n",
        i, top_name, degrees[top_node], max(comp$csize), n_orig,
        100 * max(comp$csize) / n_orig, comp$no
    ))
    g_test <- delete_vertices(g_test, top_node)
}

# TEST 12: Verify efficiency drops when we remove hubs
cat("\nTEST 12: Efficiency degradation\n")
calc_efficiency <- function(g) {
    if (vcount(g) == 0) {
        return(0)
    }
    d <- distances(g)
    diag(d) <- NA
    finite_mask <- is.finite(d)
    if (sum(finite_mask) > 0) {
        return(mean(1 / d[finite_mask]))
    } else {
        return(0)
    }
}

g_test <- g_2019
baseline_eff <- calc_efficiency(g_test)
cat(sprintf("  Baseline efficiency: %.4f\n", baseline_eff))

for (pct in c(5, 10, 15, 20, 25, 30)) {
    g_test <- g_2019
    n_to_remove <- round(vcount(g_2019) * pct / 100)
    for (i in 1:n_to_remove) {
        if (vcount(g_test) == 0) break
        degrees <- degree(g_test)
        top_node <- which.max(degrees)
        g_test <- delete_vertices(g_test, top_node)
    }
    eff <- calc_efficiency(g_test)
    cat(sprintf(
        "  After removing %d%% (%d nodes): efficiency=%.4f (%.1f%% of baseline)\n",
        pct, n_to_remove, eff, 100 * eff / baseline_eff
    ))
}

# TEST 13: Verify random removal is indeed random (check variance)
cat("\nTEST 13: Random removal variance\n")
set.seed(42)
efficiencies <- c()
for (sim in 1:20) {
    g_test <- g_2019
    n_to_remove <- round(vcount(g_2019) * 0.2)
    nodes_to_remove <- sample(1:vcount(g_test), n_to_remove)
    # Remove one by one (indices shift!)
    for (i in 1:n_to_remove) {
        if (vcount(g_test) == 0) break
        # Sample new index each time since graph shrinks
        if (vcount(g_test) > 0) {
            idx <- sample(1:vcount(g_test), 1)
            g_test <- delete_vertices(g_test, idx)
        }
    }
    eff <- calc_efficiency(g_test)
    efficiencies <- c(efficiencies, eff)
}
cat(sprintf("  Mean efficiency after 20%% random removal: %.4f\n", mean(efficiencies)))
cat(sprintf("  Std deviation: %.6f\n", sd(efficiencies)))
cat(sprintf("  Range: %.4f - %.4f\n", min(efficiencies), max(efficiencies)))

# TEST 14: Check components behavior
cat("\nTEST 14: Network connectivity\n")
comp <- components(g_2019)
cat(sprintf("  2019 network is fully connected: %s\n", comp$no == 1))
cat(sprintf("  Number of components: %d\n", comp$no))

data_2025 <- read_delim("2025fix.csv", delim = ";", show_col_types = FALSE)
adj_2025 <- create_adjacency_matrix(data_2025)
g_2025 <- graph_from_adjacency_matrix(adj_2025, mode = "undirected")
comp <- components(g_2025)
cat(sprintf("  2025 network is fully connected: %s\n", comp$no == 1))
cat(sprintf("  Number of components: %d\n", comp$no))

# TEST 15: Verify the 50% threshold calculation
cat("\nTEST 15: 50% LCC threshold verification\n")
g_test <- g_2019
n_orig <- vcount(g_test)
step <- 0
nodes_per_step <- floor(n_orig / 15)
cat(sprintf("  Nodes per step: %d (total nodes: %d)\n", nodes_per_step, n_orig))

while (vcount(g_test) > 0) {
    comp <- components(g_test)
    lcc_pct <- 100 * max(comp$csize) / n_orig
    removed_pct <- 100 * (n_orig - vcount(g_test)) / n_orig

    if (lcc_pct < 50 || step == 0 || step %% 2 == 0) {
        cat(sprintf("  Step %d: Removed %.0f%%, LCC=%.1f%%\n", step, removed_pct, lcc_pct))
    }

    if (lcc_pct < 50) {
        cat(sprintf("  --> First drops below 50%% at %.0f%% removal\n", removed_pct))
        break
    }

    # Remove highest degree nodes
    degrees <- degree(g_test)
    to_remove <- min(nodes_per_step, vcount(g_test))
    top_nodes <- order(degrees, decreasing = TRUE)[1:to_remove]
    g_test <- delete_vertices(g_test, top_nodes)
    step <- step + 1
}

# TEST 16: Check for any NA/NaN in critical calculations
cat("\nTEST 16: NA/NaN checks\n")
distances_2019 <- distances(g_2019)
distances_2025 <- distances(g_2025)
cat(sprintf(
    "  2019 distances - any NA: %s, any NaN: %s\n",
    any(is.na(distances_2019)), any(is.nan(distances_2019))
))
cat(sprintf(
    "  2025 distances - any NA: %s, any NaN: %s\n",
    any(is.na(distances_2025)), any(is.nan(distances_2025))
))

cat("\n=== ALL DEEP TESTS COMPLETE ===\n")
