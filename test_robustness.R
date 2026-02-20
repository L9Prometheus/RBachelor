# Systematic tests for robustness analysis
library(tidyverse)
library(igraph)

source("robustness_analysis.R", local = FALSE, echo = FALSE)

# Actually, let's just test the functions directly
cat("=== ROBUSTNESS ANALYSIS BUG TESTS ===\n\n")

# Load data
data_2019 <- read_delim("2019fix.csv", delim = ";", show_col_types = FALSE)
data_2025 <- read_delim("2025fix.csv", delim = ";", show_col_types = FALSE)

# Recreate functions for testing
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

adj_2019 <- create_adjacency_matrix(data_2019)
adj_2025 <- create_adjacency_matrix(data_2025)

# TEST 1: Matrix symmetry
cat("TEST 1: Matrix symmetry\n")
cat("  2019 symmetric:", all(adj_2019 == t(adj_2019)), "\n")
cat("  2025 symmetric:", all(adj_2025 == t(adj_2025)), "\n")

# TEST 2: Diagonal all zeros
cat("\nTEST 2: Diagonal (self-loops)\n")
cat("  2019 diagonal sum:", sum(diag(adj_2019)), "(should be 0)\n")
cat("  2025 diagonal sum:", sum(diag(adj_2025)), "(should be 0)\n")

# TEST 3: Graph creation
g_2019 <- graph_from_adjacency_matrix(adj_2019, mode = "undirected")
g_2025 <- graph_from_adjacency_matrix(adj_2025, mode = "undirected")

cat("\nTEST 3: Graph node count\n")
cat("  2019 nodes:", vcount(g_2019), "expected:", length(unique(c(data_2019$Von_Haltestelle, data_2019$Nach_Haltestelle))), "\n")
cat("  2025 nodes:", vcount(g_2025), "expected:", length(unique(c(data_2025$Von_Haltestelle, data_2025$Nach_Haltestelle))), "\n")

# TEST 4: Distance matrix
cat("\nTEST 4: Distance matrix sanity\n")
dist_2019 <- distances(g_2019)
dist_2025 <- distances(g_2025)
cat("  2019 - min dist (excl diag):", min(dist_2019[row(dist_2019) != col(dist_2019)]), "\n")
cat("  2019 - max dist:", max(dist_2019[is.finite(dist_2019)]), "\n")
cat("  2019 - any Inf:", any(is.infinite(dist_2019[row(dist_2019) != col(dist_2019)])), "\n")
cat("  2025 - min dist (excl diag):", min(dist_2025[row(dist_2025) != col(dist_2025)]), "\n")
cat("  2025 - max dist:", max(dist_2025[is.finite(dist_2025)]), "\n")
cat("  2025 - any Inf:", any(is.infinite(dist_2025[row(dist_2025) != col(dist_2025)])), "\n")

# TEST 5: Efficiency calculation
cat("\nTEST 5: Efficiency calculation\n")
calc_efficiency <- function(g) {
    d <- distances(g)
    diag(d) <- NA
    finite_mask <- is.finite(d)
    if (sum(finite_mask) > 0) {
        return(mean(1 / d[finite_mask]))
    } else {
        return(0)
    }
}
eff_2019 <- calc_efficiency(g_2019)
eff_2025 <- calc_efficiency(g_2025)
cat("  2019 efficiency:", eff_2019, "\n")
cat("  2025 efficiency:", eff_2025, "\n")
cat("  Valid range (0-1):", eff_2019 >= 0 && eff_2019 <= 1 && eff_2025 >= 0 && eff_2025 <= 1, "\n")

# TEST 6: Check if deleting a node works correctly
cat("\nTEST 6: Node deletion\n")
g_test <- g_2019
original_nodes <- vcount(g_test)
g_reduced <- delete_vertices(g_test, 1)
cat("  Original nodes:", original_nodes, "\n")
cat("  After deletion:", vcount(g_reduced), "\n")
cat("  Difference:", original_nodes - vcount(g_reduced), "(should be 1)\n")

# TEST 7: LCC calculation after deletion
cat("\nTEST 7: LCC calculation\n")
comp_orig <- components(g_2019)
comp_reduced <- components(g_reduced)
cat("  Original LCC size:", max(comp_orig$csize), "\n")
cat("  After removing 1 node LCC:", max(comp_reduced$csize), "\n")
cat("  Network stays connected:", comp_reduced$no == 1, "\n")

# TEST 8: Targeted attack - verify degree ordering
cat("\nTEST 8: Targeted attack degree ordering\n")
degrees <- degree(g_2019)
top_5 <- names(sort(degrees, decreasing = TRUE)[1:5])
cat("  Top 5 highest degree nodes:\n")
for (i in 1:5) {
    cat(sprintf("    %d. %s (degree: %d)\n", i, top_5[i], degrees[top_5[i]]))
}

# TEST 9: Verify LCC fraction calculation
cat("\nTEST 9: LCC fraction calculation\n")
n_nodes <- vcount(g_2019)
# Remove first 7 nodes (about 6.5%)
g_partial <- delete_vertices(g_2019, 1:7)
comp_partial <- components(g_partial)
lcc_size <- max(comp_partial$csize)
lcc_frac_of_remaining <- lcc_size / vcount(g_partial)
lcc_frac_of_original <- lcc_size / n_nodes
cat("  Removed 7 of", n_nodes, "nodes\n")
cat("  Remaining nodes:", vcount(g_partial), "\n")
cat("  LCC size:", lcc_size, "\n")
cat("  LCC as % of remaining:", round(100 * lcc_frac_of_remaining, 1), "%\n")
cat("  LCC as % of original:", round(100 * lcc_frac_of_original, 1), "%\n")
cat("  (We use % of original in the analysis)\n")

# TEST 10: Edge count verification
cat("\nTEST 10: Edge count\n")
cat("  2019 edges:", ecount(g_2019), "\n")
cat("  2025 edges:", ecount(g_2025), "\n")
cat("  Expected 2019 (from adj matrix):", sum(adj_2019) / 2, "\n")
cat("  Expected 2025 (from adj matrix):", sum(adj_2025) / 2, "\n")

cat("\n=== ALL TESTS COMPLETE ===\n")
