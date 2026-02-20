# Bus Network Robustness Analysis
# Tests how resilient the 2019 and 2025 networks are to node/edge removal

library(tidyverse)
library(igraph)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

create_adjacency_matrix <- function(data) {
    all_stops <- unique(c(data$Von_Haltestelle, data$Nach_Haltestelle))
    n_stops <- length(all_stops)
    adj_matrix <- matrix(0,
        nrow = n_stops, ncol = n_stops,
        dimnames = list(all_stops, all_stops)
    )
    lines <- unique(data$Linie)
    for (line in lines) {
        line_data <- data %>% filter(Linie == line)
        stops_on_line <- unique(c(line_data$Von_Haltestelle, line_data$Nach_Haltestelle))
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

create_graph <- function(adj_matrix) {
    undirected_matrix <- 1 * (adj_matrix | t(adj_matrix))
    g <- graph_from_adjacency_matrix(undirected_matrix, mode = "undirected")
    return(g)
}

# =============================================================================
# ROBUSTNESS METRICS
# =============================================================================

# Calculate network metrics for a given graph
get_network_health <- function(g) {
    if (vcount(g) == 0) {
        return(list(
            largest_component_fraction = 0,
            avg_path_length = Inf,
            efficiency = 0,
            n_components = 0
        ))
    }

    comp <- components(g)
    lcc_size <- max(comp$csize)
    lcc_fraction <- lcc_size / vcount(g)

    # Get largest connected component for path calculations
    if (lcc_size > 1) {
        lcc_nodes <- which(comp$membership == which.max(comp$csize))
        g_lcc <- induced_subgraph(g, lcc_nodes)
        avg_path <- mean_distance(g_lcc, directed = FALSE)
        if (is.na(avg_path) || is.infinite(avg_path)) avg_path <- Inf

        # Global efficiency: mean(1 / d_ij) over all pairs i != j
        # Disconnected pairs have d_ij = Inf and contribute 0 because 1/Inf = 0
        d <- distances(g)
        d <- d[upper.tri(d)] # unordered pairs only, no diagonal
        efficiency <- mean(1 / d) # includes zeros from Inf distances
        if (is.na(efficiency) || is.infinite(efficiency)) efficiency <- 0
    } else {
        avg_path <- Inf
        efficiency <- 0
    }

    return(list(
        largest_component_fraction = lcc_fraction,
        avg_path_length = avg_path,
        efficiency = efficiency,
        n_components = comp$no
    ))
}

# =============================================================================
# TARGETED ATTACK SIMULATION
# Remove nodes by highest degree (simulates attack on hubs)
# =============================================================================

simulate_targeted_attack <- function(g, steps = 20) {
    n_nodes <- vcount(g)
    nodes_to_remove_per_step <- max(1, floor(n_nodes / steps))

    results <- data.frame(
        fraction_removed = numeric(),
        lcc_fraction = numeric(),
        efficiency = numeric(),
        n_components = numeric()
    )

    g_current <- g
    removed <- 0

    for (step in 0:steps) {
        if (vcount(g_current) == 0) {
            # Network completely destroyed
            results <- rbind(results, data.frame(
                fraction_removed = removed / n_nodes,
                lcc_fraction = 0,
                efficiency = 0,
                n_components = 0
            ))
            break
        }

        health <- get_network_health(g_current)
        # LCC as fraction of ORIGINAL network size
        current_lcc_size <- health$largest_component_fraction * vcount(g_current)
        results <- rbind(results, data.frame(
            fraction_removed = removed / n_nodes,
            lcc_fraction = current_lcc_size / n_nodes,
            efficiency = health$efficiency,
            n_components = health$n_components
        ))

        if (step == steps) break

        # Remove highest degree nodes
        degrees <- degree(g_current)
        if (length(degrees) == 0) break

        # Sort by degree, remove top nodes
        to_remove <- min(nodes_to_remove_per_step, vcount(g_current))
        top_nodes <- order(degrees, decreasing = TRUE)[1:to_remove]
        g_current <- delete_vertices(g_current, top_nodes)
        removed <- removed + to_remove
    }

    return(results)
}

# =============================================================================
# RANDOM FAILURE SIMULATION
# Remove nodes randomly (simulates random failures)
# =============================================================================

simulate_random_failure <- function(g, steps = 20, n_simulations = 10) {
    n_nodes <- vcount(g)
    nodes_to_remove_per_step <- max(1, floor(n_nodes / steps))

    all_results <- list()

    for (sim in 1:n_simulations) {
        results <- data.frame(
            fraction_removed = numeric(),
            lcc_fraction = numeric(),
            efficiency = numeric()
        )

        g_current <- g
        removed <- 0

        for (step in 0:steps) {
            if (vcount(g_current) == 0) {
                # Network completely destroyed
                results <- rbind(results, data.frame(
                    fraction_removed = removed / n_nodes,
                    lcc_fraction = 0,
                    efficiency = 0
                ))
                break
            }

            health <- get_network_health(g_current)
            # LCC as fraction of ORIGINAL network size
            current_lcc_size <- health$largest_component_fraction * vcount(g_current)
            results <- rbind(results, data.frame(
                fraction_removed = removed / n_nodes,
                lcc_fraction = current_lcc_size / n_nodes,
                efficiency = health$efficiency
            ))

            if (step == steps) break

            # Remove random nodes
            to_remove <- min(nodes_to_remove_per_step, vcount(g_current))
            random_nodes <- sample(1:vcount(g_current), to_remove)
            g_current <- delete_vertices(g_current, random_nodes)
            removed <- removed + to_remove
        }

        all_results[[sim]] <- results
    }

    # Average across simulations - ensure same length by padding
    max_rows <- max(sapply(all_results, nrow))

    # Pad shorter results with final values (network destroyed)
    for (sim in 1:n_simulations) {
        current_rows <- nrow(all_results[[sim]])
        if (current_rows < max_rows) {
            last_row <- all_results[[sim]][current_rows, ]
            for (r in (current_rows + 1):max_rows) {
                # Extrapolate fraction_removed, keep lcc and efficiency at 0
                new_row <- data.frame(
                    fraction_removed = min(1, last_row$fraction_removed +
                        (r - current_rows) * (1 / steps)),
                    lcc_fraction = 0,
                    efficiency = 0
                )
                all_results[[sim]] <- rbind(all_results[[sim]], new_row)
            }
        }
    }

    # Now average across simulations
    avg_results <- all_results[[1]]
    avg_results$lcc_fraction <- 0
    avg_results$efficiency <- 0

    for (sim in 1:n_simulations) {
        avg_results$lcc_fraction <- avg_results$lcc_fraction + all_results[[sim]]$lcc_fraction
        avg_results$efficiency <- avg_results$efficiency + all_results[[sim]]$efficiency
    }
    avg_results$lcc_fraction <- avg_results$lcc_fraction / n_simulations
    avg_results$efficiency <- avg_results$efficiency / n_simulations

    return(avg_results)
}

# =============================================================================
# CRITICAL NODE ANALYSIS
# Which single node removal causes the biggest damage?
# =============================================================================

find_critical_nodes <- function(g, top_n = 10) {
    baseline <- get_network_health(g)
    baseline_efficiency <- baseline$efficiency

    node_names <- V(g)$name
    n_nodes <- length(node_names)

    impacts <- data.frame(
        node = character(),
        efficiency_loss = numeric(),
        efficiency_loss_pct = numeric(),
        components_added = numeric()
    )

    for (i in 1:n_nodes) {
        g_reduced <- delete_vertices(g, i)
        health <- get_network_health(g_reduced)

        eff_loss <- baseline_efficiency - health$efficiency
        eff_loss_pct <- 100 * eff_loss / baseline_efficiency
        comp_added <- health$n_components - baseline$n_components

        impacts <- rbind(impacts, data.frame(
            node = node_names[i],
            efficiency_loss = eff_loss,
            efficiency_loss_pct = eff_loss_pct,
            components_added = comp_added
        ))
    }

    # Sort by efficiency loss
    impacts <- impacts %>% arrange(desc(efficiency_loss_pct))

    return(head(impacts, top_n))
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

cat("Loading data...\n")

data_2019 <- read_delim("2019fix.csv", delim = ";", show_col_types = FALSE)
data_2025 <- read_delim("2025fix.csv", delim = ";", show_col_types = FALSE)

adj_2019 <- create_adjacency_matrix(data_2019)
adj_2025 <- create_adjacency_matrix(data_2025)

g_2019 <- create_graph(adj_2019)
g_2025 <- create_graph(adj_2025)

cat("\n=======================================================\n")
cat("  ROBUSTNESS ANALYSIS: KONSTANZ BUS NETWORK\n")
cat("=======================================================\n")

# -----------------------------------------------------------------------------
# 1. CRITICAL NODE ANALYSIS
# -----------------------------------------------------------------------------

cat("\n--- CRITICAL NODE ANALYSIS ---\n")
cat("Which stops, if removed, cause the biggest efficiency loss?\n\n")

cat("2019 - TOP 10 CRITICAL STOPS:\n")
critical_2019 <- find_critical_nodes(g_2019, 10)
for (i in 1:nrow(critical_2019)) {
    cat(sprintf(
        "  %2d. %-35s  -%.2f%% efficiency\n",
        i, critical_2019$node[i], critical_2019$efficiency_loss_pct[i]
    ))
}

cat("\n2025 - TOP 10 CRITICAL STOPS:\n")
critical_2025 <- find_critical_nodes(g_2025, 10)
for (i in 1:nrow(critical_2025)) {
    cat(sprintf(
        "  %2d. %-35s  -%.2f%% efficiency\n",
        i, critical_2025$node[i], critical_2025$efficiency_loss_pct[i]
    ))
}

# Compare criticality
cat("\n--- CRITICALITY COMPARISON ---\n")
max_crit_2019 <- max(critical_2019$efficiency_loss_pct, na.rm = TRUE)
max_crit_2025 <- max(critical_2025$efficiency_loss_pct, na.rm = TRUE)
if (is.na(max_crit_2019) || is.infinite(max_crit_2019)) max_crit_2019 <- 0
if (is.na(max_crit_2025) || is.infinite(max_crit_2025)) max_crit_2025 <- 0
cat(sprintf("  Max single-node impact 2019: %.2f%%\n", max_crit_2019))
cat(sprintf("  Max single-node impact 2025: %.2f%%\n", max_crit_2025))

if (max_crit_2025 < max_crit_2019) {
    cat("  -> 2025 is MORE ROBUST (less dependent on single nodes)\n")
} else {
    cat("  -> 2019 is MORE ROBUST (less dependent on single nodes)\n")
}

# -----------------------------------------------------------------------------
# 2. TARGETED ATTACK SIMULATION
# -----------------------------------------------------------------------------

cat("\n--- TARGETED ATTACK SIMULATION ---\n")
cat("Removing nodes by highest degree (attack on hubs)...\n")

attack_2019 <- simulate_targeted_attack(g_2019, steps = 15)
attack_2025 <- simulate_targeted_attack(g_2025, steps = 15)

# Find when LCC drops below 50%
threshold_2019 <- attack_2019 %>%
    filter(lcc_fraction < 0.5) %>%
    slice(1) %>%
    pull(fraction_removed)
threshold_2025 <- attack_2025 %>%
    filter(lcc_fraction < 0.5) %>%
    slice(1) %>%
    pull(fraction_removed)

if (length(threshold_2019) == 0) threshold_2019 <- 1.0
if (length(threshold_2025) == 0) threshold_2025 <- 1.0

cat(sprintf(
    "\n  2019: Network fragments (LCC < 50%%) after removing %.0f%% of nodes\n",
    100 * threshold_2019
))
cat(sprintf(
    "  2025: Network fragments (LCC < 50%%) after removing %.0f%% of nodes\n",
    100 * threshold_2025
))

# -----------------------------------------------------------------------------
# 3. RANDOM FAILURE SIMULATION
# -----------------------------------------------------------------------------

cat("\n--- RANDOM FAILURE SIMULATION ---\n")
cat("Removing nodes randomly (averaged over 10 simulations)...\n")

set.seed(42)
random_2019 <- simulate_random_failure(g_2019, steps = 15, n_simulations = 10)
random_2025 <- simulate_random_failure(g_2025, steps = 15, n_simulations = 10)

# Compare efficiency at 20% removal (find closest value >= 0.2)
eff_20_2019_df <- random_2019 %>% filter(fraction_removed >= 0.19)
eff_20_2025_df <- random_2025 %>% filter(fraction_removed >= 0.19)

if (nrow(eff_20_2019_df) > 0) {
    eff_20_2019 <- eff_20_2019_df$efficiency[1]
} else {
    eff_20_2019 <- tail(random_2019$efficiency, 1)
}

if (nrow(eff_20_2025_df) > 0) {
    eff_20_2025 <- eff_20_2025_df$efficiency[1]
} else {
    eff_20_2025 <- tail(random_2025$efficiency, 1)
}

# Get baseline efficiency
baseline_2019 <- random_2019$efficiency[1]
baseline_2025 <- random_2025$efficiency[1]

# Handle edge cases
if (is.na(baseline_2019) || baseline_2019 == 0) baseline_2019 <- 1
if (is.na(baseline_2025) || baseline_2025 == 0) baseline_2025 <- 1

cat(sprintf("\n  Efficiency at 20%% random removal:\n"))
cat(sprintf(
    "    2019: %.4f (%.1f%% of original)\n",
    eff_20_2019, 100 * eff_20_2019 / baseline_2019
))
cat(sprintf(
    "    2025: %.4f (%.1f%% of original)\n",
    eff_20_2025, 100 * eff_20_2025 / baseline_2025
))

# -----------------------------------------------------------------------------
# 4. SUMMARY
# -----------------------------------------------------------------------------

cat("\n=======================================================\n")
cat("  ROBUSTNESS SUMMARY\n")
cat("=======================================================\n")

cat("\n  Metric                              2019      2025\n")
cat("  -------------------------------------------------\n")
cat(sprintf(
    "  Baseline efficiency:               %.4f    %.4f\n",
    baseline_2019, baseline_2025
))
cat(sprintf(
    "  Max single-node impact:            %.2f%%    %.2f%%\n",
    max_crit_2019, max_crit_2025
))
cat(sprintf(
    "  Attack threshold (50%% LCC):        %.0f%%      %.0f%%\n",
    100 * threshold_2019, 100 * threshold_2025
))
cat(sprintf(
    "  Efficiency after 20%% random loss:  %.1f%%     %.1f%%\n",
    100 * eff_20_2019 / baseline_2019, 100 * eff_20_2025 / baseline_2025
))

cat("\n  INTERPRETATION:\n")
if (max_crit_2025 < max_crit_2019) {
    cat("  - 2025 has lower single-node vulnerability\n")
}
if (threshold_2025 > threshold_2019) {
    cat("  - 2025 survives targeted attacks longer\n")
}
if (eff_20_2025 / baseline_2025 > eff_20_2019 / baseline_2019) {
    cat("  - 2025 retains more efficiency under random failures\n")
}

cat("=======================================================\n")

# -----------------------------------------------------------------------------
# 5. SAVE PLOTS
# -----------------------------------------------------------------------------

cat("\nGenerating plots...\n")

# Combine data for plotting
attack_2019$year <- "2019"
attack_2025$year <- "2025"
attack_combined <- rbind(attack_2019, attack_2025)

random_2019$year <- "2019"
random_2025$year <- "2025"
random_combined <- rbind(random_2019, random_2025)

# Plot 1: Targeted Attack
# Find exact threshold points for annotation
threshold_point_2019 <- attack_2019 %>%
    filter(lcc_fraction < 0.5) %>%
    slice(1)
threshold_point_2025 <- attack_2025 %>%
    filter(lcc_fraction < 0.5) %>%
    slice(1)

p1 <- ggplot(attack_combined, aes(x = fraction_removed * 100, y = lcc_fraction * 100, color = year)) +
    # Step plot for discrete removal process
    geom_step(linewidth = 1.2, direction = "hv") +
    geom_point(size = 2.5) +
    # 50% threshold line
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    # Vertical lines at fragmentation thresholds
    geom_vline(
        xintercept = threshold_point_2019$fraction_removed * 100,
        linetype = "dotted", color = "#E69F00", linewidth = 0.8
    ) +
    geom_vline(
        xintercept = threshold_point_2025$fraction_removed * 100,
        linetype = "dotted", color = "#0072B2", linewidth = 0.8
    ) +
    # Annotate threshold points
    annotate("text",
        x = threshold_point_2019$fraction_removed * 100 + 2, y = 55,
        label = paste0(round(threshold_point_2019$fraction_removed * 100), "%"),
        color = "#E69F00", size = 3.5, hjust = 0
    ) +
    annotate("text",
        x = threshold_point_2025$fraction_removed * 100 + 2, y = 55,
        label = paste0(round(threshold_point_2025$fraction_removed * 100), "%"),
        color = "#0072B2", size = 3.5, hjust = 0
    ) +
    labs(
        title = "Targeted Attack Robustness",
        subtitle = "Sequential removal of highest-degree nodes (hub attack)",
        x = "Share of Nodes Removed (%)",
        y = "Largest Connected Component Size (% of nodes)",
        color = "Year"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90")
    ) +
    scale_color_manual(
        values = c("2019" = "#E69F00", "2025" = "#0072B2"),
        labels = c("2019" = "2019", "2025" = "2025")
    ) +
    scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100))

ggsave("robustness_targeted_attack.png", p1, width = 9, height = 6, dpi = 150)

# Plot 2: Random Failure
p2 <- ggplot(random_combined, aes(x = fraction_removed * 100, y = efficiency, color = year)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    labs(
        title = "Random Failure Robustness",
        subtitle = "Average efficiency under random node removal (10 simulations)",
        x = "% of Nodes Removed",
        y = "Global Efficiency",
        color = "Year"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("2019" = "#E69F00", "2025" = "#0072B2"))

ggsave("robustness_random_failure.png", p2, width = 8, height = 6, dpi = 150)

# Plot 3: Critical Nodes Comparison (only stops in both years' top 10)
critical_2019$year <- "2019"
critical_2025$year <- "2025"

# Find stops that appear in both top 10 lists
common_critical_stops <- intersect(critical_2019$node, critical_2025$node)

# Filter to only common stops and combine
critical_2019_common <- critical_2019 %>%
    filter(node %in% common_critical_stops)
critical_2025_common <- critical_2025 %>%
    filter(node %in% common_critical_stops)

critical_combined <- rbind(critical_2019_common, critical_2025_common)

# Order stops by average criticality (most critical at top)
stop_order <- critical_combined %>%
    group_by(node) %>%
    summarise(avg_loss = mean(efficiency_loss_pct)) %>%
    arrange(desc(avg_loss)) %>%
    pull(node)

critical_combined$node <- factor(critical_combined$node, levels = rev(stop_order))

p3 <- ggplot(critical_combined, aes(x = efficiency_loss_pct, y = node, fill = year)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.2f%%", efficiency_loss_pct)),
        position = position_dodge(width = 0.8),
        hjust = -0.1, size = 3
    ) +
    labs(
        title = "Critical Node Impact Comparison",
        subtitle = paste0(
            "Efficiency loss when single stop is removed (",
            length(common_critical_stops), " stops in both top 10)"
        ),
        x = "Efficiency Loss (%)",
        y = NULL,
        fill = "Year"
    ) +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10)
    ) +
    scale_fill_manual(values = c("2019" = "#E69F00", "2025" = "#0072B2")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("robustness_critical_nodes.png", p3, width = 9, height = 6, dpi = 150)

# Print info about stops only in one year's top 10
only_2019 <- setdiff(critical_2019$node, critical_2025$node)
only_2025 <- setdiff(critical_2025$node, critical_2019$node)

if (length(only_2019) > 0 || length(only_2025) > 0) {
    cat("\n  Note: Some stops appear only in one year's top 10:\n")
    if (length(only_2019) > 0) {
        cat("    Only in 2019:", paste(only_2019, collapse = ", "), "\n")
    }
    if (length(only_2025) > 0) {
        cat("    Only in 2025:", paste(only_2025, collapse = ", "), "\n")
    }
}

cat("\nPlots saved:\n")
cat("  - robustness_targeted_attack.png\n")
cat("  - robustness_random_failure.png\n")
cat("  - robustness_critical_nodes.png\n")

cat("\n\nRobustness analysis complete!\n")
