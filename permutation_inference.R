# =============================================================================
# Permutation-based Inference for Network Centrality Analysis
# =============================================================================
#
# Tests whether observed associations between network metrics are stronger
# than expected under random assignment of centrality values.
#
# Methodology:
# - All predictors (degree, closeness, n_lines) are z-standardized
# - Betweenness is normalized by igraph, then z-standardized
# - P-values: p = (1 + count) / (B + 1), one-sided tests
# - B = 10,000 permutations
#
# =============================================================================

library(tidyverse)
library(igraph)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Create adjacency matrix from edge data
#' All stops on the same line are connected to each other
create_adjacency_matrix <- function(data) {
    all_stops <- unique(c(data$Von_Haltestelle, data$Nach_Haltestelle))
    n_stops <- length(all_stops)

    # Initialize empty matrix
    adj_matrix <- matrix(0,
        nrow = n_stops, ncol = n_stops,
        dimnames = list(all_stops, all_stops)
    )

    # For each line, connect all stops on that line
    for (line in unique(data$Linie)) {
        line_stops <- unique(c(
            data$Von_Haltestelle[data$Linie == line],
            data$Nach_Haltestelle[data$Linie == line]
        ))
        # Make all stops on this line adjacent to each other
        for (s1 in line_stops) {
            for (s2 in line_stops) {
                if (s1 != s2) {
                    adj_matrix[s1, s2] <- 1
                }
            }
        }
    }
    return(adj_matrix)
}

#' Create undirected graph from adjacency matrix
create_graph <- function(adj_matrix) {
    # Make symmetric (undirected)
    sym_matrix <- pmax(adj_matrix, t(adj_matrix))
    g <- graph_from_adjacency_matrix(sym_matrix, mode = "undirected", diag = FALSE)
    return(g)
}

#' Calculate centrality measures for all nodes
calculate_centralities <- function(g, data) {
    stops <- V(g)$name
    n_stops <- length(stops)

    # Degree centrality
    deg <- degree(g)

    # Betweenness centrality (normalized)
    bet <- betweenness(g, normalized = TRUE)

    # Closeness centrality - only meaningful for connected components
    # Use largest connected component
    comp <- components(g)
    largest_comp <- which.max(comp$csize)
    in_largest <- comp$membership == largest_comp

    close <- rep(NA_real_, n_stops)
    names(close) <- stops
    if (sum(in_largest) > 1) {
        g_sub <- induced_subgraph(g, which(in_largest))
        close[in_largest] <- closeness(g_sub, normalized = TRUE)
    }

    # Number of lines per stop
    stop_lines <- data %>%
        pivot_longer(
            cols = c(Von_Haltestelle, Nach_Haltestelle),
            values_to = "stop", names_to = "type"
        ) %>%
        group_by(stop) %>%
        summarise(n_lines = n_distinct(Linie), .groups = "drop")

    # Build result data frame
    result <- data.frame(
        stop = stops,
        degree = as.numeric(deg[stops]),
        betweenness = as.numeric(bet[stops]),
        closeness = as.numeric(close[stops]),
        stringsAsFactors = FALSE
    )

    # Add n_lines by matching
    result$n_lines <- stop_lines$n_lines[match(result$stop, stop_lines$stop)]

    return(result)
}

#' Run permutation test for association between predictor and outcome
#' @param df Data frame with columns for predictor and outcome
#' @param predictor Name of predictor column
#' @param outcome Name of outcome column
#' @param n_perms Number of permutations
#' @param seed Random seed
#' @return List with observed coefficient, R-squared, p-value, and permutation distribution
permutation_test <- function(df, predictor, outcome, n_perms = 10000, seed = 42) {
    set.seed(seed)

    # Remove NA values
    df_clean <- df[complete.cases(df[, c(predictor, outcome)]), ]
    n <- nrow(df_clean)

    if (n < 3) {
        warning("Not enough observations for permutation test")
        return(list(
            observed_coef = NA, r_squared = NA, p_value = NA,
            perm_distribution = NA, n_obs = n
        ))
    }

    # Z-standardize both variables
    x <- df_clean[[predictor]]
    y <- df_clean[[outcome]]

    x_std <- (x - mean(x)) / sd(x)
    y_std <- (y - mean(y)) / sd(y)

    # Observed coefficient (equals correlation for standardized variables)
    observed_coef <- sum(x_std * y_std) / (n - 1)
    observed_r_sq <- observed_coef^2

    # Permutation distribution
    perm_coefs <- numeric(n_perms)
    for (i in seq_len(n_perms)) {
        y_perm <- sample(y_std)
        perm_coefs[i] <- sum(x_std * y_perm) / (n - 1)
    }

    # One-sided p-value (testing if observed is greater than expected)
    count <- sum(perm_coefs >= observed_coef)
    p_value <- (1 + count) / (n_perms + 1)

    return(list(
        observed_coef = observed_coef,
        r_squared = observed_r_sq,
        p_value = p_value,
        perm_distribution = perm_coefs,
        n_obs = n
    ))
}

#' Format p-value for display
format_pvalue <- function(p) {
    if (is.na(p)) {
        return("NA")
    }
    if (p <= 0.0001) {
        return("p <= 0.0001")
    }
    return(sprintf("p = %.4f", p))
}

#' Print test results
print_test_result <- function(result, year) {
    cat(sprintf("%s:\n", year))
    cat(sprintf("  Observed coefficient: %.4f\n", result$observed_coef))
    cat(sprintf("  R-squared: %.4f\n", result$r_squared))
    cat(sprintf("  %s\n", format_pvalue(result$p_value)))
    sig <- if (result$p_value < 0.05) "SIGNIFICANT" else "NOT significant"
    cat(sprintf("  Interpretation: %s\n", sig))
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

cat("Loading data...\n\n")

# Read data
data_2019 <- read_delim("2019fix.csv", delim = ";", show_col_types = FALSE)
data_2025 <- read_delim("2025fix.csv", delim = ";", show_col_types = FALSE)

# Build networks
cat("Building networks...\n")
adj_2019 <- create_adjacency_matrix(data_2019)
adj_2025 <- create_adjacency_matrix(data_2025)

g_2019 <- create_graph(adj_2019)
g_2025 <- create_graph(adj_2025)

# Calculate centralities
cat("Calculating centralities...\n\n")
cent_2019 <- calculate_centralities(g_2019, data_2019)
cent_2025 <- calculate_centralities(g_2025, data_2025)

# Print header
cat("=======================================================\n")
cat("  PERMUTATION-BASED INFERENCE FOR CENTRALITY METRICS\n")
cat("=======================================================\n\n")

cat("Methodology:\n")
cat("- All predictors are z-standardized\n")
cat("- Betweenness is normalized [0,1] then z-standardized\n")
cat("- P-values: p = (1 + count) / (B + 1), one-sided\n")
cat("- B = 10,000 permutations\n\n")

# -----------------------------------------------------------------------------
# TEST 1: Degree -> Betweenness
# -----------------------------------------------------------------------------
cat("TEST 1: Degree -> Betweenness\n")
cat("-----------------------------------------------------------\n")

result_deg_2019 <- permutation_test(cent_2019, "degree", "betweenness")
result_deg_2025 <- permutation_test(cent_2025, "degree", "betweenness")

print_test_result(result_deg_2019, "2019")
cat("\n")
print_test_result(result_deg_2025, "2025")
cat("\n")

# -----------------------------------------------------------------------------
# TEST 2: Closeness -> Betweenness
# -----------------------------------------------------------------------------
cat("TEST 2: Closeness -> Betweenness\n")
cat("-----------------------------------------------------------\n")

result_close_2019 <- permutation_test(cent_2019, "closeness", "betweenness")
result_close_2025 <- permutation_test(cent_2025, "closeness", "betweenness")

print_test_result(result_close_2019, "2019")
cat("\n")
print_test_result(result_close_2025, "2025")
cat("\n")

# -----------------------------------------------------------------------------
# TEST 3: N.Lines -> Betweenness
# -----------------------------------------------------------------------------
cat("TEST 3: N.Lines -> Betweenness\n")
cat("-----------------------------------------------------------\n")

result_lines_2019 <- permutation_test(cent_2019, "n_lines", "betweenness")
result_lines_2025 <- permutation_test(cent_2025, "n_lines", "betweenness")

print_test_result(result_lines_2019, "2019")
cat("\n")
print_test_result(result_lines_2025, "2025")
cat("\n")

# =============================================================================
# VISUALIZATION
# =============================================================================

cat("Generating visualizations...\n")

# Helper function to create a permutation plot for both years
create_perm_plot <- function(result_2019, result_2025, title) {
    # Combine data for both years
    plot_data <- data.frame(
        coef = c(result_2019$perm_distribution, result_2025$perm_distribution),
        year = factor(rep(
            c("2019", "2025"),
            c(
                length(result_2019$perm_distribution),
                length(result_2025$perm_distribution)
            )
        ))
    )

    # Observed values for vertical lines
    observed_data <- data.frame(
        year = factor(c("2019", "2025")),
        observed = c(result_2019$observed_coef, result_2025$observed_coef)
    )

    ggplot(plot_data, aes(x = coef)) +
        geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
        geom_vline(
            data = observed_data, aes(xintercept = observed),
            color = "red", linewidth = 1.2, linetype = "solid"
        ) +
        facet_wrap(~year, scales = "free_x") +
        labs(
            title = title,
            subtitle = sprintf(
                "2019: coef=%.3f, %s | 2025: coef=%.3f, %s",
                result_2019$observed_coef, format_pvalue(result_2019$p_value),
                result_2025$observed_coef, format_pvalue(result_2025$p_value)
            ),
            x = "Permutation Coefficient",
            y = "Frequency"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10),
            strip.text = element_text(face = "bold", size = 12)
        )
}

# Create and save plots with error handling
tryCatch(
    {
        p1 <- create_perm_plot(
            result_deg_2019, result_deg_2025,
            "Permutation Test: Degree -> Betweenness"
        )
        ggsave("permutation_degree_betweenness.png", p1, width = 10, height = 5, dpi = 300)

        p2 <- create_perm_plot(
            result_close_2019, result_close_2025,
            "Permutation Test: Closeness -> Betweenness"
        )
        ggsave("permutation_closeness_betweenness.png", p2, width = 10, height = 5, dpi = 300)

        p3 <- create_perm_plot(
            result_lines_2019, result_lines_2025,
            "Permutation Test: N.Lines -> Betweenness"
        )
        ggsave("permutation_lines_betweenness.png", p3, width = 10, height = 5, dpi = 300)

        cat("Plots saved.\n\n")
    },
    error = function(e) {
        cat("\nWarning: Could not generate plots due to error:\n")
        cat(paste0("  ", e$message, "\n"))
        cat("  Try running: remove.packages('gtable'); install.packages('gtable')\n")
        cat("  Then restart R and run the script again.\n\n")
    }
)

# =============================================================================
# SUMMARY TABLE
# =============================================================================

summary_table <- data.frame(
    Predictor = rep(c("Degree", "Closeness", "N.Lines"), each = 2),
    Year = rep(c(2019, 2025), 3),
    Coefficient = c(
        result_deg_2019$observed_coef, result_deg_2025$observed_coef,
        result_close_2019$observed_coef, result_close_2025$observed_coef,
        result_lines_2019$observed_coef, result_lines_2025$observed_coef
    ),
    R_squared = c(
        result_deg_2019$r_squared, result_deg_2025$r_squared,
        result_close_2019$r_squared, result_close_2025$r_squared,
        result_lines_2019$r_squared, result_lines_2025$r_squared
    ),
    P_value = c(
        result_deg_2019$p_value, result_deg_2025$p_value,
        result_close_2019$p_value, result_close_2025$p_value,
        result_lines_2019$p_value, result_lines_2025$p_value
    ),
    stringsAsFactors = FALSE
)

summary_table$Significant <- ifelse(summary_table$P_value < 0.05, "YES", "NO")
summary_table$P_value <- sapply(summary_table$P_value, function(p) sprintf("%.4e", p))

write.csv(summary_table, "permutation_test_summary.csv", row.names = FALSE)

cat("=======================================================\n")
cat("SUMMARY TABLE\n")
cat("=======================================================\n\n")
print(summary_table)

cat("\n\nAnalysis complete. Results saved.\n")
