# Load libraries
library(igraph)   # network analysis
library(dplyr)    # data manipulation
library(tidyr)    # data manipulation
library(readr)    # reading/writing delimited files
library(stringr)  # text handling
library(scales)   # rescaling for plot aesthetics

# Expected columns in node_attributes_data.csv:
# patient_id, patient_hrr, npi, provider_hrr, provider_specialty,
# cancer_site, pc_flag, pc_specialist_flag, encounter_year

### Part 1: Load data
df_attr <- read.csv("node_attributes_data.csv", stringsAsFactors = FALSE)

### Part 2: Build provider network
# Keep unique patient–provider pairs (unweighted bipartite incidence)
benephys_pairs <- df_attr %>%
  select(patient_id, npi) %>%
  distinct()

# Build bipartite graph: patients + providers
g_bip <- graph_from_data_frame(benephys_pairs, directed = FALSE)
# Assign bipartite types
V(g_bip)$type <- bipartite_mapping(g_bip)$type

# Project to provider–provider (unipartite) with multiplicity = shared patients
proj <- bipartite_projection(g_bip, multiplicity = TRUE)

# choose the projection that matches NPIs
g_prov <- if (sum(V(proj$proj1)$name %in% df_attr$npi) >=
              sum(V(proj$proj2)$name %in% df_attr$npi)) proj$proj1 else proj$proj2

cat("Providers in projected network: ", vcount(g_prov), "\n")

# Add physician specialty as node attribute to annual network
# Build an attribute frame per NPI (one row per provider)
provider_attributes <- df_attr %>%
  select(npi, provider_hrr, provider_specialty, pc_specialist_flag) %>%
  mutate(
    provider_specialty = ifelse(is.na(provider_specialty) | provider_specialty=="",
                                "Unknown", provider_specialty)
  ) %>%
  distinct(npi, .keep_all = TRUE)

# Compute patient volume (unique patients per provider)
patient_volume <- df_attr %>%
  group_by(npi) %>%
  summarise(patient_volume = n_distinct(patient_id), .groups = "drop")

# Merge patient volume into provider_attributes
provider_attributes <- provider_attributes %>%
  left_join(patient_volume, by = "npi")

# Restrict attribute table to nodes present in the graph
npi_nodes <- tibble(npi = V(g_prov)$name)
provider_attributes <- npi_nodes %>%
  left_join(provider_attributes, by = "npi")

# Add attributes to graph
g_prov <- set_vertex_attr(g_prov, name = "spec", value = provider_attributes$provider_specialty)

# Edge weights already exist from projection; compute node strength (sum of incident weights)
V(g_prov)$strength <- strength(g_prov, mode = "all", weights = E(g_prov)$weight)

# Patient volume as vertex attribute
V(g_prov)$patient_volume <- provider_attributes$patient_volume

# Save objects
saveRDS(g_prov, file = "provider_network.rds")
saveRDS(provider_attributes, file = "node_attributes.rds")

### Part 3: Prepare for HRR visualization
# Define groupings to collapse granular specialties
pcp_spec <- c("Family Medicine","Family Practice","General Practice","Internal Medicine",
              "Geriatric Medicine","Certified Clinical Nurse Specialist",
              "Nurse Practitioner","Physician Assistant")

surgeon_spec <- c("General Surgery","Orthopedic Surgery","Cardiac Surgery","Thoracic Surgery",
                  "Vascular Surgery","Neurosurgery","Urology","Ophthalmology",
                  "Surgical Oncology","Colorectal Surgery (Proctology)")

medonc_spec <- c("Hematology-Oncology","Medical Oncology","Hematology")
pal_taxonomy <- "Hospice and Palliative Care"

term_spec <- function(spec) {
  case_when(
    spec %in% pcp_spec                     ~ "PCP",
    spec %in% medonc_spec                  ~ "Medical Oncologist",
    spec %in% surgeon_spec                 ~ "Surgeon",
    spec == "Radiation Oncology"           ~ "Radiation Oncologist",
    spec == "Hospitalist"                  ~ "Hospitalist",
    spec == pal_taxonomy                   ~ "Formally-Trained PC Specialist",
    TRUE                                   ~ "Others"
  )
}

# Color palette
specialty_color <- c(
  "Formally-Trained PC Specialist" = "red",
  "non-specialist PC"              = "#4CBB17",   # green for non-spec PC PCPs
  "Non-PC"                         = "lightgrey",
  "Medical Oncologist"             = "blue",
  "Radiation Oncologist"           = "orange",
  "Surgeon"                        = "purple",
  "PCP"                            = "#4CBB17",
  "Hospitalist"                    = "saddlebrown",
  "Others"                         = "lightgrey"
)

df_network_fig <- provider_attributes %>%
  mutate(
    specialty_group = term_spec(provider_specialty),
    
    # type of PC provision
    pc_type = case_when(
      pc_specialist_flag == 1 & provider_specialty == pal_taxonomy ~ "Formally-Trained PC Specialist",
      pc_specialist_flag == 1 & provider_specialty != pal_taxonomy ~ "non-specialist PC",
      TRUE                                                         ~ "Non-PC"
    ),
    
    # SHAPE indicates PC provision (square = provided PC; circle = non-PC)
    shape = ifelse(pc_type %in% c("Formally-Trained PC Specialist","non-specialist PC"), "square", "circle"),
    
    # COLOR indicates SPECIALTY
    color = specialty_color[specialty_group]
  )

### Part 4: Visualize the subnetwork
# Output dir
out_dir <- "network_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Plot once per HRR (unique)
for (hrr in unique(df_network_fig$provider_hrr)) {
  
  cat("Plotting subnetwork for HRR: ", hrr)
  
  # Providers to include: this HRR; drop non-PC with 'Others' specialty to declutter
  npi_list <- df_network_fig %>%
    filter(provider_hrr == hrr) %>%
    filter(!(pc_type == "Non-PC" & specialty_group == "Others")) %>%
    pull(npi)
  
  # Subgraph induced by those providers
  g_sub <- induced_subgraph(g_prov, vids = V(g_prov)[name %in% npi_list])
  df_sub <- df_network_fig %>% filter(npi %in% npi_list)
  
  # Attach visual attributes from df_network_fig
  idx <- match(V(g_sub)$name, df_sub$npi)
  V(g_sub)$vertex.shape       <- df_sub$shape[idx]
  V(g_sub)$vertex.color       <- df_sub$color[idx]
  V(g_sub)$vertex.frame.color <- ifelse(df_sub$pc_type[idx] != "Non-PC", "black", NA)
  V(g_sub)$vertex.size        <- rescale(df_sub$patient_volume[idx], to = c(3, 6), na.rm = TRUE)
  E(g_sub)$width              <- rescale(E(g_sub)$weight, to = c(0.5, 2.0), na.rm = TRUE)
  
  # Remove isolcates
  g_sub <- igraph::simplify(g_sub, remove.loops = TRUE)
  g_sub <- delete_vertices(g_sub, which(degree(g_sub) == 0))
  
  # Layout (Kamada–Kawai is a good default for smaller subnetworks)
  set.seed(139)
  layout_kk <- layout_with_kk(g_sub)
  
  # File name
  file_name <- file.path(out_dir, paste0("network_figure_", hrr, ".png"))
  
  # Plot
  png(file_name, width = 1600, height = 1200, res = 150)
  par(mar = c(1,1,2,1))
  plot(
    g_sub,
    layout = layout_kk,
    vertex.label = NA,
    vertex.shape = V(g_sub)$vertex.shape,
    vertex.color = V(g_sub)$vertex.color,
    vertex.frame.color = V(g_sub)$vertex.frame.color,
    vertex.size = V(g_sub)$vertex.size,
    edge.color = "gray66",
    edge.width = E(g_sub)$width,
    main = paste0("HRR: ", hrr, " — Palliative Care (PC) Subnetwork")
  )
  
  # Legend (compact, grouped)
  legend_items <- c(
    "Provided PC (shape)", "Yes", "No", "",
    "Provider Specialty (color)",
    "PC Specialist","Medical Oncologist","Radiation Oncologist",
    "Surgeon","PCP","Hospitalist","Others"
  )
  legend_shapes <- c(NA, 22, 21, NA, NA, rep(21, 7))
  legend_colors <- c(NA, "black", "black", NA, NA, rep(NA, 7))
  legend_fill   <- c(NA, NA, NA, NA, NA,
                     "red","blue","orange","purple","#4CBB17","saddlebrown","lightgrey")
  
  legend("topright",
         legend = legend_items,
         pch = legend_shapes,
         pt.bg = legend_fill,
         col = legend_colors,
         pt.cex = 1.0,
         cex = 0.8,
         bty = "n",
         text.font = c(2,1,1,1,2, rep(1,7)),
         x.intersp = 0.7,
         y.intersp = 1.1)
  
  dev.off()
}