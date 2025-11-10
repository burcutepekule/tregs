div0 = function(x, y) {
  ifelse(y == 0, 0, x / y)
}

vectors_to_treg_df = function() {
  data.frame(
    x = treg_x,
    y = treg_y,
    active_age = treg_active_age,
    phenotype = treg_phenotype,
    activity_SAMPs_binary = treg_activity_SAMPs_binary,
    id = treg_id
  )
}

matrix_to_pathogen_df = function() {
  if(nrow(pathogen_coords) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), id = numeric(0)))
  }
  data.frame(
    x = pathogen_coords[, "x"],
    y = pathogen_coords[, "y"], 
    id = pathogen_coords[, "id"]
  )
}

matrix_to_commensal_df = function() {
  if(nrow(commensal_coords) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), id = numeric(0)))
  }
  data.frame(
    x = commensal_coords[, "x"],
    y = commensal_coords[, "y"],
    id = commensal_coords[, "id"]
  )
}

agent_colors = c(
  epithelial_healthy = "#B0E2FF",
  epithelial_inj_1   = "#8CB4E5",
  epithelial_inj_2   = "#6987CC",
  epithelial_inj_3   = "#465AB2",
  epithelial_inj_4   = "#232D99",
  epithelial_inj_5   = "#000080",
  phagocyte_M0       = "grey70",
  phagocyte_M1_L_0   = "#F8C8E8",
  phagocyte_M1_L_1   = "#F397D6",
  phagocyte_M1_L_2   = "#E754C4",
  phagocyte_M1_L_3   = "#D12CA0",
  phagocyte_M1_L_4   = "#A5177A",
  phagocyte_M1_L_5   = "#6B0C4F",
  phagocyte_M2_L_0   = "#CDEFE3", 
  phagocyte_M2_L_1   = "#97D6BC",  
  phagocyte_M2_L_2   = "#61BD96",  
  phagocyte_M2_L_3   = "#3BA578",  
  phagocyte_M2_L_4   = "#2E8B57",  
  phagocyte_M2_L_5   = "#1F5C3B",  
  phagocyte_M0   = "grey70",
  phagocyte_M1   = "#E754C4",
  phagocyte_M2   = "#3BA578",
  treg_resting = "#D8BFD8",
  treg_active  = "#967BB6",
  commensal    = "turquoise2",
  pathogen     = "firebrick1",
  C_ROS = "black",
  C_M0  = "grey70",
  C_M1  = "#D12CA0",
  C_M2  = "#3BA578",
  P_ROS = "black",
  P_M0  = "grey70",
  P_M1  = "#D12CA0",
  P_M2  = "#3BA578"
)


plot_faceted = function(data, variables, title) {
  data_long = data %>%
    dplyr::select(t, sterile, allow_tregs, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
    geom_line(alpha = 1, linewidth = 1) +
    facet_grid(sterile ~ allow_tregs, labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = title, x = "Time", y = "Count", color = "Agent")
  
  return(p)
}

plot_faceted_stationary = function(data, variables, title) {
  data_long = data %>%
    dplyr::select(t, sterile, allow_tregs, first_stable_t, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  # Summarize the vertical line position per facet
  vline_data = data_long %>%
    distinct(sterile, allow_tregs, first_stable_t)
  
  p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
    geom_line(alpha = 1, linewidth = 1) +
    geom_vline(data = vline_data, 
               aes(xintercept = first_stable_t), 
               linetype = "dashed", color = "black") +
    facet_grid(sterile ~ allow_tregs, labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = title, x = "Time", y = "Count", color = "Agent")
  
  return(p)
}


plot_faceted_8 = function(data, variables, title) {
  data_long = data %>%
    dplyr::select(t, sterile, allow_tregs, randomize_tregs, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p = ggplot(data_long, aes(x = t, y = value, color = variable)) +
    geom_line(alpha = 1, linewidth = 1) +
    # facet_grid(sterile ~ allow_tregs + randomize_tregs, labeller = label_both) +
    facet_grid(randomize_tregs ~ allow_tregs + sterile, labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = title, x = "Time", y = "Count", color = "Agent")
  
  return(p)
}

plot_grid_DAMPs = function() {
  
  # Create full grid background (invisible or white)
  full_grid = expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  # --- Convert DAMP matrix to long format ---
  damps_df = as.data.frame(DAMPs)
  colnames(damps_df) = paste0("x", 1:ncol(damps_df))
  damps_df = damps_df %>%
    mutate(y = nrow(damps_df):1) %>%
    pivot_longer(cols = starts_with("x"), names_to = "x", values_to = "value") %>%
    mutate(x = as.integer(gsub("x", "", x)))
  
  # --- Epithelial layer: blue/orange by health ---
  epithelial_df = epithelium %>%
    mutate(
      type = ifelse(level_injury == 0, "epithelial_healthy",
                    paste0("epithelial_inj_", level_injury)),
      fill = agent_colors[type],
      y = grid_size+ 1 - y
    ) %>%
    dplyr::select(x, y, fill)
  
  # --- Plot ---
  p_damps = ggplot() +
    geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
    # DAMP heatmap with colorbar
    geom_tile(data = damps_df, aes(x = x, y = y, fill = value)) +
    scale_fill_gradient(low = "white", high = "black", name = "DAMPs Density", limits = c(0, lim_DAMP)) +
    
    # Epithelium on top with inline fill colors (no legend)
    geom_tile(data = epithelial_df, aes(x = x, y = y), fill = epithelial_df$fill, width = 1, height = 1) +
    
    coord_fixed() +
    # labs(title = "DAMP Heatmap with Epithelium (No Legend for Epithelium)") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )
  return(p_damps)
  
}

plot_grid_SAMPs = function() {
  
  # Create full grid background (invisible or white)
  full_grid = expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  # --- Convert DAMP matrix to long format ---
  samps_df = as.data.frame(SAMPs)
  colnames(samps_df) = paste0("x", 1:ncol(samps_df))
  samps_df = samps_df %>%
    mutate(y = nrow(samps_df):1) %>%
    pivot_longer(cols = starts_with("x"), names_to = "x", values_to = "value") %>%
    mutate(x = as.integer(gsub("x", "", x)))
  
  # --- Epithelial layer: blue/orange by health ---
  epithelial_df = epithelium %>%
    mutate(
      type = ifelse(level_injury == 0, "epithelial_healthy",
                    paste0("epithelial_inj_", level_injury)),
      fill = agent_colors[type],
      y = grid_size+ 1 - y
    ) %>%
    dplyr::select(x, y, fill)
  
  # --- Plot ---
  p_samps = ggplot() +
    geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
    # DAMP heatmap with colorbar
    geom_tile(data = samps_df, aes(x = x, y = y, fill = value)) +
    scale_fill_gradient(low = "white", high = "black", name = "SAMPs Density", limits = c(0, lim_DAMP)) +
    
    # Epithelium on top with inline fill colors (no legend)
    geom_tile(data = epithelial_df, aes(x = x, y = y), fill = epithelial_df$fill, width = 1, height = 1) +
    
    coord_fixed() +
    # labs(title = "DAMP Heatmap with Epithelium (No Legend for Epithelium)") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )
  return(p_samps)
  
}

plot_grid_ROS = function() {
  
  # Create full grid background (invisible or white)
  full_grid = expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  # --- Convert DAMP matrix to long format ---
  ros_df = as.data.frame(ROS)
  colnames(ros_df) = paste0("x", 1:ncol(ros_df))
  ros_df = ros_df %>%
    mutate(y = nrow(ros_df):1) %>%
    pivot_longer(cols = starts_with("x"), names_to = "x", values_to = "value") %>%
    mutate(x = as.integer(gsub("x", "", x)))
  
  # --- Epithelial layer: blue/orange by health ---
  epithelial_df = epithelium %>%
    mutate(
      type = ifelse(level_injury == 0, "epithelial_healthy",
                    paste0("epithelial_inj_", level_injury)),
      fill = agent_colors[type],
      y = grid_size+ 1 - y
    ) %>%
    dplyr::select(x, y, fill)
  
  # --- Plot ---
  p_ros = ggplot() +
    geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
    # DAMP heatmap with colorbar
    geom_tile(data = ros_df, aes(x = x, y = y, fill = value)) +
    scale_fill_gradient(low = "white", high = "black", name = "ROS Density", limits = c(0, lim_ROS)) +
    
    # Epithelium on top with inline fill colors (no legend)
    geom_tile(data = epithelial_df, aes(x = x, y = y), fill = epithelial_df$fill, width = 1, height = 1) +
    
    coord_fixed() +
    # labs(title = "ROS Heatmap with Epithelium (No Legend for Epithelium)") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )
  return(p_ros)
}

plot_grid_antiinf = function() {
  
  # Create a base dataframe for the grid
  grid = expand.grid(x = 1:grid_size, y = 0:grid_size)
  
  # 1. Epithelial cell state (y = 0 only)
  # 1. Epithelial cell state (y = 0 only)
  epithelial_layer = epithelium %>%
    mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                         paste0("epithelial_inj_", level_injury)))
  
  # 1. Phagocyte state 
  phagocytes_plot = phagocytes %>%
    mutate(type = ifelse(phenotype==0, "phagocyte_M0", 
                         ifelse(phenotype==1, "phagocyte_M1","phagocyte_M2")))
  
  phagocytes_plot = phagocytes_plot %>% filter(type %in% c("phagocyte_M2"))
  
  tregs_plot = tregs %>%
    mutate(type = ifelse(phenotype==0, "treg_resting","treg_active"))
  
  tregs_plot = tregs_plot %>% filter(type=="treg_active")
  
  # Create full grid background (invisible or white)
  full_grid = expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  # 2. Combine all agents into one dataframe with their type
  all_types = c(
    "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
    "phagocyte_M2","treg_active")
  
  agent_plot_df = bind_rows(
    epithelial_layer %>% dplyr::select(x, y, type),
    phagocytes_plot %>% dplyr::select(x, y, type),
    tregs_plot  %>% dplyr::select(x, y, type)
  ) %>%
    mutate(type = factor(type, levels = all_types))
  
  p_lym = ggplot() +
    geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
    geom_tile(data = agent_plot_df, aes(x = x, y = y, fill = type), width = 1, height = 1) +
    scale_fill_manual(
      values = agent_colors,
      name = "Cell Type",
      drop = FALSE  # =- ensures unused levels are shown
    ) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size+1)) +
    scale_y_reverse(expand = c(0, 0)) +  # Flip Y so y=0 is on top
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      legend.position = "right"
    )
  
  return(p_lym)
}

plot_grid_phagocyte_M1 = function() {
  
  # Create a base dataframe for the grid
  grid = expand.grid(x = 1:grid_size, y = 0:grid_size)
  
  # 1. Epithelial cell state (y = 0 only)
  epithelial_layer = epithelium %>%
    mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                         paste0("epithelial_inj_", level_injury)))
  # 1. Phagocyte state
  phagocytes_plot = phagocytes %>%
    mutate(type = ifelse(phenotype==0, "phagocyte_M0", 
                         ifelse(phenotype==1, "phagocyte_M1","phagocyte_M2")))
  
  phagocytes_plot = phagocytes_plot %>% filter(type=="phagocyte_M1")
  # Create full grid background (invisible or white)
  full_grid = expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  phagocytes_plot$level = sapply(phagocytes_plot$bacteria_registry, sum)
  
  if(dim(phagocytes_plot)[1]>0){
    phagocytes_plot$type = paste0("phagocyte_M1_L_", phagocytes_plot$level)
  }
  
  # Combine all agents into one dataframe with their type
  agent_plot_df = bind_rows(
    epithelial_layer %>% dplyr::select(x, y, type),
    phagocytes_plot %>% dplyr::select(x, y, type)
  )
  
  # 2. Combine all agents into one dataframe with their type
  # all_types = c(
  #   "epithelial_healthy", "epithelial_unhealthy","phagocyte_M1")
  
  all_types = c(
    "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
    "phagocyte_M1_L_0","phagocyte_M1_L_1","phagocyte_M1_L_2","phagocyte_M1_L_3","phagocyte_M1_L_4","phagocyte_M1_L_5")
  
  agent_plot_df = bind_rows(
    epithelial_layer %>% dplyr::select(x, y, type),
    phagocytes_plot %>% dplyr::select(x, y, type)
  ) %>%
    mutate(type = factor(type, levels = all_types))
  
  p_lym = ggplot() +
    geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
    geom_tile(data = agent_plot_df, aes(x = x, y = y, fill = type), width = 1, height = 1) +
    scale_fill_manual(
      values = agent_colors,
      name = "Cell Type",
      drop = FALSE  # =- ensures unused levels are shown
    ) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size+1)) +
    scale_y_reverse(expand = c(0, 0)) +  # Flip Y so y=0 is on top
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      legend.position = "right"
    )
  
  return(p_lym)
}

plot_grid_mat = function(mat, grid_size, lim_mat){
  # Create full grid background (invisible or white)
  full_grid = expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  # --- Convert DAMP matrix to long forDAMPs ---
  mat_df = as.data.frame(mat)
  colnames(mat_df) = paste0("x", 1:ncol(mat_df))
  mat_df = mat_df %>%
    mutate(y = nrow(mat_df):1) %>%
    pivot_longer(cols = starts_with("x"), names_to = "x", values_to = "value") %>%
    mutate(x = as.integer(gsub("x", "", x)))
  
  # --- Plot ---
  p = ggplot() +
    geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
    geom_tile(data = mat_df, aes(x = x, y = y, fill = value)) +
    scale_fill_gradient(low = "white", high = "black", name = "DAMP Density", limits = c(0, lim_mat)) +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )
  return(p)
}

plot_microbes_nomem = function() {
  # Count microbes
  num_pat = nrow(pathogens_lp)
  num_com = nrow(commensals_lp)
  
  # Create a new row of data
  new_data = data.frame(
    time = t,
    count = c(num_pat, num_com),
    type = c("pathogen", "commensal")
  )
  
  microbe_colors = agent_colors[c("pathogen", "commensal")]
  
  
  # Horizontal bar plot
  p = ggplot(new_data, aes(x = type, y = count, fill = type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = microbe_colors, name = "Microbe Type") +
    coord_flip() +  # Make bars horizontal
    ylim(0, max(200, num_pat, num_com)) +
    # labs(x = "", y = "Microbe Count", title = paste("Time =", t)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  
  return(p)
}

plot_phagocyte_nomem = function() {
  # Count lymph
  # 1. Phagocyte state 
  phagocytes_plot = phagocytes %>%
    mutate(type = ifelse(phenotype==0, "phagocyte_M0", 
                         ifelse(phenotype==1, "phagocyte_M1","phagocyte_M2")))
  
  num_M0 = dim(phagocytes_plot %>% filter(type=="phagocyte_M0"))[1]
  num_M1 = dim(phagocytes_plot %>% filter(type=="phagocyte_M1"))[1]
  num_M2 = dim(phagocytes_plot %>% filter(type=="phagocyte_M2"))[1]
  
  # Create a new row of data
  new_data = data.frame(
    time = t,
    count = c(num_M0, num_M1, num_M2),
    type = c("phagocyte_M0", "phagocyte_M1", "phagocyte_M2")
  )
  
  cell_colors = agent_colors[c("phagocyte_M0", "phagocyte_M1", "phagocyte_M2")]
  
  
  # Horizontal bar plot
  p = ggplot(new_data, aes(x = type, y = count, fill = type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cell_colors, name = "Cell Type") +
    coord_flip() +  # Make bars horizontal
    ylim(0, max(30, num_M0, num_M1, num_M2)) +
    theme_minimal() +
    theme(legend.position = "none")
  return(p)
}

plot_treg_nomem = function() {
  # Count lymph
  # Treg state 
  tregs_plot = tregs %>%
    mutate(type = ifelse(phenotype==0, "treg_resting","treg_active"))
  
  num_treg_0 = dim(tregs_plot %>% filter(type=="treg_resting"))[1]
  num_treg_1 = dim(tregs_plot %>% filter(type=="treg_active"))[1]
  
  # Create a new row of data
  new_data = data.frame(
    time = t,
    count = c(num_treg_0, num_treg_1),
    type = c("treg_resting", "treg_active")
  )
  
  cell_colors = agent_colors[c("treg_resting", "treg_active")]
  
  
  # Horizontal bar plot
  p = ggplot(new_data, aes(x = type, y = count, fill = type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cell_colors, name = "Cell Type") +
    coord_flip() +  # Make bars horizontal
    ylim(0, max(30, num_treg_0, num_treg_1)) +
    theme_minimal() +
    theme(legend.position = "none")
  return(p)
}

plot_cumdeath_nomem = function() {
  # Count microbes
  
  # Create a new row of data
  new_data = data.frame(
    time = t,
    count = c(pathogens_killed_by_ROS, pathogens_killed_by_Mac,commensals_killed_by_ROS, commensals_killed_by_Mac),
    type = c("P_ROS","P_M0","P_M1","P_M2","C_ROS","C_M0","C_M1","C_M2")
  )
  
  microbe_colors = agent_colors[c("P_ROS","P_M0","P_M1","P_M2","C_ROS","C_M0","C_M1","C_M2")]
  
  
  # Horizontal bar plot
  p = ggplot(new_data, aes(x = type, y = count, fill = type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = microbe_colors, name = "Microbe Type") +
    coord_flip() +  # Make bars horizontal
    ylim(0, max(200, pathogens_killed_by_ROS, pathogens_killed_by_Mac,commensals_killed_by_ROS, commensals_killed_by_Mac)) +
    labs(x = "", y = "Cum. death") +
    theme_minimal() +
    theme(legend.position = "none")
  
  
  return(p)
}

plot_grid_pathogens = function() {
  
  # Create a base dataframe for the grid
  grid = expand.grid(x = 1:grid_size, y = 0:grid_size)
  
  # 1. Epithelial cell state (y = 0 only)
  # 1. Epithelial cell state (y = 0 only)
  epithelial_layer = epithelium %>%
    mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                         paste0("epithelial_inj_", level_injury)))
  
  # Create full grid background (invisible or white)
  full_grid = expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  # 2. Combine all agents into one dataframe with their type
  all_types = c(
    "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
    "pathogen"
  )
  
  pathogen_density = pathogens_lp %>%
    dplyr::group_by(x, y) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    mutate(
      type = "pathogen"
    )
  
  # Check if there is variation in counts
  if (nrow(pathogen_density) > 0) {
    if (length(unique(pathogen_density$count)) == 1) {
      # All counts are the same: use fixed alpha
      pathogen_density$alpha_val = 0.8
    } else {
      # Vary alpha by density
      pathogen_density$alpha_val = pmin(1, pathogen_density$count / max(pathogen_density$count))
    }
  }
  
  epithelial_layer_plot = epithelial_layer %>%
    dplyr::select(x, y, type) %>%
    mutate(alpha_val = 1)
  
  agent_plot_df = bind_rows(
    epithelial_layer_plot,
    pathogen_density
  ) %>%
    mutate(type = factor(type, levels = all_types))
  
  
  p_mic = ggplot() +
    geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
    geom_tile(data = agent_plot_df,
              aes(x = x, y = y, fill = type, alpha = alpha_val),
              width = 1, height = 1) +
    scale_fill_manual(
      values = agent_colors,
      name = "Cell Type",
      drop = FALSE
    ) +
    scale_alpha_continuous(range = c(0.5, 1), guide = FALSE) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size + 1)) +
    scale_y_reverse(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      legend.position = "right"
    )
  
  return(p_mic)
}

plot_grid_commensals = function() {
  
  # Create a base dataframe for the grid
  grid = expand.grid(x = 1:grid_size, y = 0:grid_size)
  
  # 1. Epithelial cell state (y = 0 only)
  epithelial_layer = epithelium %>%
    mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                         paste0("epithelial_inj_", level_injury)))
  
  # Create full grid background (invisible or white)
  full_grid = expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  # 2. Combine all agents into one dataframe with their type
  all_types = c(
    "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
    "commensal"
  )
  
  commensal_density = commensals_lp %>%
    dplyr::group_by(x, y) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    mutate(
      type = "commensal"
    )
  
  # Check if there is variation in counts
  if (nrow(commensal_density) > 0) {
    if (length(unique(commensal_density$count)) == 1) {
      # All counts are the same: use fixed alpha
      commensal_density$alpha_val = 0.8
    } else {
      # Vary alpha by density
      commensal_density$alpha_val = pmin(1, commensal_density$count / max(commensal_density$count))
    }
  }
  
  epithelial_layer_plot = epithelial_layer %>%
    dplyr::select(x, y, type) %>%
    mutate(alpha_val = 1)
  
  agent_plot_df = bind_rows(
    epithelial_layer_plot,
    commensal_density
  ) %>%
    mutate(type = factor(type, levels = all_types))
  
  
  p_mic = ggplot() +
    geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
    geom_tile(data = agent_plot_df,
              aes(x = x, y = y, fill = type, alpha = alpha_val),
              width = 1, height = 1) +
    scale_fill_manual(
      values = agent_colors,
      name = "Cell Type",
      drop = FALSE
    ) +
    scale_alpha_continuous(range = c(0.5, 1), guide = FALSE) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size + 1)) +
    scale_y_reverse(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      legend.position = "right"
    )
  
  return(p_mic)
}

plot_grid_resting = function() {
  
  # Create a base dataframe for the grid
  grid = expand.grid(x = 1:grid_size, y = 0:grid_size)
  
  # 1. Epithelial cell state (y = 0 only)
  # 1. Epithelial cell state (y = 0 only)
  epithelial_layer = epithelium %>%
    mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                         paste0("epithelial_inj_", level_injury)))
  # 1. Phagocyte state 
  phagocytes_plot = phagocytes %>%
    mutate(type = ifelse(phenotype==0, "phagocyte_M0", 
                         ifelse(phenotype==1, "phagocyte_M1","phagocyte_M2")))
  
  phagocytes_plot = phagocytes_plot %>% filter(type=="phagocyte_M0")
  
  tregs_plot = tregs %>%
    mutate(type = ifelse(phenotype==0, "treg_resting","treg_active"))
  tregs_plot = tregs_plot %>% filter(type=="treg_resting")
  
  # Create full grid background (invisible or white)
  full_grid = expand.grid(x = 1:grid_size, y = 1:grid_size)
  
  # 2. Combine all agents into one dataframe with their type
  all_types = c(
    "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
    "phagocyte_M0","treg_resting"
  )
  
  agent_plot_df = bind_rows(
    epithelial_layer %>% dplyr::select(x, y, type),
    phagocytes_plot %>% dplyr::select(x, y, type),
    tregs_plot %>% dplyr::select(x, y, type)
  ) %>%
    mutate(type = factor(type, levels = all_types))
  
  
  
  p_lym = ggplot() +
    geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
    geom_tile(data = agent_plot_df, aes(x = x, y = y, fill = type), width = 1, height = 1) +
    scale_fill_manual(
      values = agent_colors,
      name = "Cell Type",
      drop = FALSE  # =- ensures unused levels are shown
    ) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size+1)) +
    scale_y_reverse(expand = c(0, 0)) +  # Flip Y so y=0 is on top
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      legend.position = "right"
    )
  
  return(p_lym)
}

plot_simtime_simple = function(){
  
  p_DAMPs  = plot_grid_DAMPs()
  p_SAMPs  = plot_grid_SAMPs()
  p_ROS    = plot_grid_ROS()
  p_d      = plot_grid_antiinf()
  p_a      = plot_grid_phagocyte_M1()
  p_lymp   = plot_grid_resting()
  p_com    = plot_grid_commensals()
  p_pat    = plot_grid_pathogens()
  
  p_mic_bar  = plot_microbes_nomem()
  p_lym_bar  = plot_phagocyte_nomem()
  p_treg_bar = plot_treg_nomem()
  p_cumdeath = plot_cumdeath_nomem()
  
  if(sterile==0){
    injury_type= 'pathogenic'
  }else{
    injury_type= 'sterile'
  }
  
  tit_add = paste0('Sim. time : ',t,', Injury: ',injury_type,', Tregs allowed : ',allow_tregs)
  row_0 = plot_grid(p_mic_bar,p_lym_bar,p_treg_bar, p_cumdeath, ncol = 4, rel_widths =  c(1,1,1,1))
  row_1 = plot_grid(p_com,p_pat,p_a,p_d,ncol = 4, rel_widths =  c(1,1,1.1,1.05))
  row_2 = plot_grid(p_DAMPs,p_SAMPs,p_ROS,p_lymp,ncol = 4, rel_widths =  c(1,1,0.98,1.1))
  
  combined_plot = plot_grid(row_0,row_1,row_2,align='v',nrow = 3, rel_heights =  c(0.35,1,1))
  # combined_plot = plot_grid(row_1,row_2,align='v',nrow = 2, rel_heights =  c(1,1))
  
  # Add a title
  combined_plot = plot_grid(
    ggdraw() + draw_label(tit_add, size = 14, hjust = 0.5),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.05, 1)  # Adjust title vs. plot height ratio
  )
  
  # print(combined_plot)
  return(combined_plot)
}
