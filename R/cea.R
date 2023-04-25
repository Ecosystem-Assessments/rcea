#' Cumulative effects assessments
#'
#' Series of functions to assess cumulative effects and related metrics using the Halpern et al. 2008 or the Beauchesne et al. 2021 methodologies. 
#'
#' @eval arguments(c("drivers","vc","sensitivity","metaweb","trophic_sensitivity"))
#' @param motif_summary list of parameters obtained from \link{ncea_motifs} to allow for faster assessments, since that part of the assessment is the longest to run. 
#' @param exportAs string, the type of object that should be created by \link{exposure} ("list" | "stars"), \link{cea} ("list" | "stars") and \link{ncea} ("data.frame" | "stars")
#'
#' @describeIn cea evaluate exposure of valued components to environmental drivers
#' @export
#'
#' @examples
#' # Data
#' drivers <- rcea:::drivers 
#' vc <- rcea:::vc
#' sensitivity <- rcea:::sensitivity
#' metaweb <- rcea:::metaweb
#' data(trophic_sensitivity)
#'
#' # Exposure
#' (expo <- exposure(drivers, vc, "stars"))
#' plot(expo)
#' expo <- merge(expo, name = "vc") |>
#'         split("drivers")
#' plot(expo)
#'
#' # Species-scale effects
#' (halpern <- cea(drivers, vc, sensitivity, "stars"))
#' plot(halpern)
#' halpern <- merge(halpern, name = "vc") |>
#'           split("drivers")
#' plot(halpern)
#'
#' \dontrun{
#' # Network-scale effects 
#' beauchesne <- ncea(drivers, vc, sensitivity, metaweb, trophic_sensitivity)
#' plot(beauchesne$net)
#' plot(beauchesne$direct)
#' plot(beauchesne$indirect)
#' }
exposure <- function(drivers, vc, exportAs = "list") {
  # Drivers
  dr_df <- as.data.frame(drivers) |>
        dplyr::select(-x,-y)
  
  # Valued components
  vc_df <- as.data.frame(vc) |>
        dplyr::select(-x,-y)
        
  # Exposure of valued components to drivers (Dj * VCi)
  dat <- apply(
    dr_df,
    MARGIN = 2,
    function(x) {
      sweep(vc_df, MARGIN=1, x, `*`)
    }
  ) 
  
  # Return
  if (exportAs == "list") {
    dat
  } else if (exportAs == "stars") {
    cea_stars(dat, drivers, vc)
  }
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea effects assessment using the Halpern approach, i.e. assessment of direct effects of drivers on valued components
#' @export
cea <- function(drivers, vc, sensitivity, exportAs = "list") {
  # Exposure
  dat <- exposure(drivers, vc)

  # Sensitivity
  nmDr <- names(dat)
  nmVC <- names(dat[[1]])
  sensitivity <- sensitivity[nmVC, nmDr] 
  
  # Effect of drivers on valued components (D * VC * u)
  for(i in 1:length(dat)) {
    dat[[i]] <- sweep(dat[[i]], MARGIN = 2, sensitivity[,i], `*`)
  }

  # Return
  if (exportAs == "list") {
    dat
  } else if (exportAs == "stars") {
    xy <- sf::st_coordinates(vc)
    drNames <- names(drivers)
    cea_stars(dat, drivers, vc)
  }

}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea assessment of effects across all motifs 
#' @export
ncea <- function(drivers, vc, sensitivity, metaweb, trophic_sensitivity, w_d = 0.5, w_i = 0.25, motif_summary = NULL) {
  # NOTE: `motif_summary` Allows a user to directly provide the motif_summary object if it's 
  #       been run before. This saves computation time if you wish to play around with weights 
  #       or other method parameters.
  if (is.null(motif_summary)) {
    # 3-species motifs for full metaweb
    motifs <- triads(metaweb)

    # Direct effects, i.e. Halpern approach
    direct_effect <- cea(drivers, vc, sensitivity) |>
                     cea_array()

    # Pathways of direct effect
    direct_pathways <- cea_pathways(direct_effect, vc)

    # Pathways of indirect effect
    indirect_pathways <- ncea_pathways_(direct_pathways, motifs)

    # Effects for motifs in each cell 
    motif_summary <- ncea_motifs(direct_effect, indirect_pathways)
  }
  
  # Measure effects on each motif
  motif_effects <- ncea_effects(motif_summary, w_d, w_i)
  
  # params for stars object creation
  xy <- sf::st_coordinates(drivers)
  drNames <- names(drivers)
  
  # Species contribution to indirect effects 
  species_contribution <- get_species_contribution(motif_effects) |>
                          dplyr::rename(vc_id = interaction) |>
                          cea_stars(drivers, vc)
  
  # Direct & indirect effects
  direct_indirect <- get_direct_indirect(motif_effects)
  direct <- dplyr::filter(direct_indirect, direct) |>
            cea_stars(drivers, vc)
  indirect <- dplyr::filter(direct_indirect, !direct) |>
              cea_stars(drivers, vc)

  # Net effects
  net <- get_net(motif_effects) |>
         cea_stars(drivers, vc)
  
  # Effects / km2 
  cekm <- get_cekm(motif_effects, vc)
  
  # Return
  list(
    motif_effects = motif_effects,
    net = net,
    direct = direct, 
    indirect = indirect,
    species_contribution = species_contribution,
    cekm = cekm
  )
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea create array from list of cea matrices
#' @export
cea_array <- function(dat) {
  unlist(dat) |>
  array(
    dim = c(
      nrow(dat[[1]]), 
      ncol(dat[[1]]), 
      length(dat)
    ),
    dimnames = list(
      c(), 
      names(dat[[1]]), 
      names(dat)
    )
  )  
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea create stars object from list of cea matrices or cea array
#' @export
cea_stars <- function(dat, drivers, vc) {
  xy <- sf::st_coordinates(drivers)
  drNames <- names(drivers)
  vc_index <- data.frame(
    vc = names(vc),
    vc_id = 1:length(vc)
  )
    
  # For data exported as list format (exposure and cea assessments)
  if ("list" %in% class(dat)) {
    for(i in 1:length(dat)) {
      dat[[i]] <- cbind(xy, dat[[i]]) |>
                  dplyr::mutate(drivers = drNames[i])
    }
    dat <- dplyr::bind_rows(dat) |>
           stars::st_as_stars(dims = c("x","y","drivers"))          
  }
  
  # For data exported in long format, from the ncea assessment
  if ("data.frame" %in% class(dat)) {
    # Get coordinates with repeated driver names
    xy$id_cell <- 1:nrow(xy)
    xyd <- cbind(
      dplyr::slice(xy, rep(1:dplyr::n(), each = length(drNames))),
      drivers = rep(drNames, nrow(xy))
    )
    
    # Pivot data to get drivers as lines and add vc names
    dat <- tidyr::pivot_longer(
      dat,
      cols = -c(id_cell, vc_id), 
      names_to = "drivers", 
      values_to = "effect"
    ) |>
    dplyr::left_join(vc_index, by = "vc_id") |>
    dplyr::select(-vc_id) |>
    tidyr::pivot_wider(names_from = vc, values_from = effect) 
    
    # Join with xy data and transform into stars object
    dat <- dplyr::left_join(xyd, dat, by = c("id_cell", "drivers")) |>
           dplyr::select(-id_cell) |>
           dplyr::select(x, y, drivers, dplyr::all_of(vc_index$vc)) |>
           stars::st_as_stars(dims = c("x","y","drivers"))          
  }
  
  # Return 
  dat
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea transform effects assessment into binary 2D matrix to assess the presence of an effect to a valued component in a specific grid cell
#' @export
cea_binary <- function(effect) {
  (effect / effect) |>
  apply(c(1,2), sum, na.rm = TRUE) |>
  as.data.frame() |>
  dplyr::mutate(id_cell = 1:dplyr::n()) |>
  tidyr::pivot_longer(cols = -c(id_cell), names_to = "vc", values_to = "effect") |>
  dplyr::mutate(effect = as.logical(effect))  
}


#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea assess all triads of interest from metaweb
#' @export
triads <- function(metaweb) {
  dat <- as.matrix(metaweb) |>
         motifcensus::motif_census_triplet() |>
         dplyr::mutate(motif_id = 1:dplyr::n())

  # Select only motifs of interest
  motifs <- c('exploitative competition','linear chain','apparent competition','omnivory')
  uid <- dat$name_uni %in% motifs
  dat <- dat[uid, ]

  # Modify index for positions i,j,k (currently starts at 0)
  dat$i_id <- dat$i + 1
  dat$j_id <- dat$j + 1
  dat$k_id <- dat$k + 1

  # Modify position numbers
  # Species:
  #   1: exploitative competition bottom
  #   2: exploitative competition top
  #   3: linear chains bottom
  #   4: linear chains middle
  #   5: linear chains top
  #   9: apparent competition bottom
  #  10: apparent competition top
  #  11: omnivory bottom
  #  12: omnivory middle
  #  13: omnivory top
  mid <- data.frame(motifcensus = c(1,2,3,4,5,9,10,11,12,13),
                    id = c(1,3,1,2,3,1,3,1,2,3))

  # Import in the table
  dat <- dat |>
         dplyr::left_join(mid, by = c("pos_i" = "motifcensus")) |>
         dplyr::rename(i_num = id) |>
         dplyr::left_join(mid, by = c("pos_j" = "motifcensus")) |>
         dplyr::rename(j_num = id) |>
         dplyr::left_join(mid, by = c("pos_k" = "motifcensus")) |>
         dplyr::rename(k_num = id)

  # Transform duplicated positions
  dat$j_num[dat$j_num == dat$i_num] <- 2
  dat$j_num[dat$j_num == dat$k_num] <- 2
  dat$k_num[dat$k_num == dat$i_num] <- 2  
           
  # Return 
  dat
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea pathways of direct effect 
#' @export
cea_pathways <- function(effect, vc) {
  # Binary effects (for pathways)
  bin <- cea_binary(effect)

  # vc as data.frame
  vc_df <- as.data.frame(vc) |>
           dplyr::select(-x, -y)

  # Index of vc
  vc_index <- data.frame(
    vc = colnames(vc_df), 
    vc_id = 1:ncol(vc_df)
  )              

  # Species in each cell and presence of direct effect
  vc_df <- vc_df |>
    dplyr::mutate(id_cell = 1:dplyr::n()) |>
    tidyr::pivot_longer(cols = -c(id_cell), names_to = "vc", values_to = "presence") |>
    tidyr::drop_na() |>
    dplyr::left_join(vc_index, by = c("vc")) |>
    dplyr::left_join(bin, by = c("id_cell","vc")) |>
    dplyr::select(-vc, -presence) 
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea pathways of indirect effect and trophic sensitivity
#' @export
ncea_pathways <- function(vc_id, motifs) {
  # Select triads involving species found in cell j
  uid <- motifs$i_id %in% vc_id$vc_id & 
         motifs$j_id %in% vc_id$vc_id & 
         motifs$k_id %in% vc_id$vc_id
  dat <- motifs[uid, ] 

  if (nrow(dat) > 0) {
    # Identify pathways of effect
    dat <- dat |>
      dplyr::mutate(i_focus = i_id, j_focus = j_id, k_focus = k_id) |>
      tidyr::pivot_longer(cols = c("i_focus","j_focus","k_focus"), values_to = "focus") |>
      dplyr::left_join(vc_id, by = c("focus" = "vc_id")) |>
      dplyr::select(-focus) |>
      tidyr::pivot_wider(names_from = name, values_from = effect) |>
      dplyr::mutate(path = 
       (i_focus + j_focus + k_focus) * 
       (i_focus*i_num + j_focus*j_num + k_focus*k_num)
      )
    
    # Extract trophic sensitivity per motif per species
    # NOTE: Trophic sensitivity values come from Beauchesne et al. 2021 and are available as data in
    #       the `rcea` package
    ts <- dat |>
      dplyr::mutate(i2 = i_id, j2 = j_id, k2 = k_id) |>
      tidyr::pivot_longer(
        cols = c("i_id","j_id","k_id","i_num","j_num","k_num"), 
        names_to = c("pos",".value"), 
        names_sep = "_"
      ) |>
      dplyr::left_join(
        trophic_sensitivity, 
        by = c("sum_pos" = "motifID", "num" = "speciesID", "path" = "pathID")
      ) |>
      tidyr::pivot_longer(cols = c("i2","j2","k2"), values_to = "interaction")
    
    # Join with direct effects table containing list of species 
    dat <- dplyr::full_join(vc_id, ts, by = c("vc_id" = "id")) |>
           dplyr::select(vc_id, interaction, Sensitivity)
    
    # Modify NAs for the name of the species, and a sensitivity of 1, meaning that species that are not involved in any motif should have a trophic sensitivity = to 1 
    uid <- is.na(dat$interaction)
    dat$interaction[uid] <- dat$vc_id[uid]
    dat$Sensitivity[uid] <- 1    
  } else {
    dat <- data.frame(
      vc_id = vc_id$vc_id,
      interaction = vc_id$vc_id,
      Sensitivity = 1
    )
  }
  
  # Return
  dat
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea apply ncea_pathways to get pathways of indirect effect and trophic sensitivity for all cells
#' @export
ncea_pathways_ <- function(direct_pathways, motifs) {
  dat <- dplyr::group_by(direct_pathways, id_cell) |>
         dplyr::group_split() |>
         lapply(function(x) ncea_pathways(x, motifs))
  for(i in 1:length(dat)) {
    dat[[i]] <- dplyr::mutate(
      dat[[i]],
      id_cell = i
    ) 
  }
  dat
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea get effects of drivers for species in all motifs in each cell
#' @export
ncea_motifs <- function(direct_effect, indirect_pathways) {
  drNames <- dimnames(direct_effect)[3][[1]]
  lapply(
    indirect_pathways,
    function(x) {
      cbind(x,direct_effect[x$id_cell[1], x$interaction, ])
    }
  ) |>
  dplyr::bind_rows() |>
  dplyr::select(
    id_cell, 
    vc_id, 
    interaction, 
    Sensitivity, 
    dplyr::all_of(drNames)
  ) |>
  dplyr::mutate(direct = vc_id == interaction) |>
  dplyr::group_by(id_cell, vc_id) |>
  dplyr::mutate(M = sum(direct)) |>
  dplyr::ungroup()  
}


#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea evaluate effects for all motifs using trophic sensitivity and effect weights
#' @export
ncea_effects <- function(motif_summary, w_d = 0.5, w_i = 0.25) {
  # w_d + 2*w_i = 1
  stopifnot(w_d + 2*w_i == 1)
  notDr <- c("id_cell","vc_id","interaction","Sensitivity","direct","M","weight")
  drNames <- colnames(motif_summary)
  drNames <- drNames[!drNames %in% notDr]
  
  # Direct & indirect weights 
  w <- data.frame(
    direct = c(TRUE, FALSE), 
    weight = c(w_d, w_i)
  )
  motif_summary <- dplyr::left_join(motif_summary, w, by = "direct")
  
  # Effects, i.e. multiply driver columns by trophic sensitivity and weights
  motif_effects <- motif_summary |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(drNames), 
        \(x) (x * weight * Sensitivity) / M
      )    
    ) |>
    dplyr::select(
      id_cell, M, vc_id, interaction, direct, dplyr::all_of(drNames), 
      -weight, -Sensitivity, -M
    )
  
  # Return
  motif_effects
}


#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea get contribution of species to indirect effects 
#' @export
get_species_contribution <- function(motif_effects) {
  notDr <- c("id_cell","vc_id","interaction","Sensitivity","direct","M","weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]
  dat = dplyr::filter(motif_effects, !direct) |>
  dplyr::group_by(id_cell, interaction) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(drNames),
      \(x) sum(x, na.rm = TRUE)
    )
  )
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea get direct and indirect effects of drivers
#' @export
get_direct_indirect <- function (motif_effects) {
  notDr <- c("id_cell","vc_id","interaction","Sensitivity","direct","M","weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]
  dplyr::group_by(motif_effects, id_cell, vc_id, direct) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(drNames), 
      \(x) sum(x, na.rm = TRUE)
    )
  )
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea get net effects of drivers
#' @export
get_net <- function (motif_effects) {
  notDr <- c("id_cell","vc_id","interaction","Sensitivity","direct","M","weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]
  dplyr::group_by(motif_effects, id_cell, vc_id) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(drNames), 
      \(x) sum(x, na.rm = TRUE)
    )
  )
}  

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn cea get effects per km2
#' @export
get_cekm <- function (motif_effects, vc) {
  # Driver names
  notDr <- c("id_cell","vc_id","interaction","Sensitivity","direct","M","weight")
  drNames <- colnames(motif_effects)
  drNames <- drNames[!drNames %in% notDr]

  # vc as data.frame
  vc_df <- as.data.frame(vc) |>
           dplyr::select(-x, -y)

  # Index of vc
  vc_index <- data.frame(
    vc = colnames(vc_df), 
    vc_id = 1:ncol(vc_df)
  )              
  
  # Calculate area, i.e. number of cells (assuming 1km2 grid cells)
  vc_df <- vc_df |>
           dplyr::mutate(id_cell = 1:dplyr::n()) |>
           tidyr::pivot_longer(cols = -c(id_cell), names_to = "vc", values_to = "presence") |>
           dplyr::group_by(vc) |>
           dplyr::summarise(km2 = sum(presence, na.rm = TRUE)) |>
           dplyr::left_join(vc_index, by = "vc") |>
           dplyr::select(-vc)
        
  # Direct & indirect effects      
  direct_indirect <- get_direct_indirect(motif_effects) |>
                     dplyr::left_join(vc_df, by = "vc_id") |>
                     dplyr::mutate(
                       dplyr::across(
                         dplyr::all_of(drNames),
                         \(x) x / km2
                       )
                     ) |>
                     dplyr::group_by(vc_id, direct) |>
                     dplyr::summarise(
                       dplyr::across(
                         dplyr::all_of(drNames),
                         \(x) sum(x, na.rm = TRUE)
                       )
                     ) |>
                     dplyr::mutate(type = 
                       dplyr::case_when(
                         direct ~ "direct",
                         !direct ~ "indirect"
                       ) 
                     ) |>
                     dplyr::select(-direct)
   
   # Net effects 
   net <- direct_indirect |>
          dplyr::select(-type) |>
          dplyr::group_by(vc_id) |>
          dplyr::summarise(
            dplyr::across(
              dplyr::all_of(drNames),
              \(x) sum(x, na.rm = TRUE)
            )
          ) |>
          dplyr::mutate(type = "net")
  
  
  # Summary table
  cekm <- dplyr::bind_rows(direct_indirect, net) |>
          dplyr::arrange(vc_id) |>
          dplyr::select(vc_id, type, dplyr::all_of(drNames))
          
  # Return
  cekm
}  
