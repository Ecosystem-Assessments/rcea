#' Network-scale cumulative effects assessments
#'
#' Assessment of cumulative effects and related metrics using the Beauchesne et al. 2021 method. 
#'
#' @eval arguments(c("drivers","vc","sensitivity","metaweb","trophic_sensitivity"))
#' @param motif_summary list of parameters obtained from \link{ncea_motifs} to allow for faster assessments, since that part of the assessment is the longest to run. 
#' @param exportAs string, the type of object that should be created, either a "data.frame" or a "stars" object
#'
#' @examples
#' # Data
#' drivers <- rcea:::drivers 
#' vc <- rcea:::vc
#' sensitivity <- rcea:::sensitivity
#' metaweb <- rcea:::metaweb
#' data(trophic_sensitivity)
#'
#' \dontrun{
#' # Network-scale effects 
#' beauchesne <- ncea(drivers, vc, sensitivity, metaweb, trophic_sensitivity)
#' plot(beauchesne$net)
#' plot(beauchesne$direct)
#' plot(beauchesne$indirect)
#' }
#' @export
#' @describeIn ncea_ effects assessment using the Beauchesne approach, i.e. assessment of direct and indirect effects of drivers on valued components (species), with method applied on a single species at a time for ease of parallelisation 
#' @export
ncea_ <- function(drivers, vc, sensitivity, metaweb, trophic_sensitivity, w_d = 0.5, w_i = 0.25, vc_focus, direct_effect = NULL, motif_summary = NULL) {
  # NOTE: `motif_summary` Allows a user to directly provide the motif_summary object if it's 
  #       been run before. This saves computation time if you wish to play around with weights 
  #       or other method parameters.
  if (is.null(motif_summary)) {
    # 3-species motifs for full metaweb
    motifs <- triads(metaweb)

    # Direct effects, i.e. Halpern approach
    direct_effect <- cea(drivers, vc, sensitivity) |>
                     make_array()

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
                          make_stars(drivers, vc)
  
  # Direct & indirect effects
  direct_indirect <- get_direct_indirect(motif_effects)
  direct <- dplyr::filter(direct_indirect, direct) |>
            make_stars(drivers, vc)
  indirect <- dplyr::filter(direct_indirect, !direct) |>
              make_stars(drivers, vc)

  # Net effects
  net <- get_net(motif_effects) |>
         make_stars(drivers, vc)
  
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


