% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncea_species.R
\name{ncea_species}
\alias{ncea_species}
\title{Network-scale cumulative effects assessments for a single focal species}
\usage{
ncea_species(
  focus,
  drivers,
  vc,
  sensitivity,
  metaweb,
  trophic_sensitivity,
  w_d = 0.5,
  w_i = 0.25,
  output = "output/ncea",
  output_format = "geotiff"
)
}
\arguments{
\item{focus}{named argument, string with name of valued component on which to perform the assessment.}

\item{drivers}{distribution and intensity of environmental drivers as stars object}

\item{vc}{distribution of valued components as stars object}

\item{sensitivity}{matrix of environmental drivers and valued component, with same name as those used in \code{drivers} and \code{vc}}

\item{metaweb}{matrix of valued component by valued component describing the binary interations structuring the network of valued components}

\item{trophic_sensitivity}{data.frame of trophic sensitivities, default from Beauchesne. Available as data package with \code{data(trophic_sensitivity)}}

\item{w_d, w_i}{weight for the direct (\code{w_d}) and indirect (\code{w_i}) modules when calculating network-scale cea scores; w_d + 2*w_i should be equal to 1.}

\item{output}{relative path to export results of assessment.}

\item{output_format}{output format, one of \code{geotiff} or \code{COG}.}
}
\description{
Assessment of cumulative effects and related metrics using the Beauchesne et al. 2021 method.
}
\examples{
# Data
drivers <- rcea:::drivers
vc <- rcea:::vc
sensitivity <- rcea:::sensitivity
metaweb <- rcea:::metaweb
trophic_sensitivity <- rcea::trophic_sensitivity

\dontrun{
# Network-scale effects for individual species
ncea_species(focus = "vc1", drivers, vc, sensitivity, metaweb, trophic_sensitivity)
ncea_species(focus = "vc1", drivers, vc, sensitivity, metaweb, trophic_sensitivity, output_format = "COG")
}
}
