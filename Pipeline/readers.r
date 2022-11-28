#################################################################
# Thomas Lamot, 2022
# Used in Santin et al, 2022

## The functions below depend on the following functions in extract_from.R: extract_oufti_celllist(), extract_microbeJ_mesh() & extract_microbeJ_cell()

## Dependencies: R.matlab, tidyverse

## Oufti

read_oufti <- function(oufti_matfile) {
  oufti_matfile |>
    R.matlab::readMat(oufti_matfile) |>
    extract_oufti_celllist() -> cell_list

    
  return(cell_list)
}

read_oufti_parameters <- function(oufti_matfile) {
  oufti_matfile |>
    R.matlab::readMat(oufti_matfile) -> raw
  raw$p[,,1] |>  unlist()  |> as.list() -> param
  return(param)
}

##MicrobeJ

read_microbeJ_mesh <- function(microbeJ_csv){
  microbeJ_df <- read_csv(microbeJ_csv)
  microbeJ_df |> 
    extract_microbeJ_mesh() -> mesh
  return(mesh)
}

read_microbeJ_cell <- function(microbeJ_csv){
  microbeJ_df <- read_csv(microbeJ_csv) |>
    extract_microbeJ_cell()
  return(microbeJ_df)
}

##################################################################

#Functions for other input programs; not tested in this context:

#Supersegger

read_supersegger_cell <- function(supersegger_matfile) {
  supersegger_matfile |> 
    R.matlab::readMat() |>
    extract_supersegger_celllist() -> cell_list
  
  return(cell_list)
}

read_supersegger_mesh_ <- function(supersegger_matfile) {
  supersegger_matfile |>
    R.matlab::readMat() |>
    extract_supersegger_mesh() -> mesh
  return(mesh)
}

read_supersegger_mesh <- function(supersegger_path) {
  supersegger_path <- fs::path(supersegger_path)
  if (fs::is_dir(supersegger_path)) {
    mesh_files <- fs::dir_ls(supersegger_path, glob = "*.mat")
    purrr::map(mesh_files, \(x) read_supersegger_mesh_(x)) |>
      bind_rows() -> meshes
    return(meshes)
  }
  else if (fs::is_file(supersegger_path) && fs::path_ext(supersegger_path) == "mat") {
    return(read_supersegger_mesh_(supersegger_path))
  }
  stop("Invalid file format")
}

#Morphometrics

read_morphometrics <- function(morphometrics_matfile){
  morphometrics_list <- R.matlab::readMat(morphometrics_matfile)
  cell_list <- extract_morphometrics(morphometrics_list)
  return(cell_list)
}
