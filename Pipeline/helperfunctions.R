

get_points_inside_polygons <- function(raster, polygons_sf, join = st_within, left = FALSE) {
  names(raster) <- "value"
  points <- st_as_sf(as.points(raster))
  points_inside <- st_join(points, polygons_sf, join = join, left = left) |> as_tibble() |> rename(point = geometry)
  return(points_inside)
}



compute_buffer_pix <- function(um, pixel_um_ratio) {
  return(um / pixel_um_ratio) 
}


weigthed_sd <- function(value, weight, na.rm = TRUE) {
  sum(weight * ((value - weighted.mean(value, weight, na.rm = na.rm))^2)) / (length(weight)-1/length(weight)) * sum(weight)
}

median_mad_threshold <- function(x, w, n, constant = 1.482602, na.rm = TRUE){
  w_median  <- robsurvey::weighted_median(x, w, na.rm = na.rm)
  w_mad <- robsurvey::weighted_mad(x, w, constant = constant, na.rm = na.rm)
  return(w_median + n * w_mad)
}

image_pixels_in_mesh <- function(image, mesh, include_cols = NULL) {
  
  img_frames <- imager::depth(image)
  if(img_frames == 1){
    
    
    return(
      terra::rast(t(as.matrix(image))) |>
        terra::flip(direction = "vertical") |>
        exactextractr::exact_extract(
          mesh,
          include_xy = TRUE,
          force_df = TRUE,
          progress = FALSE,
          include_cols = include_cols) |> bind_rows()  
    )
  }
  
  return(
    purrr::map2_dfr(
      1:img_frames,
      imager::imsplit(image, "z"),
      \(i, img) {image_pixels_in_mesh(img, mesh, include_cols = include_cols)},
      .id = "frame") 
  )
}

points_in_mesh <- function(point, mesh) {
  st_intersection(point, mesh) -> intersection
  if (dim(intersection)[1]>0){
    return( 
      intersection |>
        sfheaders::sf_to_df(fill = TRUE) |> 
        as_tibble())
  }
  NULL
}




extract_image_pixels_from_meshes <- function(images_paths, meshes, include_cols = NULL, parallel = TRUE) {
  
  if(length(images_paths) != length(meshes)){
    stop("length(images_paths) must be equal to (==) length(meshes)")
  }
  

  if(length(images_paths) == 1){
    
    print(glue::glue("Reading {fs::path_file(file)}..."))
    image <- imager::load.image(images_paths) 
    return(
      image_pixels_in_mesh(image, meshes, include_cols = c("id")) |>
        mutate(frame = as.integer(frame))
    )
  }
  
  if (is.numeric(parallel)) {
    plan(multisession, workers = parallel)
    map_function <- furrr::future_map2_dfr
  }
  else if (isTRUE(parallel)) {
    plan(multisession, workers = future::availableCores()-2)
    map_function <- furrr::future_map2_dfr
  }
  else {
    map_function <- purrr::map2_dfr
  }
  
  map_function(images_paths, meshes, \(x, y) {extract_image_pixels_from_meshes(x, y, include_cols, parallel = FALSE)})

}


running_means_difference <- function(x, span) {
  runner::runner(x, k = span, mean, na_pad = F, lag = 0) - runner::runner(x, k = span, mean, na_pad = F, lag = -span+1)
}


circle_area <- function(radius) { pi * radius^2}
sphere_volume <-  function(radius) {4/3 * pi * radius^3}

find_min <- function(vector, for_which) {min(vector[which(for_which)])}
find_max <- function(vector, for_which) {max(vector[which(for_which)])}


get_mesh_list <- function(meshes, mesh_col, mesh_id_col, split_by){
   active_geometry <- st_geometry(meshes)
   st_geometry(meshes) <- deparse(substitute(mesh_id_col))
   
   meshes |>
   select(all_of(mesh_col, mesh_id_col, split_by)) |>
   arrange(split_by) |>
   group_by(split_by) |>
   group_split() -> l_meshes
   
   st_geometry(meshes) <- active_geometry
   return(l_meshes)
   
  
}

get_image_list <- function(meshes, image_path_col) {
  (meshes |> arrange(image_path_col))[[deparse(substitute(image_path_col))]] |> unique()
}


