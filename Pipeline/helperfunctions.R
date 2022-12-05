

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
      .id = "frame")  |> mutate(frame =as.integer(frame))
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


st_polygon_autoclose <- function(x, y) {
  if (first(x) != last(x) || first(y) != last(y)) {
    x <- append(x, first(x))
    y <- append(y, first(y))
  }
  return(sf::st_polygon(list(cbind(x = x, y = y))))
}
 
angle_from_horizontal <- function(..., type = "rad") {
  result <- switch(type,
                   "rad" = pi,
                   "deg" = 180,
                   "hour" = 12
  )
  angles <- c(...)
  for (theta in angles) {
    result <- result - theta
  }
  return(result)
}

angle_from_vertical <- function(..., type = "rad") {
  result <- switch(type,
                   "rad" = pi / 2,
                   "deg" = 90,
                   "hour" = 6
  )
  angles <- c(...)
  for (theta in angles) {
    result <- result - theta
  }
  return(result)
}

polar_distance <- function(x, y) {
  return(sqrt(x^2 + y^2))
}

polar_angle <- function(x, y) {
  return(atan(x / y))
}

st_box <- function(left, right, top, bot, center, width, height, angle = 0) {
  if (!missing(center)) {
    if (!missing(width) & !missing(height)) {
      left <- center[1] - width/2
      right <- center[1] +  width/2
      bot <- center[2] - height/2
      top <- center[2] + height/2
    }
    else if (!missing(left) & !missing(right)) {
      right <- center[1] + (center-left)/2
      top <- center[2] + (center-bot)/2
    }
    else {
      stop("you need to provide 'width' and 'height' OR 'left' and 'bot' arguments when using 'center' argument")
    }
  }
  else {
    
    if (!xor(missing(top) | missing(right), missing(width) | missing(height))) {
      stop("you need to provide 'width' and 'height' OR 'top' arguments when using 'left' and bot argument")
    }
    if (missing(top) | missing(right)) {
      right <- left +  width/2
      top <- bot +  height/2
    }
  }
  
  points <- matrix(c(right, top, left, top, left, bot, right, bot, right, top), ncol = 2, byrow = TRUE)
  ox <- (right + left)/2
  oy <- (top + bot)/2
  points <- matrix(c(
    cos(-angle) * (points[,1] - ox)  - sin(-angle) * (points[,2] - oy) + ox,
    sin(-angle) * (points[,1] - ox)  + cos(-angle) * (points[,2] - oy) + oy), ncol = 2)
  box <- st_polygon(list(points))
  return(box)
}


extract_image_pixels_from_meshes <- function(images_paths, meshes, include_cols = NULL, parallel = TRUE) {

  if(length(images_paths) == 1){
    
    print(glue::glue("Reading {fs::path_file(images_paths)}..."))
    
    image <- imager::load.image(images_paths) 
    return(
      image_pixels_in_mesh(image, meshes, include_cols = c("id")) |>
        mutate(frame = as.integer(frame))
    )
  }
  
  if(length(images_paths) != length(meshes)){
    stop("length(images_paths) must be equal to (==) length(meshes)")
  }
  
  if (is.numeric(parallel)) {
    plan(multisession, workers = parallel)
    map_function <- furrr::future_map2_dfr
  }
  else if (isTRUE(parallel)) {
    plan(multisession, workers = future::availableCores()-4)
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
  active_geometry <- attr(meshes, "sf_column")
  st_geometry(meshes) <- deparse(substitute(mesh_id_col))
   
  vars <- c(
    deparse(substitute(mesh_col)),
    deparse(substitute(mesh_id_col)),
    deparse(substitute(split_by)))
     
     meshes |>
     as_tibble() |>
     select(all_of(vars)) |>
     st_sf() |>
     arrange({{split_by}}) |>
     group_by({{split_by}}) |>
     group_split() -> l_meshes
     
     st_geometry(meshes) <- active_geometry
     return(l_meshes)
     
}

get_image_list <- function(meshes, image_path_col) {
  meshes[[deparse(substitute(image_path_col))]] |> unique() |> sort()
}


