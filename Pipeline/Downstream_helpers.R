

plot_pixel_intensities <- function(df, geometry = "mesh_buffered", channel = "phase", type = "wrap", threshold = FALSE) {
  channel_filename <- switch (channel,
                              "phase" = "phase_image_path",
                              "focis" = "focis_image_path",
                              "diffuse" = "diffuse_image_path"
  )
  low <- "#000000"
  high <- "#FFFFFF"
  st_geometry(df) <- geometry
  
  df |>
    mutate(thresh = threshold) |>
    group_by(id) |>
    select(all_of(geometry), all_of(channel_filename), id, thresh) |> group_split() -> df_splitted
  purrr::map_dfr(df_splitted, \(mesh){
    print(glue::glue("Loading {mesh[[channel_filename]]}"))
    img <- imager::load.image(mesh[[channel_filename]])
    image_pixels_in_mesh(img, mesh, include_cols = c("id", "thresh")) |>
      mutate(frame = as.numeric(frame)) 
  }) -> intensity
  
  if(!isFALSE(threshold)){
    intensity |> 
      ungroup() |>
      mutate(value = if_else(value >= thresh,1,0)) -> intensity
    low <- "#f48579"
    high <- "#33c4cb"
  }
  
  int_tmp <<- intensity
  
  if (type == "movie_frame"){
    intensity |>
      group_by(id, frame) |>
      group_split() |>
      purrr::map(\(df) {
        df |>  ggplot(aes(x = x, y = y, color = value, fill = value), width = 1, height = 1) +
          geom_tile() +
          scale_fill_gradient(low = low, high = high) +
          scale_color_gradient(low = low, high = high) +
          theme_void() +
          labs(title = df$id[[1]]) + 
          coord_fixed() + theme(legend.position = "none")
      })
  }
  
  
  else if (type == "movie") {
    intensity |>
      group_by(id) |>
      group_split() |>
      purrr::map(\(df) {
        df |>  ggplot(aes(x = x, y = y, color = value, fill = value), width = 1, height = 1) +
          geom_tile() +
          gganimate::transition_time(frame) +
          scale_fill_gradient(low = low, high = high) +
          scale_color_gradient(low = low, high = high) +
          theme_void() +
          ggtitle("{formatC((frame*8)%/%60, width = 2, format = 'd', flag = '0')}:{formatC((frame*8)%%60, width = 2, format = 'd', flag = '0')}") + 
          coord_fixed() +theme(legend.position = "none")
      })
  }
  
  
  else if (type == "wrap") {
    intensity |>
      group_by(id) |>
      group_split() |>
      purrr::map(\(df) {
        df |>  ggplot(aes(x = x, y = y, color = value, fill = value), width = 1, height = 1) +
          geom_tile() +
          gganimate::transition_time(frame) +
          scale_fill_gradient(low = low, high = high) +
          scale_color_gradient(low = low, high = high) +
          theme_void() +
          labs(title = df$id[[1]]) + 
          coord_fixed() + theme(legend.position = "none")
      })
  }
  
}
