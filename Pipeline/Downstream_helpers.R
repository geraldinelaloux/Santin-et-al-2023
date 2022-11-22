

plot_pixel_intensities <- function(df, geometry = "mesh_buffered", channel = "phase", type = "wrap", threshold = FALSE) {
  channel_filename <- switch (channel,
                              "phase" = "phase_image_path",
                              "focis" = "focis_image_path",
                              "diffuse" = "diffuse_image_path"
  )
  
  st_geometry(df) <- geometry
  
  df |>
    group_by(id) |>
    select(all_of(geometry), all_of(channel_filename), id) |> group_split() -> df_splitted
  purrr::map_dfr(df_splitted, \(mesh){
    print(glue::glue("Loading {mesh[[channel_filename]]}"))
    img <- imager::load.image(mesh[[channel_filename]])
    image_pixels_in_mesh(img, mesh, include_cols = c("id")) |>
      mutate(frame = as.numeric(frame))
  }) -> intensity
  
  if(!isFALSE(threshold)){
    intensity |> 
      mutate(value = if_else(value >= threshold,0,1)) -> intensity
  }
  
  int_tmp <<- intensity
  
  if (type == "movie_frame"){
    intensity |>
      group_by(id, frame) |>
      group_split() |>
      purrr::map(\(df) {
        df |>  ggplot(aes(x = x, y = y, color = value, fill = value), width = 1, height = 1) +
          geom_tile() +
          scale_fill_gradient(low = "#000000", high = "#FFFFFF") +
          scale_color_gradient(low = "#000000", high = "#FFFFFF") +
          theme_classic() +
          labs(title = df$id[[1]]) + 
          coord_fixed()
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
          scale_fill_gradient(low = "#000000", high = "#FFFFFF") +
          scale_color_gradient(low = "#000000", high = "#FFFFFF") +
          theme_classic() +
          ggtitle(" {frame*8} / {nframes*8}") + 
          coord_fixed()
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
          scale_fill_gradient(low = "#000000", high = "#FFFFFF") +
          scale_color_gradient(low = "#000000", high = "#FFFFFF") +
          theme_classic() +
          labs(title = df$id[[1]]) + 
          coord_fixed()
      })
  }
  
  else {
    intensity |>
      group_by(id) |>
      group_split() |>
      purrr::map(\(df) {
        df |>  ggplot(aes(x = x, y = y, color = value, fill = value), width = 1, height = 1) +
          geom_tile() +
          facet_wrap(~frame) +
          scale_fill_gradient(low = "#000000", high = "#FFFFFF") +
          scale_color_gradient(low = "#000000", high = "#FFFFFF") +
          theme_classic() +
          labs(title = {df$id[[1]]}) + 
          coord_fixed()
      })
  }
}