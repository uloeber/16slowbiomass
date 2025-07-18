# ggsave in different formats with current date at once

ggsaveAll <- function(filepath, 
                      filename, 
                      plot, 
                      devices = "svg", 
                      date = sub("-", "", Sys.Date()), 
                      width = NA,
                      height = NA) {
  for (dev in devices) {
    ggsave(
      filename = paste0(filepath, date, "_", filename, ".", dev),
      plot = plot,
      device = dev,
      width = width,
      height = height
    )
  }
}
