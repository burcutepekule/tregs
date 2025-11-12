library(magick)
library(stringr)

folder <- "frames_5091_sterile_1_tregs_0"

# List all PNG files
png_files <- list.files(folder, pattern = "\\.png$", full.names = TRUE)

# Extract numeric index from filename (the number before .png)
frame_numbers <- as.numeric(str_extract(png_files, "(?<=_trnd_0_)\\d+"))

# Sort by numeric order
png_files_sorted <- png_files[order(frame_numbers)]

png_files_sorted = png_files_sorted[1:25]

# Read and combine into a PDF
imgs <- image_read(png_files_sorted)
pdf_path <- file.path(folder, "combined_frames.pdf")
image_write(image_join(imgs), path = pdf_path, format = "pdf")

cat("✅ Saved combined PDF at:", pdf_path, "\n")
