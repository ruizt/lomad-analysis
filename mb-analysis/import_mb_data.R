# data-raw/import_mb_data.R
#
# Copies clean Morro Bay data files from the mb-qartod repository into
# _mb-data/ in this repository. Re-run whenever mb-qartod data is updated.
#
# Assumes mb-qartod and lomad are siblings under the same parent directory.

mb_clean <- fs::path(here::here(), "../mb-qartod/_data/clean")

if (!fs::dir_exists(mb_clean)) {
  stop("mb-qartod clean data directory not found at: ", mb_clean,
       "\nCheck that the mb-qartod repo is a sibling of this one.")
}

out_dir <- here::here("_mb-data")
fs::dir_create(out_dir)

fs::file_copy(
  fs::dir_ls(mb_clean),
  fs::path(out_dir, fs::path_file(fs::dir_ls(mb_clean))),
  overwrite = TRUE
)

message("Copied ", length(fs::dir_ls(mb_clean)), " file(s) to ", out_dir)
