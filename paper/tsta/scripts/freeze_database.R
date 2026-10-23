args <- commandArgs(TRUE)
root <- normalizePath(args[1])
.libPaths(c(file.path(root, "lib"), .libPaths()))
suppressPackageStartupMessages({library(SpatialESS); library(SpatialCellChat); library(jsonlite)})
stopifnot(as.character(packageVersion("SpatialESS")) == "0.1.4",
          normalizePath(find.package("SpatialESS")) == normalizePath(file.path(root, "lib/SpatialESS")))
e <- new.env()
data("CellChatDB.human", package="SpatialCellChat", envir=e)
saveRDS(e$CellChatDB.human, file.path(root, "inputs/database.rds"))
writeLines(capture.output(sessionInfo()), file.path(root, "reports/setup_sessionInfo.txt"))
write_json(list(SpatialESS_version=as.character(packageVersion("SpatialESS")),
                SpatialESS_path=find.package("SpatialESS"),
                SpatialCellChat_version=as.character(packageVersion("SpatialCellChat")),
                interaction_rows=nrow(e$CellChatDB.human$interaction),
                database_md5=unname(tools::md5sum(file.path(root, "inputs/database.rds")))),
           file.path(root, "reports/database.json"), pretty=TRUE, auto_unbox=TRUE)
cat("FROZEN database", nrow(e$CellChatDB.human$interaction), "LR rows; SpatialESS", as.character(packageVersion("SpatialESS")), "\n")
