rmp_2 <- function (x=dat, minimum.coverage = 0.5, score.estimator = "median", 
          annotation = 2, module.db = NULL, threads = 1, normalize.by.length = FALSE, 
          distribute = FALSE, java.mem = NULL, out.dir=NULL) 
{
  rpm.exec <- system.file("java", "omixer-rpm.jar", package = "omixerRpm")
  if (is.null(module.db)) {
    module.db <- loadDefaultDB()
  }
  out.dir <- out.dir
  dir.create(out.dir)
  if (!is.character(x)) {
    in.dir <- out.dir
    dir.create(in.dir)
    input <- file.path(in.dir, "input.tsv")
    write.table(x, input, col.names = T, row.names = F, 
                quote = F, sep = "\t")
  }
  else {
    input <- x
  }
  command <- "java -server"
  if (!is.null(java.mem)) {
    command <- paste(command, paste0("-Xmx", java.mem, "G"))
  }
  command <- paste(command, "-jar", rpm.exec, "-c", minimum.coverage, 
                   "-s", score.estimator, "-d", file.path(module.db@directory, 
                                                          module.db@modules), "-i", input, "-o", out.dir, 
                   "-a", annotation, "-t", threads, "-e", 2)
  if (normalize.by.length == TRUE) {
    command <- paste(command, "-n")
  }
  if (distribute == TRUE) {
    command <- paste(command, "--Xdistribute")
  }
  tryCatch({
    system(command)
  }, warning = function(e) {
    print(e)
    if (e$message == "error in running command") {
      stop(e)
    }
  }, error = function(e) {
    print(geterrmessage())
    stop(e)
  })
  abundance <- read.table(file.path(out.dir, "modules.tsv"), 
                          sep = "\t", header = TRUE)
  coverage <- read.table(file.path(out.dir, "modules-coverage.tsv"), 
                         sep = "\t", header = TRUE)
  annotation.df <- NULL
  if (annotation == 1) {
    abundance.colnames <- colnames(abundance)
    annotation.df <- as.data.frame(abundance[, 1])
    abundance <- as.data.frame(abundance[, -c(1)])
    coverage <- as.data.frame(coverage[, -c(1)])
    colnames(annotation.df) <- abundance.colnames[1]
    colnames(abundance) <- abundance.colnames[2:length(abundance.colnames)]
    colnames(coverage) <- abundance.colnames[2:length(abundance.colnames)]
  }
  else if (annotation == 2) {
    abundance.colnames <- colnames(abundance)
    annotation.df <- as.data.frame(abundance[, c(1, 2)])
    abundance <- as.data.frame(abundance[, -c(1, 2)])
    coverage <- as.data.frame(coverage[, -c(1, 2)])
    colnames(annotation.df) <- abundance.colnames[1:2]
    colnames(abundance) <- abundance.colnames[3:length(abundance.colnames)]
    colnames(coverage) <- abundance.colnames[3:length(abundance.colnames)]
  }
  Modules(abundance = abundance, coverage = coverage, annotation = annotation.df, 
          db = module.db)
}
