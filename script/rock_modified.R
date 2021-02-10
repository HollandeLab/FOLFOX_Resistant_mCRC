
get_caller <-  function (cnv) {
  stopifnot(file.exists(cnv))
  h <- readr::read_lines(cnv, n_max = 1)
  caller <- NULL
  if (grepl("probes", h)) {
    caller <- "cnvkit"
  }
  else if (grepl("tcn\\.em", h)) {
    caller <- "facets"
  }
  else if (grepl("segmentStartSupport", h)) {
    caller <- "purple"
  }
  else if (grepl("TITAN", h)) {
    caller <- "titan"
  }
  else if (grepl("ASCAT", h)) {
    caller <- "ascat"
  }
  else {
    message(
      "Unknown caller. Make sure you're reading CNV segments ",
      "from\nCNVkit (cnvkit-call.cns),\nFACETS (facets_segs.tsv),\n",
      "PURPLE (purple.cnv), \nTITAN (titan_segs.tsv), or \nASCAT (ascat.genes.known.txt)."
    )
  }
  return(caller)
}





prep_cnv_circos <- function (cnv)
{
  cnv <- read_cnv(cnv)$cnv
  cnv %>% dplyr::mutate(chrom = paste0("hs", .data$chrom),
    tot_cn = .data$tot_cn - 2) %>% dplyr::rename(value = .data$tot_cn)
}




read_cnv <- function (cnv){
  caller <- get_caller(cnv)
  stopifnot(!is.null(caller))
  fl <- list(
    facets = prep_facets_seg,
    ascat = prep_ascat_seg
    # ,
    # cnvkit = prep_cnvkit_seg,
    # purple = prep_purple_seg,
    # titan = prep_titan_seg
  )
  fl[[caller]](cnv)
}





circos_prep <-  function (outdir = "circos", manta = NULL, cnv = NULL) 
{
  template <- NULL
  stopifnot((!is.null(manta) && file.exists(manta)) || (!is.null(cnv) && 
      file.exists(cnv)))
  dir.create(outdir, recursive = TRUE)
  message(glue::glue("Exporting Manta and/or CNV circos files to '{outdir}'."))
  if (!is.null(manta)) {
    template <- "sv"
    manta <- prep_manta_vcf2(manta)
    readr::write_tsv(manta, file.path(outdir, "SAMPLE.link.circos"), 
      col_names = FALSE)
  }
  if (!is.null(cnv)) {
    template <- paste0("cnv", template)
    cnv <- prep_cnv_circos(cnv)
    readr::write_tsv(cnv, file.path(outdir, "SAMPLE.cnv.circos"), 
      col_names = FALSE)
  }
  message(glue::glue("Copying circos templates to '{outdir}'. 'template' is {template}."))
  file.copy(from = system.file("templates/circos", glue::glue("circos_{template}.conf"), 
    package = "pebbles"), to = file.path(outdir, "circos.conf"), 
    overwrite = TRUE)
  file.copy(system.file("templates/circos", "gaps.txt", package = "pebbles"), 
    outdir, overwrite = TRUE)
  file.copy(system.file("templates/circos", "ideogram.conf", 
    package = "pebbles"), outdir, overwrite = TRUE)
  invisible(list(sv = manta, cnv = cnv))
}








