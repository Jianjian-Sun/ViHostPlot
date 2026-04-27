#' Read Host Chromosome Sizes
#'
#' @param path Path to a tab-delimited host chromosome-size table.
#' @return A data frame with columns \code{chr}, \code{start}, and \code{end}.
#' @keywords internal
read_chrom_sizes <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  df <- utils::read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c("chr", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Host chromosome file must contain columns: chr, start, end")
  }

  df$chr <- trimws(gsub("(?i)^chr", "", as.character(df$chr)))
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  df <- df[!is.na(df$chr) & !is.na(df$start) & !is.na(df$end), , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No valid chromosome sizes were found.")
  }

  df
}

resolve_ucsc_assembly <- function(host) {
  host <- tolower(trimws(as.character(host)))

  if (length(host) != 1 || is.na(host) || !nzchar(host)) {
    stop("host must be a single non-empty string.")
  }

  aliases <- c(
    human = "hg38",
    homo_sapiens = "hg38",
    mouse = "mm39",
    mus_musculus = "mm39",
    rat = "rn7",
    rattus_norvegicus = "rn7",
    zebrafish = "danrer11",
    danio_rerio = "danRer11",
    fruitfly = "dm6",
    drosophila = "dm6",
    yeast = "sacCer3"
  )

  if (host %in% names(aliases)) {
    return(unname(aliases[[host]]))
  }

  host
}

filter_primary_chromosomes <- function(df, host, assembly = NULL) {
  host_key <- tolower(trimws(as.character(host)))
  assembly_key <- if (is.null(assembly)) "" else tolower(trimws(as.character(assembly)))

  if (host_key %in% c("human", "homo_sapiens") || assembly_key == "hg38") {
    keep_chr <- c(as.character(1:22), "X", "Y")
    return(df[df$chr %in% keep_chr, , drop = FALSE])
  }
  # 后续这里是否可以添加更多的宿主
  
  keep_generic <- grepl("^[0-9]+$|^X$|^Y$", df$chr, ignore.case = TRUE)
  df <- df[keep_generic, , drop = FALSE]

  if (nrow(df) == 0) {
    stop(
      "No primary chromosomes remained after filtering UCSC chromosome sizes for host: ",
      host,
      call. = FALSE
    )
  }

  df
}

fetch_ucsc_chrom_sizes <- function(host) {
  assembly <- resolve_ucsc_assembly(host)
  urls <- c(
    sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes", assembly, assembly),
    sprintf("https://hgdownload.soe.ucsc.edu/goldenpath/%s/bigZips/%s.chrom.sizes", assembly, assembly),
    sprintf("https://hgdownload.soe.ucsc.edu/goldenpath/%s/bigZips/latest/%s.chrom.sizes", assembly, assembly)
  )

  last_error <- NULL
  for (u in urls) {
    fetched <- tryCatch(
      utils::read.table(
        u,
        header = FALSE,
        sep = "\t",
        stringsAsFactors = FALSE,
        comment.char = "",
        quote = ""
      ),
      error = function(e) {
        last_error <<- e$message
        NULL
      }
    )

    if (!is.null(fetched) && ncol(fetched) >= 2) {
      df <- data.frame(
        chr = trimws(gsub("(?i)^chr", "", as.character(fetched[[1]]), perl = TRUE)),
        start = 0,
        end = as.numeric(fetched[[2]]),
        stringsAsFactors = FALSE
      )
      df <- df[!is.na(df$chr) & nzchar(df$chr) & !is.na(df$end), , drop = FALSE]
      df <- filter_primary_chromosomes(df, host = host, assembly = assembly)

      if (nrow(df) > 0) {
        return(df)
      }
    }
  }

  stop(
    "Unable to fetch chromosome sizes online for host/assembly: ", host,
    ". Tried UCSC assembly: ", assembly,
    if (!is.null(last_error)) paste0(". Last error: ", last_error) else "",
    call. = FALSE
  )
}

resolve_host_chrom_sizes <- function(host = NULL, chrom_file = NULL) {
  if (!is.null(host) && !is.null(chrom_file)) {
    stop("Please supply either host or chrom_file, not both.", call. = FALSE)
  }

  if (!is.null(host)) {
    return(fetch_ucsc_chrom_sizes(host))
  }

  if (!is.null(chrom_file)) {
    return(read_chrom_sizes(chrom_file))
  }

  stop("Please supply either a built-in host name or a host chromosome file.", call. = FALSE)
}
