#' Run BLASTP with explicit binary and DB paths; hermetically scoped
#'
#' @param query_faa Path to query FASTA/FAA file
#' @param db BLAST database prefix (e.g., "swissprot" or absolute "/path/to/swissprot")
#' @param dbdir Directory containing BLAST DBs; if db is relative, dbdir is prepended
#' @param blastp_bin Full path to blastp binary (preferred explicit control)
#' @param blastp_dir Directory containing blastp; used if blastp_bin is NULL
#' @param evalue E-value threshold
#' @param threads Integer number of threads
#' @param fields Vector of outfmt 6 columns
#' @param validate_db Logical; check DB via blastdbcmd before running
#' @param more_args Character vector of additional BLASTP args (e.g., c("-max_target_seqs","10"))
#' @return tibble of hits; aborts with rich diagnostics on failure
#' @export
blastp_capture <- function(
    query_faa,
    db = "swissprot",
    dbdir = Sys.getenv("BLASTDB", ""),
    blastp_bin = Sys.getenv("BLASTP_BIN", ""),
    blastp_dir = Sys.getenv("BLASTP_PATH", ""),
    evalue = 1e-10,
    threads = parallel::detectCores(),
    fields = c("qseqid","sacc","stitle","pident","length","qcovs","bitscore","evalue","staxids"),
    validate_db = TRUE,
    more_args = character()
) {
  # Dependencies: avoid attachment
  if (!requireNamespace("cli", quietly = TRUE)) stop("Package 'cli' is required.")
  if (!requireNamespace("readr", quietly = TRUE)) stop("Package 'readr' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")

  # --- Resolve blastp binary deterministically ---
  resolve_blastp <- function(bin, dir) {
    if (nzchar(bin)) {
      if (file.exists(bin) && isTRUE(file.info(bin)$isdir)) file.path(bin, "blastp") else bin
    } else if (nzchar(dir)) {
      if (file.exists(dir) && isTRUE(file.info(dir)$isdir)) file.path(dir, "blastp") else dir
    } else {
      Sys.which("blastp")
    }
  }
  candidate_bin <- resolve_blastp(blastp_bin, blastp_dir)

  if (!nzchar(candidate_bin) || !file.exists(candidate_bin)) {
    cli::cli_abort(c(
      "Could not locate {.val blastp}.",
      i = "Set {.env BLASTP_BIN} to the full binary or {.env BLASTP_PATH} to its directory, ",
      i = "or ensure 'blastp' is on PATH."
    ))
  }
  if (isTRUE(file.info(candidate_bin)$isdir)) {
    cli::cli_abort(c(
      "Provided path is a directory, not a binary: {.path {candidate_bin}}",
      i = "The wrapper accepts a directory, but it must contain a {.file blastp} executable."
    ))
  }
  if (.Platform$OS.type == "unix" && file.access(candidate_bin, 1) != 0) {
    cli::cli_abort("blastp binary is not executable: {.path {candidate_bin}}")
  }
  blastp <- candidate_bin

  # Resolve blastdbcmd alongside blastp; fallback to PATH
  blastdbcmd <- file.path(dirname(blastp), "blastdbcmd")
  if (!file.exists(blastdbcmd)) {
    blastdbcmd <- Sys.which("blastdbcmd")
  }

  # --- Preflight inputs ---
  if (!file.exists(query_faa)) {
    cli::cli_abort("Query file not found: {.path {query_faa}}")
  }

  # Derive db_path: if db is relative and dbdir provided, join; otherwise absolute db is used verbatim.
  db_path <- if (nzchar(dbdir) && !grepl("^/", db)) file.path(dbdir, db) else db

  # Warn if using relative path without explicit dbdir
  if (!nzchar(dbdir) && !grepl("^/", db)) {
    cli::cli_warn(c(
      "Using relative database path {.val {db}} without explicit {.arg dbdir}.",
      i = "Set {.env BLASTDB} or provide {.arg dbdir} parameter for reliable database location.",
      i = "Current working directory: {.path {getwd()}}"
    ))
  }

  # --- Optional DB sanity check (quote the path; out-of-band shell may split otherwise) ---
  if (validate_db && nzchar(blastdbcmd)) {
    db_path_q <- shQuote(db_path, type = "sh")
    info <- tryCatch(
      system2(blastdbcmd, args = c("-db", db_path_q, "-info"), stdout = TRUE, stderr = TRUE),
      error = function(e) character()
    )
    if (length(info) == 0 || any(grepl("\\bError\\b", info))) {
      cli::cli_abort(c(
        "Could not open BLAST DB at {.path {db_path}}.",
        ">" = paste(info, collapse = "\n"),
        i = "Provide an absolute DB prefix or set {.env BLASTDB} correctly."
      ))
    }
  } else if (validate_db) {
    # Fallback validation: check for common DB file extensions without blastdbcmd
    cli::cli_warn("blastdbcmd not found; attempting manual database file check.")

    # Determine search directory and database basename
    db_dir <- dirname(db_path)
    if (db_dir == ".") db_dir <- getwd()
    db_base <- basename(db_path)

    # Check for protein database files (.phr, .pin, .psq) or nucleotide (.nhr, .nin, .nsq)
    db_files <- list.files(
      path = db_dir,
      pattern = paste0("^", gsub("([.+*?^$(){}|\\[\\]\\\\])", "\\\\\\1", db_base),
                       "\\.(phr|pin|psq|nhr|nin|nsq)$"),
      full.names = FALSE
    )

    if (length(db_files) == 0) {
      # Provide helpful suggestions
      common_location <- "/home/william-ackerman/Desktop/Link to Desktop/tmp_ncbi_blast/"
      common_dbs <- character()
      if (file.exists(common_location)) {
        common_dbs <- list.files(
          path = common_location,
          pattern = "\\.(phr|nhr)$"
        )
        common_dbs <- unique(sub("\\.(phr|nhr)$", "", common_dbs))
      }

      cli::cli_abort(c(
        "Cannot find BLAST database files at {.path {db_path}}.",
        i = "No .phr/.pin/.psq (protein) or .nhr/.nin/.nsq (nucleotide) files found for database prefix {.val {db_base}}.",
        if (length(common_dbs) > 0)
          c(i = "Common location: {.path {common_location}}",
            i = "Available databases: {.val {head(common_dbs, 5)}}{if (length(common_dbs) > 5) ', ...'}"),
        i = "Set {.env BLASTDB} environment variable or provide absolute path via {.arg dbdir}.",
        i = "Or set {.code validate_db = FALSE} to skip this check (not recommended)."
      ))
    }

    cli::cli_alert_success(
      "Found {length(db_files)} database file{?s} for {.val {db_base}} in {.path {db_dir}}"
    )
  }

  # --- Compose outfmt 6 (must be quoted if a shell is involved) ---
  fmt <- paste("6", paste(fields, collapse = " "))
  fmt_q <- shQuote(fmt, type = "sh")

  # --- Build BLASTP args ---
  # Quote any path-like arguments that may contain spaces if we fall back to system2.
  query_q <- shQuote(query_faa, type = "sh")
  db_q    <- shQuote(db_path,   type = "sh")

  args <- c(
    "-query", query_faa,            # raw for processx
    "-db",    db_path,              # raw for processx
    "-evalue", as.character(evalue),
    "-seg",   "yes",
    "-comp_based_stats", "2",
    "-num_threads", as.integer(threads),
    "-outfmt", fmt
  )
  args_sh <- c(                      # quoted for system2 with stderr/stdout capture
    "-query", query_q,
    "-db",    db_q,
    "-evalue", as.character(evalue),
    "-seg",   "yes",
    "-comp_based_stats", "2",
    "-num_threads", as.integer(threads),
    "-outfmt", fmt_q
  )
  if (length(more_args)) {
    # preserve extra args; quote any that have spaces to be safe under a shell
    more_args_sh <- ifelse(grepl("\\s", more_args), shQuote(more_args, type = "sh"), more_args)
    args    <- c(args,    more_args)
    args_sh <- c(args_sh, more_args_sh)
  }

  cli::cli_alert_info("Running: {blastp} -db {db_path} -query {query_faa} -outfmt {fmt}")

  # --- Hermetically set BLASTDB just for the call ---
  old_BLASTDB <- Sys.getenv("BLASTDB", unset = NA_character_)
  on.exit({
    if (!is.na(old_BLASTDB)) Sys.setenv(BLASTDB = old_BLASTDB) else Sys.unsetenv("BLASTDB")
  }, add = TRUE)
  if (nzchar(dbdir)) Sys.setenv(BLASTDB = dbdir)

  # --- Execute and capture (prefer processx if available to avoid shell quoting entirely) ---
  use_px <- requireNamespace("processx", quietly = TRUE)
  if (isTRUE(use_px)) {
    res <- processx::run(blastp, args = args, error_on_status = FALSE, echo_cmd = FALSE)
    status <- res$status
    out    <- res$stdout
    err    <- res$stderr
  } else {
    out <- suppressWarnings(system2(command = blastp, args = args_sh, stdout = TRUE, stderr = TRUE))
    status <- attr(out, "status")
    err    <- attr(out, "stderr")
  }

  if (!is.null(status) && status != 0) {
    msg <- if (!is.null(err) && length(err)) err else out
    cli::cli_abort(c(
      "blastp exited with status {status}.",
      ">" = paste(msg, collapse = "\n"),
      i = "Check binary path, DB prefix (quoted if it contains spaces), permissions, and device mounts."
    ))
  }

  # --- Parse and rank ---
  df <- readr::read_tsv(I(out), col_names = fields, show_col_types = FALSE)
  df |>
    dplyr::mutate(dplyr::across(c(pident, qcovs, bitscore, evalue), suppressWarnings(as.numeric))) |>
    dplyr::arrange(evalue, dplyr::desc(bitscore), dplyr::desc(qcovs))
}


#' Run BLASTP on one or more query FAA files and return tidy hits
#'
#' Inherits binary/db controls from blastp_capture; no reticulate involved.
#'
#' @inheritParams blastp_capture
#' @param faa_path Character vector of one or more query FASTA/FAA file paths
#' @return tibble of combined BLAST hits from all queries
#' @export
blastp_roi <- function(
    faa_path,
    db = "swissprot",
    dbdir = Sys.getenv("BLASTDB", ""),
    blastp_bin = Sys.getenv("BLASTP_BIN", ""),
    blastp_dir = Sys.getenv("BLASTP_PATH", ""),
    evalue = 1e-10,
    threads = parallel::detectCores(),
    fields = c("qseqid","sacc","stitle","pident","length","qcovs","bitscore","evalue","staxids"),
    validate_db = TRUE,
    more_args = character()
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("cli", quietly = TRUE)) stop("Package 'cli' is required.")

  # Allow vector input; run sequentially unless you later want parallelization
  faa_path <- unique(faa_path)
  res_list <- vector("list", length(faa_path))

  for (i in seq_along(faa_path)) {
    q <- faa_path[i]
    cli::cli_alert_info("BLASTP ROI: {i}/{length(faa_path)} - {.path {q}}")

    res_list[[i]] <- blastp_capture(
      query_faa   = q,
      db          = db,
      dbdir       = dbdir,
      blastp_bin  = blastp_bin,
      blastp_dir  = blastp_dir,
      evalue      = evalue,
      threads     = threads,
      fields      = fields,
      validate_db = validate_db,
      more_args   = more_args
    )
  }

  # Bind rows; ensure required columns exist
  out <- dplyr::bind_rows(res_list)
  needed <- c("qseqid","sacc","stitle","pident","length","qcovs","bitscore","evalue")
  missing <- setdiff(needed, colnames(out))
  if (length(missing)) {
    stop("Missing columns in BLAST output: ", paste(missing, collapse = ", "))
  }

  out
}


#' Reduce BLAST hits with cov/identity thresholds; select best hit or top-N per query
#'
#' @param hits Data frame of BLAST hits (from blastp_capture or blastp_roi)
#' @param min_qcov Minimum query coverage percentage (default 40)
#' @param min_pident Minimum percent identity (default 25)
#' @param besthit Logical; if TRUE, return only best hit per query (default TRUE)
#' @param max_per_query If not NULL and besthit=FALSE, return top N hits per query
#' @return Filtered tibble of BLAST hits
#' @export
reduce_hits <- function(
    hits,
    min_qcov = 40,
    min_pident = 25,
    besthit = TRUE,
    max_per_query = NULL
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!is.data.frame(hits) || nrow(hits) == 0) return(hits)

  h <- hits

  # Coerce numerics once, defensively; preserve NA semantics
  num_cols <- intersect(c("pident","qcovs","bitscore","evalue"), colnames(h))
  h <- dplyr::mutate(h, dplyr::across(dplyr::all_of(num_cols), suppressWarnings(as.numeric)))

  # Thresholds: allow NA to pass (as you did), which is reasonable for partial outputs
  if (!is.null(min_qcov) && "qcovs" %in% colnames(h)) {
    h <- dplyr::filter(h, is.na(qcovs) | qcovs >= min_qcov)
  }
  if (!is.null(min_pident) && "pident" %in% colnames(h)) {
    h <- dplyr::filter(h, is.na(pident) | pident >= min_pident)
  }

  h <- dplyr::group_by(h, qseqid)

  if (isTRUE(besthit)) {
    # Three-step tie-break: evalue -> bitscore -> qcovs
    h <- dplyr::slice_min(h, order_by = evalue, n = 1, with_ties = TRUE)
    h <- dplyr::slice_max(h, order_by = bitscore, n = 1, with_ties = TRUE)
    h <- dplyr::slice_max(h, order_by = qcovs, n = 1, with_ties = FALSE)
  } else {
    h <- dplyr::arrange(h, evalue, dplyr::desc(bitscore), dplyr::desc(qcovs))
    if (!is.null(max_per_query) && is.finite(max_per_query)) {
      h <- dplyr::slice_head(h, n = max_per_query)
    }
  }

  dplyr::ungroup(h)
}
