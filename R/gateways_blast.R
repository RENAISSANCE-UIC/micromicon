# Suppress R CMD check NOTEs for dplyr NSE variables used in gateway
utils::globalVariables(c("pident", "qcovs", "bitscore", "evalue", "qseqid"))

#' Create BLAST Gateway
#'
#' @description
#' Creates a gateway for running BLAST searches against local databases.
#' This gateway provides methods for protein BLAST searches (BLASTP).
#'
#' Designed for extensibility - future versions may support BLASTN, BLASTX, etc.
#' by adding corresponding methods to the gateway.
#'
#' Gateway Contract (for use cases):
#' - blastp_capture(query_faa, ...): Run BLASTP on single query file
#' - blastp_roi(faa_path, ...): Run BLASTP on multiple query files
#' - reduce_hits(hits, ...): Filter and rank BLAST results
#'
#' @param blast_type Type of BLAST to configure (default "blastp"). Reserved
#'   for future extensibility (blastn, blastx, etc.)
#'
#' @return Gateway object (list with BLAST methods)
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' gateway <- create_blast_gateway()
#' results <- gateway$blastp_capture(
#'   query_faa = "proteins.faa",
#'   db = "swissprot"
#' )
#' filtered <- gateway$reduce_hits(results, min_qcov = 80)
#' }
create_blast_gateway <- function(blast_type = "blastp") {

  # ---- Internal helper functions (shared across BLAST types) ----

  # Resolve BLAST binary path deterministically
  resolve_blast_binary <- function(bin, dir, program = "blastp") {
    if (nzchar(bin)) {
      if (file.exists(bin) && isTRUE(file.info(bin)$isdir)) {
        file.path(bin, program)
      } else {
        bin
      }
    } else if (nzchar(dir)) {
      if (file.exists(dir) && isTRUE(file.info(dir)$isdir)) {
        file.path(dir, program)
      } else {
        dir
      }
    } else {
      Sys.which(program)
    }
  }

  # Validate BLAST database exists
  validate_blast_db <- function(db_path, blastdbcmd) {
    if (nzchar(blastdbcmd)) {
      # Use blastdbcmd for validation
      db_path_q <- shQuote(db_path, type = "sh")
      info <- tryCatch(
        system2(blastdbcmd, args = c("-db", db_path_q, "-info"),
                stdout = TRUE, stderr = TRUE),
        error = function(e) character()
      )

      if (length(info) == 0 || any(grepl("\\bError\\b", info))) {
        cli::cli_abort(c(
          "Could not open BLAST DB at {.path {db_path}}.",
          ">" = paste(info, collapse = "\n"),
          i = "Provide an absolute DB prefix or set {.env BLASTDB} correctly."
        ))
      }
    } else {
      # Fallback: manual file check
      cli::cli_warn("blastdbcmd not found; attempting manual database file check.")

      db_dir <- dirname(db_path)
      if (db_dir == ".") db_dir <- getwd()
      db_base <- basename(db_path)

      # Check for protein (.phr, .pin, .psq) or nucleotide (.nhr, .nin, .nsq) files
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
  }

  # ---- Gateway Methods ----

  list(
    # Run BLASTP with explicit binary and DB paths
    #
    # @param query_faa Path to query FASTA/FAA file
    # @param db BLAST database prefix (e.g., "swissprot" or absolute "/path/to/swissprot")
    # @param dbdir Directory containing BLAST DBs; if db is relative, dbdir is prepended
    # @param blastp_bin Full path to blastp binary (preferred explicit control)
    # @param blastp_dir Directory containing blastp; used if blastp_bin is NULL
    # @param evalue E-value threshold
    # @param threads Integer number of threads
    # @param fields Vector of outfmt 6 columns
    # @param validate_db Logical; check DB via blastdbcmd before running
    # @param more_args Character vector of additional BLASTP args
    # @return tibble of hits; aborts with rich diagnostics on failure
    blastp_capture = function(
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
      # Dependencies
      if (!requireNamespace("cli", quietly = TRUE)) {
        cli::cli_abort("Package 'cli' is required")
      }
      if (!requireNamespace("readr", quietly = TRUE)) {
        cli::cli_abort("Package 'readr' is required")
      }
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        cli::cli_abort("Package 'dplyr' is required")
      }

      # Resolve blastp binary
      candidate_bin <- resolve_blast_binary(blastp_bin, blastp_dir, "blastp")

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

      # Preflight: check query file
      if (!file.exists(query_faa)) {
        cli::cli_abort("Query file not found: {.path {query_faa}}")
      }

      # Derive db_path
      db_path <- if (nzchar(dbdir) && !grepl("^/", db)) {
        file.path(dbdir, db)
      } else {
        db
      }

      # Warn if using relative path without explicit dbdir
      if (!nzchar(dbdir) && !grepl("^/", db)) {
        cli::cli_warn(c(
          "Using relative database path {.val {db}} without explicit {.arg dbdir}.",
          i = "Set {.env BLASTDB} or provide {.arg dbdir} parameter for reliable database location.",
          i = "Current working directory: {.path {getwd()}}"
        ))
      }

      # Validate DB if requested
      if (validate_db) {
        validate_blast_db(db_path, blastdbcmd)
      }

      # Compose outfmt 6
      fmt <- paste("6", paste(fields, collapse = " "))
      fmt_q <- shQuote(fmt, type = "sh")

      # Build BLASTP args
      query_q <- shQuote(query_faa, type = "sh")
      db_q    <- shQuote(db_path,   type = "sh")

      args <- c(
        "-query", query_faa,
        "-db",    db_path,
        "-evalue", as.character(evalue),
        "-seg",   "yes",
        "-comp_based_stats", "2",
        "-num_threads", as.integer(threads),
        "-outfmt", fmt
      )

      args_sh <- c(
        "-query", query_q,
        "-db",    db_q,
        "-evalue", as.character(evalue),
        "-seg",   "yes",
        "-comp_based_stats", "2",
        "-num_threads", as.integer(threads),
        "-outfmt", fmt_q
      )

      if (length(more_args)) {
        more_args_sh <- ifelse(grepl("\\s", more_args),
                               shQuote(more_args, type = "sh"),
                               more_args)
        args    <- c(args,    more_args)
        args_sh <- c(args_sh, more_args_sh)
      }

      cli::cli_alert_info("Running: {blastp} -db {db_path} -query {query_faa} -outfmt {fmt}")

      # Hermetically set BLASTDB just for the call
      old_BLASTDB <- Sys.getenv("BLASTDB", unset = NA_character_)
      on.exit({
        if (!is.na(old_BLASTDB)) {
          Sys.setenv(BLASTDB = old_BLASTDB)
        } else {
          Sys.unsetenv("BLASTDB")
        }
      }, add = TRUE)
      if (nzchar(dbdir)) Sys.setenv(BLASTDB = dbdir)

      # Execute and capture (prefer processx to avoid shell quoting)
      use_px <- requireNamespace("processx", quietly = TRUE)
      if (isTRUE(use_px)) {
        res <- processx::run(blastp, args = args, error_on_status = FALSE, echo_cmd = FALSE)
        status <- res$status
        out    <- res$stdout
        err    <- res$stderr
      } else {
        out <- suppressWarnings(system2(command = blastp, args = args_sh,
                                        stdout = TRUE, stderr = TRUE))
        status <- attr(out, "status")
        err    <- attr(out, "stderr")
      }

      if (!is.null(status) && status != 0) {
        msg <- if (!is.null(err) && length(err)) err else out
        cli::cli_abort(c(
          "blastp exited with status {status}.",
          ">" = paste(msg, collapse = "\n"),
          i = "Check binary path, DB prefix, permissions, and device mounts."
        ))
      }

      # Parse and rank
      df <- readr::read_tsv(I(out), col_names = fields, show_col_types = FALSE)
      df |>
        dplyr::mutate(dplyr::across(c(pident, qcovs, bitscore, evalue),
                                     suppressWarnings(as.numeric))) |>
        dplyr::arrange(evalue, dplyr::desc(bitscore), dplyr::desc(qcovs))
    },

    # Run BLASTP on one or more query FAA files and return tidy hits
    #
    # @inheritParams blastp_capture
    # @param faa_path Character vector of one or more query FASTA/FAA file paths
    # @return tibble of combined BLAST hits from all queries
    blastp_roi = function(
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
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        cli::cli_abort("Package 'dplyr' is required")
      }
      if (!requireNamespace("cli", quietly = TRUE)) {
        cli::cli_abort("Package 'cli' is required")
      }

      # Allow vector input; run sequentially
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
        cli::cli_abort("Missing columns in BLAST output: {paste(missing, collapse = ', ')}")
      }

      out
    },

    # Reduce BLAST hits with coverage/identity thresholds
    #
    # @param hits Data frame of BLAST hits (from blastp_capture or blastp_roi)
    # @param min_qcov Minimum query coverage percentage (default 40)
    # @param min_pident Minimum percent identity (default 25)
    # @param besthit Logical; if TRUE, return only best hit per query (default TRUE)
    # @param max_per_query If not NULL and besthit=FALSE, return top N hits per query
    # @return Filtered tibble of BLAST hits
    reduce_hits = function(
        hits,
        min_qcov = 40,
        min_pident = 25,
        besthit = TRUE,
        max_per_query = NULL
    ) {
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        cli::cli_abort("Package 'dplyr' is required")
      }

      if (!is.data.frame(hits) || nrow(hits) == 0) return(hits)

      h <- hits

      # Coerce numerics defensively
      num_cols <- intersect(c("pident","qcovs","bitscore","evalue"), colnames(h))
      h <- dplyr::mutate(h, dplyr::across(dplyr::all_of(num_cols),
                                          suppressWarnings(as.numeric)))

      # Apply thresholds
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
  )
}
