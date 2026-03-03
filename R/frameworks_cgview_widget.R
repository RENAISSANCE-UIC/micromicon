# cgview_widget.R
# htmlwidgets binding for CGView.js circular genome viewer.

#' Render a circular genome map with CGView.js
#'
#' @description
#' Creates an interactive CGView.js htmlwidget from a CGView JSON payload (as
#' produced by \code{build_cgview_json()} or \code{beta_cgview_from_entity()}).
#'
#' @param json A CGView JSON structure. Accepts either:
#'   - An R list matching the CGView JSON schema (serialised internally via
#'     [jsonlite::toJSON()])
#'   - A character string of already-serialised JSON
#' @param width,height Widget dimensions. `NULL` lets the browser decide.
#'   `height` defaults to `700` px.
#' @param elementId Optional fixed HTML id for the container div.
#'
#' @return An htmlwidget.
#'
#' @seealso \code{build_cgview_json()}, \code{beta_cgview_from_entity()}
#'
#' @keywords internal
#' @noRd
beta_cgview <- function(json, width = NULL, height = 2000, elementId = NULL) {

  # htmlwidgets calls jsonlite::toJSON() internally on the x list.
  # Hand it a plain R list, never a pre-serialised json-class string.
  if (is.character(json)) {
    json_list <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  } else if (is.list(json)) {
    json_list <- json
  } else {
    stop("`json` must be a character string or a list.", call. = FALSE)
  }

  x <- list(cgview_json = json_list)

  htmlwidgets::createWidget(
    name      = 'cgview',
    x         = x,
    width     = width,
    height    = height,
    package   = 'micromicon',
    elementId = elementId
  )
}

#' Side-by-side pair of CGView circular genome maps
#'
#' @description
#' Renders two \code{beta_cgview()} widgets in a flex row. Typical use: left panel
#' for reference features, right panel for variants.
#'
#' @param left_json,right_json CGView JSON payloads (R lists or character
#'   strings) as produced by \code{build_cgview_json()}.
#' @param height Widget height in pixels (applied to both panels).
#' @param left_title,right_title Optional character title rendered above each
#'   panel. `NULL` suppresses the title.
#' @param gap CSS gap between the two panels (default `"12px"`).
#'
#' @return An `htmltools::browsable` tagList.
#'
#' @seealso \code{build_cgview_json()}, \code{beta_cgview()}
#'
#' @keywords internal
#' @noRd
beta_cgview_pair <- function(left_json, right_json,
                              height      = 2000,
                              left_title  = NULL,
                              right_title = NULL,
                              gap         = "12px") {

  left_w  <- beta_cgview(left_json,  height = height)
  right_w <- beta_cgview(right_json, height = height)

  title_style <- paste0(
    "margin:0 0 4px 0;",
    "font:bold 13px/1.2 sans-serif;",
    "color:#444;",
    "text-align:center;"
  )

  make_panel <- function(widget, title) {
    kids <- if (!is.null(title))
      list(htmltools::tags$p(title, style = title_style), widget)
    else
      list(widget)
    htmltools::tags$div(style = "flex:1; min-width:0;", kids)
  }

  htmltools::browsable(
    htmltools::tags$div(
      style = paste0("display:flex; gap:", gap, "; align-items:flex-start;"),
      make_panel(left_w,  left_title),
      make_panel(right_w, right_title)
    )
  )
}

#' Shiny output for beta_cgview
#'
#' @param outputId Shiny output id.
#' @param width,height CSS dimensions (default `"100%"` / `"700px"`).
#'
#' @keywords internal
#' @noRd
beta_cgviewOutput <- function(outputId, width = "100%", height = "700px") {
  htmlwidgets::shinyWidgetOutput(
    outputId = outputId,
    name     = 'cgview',
    package  = 'micromicon',
    width    = width,
    height   = height
  )
}

#' Shiny render function for beta_cgview
#'
#' @param expr An expression that returns a \code{beta_cgview()} widget.
#' @param env,quoted Passed to [htmlwidgets::shinyRenderWidget()].
#'
#' @keywords internal
#' @noRd
render_beta_cgview <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) expr <- substitute(expr)
  htmlwidgets::shinyRenderWidget(expr, beta_cgviewOutput, env, quoted = TRUE)
}

#' @keywords internal
.cgview_open <- function(x, viewer) {
  viewer <- match.arg(viewer, c("pane", "browser"))
  if (identical(viewer, "browser")) {
    tmp <- tempfile(fileext = ".html")
    htmltools::save_html(x, tmp)
    utils::browseURL(tmp)
    return(invisible(x))
  }
  x
}
