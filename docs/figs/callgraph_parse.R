# docs/figs/callgraph_parse.R -------------------------------------------------
# Exact call graph of the R sources, from R's own parser.
#
# Parses every file under R/, run/ and tests/ plus _targets.R with
# parse(keep.source = TRUE) and reads getParseData(). No regular expressions
# touch the code. Writes two tables under docs/figs/callgraph/:
#
#   definitions.csv  one row per function definition: name, file, the line
#                    span of the function expression, the enclosing function
#                    when the definition is nested, and whether it is top level.
#   edges.csv        one row per reference from a caller to a callee. The
#                    caller is the innermost function whose span contains the
#                    token, or `<script:FILE>` for top-level code. `kind` says
#                    how the callee was referenced:
#                      call      name(...)                 SYMBOL_FUNCTION_CALL
#                      ref       bare name as an argument  SYMBOL matching a definition
#                      string    "name" as a string        STR_CONST matching a definition
#                      external  pkg::name(...)            SYMBOL_PACKAGE + call
#                      target    a tar_target body (or script statement) reads
#                                another target's value (callee is its name)
#                      target_write  a script statement assigns a script
#                                value inside itself (`if (...) x <- f()`),
#                                so it is a producer of x
#                    Only names that are defined somewhere in the parsed files
#                    become `call`, `ref` or `string` edges; base R and package
#                    calls without a namespace prefix are not recorded.
#                    `target` names the tar_target() whose body holds the
#                    token, for tokens inside _targets.R; empty elsewhere.
#   targets.csv      one row per tar_target() in _targets.R and per top-level
#                    statement of every other script (named by the assigned
#                    symbol, or <expr:FILE:LINE>): file and line span. A
#                    script's statements are targets in all but name: a value
#                    assigned by one is read by the next.
#
# The target column and the target edges are what lets a chart's arrows be
# checked: the pipeline's data flow between two functions that never call each
# other runs through target values, and those dependencies are here.
#
# Reads source text only. Nothing under data/ is opened (hard rule 1).
#
#   Rscript docs/figs/callgraph_parse.R
# ------------------------------------------------------------------------------

# Run from the project root (the directory holding _targets.R).
if (!file.exists("_targets.R")) stop("run from the project root: Rscript docs/figs/callgraph_parse.R")

files <- c(list.files("R", "\\.R$", full.names = TRUE),
           list.files("run", "\\.R$", full.names = TRUE),
           list.files("tests", "\\.R$", full.names = TRUE),
           "_targets.R")
files <- files[file.exists(files)]

parsed <- lapply(files, function(f) {
  p <- tryCatch(parse(f, keep.source = TRUE), error = function(e) NULL)
  if (is.null(p)) { message("parse error, skipped: ", f); return(NULL) }
  pd <- getParseData(p, includeText = TRUE)
  pd$file <- f
  pd
})
names(parsed) <- files
parsed <- Filter(Negate(is.null), parsed)

# --- definitions ---------------------------------------------------------------
# A definition is  SYMBOL <- FUNCTION ...  (or =, or <<-): the FUNCTION token's
# parent is the function expr; that expr's parent is the assignment expr whose
# first child is an expr wrapping a SYMBOL and whose second child is an
# assignment operator.
def_rows <- list()
for (f in names(parsed)) {
  pd <- parsed[[f]]
  fun_tok <- pd[pd$token == "FUNCTION", , drop = FALSE]
  for (i in seq_len(nrow(fun_tok))) {
    fun_expr <- pd[pd$id == fun_tok$parent[i], , drop = FALSE]
    assign_expr <- pd[pd$id == fun_expr$parent, , drop = FALSE]
    nm <- NA_character_
    if (nrow(assign_expr) == 1L) {
      kids <- pd[pd$parent == assign_expr$id, , drop = FALSE]
      kids <- kids[order(kids$line1, kids$col1), , drop = FALSE]
      if (nrow(kids) >= 3L && kids$token[2] %in% c("LEFT_ASSIGN", "EQ_ASSIGN") &&
          kids$id[3] == fun_expr$id) {
        lhs <- pd[pd$parent == kids$id[1] & pd$token == "SYMBOL", , drop = FALSE]
        if (nrow(lhs) == 1L) nm <- lhs$text
      }
    }
    def_rows[[length(def_rows) + 1L]] <- data.frame(
      name = nm, file = f,
      line_start = fun_expr$line1, col_start = fun_expr$col1,
      line_end = fun_expr$line2, col_end = fun_expr$col2,
      fun_id = fun_expr$id, stringsAsFactors = FALSE)
  }
}
defs <- do.call(rbind, def_rows)
defs$anonymous <- is.na(defs$name)
defs$name[defs$anonymous] <- sprintf("<anon:%s:%d>", basename(defs$file[defs$anonymous]),
                                     defs$line_start[defs$anonymous])

# Enclosing definition: the smallest other definition in the same file whose
# span strictly contains this one.
.contains <- function(a, b) {  # does span a contain span b
  (a$line_start < b$line_start | (a$line_start == b$line_start & a$col_start <= b$col_start)) &
  (a$line_end   > b$line_end   | (a$line_end   == b$line_end   & a$col_end   >= b$col_end))
}
defs$enclosing <- ""
for (i in seq_len(nrow(defs))) {
  same <- defs[defs$file == defs$file[i] & seq_len(nrow(defs)) != i, , drop = FALSE]
  if (!nrow(same)) next
  inside <- same[.contains(same, defs[i, ]), , drop = FALSE]
  if (nrow(inside)) {
    span <- (inside$line_end - inside$line_start)
    defs$enclosing[i] <- inside$name[which.min(span)]
  }
}
defs$top_level <- defs$enclosing == ""
named_defs <- unique(defs$name[!defs$anonymous])
# Names visible everywhere: top-level definitions under R/, which the runners
# source as a whole. A bare-name or string reference is resolved only against
# these plus the top-level definitions of the file it sits in, so that a local
# helper called `p` in one test script does not turn every `p` into an edge.
global_defs <- unique(defs$name[!defs$anonymous & defs$top_level & startsWith(defs$file, "R/")])
visible_in <- function(f) unique(c(global_defs,
  defs$name[!defs$anonymous & defs$top_level & defs$file == f]))

# --- edges ---------------------------------------------------------------------
# Innermost enclosing definition of a token position within a file.
# Anonymous closures (`function(i) ...` inside lapply and the like) are
# transparent: a call inside one belongs to the named function around it.
owner_of <- function(f, line, col) {
  d <- defs[defs$file == f & !defs$anonymous, , drop = FALSE]
  if (!nrow(d)) return(sprintf("<script:%s>", f))
  hit <- d[(d$line_start < line | (d$line_start == line & d$col_start <= col)) &
           (d$line_end   > line | (d$line_end   == line & d$col_end   >= col)), , drop = FALSE]
  if (!nrow(hit)) return(sprintf("<script:%s>", f))
  hit$name[which.min(hit$line_end - hit$line_start)]
}

# --- targets ---------------------------------------------------------------------
# Every tar_target(NAME, BODY, ...) call: its name and the span of the call
# expression, so a token can be assigned to the target whose body holds it.
tgt_rows <- list()
for (f in names(parsed)) {
  pd <- parsed[[f]]
  tt <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "tar_target", , drop = FALSE]
  for (i in seq_len(nrow(tt))) {
    call_expr <- pd[pd$id == pd$parent[pd$id == tt$parent[i]], , drop = FALSE]
    kids <- pd[pd$parent == call_expr$id & pd$token == "expr", , drop = FALSE]
    kids <- kids[order(kids$line1, kids$col1), , drop = FALSE]
    if (nrow(kids) < 2L) next
    nm <- pd[pd$parent == kids$id[2] & pd$token == "SYMBOL", , drop = FALSE]
    if (nrow(nm) != 1L) next
    tgt_rows[[length(tgt_rows) + 1L]] <- data.frame(
      name = nm$text, file = f, line_start = call_expr$line1, col_start = call_expr$col1,
      line_end = call_expr$line2, col_end = call_expr$col2, name_id = nm$id,
      stringsAsFactors = FALSE)
  }
}
# Top-level statements of the run and test scripts are the same kind of node:
# `x <- f(...)` defines a script value that later statements read. Each
# top-level expression that is not a function definition becomes a pseudo-
# target named by the assigned symbol, or `<expr:FILE:LINE>` when nothing is
# assigned. _targets.R is excluded here because its statements are the
# tar_target() calls above.
for (f in setdiff(names(parsed), "_targets.R")) {
  pd <- parsed[[f]]
  tops <- pd[pd$parent == 0 & pd$token == "expr", , drop = FALSE]
  fn_names <- defs$name[defs$file == f & !defs$anonymous & defs$top_level]
  for (i in seq_len(nrow(tops))) {
    kids <- pd[pd$parent == tops$id[i], , drop = FALSE]
    kids <- kids[order(kids$line1, kids$col1), , drop = FALSE]
    nm <- sprintf("<expr:%s:%d>", f, tops$line1[i]); name_id <- NA_integer_
    if (nrow(kids) >= 3L && kids$token[2] %in% c("LEFT_ASSIGN", "EQ_ASSIGN")) {
      lhs <- pd[pd$parent == kids$id[1] & pd$token == "SYMBOL", , drop = FALSE]
      if (nrow(lhs) == 1L) {
        if (lhs$text %in% fn_names) next
        nm <- lhs$text; name_id <- lhs$id
      }
    }
    tgt_rows[[length(tgt_rows) + 1L]] <- data.frame(
      name = nm, file = f, line_start = tops$line1[i], col_start = tops$col1[i],
      line_end = tops$line2[i], col_end = tops$col2[i], name_id = name_id,
      stringsAsFactors = FALSE)
  }
}
targets <- if (length(tgt_rows)) do.call(rbind, tgt_rows) else
  data.frame(name = character(0), file = character(0), line_start = integer(0),
             col_start = integer(0), line_end = integer(0), col_end = integer(0),
             name_id = integer(0))
target_of <- function(f, line, col) {
  t <- targets[targets$file == f, , drop = FALSE]
  if (!nrow(t)) return("")
  hit <- t[(t$line_start < line | (t$line_start == line & t$col_start <= col)) &
           (t$line_end   > line | (t$line_end   == line & t$col_end   >= col)), , drop = FALSE]
  if (!nrow(hit)) return("")
  hit$name[which.min(hit$line_end - hit$line_start)]
}

edge_rows <- list()
for (f in names(parsed)) {
  pd <- parsed[[f]]
  # target-to-target dependencies: a bare SYMBOL inside one target's body that
  # names another target (never the target's own name token)
  tnames <- targets$name[targets$file == f]
  if (length(tnames)) {
    # parse ids are unique within a file only, so the name tokens to skip are
    # this file's
    own_name_ids <- targets$name_id[targets$file == f & !is.na(targets$name_id)]
    syms <- pd[pd$token == "SYMBOL" & pd$text %in% tnames & !pd$id %in% own_name_ids, , drop = FALSE]
    for (i in seq_len(nrow(syms))) {
      own <- target_of(f, syms$line1[i], syms$col1[i])
      if (own == "" || own == syms$text[i]) next
      # `x <- ...`, `x[[k]] <- ...` or `x$a <- ...` inside a statement (an
      # if-branch or a loop, say) writes x rather than reading it: the
      # statement is then a producer of x. Walk up through indexing while the
      # symbol stays the leftmost operand.
      writes <- FALSE
      node <- syms$parent[i]
      repeat {
        par <- pd[pd$id == pd$parent[pd$id == node], , drop = FALSE]
        if (!nrow(par)) break
        kids <- pd[pd$parent == par$id, , drop = FALSE]
        kids <- kids[order(kids$line1, kids$col1), , drop = FALSE]
        if (nrow(kids) < 2L || kids$id[1] != node) break
        if (kids$token[2] %in% c("LEFT_ASSIGN", "EQ_ASSIGN")) { writes <- TRUE; break }
        if (kids$token[2] %in% c("'['", "LBB", "'$'", "'@'")) { node <- par$id; next }
        break
      }
      edge_rows[[length(edge_rows) + 1L]] <- data.frame(
        caller = sprintf("<script:%s>", f), caller_file = f, callee = syms$text[i],
        kind = if (writes) "target_write" else "target", line = syms$line1[i],
        target = own, stringsAsFactors = FALSE)
    }
  }
  # name tokens of definitions in this file (the LHS symbols), to exclude from refs
  lhs_ids <- integer(0)
  dd <- defs[defs$file == f & !defs$anonymous, , drop = FALSE]
  for (i in seq_len(nrow(dd))) {
    fun_expr <- pd[pd$id == dd$fun_id[i], , drop = FALSE]
    assign_expr <- pd[pd$id == fun_expr$parent, , drop = FALSE]
    kids <- pd[pd$parent == assign_expr$id, , drop = FALSE]
    kids <- kids[order(kids$line1, kids$col1), , drop = FALSE]
    lhs <- pd[pd$parent == kids$id[1] & pd$token == "SYMBOL", , drop = FALSE]
    lhs_ids <- c(lhs_ids, lhs$id)
  }
  calls <- pd[pd$token == "SYMBOL_FUNCTION_CALL", , drop = FALSE]
  for (i in seq_len(nrow(calls))) {
    sib <- pd[pd$parent == calls$parent[i], , drop = FALSE]
    pkg <- sib$text[sib$token == "SYMBOL_PACKAGE"]
    callee <- calls$text[i]
    kind <- "call"
    if (length(pkg)) { callee <- paste0(pkg, "::", callee); kind <- "external" }
    else if (!callee %in% named_defs) next
    edge_rows[[length(edge_rows) + 1L]] <- data.frame(
      caller = owner_of(f, calls$line1[i], calls$col1[i]), caller_file = f,
      callee = callee, kind = kind, line = calls$line1[i],
      target = target_of(f, calls$line1[i], calls$col1[i]), stringsAsFactors = FALSE)
  }
  vis <- visible_in(f)
  refs <- pd[pd$token == "SYMBOL" & pd$text %in% vis & !pd$id %in% lhs_ids, , drop = FALSE]
  for (i in seq_len(nrow(refs))) {
    edge_rows[[length(edge_rows) + 1L]] <- data.frame(
      caller = owner_of(f, refs$line1[i], refs$col1[i]), caller_file = f,
      callee = refs$text[i], kind = "ref", line = refs$line1[i],
      target = target_of(f, refs$line1[i], refs$col1[i]), stringsAsFactors = FALSE)
  }
  strs <- pd[pd$token == "STR_CONST", , drop = FALSE]
  if (nrow(strs)) {
    strs$bare <- gsub("^[\"']|[\"']$", "", strs$text)
    strs <- strs[strs$bare %in% vis, , drop = FALSE]
    for (i in seq_len(nrow(strs))) {
      edge_rows[[length(edge_rows) + 1L]] <- data.frame(
        caller = owner_of(f, strs$line1[i], strs$col1[i]), caller_file = f,
        callee = strs$bare[i], kind = "string", line = strs$line1[i],
        target = target_of(f, strs$line1[i], strs$col1[i]), stringsAsFactors = FALSE)
    }
  }
}
edges <- do.call(rbind, edge_rows)
edges <- edges[edges$caller != edges$callee, , drop = FALSE]
edges <- unique(edges)

dir.create("docs/figs/callgraph", showWarnings = FALSE, recursive = TRUE)
write.csv(defs[, c("name", "file", "line_start", "line_end", "enclosing", "top_level", "anonymous")],
          "docs/figs/callgraph/definitions.csv", row.names = FALSE)
write.csv(edges, "docs/figs/callgraph/edges.csv", row.names = FALSE)
write.csv(targets[, c("name", "file", "line_start", "line_end")],
          "docs/figs/callgraph/targets.csv", row.names = FALSE)
cat(sprintf("parsed %d files; %d function definitions (%d named, %d nested); %d targets; %d edges (%s)\n",
            length(parsed), nrow(defs), sum(!defs$anonymous), sum(!defs$top_level), nrow(targets), nrow(edges),
            paste(sprintf("%s %d", names(table(edges$kind)), table(edges$kind)), collapse = ", ")))
