r"""Reachability per arm of the pipeline, and the check of a chart against it.

Consumes the exact parse tables written by docs/figs/callgraph_parse.R
(definitions.csv, edges.csv) and, for each arm declared in ARMS below,

  * walks the call graph from the arm's entry scripts, expanding functions
    that live in the arm's own files and stopping at functions that live
    elsewhere (recorded as the arm's boundary: what it borrows);
  * writes docs/figs/callgraph/<arm>.md with the reachable functions, their
    callers and callees, the boundary, the external package calls, and the
    owned functions that were NOT reached (the excluded set, listed by name);
  * lists, in the same file, every tar_target whose body calls an owned
    function: what it calls and which other targets it reads. That is the
    data flow between functions that never call each other;
  * if the arm has a chart, checks that every reachable owned function is
    named on it or is in the arm's allow-list with a reason, and that every
    arrow on the chart is backed by an edge in the tables (see check_arrows),
    and exits 1 otherwise.

Arrows. A chart arrow from box u to box v is accepted when the tables hold
one of: a function of v calls a function of u (v consumes u); a function of u
calls a function of v (u drives v); a target of v reads a target of u; a
target of v calls a function of u; a function of v is called by a target that
reads a target of u; the targets calling the two functions are linked by a
read; one target or function body calls a function of each box; or a function
of v is called by a function that a target reading a target of u calls (the
bundle a runner reads, applied through the engine's apply function). A box's functions are the names on its function line (the third
argument of \B / \Bs); an object box (\Obj) names targets. Arrows drawn with
the `ex` style cross a process boundary (a bundle written by one script and
read by another) and cannot be verified from source: they are counted and
listed, not checked. Arrows to note boxes (`en`) are ignored. Pairs of boxes
whose functions call each other but have no arrow are reported as warnings.

    python docs/figs/callgraph_check.py            # parse, then check every arm
    python docs/figs/callgraph_check.py --no-parse # reuse the last parse
    python docs/figs/callgraph_check.py engine     # one arm

Reads source-derived tables only; nothing under data/ is opened.
"""
import csv, os, re, subprocess, sys
from collections import defaultdict
from fnmatch import fnmatch

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
CG = os.path.join(ROOT, "docs", "figs", "callgraph")

ENGINE_FILES = ["R/00_utils.R", "R/01_load.R", "R/02_validate.R", "R/03_folds.R",
                "R/04_features.R", "R/04b_conditional.R", "R/04c_prior_diagnostics.R",
                "R/05_formula.R", "R/06_layer1.R", "R/07_layer2.R", "R/08_diagnostics.R",
                "R/09_metrics.R", "R/10_bundle.R", "R/11_run.R"]

ARMS = {
    "engine": dict(
        # tests/gam_qc.R is the diagnostics runner of record: the permutation
        # null and the final-model diagnostics table come from it, not from
        # the targets graph.
        entries=["_targets.R", "run/internal.R", "run/test_look.R", "run/external.R",
                 "tests/gam_qc.R", "tests/concurvity_null.R"],
        owned=ENGINE_FILES + ["tests/gam_qc.R", "tests/gam_qc_common.R", "tests/concurvity_null.R"],
        chart="docs/figs/pipeline_implementation.tex",
        # WHY THERE IS AN ALLOW-LIST, AND WHY IT NAMES FUNCTIONS ONE BY ONE.
        #
        # The check asks that every top-level function reachable from the
        # arm's entry scripts be named on the chart. A chart that named all of
        # them would be unreadable: a handful of reachable functions are
        # numeric internals of one estimator (the proportional-odds and
        # variance-component code inside the ordinal delta fit), constructors
        # that return an empty frame for a design that does not use the
        # construct, or three-line helpers local to a QC script. Those are
        # left off the chart, and this list is where that decision lives.
        #
        # It is a list of names rather than a rule (such as "ignore names that
        # start with a dot") because a rule hides. The regex version of this
        # script used exactly that kind of rule, and it is how `.prior_jobs()`
        # and `.input_fingerprint()`, both essential to the fold logic and the
        # apply provenance, stayed off the first chart without anyone noticing.
        # Each exclusion here carries its reason, the count of exclusions is
        # printed beside the coverage on every run, and the check warns if an
        # allow-listed function later appears on the chart so the entry can be
        # removed. A reviewer can read the list and disagree with any line of
        # it; a rule would give them nothing to disagree with.
        #
        # What belongs here: functions whose omission changes nothing a reader
        # of the chart needs to know. What does not belong here: anything that
        # selects rows, resolves a role or fold, builds a formula, reads
        # config, or decides what enters the bundle. If such a function ever
        # shows up in the failure list, the answer is to draw it, not to add it.
        allow={
            ".resolve_resample_col": "inside resample_cols(): resolves one declared column name",
            ".empty_intervention_priors": "empty-frame constructor for a design with no lambda",
            ".empty_magnitude_priors": "empty-frame constructor for a design with no delta",
            ".delta_raw_resid": "numeric internals of delta_value() / fit_delta_ordinal()",
            ".fit_var_components": "variance-component estimator inside fit_delta() / fit_lambda_ln()",
            ".fit_var_components_rep": "replicate-statistics variant of the same",
            ".mid_pit_logit": "ordinal residual internals",
            ".par_levels": "ordinal parameter unpacking",
            ".par_theta": "ordinal parameter unpacking",
            ".polr_cum_below": "proportional-odds internals",
            ".polr_probs": "proportional-odds internals",
            ".theta_from_par": "proportional-odds internals",
            ".value_scale_resid": "ordinal residual internals",
            "layer1_nested_l": "serves the stacked booster only; drawn on the boosters chart",
            ".opt": "tests/gam_qc.R: command-line option reader",
            ".run_tab": "tests/gam_qc.R: reads one table of a run directory",
            "ll": "tests/gam_qc.R: one-line log helper",
            "win2": "tests/gam_qc.R: two-sided window helper for the support grid",
        },
        # utilities named on the chart only when they carry meaning
        utility_files=["R/00_utils.R"]),
    "boosters": dict(
        entries=["_targets.R", "run/test_look.R", "run/external.R"],
        owned=["R/09b_xgboost.R"],
        chart="docs/figs/pipeline_boosters.tex", allow={}, utility_files=[]),
    "severity": dict(
        entries=["_targets.R", "run/test_look.R", "run/external.R"],
        owned=["R/09c_apache.R", "R/09d_sofa.R", "R/10b_severity.R"],
        chart="docs/figs/pipeline_severity.tex", allow={}, utility_files=[]),
    "attribution": dict(
        # The scripts of record, confirmed from the run manifests (2026-09-25):
        # attr_replicates.R writes the attrgen store, attr_metrics.R reads it
        # (attrmetrics, and the attrgate check of the frozen gate runs),
        # attr_external.R writes attrext, attr_bands.R writes attrbands, and
        # run/patient_card.R writes card. run/attribution.R (R/13b floor) and
        # the coupattr / attrties / shapfloor scripts produced the frozen gate
        # inputs on 2026-09-05/06 and are not re-run; they are left out so the
        # walk shows what the current design reaches.
        entries=["tests/attr_replicates.R", "tests/attr_metrics.R", "tests/attr_external.R",
                 "tests/attr_bands.R", "run/patient_card.R"],
        owned=["R/13_attribution.R", "R/13b_attribution_floor.R", "R/14_attribution_eval.R",
               "tests/attr_replicates.R", "tests/attr_metrics.R", "tests/attr_external.R",
               "tests/attr_bands.R", "run/patient_card.R"],
        chart="docs/figs/pipeline_attribution.tex",
        # same rule as the engine's allow-list (see the comment there)
        allow={
            ".opt": "command-line option reader, one copy per script",
        },
        utility_files=[]),
}

def match_any(path, patterns):
    return any(fnmatch(path, p) for p in patterns)

def load_tables():
    defs = list(csv.DictReader(open(os.path.join(CG, "definitions.csv"), encoding="utf-8")))
    edges = list(csv.DictReader(open(os.path.join(CG, "edges.csv"), encoding="utf-8")))
    targets = list(csv.DictReader(open(os.path.join(CG, "targets.csv"), encoding="utf-8")))
    return defs, edges, targets

def target_tables(edges, targets):
    """TCALL[target] = functions its body calls; TDEP[target] = targets it reads."""
    # A tar_target is keyed by its name; a script statement by name@file, since
    # the run scripts reuse the same local names.
    def tkey(name, file):
        return name if file == "_targets.R" else f"{name}@{file}"
    tcall, tdep, writes = defaultdict(set), defaultdict(set), defaultdict(set)
    for e in edges:
        if not e["target"]:
            continue
        t = tkey(e["target"], e["caller_file"])
        if e["kind"] == "target":
            tdep[t].add(tkey(e["callee"], e["caller_file"]))
        elif e["kind"] == "target_write":
            writes[tkey(e["callee"], e["caller_file"])].add(t)
        elif e["kind"] in ("call", "ref", "string"):
            tcall[t].add(e["callee"])
    # a statement that writes x is read by whatever reads x
    for t in list(tdep):
        for x in list(tdep[t]):
            tdep[t] |= writes.get(x, set()) - {t}
    names = []
    for t in targets:
        k = tkey(t["name"], t["file"])
        if k not in names: names.append(k)
    return tcall, tdep, names

def resolve(callee, caller_file, by_name):
    """Pick the definition a name refers to: the same file first; else a
    top-level definition under R/ (sourced everywhere); else a top-level
    definition in a `*_common.R` file of the caller's directory (sourced by
    the scripts beside it). Anything else is a name clash with some other
    script's local helper and is not an edge."""
    cands = by_name.get(callee, [])
    if not cands:
        return None
    same = [d for d in cands if d["file"] == caller_file]
    if same:
        return same[0]
    glob = [d for d in cands if d["file"].startswith("R/") and d["top_level"] == "TRUE"]
    if glob:
        return glob[0]
    d0 = os.path.dirname(caller_file)
    common = [d for d in cands if d["top_level"] == "TRUE" and d["file"].endswith("_common.R")
              and os.path.dirname(d["file"]) == d0]
    if common:
        return common[0]
    return None

def analyse(arm_name, arm, defs, edges, targets):
    by_name = defaultdict(list)
    for d in defs:
        if d["anonymous"] == "FALSE":
            by_name[d["name"]].append(d)
    key = lambda d: (d["name"], d["file"])
    owned_defs = {key(d): d for d in defs if d["anonymous"] == "FALSE" and match_any(d["file"], arm["owned"])}

    # adjacency on (name, file) keys; script nodes are ("<script:FILE>", FILE)
    adj = defaultdict(list)     # node -> list of (callee_key, kind)
    ext = defaultdict(set)      # node -> external calls
    for e in edges:
        caller_key = (e["caller"], e["caller_file"])
        if e["kind"] == "external":
            ext[caller_key].add(e["callee"]); continue
        tgt = resolve(e["callee"], e["caller_file"], by_name)
        if tgt is None:
            continue
        adj[caller_key].append((key(tgt), e["kind"]))

    entry_nodes = [n for n in adj if n[0].startswith("<script:") and match_any(n[1], arm["entries"])]
    reach, boundary, frontier = {}, {}, []
    for n in entry_nodes:
        for tgt, kind in adj[n]:
            frontier.append((tgt, kind, n))
    while frontier:
        tgt, kind, src = frontier.pop()
        if tgt in owned_defs:
            if tgt in reach:
                reach[tgt]["callers"].add((src, kind)); continue
            reach[tgt] = {"callers": {(src, kind)}}
            for t2, k2 in adj.get(tgt, []):
                frontier.append((t2, k2, tgt))
        else:
            boundary.setdefault(tgt, set()).add((src, kind))

    # Nested definitions (closures inside a function) are part of their
    # enclosing function's box and are reported under it, not required on the
    # chart. Coverage is asked of top-level definitions only.
    util = {k for k in reach if match_any(k[1], arm["utility_files"])}
    nested = {k: owned_defs[k]["enclosing"] for k in reach if owned_defs[k]["top_level"] == "FALSE"}
    main = sorted((k for k in reach if k not in util and k not in nested), key=lambda k: (k[1], k[0]))
    unreached = sorted((k for k in owned_defs if k not in reach and k not in util
                        and owned_defs[k]["top_level"] == "TRUE"), key=lambda k: (k[1], k[0]))
    externals = sorted({c for n in list(reach) + entry_nodes for c in ext.get(n, [])})
    return dict(main=main, util=sorted(util), nested=nested, unreached=unreached, boundary=boundary,
                reach=reach, adj=adj, entries=entry_nodes, externals=externals, owned=owned_defs,
                targets=target_tables(edges, targets))

def write_md(arm_name, arm, r):
    p = os.path.join(CG, f"{arm_name}.md")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(f"# Call graph of the {arm_name} arm (from R's parser)\n\n")
        fh.write("Generated by `docs/figs/callgraph_check.py` from the tables of "
                 "`docs/figs/callgraph_parse.R`. Owned files: " +
                 ", ".join(f"`{o}`" for o in arm["owned"]) + ". Entry scripts: " +
                 ", ".join(f"`{n[1]}`" for n in r["entries"]) + ".\n\n")
        fh.write("Edge kinds: call = `name(...)`, ref = bare name passed as an argument, "
                 "string = name given as a string.\n\n")
        fh.write(f"## Reachable functions of the arm ({len(r['main'])})\n\n")
        fh.write("| Function | File | Called by | Calls |\n|---|---|---|---|\n")
        for k in r["main"]:
            cb = ", ".join(sorted({f"`{s[0]}`" + ("" if kd == "call" else f" ({kd})")
                                   for s, kd in r["reach"][k]["callers"]}))
            cs = ", ".join(sorted({f"`{t[0]}`" + ("" if kd == "call" else f" ({kd})")
                                   for t, kd in r["adj"].get(k, []) if t in r["reach"] and t not in set(r["util"])}))
            fh.write(f"| `{k[0]}` | {k[1]} | {cb} | {cs} |\n")
        if r["util"]:
            fh.write("\nUtilities reached: " + ", ".join(f"`{k[0]}`" for k in r["util"]) + "\n")
        if r["nested"]:
            fh.write(f"\n## Nested helpers, folded into their enclosing function ({len(r['nested'])})\n\n")
            by_parent = defaultdict(list)
            for k, parent in r["nested"].items():
                by_parent[parent].append(k[0])
            for parent in sorted(by_parent):
                fh.write(f"- `{parent}`: " + ", ".join(f"`{n}`" for n in sorted(by_parent[parent])) + "\n")
        owned_names = {k[0] for k in r["main"]} | {k[0] for k in r["util"]}
        tcall, tdep, tnames = r["targets"]
        in_entry = lambda t: match_any(t.split("@", 1)[1] if "@" in t else "_targets.R", arm["entries"])
        touching = [t for t in tnames if tcall.get(t, set()) & owned_names and in_entry(t)]
        if touching:
            fh.write(f"\n## Targets and entry-script statements whose body calls a function of the arm ({len(touching)})\n\n")
            fh.write("Data flow between functions that never call each other runs through these. "
                     "A tar_target is keyed by its name, a script statement by the symbol it assigns "
                     "(or `<expr:FILE:LINE>`) followed by `@FILE`. `reads` lists the other targets the body "
                     "uses; `read by` lists the targets that use this one.\n\n")
            fh.write("| Target | Calls | Reads | Read by |\n|---|---|---|---|\n")
            for t in touching:
                rb = sorted(t2 for t2 in tnames if t in tdep.get(t2, set()))
                fh.write(f"| `{t}` | " + ", ".join(f"`{c}`" for c in sorted(tcall[t] & owned_names)) + " | " +
                         ", ".join(f"`{c}`" for c in sorted(tdep.get(t, set()))) + " | " +
                         ", ".join(f"`{c}`" for c in rb) + " |\n")
        fh.write(f"\n## Boundary: functions of other arms that this arm calls ({len(r['boundary'])})\n\n")
        fh.write("| Function | File | Called from |\n|---|---|---|\n")
        for k in sorted(r["boundary"], key=lambda k: (k[1], k[0])):
            fh.write(f"| `{k[0]}` | {k[1]} | " + ", ".join(sorted({f"`{s[0]}`" for s, _ in r["boundary"][k]})) + " |\n")
        fh.write(f"\n## Owned but not reached from the entry scripts ({len(r['unreached'])})\n\n")
        fh.write(", ".join(f"`{k[0]}` ({k[1]})" for k in r["unreached"]) + "\n" if r["unreached"] else "none\n")
        fh.write(f"\n## External package calls ({len(r['externals'])})\n\n")
        fh.write(", ".join(f"`{c}`" for c in r["externals"]) + "\n")
    return p

def check_chart(arm_name, arm, r):
    tex = open(os.path.join(ROOT, arm["chart"]), encoding="utf-8").read()
    names = {t.replace("\\_", "_") for t in re.findall(r"([A-Za-z_.][A-Za-z0-9_.\\]*)\(", tex)}
    missing = [k for k in r["main"] if k[0] not in names and k[0] not in arm["allow"]]
    stale = [n for n in arm["allow"] if n in names]
    print(f"[{arm_name}] reachable owned functions: {len(r['main'])} (+{len(r['util'])} utilities); "
          f"named on the chart: {sum(1 for k in r['main'] if k[0] in names)}; allow-listed: {len(arm['allow'])}; "
          f"boundary: {len(r['boundary'])}; owned but unreached: {len(r['unreached'])}")
    if stale:
        print(f"[{arm_name}] allow-listed but now drawn (remove from allow): " + ", ".join(stale))
    if missing:
        print(f"[{arm_name}] NOT on the chart and NOT allow-listed:")
        for k in missing:
            print(f"    {k[0]:34s} {k[1]:28s} called by " +
                  ", ".join(sorted({s[0] for s, _ in r['reach'][k]['callers']})))
    else:
        print(f"[{arm_name}] the chart names every reachable owned function outside the allow-list")
    return not missing

# --- the chart's arrows against the tables -----------------------------------
NAME_RE = re.compile(r"([A-Za-z_.][A-Za-z0-9_.\\]*)\(")

def _unesc(s):
    return s.replace("\\_", "_")

def _brace_args(s, start):
    """Arguments of a macro call whose first '{' is at s[start]; returns (args, end)."""
    args, i = [], start
    while i < len(s) and s[i] == "{":
        depth, j = 0, i
        while j < len(s):
            if s[j] == "{": depth += 1
            elif s[j] == "}":
                depth -= 1
                if depth == 0: break
            j += 1
        args.append(s[i + 1:j]); i = j + 1
        while i < len(s) and s[i] in " \n\t": i += 1
    return args, i

def parse_chart(tex):
    """Boxes and arrows of a chart.  boxes[id] = dict(num, fns, inner, objs);
    arrows = list of (style, u, v)."""
    boxes = {}
    for m in re.finditer(r"\\node\[[^\]]*\]\s*\((\w+)\)\s*(?:at\s*\([^)]*\)\s*)?\{", tex):
        nid, i = m.group(1), m.end() - 1
        body, _ = _brace_args(tex, i)
        body = body[0]
        b = dict(num="", fns=set(), inner=set(), objs=set())
        mm = re.match(r"\\(B|Bs|Obj)\{", body)
        if mm:
            args, _ = _brace_args(body, mm.end() - 1)
            b["num"] = args[0].strip()
            if mm.group(1) == "Obj" and len(args) >= 2:
                b["objs"] = {_unesc(x.strip()) for x in args[1].split(",")}
                b["fns"] = {_unesc(x) for x in NAME_RE.findall(args[2])} if len(args) > 2 else set()
            elif len(args) >= 3:
                b["fns"] = {_unesc(x) for x in NAME_RE.findall(args[2])}
                if mm.group(1) == "B" and len(args) >= 4:
                    b["inner"] = {_unesc(x) for x in NAME_RE.findall(args[3])}
        boxes[nid] = b
    arrows = []
    for m in re.finditer(r"\\draw\[([^\]]*)\]([^;]*);", tex):
        style = m.group(1).split(",")[0].strip()
        refs = [x for x in re.findall(r"\((\w+)(?:\.[a-z ]+)?\)", m.group(2)) if x in boxes]
        if len(refs) >= 2:
            arrows.append((style, refs[0], refs[-1]))
    return boxes, arrows

def check_arrows(arm_name, arm, r, defs, edges):
    tex = open(os.path.join(ROOT, arm["chart"]), encoding="utf-8").read()
    boxes, arrows = parse_chart(tex)
    tcall, tdep, tnames = r["targets"]
    tset = set(tnames)
    # name-level call relation, nested helpers folded into their top-level function
    encl = {}
    for d in defs:
        if d["anonymous"] == "FALSE" and d["top_level"] == "FALSE":
            encl[d["name"]] = d["enclosing"]
    def top(n):
        seen = set()
        while n in encl and n not in seen:
            seen.add(n); n = encl[n]
        return n
    calls = defaultdict(set)
    for e in edges:
        if e["kind"] in ("call", "ref", "string") and not e["caller"].startswith("<script:"):
            calls[top(e["caller"])].add(e["callee"])
    ft = defaultdict(set)                       # function -> targets whose body calls it
    for t, fs in tcall.items():
        for f in fs: ft[f].add(t)
    def reads(t):
        return tdep.get(t, set())
    def objs_of(box):
        """An object box names tar_targets, or script values of an entry
        script (`bundle` in run/test_look.R is the target `bundle@run/test_look.R`)."""
        out = set()
        for o in box["objs"]:
            if o in tset: out.add(o)
            out |= {t for t in tset if t.startswith(o + "@") and match_any(t.split("@", 1)[1], arm["entries"])}
        return out

    def why(u, v):
        U, V = boxes[u], boxes[v]
        for fv in V["fns"]:
            if calls[fv] & U["fns"]: return f"{fv}() calls " + ", ".join(sorted(calls[fv] & U["fns"]))
        for fu in U["fns"]:
            if calls[fu] & V["fns"]: return f"{fu}() calls " + ", ".join(sorted(calls[fu] & V["fns"]))
        TU = objs_of(U) | {t for f in U["fns"] for t in ft[f]}
        TV = objs_of(V) | {t for f in V["fns"] for t in ft[f]}
        for tv in TV:
            hit = reads(tv) & TU
            if hit: return f"target {tv} reads " + ", ".join(sorted(hit))
        for tv in objs_of(V):
            if tcall.get(tv, set()) & U["fns"]: return f"target {tv} calls " + ", ".join(sorted(tcall[tv] & U["fns"]))
        for t in TU & TV:                       # one target body calls a function of each box
            if tcall.get(t, set()) & U["fns"] and tcall.get(t, set()) & V["fns"]:
                return f"both called within target {t}"
        for f, cs in calls.items():             # one function body calls a function of each box
            if cs & U["fns"] and cs & V["fns"]:
                return f"both called within {f}()"
        for fv in V["fns"]:                     # v reached through one intermediate function
            for f2, cs in calls.items():
                if fv in cs:
                    for t in ft[f2]:
                        if reads(t) & TU:
                            return f"{fv}() called by {f2}(), called by target {t}, which reads " + ", ".join(sorted(reads(t) & TU))
        for fu in U["fns"]:                     # an object produced by a target through one intermediate function
            for f2, cs in calls.items():
                if fu in cs and ft[f2] & objs_of(V):
                    return f"target {sorted(ft[f2] & objs_of(V))[0]} calls {f2}(), which calls {fu}()"
        return None

    if PAIRS is not None:                      # --pairs: query mode, no verdict
        ids = list(boxes)
        for u in ids:
            for v in ids:
                if u == v or (PAIRS and u not in PAIRS and v not in PAIRS): continue
                w = why(u, v)
                if w: print(f"    {boxes[u]['num'] or u:>5} -> {boxes[v]['num'] or v:<5} ({u} -> {v}): {w}")
        return True
    verified, declared, bad = [], [], []
    for style, u, v in arrows:
        if style == "en":
            continue
        if style == "ex":
            declared.append((u, v)); continue
        w = why(u, v)
        (verified if w else bad).append((u, v, w))
    # call relations between boxes with no arrow and no textual containment
    missing = []
    ids = list(boxes)
    drawn = {(u, v) for _, u, v in arrows} | {(v, u) for _, u, v in arrows}
    for u in ids:
        for v in ids:
            if u >= v or (u, v) in drawn: continue
            U, V = boxes[u], boxes[v]
            rel = {f"{a}()->{b}()" for a in U["fns"] for b in calls[a] & V["fns"]} | \
                  {f"{a}()->{b}()" for a in V["fns"] for b in calls[a] & U["fns"]}
            if not rel: continue
            if (U["fns"] & (V["inner"] | V["fns"])) or (V["fns"] & (U["inner"] | U["fns"])): continue
            missing.append((boxes[u]["num"] or u, boxes[v]["num"] or v, sorted(rel)))
    num = lambda x: boxes[x]["num"] or x
    print(f"[{arm_name}] arrows: {len(verified)} backed by the tables, {len(declared)} declared cross-process (ex), "
          f"{len(bad)} NOT backed")
    for u, v, _ in bad:
        print(f"    {num(u)} -> {num(v)}   ({u} -> {v}): no call, target read or target call links these boxes")
    for u, v in declared:
        print(f"    declared: {num(u)} -> {num(v)}")
    for u, v, rel in missing:
        print(f"    warning: {u} and {v} are linked by " + "; ".join(rel) + " but no arrow joins them")
    return not bad

PAIRS = None   # --pairs[=id1,id2,...]: print every backed pair among the chart's boxes (or those touching the ids)

if __name__ == "__main__":
    for a in sys.argv[1:]:
        if a.startswith("--pairs"):
            PAIRS = [x for x in a.partition("=")[2].split(",") if x]
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--no-parse" not in sys.argv:
        subprocess.run(["Rscript", "docs/figs/callgraph_parse.R"], cwd=ROOT, check=True)
    defs, edges, targets = load_tables()
    ok = True
    for arm_name, arm in ARMS.items():
        if args and arm_name not in args:
            continue
        r = analyse(arm_name, arm, defs, edges, targets)
        p = write_md(arm_name, arm, r)
        if arm["chart"]:
            ok = check_chart(arm_name, arm, r) and ok
            ok = check_arrows(arm_name, arm, r, defs, edges) and ok
        else:
            print(f"[{arm_name}] reachable owned functions: {len(r['main'])}; boundary: {len(r['boundary'])}; "
                  f"owned but unreached: {len(r['unreached'])}; no chart yet -> {os.path.relpath(p, ROOT)}")
    sys.exit(0 if ok else 1)
