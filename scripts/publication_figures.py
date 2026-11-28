"""Generate aggregate-only draft scientific figures without new statistical inference."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile

TYPES = {"forest", "volcano", "heatmap", "locuscompare", "schematic"}
EDGE_TYPES = {"workflow", "documented", "hypothesis"}
BLUE, ORANGE, GRAY = "#0072B2", "#D55E00", "#757575"


def check(condition, message):
    if not condition:
        raise ValueError(message)


def nonempty(value):
    return isinstance(value, str) and bool(value.strip())


def num(value):
    value = float(value)
    check(math.isfinite(value), "Non-finite numeric input")
    return value


def prob(value):
    value = num(value)
    check(0 < value <= 1, "P value must be in (0, 1]; no zero imputation")
    return value


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def source_path(root, relative):
    check(nonempty(relative) and not Path(relative).is_absolute(), "Relative source path required")
    path = (root / relative).resolve()
    check(path.is_relative_to(root) and path.is_file(), "Source outside input root or missing")
    return path


def load(manifest, root):
    spec = json.loads(manifest.read_text(encoding="utf-8"))
    check(spec.get("schema_version") == 1, "Unsupported schema")
    check(spec.get("data_status") in {"synthetic", "aggregate"}, "Only synthetic or aggregate inputs")
    check(nonempty(spec.get("title")), "Title required")
    w, h = num(spec.get("width_mm", 180)), num(spec.get("height_mm", 180))
    check(80 <= w <= 220 and 60 <= h <= 300, "Unsupported figure dimensions")
    panels = spec.get("panels")
    check(isinstance(panels, list) and 1 <= len(panels) <= 6, "Use 1-6 panels")
    cols = spec.get("columns", 2)
    check(type(cols) is int and 1 <= cols <= min(3, len(panels)), "Invalid columns")
    records = []
    for p in panels:
        k = p.get("type")
        check(k in TYPES and nonempty(p.get("title")) and nonempty(p.get("reference")),
              "Panel needs valid type, title and source reference")
        record = {"type": k, "title": p["title"], "reference": p["reference"]}
        if k == "schematic":
            nodes, edges = p.get("nodes", []), p.get("edges", [])
            check(1 <= len(nodes) <= 12, "Use 1-12 schematic nodes")
            ids = set()
            for n in nodes:
                check(nonempty(n.get("id")) and n["id"] not in ids and nonempty(n.get("label")),
                      "Invalid/duplicate node")
                ids.add(n["id"])
                check(.12 <= num(n["x"]) <= .88 and .14 <= num(n["y"]) <= .88,
                      "Node outside safe canvas")
            for e in edges:
                check(e.get("from") in ids and e.get("to") in ids and e["from"] != e["to"],
                      "Invalid schematic edge")
                check(e.get("kind") in EDGE_TYPES, "Edge evidence class required")
                if e["kind"] == "documented":
                    check(nonempty(e.get("reference")), "Documented edge reference required")
            record["edges"] = edges
        else:
            path = source_path(root, p.get("source"))
            required = {"forest": ["label", "effect", "lower", "upper"],
                        "volcano": ["label", "log2fc", "p"], "heatmap": ["label"],
                        "locuscompare": ["variant", "p_gwas", "p_qtl"]}[k]
            with path.open(encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream, delimiter="\t")
                header = reader.fieldnames or []
                check(len(set(header)) == len(header) and set(required).issubset(header),
                      "Missing or duplicate TSV columns")
                rows = list(reader)
            check(rows and all(None not in row and all(nonempty(v) for v in row.values())
                               for row in rows), "Empty or malformed rows")
            if k == "forest":
                check(p.get("scale") in {"ratio", "difference"} and nonempty(p.get("effect_label")),
                      "Forest effect scale and label required")
                check(len(rows) <= 25, "Forest limited to 25 rows")
                for row in rows:
                    lo, est, hi = [num(row[f]) for f in ("lower", "effect", "upper")]
                    check(lo <= est <= hi and lo < hi, "Invalid interval")
                    check(p["scale"] != "ratio" or lo > 0, "Ratio interval must be positive")
            elif k == "volcano":
                check(p.get("p_label") in {"P value", "Adjusted P value", "FDR"},
                      "P value semantics required")
                check(prob(p.get("p_threshold", .05)) < 1 and num(p.get("fc_threshold", 1)) > 0,
                      "Invalid display thresholds")
                for row in rows:
                    num(row["log2fc"])
                    prob(row["p"])
            elif k == "heatmap":
                check(header[0] == "label" and 2 <= len(header) <= 16 and len(rows) <= 30
                      and nonempty(p.get("units")), "Heatmap needs label, numeric columns and units")
                for row in rows:
                    for col in header[1:]:
                        num(row[col])
            else:
                check(len({row["variant"] for row in rows}) == len(rows), "Duplicate variant")
                for row in rows:
                    prob(row["p_gwas"])
                    prob(row["p_qtl"])
            p["_header"], p["_rows"] = header, rows
            record.update(source=p["source"], sha256=digest(path), rows=len(rows))
        records.append(record)
    return spec, records


def plot(ax, p):
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    k = p["type"]
    ax.set_title(p["title"], loc="left", fontsize=9, fontweight="bold")
    if k == "schematic":
        ax.set(xlim=(0, 1), ylim=(0, 1))
        ax.axis("off")
        nodes = {n["id"]: n for n in p["nodes"]}
        for e in p.get("edges", []):
            a, b = nodes[e["from"]], nodes[e["to"]]
            ax.add_patch(FancyArrowPatch((num(a["x"]), num(a["y"])), (num(b["x"]), num(b["y"])),
                         arrowstyle="-|>", mutation_scale=10, shrinkA=24, shrinkB=24,
                         linestyle="--" if e["kind"] == "hypothesis" else "-",
                         color=ORANGE if e["kind"] == "hypothesis" else "#263746",
                         linewidth=1.1, zorder=1))
        for n in nodes.values():
            x, y = num(n["x"]), num(n["y"])
            ax.add_patch(FancyBboxPatch((x-.105, y-.07), .21, .14, boxstyle="round,pad=0.008",
                                        fc="#E8F2F7", ec=BLUE, lw=1, zorder=2))
            ax.text(x, y, n["label"], ha="center", va="center", fontsize=7, zorder=3)
        ax.text(.5, .015, "Dashed: hypothesis; solid: workflow/documented",
                ha="center", va="bottom", fontsize=6, transform=ax.transAxes)
        return
    rows = p["_rows"]
    if k == "forest":
        effects = [num(r["effect"]) for r in rows]
        errors = [[v-num(r["lower"]) for v, r in zip(effects, rows)],
                  [num(r["upper"])-v for v, r in zip(effects, rows)]]
        ax.errorbar(effects, range(len(rows)), xerr=errors, fmt="o", color=BLUE,
                    markersize=3.5, lw=1, capsize=2)
        ax.axvline(1 if p["scale"] == "ratio" else 0, color=GRAY, ls="--", lw=.8)
        if p["scale"] == "ratio":
            ax.set_xscale("log")
            from matplotlib.ticker import MaxNLocator, FixedLocator, FuncFormatter, NullLocator
            # Select a small number of readable ratio ticks on the log axis.
            lo, hi = ax.get_xlim()
            ticks = [float(t) for t in MaxNLocator(nbins=4).tick_values(lo, hi) if lo <= t <= hi and t > 0]
            if lo <= 1 <= hi and not any(abs(t - 1) < 1e-9 for t in ticks):
                ticks.append(1.0)
            ax.xaxis.set_major_locator(FixedLocator(sorted(ticks)))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.3g}"))
            ax.xaxis.set_minor_locator(NullLocator())
        ax.set_yticks(range(len(rows)), [r["label"] for r in rows])
        ax.invert_yaxis()
        ax.set_xlabel(p["effect_label"])
    elif k == "volcano":
        x = [num(r["log2fc"]) for r in rows]
        y = [-math.log10(prob(r["p"])) for r in rows]
        pt, fc = -math.log10(prob(p.get("p_threshold", .05))), num(p.get("fc_threshold", 1))
        colors = [ORANGE if a >= fc and b >= pt else BLUE if a <= -fc and b >= pt
                  else "#A4A4A4" for a, b in zip(x, y)]
        ax.scatter(x, y, c=colors, s=10, alpha=.8, linewidths=0)
        ax.axhline(pt, ls=":", lw=.8, color=GRAY)
        for v in (-fc, fc):
            ax.axvline(v, ls=":", lw=.8, color=GRAY)
        ax.set(xlabel="log2 fold change", ylabel="-log10(" + p["p_label"] + ")")
    elif k == "heatmap":
        columns = p["_header"][1:]
        matrix = [[num(r[c]) for c in columns] for r in rows]
        limit = max(abs(v) for row in matrix for v in row) or 1
        im = ax.imshow(matrix, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
        ax.set_xticks(range(len(columns)), columns, rotation=40, ha="right")
        ax.set_yticks(range(len(rows)), [r["label"] for r in rows])
        ax.figure.colorbar(im, ax=ax, fraction=.05, pad=.04, label=p["units"])
    else:
        ax.scatter([-math.log10(prob(r["p_gwas"])) for r in rows],
                   [-math.log10(prob(r["p_qtl"])) for r in rows], s=10, color=BLUE, alpha=.75)
        ax.set(xlabel="-log10(GWAS P)", ylabel="-log10(QTL P)")
        ax.text(.02, .98, "Descriptive only; no PP4 inference", transform=ax.transAxes,
                ha="left", va="top", fontsize=6)
    for name in ("top", "right"):
        ax.spines[name].set_visible(False)
    ax.tick_params(labelsize=7, length=3)


def render(manifest, root, output):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    import PIL
    from PIL import Image
    spec, records = load(manifest, root)
    check(not output.exists(), "Output exists: refusing overwrite")
    output.parent.mkdir(parents=True, exist_ok=True)
    cols = spec.get("columns", 2)
    rows = math.ceil(len(spec["panels"])/cols)
    banner = ("SYNTHETIC DEMO - NOT RESEARCH DATA" if spec["data_status"] == "synthetic"
              else "DRAFT - SCIENTIFIC REVIEW REQUIRED")
    matplotlib.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                                "svg.fonttype": "none", "pdf.fonttype": 42,
                                "axes.labelsize": 8, "svg.hashsalt": "publication-figures-v1"})
    size = [num(spec.get(k, 180))/25.4 for k in ("width_mm", "height_mm")]
    fig, axes = plt.subplots(rows, cols, squeeze=False, layout="constrained", figsize=size)
    try:
        fig.suptitle(spec["title"] + "\n" + banner, fontsize=9, fontweight="bold")
        for i, ax in enumerate(axes.flat):
            if i >= len(spec["panels"]):
                ax.axis("off")
            else:
                plot(ax, spec["panels"][i])
                ax.text(-.07, 1.04, chr(65+i), transform=ax.transAxes,
                        fontsize=10, fontweight="bold", va="bottom")
        with tempfile.TemporaryDirectory(prefix=".figures-", dir=output.parent) as tmp:
            folder = Path(tmp)
            exports = {}
            for ext in ("svg", "pdf", "png", "tiff"):
                path = folder / ("figure." + ext)
                kwargs = {"dpi": 600 if ext in {"png", "tiff"} else 300, "facecolor": "white"}
                if ext == "tiff":
                    kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
                if ext == "svg":
                    kwargs["metadata"] = {"Date": None}
                if ext == "pdf":
                    kwargs["metadata"] = {"CreationDate": None, "ModDate": None}
                fig.savefig(path, **kwargs)
                check(path.stat().st_size > 0, "Empty export")
                exports[path.name] = {"sha256": digest(path), "bytes": path.stat().st_size}
                if ext in {"png", "tiff"}:
                    with Image.open(path) as im:
                        check(im.format == ("PNG" if ext == "png" else "TIFF"), "Raster format mismatch")
                        check(all(abs(a-b*600) <= 2 for a, b in zip(im.size, size)), "Raster size mismatch")
                        dpi = im.info.get("dpi")
                        check(dpi and all(abs(float(v)-600) < 1 for v in dpi), "Raster DPI mismatch")
                        # TIFF metadata may contain Pillow IFDRational values, not JSON primitives.
                        exports[path.name].update(pixels=list(im.size), dpi=[float(v) for v in dpi])
            report = {"status": "DRAFT_REVIEW_REQUIRED", "schema_version": 1,
                      "data_status": spec["data_status"], "manifest_sha256": digest(manifest),
                      "code_sha256": digest(Path(__file__)), "matplotlib": matplotlib.__version__,
                      "pillow": PIL.__version__, "panels": records, "exports": exports,
                      "limitations": ["No new statistics or causal claim", "Scientific, visual and journal review required"]}
            (folder / "provenance.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
            folder.rename(output)
    finally:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    render(args.manifest.resolve(), args.input_root.resolve(), args.out.resolve())


if __name__ == "__main__":
    main()
