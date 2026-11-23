#!/usr/bin/env python3
"""Build the external Acme DAG with Graphviz, dot2tex, and LaTeX.

Requires Python's standard library, dot, dot2tex, and pdflatex with TikZ.
Run from any directory; --document also renders the PDF Quarto walkthrough.
Generated TeX, layout JSON and logs stay in ../tmp/pdfs/acme-dag.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent
WORK = ROOT.parent / "tmp/pdfs/acme-dag"


def command(name):
    """Resolve a build dependency with a clear missing-tool message."""
    path = shutil.which(name)
    if path is None:
        raise SystemExit(f"Missing {name}; install it and put it on PATH.")
    return path


def make_figure():
    """Compile DOT through dot2tex; verify the requested top-row ordering."""
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / ".gitignore").write_text("*\n")
    layout = subprocess.run([command("dot"), "-Tjson", str(ROOT / "acme.dot")],
                            check=True, capture_output=True, text=True)
    (WORK / "acme-layout.json").write_text(layout.stdout)
    nodes = {n["name"]: n for n in json.loads(layout.stdout)["objects"] if "pos" in n}
    row = [list(map(float, nodes[k]["pos"].split(",")))
           for k in ("Aju", "Yju", "Aad", "Yad")]
    assert max(p[1] for p in row) - min(p[1] for p in row) < .01
    assert all(row[i][0] < row[i+1][0] for i in range(3))
    picture = WORK / "acme-picture.tex"
    subprocess.run([command("dot2tex"), "--format", "tikz", "--figonly",
                    "--nodeoptions", r"minimum size=51.84bp,inner sep=0pt,font=\large",
                    "--output", str(picture), str(ROOT / "acme.dot")], check=True)
    # dot2tex emits 'filled' when fill and stroke colors match. Define that
    # shorthand explicitly rather than changing the generated TikZ by hand.
    tex = WORK / "acme-dag.tex"
    tex.write_text(r"""\documentclass[tikz,border=9pt]{standalone}
\usepackage{amsmath,amssymb,lmodern}
\usetikzlibrary{arrows.meta,shapes.geometric}
\tikzset{filled/.style={fill=strokecolor}}
\begin{document}
""" + picture.read_text() + "\n\\end{document}\n")
    proc = subprocess.run([command("pdflatex"), "-interaction=nonstopmode", "-halt-on-error",
                           f"-output-directory={WORK}", str(tex)], cwd=ROOT,
                          capture_output=True, text=True)
    (WORK / "figure-build.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode:
        raise SystemExit(f"LaTeX figure failed; see {WORK / 'figure-build.log'}")
    for directory in (ROOT / "figures", ROOT / "output/pdf"):
        directory.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(WORK / "acme-dag.pdf", directory / "acme-dag.pdf")
    print(ROOT / "output/pdf/acme-dag.pdf")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--document", action="store_true")
    args = parser.parse_args()
    make_figure()
    if args.document:
        subprocess.run([command("quarto"), "render", "acme.qmd", "--to", "pdf",
                        "--output-dir", "output/pdf"], cwd=ROOT, check=True)
        # Quarto refreshes its output directory; restore the standalone figure.
        shutil.copyfile(ROOT / "figures/acme-dag.pdf", ROOT / "output/pdf/acme-dag.pdf")


if __name__ == "__main__":
    main()
