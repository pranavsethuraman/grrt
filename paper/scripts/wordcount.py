#!/usr/bin/env python3
"""Word-count / figure / reference gate for the RNAAS note (gate 6).

Usage: ``python wordcount.py manuscript.tex``

Counts the manuscript *body* and checks it against the limits in
docs/phase-3d-plan.md Sec. 2.5 / Sec. 5:

  * body words  <= 1000   (excludes the title, section headings, the
                           acknowledgments block, math, citation/label/
                           includegraphics arguments, and comments;
                           abstract and figure captions are included)
  * figures     <= 2       (``\\begin{figure}`` floats)
  * references  <= 15      (``\\bibitem`` entries, or @-entries in the
                           sibling refs.bib when BibTeX is used)

Exits non-zero if any limit is exceeded.
"""

import re
import sys
from pathlib import Path

BODY_LIMIT = 1000
FIGURE_LIMIT = 2
REF_LIMIT = 15

# Commands whose bracketed/braced arguments must not be counted.
DROP_WITH_ARGS = (
    "title|author|affiliation|section|subsection|"
    "citep|citet|cite|ref|eqref|label|includegraphics|"
    "bibliography|bibliographystyle|documentclass|usepackage"
)


def strip_comments(text: str) -> str:
    return re.sub(r"(?<!\\)%.*", "", text)


def body_words(tex: str) -> int:
    tex = strip_comments(tex)
    match = re.search(r"\\begin\{document\}(.*)\\end\{document\}", tex, re.S)
    body = match.group(1) if match else tex

    # Drop the acknowledgments block entirely.
    body = re.sub(r"\\begin\{acknowledgments\}.*?\\end\{acknowledgments\}",
                  " ", body, flags=re.S)
    # Drop commands together with their optional/mandatory arguments.
    body = re.sub(
        r"\\(?:" + DROP_WITH_ARGS + r")\*?\s*(?:\[[^\]]*\])*\{[^{}]*\}",
        " ", body, flags=re.S)
    # Drop math.
    body = re.sub(r"\$[^$]*\$", " ", body)
    body = re.sub(r"\\\[.*?\\\]", " ", body, flags=re.S)
    # Drop environment markers, remaining commands, and braces.
    body = re.sub(r"\\(?:begin|end)\{[^}]*\}(?:\[[^\]]*\])?", " ", body)
    body = re.sub(r"\\[a-zA-Z@]+\*?", " ", body)
    body = re.sub(r"[{}]", " ", body)

    return sum(1 for tok in body.split() if re.search(r"[A-Za-z]", tok))


def count_figures(tex: str) -> int:
    return len(re.findall(r"\\begin\{figure\*?\}", strip_comments(tex)))


def count_refs(tex_path: Path, tex: str) -> int:
    bibitems = len(re.findall(r"\\bibitem\b", tex))
    if bibitems:
        return bibitems
    bib = tex_path.parent / "refs.bib"
    if bib.exists():
        return len(re.findall(r"^\s*@\w+\s*\{", bib.read_text(), flags=re.M))
    return 0


def main(argv) -> int:
    if len(argv) != 2:
        print("usage: wordcount.py manuscript.tex", file=sys.stderr)
        return 2
    path = Path(argv[1])
    tex = path.read_text(encoding="utf-8")

    words = body_words(tex)
    figures = count_figures(tex)
    refs = count_refs(path, tex)

    print(f"body words : {words:4d}  (limit {BODY_LIMIT})")
    print(f"figures    : {figures:4d}  (limit {FIGURE_LIMIT})")
    print(f"references : {refs:4d}  (limit {REF_LIMIT})")

    failures = []
    if words > BODY_LIMIT:
        failures.append(f"body words {words} > {BODY_LIMIT}")
    if figures > FIGURE_LIMIT:
        failures.append(f"figures {figures} > {FIGURE_LIMIT}")
    if refs > REF_LIMIT:
        failures.append(f"references {refs} > {REF_LIMIT}")

    if failures:
        print("FAIL: " + "; ".join(failures), file=sys.stderr)
        return 1
    print("PASS: within RNAAS limits")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
