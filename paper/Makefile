# grrt RNAAS paper build (Phase 3D).
#
#   make figures   regenerate both figures into figures/  (Phase 3D.2; gate 4)
#   make deps      create paper/.venv and install figure dependencies
#   make paper     compile manuscript.pdf                  (Phase 3D.3)
#   make wordcount body word-count gate                    (Phase 3D.3)
#   make clean     remove LaTeX/build artifacts (keeps figures/.gitkeep)
#   make veryclean also remove generated PDFs and the venv
#
# Figures and the paper build read only ../output/sweep.csv and the
# ../output/sweep_*.fits frames; no hand-edited plot data.

PYTHON ?= python3
VENV := .venv
VENV_PY := $(VENV)/bin/python
REQS := scripts/requirements.txt
FIGURES := figures/asymmetry_vs_epsilon3.pdf figures/image_gallery.pdf

# LaTeX engine autodetect (Phase 3D.3). Empty if none installed.
LATEXMK := $(shell command -v latexmk 2>/dev/null)
PDFLATEX := $(shell command -v pdflatex 2>/dev/null)
TECTONIC := $(shell command -v tectonic 2>/dev/null)

.PHONY: all figures deps paper wordcount clean veryclean

all: paper

# --- Phase 3D.2: figures -------------------------------------------------

$(VENV_PY): $(REQS)
	$(PYTHON) -m venv $(VENV)
	$(VENV_PY) -m pip install --quiet --upgrade pip
	$(VENV_PY) -m pip install --quiet -r $(REQS)
	@touch $(VENV_PY)

deps: $(VENV_PY)

figures: deps
	$(VENV_PY) scripts/make_figures.py

# --- Phase 3D.3: manuscript ---------------------------------------------

paper: figures manuscript.tex
	@if [ -n "$(LATEXMK)" ]; then \
		latexmk -pdf -interaction=nonstopmode manuscript.tex; \
	elif [ -n "$(TECTONIC)" ]; then \
		tectonic manuscript.tex; \
	elif [ -n "$(PDFLATEX)" ]; then \
		pdflatex -interaction=nonstopmode manuscript.tex && \
		(bibtex manuscript || true) && \
		pdflatex -interaction=nonstopmode manuscript.tex && \
		pdflatex -interaction=nonstopmode manuscript.tex; \
	else \
		echo "make: *** no LaTeX engine found (need latexmk, pdflatex, or tectonic)" >&2; \
		exit 1; \
	fi

wordcount: manuscript.tex
	$(PYTHON) scripts/wordcount.py manuscript.tex

# --- housekeeping --------------------------------------------------------

clean:
	rm -f *.aux *.log *.out *.bbl *.blg *.fls *.fdb_latexmk *.toc manuscript.pdf
	rm -rf scripts/__pycache__

veryclean: clean
	rm -f $(FIGURES)
	rm -rf $(VENV)
