PY ?= python
SEEDS ?= 10

.PHONY: all env data bde smartcyp curate evaluate analyse figures report test clean clean-results help

help:
	@echo "one command per stage:"
	@echo "  make env       create the pinned conda environment (som-bench)"
	@echo "  make data      download the raw SDFs from the pinned upstream commit"
	@echo "  make bde       precompute ALFABET C-H bond dissociation energies"
	@echo "  make smartcyp  run the SMARTCyp baseline (needs SMARTCYP_JAR=...)"
	@echo "  make curate    standardize, relabel, deduplicate; write the report"
	@echo "  make evaluate  train and evaluate all models x splits x seeds"
	@echo "  make analyse   failure analysis from the saved predictions"
	@echo "  make figures   figures from the result tables"
	@echo "  make report    render the result tables into results/summary.md"
	@echo "  make test      run the test suite"
	@echo "  make all       data -> curate -> bde -> evaluate -> analyse -> figures -> report"
	@echo ""
	@echo "override the repeat count with e.g.  make evaluate SEEDS=3"

all: data curate bde evaluate analyse figures report

env:
	mamba env create -f environment.yml || conda env create -f environment.yml
	conda run -n som-bench pip install -e .
	@echo "activate with: conda activate som-bench"

data:
	$(PY) scripts/01_download.py

curate:
	$(PY) scripts/02_curate.py

# Runs in its own interpreter: ALFABET needs TensorFlow and an ancient
# scikit-learn pin, neither of which belongs in som-bench. See README.
ALFABET_PY ?= conda run -n alfabet-env python
bde:
	TF_USE_LEGACY_KERAS=1 $(ALFABET_PY) scripts/07_bde.py

# External published baseline. Needs a SMARTCyp jar, which this repo does not
# ship (LGPL-3.0): build it, then
#   make smartcyp SMARTCYP_JAR=/path/to/smartcyp.jar
SMARTCYP_JAR ?=
smartcyp:
	@test -n "$(SMARTCYP_JAR)" || { \
	  echo "set SMARTCYP_JAR=/path/to/smartcyp.jar (see scripts/08_smartcyp.py)"; \
	  exit 1; }
	$(PY) scripts/08_smartcyp.py --jar $(SMARTCYP_JAR)

evaluate:
	$(PY) scripts/03_evaluate.py --seeds $(SEEDS)

analyse:
	$(PY) scripts/04_analyse.py --seeds $(SEEDS)

figures:
	$(PY) scripts/05_figures.py

report:
	$(PY) scripts/06_report.py

test:
	$(PY) -m pytest -q

clean-results:
	rm -rf results figures data/interim/predictions

clean: clean-results
	rm -rf data
