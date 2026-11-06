PY ?= python
SEEDS ?= 10

.PHONY: all env data curate evaluate analyse figures test clean clean-results help

help:
	@echo "one command per stage:"
	@echo "  make env       create the pinned conda environment (som-bench)"
	@echo "  make data      download the raw SDFs from the pinned upstream commit"
	@echo "  make curate    standardize, relabel, deduplicate; write the report"
	@echo "  make evaluate  train and evaluate all models x splits x seeds"
	@echo "  make analyse   failure analysis from the saved predictions"
	@echo "  make figures   figures from the result tables"
	@echo "  make report    render the result tables into results/summary.md"
	@echo "  make test      run the test suite"
	@echo "  make all       data -> curate -> evaluate -> analyse -> figures -> report"
	@echo ""
	@echo "override the repeat count with e.g.  make evaluate SEEDS=3"

all: data curate evaluate analyse figures report

env:
	mamba env create -f environment.yml || conda env create -f environment.yml
	conda run -n som-bench pip install -e .
	@echo "activate with: conda activate som-bench"

data:
	$(PY) scripts/01_download.py

curate:
	$(PY) scripts/02_curate.py

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
