SHELL := /usr/bin/env bash
PYTHON_VERSION := 3.11
PYTHON_VERSION_CONDENSED := 311
PACKAGE_NAME := lignova
REPO_PATH := $(shell git rev-parse --show-toplevel)
PACKAGE_PATH := $(REPO_PATH)/$(PACKAGE_NAME)
TESTS_PATH := $(REPO_PATH)/tests
CONDA_NAME := $(PACKAGE_NAME)-dev
CONDA := conda run -n $(CONDA_NAME)
CONDA_LOCK_OPTIONS := -p linux-64 -p osx-64 -p win-64 --channel conda-forge
DOCS_URL := https://durrantlab.github.io/lignova



###   ENVIRONMENT   ###

# See https://github.com/pypa/pip/issues/7883#issuecomment-643319919
export PYTHON_KEYRING_BACKEND := keyring.backends.null.Keyring

.PHONY: conda-create
conda-create:
	- conda deactivate
	conda remove -y -n $(CONDA_NAME) --all
	conda create -y -n $(CONDA_NAME)
	$(CONDA) conda install -y python=$(PYTHON_VERSION)
	$(CONDA) conda install -y conda-lock

# Default packages that we always need.
.PHONY: conda-setup
conda-setup:
	$(CONDA) conda install -y -c conda-forge poetry
	$(CONDA) conda install -y -c conda-forge pre-commit
	$(CONDA) conda install -y -c conda-forge tomli tomli-w
	$(CONDA) conda install -y -c conda-forge conda-poetry-liaison
	$(CONDA) conda install -c conda-forge openbabel


# Conda-only packages specific to this project.
.PHONY: conda-dependencies
conda-dependencies:
	$(CONDA) conda install -y -c conda-forge -c bioconda mmseqs2
	
.PHONY: conda-lock
conda-lock:
	- rm $(REPO_PATH)/conda-lock.yml
	$(CONDA) conda env export --from-history | grep -v "^prefix" > environment.yml
	$(CONDA) conda-lock -f environment.yml $(CONDA_LOCK_OPTIONS)
	rm $(REPO_PATH)/environment.yml
	$(CONDA) cpl-deps $(REPO_PATH)/pyproject.toml --env_name $(CONDA_NAME)
	$(CONDA) cpl-clean --env_name $(CONDA_NAME)

.PHONY: from-conda-lock
from-conda-lock:
	$(CONDA) conda-lock install -n $(CONDA_NAME) $(REPO_PATH)/conda-lock.yml
	$(CONDA) cpl-clean --env_name $(CONDA_NAME)

.PHONY: pre-commit-install
pre-commit-install:
	$(CONDA) pre-commit install

# Reads `pyproject.toml`, solves environment, then writes lock file.
.PHONY: poetry-lock
poetry-lock:
	$(CONDA) poetry lock --no-interaction

.PHONY: install
install:
	$(CONDA) poetry install --no-interaction
	- mkdir .mypy_cache
	- $(CONDA) poetry run mypy --install-types --non-interactive --explicit-package-bases $(PACKAGE_NAME)

.PHONY: environment
environment: conda-create from-conda-lock pre-commit-install install

.PHONY: refresh-locks
refresh-locks: conda-create conda-setup conda-lock pre-commit-install poetry-lock install



###   FORMATTING   ###

.PHONY: validate
validate:
	- $(CONDA) pre-commit run --all-files

.PHONY: formatting
formatting:
	- $(CONDA) isort --settings-path pyproject.toml ./
	- $(CONDA) black --config pyproject.toml ./




###   LINTING   ###


.PHONY: test
test:
	$(CONDA) pytest -c pyproject.toml --cov=$(PACKAGE_PATH) --cov-report=xml $(TESTS_PATH)

.PHONY: check-codestyle
check-codestyle:
	$(CONDA) isort --diff --check-only $(PACKAGE_PATH)
	$(CONDA) black --diff --check --config pyproject.toml $(PACKAGE_PATH)
	- $(CONDA) pylint --rcfile pyproject.toml $(PACKAGE_PATH)

.PHONY: mypy
mypy:
	-$(CONDA) mypy --config-file pyproject.toml $(PACKAGE_PATH)

.PHONY: lint
lint: check-codestyle mypy



###   CLEANING   ###

.PHONY: pycache-remove
pycache-remove:
	find . | grep -E "(__pycache__|\.pyc|\.pyo$$)" | xargs rm -rf

.PHONY: dsstore-remove
dsstore-remove:
	find . | grep -E ".DS_Store" | xargs rm -rf

.PHONY: mypycache-remove
mypycache-remove:
	find . | grep -E ".mypy_cache" | xargs rm -rf

.PHONY: ipynbcheckpoints-remove
ipynbcheckpoints-remove:
	find . | grep -E ".ipynb_checkpoints" | xargs rm -rf

.PHONY: pytestcache-remove
pytestcache-remove:
	find . | grep -E ".pytest_cache" | xargs rm -rf

.PHONY: build-remove
build-remove:
	rm -rf build/

.PHONY: cleanup
cleanup: pycache-remove dsstore-remove mypycache-remove ipynbcheckpoints-remove pytestcache-remove



###   BUILDING   ###

.PHONY: build
build:
	$(CONDA) poetry build


###   DOCS   ###

.PHONY: docs
docs:
	rm -rf ./docs/html/
	$(CONDA) sphinx-build -nT ./docs/source/ ./docs/html/
	touch ./docs/html/.nojekyll

.PHONY: docs-versioned
docs-versioned:
	rm -rf ./docs/html/
	$(CONDA) sphinx-multiversion -nT ./docs/source/ ./docs/html/
	touch ./docs/html/.nojekyll

	# Create html redirect to main
	echo "<head>" > ./docs/html/index.html
	echo "  <meta http-equiv='refresh' content='0; URL=$(DOCS_URL)/main/index.html'>" >> ./docs/html/index.html
	echo "</head>" >> ./docs/html/index.html

.PHONY: open-docs
open-docs:
	xdg-open ./docs/html/index.html 2>/dev/null
