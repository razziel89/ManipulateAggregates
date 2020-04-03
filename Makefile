# makefile for ManipulateAggregates

SHELL:=bash

default: clean
	$(MAKE) pydoc
# ===== distribution rules =====

.PHONY : dist
dist: clean_dist
	@echo "Creating source distribution"
	@python setup.py sdist

.PHONY : publish
publish: dist
	@echo "Publishing to PyPI"
	@python -m twine upload dist/*

# ===== clean rules =====

.PHONY : clean
clean : clean_doc clean_dist

.PHONY : clean_doc
clean_doc : 
	rm -rf docs

.PHONY : clean_dist
clean_dist: 
	rm -rf build dist

# ===== Python documetnation =====
# Sphinx documentation
# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = ManipulateAggregates
SOURCEDIR     = .
BUILDDIR      = docs/tmp_python

spinx-help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
pydoc: Makefile manipagg.rst energyscan.rst
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	mv -t docs/ $(BUILDDIR)/html/*
	rm -r $(BUILDDIR)
	touch docs/.nojekyll
	rm manipagg.rst energyscan.rst

manipagg.rst: Makefile ManipulateAggregates/bin/manipagg.py
	set -o pipefail && \
	manipagg --full-help | \
	sed -r 's/^ +/ /' | \
	awk '$$1 ~ /^--/{print ""; $$1=$$1 "\n"}{print}' | \
	sed 's/-/\\\\-/g' | \
	sed 's/~/-/g' | \
	sed -r 's/ $$/\n/' \
	> $@

energyscan.rst: Makefile ManipulateAggregates/bin/energyscan.py
	set -o pipefail&& \
	energyscan --help | \
	sed -r 's/^ +/ /' | \
	awk '$$1 ~ /^--/{print ""; $$1=$$1 "\n"}{print}' | \
	sed 's/-/\\\\-/g' | \
	sed 's/~/-/g' | \
	sed -r 's/ $$/\n/' \
	> $@
	echo >> $@
	echo "Default config file:" >> $@
	echo "--------------------" >> $@
	echo >> $@
	echo ".. code-block::" >> $@
	echo >> $@
	set -o pipefail && \
	energyscan --longhelp | \
	sed 's/^/    /' \
	>> $@
