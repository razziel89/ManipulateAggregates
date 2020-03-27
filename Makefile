# makefile for libFireDeamon

default: clean
	$(MAKE) pydoc

.PHONY : clean
clean : clean_doc

.PHONY : clean_doc
clean_doc : 
	rm -rf docs

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
pydoc: Makefile
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	mv -t docs/ $(BUILDDIR)/html/*
	rm -r $(BUILDDIR)
