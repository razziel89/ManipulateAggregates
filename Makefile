# makefile for libFireDeamon

#import build parameters
include make.vars

.PHONY : doc
doc:
	@which doxygen
	@echo "Creating documentation using doxygen..."
	export PATH=$(THISDIR)/doc_bin:$(PATH) && doxygen --verbose Doxyfile

.PHONY : clean
clean : clean_doc

.PHONY : clean_doc
clean_doc : 
	rm -rf docs

.PHONY : uninstall
uninstall : 
	@test -d $(PYTHONPREFIX)/$(MODULENAME)
	rm -r $(PYTHONPREFIX)/$(MODULENAME)
	@test -f $(PREFIX)/bin/hashsort
	rm $(PREFIX)/bin/hashsort
	@test -f $(PREFIX)/bin/energyscan
	rm $(PREFIX)/bin/energyscan

.PHONY : install
install:
	@cd $(THISDIR)
	mkdir -p $(PYTHONPREFIX)/$(MODULENAME)
	@test -d $(PYTHONPREFIX)/$(MODULENAME)
	@cp -av __init__.py $(PYTHONPREFIX)/$(MODULENAME)
	@cp -av collection $(PYTHONPREFIX)/$(MODULENAME)
	@cp -av energyscan $(PYTHONPREFIX)/$(MODULENAME)
	@cp -av manipulation $(PYTHONPREFIX)/$(MODULENAME)
	@cp -av orbitalcharacter $(PYTHONPREFIX)/$(MODULENAME)
	@echo "Deleting VIM swapfiles"
	@find $(PYTHONPREFIX)/$(MODULENAME) -type f -name ".*.swp" -exec rm {} +
	mkdir -p $(PREFIX)/bin
	@cp -av bin/hashsort.py $(PREFIX)/bin/hashsort
	sed -i 's/REPLACEMODULENAME/$(MODULENAME)/' $(PREFIX)/bin/hashsort
	@cp -av bin/energyscan.py $(PREFIX)/bin/energyscan
	sed -i 's/REPLACEMODULENAME/$(MODULENAME)/' $(PREFIX)/bin/energyscan
