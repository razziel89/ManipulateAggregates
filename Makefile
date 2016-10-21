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
	@test -d $(PREFIX)/$(MODULENAME)
	rm -rf $(PREFIX)/$(MODULENAME)

.PHONY : install
install:
	@mkdir -p $(PREFIX)/$(MODULENAME)
	@test -d $(PREFIX)/$(MODULENAME)
	@cp -rv __init__.py $(PREFIX)/$(MODULENAME)
	@cp -rv collection $(PREFIX)/$(MODULENAME)
	@cp -rv energyscan $(PREFIX)/$(MODULENAME)
	@cp -rv manipulation $(PREFIX)/$(MODULENAME)
	@cp -rv orbitalcharacter $(PREFIX)/$(MODULENAME)
	@find $(PREFIX)/$(MODULENAME) -type f -name ".*.swp" -exec rm {} +
