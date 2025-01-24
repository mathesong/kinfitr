# Makefile
.PHONY: reload clean install

reload:
	Rscript -e "if('kinfitr' %in% (.packages())) { detach('package:kinfitr', unload=TRUE); try(unloadNamespace('kinfitr'), silent=TRUE) }; devtools::install('.', force=TRUE); library(kinfitr)"

clean:
	Rscript -e "remove.packages('kinfitr'); devtools::clean_dll()"

install:
	Rscript -e "devtools::install('.', force=TRUE)"

interactive:
	R

debug: reload interactive
