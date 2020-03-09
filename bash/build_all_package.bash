#!/bin/bash

R --vanilla --slave -e "devtools::test()"
R --vanilla --slave -e "devtools::run_examples()"
R --vanilla --slave -e "devtools::document()"
R --vanilla --slave -e "pkgdown::build_site()"
