About
-----

rFerns is an extended [random ferns](http://cvlab.epfl.ch/alumni/oezuysal/ferns.html) implementation for [R](http://r-project.org); in comparison to original, it can handle standard information system containing both categorical and continuous attributes. Moreover, it generates OOB error approximation and permutation-based attribute importance measure similar to [randomForest](http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm). Here is a [preprint of a paper with all the details](http://arxiv.org/abs/1202.1121).

How to use
---------

Quite fresh version should be on [CRAN](http://cran.r-project.org/web/packages/rFerns/index.html).

If you want to compile R package, fire `updpak.sh` and then execute `R CMD build rFerns`.

If you want to use it / test it apart from R, it is quite possible -- consult `side_src/test.c` to see how this may work.
Yet don't expect that this will ever become a standalone library.
