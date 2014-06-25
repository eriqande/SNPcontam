# SNPcontam

This is an R package with some code implementing methods under development
to detect contaminated samples from SNP genotyes sampled in molecular
ecology contexts.  

Currently, the directory `not-package` contains things that will not go into
the final R package.  This includes the manuscript that we will be preparing
and several readings and resources.  

To typeset the manuscript (which is really just a stub at this point)
`cd` into the `not-package/manuscript` directory and issue this
command (note that you have to have a TeX system installed on your
system, as well as R, with the `knitr` package installed):
```
Rscript -e "knitr::knit2pdf('snp-contam-main.Rnw')"
```
The resulting PDF file, `snp-contam-main.pdf` should then be complete.

## Compiling manuscript with TeXShop
To compile the manuscript in TeXShop you have to add a knitr engine and do a
few things to try to make syncing work better:

1. First, copy the file in the repo `./not-package/extras/Knitr.engine` to
`~/Library/TeXShop/Engines/`.
2. Then make the directory `~/Library/TeXShop/Rscripts` and copy the file
`./not-package/extras/patchKnitrSynctex.R` from the repo to `~/Library/TeXShop/Rscripts`.
Note that this file is from https://github.com/jan-glx/patchKnitrSynctex.
3. Once you have done that you will have to restart TeXShop.  Then you ought to
choose TeXShop Preferences -> Typesetting -> Sync Method -> Both PDF sync and Search.

Note that synctex and knitr and TeXShop don't play very well with multiple Rnw files in a project.
At some point that might get figured out.  Until then we just use `\input{}` to input the 
child .Rnw files and don't have any R code chunks in the child files.  Instead, we do
all the R analysis in separate scripts and will save the output in an .rda  file and just
slurp the results out via, perhaps, a single .Rnw file that we input via knitr's child
functions, but which does not create any text output.

## Terms 

As a work partially of the United States Government, this package is in the
public domain within the United States. 

Of course, it would be awfully lame of anyone to take anything found within
here and use it as their own before we tried publishing any of this, should
we choose to do that.

See TERMS.md for more information.

