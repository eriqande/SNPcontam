library(rmarkdown)

# script to compile some of the documents and reports here

# check that you are in correct directory
if(any(!file.exists("not-package", "SNPcontam.Rproj"))) stop("It appears you are not running Compile.R in the correct directory")



# get the absolute path for the output directory
OUTD <- file.path(getwd(), "Compiled_Documents")


# render elena's notes and stuff in HTML and PDF formats
INPUTS <- c("not-package/notes/elena-work-log.md", "not-package/notes/work-plan.md")
lapply(INPUTS, function(x) {
  render(x,  output_format = "html_document", output_dir = OUTD)
  render(x,  output_format = "pdf_document", output_dir = OUTD)
})



