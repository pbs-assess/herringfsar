export PATH := /Library/TeX/texbin:$(PATH)    # add LaTeX path
export PATH := /usr/local/bin:~/.local/bin:$(PATH) # add pandoc-citeproc-preamble

# Cluster targets
html:  SOG_MPevalReport.html

all: html



# html version (all sections)
SOG_MPevalReport.html: SOG_MPevalReport.Rmd makefile
	Rscript -e 'bookdown::render_book(input = "SOG_MPevalReport.Rmd", clean = TRUE, output_format = "bookdown::html_document2")'	
