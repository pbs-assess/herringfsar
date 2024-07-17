export PATH := /Library/TeX/texbin:$(PATH)    # add LaTeX path
export PATH := /usr/local/bin:~/.local/bin:$(PATH) # add pandoc-citeproc-preamble

# Cluster targets
html:  index.html
#pdf:  index.pdf
docx:  index.docx

all: html docx

# pdf version (all sections)
#index.pdf: index.Rmd makefile
#	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, config_file = "_bookdown.yml", output_format = "bookdown::pdf_document2")'	

# html version (all sections)
index.html: index.Rmd makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, config_file = "_bookdown.yml", output_format = "bookdown::html_document2")'	

# docx version (all sections)
index.docx: index.Rmd makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, config_file = "_bookdown.yml", output_format = "csasdown::fsar_word")'
