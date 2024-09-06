# Clean directory after knitting the document (Windows version)

del fsar.aux
del fsar.docx
del fsar.html
del fsar.knit.md
del fsar.log
del fsar.pdf
del fsar.Rmd
del fsar.tex
del fsar.upa
del fsar.upb
del fsar.utf8.md
del texput.log

rmdir /S /Q csas-style
rmdir /S /Q knitr-figs-pdf
rmdir /S /Q knitr-figs-docx
rmdir /S /Q _book
rmdir /S /Q _bookdown_files
