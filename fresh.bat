# Clean up the directory after building the document (Windows version)

del HerringSR.aux
del HerringSR.knit.md
del HerringSR.log
del HerringSR.pdf
del HerringSR.Rmd
del HerringSR.tex
del HerringSR.upa
del HerringSR.upb
del HerringSR.utf8.md
del texput.log

rmdir /S /Q csas-style
rmdir /S /Q knitr-figs-pdf
rmdir /S /Q knitr-figs-docx
rmdir /S /Q _book
rmdir /S /Q lrp_ms
