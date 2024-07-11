# Render tech report for IPHC survey analysis
# See more info on csasdown at:
# https://github.com/pbs-assess/csasdown

library(bookdown)

# Render html
bookdown::render_book( 	input = "SOG_MPevalReport.Rmd",
                        clean = TRUE,
                        output_format = "bookdown::html_document2" )

