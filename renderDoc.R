# Render tech report for IPHC survey analysis
# See more info on csasdown at:
# https://github.com/pbs-assess/csasdown

library(bookdown)
library(csasdown)


csasdown::render_sar(  input = "index.rmd",
                        clean = TRUE,
                        output_format = "csasdown::fsar_word" )

                        # # config_file = "_bookdown.yml",
                        

render_book(  input = "index.rmd",
                        clean = TRUE,
                        config_file = "_bookdown.yml",
                        output_format = "csasdown::fsar_word" )

# bookdown::render_book(  input = "index.rmd",
#                         clean = TRUE,
#                         config_file = "_bookdown.yml",
#                         output_format = "bookdown::html_document2" )

