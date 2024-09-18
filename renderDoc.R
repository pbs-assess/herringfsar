# Render tech report for IPHC survey analysis
# See more info on csasdown at:
# https://github.com/pbs-assess/csasdown

library(bookdown)
library(csasdown)

render_book(  input = "index.Rmd",
              clean = FALSE,
              config_file = "_bookdown.yml",
              output_format = "csasdown::fsar_word" )


# csasdown::render_sar(  input = "index.Rmd",
#                         clean = TRUE,
#                         output_format = "csasdown::fsar_word" )

                        # # config_file = "_bookdown.yml",
                        



# bookdown::render_book(  input = "index.rmd",
#                         clean = TRUE,
#                         config_file = "_bookdown.yml",
#                         output_format = "bookdown::html_document2" )

