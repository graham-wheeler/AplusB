# AplusB
A + B Design Investigator for phase I dose-escalation studies

To install and run the AplusB application from your machine, you will require R (available for download free from https://cran.rstudio.com/) or RStudio (available for download free from https://www.rstudio.com/products/rstudio/download/).

1) Download the .zip file for AplusB using the "Download ZIP" button.

2) Extract the folder "AplusB-master" to a directory of your choosing.

3) Open R/RStudio and run the following code:

```
setwd(<directory where "AplusB-master" folder is stored>)
install.packages("shiny")
library(shiny)
runApp("AplusB-master")
```
