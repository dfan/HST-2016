library(RMySQL)
library(shiny)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/Week_1/")
all_cons <- dbListConnections(MySQL())
for (con in all_cons)
  dbDisconnect(con)
runApp(paste(getwd(),"Incidentalome_Figure",sep="/"))
runApp(paste(getwd(),"ExAC",sep="/"))


#library(RMySQL)
#con <- dbConnect(MySQL(), user = 'root', unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock"),
#    password = 'root', dbname = 'kohane_lab', host = 'localhost')
#dbListTables(con)
