# week_1
HST Summer Institute 2016

This week you will work with large-scale genomic datasets to help us better understand incidental findings in diverse populations, especially those that may be obtained upon routine clinical sequencing (e.g. for HCM). In addition to informing clinical interpretation of sequence variants (a huge issue right now), this will give you a chance to build your R/SQL skills, which will be useful throughout the summer and beyond. 

The goal is to create an interactive ‘Incidentalome App’ using the R/Shiny framework. You will use the app to grapple with the Incidentalome generally and understand genetic diversity, as well as a lens to understand the recently reported ‘Resilience’ finding (see literature below).  Please create and share a git repo to save/version all code you write (commit regularly!). Below are specific steps to get you started (feel free to skip the R/SQL steps you’ve already completed). After you’ve had a chance to build a basic app by completing the steps below, we’ll map out how to study the connections between the Incidentalome, genetic diversity, and ‘Resilience’.
1.	For your Science Communication class, summarize the following three papers (attached):
a. "The incidentalome: a threat to genomic medicine" by Kohane et al. JAMA 2006
b. “Genetic misdiagnoses and the potential for health disparities” by Manrai et al. NEJM 2016
c. “Analysis of 589,306 genomes identifies individuals resilient to severe Mendelian childhood diseases” by Chen et al. Nature Biotechnology 2016
2.	Install R [https://www.r-project.org/] and RStudio [https://www.rstudio.com/], and complete a basic tutorial to get comfortable with file input/output, manipulating data frames, and basic statistical operations in R. The following R tutorials are quite good: https://www.codeschool.com/courses/try-r and http://www.cyclismo.org/tutorial/R/
3.	If you don't already have one, download and install a MySQL/Apache/PHP stack e.g. MAMP (https://www.mamp.info/en/) or WAMP (http://www.wampserver.com/en/). The following SQL tutorials are great ways to get started with SQL: https://www.codecademy.com/learn/learn-sql and http://dev.mysql.com/doc/refman/5.7/en/tutorial.html
4.	Write a single paragraph on the following — What assumptions do Kohane et al. make to create Figure 1 in their JAMA 2006 paper? What parameters control the slope/concavity of the curve in Figure 1?
5.	Write an R function that reproduces Figure 1 in Kohane et al. Function inputs should include all parameters identified in 4.
6.	Familiarize yourself with R/Shiny (http://shiny.rstudio.com/tutorial/).
7.	Build an interactive app with R/Shiny to plot Kohane et al.’s Figure 1 as a function of the parameters you identified in 4 above.
8.	Download the ExAC allele frequency data (http://exac.broadinstitute.org/) for genes MYBPC3 and MYH7 (two key HCM genes)
9.	Stage the allele frequency data downloaded in 8 in a MySQL database locally.
10.	Write a SQL script to query the top variants in MYBPC3 and MYH7 by the absolute allele frequency difference between the AFR and NFE populations. Plot an empirical CDF of all allele frequency differences using R.
11.	Use RMySQL to connect your R/Shiny app to your MySQL db, and create basic UI functionality to allow the user to input an allele frequency threshold, and have your app report the fraction of variants in MYBPC3 and MYH7 separately that are differentiated at that threshold between AFR & NFE.