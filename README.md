# Individual-level analyses as published in:  
***Barneche DR, Kulbicki M, Floeter SR, Friedlander AM, Maina J, Allen AP. Scaling metabolism from individuals to reef-fish communities at broad spatial scales. Ecology Letters 10.1111/ele.12309***  

One can click [here](http://doi.wiley.com/10.1111/ele.12309) for the article URL  

### Overview  
This repository contains all analyses, figures and tables presented at the above-mentioned paper. It deals only with the individual-level metabolic rates analyses, both in the [main text]() and the [SI]().  

* first I advise that you open your R GUI (R, RStudio, or starting on terminal) from the ELEbarneche2014.RData file...
* ...this is an empty file that sets the absolute path to your project so everything will work independently on any machine; 
* pay attention to the required packages at the beginning of each file;
* in the project root folder, some summary statistics are provided in the file analyses-00-all.R;
* a file called figures.R also reproduces all the figures 2, S1-S5 exactly as they are shown in the paper and online SI;
* the figures will be automatically placed in a folder called output (it is going to be automatically created for you);
* I have also provided a tables.R file that creates Tables 1, S1-S3 as presented in the paper and online SI; 
* analyses are broken down into individual files and are located in the re-run folder; they should all run independently;
* files that start with analyses-* followed by a number relate to analyses presented on the main text; 
* analyses presented at the SI are found within files that follow the pattern analyses-a*;
* I have provided .RData files in re-run to facilitate visualization of summary statistics, but feel free to re-run all;
* notice, though, that some of these analyses may take up to 2 days to run on a regular computer;
* finally, I created a file with hypothetical data demonstrating how one propagates model uncertainties from JAGS to the community level. The file is called hardwire-example.R.

### Important remarks: these analyses require specific versions of:
* R: version 2.15.3 (for lme4 analyses) or R version 3.0.2 (for R2jags analyses);
* package lme4: version 0.999999-0;
* package R2jags: version 0.04-01;
* JAGS 3.3.0 (for R2jags analyses)

### Please report if you run into problems or spot a bug on the code:
d13g0 DOT b4rn3ch3 AT mq DOT 3du DOT 4u (replace the 0 for o, 1 for i, 3 for e, 4 for a)  

### How to download this project for people not familiar with GitHub:  
* on the project main page on GitHub, click on the right-bottom corner where it says Download ZIP;  
