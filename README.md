# tmbpop

[![Travis build status](https://travis-ci.org/pbs-assess/tmbpop.svg?branch=master)](https://travis-ci.org/pbs-assess/tmbpop)
[![Project Status: WIP - Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

**In development. Not for use yet.**

R package for implementing TMB age-structured POP model

Install with:

```r
# install.packages("remotes")
remotes::install_github("pbs-assess/tmbpop")
```

See example model fitting:

```
library(tmbpop)
?tmbpop::fit
```
# TMB-POP
Template Model Builder Code for Pacific Ocean Perch

## Background

This is initially the TMB and R code that William Aeberhard wrote (as a contract) based on the equations from our [Pacific Ocean Perch stock assessment for the west coast of Vancouver Island](http://waves-vagues.dfo-mpo.gc.ca/Library/361330.pdf). 

Initially this is to easily share with fellow Pacific participants of the TMB course in Halifax, lest they go their own way and spend time starting from scratch. I ran this code on Friday in the workshop and it worked fine and was very fast. It will be a good starting point for future assessments.

I haven't had time (or expertise) to check it, test it, or think about it in detail. But that can now happen given our newly acquired TMB knowledge.

An ideal plan would be something like

1. Test and understand the code.
1. Compare to previous results.
1. If satisfied then this could be used for the next POP assessment, or in conjunction with exiting (Awatea) code. Using the TMB code for a full assessment will initially require a lot of setting up of directories and structure and writing code for figures and tables, but we can utilise the structures that Chris developed for hake and then herring, and have code we can adapt for the figures and tables.
1. Maybe there should be a knitr file that plots results for a single model run (for quick understanding, and checking of convergence etc.), and also an assessment file that builds the assessment (with sensitivities and bridging analyses as necessary).
1. I think I can even automate the writing of the equations (based on the input file). 

Thinking and discussing this stuff before starting will help a lot. 

If you want to contribute to this repository (in a way that may help us) then please **fork** (on GitHub) and then **clone** to your machine. If you don't know what that means, or you just want to take the code and use it yourself (i.e. you may make changes but don't anticipate they will be useful for anyone else), then you can just download the repository (see the green 'Clone or download' button and download the .zip file).  Ideally I think we'd want to keep this repository fairly 'clean', in that changes are expected to be shared and merged into the main repository.
