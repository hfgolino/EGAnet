---
title: "cran-comments"
author: "Hudson Golino"
date: "5/07/2019"
output: pdf_document
---

## 06/03/2019
# Updated version of the EGAnet package
- This is an updated version of the EGAnet package. A new function was added (vn.entropy) and another was modified (entropyFit). 

- The EGA function was also updated to accomodate a new check for unidimensionality.

## 04/30/2019
# Resubmission
This is a resubmission. In this version I have:

* Replaced \dontrun{} to \donttest{}

* Did not include toy examples because this is not feasible. Even with simulated data, in the most simple conditions possible, the functions are taking > 5 sec to run. 

* Also, it is important to note that ALL authors and contributors are listed in the DESCRIPTION file. If someone was a co-author in a paper cited in the DESCRIPTION file, but did not contribute to the development of the package itself, this person was not included as an author/contributor. The reason for this is quite simple: some people had a role in the papers that were more substantive (i.e. helped in the theoretical background of the papers instead of helping with the package development).

* Corrected the URL below (From: man/itemConfirm.Rd)
     URL:
https://iopscience.iop.org/article/10.1088/1742-5468/2005/09/P09008/meta

 To:
 
<doi:10.1088/1742-5468/2005/09/P09008>

* Explained all acronyms (e.g. CFA, TMFG) in the Description field
to avoid misunderstandings.

* Also added references describing the methods in Description field.

* I haven't unwrapped the examples in \dontrun{}, because this is not feasible, because the simplest example in the package takes 0.88 sec to be executed. Also, small toy examples are not informative.  

* Corrected the DOIs in the Description field of the DESCRIPTION file.

* Updated the The Date field in the DESCRIPTION file.
