# muscle/ligament wrapping
This repository contains a subroutine for MSC ADAMS to implement the wrapping behaviour of tendons and ligaments. Usually, the ligaments are modelled as lines connecting the two attachment points, however, this causes inaccuracies if this line penetrates the bone surfaces. This subroutine uses a torus geometry to model the bone surfaces that the ligament might contact during the range of motion in order to calculate the ligament length change and the force in it. 
## files
* torus_wrap.f  contains the Fortran source file of the ADAMS GFORCE subroutine
* torus_req.f contains the Fortran source file of the ADAMS request subroutine
* torus_wrapping_subroutine.dll is the subroutine file that you can use in your ADAMS models. 
