## Test environments

* Tested on Local Windows machine with R 4.4.4
* There were no ERRORs or WARNINGs or NOTES. 

* RHub: Windows Server 2022, R-devel, 64 bit
* There were no ERRORs or WARNINGs. 3 Notes.

❯ checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

I think this is related to an issue in RHub.

* RHub: Ubuntu Linux 20.04.1 LTS, R-release, GCC
* There were no ERRORs or WARNINGs. 2 NOTES.  

* checking examples ... [18s/35s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
     user system elapsed
IAB 2.779   0.02   5.518
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable


* RHub: Fedora Linux, R-devel, clang, gfortran
* There were no ERRORs or WARNINGs. 2 Notes.

* checking examples ... [18s/34s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
     user system elapsed
IAB 3.001  0.017   6.266
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable


## Downstream dependencies

* There are currently no downstream dependencies for this package.
