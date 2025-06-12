GWAsimulator Windows package

The executable file is standalone.  However, the feature of output file compression relies on an external program, gzip, that may not be available in many Windows systems.  If you need this feature and don't have gzip, you can install gzip in one of the following ways:

(1) Double click the copygzip.bat file included in the package to copy the gzip.exe file into the folder c:\windows\system32\ .  (This gzip program is version 1.3.12 and was obtained from the website in (2)).

(2) Install the gzip program from http://gnuwin32.sourceforge.net/packages/gzip.htm

(3) Install Cygwin with gzip included.  See http://www.cygwin.com for details.
