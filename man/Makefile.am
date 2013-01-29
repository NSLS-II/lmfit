man_MANS = lmfit.3 lmmin.3 lmcurve_fit.3
auxdir = $(datadir)/man/html
aux_DATA = lmfit.html

lmfit.3 lmmin.3 lmcurve_fit.3 : lmfit.pod
	pod2man -s 3 -c "lmfit manual" lmfit.pod > lmfit.3
	pod2man -s 3 -c "lmfit manual" -n "lmmin" lmfit.pod > lmmin.3
	pod2man -s 3 -c "lmfit manual" -n "lmcurve_fit" lmfit.pod > lmcurve_fit.3

lmfit.html : lmfit.pod
	pod2html --title="Levenberg-Marquardt least-squares minimization and curve fitting in C/C++" --noindex lmfit.pod > lmfit.html

