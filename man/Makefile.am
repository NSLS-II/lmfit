man_MANS = lmmin.3 lmcurve.3
noinst_DATA = lmmin.html lmcurve.html

%.3 : %.pod
	pod2man -s 3 -c "lmfit manual" -n "$*"      $< > $@

%.html : %.pod
	pod2html --title="lmfit: a self-contained C library for Levenberg-Marquardt least-squares minimization and curve fitting" --noindex $< > $@
