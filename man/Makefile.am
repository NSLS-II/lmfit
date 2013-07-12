man_MANS = lmfit.7 lmmin.3 lmcurve.3
noinst_DATA = lmfit.html lmmin.html lmcurve.html

lmfit.7 : lmfit.pod
	pod2man -s 7 -c "lmfit manual"              lmfit.pod   > lmfit.7

%.3 : %.pod
	pod2man -s 3 -c "lmfit manual" -n "$*"      $< > $@

%.html : %.pod
	pod2html --title="Levenberg-Marquardt least-squares minimization and curve fitting: a C library" --noindex $< > $@
