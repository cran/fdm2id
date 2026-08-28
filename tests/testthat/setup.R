# Everything the tests draw goes to a null device.
#
# Without one, a plot with no device open makes R open its fallback, Rplots.pdf, in the working
# directory -- and testthat sets that to this folder, so the file lands in the sources. It was
# picked up by a commit once already.
#
# Tests that open their own device and close it on exit are unaffected: this one sits
# underneath them, and the session closes it on the way out.
grDevices::pdf (NULL)
