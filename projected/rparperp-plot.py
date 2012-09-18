import pyx
import numpy
pyx.text.set(mode="latex")
pyx.text.preamble("\usepackage{txfonts,color}")

basename = __file__.split('.')[0]

betalist = [0.3, 0.1, 0.01, 0.001, 1.e-4, 1.e-5, 1.e-7]
ratiolines = []
rparlines = []
rperplines = []
thparlines = []
thperplines = []
for beta in betalist:
    betastring = "%.2e" % (beta)
    betamantissa = betastring[:4]
    betaexponent = betastring[-3:]
    betatitle = r"\(\beta = %s \times 10^{%s}\)" % (betamantissa, betaexponent)

    betas, incs, thpars, thperps, rpars, rperps = \
	numpy.loadtxt("rparperp-B%.2e.dat" % (beta), skiprows=1, unpack=True)

    ratios = rperps/rpars
    incs *= 180.0/numpy.pi
    thpars *= 180.0/numpy.pi
    thperps *= 180.0/numpy.pi

    ratiolines.append(
	pyx.graph.data.values( x=incs, y = ratios,
			       title=betatitle
			       )
	)
    rparlines.append(
	pyx.graph.data.values( x=incs, y = rpars,
			       title=betatitle
			       )
	)
    rperplines.append(
	pyx.graph.data.values( x=incs, y = rperps,
			       title=None
			       )
	)
    thparlines.append(
	pyx.graph.data.values( x=incs, y = thpars,
			       title=betatitle
			       )
	)
    thperplines.append(
	pyx.graph.data.values( x=incs, y = thperps,
			       title=None
			       )
	)

g = pyx.graph.graphxy(
    width=10,
    x = pyx.graph.axis.linear(title=r'Inclination, \(i\)'),
    y = pyx.graph.axis.linear(min = 0.0, max=6.0, title=r'\(R_\perp/R_\parallel\)'),
    key = pyx.graph.key.key(pos="tr", textattrs=[pyx.trafo.scale(0.7)]),
    )

# g.plot(rparlines, [pyx.graph.style.line([pyx.style.linestyle.solid, pyx.color.gradient.Rainbow])])
# g.plot(rperplines, [pyx.graph.style.line([pyx.style.linestyle.dashed, pyx.color.gradient.Rainbow])])
g.plot(ratiolines, [pyx.graph.style.line([pyx.style.linestyle.solid, pyx.color.gradient.Rainbow])])
g.writePDFfile(basename)

g2 = pyx.graph.graphxy(
    width=10,
    x = pyx.graph.axis.linear(title=r'Inclination, \(i\)'),
    y = pyx.graph.axis.linear(min=0, max=180, title=r'\(\theta_\parallel, \theta_\perp\)'),
    key = pyx.graph.key.key(pos="tr", textattrs=[pyx.trafo.scale(0.7)]),
    )

g2.plot(thparlines, [pyx.graph.style.line([pyx.style.linestyle.solid, pyx.color.gradient.Rainbow])])
g2.plot(thperplines, [pyx.graph.style.line([pyx.style.linestyle.dashed, pyx.color.gradient.Rainbow])])
# g2.plot(thparapproxlines, [pyx.graph.style.line([pyx.style.linestyle.dashed, pyx.color.gradient.Rainbow])])

g2.writePDFfile(basename + "-thpar")
