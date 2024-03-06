# pafplot

I use dotplots to compare genomes and genome assemblies, but there aren't really any good tools to plot these from the command line or REPLs.

This is an early attempt to remedy this gap.

Existing tools:

- [mummerplot](https://mummer4.github.io/tutorial/tutorial.html) is the original and has a great layout algorithm, but it constantly breaks (in my hands) and looks pretty dated.
- [generateDotPlot from MashMap](https://github.com/marbl/MashMap) is essentially a quick clone of mummerplot, but with a parser for their own alignment format.
- [dgenies](https://dgenies.toulouse.inra.fr/) is a web tool to align genomes and plot the outputs. It's pretty nice, but there is no scripting interface so comparing many genomes or for general quality control this is pretty tedious.
  Even the local version has to spin up a web-server.
- [minidot](https://github.com/thackl/minidot) is quite pretty, but doesn't work for multi-sequence comparisons (e.g. a genome with multiple chromosomes) and it hasn't been updated in a long time.
- [pafr](https://dwinter.github.io/pafr/reference/dotplot.html) is probably the closest thing to what I'm after.
  It also looks a bit ugly, but that could be styled with some customisation. It also hasn't been updated in a long time.



What I hope to achieve here is a tool that:

- Can parse multiple input formats (mummer.coords/.delta, PAF, MashMap).
- Easily plot pairwise dotplots comparing genomes (two initially, eventually multiple genomes) from a command line interface and an easily hackable python interface.
  - These will include properly labelled scaffolds, facets etc.
- Simple installation with conda or pip, no more weird GNUplot backwards compatibility issues.
- Options to include marginal information (e.g. read coverage, genes, gene density, GC content).
- [Eventually] the ability to generate a static interactive HTML plot.

This covers my immediate use case.
In the future a more complete synteny/genome alignment plotting might be included, but it isn't a priority for me at the moment.


Cheers,
Darcy
