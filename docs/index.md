# Welcome to the EVN Observation Planner (PlanObs)

The EVN Observation Planner is a tool to determine the visibility of a given astronomical source when planning very long baseline interferometry (VLBI) observations. The tool is specially written for the preparation of observations with the [European VLBI Network (EVN)](https://www.evlbi.org), but it can be used for any kind of VLBI observations than can be currently arranged (e.g. with the [Very Long Baseline Array, VLBA](https://public.nrao.edu/telescopes/vlba/); [the Australian Long Baseline Array, LBA](https://www.atnf.csiro.au/vlbi/overview/index.html); [eMERLIN](http://www.merlin.ac.uk/e-merlin/index.html); or [the global mm-VLBI array](https://www3.mpifr-bonn.mpg.de/div/vlbi/globalmm/), for example). An ad-doc VLBI array can also be quickly configured.

In addition to the determination of the source visibility by the different antennas, the EVN Observation Planner would provide an estimation of the expected rms noise level (sensitivity) reached during the planned observations, and an estimation of the resolution. The EVN Observation Planner can thus be used while [preparing an observing proposal](https://www.evlbi.org/using-evn).
Note that the EVN Observation Planner has been designed as a more complete version of the previous [EVN Calculator](http://old.evlbi.org/cgi-bin/EVNcalc.pl).


## Different ways to use PlanObs

You can use PlanObs just by directly going to [the online version hosted at JIVE](https://planobs.jive.eu), without the need of installing anything.

However, you can also install it locally at your computer and then use it in three different ways:

* `planobs-server` - Launches the same GUI as the one in the JIVE servers, so you can use it graphically.
* `planobs` - a command-line interface (CLI) tool that allows you a bit more flexibility, despite not having the graphical interface.
* By importing `vlbiplanobs` in your Python scripts.

[See getting started](getting-started.md) for a detailed explanation of all these three ways to use PlanObs.

