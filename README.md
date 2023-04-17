[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![NSF-1504267](https://img.shields.io/badge/NSF-1504267-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1504267) [![NSF-1503959](https://img.shields.io/badge/NSF-1503959-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1503959)

# Cosmo3D

A 3D cosmogenic nuclide production model through complex (i.e. non-infinite half space) to estimate nuclide concentrations at depth through a surface. The original project is to determine the profile of a quarried zone based on collected samples from inside the quarried outline, as observed in the field. Input data includes: Cosmogenic Nuclide sample information, and it could include a surface profile, but initially that is hard coded for my purposes. Outputs would include figures, and minimized profiles.

A short, one paragraph description of the project goes here.

  * Please add all project Award numbers to the Badges above.  This will allow us to crawl and add these repositories to the graph, and also show NSF that we are being good citizens.
  * Please ensure that all contributors from Throughput have ORCIDs and that these are linked below in the contributors section.
  * Please consider a clear directory structure early on, and report it below in the "How to use this repository".
  * Please define which data products

## Contributors

This project is an open project, and contributions are welcome from any individual.  All contributors to this project are bound by a [code of conduct](CODE_OF_CONDUCT.md).  Please review and follow this code of conduct as part of your contribution.

  * [Brandon Graham](https://www.usgs.gov/staff-profiles/brandon-lars-graham) [![orcid](https://img.shields.io/badge/orcid-0000--0002--7197--0413-brightgreen.svg)](https://orcid.org/0000-0002-7197-0413)


### Tips for Contributing

Issues and bug reports are always welcome.  Code clean-up, and feature additions can be done either through pull requests to [project forks]() or branches.

All products of the Throughput Annotation Project are licensed under an [MIT License](LICENSE) unless otherwise noted.

## How to use this repository

The primary Matlab scripts to call directly are:  
[JakMCMCCall.m](src/JakMCMCCall.m) This is the main Markov Chain Monte Carlo simulation to do the 3D model. This call the associated data, scripts, and functions:  
    [JakMCMC.m](src/JakMCMC.m)  
    [Camp3Samples](data/Input/Camp3Samples.csv)  
    [surfaceModelSingle.m](src/surfaceModelSingle.m) This is used to take inflection point inputs for the surface model and simulated abrasion depth and create a surface model for forward modeling of cosmic ray bombardment.  
    [cosmoDistribution.m](src/cosmoDistribution.m) This generates the randomized incoming rays based on the distribution described in Gosse and Philips (2001).


It also relies on the MCMCStat model which is maintained at :[MJLaine](https://mjlaine.github.io/mcmcstat/)  

[JakCosmo3DFull.m](src/JakCosmo3DFull.m) Is designed to test and specify initial parameters, and was used to develop the MCMC model simulation and inversion.

[JakBinaryTestShielding](src/JakBinaryTestShielding.m) is written to determine if the estimated shielding factor from the measured concentrations and known exposure history matches the measured shielding in the field.

[individualShieldingDepth](src/individualShieldingDepth.m) is a stand alone version that can be used to directly solve the depth of erosion of a measured cosmogenic nuclide sample with a known exposure history.

[RoseDiagram](src/RoseDiagram.m) is used to generate the Rose Diagram figures based on the data collected in [JAK2018Camp3Flow](data/Input/JAK2018Camp3Flow.csv). 


[skyline](src/skyline.m) This is the shielding factor matlab script written by Greg Balco (2006) and is distributed under the GNU General Public License, version 2, and is used in determining shielding factor in [JakBinaryTestShielding](src/JakBinaryTestShielding.m).



### Workflow Overview

Th project uses X core information, manages it and passes our some stuff.

### System Requirements

This project is developed using Matlab 2020a.  It runs on Linux, but should also run on Windows and Mac. Older version of Matlab will likely work, however there are no guarantees.

### Data Requirements

The project pulls data from (?).

### Key Outputs

This project generates (an API, some log files, what?)

## Metrics

This project is to be evaluated using the following metrics. . .

