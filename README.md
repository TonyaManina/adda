<img src='https://raw.githubusercontent.com/wiki/adda-team/adda/img/adda.svg?sanitize=true' align='right'>

test change

ADDA (pronounced /ˈɑddɑ/) is a C software package (console application) to calculate scattering and absorption of electromagnetic waves by particles of arbitrary shape and composition using the discrete dipole approximation (DDA). The particles can be located in a homogeneous medium or near a plane substrate; emission (decay-rate) enhancement of point emitters can also be calculated. The main feature of ADDA is the ability to run on a multiprocessor system or multicore processors (parallelizing a _single_ DDA simulation). It can also employ modern GPUs to accelerate computations. ADDA is intended to be a [versatile tool](https://github.com/adda-team/adda/wiki/Features), suitable for a wide variety of applications ranging from interstellar dust and atmospheric aerosols to metallic nanoparticles and biological cells. Its applicability is limited only by [available computer resources](https://github.com/adda-team/adda/wiki/LargestSimulations).

ADDA [originated](https://github.com/adda-team/adda/wiki/EarlyHistory) at the University of Amsterdam but has then evolved into an open-source international project, used in [about 30 countries](https://github.com/adda-team/adda/wiki/Publications). We recommend to start with one of the following:
* [Getting started (wikis)](https://github.com/adda-team/adda/wiki)
* [Manual](doc/manual.pdf)
* [FAQ](https://github.com/adda-team/adda/wiki/FAQ)
* [Features](https://github.com/adda-team/adda/wiki/Features)
* [Releases/Downloads](https://github.com/adda-team/adda/releases) (source code and 64-bit Windows executables)
* [Compiling ADDA](https://github.com/adda-team/adda/wiki/CompilingADDA)

If you choose to use ADDA, please [subscribe](mailto:adda-announce+subscribe@googlegroups.com) to [announcement mailing list](http://groups.google.com/group/adda-announce); "registered" users of ADDA will be notified when updates to the code are made. If you publish results obtained using ADDA, you should acknowledge the source of the code. The general reference is - Yurkin M.A. and Hoekstra A.G. The discrete-dipole-approximation code ADDA: capabilities and known limitations, [_J. Quant. Spectrosc. Radiat. Transfer_ **112**, 2234–2247](http://doi.org/10.1016/j.jqsrt.2011.01.031) (2011).
Please also look at [a list of more specific references](https://github.com/adda-team/adda/wiki/References).

We encourage users to provide feedback in any possible way. If you have any questions, please [write](mailto:adda-discuss@googlegroups.com) to the [discussion group](http://groups.google.com/group/adda-discuss). If you have any suggestions or a bug report, submit it directly to the [issue tracker](https://github.com/adda-team/adda/issues), taking advantage of the open development process. Please also provide feedback on the existing issues - you may upvote issues that you find important (GitHub account is required) and/or provide some meaningful comments.

### Directory structure

* `.github/` - GitHub-related settings
* `devtools/` - tools used by ADDA developers
* `doc/` – documentation
* `input/` – default input files
* `misc/` – additional tools, relatively independent from the main part of ADDA
* `sample/` – sample output and other files
* `src/` - source files and makefiles
* `tests/2exec` - extensive test of two ADDA executables (in terms of output) against each other
* `win64/` - executables for 64-bit (x86\_64) Windows.
* `LICENSE` - GNU General Public License v3.0
