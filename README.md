# TOPAS-nBio (Version 4.1.0)
This is the TOPAS-nBio extension repository, a Monte Carlo simulation framework for (sub-) cellular radiobiology.

> [!WARNING]
> TOPAS-nBio 4.1.0 is not compatible with Geant4 11.4.2. In this combination,
> Geant4-DNA chemistry may produce no molecules, so the TOPAS-nBio build is
> intentionally stopped with a compatibility error. Use Geant4 11.3.2 when
> building TOPAS-nBio. This restriction concerns the TOPAS-nBio and Geant4
> combination; it does not indicate a known incompatibility with the TOPAS
> 4.3.0 candidate release itself.

TOPAS-nBio is described here: https://topas-nbio.readthedocs.io/. 
This page includes a class documentation and the license.

TOPAS-nBio is an extension of OpenTOPAS (TOol for PArticle Simulations), which can be obtained from https://OpenTOPAS.github.io. The TOPAS documentation can be found at https://opentopas.readthedocs.io/. 

The TOPAS-nBio package was described in [Schuemann et al., Radiation Research, 2019, 191(2), p.125](http://www.rrjournal.org/doi/10.1667/RR15226.1). This reference should be cited for all work using the TOPAS-nBio package.

We encourage contribution of user-developed extensions that fit within the TOPAS-nBio scheme. If you would like to contribute code, please contact the developers.

---

## General information

> [!NOTE]
> The main difference in installing for Debian/Ubuntu and MacOS is the installation paths. 

> [!TIP]
> Users can install multiple extensions by either:
> - Copy all the Extension directories into a common global directory, or
> - Using the `\;` between directory paths e.g., (and see below):
>
>      `-DTOPAS_EXTENSIONS_DIR=/MyProject1Extensions/\;/MyProject2Extensions/`

### Debian

Follow [installation for Debian/Ubuntu](https://opentopas.readthedocs.io/en/latest/getting-started/Debian.html) procedures after and including **step 8.4**. Then:

    $ cd $HOME/Applications/

Clone the repository:

    $ git clone https://github.com/topas-nbio/TOPAS-nBio.git

Return to the OpenTopas-build directory:

    $ cd $HOME/Applications/TOPAS/OpenTOPAS-build

Repeat `cmake` command by including the `TOPAS_EXTENSIONS_DIR` variable as follows (don't miss the `-D`):

    $ cmake ../OpenTOPAS -DCMAKE_INSTALL_PREFIX=../OpenTOPAS-install -DTOPAS_EXTENSIONS_DIR=$HOME/Applications/TOPAS-nBio

Build:

    $ make -j8 install

### MacOS

Follow [installation for MacOS](https://opentopas.readthedocs.io/en/latest/getting-started/MacOS.html) procedures after and including **step 7.4**. Then:

    $ cd /Applications

Clone the repository:

    $ git clone https://github.com/topas-nbio/TOPAS-nBio.git

Return to the OpenTopas-build directory:

    $ cd /Applications/TOPAS/OpenTOPAS-build

Repeat `cmake` command by including the `TOPAS_EXTENSIONS_DIR` variable as follows (don't miss the `-D`):

    $ cmake ../OpenTOPAS -DCMAKE_INSTALL_PREFIX=../OpenTOPAS-install -DTOPAS_EXTENSIONS_DIR=/Applications/TOPAS-nBio

Build:

    $ make -j8 install
    
---

Run the demos. For some demos, a pause before quit is enabled, then, write `exit` at the terminal prompt.

   Linux:
        
        $ source rundemos.csh

   Mac:
        
        $ source rundemos.csh

---
