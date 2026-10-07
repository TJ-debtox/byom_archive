# byom_archive
This repository is the archive for the BYOM platform in Matlab (Bring Your Own Model). It archives the last version (v69a of 11 April 2025), with only a few small changes to the header texts of the m files (clarifying the licensing per file).

## Licensing

The BYOM project as a whole is licensed under the **GNU General Public License v3.0** (GPLv3) — see the [LICENSE](LICENSE) file in the root directory for the full terms. 

Because the core system integrates and depends on GPLv3-licensed components, any interactive or combined execution of this software falls under the copyleft protections of the GPLv3.

### Third-Party & Component Licenses
While the combined system is distributed under the GPLv3, individual files and sub-folders within this repository retain their original, compatible open-source licenses:

* **BYOM basics:** The main files of BYOM (apart from most files in the `/BYOM/engine*/parspace/` and `/BYOM/engine*/external/` folders) are originally written under the **MIT License**. You can find the full license text in [LICENSE-MIT.txt](LICENSE-MIT.txt) in the `/BYOM/` folder.
* **Parameter space explorer:** Most of the files in the `/BYOM/engine*/parspace/` folders are slightly modified from the code that is distributed as part of the Matlab version of openGUTS (see [http://www.openguts.info](http://www.openguts.info)). Therefore, this code is distributed under the same license as openGUTS (GPLv3). The modifications are not in the algorithm itself but only to ensure that the code operates in the general BYOM framework.
* **External Dependencies:** All third-party scripts and tools that were downloaded from other authors are isolated in the `/BYOM/engine*/external/` sub-folders ((exception is the file TSK.m, which is in the folder `BYOM/packages_GUTS/GUTS/sim_guts/` for optional LC50 calculation on simulated data sets). Each third-party file retains its original author's copyright headers including specific license information.
* Please check the individual file headers for details before reusing any of the specific components independently.

## What is BYOM?

Bring Your Own Model (BYOM) is a flexible set of Matlab scripts and functions to help you build, simulate, fit and analyse your own models. Not just DEB models or TKTD models; any model that you want (well, as long as it can be expressed in terms of ordinary differential equations, or as explicit function). Using these files requires Matlab and a working knowledge of Matlab, as there is no nice GUI. 

Downside of the enormous flexibility is a lack of user-friendliness. There is quite a learning curve to BYOM. Start by examining and running the files in the `\BYOM\examples\` directory, and start with simple models. The files in the included packages are already set up for various more complex analyses.

### End of life

As of 1 December 2026, I have retired my business DEBtox Research, and BYOM is retired as well. This archive on GitHub contains the last release: version 6.9a (11 April 2025). There will be no more updates, bug fixes or support. BYOM was last tested on Matlab R2024b.

This archive is particularly useful for those that want to redo a BYOM analysis, or those that are already familiar with it and want to update to the last version. Also, there may be someone interested in developing (parts of) BYOM further. If you are looking to start with modelling, using BYOM is probably not a good idea as it is EOL. However, I am not aware of any platform with similar features.

### History of BYOM

Development of BYOM started in 2012, as I wanted to have a platform for teaching purposes in the TKTD summercourses. I wanted students to build and fit their own models (as systems of ODE's), but I wanted them to focus on the modelling and not on the coding (BYOM takes most of the difficulties away from the user). However, I quickly found myself using it almost exclusively, as it was a perfect platform (for me) for both standard and especially for non-standard (one of a kind) modelling. 

### Supporting web page

The supporting web page for BYOM is [https://www.debtox.info/byom.html](https://www.debtox.info/byom.html). This page no longer contains the downloads for the Matlab files, but still hosts manuals, code 'walk-throughs', short description of the packages, etc. 

## Installation and folder structure

To install BYOM, unpack the downloaded ZIP file to a location of your choice. Do NOT rename the directories `/BYOM/`, `/BYOM/engine/` (and its subfolders), `/BYOM/engine_par/` (and its subfolders). You can move the complete folder `/BYOM/` (with all its subfolders) out of the main folder `/byom_archive-main/` if you like. Do NOT set the path manually to any of the BYOM folders, and please make sure that the directory where your main script is located is the active directory (if you run a script when it is not in the active directory please select <Change folder> and NOT <Add to Path>).

The reason for not changing the folder names and modifying the Matlab path is that this is done automatically in `pathdefine.m`, which should be present in any folder where you run a script file. The folder names are hard-coded in `pathdefine.m`. This function makes sure that the path is set correctly. It will also remove redundant path settings from previous runs. The changes to the path are temporary. Also note that the `pathdefine.m` of BYOM will remove any path settings made by the openGUTS Matlab version. This allows both platforms to be used in the same Matlab session with path conflicts.

The `/engine/` and `/engine_par/` directories contains all the Matlab functions for fitting the model to data, plotting results, statistical analysis etc. There should be no need to modify these files. The folders `/engine/` and `/engine_par/` are almost identical. The difference is that `/engine_par/` is set up for use with Matlab's parallel computing toolbox. A number of analyses (particularly with the parameter-space explorer) can be sped up with this toolbox (if your computer has at least 4 physical cores, preferably more). Calling `pathdefine.m` from your script with the option 1 will automatically place the `/engine_par/` directories in the path, but only if the parallel toolbox is installed in your Matlab system.

## Getting started

The first files to study are in the `/BYOM/examples/` directory. It is best not to modify these files, so you can always return to a working model. Run the `byom_...` scripts to see what they do. Check out the other functions in each folder to see if you can follow the logic. The code is heavily commented, so I hope this helps.

To create your own model, copy the files in examples to a new directory, somewhere as sub- or sub-sub-directory of the BYOM directory (do NOT rename the `/BYOM/` directory, put all new directories BELOW `/BYOM/`, and don't start any of the directory or file names with BYOM (uppercase).

Make sure that the new directory includes a copy of the functions `derivatives.m` and/or `simplefun.m` (which hold your model equations as ODE's or as explicit functions, respectively), `call_deri.m` (calls the derivatives functions with initial state values and includes the events function), and `pathdefine.m`. Modify `derivatives.m` to represent your own models. In several cases, (slight) modifications of `call_deri.m` will be needed. 

There is no need to change the functions in the `/engine/` or `/engine_par/` folders, and in general I would advise to leave them alone unless you REALLY know what you're doing. I suggest naming your script files starting with `byom_` to clearly distinguish them from the functions. Do not include spaces in script or directory names (Matlab does not like that), and don't use capitals for byom in your script names. 

### Packages

The folders `/packages_.../` contain all the BYOM packages that have been developed by me. With a package I mean a set of scripts and `derivatives.m`, `simplefun.m`, `call_deri.m` files that have already been designed for a specific purpose. This includes various DEB-based calculations, GUTS models, dose-response curve fitting, and a series of packages to redo some of the calculations in published papers or book chapters (`/packages_special_support/`). Many of the packages have a text file in their folder that explains some of the features and how to use it.

Originally, there were separate downloads on the BYOM support page, but in this archive I decided to include them all. Make sure they remain below the `/BYOM/` folder, otherwise `pathdefine.m` won't function properly. You can always remove the packages that you are not interested in, if you need the disk space.

Note that for more standard GUTS analyses, I would suggest using openGUTS, either the standalone PC software or the Matlab version, which is optimised for such analyses (and much more user-friendly). OpenGUTS can be freely downloaded from [http://www.openguts.info](http://www.openguts.info).
