# byom_archive
This repository is the archive for the BYOM platform in Matlab (Bring Your Own Model). It archives the last version (v69a of 11 April 2025), with only a few small changes to the header text (clarifying the licensing per file).

## Licensing

The BYOM project as a whole is licensed under the **GNU General Public License v3.0** (GPLv3) — see the [LICENSE](LICENSE) file in the root directory for the full terms. 

Because the core system integrates and depends on GPLv3-licensed components, any interactive or combined execution of this software falls under the copyleft protections of the GPLv3.

### Third-Party & Component Licenses
While the combined system is distributed under the GPLv3, individual files and sub-folders within this repository retain their original, compatible open-source licenses:

* **BYOM basics:** The main files of BYOM (apart from most files in the `/engine_*/parspace/` and `/engine_*/external/` folders) are originally written under the **MIT License**. You can find the full license text in [LICENSE-MIT.txt](LICENSE-MIT.txt).
* **Parameter space explorer:** Most of the files in the `/engine_*/parspace/` folders are slightly modified from the code that is distributed as part of the Matlab version of openGUTS (see [http://www.openguts.info](http://www.openguts.info)). Therefore, this code is distributed under the same license as openGUTS (GPLv3). The modifications are not in the algorithm itself but only to ensure that the code operates in the general BYOM framework.
* **External Dependencies:** All third-party scripts and tools that were downloaded from other authors are isolated in the `/engine/external/` and `/engine_par/external/` sub-folders. 
  * Each file inside `external/` retains its original author's copyright headers including specific license information.
  * Please check the individual file headers for details before reusing those specific components independently.
