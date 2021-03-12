# CROSSINN - Cross-valley flow in the Inn Valley investigated by dual-Doppler lidar measurements
## Abstract 
Script collection for data processing and scan file generation of [HALO Photonics](https://halo-photonics.com/) Doppler wind lidars operated during the [CROSSINN](https://www.imk-tro.kit.edu/english/844_8306.php) (Cross-valley flow in the Inn Valley investigated by dual-Doppler lidar measurements) measurement campaign. From 14 June 2019 to 10 October 2019, two Doppler wind lidars, model Stream Line (SL_88) and Stream Line XR (SLXR_142), were operated by the [ACINN](https://www.uibk.ac.at/acinn/index.html.en) (Department of Atmospheric and Cryospheric Sciences, [University of Innsbruck](https://www.uibk.ac.at/index.html.en) at the [i-Box](https://www.uibk.ac.at/acinn/research/atmospheric-dynamics/projects/innsbruck-box-i-box.html.en) site in [Kolsass](https://acinn-data.uibk.ac.at/pages/i-box-kolsass.html) (47.305341°N, 11.62219°E). The lidars have an roating scanner head which allows for a variety of scan pattern. The Stream Line software provides profiles of radial velocity stored in txt files. 

## `crossinn_netCDF_l1.py`
**Generation of .nc files of corrected data.** 

The StreamLine software generates .hpl files (text files) containing time-distance resolved data of radial velocity, Signal-to-Noise (SNR) and attenuated backscatter (newer software also delivers spectral width). Addtionally, for each time stamp, the scanner position is given as azimuth and elevation angle. This data is converted into matrices (time x distance) and converted into netCDF (.nc) level0 (l0) data (e.g., `hpl_to_netcdf()` in marenha/doppler_wind_lidar_toolbox/2NetCDF/hpl2NetCDF.py). Level 0 data contains the same data as the original .hpl files. In `crossinn_netCDF_l1.py`, the data of these l0 .nc files is corrected and stored as level 1 (l1) netCDF data. The corrected variables are the azimuth angle and the range gate centre. The azimuth angle needs to be corrected due to a misalignment of the lidar scanner to geographical North in the field. Further, the range gates centres of the SLXR_142 show to have an offset which must be corrected. Additionally, information about the lidar location is added to the l1 netCF files. 

## `crossinn_vad2netCDF.py`
**Generation of daily .nc files of vertical profiles of horizontal wind.** 

During the campaign, conical PPI scans at an elevation angle of 70 deg were performed regularly. The data of these scans is used to derive vertical profiles of horizontal wind using the VAD algorithm. For each day one data file is created containing the retrieved wind profiles. Additionally, quicklooks of the retrieved data are created. Therefore, the modules `calc_vad.py` (marenha/doppler_wind_lidar_toolbox/VAD_retrieval) and `plot_vad.py` (marenha/doppler_wind_lidar_toolbox/quicklooks) are used. 

**To create these quicklooks, the `quiver.py` module of `matplotlib` (matplotlib/matplotlib) was modified. The modified `quiver.py` module depicts an empty barb for wind speeds between 1 and 5 kn.**

## `crossinn_write_dss.py`
**Generation of a daily scan schedule (dss) and the involved scan files for a certain scan scenario.** 

To perform different scan patterns, the StreamLine software needs a .txt files in a specific format containing the information about the scan pattern. For CROSSINN, these scan files were created with `write_scan_file.py` (marenha/doppler_wind_lidar_toolbox/SL_scan_files). The daily scan schedule created with `crossinn_write_dss.py` starts the involved scan files at a certain time of the day.

## References

Adler, B., A. Gohm, N. Kalthoff, N. Babić, U. Corsmeier, M. Lehner, M. W. Rotach, M. Haid, P. Markmann, E. Gast, G. Tsaknakis, G. Georgoussis, 2021a: CROSSINN: A Field Experiment to Study the Three-Dimensional Flow Structure in the Inn Valley, Austria. Bulletin of the American Meteorological Society, 102, E38-E60, https://doi.org/10.1175/BAMS-D-19-0283.1

Adler, B., N. Babić, N. Kalthoff, A. Wieser, 2021b: CROSSINN (Cross-valley flow in the Inn Valley investigated by dual-Doppler lidar measurements) - KITcube data sets [WLS200s]. Data set. KITopen. https://doi.org/10.5445/IR/1000127847

Adler, B., N. Babić, N. Kalthoff, A. Wieser, 2021c: CROSSINN (Cross-valley flow in the Inn Valley investigated by dual-Doppler lidar measurements) - KITcube data sets [CHM 15k, GRAW, HATPRO2, Mobotix, Photos]. Data set. KITopen. https://doi.org/10.5445/IR/1000127577

Adler, B., N. Babić, N. Kalthoff, U. Corsmeier, C. Kottmeier, C. Mallaun, 2021d: CROSSINN (Cross-valley flow in the Inn Valley investigated by dual-Doppler lidar measurements) - Aircraft data set [FDLR]. Data set. KITopen. https://doi.org/10.5445/IR/1000127862

Ladstätter, P. J., 2020: Vertical structure of the atmospheric boundary layer in the Inn Valley during CROSSINN. Master's thesis, University of Innsbruck, 114 pp. https://resolver.obvsg.at/urn:nbn:at:at-ubi:1-68131
