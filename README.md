# HiRes2LowRes_Converter.py

## Description

`HiRes2LowRes_Converter.py` is used to convert high-resolution MS data to low-resolution MS data (unit resolution) for further alignment using Guineu software (v1.0.3; https://doi.org/10.1021/ac103308x) as Guineu can handle only unit-resolution MS data. In addition, several other important amends are made by the script (see below), in a way to fit the file format accepted by Guineu (v1.0.3).

## What does the script do

The script takes tab-separated `TXT` files as input and does the following:
1.	Checks the input filenames and names output files in a format `“YYMMDD_LR_InitialFilename_Guineu”`, where `“YYMMDD”` is today’s date and `“LR”` stands for `Low Resolution”`. In case the initial filename already contains date in a format `"YYMMDD"`, it will be replaced with today’s date. Also, a suffix `“_Result”` (ChromaTOF export artifact) will be removed in case it was included in the input filenames. 
2.	Checks whether `"R.T. (s)"` column is present in the input file and, if yes, proceeds with further processing. Otherwise, the script returns error. 
3.	Deletes features (rows) for which there are empty cells in `"Area"` column.
4.	Removes all features (rows), for which there is a value `“Bleed”`, `“Toluene”`, or `“DCM”` in `“Classifications”` column. The classification is done in ChromaTOF using spectral peak filter and region functionalities; see ChromaTOF manual for more information.
5.	Removes all features (rows) where `"Name"` column contains unwanted keywords, e.g., `“silo”` for any `“siloxane”` or `“syloxy”` compounds or `“bora”` for any `“borane”` or `“borate”` features. 
6.	Removes all features (rows) where `"Name"` column contains unwanted keywords and does not contain specific m/z, e.g., compounds with spectra similar to toluene's (m/z 91 or 92) if their spectra do not contain other higher m/z values, like molecular ions.
7.	Rounds RT1 values of the column `“R.T. (s)”`. 
8.	Renames “non-hit” entries (unknown features; no match in the NIST library) from `“Peak #”` to `“Peak_RT1 (rounded)_RT2”`, where RT1 and RT2 are retention times in the 1st and 2nd dimensions, respectively.
9.	Populates empty `“Similarity”` column cells for “non-hits” with a value of `800` as it is done by the older version of Guineu (v0.9) to enable simultaneous processing of both “hits” and “non-hits” in Guineu (v1.0.3); otherwise, “non-hits” are skipped during the inter-sample data alignment in Guineu (v1.0.3) (`“Score Alignment”` feature).
10.	Rounds `“Spectrum”` column values (m/z & intensity). In case duplicate m/z's appear in a given spectrum after rounding, only one m/z is stored and the respective intensity values are added up; otherwise, one intensity value is assigned to one m/z value.
11.	Rescales `“Spectrum”` column values to max intensity of `1000`. 
12.	Removes all fragments with intensities <1 from `“Spectrum”` column values.
13.	Removes all fragments with intensities <3% from `“Spectrum”` column values.

## Prerequisites

Before using the script, several applications/tools have to be installed:
1.	Visual Studio Code; https://code.visualstudio.com/download.
2.	Python 3; https://www.python.org/downloads/windows/.
3.	Python Extension in Visual Studio Code > Extensions (`Ctrl + Shift + X`) > Search “python” > Press `Install`.

Then, the packages must be installed as follows:
Visual Studio Code > Terminal > New Terminal > In terminal, type `pip install package_name`, where `package_name` is a desired package name, e.g., `numpy` > Press `Enter`.

## How to use the script

To use the script, the following steps must be executed:
1.	Specify required information:
- Specify keywords to be filtered out in **lines 198-210** as, for example: `keywords = ["sila", "silic", "silo", "TMS"]`
- Specify filters as a list of tuples (`keyword, specific integer m/z in spectrum`) in **lines 212-340** as, for example: `filters = [("(Benzyloxy)(methyl)amine", 105), ("(Isopropoxymethy)lbenzene", 107)]`
- Specify filters for `"Peak"` (unknown features) and bleed m/z as a list of tuples (`keyword, specific integer m/z in spectrum`) in **lines 342-358** as, for example: `bleed_filters = [("Peak", 73), ("Peak", 147), ("Peak", 207)]`. **NB!** All three filters can be left blank in case no filtration is needed as follows: `keywords = []`, `filters = []`, `bleed_filters = []`
- Save the updated script by pressing `Ctrl + S`.

2.	To run the script:
- Right mouse click anywhere in Visual Studio Code script file > Run Python > Run Python File in Terminal or press `play` button in the top-right corner.
- Choose the files for processing in the new pop-up window and press `Open`.

## Notes and recommendations

The input files must contain at least the following columns to be processed: 
`"Name" "R.T. (s)" "Similarity" "Area" "Spectrum" "Classifications"`

**NB!** `"Spectrum"` values must be in LECO ChromaTOF format: `39:4500 52:220 67:9999`.

## License
[![MIT License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit)

Intended for academic and research use.
