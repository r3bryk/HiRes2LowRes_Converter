# importing functions
import csv
import sys
import re
import tkinter.filedialog
import os, glob, codecs, math
from os.path import dirname, basename, join
import time
from datetime import timedelta
import pandas as pd
from pandas import DataFrame
import numpy as np
import warnings
warnings.simplefilter(action="ignore", category=pd.errors.SettingWithCopyWarning)

# defyning the main function
def HiRes2LowRes(filename):
    print('')
    print('Sample start time:', time.strftime("%Y-%m-%d %H:%M:%S"))
    print('')
    print(150 * "*")
    # naming the output file in a format "YYMMDD_LR_InputFileName" & saving it to the input file directory
    # checking whether the initial filename contains date in a format "YYMMDD"
    if re.match('\d{6}', basename(filename)):
        outFilename = join(dirname(filename), time.strftime("%y%m%d") + re.sub('\d{6}', '_LR', basename(filename.replace('_Result', '').replace('.txt', '_Guineu.txt')))) # naming the output file for filenames that contain date
    else:
        outFilename = join(dirname(filename), time.strftime("%y%m%d") + "_LR_" + basename(filename.replace('_Result', '').replace('.txt', '_Guineu.txt'))) # naming the output file for filenames that do NOT contain date

    with open(filename, 'r', newline = '') as inFile: # opening the input file
        with open(outFilename, 'w', newline = '') as outFile: # opening the output file
            reader_DF = pd.read_csv(inFile, sep = '\t') # reading & storing the input file as DataFrame (DF)

            try:
                if 'R.T. (s)' in reader_DF.columns: # checking for a specific column name R.T. (s)
                    
                    reader_DF['Area'].replace('', np.nan, inplace = True) # replacing empty values in Area column with NaN
                    reader_DF.dropna(subset = ['Area'], inplace = True) # deleting all rows that contain NaN in Area column
                    init_DF = reader_DF

                    init_rows = len(reader_DF)
                    print()
                    print(f"Initial number of rows: {init_rows}.")
                    
                    # filtering away features (dropping rows) classified as 'Bleed' by ChromaTOF filters & regions
                    reader_DF = reader_DF[~reader_DF['Classifications'].str.contains('Bleed', case = False, na = False)]
                    reader_DF = reader_DF[~reader_DF['Classifications'].str.contains('Toluene', case = False, na = False)]
                    reader_DF = reader_DF[~reader_DF['Classifications'].str.contains('DCM', case = False, na = False)]
                    reader_DF.reset_index(drop = True, inplace = True) # resetting index

                    class_fltrd_rows = len(reader_DF)
                    print(f"Number of filtered rows: {init_rows - class_fltrd_rows}.")
                    print(f"Number of rows after class filtering: {class_fltrd_rows}. \n")
                    print("Done removing 'Bleed', 'Toluene' & 'DCM'", "\n")

                    # a function to escape special regex characters
                    def escape_regex(str):
                        special_chars = r"([\\.+*?^=!:${}()|\[\]\/&])"
                        return re.sub(special_chars, r"\\\1", str)

                    # a function to filter out rows where "Name" contains unwanted values
                    def filter_names(reader_DF, keywords):
                        # pre-escape keywords for regex safety
                        escaped_keywords = [re.compile("(?i)" + escape_regex(k)) for k in keywords]

                        # check if 'Name' column exists
                        if 'Name' in reader_DF.columns:
                            # replace missing values in 'Name' with "NA"
                            reader_DF['Name'] = reader_DF['Name'].fillna('NA').str.strip()

                            # print the original number of rows for debugging
                            original_rows = len(reader_DF)
                            print(f"Original number of rows before name filtering: {original_rows}.")

                            # find rows matching keywords for debugging
                            matched_rows = reader_DF[reader_DF['Name'].apply(lambda x: any(regex.search(x) for regex in escaped_keywords))]
                            print(f"Rows matching keywords: {len(matched_rows)}.")
                            print()
                            print("Matched rows preview:")
                            print(150*"-")
                            for name in matched_rows['Name']:
                                print(name)
                                print(150*"-")
                            print()

                            # delete features (rows) where 'Name' contains bleed compound name parts
                            reader_DF = reader_DF[~reader_DF['Name'].apply(lambda x: any(regex.search(x) for regex in escaped_keywords))]

                            # print the number of rows after filtering
                            filtered_rows = len(reader_DF)
                            print(f"Original number of rows before name filtering: {original_rows}.")
                            print(f"Rows matching keywords: {len(matched_rows)}.")
                            print(f"Number of filtered rows: {original_rows - filtered_rows}.")
                            print(f"Number of rows after name filtering: {filtered_rows}.", "\n")

                            return reader_DF
                        else:
                            print(f"Warning: 'Name' column not found. Skipping filtration.")
                            return reader_DF

                    # a function to filter out rows after alignment; similar to filter_names_if_missing_mz
                    def filter_names_if_missing_mz(reader_DF, filters):
                        # pre-escape keywords for regex safety
                        escaped_filters = [(re.compile("(?i)" + escape_regex(f[0])), f[1]) for f in filters]

                        # check if 'Name' and 'Spectrum' columns exist
                        if 'Name' in reader_DF.columns and 'Spectrum' in reader_DF.columns:
                            # replace missing values in 'Name' with "NA"
                            reader_DF['Name'] = reader_DF['Name'].fillna('NA').str.strip()

                            # print the original number of rows for debugging
                            original_rows = len(reader_DF)
                            print()
                            print(f"Input number of rows before name & m/z filtering: {original_rows}.")

                            # find rows matching keywords for debugging
                            matched_rows = reader_DF[reader_DF['Name'].apply(lambda x: any(regex[0].search(x) for regex in escaped_filters))]
                            print(f"Rows matching keywords: {len(matched_rows)}.", "\n")
                            
                            # print string (names) and integer (m/z's) keywords found in the input files
                            for index, row in reader_DF.iterrows():
                                for regex in escaped_filters:
                                    if regex[0].search(row['Name']):
                                        mz_values = [int(float(mz)) for mz in re.findall(r'\d+', row['Spectrum'])]
                                        print(f"Found keyword: {regex[0].pattern} in 'Name': {row['Name']}")
                                        if regex[1] in mz_values:
                                            print(f"Target m/z {regex[1]} found in 'Spectrum'")
                                        print(150*"-")

                            # delete features (rows) where 'Name' contains bleed compound name parts & DO NOT contain specific m/z's
                            def filter_condition(row):
                                return not any(regex[0].search(row['Name']) and regex[1] not in [int(float(mz)) for mz in re.findall(r'\d+', row['Spectrum'])] for regex in escaped_filters)

                            reader_DF = reader_DF[reader_DF.apply(filter_condition, axis=1)]

                            # print the number of rows after filtering
                            filtered_rows = len(reader_DF)
                            print()
                            print(f"Input number of rows after name filtering: {original_rows}.")
                            print(f"Rows matching keywords: {len(matched_rows)}.")        
                            print(f"Number of filtered rows: {original_rows - filtered_rows}.")
                            print(f"Number of rows after name & m/z filtering: {filtered_rows}.", "\n")
                            print(f"Final number of rows: {filtered_rows}.")

                            return reader_DF
                        else:
                            print(f"Warning: 'Name' or 'Spectrum' column not found. Skipping filtration.")
                            return reader_DF
                    
                    # a function to filter out rows after alignment; similar to filter_names_if_missing_mz
                    def filter_peak_names_with_bleed_mz(reader_DF, bleed_filters):
                        # pre-escape keywords for regex safety
                        escaped_filters = [(re.compile("(?i)" + escape_regex(f[0])), f[1]) for f in bleed_filters]

                        # check if 'Name' and 'Spectrum' columns exist
                        if 'Name' in reader_DF.columns and 'Spectrum' in reader_DF.columns:
                            # replace missing values in 'Name' with "NA"
                            reader_DF['Name'] = reader_DF['Name'].fillna('NA').str.strip()

                            # print the original number of rows for debugging
                            original_rows = len(reader_DF)
                            print()
                            print(f"Input number of rows before bleed filtering: {original_rows}.")

                            # find rows matching keywords for debugging
                            matched_rows = reader_DF[reader_DF['Name'].apply(lambda x: any(regex[0].search(x) for regex in escaped_filters))]
                            print(f"Rows matching keywords: {len(matched_rows)}.", "\n")
                            print(150 * "-")

                            # delete features (rows) where 'Name' contains "Peak" and the highest intensity m/z is one of the specified values
                            def filter_condition(row):
                                pairs = row['Spectrum'].split(' ')
                                intensities = [float(pair.split(':')[1]) for pair in pairs]
                                max_intensity_index = intensities.index(max(intensities))
                                max_mz = int(float(pairs[max_intensity_index].split(':')[0]))
                                condition = any(regex[0].search(row['Name']) and max_mz == regex[1] for regex in escaped_filters)
                                if condition:
                                    print(f"Removing row with 'Name': {row['Name']} and m/z: {max_mz}")
                                return not condition
                            print(150 * "-")
                                
                            reader_DF = reader_DF[reader_DF.apply(filter_condition, axis = 1)]

                            # print the number of rows after filtering
                            filtered_rows = len(reader_DF)
                            print()
                            print(f"Input number of rows for bleed filtering: {original_rows}.")
                            print(f"Rows matching keywords: {len(matched_rows)}.")        
                            print(f"Number of filtered rows: {original_rows - filtered_rows}.")
                            print(f"Number of rows after bleed filtering: {filtered_rows}.", "\n")
                            print(f"Final number of rows: {filtered_rows}.")
                            print(150 * "-")

                            return reader_DF
                        else:
                            print(f"Warning: 'Name' or 'Spectrum' column not found. Skipping filtration.")
                            return reader_DF

                    # specify the keywords to filter out
                    keywords = [
                        "1H-Tetrazol-5-amine", "1H-Tetrazole, 1-methyl-", "Acetic acid, mercapto-", 
                        "1H-Pyrrole-2-carbonitrile", "5-Diazo-1,3-cyclopentadiene", "1H-Pyrrole-3-carbonitrile", 
                        "1H-1,2,3-Triazole-4-", "1H-1,2,4-Triazole, 3-", "1-Propanamine, 3-dibenzo", 
                        "2,4,6-Cycloheptatrien", "2-Picoline, 6-nitro-", "4-Benzyloxy-3-hydroxy", 
                        "9-[(S)-2-(Hydroxymethyl)", "Arsan", "arsan", "Arsen", "arsen", 
                        "Benzyl (1,2,3-thiadiazol", "Benzyl-4,9,15-trioxa", "Bora", "bora", 
                        "Borin", "borin", "boro", "Boro", "boryl", "Chromium", "Cobalt", 
                        "Iron,", "Mercury", "sila", "Sila", "Silic", "silic", "Silo", "silo", 
                        "silver", "silyl", "TBDMS", "TMS", "triazolo[", "Tricyclo[3.2.2.0(2,4)]", 
                        "Zinc", "Zirconium"
                    ]

                    # specify the filters as a list of tuples (keyword, spectrum_number (specific integer m/z in spectrum))
                    filters = [
                        ("(4aS,4bS,10aS)-1,1,4a-Trimethyl-7-(propan-2-ylidene)-1,2,3,4,4a,4b,5,6,7,9,10,10a-dodecahydrophenanthrene", 272),
                        ("(Benzyloxy)(methyl)amine", 105),
                        ("(Isopropoxymethy)lbenzene", 107),
                        ("(R)-9-[(S)-2-(Hydroxymethyl)pyrrolidin-1-yl]-3-methyl-3,4-dihydro-2H-benzo[b][1,4,5]oxathiazepine 1,1-dioxide", 312),
                        ("(S)-9-[(S)-2-(Hydroxymethyl)pyrrolidin-1-yl]-3-methyl-3,4-dihydro-2H-benzo[b][1,4,5]oxathiazepine 1,1-dioxide", 311),
                        ("1-Benzyl-2-(trifluoromethyl)aziridine", 201),
                        ("1-Benzyl-3-hydroxypyridinium hydroxide", 185),
                        ("1-Benzylcyclopentanol-1", 129),
                        ("1H-[1,2,3]Triazole-4-carboxylic acid, 5-acetylamino-1-benzyl-, phenylamide", 335),
                        ("1-Phenyl-2-propanol", 136),
                        ("2(3H)-Benzofuranone, 6-ethenylhexahydro-6-methyl-3-methylene-7-(1-methylethenyl)-, [3aS-(3aa,6a,7ß,7aß)]-", 232),
                        ("2-(Benzyloxy)-1-chloro-3-fluorobenzene", 117),
                        ("2-(Benzyloxy)-4-methoxybenzaldehyde", 242),
                        ("2-Benzyl-3-isopropyl-cyclopentanone", 216),
                        ("2-Benzyloxyethylamine", 105),
                        ("2-Butanol, 3-benzyloxy-", 135),
                        ("3-(Iodomethyl)pyridine", 219),
                        ("3-Benzyl-4-chloro-1,2,3-triazole 1-oxide", 130),
                        ("3-Benzylsulfanyl-3-fluoro-2-trifluoromethyl-acrylonitrile", 258),
                        ("3H-Pyrazole, 5-ethynyl-3,3-dimethyl-", 120),
                        ("3-Picoline, 2-nitro-", 108),
                        ("4-(Benzyloxy)-3-fluorophenol, trifluoroacetate", 314),
                        ("4-(Benzyloxy)pyridine 1-oxide", 201),
                        ("4-Azido-2-phenylmethanesulfinyl-benzonitrile", 115),
                        ("4-Azido-2-phenylmethanesulfonyl-benzonitrile", 234),
                        ("4-Benzyloxy-2-methyl-2-buten-1-ol", 108),
                        ("4H-1,2,4-triazol-3-ol, 5-[(phenylmethyl)thio]-", 207),
                        ("5-Benzyloxy-2-nitrotoluene", 243),
                        ("5-Diazo-1,3-cyclopentadiene", 63),
                        ("7,8-Diazabicyclo[4.2.2]deca-2,4,7,9-tetraen-7-oxide", 118),
                        ("7-Chloro-2,3-dihydro-3-(4-N,N-dimethylaminobenzylidene)-5-phenyl-1H-1,4-benzodiazepin-2-one", 159),
                        ("Benzaldehyde, 4-methoxy-3-(phenylmethoxy)-", 242),
                        ("Benzene, (1,1-dimethylnonyl)-", 232),
                        ("Benzene, (1-azidoethyl)-", 147),
                        ("Benzene, (1-butylheptyl)-", 147),
                        ("Benzene, (1-butylnonyl)-", 147),
                        ("Benzene, (1-ethyldecyl)-", 246),
                        ("Benzene, (1-ethylundecyl)-", 260),
                        ("Benzene, (1-methyldecyl)-", 232),
                        ("Benzene, (1-methyldodecyl)-", 260),
                        ("Benzene, (1-methylundecyl)-", 246),
                        ("Benzene, (1-pentylheptyl)-", 246),
                        ("Benzene, (1-pentyloctyl)-", 260),
                        ("Benzene, (1-propyldecyl)-", 133),
                        ("Benzene, (1-propylnonyl)-", 133),
                        ("Benzene, (2,2-dichloroethyl)-", 174),
                        ("Benzene, (2-chloroethyl)-", 140),
                        ("Benzene, (2-chloropropyl)-", 154),
                        ("Benzene, (2-cyclohexylethyl)-", 188),
                        ("Benzene, (2-methylpropyl)-", 134),
                        ("Benzene, (3-methylpentyl)-", 162),
                        ("Benzene, (bromomethyl)-", 170),
                        ("Benzene, (butoxymethyl)-", 107),
                        ("Benzene, (ethoxymethyl)-", 135),
                        ("Benzene, (iodomethyl)-", 127),
                        ("Benzene, (phenoxymethyl)-", 184),
                        ("Benzene, (propoxymethyl)-", 107),
                        ("Benzene, [(2-propenyloxy)methyl]-", 107),
                        ("Benzene, [(methylsulfinyl)methyl]-", 154),
                        ("Benzene, [(methylsulfonyl)methyl]-", 170),
                        ("Benzene, 1,1'-(1,1,2,2-tetramethyl-1,2-ethanediyl)bis-", 119),
                        ("Benzene, 1,1'-[oxybis(methylene)]bis-", 107),
                        ("Benzene, 1-methyl-2-nitroso-", 121),
                        ("Benzene, 1-methyl-3-(1-methylethenyl)-", 132),
                        ("Benzene, n-butyl-", 134),
                        ("Benzeneacetaldehyde", 120),
                        ("Benzeneacetamide", 135),
                        ("Benzeneacetamide, N", 149),
                        ("Benzeneacetic acid", 136),
                        ("Benzeneacetic acid 1-methylethyl ester", 176),
                        ("Benzeneacetic acid, 2-propenyl ester", 103),
                        ("Benzeneethanol", 178),
                        ("Benzenemethanesulfonamide", 107),
                        ("Benzenemethanesulfonyl chloride", 126),
                        ("Benzenemethanethiol", 124),
                        ("Benzenepropanenitrile, a-phenyl-", 107),
                        ("Benzenesulfonamide, 4-methyl-", 171),
                        ("Benzonitrile, m-phenethyl-", 201),
                        ("Benzyl (1,2,3-thiadiazol-4-y)carbamate", 108),
                        ("Benzyl 2-chloroethyl sulfone", 218),
                        ("Benzyl 4-nitrophenyl carbonate", 139),
                        ("Benzyl butyl phthalate", 149),
                        ("Benzyl chloride", 126),
                        ("Benzyl chloroformate", 170),
                        ("Benzyl isopentyl ether", 107),
                        ("Benzyl lactate", 180),
                        ("Benzyl methyl disulfide", 170),
                        ("Benzyl methyl ketone", 134),
                        ("Benzyl N-[4-(4-cyano-3-fluorophenyl)phenyl]carbamate, TFA", 238),
                        ("Benzylcyclopentane", 160),
                        ("Bibenzyl", 182),
                        ("Bicyclo[2.2.2]oct-7-en-2-one, 5-methylene-", 134),
                        ("Bicyclo[3.1.1]hept-2-ene, 3,6,6-trimethyl-", 121),
                        ("Butane, 1-(benzyloxy)-2-[(benzyloxy)methyl]-", 193),
                        ("Cycloheptatrienylium, iodide", 78),
                        ("Decane, 1-chloro-", 105),
                        ("Dispiro[cyclopropane-1,3'-tricyclo[5.2.1.0(2,6)]decane-10',1''-cyclopropane]-4',8'-diene", 115),
                        ("Dodecane", 170),
                        ("Dodecane, 1-chloro-", 105),
                        ("Dodecane, 2,6,11-trimethyl-", 113),
                        ("Ethane, hexachloro-", 201),
                        ("Eucalyptol", 154),
                        ("Hydrazinecarbothioamide", 60),
                        ("Hydroxylamine, O-(phenylmethyl)-", 105),
                        ("Isophorone", 138),
                        ("MGK-264", 164),
                        ("N-(Phenylacetyl)glycine", 193),
                        ("N-Benzyloxy-2-carbomethoxyaziridine", 105),
                        ("N-Hydroxymethyl-2-phenylacetamide", 165),
                        ("Pentadecane", 212),
                        ("Pentalene, 1,2,4,5,6,6a-hexahydro-2-methylene-", 120),
                        ("Phenylacetamide", 118),
                        ("Phenylacetamide, N-propyl-", 177),
                        ("Phenylethyl Alcohol", 122),
                        ("Phosphine oxide, bis(pentamethylphenyl)-", 342),
                        ("Phthalan", 120),
                        ("Pyridine, 2-methyl-, 1-oxide", 109),
                        ("Pyridine, 4-methyl-2-nitro-", 138),
                        ("Sabinyl, 2-methylbutanoate", 119),
                        ("Spiro[4.4]non-3-en-2-one, 4-methyl-3-(1H-tetrazol-5-yl)-1-oxa-", 123),
                        ("ß-Myrcene", 136),
                        ("Sydnone, 3-(phenylmethyl)-", 176),
                        ("tert-Nonylphenol, Ac derivative", 177),
                        ("Thiocyanic acid", 149),
                        ("Thiocyanic acid, phenylmethyl ester", 149),
                        ("Tricyclo[3.2.2.0(2,4)]non-8-ene-6,6,7,7-tetracarbonitrile", 128)
                    ]

                    # specify the filters for "Peak" and bleed m/z as a list of tuples (keyword, spectrum_number (specific integer m/z in spectrum))
                    bleed_filters = [
                        ("Peak", 73),
                        ("Peak", 147),
                        ("Peak", 207),
                        ("Peak", 267),
                        ("Peak", 281),
                        ("Peak", 341),
                        ("Peak", 355),
                        ("Peak", 429),
                        ("Peak", 479),
                        ("Peak", 59), ("Peak", 78), ("Peak", 135), ("Peak", 156), ("Peak", 193), ("Peak", 195), ("Peak", 209), ("Peak", 215),
                        ("Peak", 221), ("Peak", 251), ("Peak", 253), ("Peak", 255), ("Peak", 269), ("Peak", 327), ("Peak", 329), ("Peak", 331),
                        ("Peak", 333), ("Peak", 339), ("Peak", 343), ("Peak", 377), ("Peak", 401), ("Peak", 403), ("Peak", 405), ("Peak", 415),
                        ("Peak", 417), ("Peak", 439), ("Peak", 451), ("Peak", 475), ("Peak", 477), ("Peak", 489), ("Peak", 503), ("Peak", 549),
                        ("Peak", 553), ("Peak", 563), ("Peak", 623), ("Peak", 91), ("Peak", 92)
                    ]
                    
                    # apply both filters
                    reader_DF = filter_names(reader_DF, keywords)
                    reader_DF2 = filter_names_if_missing_mz(reader_DF, filters)
                    print("DF after name & mz filtering: \n", reader_DF2)
                    reader_DF = filter_peak_names_with_bleed_mz(reader_DF2, bleed_filters)
                    print(150*"-", "\n")
                    reader_DF.reset_index(drop = True, inplace = True) # resetting index
                    
                    print(f"Number of rows after class filtering: {class_fltrd_rows}.")
                    print(f"Initial number of rows: {init_rows}. \n")
                    print("Done filtering by name and m/z", "\n")
                    print(150*"*", "\n")
                    print("Input DF: ", "\n") 
                    print(init_DF, "\n")
                    print(150*"*", "\n")
                    print("Filtered DF: ", "\n")
                    print(reader_DF, "\n")
                    print(150*"*", "\n")
                    
                    reader_DF[['RT1', 'RT2']] = reader_DF['R.T. (s)'].str.split(', ', expand = True) # splitting R.T. (s) column into RT1 & RT2
                    rt1 = reader_DF['RT1'] # storing RT1
                    rt2 = reader_DF['RT2'] # storing RT2
                    rt1_round = [] # creating an empty list for storing rounded RT1 values
                    rt2_float = [] # creating an empty list for storing float RT2 values
                    for i in range(len(rt1)):
                        rt1_round.append(round(float(rt1[i]))) # rounding & appending RT1
                        rt2_float.append(float(rt2[i])) # appending RT2
                    rt1_round = pd.DataFrame(rt1_round).set_axis(['1st Dimension Time (s)'], axis = 'columns') # converting rounded RT1 into DF & renaming the column
                    rt2_float = pd.DataFrame(rt2_float).set_axis(['2nd Dimension Time (s)'], axis = 'columns') # converting RT2 into DF & renaming the column
                    frames = [rt1_round, rt2_float] # merging rounded RT1 & RT2 into list
                    rts = pd.concat(frames, join = "outer", axis = 1).astype(str) # concatenating rounded RT1 & RT2 columns
                    rts['R.T. (s)'] = rts[['1st Dimension Time (s)', '2nd Dimension Time (s)']].agg(', '.join, axis = 1) # merging rounded RT1 & RT2 values into one column separated with commas, i.e. 123, 0.123, and renaming with R.T. (s)
                    rts = rts['R.T. (s)'] # storing just R.T. (s) column; by default, initial RT1 & RT2 columns are also present
                    print("Done RT1 and RT2 splitting and rounding", "\n")

                    for a in range(len(reader_DF.loc[:, ('Name')])): # replacing 'non-hits' Name values (Peak #) with Peak_RT1_RT2
                        if str('Peak') in str(reader_DF.loc[a, ('Name')]):
                            reader_DF.loc[a, ('Name')] = str('Peak_' + str(round(float(rt1.loc[a]))) + '_' + str(rt2.loc[a]))
                    print("Done renaming non-hits", "\n")

                    spectrum = reader_DF['Spectrum'] # storing Spectrum values
                    name = reader_DF['Name'] # storing Name values
                    other = reader_DF # storing initial table as new DF & removing Name, RTs, and Spectrum columns
                    other.drop(['Name', 'R.T. (s)', 'RT1', 'RT2', 'Spectrum'], axis = 1, inplace = True)
                    print("Done dropping unused columns", "\n")

                    other['Base Mass'] = other['Base Mass'].round(0).apply(int) # rounding Base Mass values
                    print("Done rounding 'Base Mass' values", "\n")

                    other['Similarity'] = other['Similarity'].fillna(800).round(0).apply(int) # replacing 'non-hits' Similarity values (N/A) with 800 (as it is done by Guineu 0.9) & rounding
                    print("Done replacing non-hits 'Similarity' values", "\n")
                    
                    # rounding values in Spectra (mz & intensity values)
                    # a function to remove pairs with intensity values containing '0.'
                    def remove_zero_intensity_pairs(spectrum):
                        pairs = spectrum.split(' ')
                        filtered_pairs = [pair for pair in pairs if not ':0' in pair and pair != '']
                        removed_pairs = [pair for pair in pairs if ':0' in pair or pair == '']
                        return ' '.join(filtered_pairs)
                        
                    newSpec = pd.DataFrame() # creating empty DF object for storing rounded spectra

                    for j in range(len(spectrum)):
                        # splitting mz & intensity values and rounding mz values
                        splitRow = spectrum[j].split(' ')
                        mz = [] # creating empty list for storing rounded mz's
                        intensity = [] # creating empty list for storing rounded intensities
                        # rounding values of mz & intensity pairs & storing them separately
                        for pair in splitRow:
                            mz.append(int(round(float(pair.split(':')[0]))))
                            intensity.append(int(round(float(pair.split(':')[1]))))
                        intensityUnited = {} # creating empty set for storing mz & intensity pairs
                        # checking for duplicate mz values for each spectrum after rounding
                        # if duplicate mz's exist, only one mz is stored & intensity values are added
                        # otherwise, one intensity values is assigned to one mz value
                        for n in range(len(mz)):
                            try:
                                intensityUnited[mz[n]] = intensityUnited[mz[n]] + intensity[n]
                            except:
                                intensityUnited[mz[n]] = intensity[n]
                        
                        # renormalizing intensities to a maximum of 1000
                        max_intensity = max(intensityUnited.values())
                        normalization_factor = 1000 / max_intensity
                        for key in intensityUnited.keys():
                            intensityUnited[key] = int(round(intensityUnited[key] * normalization_factor))
                        
                        newRow = str() # creating an empty string for storing rounded spectra

                        # storing spectra as a string
                        for key in intensityUnited.keys():
                            newRow += str(key) + ':' + str(intensityUnited[key]) + ' '
                        # storing rounded spectra as rows in Spectrum column
                        newSpec.at[j, 'Spectrum'] = newRow.strip()  # removing trailing space
                    
                    newSpec = pd.DataFrame(newSpec) # converting to DF for merging with other columns
                    print("Done rounding 'Spectra' values \n")

                    # removing pairs with intensity values containing '0'
                    newSpec_zero = newSpec['Spectrum'].apply(remove_zero_intensity_pairs)
                    print("Done removing pairs with intensity values containing '0' \n")

                    newSpec_zero_df = pd.DataFrame(newSpec_zero) # converting to DF for merging with other columns

                    def remove_low_intensity_pairs(spectrum):
                        pairs = spectrum.split(' ')
                        intensities = [float(pair.split(':')[1]) for pair in pairs]
                        max_intensity = max(intensities)
                        threshold = 0.03 * max_intensity
                        filtered_pairs = [pair for pair in pairs if float(pair.split(':')[1]) > threshold]
                        removed_pairs = [pair for pair in pairs if float(pair.split(':')[1]) <= threshold]
                        sorted_filtered_pairs = sorted(filtered_pairs, key=lambda pair: float(pair.split(':')[0]))
                        return ' '.join(sorted_filtered_pairs)
                    
                    newSpec_low = newSpec_zero_df['Spectrum'].apply(remove_low_intensity_pairs)
                    print("Done removing pairs with intensity values <= 3 percent of the highest peak \n")
                    newSpec_low = pd.DataFrame(newSpec_low) # converting to DF for merging with other columns
                    print(150 * "-")
                                        
                    frames_all = [name, rts, other, newSpec_low] # merging all columns into a list
                    allColsOut = pd.concat(frames_all, join = "outer", axis = 1) # concatenating all columns
                    merged_DF = pd.DataFrame(allColsOut) # converting to DF
                    out_DF = merged_DF.fillna('') # replacing poissible NaN values with empty cells
                    out_DF.to_csv(outFile, sep ='\t', index = False, quoting = csv.QUOTE_NONE, quotechar = '~') # writing DF to CSV; "quotechar = '~'" is used to avoid quoting quotes
                    print("Out DF: ", "\n")
                    print(out_DF)
                    print(150 * "*")

                else: # error is returned in case input data do not contain "R.T. (s)" column
                    sys.exit('File %s: %s' % (os.path.basename(filename), 'Oops! An error occured.'))

            except:
                print(150 * "!")
                print('File %s: %s' % (os.path.basename(filename), 'Oops! Something went wrong. \nPossible reasons: \ni) your input file does not contain "1st Dimension Time (s)"/"2nd Dimension Time (s)" or "R.T. (s)" columns, \nii) data were already converted, \niii) your input file does not contain any data (an empty file was submitted), or \niv) "For reasons unknown" (c) The Killers.'))
                print(150 * "*")
                sys.exit('File %s: %s' % (os.path.basename(filename), 'Oops! Something went wrong.')) # exiting the script in case of an error
    print('')
    print('Sample end time:', time.strftime("%Y-%m-%d %H:%M:%S"))
    print('')


# calling the input TXT file(s) loading dialog
inFilenames = tkinter.filedialog.askopenfilenames(defaultextension = '.txt', filetypes = [('TXT file', '.txt')], title = 'Select files for processing...')

print(150 * "*")
print('')
print('Overall start time:', time.strftime("%Y-%m-%d %H:%M:%S"))
startTime = time.strftime("%Y-%m-%d %H:%M:%S")
startTime4calc = time.time()
print(150 * "*")

# applying HiRes2LowRes function to the input TXT file(s)
if len(inFilenames) > 0:
    for name in inFilenames:
        HiRes2LowRes(name)
        print(150 * "*")
        print('Success! Sample:', os.path.basename(name))
        print(150 * "*")
        print(150 * "*")
else:
    print('Oops! No files selected. Exiting...')
    print(150 * "*")
    print(150 * "*", "\n")
    sys.exit() # exiting the script in case no files were selected

print('')
print('Overall start time:', startTime)
print('Overall end time:', time.strftime("%Y-%m-%d %H:%M:%S"))
endTime4calc = time.time()
print('Overall processing time:', str(timedelta(seconds=int(endTime4calc - startTime4calc))))
print('Number of samples converted:', len(inFilenames), "\n")
print(150 * "*")
print(150 * "*", "\n")