Tools for analysis of ZZ Ceti data

Written by JT Fuchs and P O'Brien.


all_solutions.py: Reads in a fitting_solutions.txt file and plots the best log(g) and Teff.

computechi.py: Reads in chi-square surfaces and fits them. 

plotfit.py: Reads in a normalized spectrum and model and makes a fitting plot.

** These should be combined into one file, but I'm not doing that now.***
Order of operations:

fit_all_surfaces.py: Fit chi-square surfaces of each WD. Save file as all_teff_logg_models.txt

get_coords.py: Finds each blue file, gets data from header. Saved as spec_details.txt

combine_tables: Combines all_teff_logg_models.txt and spec_details.txt. Output is model_results.txt

match_catalogs: Combines model_results.txt and catalog_data.tsv into catalog_models.txt. Sorts by RA. Adds column for duplicates