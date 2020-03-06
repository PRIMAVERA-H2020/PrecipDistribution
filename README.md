# Author

Ségolène Berthou, Met Office, FitzRoy Road, Exeter UK

(c) British Crown Copyright 2020, Met Office

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.

You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# Citation

Users who apply the software resulting in presentations or papers are kindly asked to cite:

    Klingaman, N. P., Martin, G. M., and Moise, A.: ASoP (v1.0): A set of methods for analyzing scales of precipitation in general circulation models, Geosci. Model Dev.,10, 57-83, doi:10.5194/gmd-10-57-2017, 2017

and 

    Berthou, S., Kendon, E., Rowell, D. P., Roberts, M. J., Tucker, S. O., & Stratton, R. A. ( 2019). Larger future intensification of rainfall in the West African Sahel in a convection‐permitting model.Geophysical Research Letters, 46, 13299–13307. https://doi.org/10.1029/2019GL083544 

as done in: 

     Demory, M.-E., S. Berthou, S. L. Sørland, M. J. Roberts, U. Beyerle, J. Seddon, R. Haarsma, C. Schär, O. B. Christensen, R. Fealy, J. Fernandez, G. Nikulin, D. Peano, D. Putrasahan, C. D. Roberts, C. Steger, C. Teichmann, R. Vautard. Can high-resolution GCMs reach the level of information provided by 12-50 km CORDEX RCMs in terms of daily precipitation distribution?. Geoscientific Model Development, submitted.

# Directions of use

PrecipDistribution calculates and plots precipitation histograms (frequency or contribution = frequency * mean bin rate)
for given bins by pooling data from a region/season.

This code is able to process a lot of datasets for different regions and seasons.
An example figure and the method explanation is given in fig.2 of this article:
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GL083544.

This code is specifically designed to compare PRIMAVERA and CORDEX precipitation distributions over Europe for the Demory et al. (submitted) GMD article.

This version of the code is set to run on the CEDA jasmin platform with python2.7 and where the PRIMAVERA and CORDEX daily precipitation datasets are located

1) create a directory where to store the pdf files (json files):
srex_pdf_json/

and an images/ directory where your code is
The images directory will store the figures created by the code

2) add the location of your datasets in european_masked_subregion.py

3) search for "sberthou" is all the code to replace all the locations to your own locations

4) now you can process the data to generate the intermediate json files (which will be used in a second step to generate the figures)
in launch_all_europe.py, set your bin options, the frequency of the data ('d' for daily), if you do 1D histograms (which is what is used with the GMD PRIMAVERA vs CORDEX paper or wet/dry spells (2D histograms used in Fig. 8 of https://link.springer.com/article/10.1007/s00382-018-4114-6 and Fig. 10 here: https://link.springer.com/content/pdf/10.1007%2Fs00382-019-04759-4.pdf).

python2.7 launch_all_europe.py

this calls launch_european_masked_subregion.py
which launches process_data_europe.py jobs on the cluster (LOTUS on jasmin)
you need to copy the mqsub file into your ~/.local/bin/

if you want to change the bin definition or the time constraint applied to the data, do so in process_data.py

process_data then calls__init__.py:
it defines the class PrecipExtremes, which will contain the histogram, centiles and metadata for the dataset (e.g. start year and end year).
Then the calculation is done and stored in a json file in srex_pdf_json/
The mask provided for the subregion in european_masked_subregion.py is linearly interpolated to the datasets. 

Note that iris is sometimes a bit harsh with metadata, which prevent concatenation of cubes: make sure all the attributes of the netCDF file are the same (e.g. if different netCDF per year or decade). If not, then you can add attributes_overwrite in callback_overwrite. This is the most common failure of the code.

5) Once you have files in srex_pdf_json/, then you can plot the figures with launch_europe_interface.py
specify a few options in that file.
It will run: 
interface_plot_european_masked_subregions.py
in series for each bigger region (e.g. "prudence", "france")

This will create a contribution histogram for a given season for different subregions (up to 6)
It will then plot the pie plot as in the GMD paper, comparing on which percentage of given intervals the distributions are different.



