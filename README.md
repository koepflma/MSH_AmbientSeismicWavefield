# Code belonging to: "Examining 22 Years of Ambient Seismic Wavefield at Mount St. Helens"

## Feature extraction
This github repository provides the code to extract the following features:
- Statistical Features: root-mean-square (RMS), root-median-square (RMeS), peak-ground-velocity (PGV), and peak-ground-acceleration (PGA)
- Real-time Seismic Amplitude Measurement (RSAM)
- Displacement Seismic Amplitude Ratio (DSAR)
- Scattering Coefficients & clustering (adapted from: [scatseisnet](https://scatseisnet.readthedocs.io/en/latest/))
- Spectral width (covariance, adapted from [covseisnet](https://covseisnet.gricad-pages.univ-grenoble-alpes.fr/covseisnet/))
- z-score normalization (pers. communication A. Ardid)

## Installation
1. Clone the github repository (MSH_AmbientSeismicWavefield).
2. Navigate into the cloned github repository (MSH_AmbientSeismicWavefield).
3. Create a conda environment form ```MSH.yml``` file
4. Activate the new created conda environment (MSH).
5. Create a kernel from the new created conda environment (MSH) to run the notebooks.
   
```python
git clone git@github.com:koepflma/MSH_AmbientSeismicWavefield.git
cd MSH_AmbientSeismicWavefield
conda env create --file MSH.yml
conda activate MSH
ipython kernel install --name "MSH" --user
```

## Structure
Mainly, the sripts (```.py```) are used to calculate the features.<br>
Beside the code used for the feature extraction, some example jupyter notebooks finalize the featyre extraction and help producing some basic figures.
Everything generated will be saved in a folder called ```output```, separated by the different methods and is therein further separated in ```data``` and ```figure```. Some example outputs are provided in the folder ```output_example```.
