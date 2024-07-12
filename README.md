# Code belonging to: "Examining 22 Years of Ambient Seismic Wavefield at Mount St. Helens"

## Paper Abstract

An increase in seismic activity precedes most volcanic eruptions. Whereas event‐based forecasting approaches have been successful, some eruptions remain unanticipated, resulting in casualties and damage. Our study leverages the recent advancements in ambient field seismology. We explore features extracted from continuous ambient fields using traditional methods, for example, peak ground velocity, peak ground acceleration, root mean square, root median square, real‐time seismic amplitude measurement, and novel methods (displacement seismic amplitude ratio and spectral width). In addition, we explore unsupervised learning of higher order wavelet features using scattering networks. We find that combining all the methods was necessary to disentangle the effects of seismic sources from structural changes at Mount St. Helens. Although the ambient wavefield‐based approach does not yield additional or more significant precursory signals than event‐based methods at Mount St. Helens, our study demonstrates that the ambient wavefield provides supplementary information, mainly about structural changes and complements traditional methods. The ambient seismic wavefield offers additional insights into long‐lasting processes. We find enhanced wave attenuation correlating with geochemical measurements. We interpret this as ongoing structural changes, such as dome growth or the evolution of the volcanic conduit system. On annual and decadal timescales, we interpret seasonal seismic attenuation in the shallow subsurface as groundwater fluctuations, corroborated by observations at the nearby Spirit Lake level. This multimethod approach at Mount St. Helens sheds light on a volcanic system’s underlying dynamics and structure.

You can find the published paper [here](https://pubs.geoscienceworld.org/ssa/srl/article/doi/10.1785/0220240079/644839/Examining-22-Years-of-Ambient-Seismic-Wavefield-at).

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
