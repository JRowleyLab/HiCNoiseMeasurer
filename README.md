# HiCNoiseMeasurer

This script is designed to approximate the noise in Hi-C data. HiCNoiseMeasurer is built on the principle that signal in neighboring bins is generally correlated, and therefore computes the auto-correlation function (acf) along rows. Essentially, this tool measures how smooth the signal is. This is probably most effective at finer-scale resolutions where the polymer most likely restricts neighboring bins from too different, but we have tested it with coarse bin sizes with success. We recommend using this tool in conjunction with read subsampling, [HiCSampler](https://github.com/JRowleyLab/HiCSampler), to estimate the effect that sequencing depth will have on the data.   

HiCNoiseMeasurer takes Hi-C files in [juicer](https://github.com/aidenlab/juicer) format (.hic). It outputs a file that summarizes the input file names, resolutions, and the corresponding noise value.

Dependencies:<br>python3<br>numpy<br>argparse<br>subprocess<br>math<br>statsmodels<br>juicer tools<br>
Usage: python noiseMeasure.py --help

