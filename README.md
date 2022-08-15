# HiCNoiseMeasurer

This script is designed to approximate the noise in Hi-C data. HiCNoiseMeasurer is built on the principle that signal in neighboring bins is generally correlated, and therefore computes the auto-correlation function (acf) along rows. Essentially, this tool measures how smooth the signal is. This is probably most effective at finer-scale resolutions where the polymer most likely restricts neighboring bins from too different, but we have tested it with coarse bin sizes with success. We recommend using this tool in conjunction with read subsampling, [HiCSampler](https://github.com/JRowleyLab/HiCSampler), to estimate the effect that sequencing depth will have on the data.   

HiCNoiseMeasurer takes Hi-C files in [juicer](https://github.com/aidenlab/juicer) format (.hic). It outputs a file that summarizes the input file names, resolutions, and the corresponding noise value.

## Dependencies:
<br>python3<br>numpy<br>argparse<br>subprocess<br>math<br>statsmodels<br>juicer tools<br>

## Usage: 

python noiseMeasure.py [-h] -i FILE --res str [-c CHROM] -j JUICERTOOLS [-o OUTPUTFILE]

-h, --help show this help message and exit

-i FILE, --hicfile FILE Input .hic file(s), can be comma separated list

--res str Resolution(s) to process the Hi-C file, can be comma separated list

-c CHROM Chromosome to examine

-j JUICERTOOLS Path to juicer tools

-o OUTPUTFILE Path to output file. If none specified, will output to terminal.

## Example:

### Single File:
python noiseMeasure.py -i InputFile.hic -j juicer_tools_1.14.08.jar --res 5000 -c 1 -o OutputFile.txt

The above command will calculate the acf at 5 kb resolution on our input .hic file and write the results to OutputFile.txt.

### Multiple Files:
python noiseMeasure.py -i InputFile1.hic,InputFile2.hic,InputFile3.hic -j juicer_tools_1.14.08.jar --res 5000,10000,25000 -c 1 -o OutputMultiple.txt

The above command will calculate the acf at 5, 10, and 25 kb for each .hic file listed.

## OutputFile Format:
InputFileName Resolution Score

Example: 

InputFile1.hic 5000 4.013<br>
InputFile2.hic 5000 3.489<br>
InputFile3.hic 5000 6.120<br>
InputFile1.hic 10000 2.014<br>
InputFile2.hic 10000 2.542<br>
InputFile3.hic 10000 2.222<br>
