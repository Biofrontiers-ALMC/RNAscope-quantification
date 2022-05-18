# RNAscope quantification

Welcome to the RNAscope quantification project! This code is being 
developed in collaboration with Connor McNulty and the Baratta lab at CU 
Boulder.

This project has two goals:
1. Identify GFP labeled cells
2. Identify and count the number of cells with RNA labeling.

## Installation and Usage

### Required MATLAB toolboxes

This project requires the following MATLAB toolboxes to be installed:
* Image Processing Toolbox

### Additional toolboxes

This project also requires the following custom toolboxes:
* [Bioformats Image Toolbox](https://github.com/Biofrontiers-ALMC/bioformats-matlab)

### Code

1. **Download** the code by clicking on the blue Code button above, then 
   selecting "Download ZIP".
2. **Extract** the ZIP file into a directory accessible by MATLAB (e.g.
   in Documents/MATLAB). Navigate to this folder in MATLAB.
3. **Add the folder to the MATLAB path** by right-clicking on
   the extracted folder in the Current Folder window in MATLAB, then select 
   "Add to Path" > "Selected Folders and Subfolders".

## Running the code

Run the file `countCells` in MATLAB. The parameters needed to run the code
 are now located at the top of the script:

* file - Path to ND2 file to process
* nuclThreshold - Greyscale to threshold nuclei in GFP channel
* maxNuclSize - Largest valid nucleus in pixels
* ch570Threshold - Threshold to count puncta in the 570 nm channel by maximum brightness in the cell
* ch650Threshold - Threshold to count puncta in the 650 nm channel by maximum brightness in the cell

## Contact

If there are any issues or bugs, please contact me at jian.tay AT colorado.edu

