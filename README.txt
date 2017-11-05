cellCycleProtocol:
Nuclear fluorescent stains such as DAPI can produce a signal that is proportional to the concentration of nuclear DNA.  Therefore, the integrated intensity of a nuclei's fluorescent signal can suggest whether a cell is in the G1, S or G2 stage of the cell cycle.  In this script I provide a workflow to process images of fluorescently statined nuclei.  The integrated intensity of segmented nuclei is then used to assign cells to specific stages of the cell cycle.   

Getting Started:

A detailed description of the workflow with reference to sample images and settings can be found on MATLAB File Exchange (search Author: Brian DuChez).
To give a brief overview of how to prepare data for using the scirpt...

Organize Images:
1. Create individual folders for each experimental group. 2. Title each folder according to the experimental conditions. 3. Place images of fluorescently stained nuclei in the appropriate folders (images should be 2-dimensional and of a file format that can be read using the imread function). 4. Ensure folders are in the current working directory of the script.

Summary Workflow:
1. Optimize three essential parameters on target image(s) to identify nuclei of interest. 2. Validate parameter settings on all images for each experimental condition (Optional). 3. Obtain integrated intensity values for all nuclei from each group. 4. Create histograms of integrated intensity for each group. 5. Use the "untreated" group histogram as a reference to determine the relative fluorescence that corresponds to cells in G1, S or G2 stages of the cell cycle. 6. Calculate proportion of cells in each stage for all groups evaluated.

Acknowlegements:
expandAxes.m and toggleFig.m were written by Brett Shoelson and are used in this workflow to visualize images.