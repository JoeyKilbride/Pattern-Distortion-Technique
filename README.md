# Pattern-Distortion-Technique
Circle detection MATLAB algorithm for extracting magnification from experimental data.

The code within contains the functions and main code to extract magnfication from the experimental images conducted in:
- Drying Dynamics publication
- Methods publication 

The algorythm is designed to extract the side of dots within an ROI, specified by the user. The dot size extraction is based on using MATLAB's regionprops() function for finding contiguous regions. 
