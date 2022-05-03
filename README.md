This github directory contains code to determine the geometric efficiency of a circular/annular detector with circular uniform/gaussian sources, considering isotropic emission.

The program uses monte carlo methods to pick a starting location and initial direction, which are than used the extrapolate the trajectories. For each detector-trajectory pair, it is determined whether or not the emission ens up in the detector. The errors can be estimated straightforwardly with Poisson statistics.

After cloning the git, execute the command "chmod +x ./build.sh" once to set up permission to use the shell script build.sh in the emission folder.

To run, from main path: "./build.isotropic.exe 'source' 'detector'"
Where 'source' can be 'uniform' or 'gaussian'; 'detector can be 'circular' or 'annular'
The program will ask for some parameters:
zmin/rd: the minimal distance for which the  geometric efficiency will be calculated (in detector radius units; for annular = outer radius);
zmax/rd: the maximal distance for which the  geometric efficiency will be calculated (in detector radius units; for annular = outer radius);
number of points: number of linspace points in which the geometric efficiency is calculated;
source/rd: source spread (in detector radius units; for annular = outer radius); for circular source = source radius, for gaussian source = sigma;
Power: 10^x monte carlo points used per distance;
Detector outer/inner: ratio of outer radius to inner radius (only for annular detector);
Filename: name of output Filename;

The ouput file is of csv type with columns showing: distance_from_source(detector radius units)    point_source_approximation   model_value relative_uncertainty(%)

For more information about the code: contact 'michael.heines@kuleuven.be'
