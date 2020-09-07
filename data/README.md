# Sample Data for 3D Fluid Flow Estimation with Integrated Particle Reconstruction

The sample data was created using the forced isotropic turbulence dataset of the [Johns Hopkins Turbulence Database (JHTDB)](http://turbulence.pha.jhu.edu/).
If you use this data please [cite](http://turbulence.pha.jhu.edu/citing.aspx) accordingly.

The setup is based on Case D of the [4th International PIV Challenge](https://www.pivchallenge.org/pivchallenge4.html). I.e., the discretization level, measurement area size and camera parameters are adopted. The measurement width was reduced to 51.2mm (originally 204.8), resulting in a smaller measurement volume size of 51.2x25.6x17.6mm and an image size of 1400x800 pixels. 

A pinhole camera is used for the particle rendering. The file calibpoints.txt contains 4000 point correspondences of the 3D world coordinates and 2D image coordinates of four cameras (X Y Z x1 y1 x2 y2 x3 y3 x4 y4). This can be used to calibrate a pinhole or polynomial camera model.

Particle images of four camera views are given for three time steps and two different particle densities (0.175ppp and 0.125ppp). All particles are rendered as Gaussian blobs with sigma=1 and varying intensity. Ground truth is given as a .mat file, containing:
- the ground truth particle point locations and intensities for all three time steps (gt_part, dimensions: 4 x numPts x 3), 
- the flow field for the second (middle/reference) time step (gt_flowgrid_t1), 
- the measurement volume dimensions and offset 
- and the pinhole camera parameters P.

Note that the grid sampling of the ground truth flow field is in general not the same as the discretization level of the estimated flow field of our algorithm. For comparison, the estimated flow field needs to be interpolated to the ground truth grid first.
