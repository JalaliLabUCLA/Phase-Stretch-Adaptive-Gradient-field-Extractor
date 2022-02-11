"""
This script is a generic test script that plots the output of the PAGE algorithm as projected into the RGB space, with
edge angle corresponding to color. This test script relies exclusively on CPU.

Implementation of Phase-stretch Adaptive Gradient-Field Extractor (PAGE) in Python
@contributors: Madhuri Suthar, Callen MacPhee, Yiming Zhou, Dale Capewell, Jalali Lab, Department of Electrical and
Computer Engineering,  UCLA

PAGE or Phase-stretch Adaptive Gradient-field Extractor is a physics-inspired algorithm for
detecting edges and their orientations in digital images at various scales. The algorithm is based on the diffraction
equations of optics.

In the original implementation published in [1], the first step is to apply an adaptive tone mapping operator (TMO) to
enhance the local contrast. Next, we reduce the noise by applying a smoothening kernel in frequency domain
(this operation can also be done in spatial domain). We then apply a spectral phase kernel that emulates the
birefringence and frequency channelized diffractive propagation. The final step of PAGE is to apply thresholding and
morphological operations on the generated feature vectors in spatial domain to produce the final output. The PAGE output
embeds the original image into a set of feature maps that select semantic information at different scale, orientation,
and spatial frequency. These feature maps include spacial frequency bins and edge angle bins.

The spectral phase operator is expressed as a product of two phase functions, ϕ_1 and ϕ_2. The first component ϕ_1 is a
symmetric gaussian filter that selects the spatial frequency range of the edges that are detected. Default center
frequency is 0, which indicates a baseband filter, the center frequency and bandwidth of which can be changed to probe
edges with different sharpness — in other words to probe edges occurring over different spatial scales. The second
component, ϕ_2, performs the edge-detection. Since the output is based on the phase, it needs to be a complex-valued
function. The PAGE operation transforms a real-value input to a complex-value quantity from which the phase is
extracted.

This code is a simplified version of the full PAGE algorithm. It does not include the adaptive tone mapping operator
introduced in [1], and it only includes one spatial frequency bin with a preselected bin central frequency, mu and bin
bandwidth, sigma. It outputs a collection of angle dependent edges, each corresponding to one angle bin. For an N x M
dimensional input, the output will be of size N x M x D, where D is the number of directional bins. This could be
repeated for each color channel of a color image.
Parameters
----------
mu_1            # Center frequency of a normal/Gaussian distributed passband filter Phi_1
mu_2            # Center frequency of log-normal distributed distributed passband filter Phi_2
sigma_1         # Standard deviation sigma of normal/Gaussian distributed passband filter Phi_1
sigma_2         # Standard deviation sigma of log-normal distributed passband filter Phi_2
S1              # Strength (Amplitude) of Phi_1 filter
S2              # Strength (Amplitude) of Phi_2 filter
Direction_bins  # Number of directional bins i.e. number of PAGE filter channels
sigma_LPF       # Standard deviation sigma of Gaussian distribution for smoothening kernel
Thresh_min      # Lower bound of bi-level (bipolar) feature thresholding for morphological operations
Thresh_max      # Upper bound of bi-level (bipolar) feature thresholding for morphological operations
Morph_flag      # Flag to choose (0) analog edge output or (1) binary edge output
See class Handles for parameter information.
Citations
---------
1. https://www.intechopen.com/chapters/70858
"""
import os
from timeit import default_timer as timer 
import mahotas as mh
import numpy as np

from page_cpu import PAGE_CPU, visualize_PAGE
from handles import Handles

#used to time blocks of code 

input_path = os.getcwd()  # This is where the code is running.
input_path = os.path.join(input_path, 'Input_Images/')

# Image and parameter Selection
Image_list=['wind_rose.png', 'jet_engine.jpeg', 'Prepona_laertes.jpeg', 'barbara.jpeg', 'X-Ray-Sunflower.jpg']  # Create list of images to run PAGE upon
# Use index below to choose test image above
image_index = 0  # Index of image in list to create filepath and run, may be looped over for many images.

################################################### PARAMETERS #########################################################
mu_1=0              # Center frequency of a normal/Gaussian distributed passband filter Phi_1
mu_2=0.35           # Center frequency of log-normal distributed distributed passband filter Phi_2
sigma_1=0.08        # Standard deviation sigma of normal/Gaussian distributed passband filter Phi_1
sigma_2=0.7         # Standard deviation sigma of log-normal distributed passband filter Phi_2
S1=0.4              # Strength (Amplitude) of Phi_1 filter
S2=0.4              # Strength (Amplitude) of Phi_2 filter
Direction_bins=30   # Number of directional bins i.e. number of PAGE filter channels, default bin width is 180/30 = 6 degrees
sigma_LPF=.1        # Standard deviation sigma of Gaussian distribution for smoothening kernel
Thresh_min=-1       # Lower bound of bi-level (bipolar) feature thresholding for morphological operations
Thresh_max=0.0003   # Upper bound of bi-level (bipolar) feature thresholding for morphological operations
Morph_flag=1        # Flag to choose (0) Analog Edge output or (1) binary edge output
########################################################################################################################

filepath = os.path.join(input_path, Image_list[image_index])  #  Create filepath string of image     
Image_orig = mh.imread(filepath)
Image_orig = mh.colors.rgb2grey(Image_orig)
Image_orig = mh.imresize(Image_orig, [640, 640])  

# Create handles instance with values as described above
handles = Handles(mu_1, mu_2, sigma_1, sigma_2, S1, S1, Direction_bins, sigma_LPF, Thresh_min, Thresh_max, Morph_flag)

# Convert image datatype to double for calculations
Image_grey_double = np.double(Image_orig);

# Run page on input image with parameters defined in Handles and time it 
start = timer()
[PAGE_output, PAGE_Kernel]=PAGE_CPU(Image_grey_double,handles);
print("Computation time on CPU:", timer()-start)

# Create a weighted color image of PAGE output to visualize directionality of edges
visualize_PAGE(PAGE_output, Image_grey_double, handles)

