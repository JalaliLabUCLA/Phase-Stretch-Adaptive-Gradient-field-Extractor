Introduction
----------
Implementation of Phase-stretch Adaptive Gradient-Field Extractor (PAGE) in Python
@contributors: Madhuri Suthar, Callen MacPhee, Yiming Zhou, Dale Capewell, Jalali Lab, Department of Electrical and
Computer Engineering,  UCLA 

Please find more information within the corresponding manuscript http://arxiv.org/abs/2202.03570 and the Wiki page "INSERT WIKI LINK HERE". 

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
