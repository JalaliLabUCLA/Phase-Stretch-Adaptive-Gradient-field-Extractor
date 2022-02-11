"""
Implementation of Phase-stretch Adaptive Gradient-Field Extractor (PAGE) in Python
@contributors: Madhuri Suthar, Callen MacPhee, Yiming Zhou, Dale Capewell, Jalali Lab, Department of Electrical and
Computer Engineering,  UCLA PAGE or Phase-stretch Adaptive Gradient-field Extractor is a physics-inspired algorithm for
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
import numpy as np
import mahotas as mh
import matplotlib
import matplotlib.pyplot as plt
import torch

tfloat = torch.float64
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

def cart2pol_torch(x, y):
    theta = torch.atan2(y, x)
    rho = torch.hypot(x, y)

    return (theta, rho)

def PAGE_torch(I,handles):
    # Define two dimensional cartesian (rectangular) vectors, X and Y
    Image_orig_size=I.shape
    L=0.5
    u = torch.linspace(-L, L, I.shape[0], dtype = tfloat, device=device)
    v = torch.linspace(-L, L, I.shape[1], dtype = tfloat, device=device)
    [U1, V1] = (torch.meshgrid(u, v))
    U = torch.transpose(U1,0,1)
    V = torch.transpose(V1,0,1)

    [THETA, RHO] = cart2pol_torch(U, V)

    # Low pass filter the original image to reduce noise
    Image_orig_f = torch.fft.fft2(I)
    expo = torch.fft.fftshift(torch.exp(-0.5*torch.pow((torch.divide(RHO, np.sqrt((handles.sigma_LPF ** 2) / np.log(2)))), 2)))
    Image_orig_filtered = torch.real(torch.fft.ifft2((torch.mul(Image_orig_f, expo))))

    # Create PAGE directional filters
    Min_Direction = 1 * np.pi / 180 # Direction in radians
    Direction_span = np.pi / handles.Direction_bins # Direction step in radians
    Direction = torch.arange(Min_Direction, np.pi, Direction_span) # Direction array in radians

    # Define the dimension of filter
    X_size = Image_orig_size[0]
    Y_size = Image_orig_size[1]
    Z_size = Direction.shape[0]
    PAGE_filter_array = torch.zeros([X_size, Y_size, Z_size], dtype = tfloat, device=device)

    # For the number of directional bins, create PAGE kernels based on spectral directionality
    for i in range(Z_size):
        tetav = Direction[i]

        # Project onto new directionality basis for PAGE filter creation
        Uprime = U*torch.cos(tetav) + V*torch.sin(tetav)
        Vprime = -U * torch.sin(tetav) + V * torch.cos(tetav)

        # Create Normal component of PAGE filter
        Phi_1 = torch.exp(-0.5 * ((abs(Uprime) - handles.mu_1) / handles.sigma_1)** 2) / (
                    1 * np.sqrt(2 * np.pi) * handles.sigma_1)
        Phi_1 = (Phi_1 / torch.max(Phi_1[:])) * handles.S1

        # Create Log-Normal component of PAGE filter ()
        Phi_2 = torch.exp(-0.5 * ((torch.log(abs(Vprime)) - handles.mu_2) / handles.sigma_2) ** 2) / (
                    abs(Vprime) * np.sqrt(2 * np.pi) * handles.sigma_2)
        Phi_2 = (Phi_2 / torch.max(Phi_2[:])) * handles.S2

        # Add overall directional filter to PAGE filter array
        PAGE_filter_array[:,:,i]= torch.mul(Phi_1,Phi_2) # log Normal distribution for X * Gaussian Distribution for Y

    # Initialized PAGE output based on original image and number of directional bins
    PAGE_output = torch.zeros((Image_orig_filtered.shape[0],Image_orig_filtered.shape[1],Direction.shape[0]), dtype=tfloat, device=device)

    # Looping over number of directional bins, apply each PAGE directional kernel and save as channel of PAGE_output
    for j in range(Direction.shape[0]):
        # Apply PAGE Kernel in spectral domain
        temp = (torch.fft.fft2(Image_orig_filtered)) * torch.fft.fftshift(torch.exp(-1j * PAGE_filter_array[:,:,j]))
        Image_orig_filtered_PAGE = torch.fft.ifft2(temp)

        # Take phase of output in spacial domain
        PAGE_Features = torch.angle(Image_orig_filtered_PAGE)

        # Morphological operations: For Morph_flag == 0, return analog edges. For Morph_flag == 1, erode edges to binary.
        if handles.Morph_flag == 0:
            out = PAGE_Features
            PAGE_output[:,:,j] = out
        else:
            features = torch.zeros(PAGE_Features.shape);
            features[torch.where(PAGE_Features > handles.Thresh_max)] = 1
            features[torch.where(PAGE_Features < handles.Thresh_min)] = 1  # Output phase has both positive and negative values
            features[torch.where(I < torch.max(torch.max(I)) / 20)] = 0  #  Ignore the features in the very dark areas of the image

            # do this as np.array because torch has no known equivalence
            out = np.array(features)
            out = mh.thin(out, 1)
            out = mh.bwperim(out, 4)
            out = mh.thin(out, 1)
            out = mh.erode(out, np.ones((1, 1)))
            
            #convert back to torch array
            out = torch.tensor(out, dtype=tfloat, device=device)
        PAGE_output[:,:,j]=out

    return PAGE_output, PAGE_filter_array #returns torch tensor


def visualize_PAGE(PAGE_output, Image_grey_double, handles):
    # Create a weighted color image of PAGE output to visualize directionality of edges
    weight_step = 255 * 3 / PAGE_output.shape[2]
    color_weight = np.arange(0, 255, weight_step)
    Edge = np.zeros((Image_grey_double.shape[0], Image_grey_double.shape[1], 3))
    step_edge = int(round(handles.Direction_bins / 3))

    # Project PAGE channels into one RGB image with color representing directionality
    for i in range(int(round(PAGE_output.shape[2] / 3))):
        Edge[:, :, 0] = color_weight[i] * PAGE_output[:, :, i] + Edge[:, :, 0]
        Edge[:, :, 1] = color_weight[i] * PAGE_output[:, :, i + step_edge] + Edge[:, :, 1]
        Edge[:, :, 2] = color_weight[i] * PAGE_output[:, :, i + (2 * step_edge)] + Edge[:, :, 2]

    Edge_normalized = (Edge - np.min(Edge)) / (np.max(Edge) - np.min(Edge))

    # Create colorbar for associated angular bins and plot the output
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]}, facecolor='black', figsize=(4, 4),
                                   dpi=200)
    ax1.imshow(Edge_normalized)
    clist = [(0, "purple"), (45. / 180., "red"), (90. / 180., 'lime'), (135. / 180., "blue"), (180. / 180., "purple")]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("name", clist)
    colors = cmap(np.arange(cmap.N))
    ax2.imshow([colors], extent=[0, 180, 0, 1], aspect=7)
    ax2.set(xticks=[0, 45, 90, 135, 180], yticks=[])
    labels = ["0\u00b0", "45\u00b0", "90\u00b0", "135\u00b0", "180\u00b0"]
    ax2.set_xticklabels(labels, color='white')
    plt.show()
