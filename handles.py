"""
This class handles the creation of an object that holds the parameter values of PAGE in one instance.

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
Copyright
---------
PAGE function  is developed in Jalali Lab at University of California,
Los Angeles (UCLA). More information about the technique can be found in our group
website: http://www.photonics.ucla.edu

Citations
---------
1. https://www.intechopen.com/chapters/70858
"""

class Handles:
    def __init__(self, mu_1=0.4, mu_2=0, sigma_1=0.7, sigma_2=0.08, S1=.3, S2=.3, Direction_bins=9, sigma_LPF=.1,
                 Thresh_min=-1, Thresh_max=0.0003, Morph_flag=1):
        self.mu_1 = mu_1
        self.mu_2 = mu_2

        self.sigma_1 = sigma_1
        self.sigma_2 = sigma_2

        self.S1 = S1
        self.S2 = S2

        self.Direction_bins = Direction_bins

        self.sigma_LPF = sigma_LPF

        self.Thresh_min = Thresh_min
        self.Thresh_max = Thresh_max

        self.Morph_flag = Morph_flag
