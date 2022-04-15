PAGE: Phase-Stretch Adaptive Gradient-field Extractor
---------------------------------------------------------------------------------
Implementation of Phase-stretch Adaptive Gradient-field Extractor (PAGE) in Python
@contributors: Madhuri Suthar, Callen MacPhee, Yiming Zhou, Dale Capewell, Jalali Lab, Department of Electrical and
Computer Engineering,  UCLA 

Please find more information within the corresponding manuscript http://arxiv.org/abs/2202.03570 [1] and the Wiki page https://en.wikipedia.org/wiki/Phase-stretch_Adaptive_Gradient-field_Extractor. 

![sunflower](https://user-images.githubusercontent.com/16159544/153650725-6e472072-4e9a-44cf-a599-02f04c0c8f31.jpeg)
_Phase-Stretch Adaptive Gradient-Field Extractor (PAGE) performed on an X-Ray of a Sunflower. The colors represent the orientation (angle) of the edge._

PAGE or Phase-stretch Adaptive Gradient-field Extractor is a physics-inspired algorithm for
detecting edges and their orientations in digital images at various scales. The algorithm is based on the diffraction
equations of optics. It builds upon the Phase Stretch Transfrom algorithm [2].

In the original implementation published in [2], the first step is to apply an adaptive tone mapping operator (TMO) to
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

![rose](https://user-images.githubusercontent.com/16159544/153651040-d32ac65c-8fc8-444f-a04c-3c641bf64bfd.gif)
_This animation demonstrates the feature channels of the PAGE algorithm. Instead of projection into the RGB space, this animation show sequentially the contents of each directional bin of the PAGE output, each of which contain edge information from the rose image._

Citations
---------
1. C. MacPhee, M. Suthar, and B. Jalali, "Phase-Stretch Adaptive Gradient-Field Extractor (PAGE)" (2022)
2. M. H. Asghari, and B. Jalali, "Edge detection in digital images using dispersive phase stretch," International Journal of Biomedical Imaging, Vol. 2015, Article ID 687819, pp. 1-6 (2015).
3. S. Madhuri and B. Jalali, “Phase-Stretch Adaptive Gradient-Field Extractor (PAGE).” (2020).

