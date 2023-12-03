import csv
import sys
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import *
from matplotlib.widgets import Slider
import numpy as np
import colour
from scipy.fft import *
from PIL import Image


# units
cm = 1e-2
mm = 1e-3
um = 1e-6
nm = 1e-9


# Loading the illuminant data into a dictionary as colour requires
data = np.loadtxt('CIE_std_illum_D65.csv', delimiter=',')
data_max = np.max(data[:,1])
data[:, 1] /= data_max
d65 = dict(data)
d65 = {float(k): float(v) for k, v in d65.items()}

# Function for a progress bar that changes in terminal for my white propagation
# (not in my report since this was just made in ChatGPT quickly to make my program nice)
def update_progress_bar(progress):
    bar_length = 35
    filled_blocks = '█'
    light_blocks = '░'
    filled_count = int(round(progress * bar_length))
    empty_count = bar_length - filled_count
    bar = filled_blocks * filled_count + light_blocks * empty_count
    sys.stdout.write(f'\rProgress: [ {bar} ] {int(progress * 100)}% ')
    sys.stdout.flush()

def imageMaker(image, dx, dy, size, width, height):

    # Create the "frame" that the aperture will be defined on (also the same size as the projection screen)
    img_back = Image.new('L', (dx, dy), 'black')

    # Determine the size of the aperture relative to the screen
    sizefactor = [width/size, height/size]
        
    # Define the size of the image and resize the image specified
    img_size = (int(dx * sizefactor[0]), int(dy * sizefactor[1]))
    resized = image.resize(img_size)

    # Paste the resized aperture image onto the black frame
    img_back.paste(resized, ((dx - img_size[0]) // 2, (dy - img_size[1]) // 2))

    # Create a numpy array from the resultant image and make it a "step function" to be used as the aperture function later
    img_arr = np.asarray(img_back)
    # Swap the last two kwargs to make it an object blocking instead, show's the Arago spot!!!
    img_arr = np.where(img_arr > 0, 1, 0)

    return img_arr

# To get accurate colors to display, use CIE color matching functions to convert from wavelength to sRGB colorspace
def wavelength_to_XYZ(wavelength):

    # The colour package uses wavelength without the nm adjustment
    wavelength = wavelength/nm
    # Use colour function to convert wavelength to XYZ
    xyz = colour.wavelength_to_XYZ(wavelength, cmfs=None)
    return xyz

# Same as above, but use illuminant
def wavelength_to_XYZ_illum(wavelength, illuminant):

    # The colour package uses wavelength without the nm adjustment
    wavelength = wavelength/nm
    # Use colour function to convert wavelength to XYZ
    xyz = colour.wavelength_to_XYZ(wavelength, cmfs=None)
    relative_power = illuminant[round(wavelength)]
    return xyz * relative_power
    
# Take the XYZ and convert to sRGB linear, then apply to intensity array specified to get the intensity in the proper color
def XYZ_to_sRGBlin(xyz, I):

    chromaticity_coordinates = [0.31272, 0.32903]
    srgb = colour.XYZ_to_sRGB(xyz, illuminant = chromaticity_coordinates)
    r = srgb[0] * I
    g = srgb[1] * I
    b = srgb[2] * I
    # imshow likes RGB color arrays to be in form NxMx3
    rgblin = np.dstack((r,g,b))
    return rgblin

# Convert linear sRGB to real sRGB space doing appropriate gamma corrections and clippings to get in range 0 to 1
def sRGBlin_to_sRGB(rgblin):

    rgb = np.where(rgblin < 0, 0, rgblin)
    rgb = np.where(rgb > 1.0, rgb / np.max(rgb), rgb)
    rgb = np.where(rgb <= 0.0031308, 12.92*rgb, 1.055*(np.sign(rgb) * (np.abs(rgb)) ** (1 / 2.4))- 0.055)
    return rgb

# Take above 3 functions and make into a single function
def XYZ_to_sRGB(wavelength, I):

    a = wavelength_to_XYZ(wavelength)
    a = XYZ_to_sRGBlin(a, I)
    a = sRGBlin_to_sRGB(a)
    return a

def XYZ_to_sRGB_illum(wavelength, I, illuminant):

    a = wavelength_to_XYZ_illum(wavelength, illuminant)
    a = XYZ_to_sRGBlin(a, I)
    a = sRGBlin_to_sRGB(a)
    return a


class Diffraction:
        
    def __init__(self, dx, dy, unit):
        
        self.dx = dx
        self.dy = dy
        self.unit = unit

        # Wavelength range is chosen to be the visible spectrum
        # Sensitivity refers to a scale from 0-1 for how much of the diffraction pattern is shown
        self.wavelength = 360
        self.sensitivity = 1
        self.z = 1

        self.frame_size = 10
        self.aperture_size = 1

    def aperture(self, path):

        # Sets the image specified by the file path, run this method to define the aperture
        self.path = path
        self.img = Image.open(self.path).convert('L')
    
    def propagate(self):

        # Separate attributes that scale to the specified units. Makes sliders cleaner
        self.wv = self.wavelength * nm
        self.size = self.frame_size * self.unit
        self.asize = self.aperture_size * self.unit

        # Sets the extent of plot, based on the dimensions given. Made to make the origin in the middle of the image
        self.extent = [-self.size/(self.unit*2), self.size/(self.unit*2), -self.size/(self.unit*2), self.size/(self.unit*2)]

        # Define the "electric field" or light field to be propagated based on the image generator from before
        self.E = imageMaker(self.img, self.dx, self.dy, self.size, self.asize, self.asize)

        # Finds the "spatial frequencies" or x-y components of the wave vector from each point on the frame using built-in FFT functions
        wx = fftshift(fftfreq(self.dx, (self.size/self.dx)))
        wy = fftshift(fftfreq(self.dy, (self.size/self.dy)))
        wx, wy = np.meshgrid(wx, wy)

        # Finds the overall wave vectors for each point on the frame using the fact that k=2π/λ
        k = (1 / self.wv)
        arg = (k**2 - wx**2 - wy**2)
        kz = 2 * np.pi * np.sqrt(np.abs(arg))
        # The square root is originally imaginary for some points, so use np.where which is a fast method I think
        kz = np.where(arg >= 0, kz, 1j*kz)

        # First FT into spatial frequencies, fftshift centers the origin properly
        U = fftshift(fft2(self.E))
        
        # Do the necessary propogation term e^(ikz), iFT to get back the spatial dimensions, and multiply by the complex conjugate 
        # to get the intensity (take real part because scipy returns imaginary even after this multiplication)
        I = ifft2(ifftshift(U * np.exp(1j * kz * self.z)))
        I = np.real(np.conjugate(I) * I)
        # Normalize the resultant intensity field at the screen so that sensitivity is easier to compute
        self.I = (I / np.max(I)) / self.sensitivity

        return self.I

    def get_variables(self):

        # Used for the sliders
        return [
            
            ('z', 0, 10, 'Screen Distance'), 
            ('aperture_size', 1, 50, 'Aperture Dimension'), 
            ('frame_size', 1, 50, 'Frame Dimension'), 
            ('sensitivity', 0.0001, 1, 'Sensitivity'),
            ('wavelength', 380, 780, 'Wavelength'), 
        
        ]
    
    # Gets extent and wavelength to be used for plot labels
    def getExtent(self):

        return self.extent
    
    def getWavelength(self):
        
        return self.wavelength
    
    def getSensitivity(self):
        
        return self.sensitivity

class Plotter:

    def __init__(self):

        self.fig, self.ax = plt.subplots()

    def plot(self, obj: object, cmap, color: bool):

        self.obj = obj
        self.I = self.obj.propagate()
        self.rgb = XYZ_to_sRGB(self.obj.getWavelength() * nm, self.I)
        self.color = color

        self.fig.set_figwidth(7)
        self.fig.set_figheight(7)

        # Default colormap if no colormap is given (grey scale linear map from black to white)
        cmp = LinearSegmentedColormap.from_list("bw", ["black", "white"], N=256, gamma=0.7)

        if cmap == None:
            self.cmp = cmp
        else:
            self.cmp = cmap
        
        # Use matplotlib imshow function and the extent defined before. Interpolate for a smoother plot and vmin vmax is 0-1 since
        # normalized. Lower origin also ensures that the image is flipped as it would in real life
        if color == False:
            self.l = plt.imshow(self.I, cmap=self.cmp, extent = self.obj.extent, interpolation="spline36", origin = "lower", vmin=0, vmax=1)
        if color == True:
            self.l = plt.imshow(self.rgb, extent = self.obj.extent, interpolation="spline36", origin = "lower", vmin=0, vmax=255)
        
        # Sets the axes labels based on the units given
        if self.obj.unit == cm:
            plt.xlabel('cm')
            plt.ylabel('cm')
        elif self.obj.unit == mm:
            plt.xlabel('mm')
            plt.ylabel('mm')
        elif self.obj.unit == um:
            plt.xlabel('μm')
            plt.ylabel('μm')
        elif self.obj.unit == nm:
            plt.xlabel('nm')
            plt.ylabel('nm')

        plt.title('λ = {} nm'.format(round(self.obj.getWavelength(), 1)))

        # Get the variables to be customized by sliders from the Diffraction object and add space between the sliders and plot
        _vars = self.obj.get_variables()
        plt.subplots_adjust(bottom=0.04*(len(_vars)+2))
        self.sliders = []

        # Creates a slider for every customizable variable
        for i, var in enumerate(_vars):
            self.add_slider(i*0.03, var[0], var[1], var[2], var[3])
        
        plt.show()

    def add_slider(self, pos, name, min, max, label):

        # Attempts to set increments for each slider
        self.allowed_vals = np.concatenate([np.linspace(min, max, 100), np.arange(min, max)])
        # Sets the slider's position
        ax = plt.axes([0.25, 0.025+pos, 0.55, 0.02])
        # Defines the slider
        slider = Slider(ax, label, min, max, valinit=getattr(self.obj, name), valstep=self.allowed_vals)
        # Adds slider to the list to be displayed
        self.sliders.append(slider)

        # Function that determines how the plot is updated when the slider value changes
        def update(val):
            setattr(self.obj, name, val)
            # Runs the propogation again when the slider value changes, using the new value
            if self.color == False:
                data = self.obj.propagate()
            if self.color == True:
                data = XYZ_to_sRGB(self.obj.getWavelength() * nm, self.obj.propagate())
            # Different behaviors when the frame changes (changes extent) and wavelength changes (changes title)
            if name == 'frame_size':
                updated_extent = self.obj.getExtent()
                self.l.set_extent(updated_extent)
            elif name == 'wavelength':
                self.ax.set_title('λ = {} nm'.format(round(self.obj.getWavelength(), 1)))
            # Updates the imshow data
            self.l.set_data(data)
            self.fig.canvas.draw()
        
        slider.on_changed(update)


class WhiteDiffraction:
    
    def __init__(self, z, dx, dy, size, width, height, unit):
        
        self.dx = dx
        self.dy = dy
        self.unit = unit
        
        self.z = z

        self.size = size * self.unit
        self.width = width * self.unit
        self.height = height * self.unit
        
        self.extent = [-self.size/(self.unit*2), self.size/(self.unit*2), -self.size/(self.unit*2), self.size/(self.unit*2)]

    def aperture(self, path):

        # Sets the image specified by the file path, run this method to define the aperture
        self.path = path
        self.img = Image.open(self.path).convert('L')
    
    def whitepropagate(self, wavelength):
        
        self.wavelength = wavelength * nm

        # Define the "electric field" or light field to be propagated based on the image generator from before
        self.E = imageMaker(self.img, self.dx, self.dy, self.size, self.width, self.height)

        # Finds the "spatial frequencies" or x-y components of the wave vector from each point on the frame using built-in FFT functions
        wx = fftshift(fftfreq(self.dx, (self.size/self.dx)))
        wy = fftshift(fftfreq(self.dy, (self.size/self.dy)))
        wx, wy = np.meshgrid(wx, wy)

        # Finds the overall wave vectors for each point on the frame using the fact that k=2π/λ
        k = (1 / self.wavelength)
        arg = (k**2 - wx**2 - wy**2)
        kz = 2 * np.pi * np.sqrt(np.abs(arg))
        # The square root is originally imaginary for some points, so use np.where which is a fast method I think
        kz = np.where(arg >= 0, kz, 1j*kz)

        # First FT into spatial frequencies, fftshift centers the origin properly
        U = fftshift(fft2(self.E))
        
        # Do the necessary propogation term e^(ikz), iFT to get back the spatial dimensions, and multiply by the complex conjugate 
        # to get the intensity (take real part because scipy returns imaginary even after this multiplication)
        I = ifft2(ifftshift(U * np.exp(1j * kz * self.z)))
        I = np.real(np.conjugate(I) * I)
        # Normalize the resultant intensity field at the screen so that sensitivity is easier to compute
        self.I = (I / np.max(I))

        return self.I

    def stack(self, iterations: int, wv_range: tuple):
    
        # Specify my iterations and wavelength range
        self.iterations = iterations
        self.wv_range = wv_range

        # If not within the colour package limitations, stop the method and say what's wrong
        if self.wv_range[0] < 360 or self.wv_range[1] > 780:
            print("Range must be between 360 and 780!")
            return

        if self.iterations > 420:
            print("Iterations cannot exceed 420!")
            return
        
        # Get the amount of wavelengths I have in my range
        self.range_difference = self.wv_range[1]-self.wv_range[0]

        # Find the change in wavelength from iterations (floor since CIE color matching functions are only in 1nm intervals)
        delta_wavelength = self.range_difference // self.iterations

        # Make a copy of my illuminant
        d65_instance = d65.copy()
        d65_instance = {key: value / self.iterations for key, value in d65_instance.items()}
        
        # Make an empty array to be added to
        self.rgb = np.zeros((self.dx, self.dy, 3))
        
        # Uses sys to make a start time
        start_time = time.time()
        for i in range(self.iterations + 1):
            # Update my progress bar as a percentage done
            progress = i / self.iterations
            update_progress_bar(progress)
            # Add my wavelength delta each iteration
            self.wv = self.wv_range[0] + (delta_wavelength*i)
            # Can't go over colour's required wavelength range so stops if it does
            if self.wv > 780 or self.wv < 360:
                break
            # Runs my propagation and applies my color with illuminant
            I = self.whitepropagate(self.wv)
            arr = XYZ_to_sRGB_illum(self.wv * nm, I, d65_instance)
            self.rgb += arr
            # Adds to my progress bar an elapsed time
            sys.stdout.write(f'{time.time()-start_time:.2f}s')
        end_time = time.time()

        total_time = end_time - start_time

        # Clears my progress bar and displays the total time to run
        sys.stdout.write(f'\r\033[KCompleted in {total_time:.2f} seconds!\n')

        # Truncates my data, since sRGB is only from a range of (0,1)
        self.rgb = np.where(self.rgb > 1, 1, self.rgb)

        return self.rgb

    # Gets extent and wavelength to be used for plot labels
    def getExtent(self):

        return self.extent
    
    def getWavelength(self):
        
        return self.wavelength

    def getZ(self):

        return self.z

class WhitePlotter:
    
    def __init__(self):

        self.fig, self.ax = plt.subplots()

    # Self explanatory!
    def plot(self, obj: object, I):

        self.obj = obj
        self.I = I
        
        self.fig.set_figwidth(7)
        self.fig.set_figheight(7)

        self.l = plt.imshow(self.I, extent = self.obj.extent, interpolation="spline36", origin = "lower", vmin=0, vmax=1)
        
        # Sets the axes labels based on the units given
        if self.obj.unit == cm:
            plt.xlabel('cm')
            plt.ylabel('cm')
        elif self.obj.unit == mm:
            plt.xlabel('mm')
            plt.ylabel('mm')
        elif self.obj.unit == um:
            plt.xlabel('μm')
            plt.ylabel('μm')
        elif self.obj.unit == nm:
            plt.xlabel('nm')
            plt.ylabel('nm')

        # Title is instead distance since this is the only major factor
        plt.title('z = {} m'.format(round(self.obj.getZ(), 1)))

        plt.show()