# Background

I created this project originally for a Python class, but expanded upon in into a full diffraction simulation! I plan on updating this overtime to optimize and add more interactivitiy options.

# How to Use

## Packages

This simulator primarily uses Pillow, NumPy, SciPy, and Matplotlib. This also uses Colour which can be installed via:

```python
pip install --user colour-science
```

## Code to Run

To use my angular spectrum diffraction simulator, run

```python
python -i angular_spectrum_diffraction.python
```

To use the interactive, single color simulator, run

```python
d = Diffraction(dx, dy, unit)
d.aperture('image path')
d.propagate()

k = Plotter()
k.plot(d, cmap, color)
```

Where dx and dy are the resolution of the simulator, unit is the units for the aperture/screen dimensions, image path is the image in the same directory to be used as an aperture, cmap is a specific color map to be used by matplotlib (None defaults to a grey map), and color is a boolean (True for accurate colors, False for grey). Use None for cmap if you want to use accurate colors.

To use the white light diffraction simulator, run

```python
w = WhiteDiffraction(z, dx, dy, screen_size, aperture_width, aperture_height, unit)
w.aperture('image path')
rgb = w.stack(iterations, wavelength range)
h = WhitePlotter()
h.plot(w, rgb)
```

These parameters are mostly the same except for a few. z is the screen distance, screen_size is the screen
dimension, aperture_width and aperture_height are the aperture dimensions, iterations is the amount of
wavelengths to be sampled, and wavelength range is the range of wavelengths you want to simulate.

## Recommendations

There are some apertures provided, the most interesting being double slit, circle, hexagon, and triangle.

For the white light simulator, I recommend starting with a resolution of 1000, a distance of 1 meter, screen size of 20 mm, aperture width and height of 1 mm, 420 iterations, and a range of (360,780) (for the whole range!). Generally, a smaller distance or bigger screen size means smaller pattern. Setting z=0 just gives you back the aperture plotted since it's a Fourier Transform! Also doubleslitthin.png is my favorite to run for the aperture.

For the units, type 1 for meters, cm for centimeters, mm for millimeters, um for micrometers, and nm for nano meters.
