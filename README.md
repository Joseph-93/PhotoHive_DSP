# PhotoHive-DSP
PhotoHive’s library of digital signal processing for image feature extraction.

## This repository is currently capable of:
* Calculating magnitude FFT of an RGB image
* measuring an approximation of blur directionality of an image, in terms of angle and radius pairs.

## Repository Components
* main.c, a C program for fast image processing
* imageDecompressor.py, a script for reading .png files to decompressed .txt tables that main.c can easily read
* imageSaver.py, a script for compressing main.c's outputted image .txt files into valid, readable png files.
