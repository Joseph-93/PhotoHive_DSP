# Create the lookup table and save it to fft_normalizer.txt

import math
import matplotlib.pyplot as plt


Y_SCALE_FACTOR = 255/math.log(2) # f(1) when M=1 is ln(2)

class LookupTable():
    input = []
    output = []

# Function used for normalization found here:
# https://www.desmos.com/calculator/fi7w3o0ced
def main(filename="lookups/fft_normalizer.txt"):
    length = 256
    lookup = LookupTable()
    with open(filename, "w") as file:
        file.write("1D\n")
        file.write(f"{length}\n")
        for index in range(length):
            output_value = index / (length - 1)
            input_value = (math.exp(index/Y_SCALE_FACTOR)-1)**2
            lookup.input.append(input_value)
            lookup.output.append(output_value)
            file.write(f"{input_value} {output_value}\n")

    plt.plot(lookup.input, lookup.output)
    plt.title("Logarithmic Transformation Lookup Table")
    plt.xlabel("Original Value")
    plt.ylabel("Transformed Value")
    plt.grid(True)
    plt.show()

# Create the lookup table and save it to "fft_normalizer.txt"
main()
