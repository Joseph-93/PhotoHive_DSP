from PIL import Image


# Open image and read out to a text file for C program to easily read
def image_to_txt(image_path, txt_path):
    with Image.open(image_path) as img:
        width, height = img.size
        pixels = list(img.getdata())

        with open(txt_path, 'w') as f:
            # Write image dimensions at the top of the file
            f.write(f"{width} {height}\n")

            # Write the RGB values of each pixel
            for y in range(height):
                for x in range(width):
                    r, g, b, a = pixels[x + y * width]
                    f.write(f"{r} {g} {b}\n")


# Use the function for an example image
filename = str(input("what is the .png file you want to convert?\n"))
png = 'images/original/' + filename + '.png'
txt = 'images/readable/' + filename + '.txt'
image_to_txt(png, txt)
