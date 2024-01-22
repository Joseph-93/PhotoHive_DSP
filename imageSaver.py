from PIL import Image


def txt_to_png(txt_path, png_path):
    with open(txt_path, 'r') as f:
        # Read the dimensions of the image
        width, height = map(int, f.readline().split())

        # Create a new image with the specified width and height
        img = Image.new('RGB', (width, height))
        pixels = img.load()

        # Read the RGB values and set the pixels
        for y in range(height):
            for x in range(width):
                r, g, b = map(int, f.readline().split())
                pixels[x, y] = (r, g, b)

    # Save the image as a PNG
    img.save(png_path)

# Use the function for an example .txt file
filename = str(input("what is the .txt file you want to convert?\n"))
txt = 'images/output/' + filename + '.txt'
png = 'images/resulting/' + filename + '.png'
txt_to_png(txt, png)
