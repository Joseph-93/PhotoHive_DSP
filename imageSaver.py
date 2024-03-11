from PIL import Image

from pathlib import Path

def list_subdirectories(base_path):
    """
    List all subdirectories in the given base path.
    """
    return [p for p in base_path.iterdir() if p.is_dir()]
    
def display_menu(subdirectories):
    """
    Display a menu of subdirectories for the user to choose from.
    """
    print("Which subfolder would you like to use?")
    for i, subdir in enumerate(subdirectories, start=1):
        print(f"{i}: {subdir.name}/")
    print("Enter the number of your choice:", end=" ")

def get_user_choice(subdirectories):
    """
    Get and validate the user's choice of subdirectory.
    """
    while True:
        user_input = input()
        try:
            choice_index = int(user_input) - 1
            # Ensure the choice is within the range of available options
            if choice_index >= 0 and choice_index < len(subdirectories):
                return subdirectories[choice_index]
            else:
                print("Invalid choice. Please enter a number from the list.")
        except ValueError:
            print("Please enter a valid number.")


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


filename = str(input("what is the .txt file you want to convert?\n"))
base_path = Path('images')
subdirectories = list_subdirectories(base_path)
display_menu(subdirectories)
chosen_subdir = get_user_choice(subdirectories)
txt = f"images/{chosen_subdir.name}/{filename}.txt"
png = 'images/resulting/' + filename + '.png'
txt_to_png(txt, png)
