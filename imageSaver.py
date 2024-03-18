from PIL import Image

from pathlib import Path

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


def navigate_directories(base_path):
    """
    Navigate through the directory structure interactively until a directory with no subdirectories is found.
    """
    current_path = base_path
    while True:
        subdirectories = list_subdirectories(current_path)
        if not subdirectories:  # If no subdirectories, break the loop
            break
        display_menu(subdirectories)
        current_path = get_user_choice(subdirectories)
    return current_path

def list_subdirectories(base_path):
    """
    List all subdirectories in the given base path.
    """
    return [p for p in base_path.iterdir() if p.is_dir()]

def display_menu(options):
    """
    Display a menu of options for the user to choose from.
    """
    print("Choose one of the following options:")
    for i, option in enumerate(options, start=1):
        print(f"{i}: {option.name}/")
    print("Enter the number of your choice:", end=" ")

def get_user_choice(options):
    """
    Get and validate the user's choice.
    """
    while True:
        user_input = input()
        try:
            choice_index = int(user_input) - 1
            if 0 <= choice_index < len(options):
                return options[choice_index]
            else:
                print("Invalid choice. Please enter a number from the list.")
        except ValueError:
            print("Please enter a valid number.")

def choose_file(directory):
    """
    Let the user choose a file from the directory.
    """
    files = [f for f in directory.iterdir() if f.is_file()]
    display_menu(files)
    chosen_file = get_user_choice(files)
    return chosen_file


def main():
    base_path = Path('images')
    chosen_dir = navigate_directories(base_path)
    chosen_file = choose_file(chosen_dir)
    
    print(f"Selected file: {chosen_file}")

    save_path = navigate_directories(Path.cwd())
    save_filename = input("Enter the desired filename (without extension): ")
    save_full_path = save_path.joinpath(f"{save_filename}.png")

    txt_to_png(chosen_file, save_full_path)
    print(f"Saved converted image to {save_full_path}")

if __name__ == "__main__":
    main()