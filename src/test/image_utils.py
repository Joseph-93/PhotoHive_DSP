import sys
import re
from PIL import Image
from pathlib import Path

def txt_to_png(txt_path, png_path):
    with open(txt_path, 'r') as f:
        width, height = map(int, f.readline().split())
        img = Image.new('RGB', (width, height))
        pixels = img.load()
        for y in range(height):
            for x in range(width):
                r, g, b = map(int, f.readline().split())
                pixels[x, y] = (r, g, b)
    img.save(png_path)

def image_to_txt(image_path, txt_path):
    with Image.open(image_path) as img:
        width, height = img.size
        pixels = list(img.getdata())
        with open(txt_path, 'w') as f:
            f.write(f"{width} {height}\n")
            for y in range(height):
                for x in range(width):
                    r, g, b = pixels[x + y * width][:3]  # Assuming image is in RGB or RGBA format
                    f.write(f"{r} {g} {b}\n")

def navigate_directories(base_path):
    current_path = base_path
    while True:
        subdirectories = list_subdirectories(current_path)
        if not subdirectories:
            break
        display_menu(subdirectories)
        current_path = get_user_choice(subdirectories)
    return current_path

def list_subdirectories(base_path):
    return [p for p in base_path.iterdir() if p.is_dir()]

def numerical_sort_key(path):
    numbers = re.findall(r'\d+', path.name)
    return int(numbers[0]) if numbers else float('inf')

def display_menu(options):
    options.sort(key=numerical_sort_key)  # Sort files numerically by the first number in their name
    print("Choose one of the following options:")
    for i, option in enumerate(options, start=1):
        print(f"{i}: {option.name}/")
    print("Enter the number of your choice:", end=" ")

def get_user_choice(options):
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

def choose_file(directory, file_extension):
    files = [f for f in directory.iterdir() if f.is_file() and f.suffix == file_extension]
    display_menu(files)
    return get_user_choice(files)

def convert_files(cmd):
    while True:
        base_path = Path('data')
        chosen_dir = navigate_directories(base_path)
        if cmd == 'txt2png':
            chosen_file = choose_file(chosen_dir, '.txt')
            save_extension = '.png'
        elif cmd == 'png2txt':
            chosen_file = choose_file(chosen_dir, '.png')
            save_extension = '.txt'
        else:
            print("Invalid command. Use 'txt2png' or 'png2txt'.")
            return

        print(f"Selected file: {chosen_file}")
        save_path = navigate_directories(base_path)
        save_filename = input("Enter the desired filename (without extension): ")
        save_full_path = save_path.joinpath(f"{save_filename}{save_extension}")

        if cmd == 'txt2png':
            txt_to_png(chosen_file, save_full_path)
        elif cmd == 'png2txt':
            image_to_txt(chosen_file, save_full_path)

        print(f"Saved converted file to {save_full_path}")
        if input("Convert another file? (y/n): ").lower() != 'y':
            break

def main():
    if len(sys.argv) > 1:
        cmd = sys.argv[1].lower()
        convert_files(cmd)
    else:
        print("Please specify the command: 'txt2png' or 'png2txt'")

if __name__ == "__main__":
    main()
