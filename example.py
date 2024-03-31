import sys
from PhotoHive_DSP_lib import set_bounding_boxes, get_report

from PIL import Image


def run_demonstration():
    global image_name
    image_path = "data/original/"
    # Check if a custom image path was provided as a command-line argument
    if len(sys.argv) > 1:
        image_name = sys.argv[1]
    else:
        image_name = "10.png"

    # Open image
    image = Image.open(image_path+image_name)

    bounds = [
        {
            "left": 61,
            "right": 383,
            "top": 212,
            "bottom": 897,
        },
        {
            "left": 363,
            "right": 591,
            "top": 130,
            "bottom": 805,
        },
        {
            "left": 467,
            "right": 944,
            "top": 94,
            "bottom": 996,
        },
    ]
    bounding_boxes = set_bounding_boxes(bounds)


    # Set up report variables, and Generate report
    report = get_report(image, salient_characters=bounding_boxes)
    # report = get_report(image)
    report.image = image

    # Display images that represent image: blur profile and color palette
    report.generate_color_palette_image()
    report.bounding_boxes = bounding_boxes
    report.display_all()
    json = report.to_json()
    print(json)


if __name__ == "__main__":
    run_demonstration()
