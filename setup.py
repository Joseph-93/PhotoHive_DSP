from setuptools import setup, find_packages

setup(
    name='PhotoHive_DSP_lib',  # Replace 'mylibrary' with the actual name of your library
    version='1.0.0',  # The current version of your library
    author='Joshua Matheson',  # Your name or your organization's name
    author_email='joshuajosephmatheson@gmail.com',  # Your email or your organization's email
    description='PhotoHive\'s library of digital signal processing for image feature extraction.',  # A short description of the library
    long_description=open('README.md').read(),  # A long description from README.md
    long_description_content_type='text/markdown',  # Specifies that the long description is in Markdown
    url='https://github.com/Joseph-93/PhotoHive-DSP',  # URL to the repository or website of the library
    packages=find_packages(),  # Automatically find and include all packages
    include_package_data=True,  # Includes package data defined in MANIFEST.in
    install_requires=[
        'numpy',  # Required dependencies for your library
        'Pillow',  # PIL Fork for image processing
        # Add other dependencies as needed
    ],
    classifiers=[  # Classifiers help users find your project by categorizing it
        'Programming Language :: Python :: 3',
        'Operating System :: Linux',
    ],
    python_requires='>=3.8',  # Specify the minimum version of Python required
)
