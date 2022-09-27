# peptoid-sequence-tools

## Description
This Python package was developed by the Knight Lab at UNC Chapel Hill to analyze the sequence of a peptoid from mass peaks in a MALDI-TOF spectrum. {more description}

## Installation Instructions
A version of Python >= 3.7 is required to use this package. We recommend using [Anaconda](https://www.anaconda.com) to install Python for those new to Python.
1. Clone this repository (requires git to be installed) or download the package from GitHub as a Zip file and unzip.


      git clone https://github.com/UNC-Knight-Lab/peptoid-sequence-tools.git

2. Open the terminal (MacOS) or Command Line (Windows) and install the package using pip. This command will install this package to your Python environment.


      pip install {path_to_this_package}

That's it!

## How to use
The peptoid sequence matching tool can be run as a Python function or from the command line terminal.
The script may prompt the user to make decisions in a tie-breaker scenario.

Data will be output as a single Excel worksheet in the output folder.

### To run from terminal:
    
    seq_match -i {path_to_input_folder} -o {path_to_output_folder}

If an input folder or output folder is not specified, the current working directory will be used for both.
Use the help `-h` tag to see more options.

### To run in Python:
Import the Python function

    from seq_match import seq_match
    seq_match(input_folder, output_folder)