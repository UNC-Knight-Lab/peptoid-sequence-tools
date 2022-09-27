# peptoid-sequence-tools

## Description
This Python package was developed by the Knight Lab at UNC Chapel Hill to analyze the sequence of a peptoid from mass peaks in a MALDI-TOF spectrum. {more description}

## Installation Instructions
A version of Python >= 3.7 is required to use this package. We recommend using [Anaconda](https://www.anaconda.com) to install Python for those new to Python.
1. Open the terminal (MacOS) or Command Prompt (Windows).
2. Download the package by either:
   1. Download the zip from GitHub (Code -> Download ZIP). Unzip the package somewhere (note the extraction path). The extracted package can be deleted after installation.
   2. Clone this repository (requires git to be installed) with:


      git clone https://github.com/UNC-Knight-Lab/peptoid-sequence-tools.git

3. Install the package using pip. This command will install this package to your Python environment.
    The package path should be `./seq_match` if cloned using git. Otherwise, replace `./seq_match` with the path to the unzipped package.


      pip install ./seq_match

That's it!

## How to use
The peptoid sequence matching tool can be run as a Python function or from the command line terminal.
The script may prompt the user to make decisions in a tie-breaker scenario.

Data will be output as a single Excel worksheet in the output folder.

### To run from terminal:
    
    seq_match -i "/path/to/input_folder" -o "/path/to/output_folder"

Instead of specifying an input or output folder, you can also navigate to your data input folder in the terminal and run the script
The current working directory will be used as default.
Use the help `-h` tag to see more options.

### To run in Python:
In a Python environment, import the Python function:

    from seq_match import seq_match
    seq_match(input_folder, output_folder)