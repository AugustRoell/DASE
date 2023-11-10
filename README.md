# DASE - Dimensionally reduced Acoustic Schrödinger Equation formalism solver

DASE is a software package used to simulate radially symmetric optical Schrödinger equation systems. 
The functionality is developed for acoustic modes inside of an HBAR device.

## File structure
### Scripts
- `DASE_example.py` contains an example of how to simulate a single HBAR design
- `DASE_SMQ.py` contains an example of how to simulate the SMQ (Single Mode Quality) values of a 2D design space
- `DASE_reverse_engineer.py` TODO. Will contain an example of how to construct and simulate the optimal cavity design for a single mode addressing scheme with a given forcing function profile

### DASE_Libraries
- `config.py` contains user-determined simulation and material parameters.
- `DASE.py` contains the core functionality of the simulation methods.
- `materials.py` contains the material specific constants, stored as dictionaries for fast lookup times.
- `profiles.py` contains different dome profiles as one-dimensional functions.

### Running a simulation
Run the desired script, either via the terminal or using your favourite source-code editor.
  
## Dependencies 
This code can require python 3.6+ and depends on the following packages:
- [matplotlib](https://matplotlib.org/)
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)