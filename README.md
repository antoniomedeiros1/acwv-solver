# acwv-solver

acwv-solver is a acoustic waves in heterogeneous medium simulator that implements a Stencil Algorithm, using a Finite Difference Method, to compute an approximation for the [Wave Equation](#wave-equation) for a given velocity field and time span, outputing the results in VTK Image Data format, which can be visualized in ParaView. *The code is parallelized using MPI and PETSc (WIP)*.

|2D Output Example|3D Output Example (WIP)|
|:-------------------------:|:-------------------------:|
|![](assets/wave.png)|![](assets/wave3d.png)|

# Instructions

## Installation

### Download the binary

You can download the available binary for Linux from the [releases page](https://github.com/antoniomedeiros1/acwv-solver/releases).

### Compiling from source

#### Dependencies

To compile the code, you need to have the following dependencies installed:
- CMake
- VTK (libvtk7-dev)
- ParaView (optional, for visualization)
- *MPI (e.g OpenMPI) - WIP*
- *PETSc - WIP*

#### Compiling

To clone the repository, run the command:

```sh
git clone https://github.com/antoniomedeiros1/acwv-solver.git
```

To compile the code, first you need to create the build directory:

```sh
mkdir build
cd build
```

Then, you can run the following command to generate the Makefile:

```sh
cmake ..
```

Finally, you can compile the code by running:

```sh
make
```

## Usage

### Input File

First, you need to generate the velocity field input file. The velocity field is a VTK Image Data file that contains the speed of sound at each point of the domain, as well as the grid spacing and size. You can use the 'generate_input.ipynb' Jupyter notebook as reference on how to generate the velocity field. 

### Running the code

To run the code, you can use the following command:

```sh
./acwv-solver <input_file> <output_folder> <number_of_steps> <dt> [number_of_frames]
```

where:
- `<input_file>` is the path to the velocity field input file.
- `<output_folder>` is the path to the folder where the output files will be saved.
- `<number_of_steps>` is the number of time steps to simulate.
- `<dt>` is the time step.
- `[number_of_frames]` is the number of frames to save. If not provided, default is 40.

# Theoretical Background

## Acoustic Waves

Acoustic waves are a type of mechanical wave that propagates through a medium as a result of the vibration of the particles of the medium (e.g sound and seismic waves). The speed of sound in a medium depends on the properties of the medium, such as density and elasticity. The speed of sound in air at room temperature is approximately 343 m/s.

## Wave Equation

The wave equation is a partial differential equation that describes the behavior of waves in a medium. The wave equation for acoustic waves is given by:

![Wave Equation](https://latex.codecogs.com/svg.latex?%5Cfrac%7B%5Cpartial%5E2%20p%7D%7B%5Cpartial%20t%5E2%7D%20%3D%20v%5E2%20%5Cnabla%5E2%20p)

where:
- ![p](https://latex.codecogs.com/svg.latex?p) is the pressure field,
- ![v](https://latex.codecogs.com/svg.latex?v) is the speed of sound in the medium,
- ![t](https://latex.codecogs.com/svg.latex?t) is time, and 
- ![%5Cnabla%5E2](https://latex.codecogs.com/svg.latex?%5Cnabla%5E2) is the Laplacian operator.

## Finite Difference Method

The Finite Difference Method is a numerical technique used to solve partial differential equations. The basic idea is to discretize the domain of the problem into a grid and approximate the derivatives in the differential equation using finite differences. 

For this project, we used a fourth-order finite difference scheme to approximate the Laplacian operator and a second-order finite difference scheme to approximate the time derivative ([more details](https://tattered-sleet-912.notion.site/2-1-M-todo-das-Diferen-as-Finitas-5fc05140011841648f7478e83fd5ffb2?pvs=4)).

### Boundary Conditions

In this project, we used the following boundary conditions:
- Non-reflecting Reynolds boundary conditions at the edges of the domain to simulate an open boundary.
- Cerjan Absorbing boundary conditions at the edges of the domain to absorb the waves.

## Stability Condition

The stability of the numerical solution is an important aspect of the Finite Difference Method. The Courant-Friedrichs-Lewy (CFL) condition is a necessary condition for the stability of the numerical solution. The CFL condition for the wave equation is given by:

![CFL Condition](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bv%20%5CDelta%20t%7D%7B%5CDelta%20x%7D%20%5Cleq%201)

where:
- ![v](https://latex.codecogs.com/svg.latex?v) is the speed of sound in the medium,
- ![%5CDelta%20t](https://latex.codecogs.com/svg.latex?%5CDelta%20t) is the time step, and
- ![%5CDelta%20x](https://latex.codecogs.com/svg.latex?%5CDelta%20x) is the grid spacing.

## Publications

- [Modelagem Acústica por Diferenças Finitas utilizando OpenMP (ERAD-RJ 2021)](https://sol.sbc.org.br/index.php/eradrj/article/view/18561)
