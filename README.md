# acwv-solver

acwv-solver is a simple acoustic waves simulator that implements a Finite Difference Method to compute the solution of the [Wave Equation](#wave-equation) for a given velocity field. It also uses MPI for parallelization. 

## Instructions

### Installation

#### Compiling from source

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

#### Download the binary

You can download the available binary for Linux from the [releases page](https://github.com/antoniomedeiros1/acwv-solver/releases).

### Usage

First, you need to generate the velocity field input file. The velocity field is a binary file that contains the speed of sound at each point of the domain. You can the 'generate_input.ipynb' Jupyter notebook to generate an example velocity field. 

To run the code, you can use the following command:

```sh
mpirun -np <number_of_processes> ./acwv-solver <input_file> <output_file> <number_of_steps> <dt>
```

where:
- `<number_of_processes>` is the number of MPI processes to use.
- `<input_file>` is the path to the velocity field input file.
- `<output_file>` is the path to the output file.
- `<number_of_steps>` is the number of time steps to simulate.
- `<dt>` is the time step.

## Theoretical Background

### Acoustic Waves

Acoustic waves are a type of mechanical wave that propagates through a medium as a result of the vibration of the particles of the medium (e.g sound and seismic waves). The speed of sound in a medium depends on the properties of the medium, such as density and elasticity. The speed of sound in air at room temperature is approximately 343 m/s.

### Wave Equation

The wave equation is a partial differential equation that describes the behavior of waves in a medium. The wave equation for acoustic waves is given by:

![Wave Equation](https://latex.codecogs.com/svg.latex?%5Cfrac%7B%5Cpartial%5E2%20p%7D%7B%5Cpartial%20t%5E2%7D%20%3D%20v%5E2%20%5Cnabla%5E2%20p)

where:
- ![p](https://latex.codecogs.com/svg.latex?p) is the pressure field,
- ![v](https://latex.codecogs.com/svg.latex?v) is the speed of sound in the medium,
- ![t](https://latex.codecogs.com/svg.latex?t) is time, and 
- ![%5Cnabla%5E2](https://latex.codecogs.com/svg.latex?%5Cnabla%5E2) is the Laplacian operator.

### Finite Difference Method

The Finite Difference Method is a numerical technique used to solve partial differential equations. The basic idea is to discretize the domain of the problem into a grid and approximate the derivatives in the differential equation using finite differences. 

For this project, we used a fourth-order finite difference scheme to approximate the Laplacian operator and a second-order finite difference scheme to approximate the time derivative ([more details](https://tattered-sleet-912.notion.site/2-1-M-todo-das-Diferen-as-Finitas-5fc05140011841648f7478e83fd5ffb2?pvs=4)).

### Boundary Conditions

In this project, we used the following boundary conditions:
- Non-reflecting Reynolds boundary conditions at the edges of the domain to simulate an open boundary.
- Damping boundary conditions at the edges of the domain to absorb the waves.

[More details](https://tattered-sleet-912.notion.site/2-3-Condi-es-de-Contorno-32e0fe6b394f49f3ab28a4b0122c3a29).

