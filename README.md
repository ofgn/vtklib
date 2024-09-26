# vtklib

vtklib is a Fortran 90 module designed for reading, writing, and processing VTK (Visualization Toolkit) files. It supports unstructured grids, scalar and vector data, and integrates seamlessly with TetGen files.

## Features

- **Unstructured Grid Handling**: Read and write VTK unstructured grids.
- **Data Processing**: Manage scalar and vector data associated with points and cells.
- **TetGen Integration**: Import TetGen `.node` and `.ele` files.
- **Data Masking**: Blank points or cells based on custom masks.
- **Flexible Output**: Export VTK files in both binary and ASCII formats.

## Installation

1. **Clone the Repository**
    ```bash
    git clone https://github.com/yourusername/vtklib.git
    cd vtklib
    ```

2. **Compile the Module**
    Ensure you have a Fortran compiler installed (e.g., `gfortran`).

    ```bash
    gfortran -c vtklib.f90
    gfortran -o vtklib_example main.f90 vtklib.o
    ```

## Usage

### **Basic Example**

```fortran
program example_usage
    use vtklib
    implicit none

    type(vtk_u_grid) :: mesh
    type(vtk_data) :: point_data, cell_data

    ! Read TetGen node and element files
    call mesh%read_tetgen_node("path/to/ex.2.node")
    call mesh%read_tetgen_ele("path/to/ex.2.ele")

    ! (Optional) Read scalar data
    call read_scalar_real64("path/to/two_blocks.sigma", cell_data)

    ! (Optional) Apply a cell mask
    call mesh%blank_on_cell_mask(cell_data%scalar_integer(1, :) .ne. 4)

    ! Write to VTK file
    call mesh%write_vtk_legacy("output/ex_2.vtk")
end program example_usage
