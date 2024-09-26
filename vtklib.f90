! -------------------------------------------------------------------------------------------------------
! @file vtklib.f90
! @brief Module for reading, writing and processing VTK files.
! @author ofgn
! @date 2024-07-01
! -------------------------------------------------------------------------------------------------------
module vtklib
    implicit none

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Type definition for unstructured grids (u_grid). Handles points, 
    ! cells, and connectivity information.
    ! ---------------------------------------------------------------------------------------------------
    type :: vtk_u_grid
        integer :: n_points                                                     ! Number of points in the grid
        integer :: n_cells                                                      ! Number of cells in the grid
        real(8), allocatable :: points(:, :)                                    ! 3D coordinates of points
        integer, allocatable :: points_per_cell(:)                              ! Number of points in each cell
        integer, allocatable :: cell_connectivity(:)                            ! Point indices for each cell
        integer, allocatable :: cell_types(:)                                   ! VTK cell types for each cell
    contains
        procedure :: read_tetgen_node                                           ! Read TetGen node file
        procedure :: read_tetgen_ele                                            ! Read TetGen element file
        procedure :: write_vtk_legacy                                           ! Write legacy VTK file in binary or ASCII
        procedure :: blank_on_point_mask                                        ! Blank points and cells based on a point mask
        procedure :: blank_on_cell_mask                                         ! Blank points and cells based on a cell mask
    end type vtk_u_grid

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Type definition for handling scalar and vector data in various 
    ! formats (integer, real32, real64) associated with points or cells.
    ! ---------------------------------------------------------------------------------------------------
    type :: vtk_data
        integer :: n = 0                                                        ! Number of points or cells
        integer, allocatable :: scalar_integer(:, :)                            ! Integer scalar data
        real(4), allocatable :: scalar_real32(:, :)                             ! Real(4) scalar data
        real(8), allocatable :: scalar_real64(:, :)                             ! Real(8) scalar data
        character(len=32), allocatable :: scalar_integer_labels(:)              ! Labels for integer scalar data
        character(len=32), allocatable :: scalar_real32_labels(:)               ! Labels for real(4) scalar data
        character(len=32), allocatable :: scalar_real64_labels(:)               ! Labels for real(8) scalar data
        integer, allocatable :: scalar_integer_components(:)                    ! Components of integer scalar data
        integer, allocatable :: scalar_real32_components(:)                     ! Components of real(4) scalar data
        integer, allocatable :: scalar_real64_components(:)                     ! Components of real(8) scalar data
        integer, allocatable :: vector_integer(:, :)                            ! Integer vector data
        real(4), allocatable :: vector_real32(:, :)                             ! Real(4) vector data
        real(8), allocatable :: vector_real64(:, :)                             ! Real(8) vector data
        character(len=32), allocatable :: vector_integer_labels(:)              ! Labels for integer vector data
        character(len=32), allocatable :: vector_real32_labels(:)               ! Labels for real(4) vector data
        character(len=32), allocatable :: vector_real64_labels(:)               ! Labels for real(8) vector data
    contains
        procedure :: info                                                       ! Print data structure info
        procedure :: add_scalar_real64                                          ! Add scalar data in real(8) format
    end type vtk_data
contains

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Read a TetGen node file and extract the node coordinates.
    ! @param[in] file_path Path to the TetGen node file.
    ! @param[out] point_data Data structure to store point information (optional).
    ! ---------------------------------------------------------------------------------------------------
    subroutine read_tetgen_node(self, file_path, point_data)
        implicit none
        class(vtk_u_grid) :: self                                               ! The unstructured grid to store node data
        character(len=*) :: file_path                                           ! Path to the TetGen node file
        type(vtk_data), optional :: point_data                                  ! Data structure for storing point info (optional)

        integer :: n_dimensions, n_scalars, id_flag                             ! Dimensions, scalars, and boundary ID flag
        integer :: point_index                                                  ! Index for the point data
        integer :: i, unit, io_status                                           ! Loop index, unit ID, I/O status
        real(8) :: start_time, end_time                                         ! Timing variables

        call cpu_time(start_time)
        call report("------------------------[  VTKLIB FILE READ    ]------------------------"//new_line('A'))

        open (newunit=unit, file=file_path, status="old", action="read", iostat=io_status)
        if (io_status .ne. 0) then
            call report("- Status:                       (ERROR) Failed to open file: "//trim(file_path), is_error=.true.)
            stop
        end if

        call report("- File:                         "//trim(file_path))

        ! Read header with number of points, dimensions, number of scalars, and id flag
        read (unit, *, iostat=io_status) self%n_points, n_dimensions, n_scalars, id_flag
        if (io_status .ne. 0) then
            call report("- Status:                       (ERROR) Error parsing file header", is_error=.true.)
            stop
        end if

        ! Validate dimensions and report
        if (n_dimensions .eq. 2) then
            call report("- Description:                  "//"2D points (Triangle)")
        else if (n_dimensions .eq. 3) then
            call report("- Description:                  "//"3D points (TetGen)")
        else
            call report("- Status:                       (ERROR) Invalid number of dimensions", is_error=.true.)
            stop
        end if

        ! Allocate memory for points
        allocate(self%points(3, self%n_points))

        ! Allocate memory for scalars and boundary ID if present
        if (present(point_data) .and. n_scalars .gt. 0) then
            allocate(point_data%scalar_real64(n_scalars, self%n_points))
            allocate(point_data%scalar_real64_labels(n_scalars))
            allocate(point_data%scalar_real64_components(n_scalars))
            do i = 1, n_scalars
                write(point_data%scalar_real64_labels(i), "(A, I0, A)") "Scalar ", i, "[double]"
            end do
            point_data%scalar_real64_components = 1
        end if

        if (present(point_data) .and. id_flag .gt. 0) then
            allocate(point_data%scalar_integer(1, self%n_points))
            allocate(point_data%scalar_integer_labels(1))
            point_data%scalar_integer_labels(1) = "Boundary Id"
            point_data%scalar_integer_components = 1
        end if

        ! Read point data directly into the appropriate arrays
        do i = 1, self%n_points
            if (present(point_data) .and. (n_scalars .gt. 0) .and. (id_flag .gt. 0)) then
                read (unit, *, iostat=io_status) point_index, self%points(1:n_dimensions, i), &
                    point_data%scalar_real64(:, i), point_data%scalar_integer(1, i)
                point_data%n = self%n_points
            else if (present(point_data) .and. (n_scalars .gt. 0)) then
                read (unit, *, iostat=io_status) point_index, self%points(1:n_dimensions, i), &
                    point_data%scalar_real64(:, i)
                point_data%n = self%n_points
            else if (present(point_data) .and. (id_flag .gt. 0)) then
                read (unit, *, iostat=io_status) point_index, self%points(1:n_dimensions, i), &
                    point_data%scalar_integer(1, i)
                point_data%n = self%n_points
            else
                read (unit, *, iostat=io_status) point_index, self%points(1:n_dimensions, i)
            end if

            ! If 2D points, set the z-coordinate to 0
            if (n_dimensions .eq. 2) then
                self%points(3, i) = 0.0d0
            end if

            if (io_status .ne. 0) then
                call report("- Status:                       (ERROR) Error reading point: "//trim(itoa(i)), is_error=.true.)
                stop
            end if
        end do

        close (unit)

        call cpu_time(end_time)
        call report("- Status:                       Completed in "&
            //trim(rtoa(end_time - start_time, decimal_places=4))//" seconds")
        call report("------------------------------------------------------------------------")
    end subroutine read_tetgen_node

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Read a TetGen element file (.ele) and extract connectivity.
    ! @param[in] file_path Path to the TetGen element file.
    ! @param[out] cell_data Data structure to store cell information (optional).
    ! ---------------------------------------------------------------------------------------------------
    subroutine read_tetgen_ele(self, file_path, cell_data)
        implicit none
        class(vtk_u_grid) :: self                                               ! The unstructured grid to store element data
        character(len=*), intent(in) :: file_path                               ! Path to the TetGen element file
        type(vtk_data), optional :: cell_data                                   ! Data structure for storing cell data (optional)

        integer :: i, unit, io_status, points_per_cell, id_flag                 ! Loop index, unit ID, I/O status
        integer :: cell_i, dummy                                                ! Element index
        real(8) :: start_time, end_time                                         ! Timing variables

        call cpu_time(start_time)
        open (newunit=unit, file=file_path, status="old", action="read", iostat=io_status)
        if (io_status .ne. 0) then
            call report("- Status:                       (ERROR) Failed to open file: "//trim(file_path), is_error=.true.)
            stop
        end if

        call report("------------------------[  VTKLIB FILE READ    ]------------------------"//new_line('A'))
        call report("- File:                         "//trim(file_path))
        ! call report("- Status:                       "//"Started")

        read (unit, *, iostat=io_status) self%n_cells, points_per_cell, id_flag
        if (io_status .ne. 0) then
            call report("- Status:                       (ERROR) Error parsing file header", is_error=.true.)
            stop
        end if

        allocate (self%cell_connectivity(self%n_cells * points_per_cell))
        allocate (self%points_per_cell(self%n_cells))
        allocate (self%cell_types(self%n_cells))

        if (points_per_cell .eq. 3) then
            self%cell_types = 5
            call report("- Description:                  "//"Unstructured grid (Triangle)")
        else if (points_per_cell .eq. 4) then
            self%cell_types = 10
            call report("- Description:                  "//"Unstructured grid (TetGen, linear)")
        else if (points_per_cell .eq. 10) then
            self%cell_types = 24
            call report("- Description:                  "//"Unstructured grid (TetGen, quadratic)")
        else
            call report("- Status:                       (ERROR) Invalid number of points per cell", is_error=.true.)
            stop
        end if

        self%points_per_cell = points_per_cell

        if (present(cell_data) .and. id_flag > 0) then
            allocate (cell_data%scalar_integer(id_flag, self%n_cells))
            allocate (cell_data%scalar_integer_labels(id_flag))
            allocate (cell_data%scalar_integer_components(id_flag))
            cell_data%scalar_integer_labels(1) = "Id"
            cell_data%scalar_integer_components = 1
            cell_data%n = self%n_cells
        end if

        if (id_flag > 0) then
            if (present(cell_data)) then
                read (unit, *, iostat=io_status) &
                    (cell_i, self%cell_connectivity((i - 1) * points_per_cell + 1:i * points_per_cell), &
                    cell_data%scalar_integer(:, i), i=1, self%n_cells)
            else 
                read (unit, *, iostat=io_status) &
                    (cell_i, self%cell_connectivity((i - 1) * points_per_cell + 1:i * points_per_cell), &
                    dummy, i=1, self%n_cells)
            end if
        else
            read (unit, *, iostat=io_status) &
                (cell_i, self%cell_connectivity((i - 1) * points_per_cell + 1:i * points_per_cell), i=1, self%n_cells)
        end if

        if (io_status .ne. 0) then
            call report("- Status:                       (ERROR) Error parsing file data", is_error=.true.)
            stop
        end if

        self%cell_connectivity = self%cell_connectivity - 1
        close (unit)

        call cpu_time(end_time)
        call report("- Status:                       Completed in "&
            //trim(rtoa(end_time - start_time, decimal_places=4))//" seconds")
        call report("------------------------------------------------------------------------")
    end subroutine read_tetgen_ele

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Write an unstructured grid to a VTK legacy file in either binary or ASCII format.
    ! @param[in] file_path Path to the output VTK file.
    ! @param[in] point_data Data associated with points (optional).
    ! @param[in] cell_data Data associated with cells (optional).
    ! @param[in] ascii Logical flag to write in ASCII format (default: false).
    ! ---------------------------------------------------------------------------------------------------
    subroutine write_vtk_legacy(self, file_path, point_data, cell_data, ascii)
        implicit none
        class(vtk_u_grid), intent(in) :: self                                   ! Unstructured grid data structure
        character(len=*), intent(in) :: file_path                               ! Path to the output VTK file
        type(vtk_data), optional, intent(in) :: point_data                      ! Data structure for point data (optional)
        type(vtk_data), optional, intent(in) :: cell_data                       ! Data structure for cell data (optional)
        logical, intent(in), optional :: ascii                                  ! Write in ASCII format (default: false)

        integer, allocatable :: full_connectivity(:)                            ! Complete cell connectivity array
        integer :: i, j, unit, io_status                                        ! Loop variables, file unit, I/O status
        real(8) :: start_time, end_time                                         ! Timing variables
        character(len=1) :: lf                                                  ! Line feed character for file format
        logical :: binary                                                       ! Flag to determine binary or ASCII format

        call cpu_time(start_time)
        
        if (present(ascii)) then
            binary = .not. ascii
        else
            binary = .true.
        end if

        if (binary) then
            open (newunit=unit, file=file_path, status="replace", access="stream", &
                action="write", form="unformatted", convert="big_endian", iostat=io_status)
        else
            open (newunit=unit, file=file_path, status="replace", access="stream", &
                action="write", form="formatted", iostat=io_status)
        end if

        if (io_status .ne. 0) then
            call report("ERROR: Failed to open file: "//trim(file_path), is_error=.true.)
            stop
        end if

        call report("------------------------[  VTKLIB FILE WRITE   ]------------------------")
        call report("- File:                         "//trim(file_path))
            ! Determine file format (ASCII or binary)
        if (binary) then
            call report("- Format:                       "//"Binary")
        else
            call report("- Format:                       "//"ASCII")
        end if

        lf = char(10)

        ! Allocate the full connectivity array
        allocate (full_connectivity(size(self%cell_connectivity) + self%n_cells))
        j = 1
        do i = 1, self%n_cells
            full_connectivity(j) = self%points_per_cell(i)
            full_connectivity(j + 1:j + self%points_per_cell(i)) = &
                self%cell_connectivity((i - 1) * self%points_per_cell(i) + 1:i * self%points_per_cell(i))
            j = j + self%points_per_cell(i) + 1
        end do

        if (binary) then
            write (unit) "# vtk DataFile Version 3.0"//lf
            write (unit) "vtk"//lf
            write (unit) "BINARY"//lf
            write (unit) "DATASET UNSTRUCTURED_GRID"//lf
        else
            write (unit, "(a)") "# vtk DataFile Version 3.0"
            write (unit, "(a)") "vtk"
            write (unit, "(a)") "ASCII"
            write (unit, "(a)") "DATASET UNSTRUCTURED_GRID"
        end if

        if (binary) then
            write (unit) "POINTS "//trim(itoa(self%n_points))//" double"//lf
            write (unit) self%points
            write (unit) lf
        else
            write (unit, "(A, I0, A)") "POINTS ", self%n_points, " double"
            do i = 1, self%n_points
                write (unit, "(3E20.8)") self%points(:, i)
            end do
        end if

        if (binary) then
            write (unit) "CELLS "//trim(itoa(self%n_cells))//" "//trim(itoa(size(full_connectivity)))//lf
            write (unit) full_connectivity
            write (unit) lf
        else
            write (unit, "(A, I0, A, I0)") "CELLS ", self%n_cells  , " ", size(full_connectivity)
            do i = 1, size(full_connectivity)
                write (unit, "(I0)") full_connectivity(i)
            end do
        end if

        if (binary) then
            write (unit) "CELL_TYPES "//trim(itoa(self%n_cells))//lf
                write (unit) self%cell_types
                write (unit) lf
        else
            write (unit, "(A, I0)") "CELL_TYPES ", self%n_cells
            do i = 1, self%n_cells
                write (unit, "(I0)") self%cell_types(i)
            end do
        end if

        call cpu_time(end_time)
        close (unit)

        call report("- Status:                       Completed in "&
            //trim(rtoa(end_time - start_time, decimal_places=4))//" seconds")
        call report("------------------------------------------------------------------------")
    end subroutine write_vtk_legacy

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Mask and remove points from the grid based on the given point mask.
    ! @param[inout] self The unstructured grid to modify.
    ! @param[in] point_mask Logical mask to determine which points to keep.
    ! ---------------------------------------------------------------------------------------------------
    subroutine blank_on_point_mask(self, point_mask)
        implicit none
        class(vtk_u_grid), intent(inout) :: self                                ! The unstructured grid to modify
        logical, intent(in) :: point_mask(:)                                    ! Logical mask for filtering points

        integer :: i, j, k, new_n_points, new_n_cells                           ! Counters for points and cells
        integer, allocatable :: new_cell_connectivity(:)                        ! New cell connectivity array
        integer, allocatable :: new_points_per_cell(:)                          ! New points-per-cell array
        integer, allocatable :: new_cell_types(:)                               ! New cell types array
        integer, allocatable :: point_map(:)                                    ! Mapping from old points to new points
        real(8), allocatable :: new_points(:, :)                                ! New point coordinates array
        real(8) :: start_time, end_time                                         ! Timing variables
        logical, allocatable :: cell_mask(:), connectivity_mask(:)              ! Masks for cells and connectivity

        call cpu_time(start_time)
        call report("------------------------[    VTKLIB BLANK     ]-------------------------"//new_line('A'))
        call report("- Status:                       Started")
        call report("- Description:                  "//"Blanking grid based on a point mask") 

        if (size(point_mask) .ne. self%n_points) then
            call report("ERROR: Mask size does not match number of points.", is_error=.true.)
            stop
        end if

        allocate(cell_mask(self%n_cells), connectivity_mask(size(self%cell_connectivity)))
        cell_mask = .true.
        connectivity_mask = .true.

        j = 1
        do i = 1, self%n_cells
            do k = 0, self%points_per_cell(i) - 1
                if (.not. point_mask(self%cell_connectivity(j + k) + 1)) then
                    cell_mask(i) = .false.
                    connectivity_mask(j:j + self%points_per_cell(i) - 1) = .false.
                    exit
                end if
            end do
            j = j + self%points_per_cell(i)
        end do

        allocate(point_map(self%n_points))
        point_map = -1
        j = 1
        do i = 1, self%n_points
            if (point_mask(i)) then
                point_map(i) = j
                j = j + 1
            end if
        end do

        new_n_points = count(point_mask)
        new_n_cells = count(cell_mask)

        allocate(new_points(3, new_n_points))
        new_points = self%points(:, pack([(i, i=1, self%n_points)], point_mask))

        new_cell_connectivity = pack(self%cell_connectivity, connectivity_mask)
        new_points_per_cell = pack(self%points_per_cell, cell_mask)
        new_cell_types = pack(self%cell_types, cell_mask)

        new_cell_connectivity = point_map(new_cell_connectivity + 1) - 1
        deallocate(cell_mask, connectivity_mask)

        self%points = new_points
        self%cell_connectivity = new_cell_connectivity
        self%points_per_cell = new_points_per_cell
        self%cell_types = new_cell_types

        call report("- Points prior:                 "//trim(itoa(self%n_points)))
        self%n_points = new_n_points
        call report("- Points after:                 "//trim(itoa(new_n_points)))
        call report("- Cells prior:                  "//trim(itoa(self%n_cells)))
        self%n_cells = new_n_cells
        call report("- Cells after:                  "//trim(itoa(new_n_cells)))
        call cpu_time(end_time)
        call report("- Status:                       Completed in "&
            //trim(rtoa(end_time - start_time, decimal_places=4))//" seconds")
        call report("------------------------------------------------------------------------")
    end subroutine blank_on_point_mask

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Mask and remove cells from the grid based on the given cell mask.
    ! @param[inout] self The unstructured grid to modify.
    ! @param[in] cell_mask Logical mask to determine which cells to keep.
    ! ---------------------------------------------------------------------------------------------------
    subroutine blank_on_cell_mask(self, cell_mask)
        implicit none
        class(vtk_u_grid), intent(inout) :: self                                ! The unstructured grid to modify
        logical, intent(in) :: cell_mask(:)                                     ! Logical mask for filtering cells

        integer :: i, j, new_n_points, new_n_cells                              ! Counters for points and cells
        integer, allocatable :: new_cell_connectivity(:)                        ! New cell connectivity array
        integer, allocatable :: new_points_per_cell(:)                          ! New points-per-cell array
        integer, allocatable :: new_cell_types(:)                               ! New cell types array
        integer, allocatable :: point_map(:)                                    ! Mapping from old points to new points
        real(8), allocatable :: new_points(:, :)                                ! New point coordinates array
        real(8) :: start_time, end_time                                         ! Timing variables
        logical, allocatable :: point_mask(:), connectivity_mask(:)             ! Masks for points and connectivity

        call cpu_time(start_time)
        call report("------------------------[    VTKLIB BLANK     ]-------------------------"//new_line('A'))
        call report("- Description:                  "//"Blanking grid based on a cell mask") 
        ! call report("- Status:                       "//"Started")

        if (size(cell_mask) .ne. self%n_cells) then
            call report("ERROR: Mask size does not match number of cells.", is_error=.true.)
            stop
        end if

        allocate(point_mask(self%n_points), connectivity_mask(size(self%cell_connectivity)))
        point_mask = .false.
        connectivity_mask = .false.

        j = 1
        do i = 1, self%n_cells
            if (cell_mask(i)) then
                point_mask(self%cell_connectivity(j:j + self%points_per_cell(i) - 1) + 1) = .true.
                connectivity_mask(j:j + self%points_per_cell(i) - 1) = .true.
            end if
            j = j + self%points_per_cell(i)
        end do

        allocate(point_map(self%n_points))
        point_map = -1
        j = 1
        do i = 1, self%n_points
            if (point_mask(i)) then
                point_map(i) = j
                j = j + 1
            end if
        end do

        new_n_points = count(point_mask)
        new_n_cells = count(cell_mask)

        allocate(new_points(3, new_n_points))
        new_points = self%points(:, pack([(i, i=1, self%n_points)], point_mask))

        new_cell_connectivity = pack(self%cell_connectivity, connectivity_mask)
        new_points_per_cell = pack(self%points_per_cell, cell_mask)
        new_cell_types = pack(self%cell_types, cell_mask)

        new_cell_connectivity = point_map(new_cell_connectivity + 1) - 1
        deallocate(point_mask, connectivity_mask)

        
        
        self%points = new_points
        self%cell_connectivity = new_cell_connectivity
        self%points_per_cell = new_points_per_cell
        self%cell_types = new_cell_types

        call report("- Points prior:                 "//trim(itoa(self%n_points)))
        self%n_points = new_n_points
        call report("- Points after:                 "//trim(itoa(new_n_points)))
        call report("- Cells prior:                  "//trim(itoa(self%n_cells)))
        self%n_cells = new_n_cells
        call report("- Cells after:                  "//trim(itoa(new_n_cells)))
        call cpu_time(end_time)
        call report("- Status:                       Completed in "&
            //trim(rtoa(end_time - start_time, decimal_places=4))//" seconds")
        call report("------------------------------------------------------------------------")
    end subroutine blank_on_cell_mask

    ! ===================================================================================================
    ! SECTION: Data Array (vtk_data) Procedures
    ! ===================================================================================================

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Report the contents of the data structure.
    ! @param[in] self The data structure to report.
    ! ---------------------------------------------------------------------------------------------------
    subroutine info(self)
        implicit none
        class(vtk_data), intent(in) :: self                                     ! Data structure to be reported
        integer :: i                                                            ! Loop index

        call report("------------------------[     VTKLIB INFO      ]------------------------"//new_line('A'))
        call report("- Number of points/cells:       "//trim(itoa(self%n)))

        ! Report integer scalar data if allocated
        if (allocated(self%scalar_integer)) then
            call report("- Scalar fields (integer):      "//trim(itoa(size(self%scalar_integer_components, 1))))
            do i = 1, size(self%scalar_integer_labels)
                call report("  * Label:                      "//trim(self%scalar_integer_labels(i)))
                call report("  * Components:                 "//trim(itoa(self%scalar_integer_components(i))))
            end do
        else
            call report("- Scalar fields (integer):      "//"0 (unallocated)")
        end if

        ! Report real32 scalar data if allocated
        if (allocated(self%scalar_real32)) then
            call report("- Scalar fields (real32):       "//trim(itoa(size(self%scalar_real32_components, 1))))
            do i = 1, size(self%scalar_real32_labels)
                call report("  * Label:                      "//trim(self%scalar_real32_labels(i)))
                call report("  * Components:                 "//trim(itoa(self%scalar_real32_components(i))))
            end do
        else
            call report("- Scalar fields (real32):       "//"0 (unallocated)")
        end if

        ! Report real64 scalar data if allocated
        if (allocated(self%scalar_real64)) then
            call report("- Scalar fields (real64):       "//trim(itoa(size(self%scalar_real64_components, 1))))
            do i = 1, size(self%scalar_real64_labels)
                call report("  * Label:                      "//trim(self%scalar_real64_labels(i)))
                call report("  * Components:                 "//trim(itoa(self%scalar_real64_components(i))))
            end do
        else
            call report("- Scalar fields (real64):       "//"0 (unallocated)")
        end if

        ! Report integer vector data if allocated
        if (allocated(self%vector_integer)) then
            call report("- Vector fields (integer):      "//trim(itoa(size(self%vector_integer_labels, 1))))
            do i = 1, size(self%vector_integer_labels)
                call report("  * Label:                      "//trim(self%vector_integer_labels(i)))
            end do
        else
            call report("- Vector fields (integer):      "//"0 (unallocated)")
        end if

        ! Report real32 vector data if allocated
        if (allocated(self%vector_real32)) then
            call report("- Vector fields (real32):       "//trim(itoa(size(self%vector_real32_labels, 1))))
            do i = 1, size(self%vector_real32_labels)
                call report("  * Label:                      "//trim(self%vector_real32_labels(i)))
            end do
        else
            call report("- Vector fields (real32):       "//"0 (unallocated)")
        end if

        ! Report real64 vector data if allocated
        if (allocated(self%vector_real64)) then
            call report("- Vector fields (real64):       "//trim(itoa(size(self%vector_real64_labels, 1))))
            do i = 1, size(self%vector_real64_labels)
                call report("  * Label:                      "//trim(self%vector_real64_labels(i)))
            end do
        else
            call report("- Vector fields (real64):       "//"0 (unallocated)")
        end if

        call report("------------------------------------------------------------------------")
    end subroutine info

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Add new real(8) scalar data to the data structure.
    ! @param[inout] self The data structure containing real(8) scalar data.
    ! @param[in] new_scalars_real64 The new real(8) scalar data to append.
    ! @param[in, optional] new_scalars_real64_labels Labels for the real(8) scalar data.
    ! @param[in, optional] new_scalars_real64_components Number of components for the new scalar.
    ! ---------------------------------------------------------------------------------------------------
    subroutine add_scalar_real64(self, new_scalars_real64, new_scalars_real64_labels, &
                                 new_scalars_real64_components)
        implicit none
        class(vtk_data), intent(inout) :: self                                  ! Data structure containing real(8) scalar data
        real(8), intent(in) :: new_scalars_real64(:, :)                         ! New real(8) scalar data to append
        character(len=*), intent(in), optional :: new_scalars_real64_labels     ! Labels for the scalar data
        integer, intent(in), optional :: new_scalars_real64_components          ! Number of components for the scalar data

        integer :: i
        integer :: n, n_new_scalars                                             ! Number of elements, number of new scalars
        integer, allocatable :: temp_components(:)                              ! Temporary array for components
        real(8), allocatable :: temp_scalars(:, :)                              ! Temporary array for scalar data
        character(len=32), allocatable :: temp_labels(:)                        ! Temporary array for scalar labels

        logical :: unallocated                                                  ! Flag to check allocation state of arrays

        n_new_scalars = size(new_scalars_real64, 1)
        n = size(new_scalars_real64, 2)

        unallocated = .not. (allocated(self%scalar_real64) .and. allocated(self%scalar_real64_labels) &
                             .and. allocated(self%scalar_real64_components))

        if (unallocated) then
            self%scalar_real64 = new_scalars_real64

            if (present(new_scalars_real64_components)) then
                self%scalar_real64_components = new_scalars_real64_components
            else
                allocate (self%scalar_real64_components(n_new_scalars))
                self%scalar_real64_components = 1
            end if

            if (present(new_scalars_real64_labels)) then
                self%scalar_real64_labels = new_scalars_real64_labels
            else
                allocate (self%scalar_real64_labels(size(self%scalar_real64_components)))
                do i = 1, size(self%scalar_real64_components)
                    write (self%scalar_real64_labels(i), "(A, I0, A)") "Scalar ", i, " [double]"
                end do
            end if
        else
            allocate (temp_scalars(size(self%scalar_real64, 1), n + n_new_scalars))
            temp_scalars(:, 1:n) = self%scalar_real64
            temp_scalars(:, n + 1:n + n_new_scalars) = new_scalars_real64
            deallocate (self%scalar_real64)
            self%scalar_real64 = temp_scalars

            allocate (temp_components(size(self%scalar_real64_components) + n_new_scalars))
            temp_components(1:size(self%scalar_real64_components)) = self%scalar_real64_components
            temp_components(size(self%scalar_real64_components) + 1:size(self%scalar_real64_components) + n_new_scalars) = 1
            deallocate (self%scalar_real64_components)
            self%scalar_real64_components = temp_components

            allocate (temp_labels(size(self%scalar_real64_labels) + n_new_scalars))
            temp_labels(1:size(self%scalar_real64_labels)) = self%scalar_real64_labels
            if (present(new_scalars_real64_labels)) then
                temp_labels(size(self%scalar_real64_labels) + 1:size(self%scalar_real64_labels) + n_new_scalars) = &
                    new_scalars_real64_labels
            else
                do i = 1, n_new_scalars
                    write (temp_labels(size(self%scalar_real64_labels) + i), "(A, I0, A)") "Scalar ", &
                        size(self%scalar_real64_labels) + i, " [double]"
                end do
            end if
            deallocate (self%scalar_real64_labels)
            self%scalar_real64_labels = temp_labels
        end if
    end subroutine add_scalar_real64

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Read scalar data from a file and add it to the cell data.
    ! @param[in] file_path Path to the scalar data file.
    ! @param[inout] geometric_data The data structure where the scalar data is stored.
    ! ---------------------------------------------------------------------------------------------------
    subroutine read_scalar_real64(file_path, geometric_data)
        implicit none
        character(len=*) :: file_path                                           ! Path to scalar data file
        type(vtk_data), intent(inout) :: geometric_data                         ! Data structure for storing scalar data
        integer :: i, unit, io_status, n, n_scalars                             ! Loop index, unit ID, I/O status, etc.
        real(8), allocatable :: scalar_real64(:, :)                             ! Array for real(8) scalar data
        real(8) :: start_time, end_time                                         ! Timing variables

        call cpu_time(start_time)

        open (newunit=unit, file=file_path, status="old", action="read", iostat=io_status)
        if (io_status .ne. 0) then
            call report("ERROR: Failed to open file: "//trim(file_path), is_error=.true.)
            stop
        end if

        call report("------------------------[   VTKLIB FILE READ   ]------------------------"//new_line('A'))
        call report("- File:                         "//trim(file_path))
        call report("- Description:                  "//"Scalar data (real64)") 
        ! call report("- Status:                       "//"Started")

        read (unit, *, iostat=io_status) n, n_scalars
        if (io_status .ne. 0) then
            call report("ERROR: Error parsing file header.", is_error=.true.)
            stop
        end if

        if (n .ne. geometric_data%n) then
            call report("ERROR: Mismatched number of points/cells.", is_error=.true.)
            stop
        end if

        allocate (scalar_real64(n_scalars, n))
        read (unit, *, iostat=io_status) (scalar_real64(1:n_scalars, i), i=1, n)
        if (io_status .ne. 0) then
            call report("ERROR: Error parsing file data.", is_error=.true.)
            stop
        end if

        close (unit)    
        call geometric_data%add_scalar_real64(scalar_real64)
        deallocate (scalar_real64)

        call cpu_time(end_time)
        call report("- Status:                       Completed in "&
            //trim(rtoa(end_time - start_time, decimal_places=4))//" seconds")
        call report("------------------------------------------------------------------------")
    end subroutine read_scalar_real64

    ! ===================================================================================================
    ! SECTION: General Procedures
    ! ===================================================================================================

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Log a message to stdout or stderr.
    ! @param[in] message The message to log.
    ! @param[in, optional] is_error Flag to log to stderr if present and true.
    ! ---------------------------------------------------------------------------------------------------
    subroutine report(message, is_error)
        use iso_fortran_env, only: output_unit, error_unit
        implicit none

        character(len=*), intent(in) :: message                         ! Message to log
        logical, intent(in), optional :: is_error                       ! Flag for logging to stderr (optional)

        if (present(is_error)) then
            if (is_error) then
                write (error_unit, "(a)") trim(message)                 ! Log to stderr
            else
                write (output_unit, "(a)") trim(message)                ! Log to stdout
            end if
        else
            write (output_unit, "(a)") trim(message)                    ! Default logging to stdout
        end if
    end subroutine report

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Convert an integer to a character string.
    ! @param[in] value The integer value to convert.
    ! @return The character string representation of the integer.
    ! ---------------------------------------------------------------------------------------------------
    function itoa(value) result(str)
        implicit none
        integer, intent(in) :: value                          ! Integer value to be converted
        character(len=32) :: str                              ! Resulting string (maximum length of 32 characters)
    
        ! Convert the integer to a string
        write (str, "(I0)") value                             ! 'I0' ensures no leading spaces in the result
    end function itoa

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Convert a real number to a character string.
    ! @param[in] value The real value to convert.
    ! @param[in, optional] decimal_places Number of decimal places to include in the output (default is 2).
    ! @return The character string representation of the real number.
    ! ---------------------------------------------------------------------------------------------------
    function rtoa(value, decimal_places) result(str)
        implicit none
        real(8), intent(in) :: value                         ! Real value to be converted
        integer, intent(in), optional :: decimal_places      ! Number of decimal places (optional, default is 2)
        character(len=64) :: str                             ! Resulting string (maximum length of 64 characters)
        integer :: dp                                        ! Final decimal places for the format string

        ! Set default decimal places to 2 if not provided
        dp = 2
        if (present(decimal_places)) dp = decimal_places

        ! Ensure that the format string uses the specified number of decimal places
        write (str, "(F0.2)") value                          ! Default format with 2 decimal places
        if (dp > 0) then
            if (abs(value) .ge. 1.0d0) then
                write (str, "(F0." // trim(itoa(dp)) // ")") value
            else
                write (str, "(A, F0." // trim(itoa(dp)) // ")") "0", value
            end if
        end if
    end function rtoa

    ! ---------------------------------------------------------------------------------------------------
    ! @brief Replace space characters from a string by replacing spaces with underscores.
    ! @param[in] input_string The input string to modify.
    ! @return The modified string with spaces replaced by underscores.
    ! ---------------------------------------------------------------------------------------------------
    function replace_space(input_string) result(output_string)
        implicit none
        character(len=*), intent(in) :: input_string                    ! Input string
        character(len=len(input_string)) :: output_string               ! Output string with spaces replaced

        integer :: i                                                    ! Loop index

        output_string = input_string                                    ! Initialise the output string

        do i = 1, len(output_string)
            if (output_string(i:i) .eq. " ") then
                output_string(i:i) = "_"                                ! Replace spaces with underscores
            end if
        end do
    end function replace_space
end module vtklib