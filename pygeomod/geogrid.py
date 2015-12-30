'''Module with classes and methods to analyse and process exported geomodel grids

Created on 21/03/2014

@author: Florian Wellmann (some parts originally developed by Erik Schaeffer)
'''
import numpy as np
import pynoddy
import subprocess
import os.path
import platform

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("\n\n\tMatplotlib not installed - plotting functions will not work!\n\n\n")

# import mpl_toolkits
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# to convert python variable types to cpp types
import ctypes
# to create array
from numpy.ctypeslib import ndpointer
# to create folder
import os
# read out and change xml file (here only used to read out model boundary information)
import geomodeller_xml_obj as GO


class GeoGrid():
    """Object definition for exported geomodel grids"""

    def __init__(self, **kwds):
        """GeoGrid contains methods to load, analyse, and process exported geomodel grids

        **Optional Keywords**:
            - *grid_filename* = string : filename of exported grid
            - *delxyz_filename* = string : file with model discretisation
            - *dimensions_filename* = string : file with model dimension (coordinates)
        """

        if kwds.has_key('grid_filename'):
            self.grid_filename = kwds['grid_filename']
        if kwds.has_key('delxyz_filename'):
            self.delxyz_filename = kwds['delxyz_filename']
        if kwds.has_key('dimensions_filename'):
            self.dimensions_filename = kwds['dimensions_filename']

    def __add__(self, G_other):
        """Combine grid with another GeoGrid if regions are overlapping"""
        # check overlap
        print self.ymin, self.ymax
        print G_other.ymin, G_other.ymax
        if (G_other.ymin < self.ymax and G_other.ymin > self.ymin):
            print("Grids overlapping in y-direction between %.0f and %.0f" %
                  (G_other.ymin, self.ymax))

    def load_grid(self):
        """Load exported grid, discretisation and dimensions from file"""
        if not hasattr(self, 'grid_filename'):
            raise AttributeError("Grid filename is not defined!")
        self.grid = np.loadtxt(self.grid_filename,
                               delimiter = ',',
                               dtype='int',
                               unpack=False)
        if hasattr(self, 'delxyz_filename'):
            self.load_delxyz(self.delxyz_filename)
            self.adjust_gridshape()
        if hasattr(self, 'dimensions_filename'):
            self.load_dimensions(self.dimensions_filename)

    def load_delxyz(self, delxyz_filename):
        """Load grid discretisation from file"""
        del_lines = open(delxyz_filename, 'r').readlines()
        d0 = del_lines[0].split("*")
        self.delx = np.array([float(d0[1]) for _ in range(int(d0[0]))])
        d1 = del_lines[1].split("*")
        self.dely = np.array([float(d1[1]) for _ in range(int(d1[0]))])
        d2 = del_lines[2].split(",")[:-1]
        self.delz = np.array([float(d) for d in d2])
        (self.nx, self.ny, self.nz) = (len(self.delx), len(self.dely), len(self.delz))
        (self.extent_x, self.extent_y, self.extent_z) = (sum(self.delx), sum(self.dely), sum(self.delz))

    def set_delxyz(self, delxyz):
        """Set delx, dely, delz arrays explicitly and update additional attributes

        **Arguments**:
            - *delxyz* = (delx-array, dely-array, delz-array): arrays with cell dimensions
        """
        self.delx, self.dely, self.delz = delxyz
        (self.nx, self.ny, self.nz) = (len(self.delx), len(self.dely), len(self.delz))
        (self.extent_x, self.extent_y, self.extent_z) = (sum(self.delx), sum(self.dely), sum(self.delz))

    def set_basename(self, name):
        """Set basename for grid exports, etc.

        **Arguments**:
            - *name* = string: basename
        """
        self.basename = name

    def load_dimensions(self, dimensions_filename):
        """Load project dimensions from file"""
        dim = [float(d) for d in open(dimensions_filename, 'r').readlines()[1].split(",")]
        (self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax) = dim
        # calculate cell centre positions in real world coordinates

    def define_regular_grid(self, nx, ny, nz):
        """Define a regular grid from defined project boundaries and given discretisations"""
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.delx = np.ones(nx) * (self.xmax - self.xmin) / nx
        self.dely = np.ones(ny) * (self.ymax - self.ymin) / ny
        self.delz = np.ones(nz) * (self.zmax - self.zmin) / nz
        # create (empty) grid object
        self.grid = np.ndarray((nx, ny, nz))
        # update model extent
        (self.extent_x, self.extent_y, self.extent_z) = (sum(self.delx), sum(self.dely), sum(self.delz))


    def define_irregular_grid(self, delx, dely, delz):
        """Set irregular grid according to delimter arrays in each direction"""
        self.delx = delx
        self.dely = dely
        self.delz = delz
        self.nx = len(delx)
        self.ny = len(dely)
        self.nz = len(delz)
        # create (empty) grid object
        self.grid = np.ndarray((self.nx, self.ny, self.nz))
        # update model extent
        (self.extent_x, self.extent_y, self.extent_z) = (sum(self.delx), sum(self.dely), sum(self.delz))

    def get_dimensions_from_geomodeller_xml_project(self, xml_filename):
        """Get grid dimensions from Geomodeller project

        **Arguments**:
            - *xml_filename* = string: filename of Geomodeller XML file
        """
        # Note: this implementation is based on the Geomodeller API
        # The boundaries could theoretically also be extracted from the XML file
        # directly, e.g. using the geomodeller_xml_obj module - but this would
        # require an additional module being loaded, so avoid here!
        filename_ctypes = ctypes.c_char_p(xml_filename)
        # get model boundaries
            #Detection of operative system:
        if platform.system() == "Linux":
            lib = ctypes.CDLL('./libgeomod.so') #linux
        elif platform.system() == "Windows":
            lib = ctypes.windll.LoadLibrary(os.path.dirname(os.path.abspath(__file__)) + os.path.sep +"libgeomodwin.dll")     #windows
        else:
            print("Your operative system is not supported")
        lib.get_model_bounds.restype = ndpointer(dtype=ctypes.c_int, shape=(6,))
        boundaries = lib.get_model_bounds(filename_ctypes)
        (self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax) = boundaries
        self.extent_x = self.xmax - self.xmin
        self.extent_y = self.ymax - self.ymin
        self.extent_z = self.zmax - self.zmin


# functions to try to improve the performance. They sucks yet.
    def update_from_geomodeller_project_OP(self, xml_filename):
        """Update grid properties directly from Geomodeller project

        **Arguments**:
            - *xml_filename* = string: filename of Geomodeller XML file
        """
        filename_ctypes = ctypes.c_char_p(xml_filename)



        #print filename_ctypes
        # create cell position list with [x0, y0, z0, ... xn, yn, zn]
        cell_position = []
        ids = []
        formations_raw = []
        # check if cell centers are defined - if not, do so!

        if platform.system() == "Linux":
            lib = ctypes.CDLL('./libgeomod.so') #linux
        elif platform.system() == "Windows":
            lib = ctypes.windll.LoadLibrary(os.path.dirname(os.path.abspath(__file__)) + os.path.sep +"libgeomodwin.dll") #windows
        else:
            print("Your operative system is not supported")

        lib.compute_irregular_grid.restype = ndpointer(dtype=ctypes.c_int, shape=(3,))
        # This is the function wich needs GPU!!!

        #lib.init_api(filename_ctypes)

        # check if cell centers are defined - if


        if not hasattr(self, 'cell_centers_x'):
            self.determine_cell_centers()
            for k in range(self.nz):
                for j in range(self.ny):
                    for i in range(self.nx):
                        cell_position.append(self.cell_centers_x[i])
                        cell_position.append(self.cell_centers_y[j])
                        cell_position.append(self.cell_centers_z[k])
                        print len(cell_position), cell_position
                        coord_ctypes = (ctypes.c_double * 3)(*cell_position)
                        formation2 = lib.compute_irregular_grid_op(coord_ctypes)
                        #formations_raw.append(lib.compute_irregular_grid(coord_ctypes))
                        print "e", formation2
                        cell_position = []
                        ids.append((i,j,k))
                        print len(formations_raw)
        # prepare variables for cpp function

        # call cpp function
                    #Detection of operative system:
        # re-sort formations into array

        for i in range(len(formations_raw)):
            self.grid[ids[i][0],ids[i][1],ids[i][2]] = formations_raw[i]


    def update_from_geomodeller_project_OP2(self,xml_filename):
        """Update grid properties directly from Geomodeller project

        **Arguments**:
            - *xml_filename* = string: filename of Geomodeller XML file
        """
        filename_ctypes = ctypes.c_char_p(xml_filename)

        if platform.system() == "Linux":
            lib = ctypes.CDLL('./libgeomod.so') #linux
        elif platform.system() == "Windows":
            lib = ctypes.windll.LoadLibrary(os.path.dirname(os.path.abspath(__file__)) + os.path.sep +"libgeomodwin_test2.dll") #windows
        else:
            print("Your operative system is not supported")

        #print filename_ctypes
        # create cell position list with [x0, y0, z0, ... xn, yn, zn]
        cell_position = []
        ids = []
        formations_sol = []
        # check if cell centers are defined - if not, do so!
        if not hasattr(self, 'cell_centers_x'):
            self.determine_cell_centers()
            for k in range(self.nz):
                for j in range(self.ny):
                    for i in range(self.nx):
                        cell_position.append(self.cell_centers_x[i])
                        cell_position.append(self.cell_centers_y[j])
                        cell_position.append(self.cell_centers_z[k])


                        coord_ctypes = (ctypes.c_double * len(cell_position))(*cell_position)
                        coord_len = len(cell_position)

                        lib.compute_irregular_grid.restype = ndpointer(dtype=ctypes.c_int, shape=(coord_len/3,))

                        formations_raw = lib.compute_irregular_grid(filename_ctypes, coord_ctypes, coord_len)
                        formations_sol.append(formations_raw)
                    #    print formations_raw
                        cell_position = []
                        ids.append((i,j,k))

        # prepare variables for cpp function

        # call cpp function
                    #Detection of operative system:

        for i in range(len(formations_sol)):
            self.grid[ids[i][0],ids[i][1],ids[i][2]] = formations_sol[i]



    def update_from_geomodeller_project(self, xml_filename):
        """Update grid properties directly from Geomodeller project

        **Arguments**:
            - *xml_filename* = string: filename of Geomodeller XML file
        """
        filename_ctypes = ctypes.c_char_p(xml_filename)
        #print filename_ctypes
        # create cell position list with [x0, y0, z0, ... xn, yn, zn]
        cell_position = []
        ids = []
        # check if cell centers are defined - if not, do so!
        if not hasattr(self, 'cell_centers_x'):
            if cell_position == []:
                self.determine_cell_centers()
                for k in range(self.nz):
                    for j in range(self.ny):
                        for i in range(self.nx):
                            cell_position.append(self.cell_centers_x[i])
                            cell_position.append(self.cell_centers_y[j])
                            cell_position.append(self.cell_centers_z[k])
                            ids.append((i,j,k))

        # prepare variables for cpp function
        #print cell_position
        coord_ctypes = (ctypes.c_double * len(cell_position))(*cell_position)

        coord_len = len(cell_position)
        # call cpp function
                    #Detection of operative system:
        if platform.system() == "Linux":
            lib = ctypes.CDLL('./libgeomod.so') #linux
        elif platform.system() == "Windows":
            lib = ctypes.windll.LoadLibrary(os.path.dirname(os.path.abspath(__file__)) + os.path.sep +"libgeomodWin.dll") #windows
        else:
            print("Your operative system is not supported")
        #print coord_len
        lib.compute_irregular_grid.restype = ndpointer(dtype=ctypes.c_int, shape=(coord_len/3,))
        # This is the function wich needs GPU!!!
        formations_raw = lib.compute_irregular_grid(filename_ctypes, coord_ctypes, coord_len)
        # re-sort formations into array
        #    print formations_raw
        for i in range(len(formations_raw)):
            self.grid[ids[i][0],ids[i][1],ids[i][2]] = formations_raw[i]

    def set_densities(self, densities):
        """Set layer densities

        **Arguments**:
            - *densities* = dictionary of floats: densities for geology ids
        """
        self.densities = densities

    def set_sus(self, sus):
        """Set layer susceptibilities

        **Arguments**:
            - *us* = dictionary of floats: magnetic susceptibilities for geology ids
        """
        self.sus = sus

    def write_noddy_files(self, **kwds):
        """Create Noddy block model files (for grav/mag calculation)

        **Optional keywords**:
            - *gps_range* = float : set GPS range (default: 1200.)

        Method generates the files required to run the forward gravity/ magnetics response
        from the block model:
        - model.g00 = file with basic model information
        - model.g12 = discretised geological (block) model
        - base.his = Noddy history file with some basic settings
        """
        self.gps_range = kwds.get("gps_range", 1200.)

        if not hasattr(self, 'basename'):
            self.basename = "geogrid"
        f_g12 = open(self.basename + ".g12", 'w')
        f_g01 = open(self.basename + ".g00", 'w')
        # method = 'numpy'    # using numpy should be faster - but it messes up the order... possible to fix?
        #         if method == 'standard':
        #             i = 0
        #             j = 0
        #             k = 0
        #             self.block = np.ndarray((self.nx,self.ny,self.nz))
        #             for line in f.readlines():
        #                 if line == '\n':
        #                     # next z-slice
        #                     k += 1
        #                     # reset x counter
        #                     i = 0
        #                     continue
        #                 l = [int(l1) for l1 in line.strip().split("\t")]
        #                 self.block[i,:,self.nz-k-1] = np.array(l)[::-1]
        #                 i += 1

        if not hasattr(self, "unit_ids"):
            self.determine_geology_ids()

        #=======================================================================
        # # create file with base settings (.g00)
        #=======================================================================
        f_g01.write("VERSION = 7.11\n")
        f_g01.write("FILE PREFIX = " + self.basename + "\r\n")
        import time
        t = time.localtime() # get current time
        f_g01.write("DATE = %d/%d/%d\n" % (t.tm_mday, t.tm_mon, t.tm_year))
        f_g01.write("TIME = %d:%d:%d\n" % (t.tm_hour, t.tm_min, t.tm_sec))
        f_g01.write("UPPER SW CORNER (X Y Z) = %.1f %.1f %.1f\n" % (self.xmin - self.gps_range,
                                                                    self.ymin - self.gps_range,
                                                                    self.zmax))
        f_g01.write("LOWER NE CORNER (X Y Z) = %.1f %.1f %.1f\n" % (self.xmax + self.gps_range,
                                                                    self.ymax + self.gps_range,
                                                                    self.zmin))
        f_g01.write("NUMBER OF LAYERS = %d\n" % self.nz)
        for k in range(self.nz):
            f_g01.write("\tLAYER %d DIMENSIONS (X Y) = %d %d\n" % (k,
                                                                   self.nx + 2 * (self.gps_range / self.delx[0]),
                                                                   self.ny + 2 * (self.gps_range / self.dely[0])))
        f_g01.write("NUMBER OF CUBE SIZES = %d\n" % self.nz)
        for k in range(self.nz):
            f_g01.write("\tCUBE SIZE FOR LAYER %d = %d\n" % (k, self.delx[0]))
        f_g01.write("CALCULATION RANGE = %d\n" % (self.gps_range / self.delx[0]))
        f_g01.write("""INCLINATION OF EARTH MAG FIELD = -67.00
        INTENSITY OF EARTH MAG FIELD = 63000.00
        DECLINATION OF VOL. WRT. MAG NORTH = 0.00
        DENSITY CALCULATED = Yes
        SUSCEPTIBILITY CALCULATED = Yes
        REMANENCE CALCULATED = No
        ANISOTROPY CALCULATED = No
        INDEXED DATA FORMAT = Yes
        """)
        f_g01.write("NUM ROCK TYPES = %d\n" % len(self.unit_ids))
        for i in self.unit_ids:
            f_g01.write("ROCK DEFINITION Layer %d = %d\n" % (i, i))
            f_g01.write("\tDensity = %f\n" % self.densities[int(i)])
            f_g01.write("\tSus = %f\n" % self.sus[int(i)])

        #=======================================================================
        # Create g12 file
        #=======================================================================

        # write geology blocks to file
        for k in range(self.nz):
            # this worked for geophysics, but not for re-import with pynoddy:
            #             for val in self.grid[:,:,k].ravel(order = 'A'):
            #                 f_g12.write("%d\t" % val)
            for i in range(self.nx):
                for val in self.grid[i,:,k]:
                    f_g12.write("%d\t" % val)
                f_g12.write("\r\n")
            # f_g12.write(['%d\t' % i for i in self.grid[:,:,k].ravel()])
            f_g12.write("\r\n")

        f_g12.close()
        f_g01.close()

        #=======================================================================
        # # create noddy history file for base settings
        #=======================================================================
        import pynoddy.history
        history = self.basename + "_base.his"
        nm = pynoddy.history.NoddyHistory('simple_two_faults.his')
        #print nm
        # add stratigraphy
        # create dummy names and layers for base stratigraphy
        layer_names = []
        layer_thicknesses = []
        for i in self.unit_ids:
            layer_names.append('Layer %d' % i)
            layer_thicknesses.append(500)
        strati_options = {'num_layers' : len(self.unit_ids),
                          'layer_names' : layer_names,
                          'layer_thickness' : layer_thicknesses}
        nm.add_event('stratigraphy', strati_options, )

        # set grid origin and extent:
        nm.set_origin(self.xmin, self.ymin, self.zmin)
        nm.set_extent(self.extent_x, self.extent_y, self.extent_z)

        nm.write_history(history)

    def analyse_geophysics(self, densities, **kwds):
        """Simulate potential-fields and use for model analysis

        It is possible to directly define filter for processing of gravity

        **Arguments**:
            - *model_dir*: directory containing sub-directories with uncertainty runs
                        (uncertainty_run_01, etc.);

        **Optional keywords**:
            - *grav_min* = float : reject models with a grav value lower than this
            - *grav_max* = float : reject models with a grav value larger than this
        """
        #os.chdir(model_dir)
        all_gravs = []
        all_gravs_filtered = []
        all_mags = []
        all_mags_filtered = []
        all_probs = {}
        all_probs_filtered = {}
        i_all = 0
        i_filtered = 0
        used_grids = []
        used_grids_filtered = []
        f = self
        #for f in os.listdir(model_dir):
        #    if os.path.splitext(f)[1] == ".pkl" and "Sandstone" in f:
        #===================================================================
        # Load grid
        #===================================================================
        #    grid_file = open(os.path.join(model_dir, f), "r")
        #    grid_ori = pickle.load(grid_file)
        #    grid_file.close()
        #===================================================================
        # Extract subgrid
        #===================================================================
        #    subrange = (40,200,30,250,0,80)
        grid = self
            # grid = grid_ori
        # substitute 0 with something else in grid ids
        tmp_grid = np.zeros_like(grid.grid)
        tmp_grid[grid.grid == 0] += 1
        #print tmp_grid.shape
    #    grid.set_basename(self.split(".")[0])
    #    print "Basename"
    #    print grid.basename
        grid.grid += tmp_grid
        #             n_cells = np.prod(grid.grid.shape)
        grid.determine_geology_ids()
        #===================================================================
        # # set densities and magnetic susceptibilities
        #===================================================================
        #densities = dens
            #densities = {0: 0.1,
            #       1: 2610,
            #       2: 2920,
            #       3: 3100,
            #       4: 2920,
            #       5: 2610}
        sus = {0: 0.001,
               1: 0.001,
               2: 0.001,
               3: 0.1,
               4: 0.001,
               5: 0.001}
        grid.set_densities(densities)
        grid.set_sus(sus)
    #    print grid
        grid.write_noddy_files(gps_range = 0.0)
    #    print grid.unit_ids
        sim_type = "ANOM_FROM_BLOCK"
        history = grid.basename + "_base.his"
        output_name = grid.basename
        #import pdb
        #pdb.set_trace()

        # save grid as vtk for testing:
        # grid_ori.export_to_vtk(vtk_filename = grid.basename)
        #===================================================================
        # Run gravity forward modeling
        #===================================================================
        out =  subprocess.Popen(['noddy.exe', history, output_name, sim_type],
                    shell=True, stderr=subprocess.PIPE,
                    stdout=subprocess.PIPE).stdout.read()

        #===================================================================
        # Initiate probability grids
        #===================================================================
        if i_all == 0:
            vals = grid.unit_ids
            for val in vals:
                if not all_probs.has_key(val):
                    all_probs[val] = np.zeros_like(grid.grid, dtype = "float")
                    all_probs_filtered[val] = np.zeros_like(grid.grid, dtype = "float")

        #===================================================================
        # Create plot and store data
        #===================================================================
        self.geophys = pynoddy.output.NoddyGeophysics(grid.basename)
        #return geophys
        """
            #=================================================
            #             Check gravity constraints
            #=================================================
        filter_val = True
        if kwds.has_key("grav_max"):
            if np.max(geophys.grv_data) > kwds['grav_max']:
                filter_val = False
        if kwds.has_key("grav_min"):
            if np.min(geophys.grv_data) < kwds['grav_min']:
                filter_val = False
        if filter_val:
            all_gravs_filtered.append(geophys.grv_data)
            all_mags_filtered.append(geophys.mag_data)
            used_grids_filtered.append("%s/%s" % (model_dir, grid.basename))
            #                 test_grid = np.zeros_like(grid.grid)
            for val in vals:
                all_probs_filtered[val] += (grid.grid == val)
                #                     test_grid += grid.grid == val
            # check probability
            #                 assert(np.sum(test_grid) == n_cells)
            i_filtered += 1
        all_gravs.append(geophys.grv_data)
        all_mags.append(geophys.mag_data)
        used_grids.append("%s/%s" % (model_dir, grid.basename))
        #             test_grid = np.zeros_like(grid.grid)
        for val in vals:
            all_probs[val] += (grid.grid == val)
            #                 test_grid += grid.grid == val
            #             assert(np.sum(test_grid) == n_cells)
        i_all += 1

        #===================================================================
        # Export to vtk for test
        #===================================================================
        #             grid_out = pynoddy.output.NoddyOutput(grid.basename)
        #             grid_out.export_to_vtk(vtk_filename = grid.basename)


        #=======================================================================
        # Analyse statistics for all simulated grids
        #=======================================================================
        # all_gravs = np.array(all_gravs)
        return all_gravs, all_mags, used_grids, all_probs, i_all,\
               all_gravs_filtered, all_mags_filtered, used_grids_filtered, all_probs_filtered, i_filtered
        #     f_all = open("all_gravs.pkl", 'w')
        #     pickle.dump(all_gravs, f_all)
        #     f_all.close()
        #     return all_gravs
        """


    def set_dimensions(self, **kwds):
        """Set model dimensions, if no argument provided: xmin = 0, max = sum(delx) and accordingly for y,z

        **Optional keywords**:
            - *dim* = (xmin, xmax, ymin, ymax, zmin, zmax) : set dimensions explicitly
        """
        if kwds.has_key("dim"):
            (self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax) = kwds['dim']
        else:
            self.xmin, self.ymin, self.zmin = (0., 0., 0.)
            self.xmax, self.ymax, self.zmax = (sum(self.delx), sum(self.dely), sum(self.delz))

    def determine_cell_centers(self):
        """Determine cell centers for all coordinate directions in "real-world" coordinates"""
        if not hasattr(self, 'xmin'):
            raise AttributeError("Please define grid dimensions first")
        sum_delx = np.cumsum(self.delx)
        sum_dely = np.cumsum(self.dely)
        sum_delz = np.cumsum(self.delz)
        self.cell_centers_x = np.array([sum_delx[i] - self.delx[i] / 2. for i in range(self.nx)]) + self.xmin
        self.cell_centers_y = np.array([sum_dely[i] - self.dely[i] / 2. for i in range(self.ny)]) + self.ymin
        self.cell_centers_z = np.array([sum_delz[i] - self.delz[i] / 2. for i in range(self.nz)]) + self.zmin

    def determine_cell_boundaries(self):
        """Determine cell boundaries for all coordinates in "real-world" coordinates"""
        if not hasattr(self, 'xmin'):
            raise AttributeError("Please define grid dimensions first")
        sum_delx = np.cumsum(self.delx)
        sum_dely = np.cumsum(self.dely)
        sum_delz = np.cumsum(self.delz)
        self.boundaries_x = np.ndarray((self.nx+1))
        self.boundaries_x[0] = 0
        self.boundaries_x[1:] = sum_delx
        self.boundaries_y = np.ndarray((self.ny+1))
        self.boundaries_y[0] = 0
        self.boundaries_y[1:] = sum_dely
        self.boundaries_z = np.ndarray((self.nz+1))
        self.boundaries_z[0] = 0
        self.boundaries_z[1:] = sum_delz

        # create a list with all bounds
        self.bounds = [self.boundaries_y[0], self.boundaries_y[-1],
                       self.boundaries_x[0], self.boundaries_x[-1],
                       self.boundaries_z[0], self.boundaries_z[-1]]


    def adjust_gridshape(self):
        """Reshape numpy array to reflect model dimensions"""
        self.grid = np.reshape(self.grid, (self.nz, self.ny, self.nx))
        self.grid = np.swapaxes(self.grid, 0, 2)
        # self.grid = np.swapaxes(self.grid, 0, 1)

    def plot_section(self, direction, cell_pos='center', **kwds):
        """Plot a section through the model in a given coordinate direction

        **Arguments**:
            - *direction* = 'x', 'y', 'z' : coordinate direction for section position
            - *cell_pos* = int/'center','min','max' : cell position, can be given as
            value of cell id, or as 'center' (default), 'min', 'max' for simplicity

        **Optional Keywords**:
            - *cmap* = mpl.colormap : define colormap for plot (default: jet)
            - *alpha* = define the transparencie for the plot (default: 1)
            - *colorbar* = bool: attach colorbar (default: True)
            - *geomod_coord* = bool:  Plotting geomodeller coordinates instead voxels (default: False)
            - *contour* = bool : Plotting contour of the layers contact
            - *plot layer* = array: Layer Number you want to plot in the contour plot
            - *rescale* = bool: rescale color bar to range of visible slice (default: False)
            - *ve* = float : vertical exageration (for plots in x,y-direction)
            - *figsize* = (x,y) : figsize settings for plot
            - *ax* = matplotlib.axis : add plot to this axis (default: new axis)
                    if axis is defined, the axis is returned and the plot not shown
                    Note: if ax is passed, colorbar is False per default!
            - *savefig* = bool : save figure to file (default: show)
            - *fig_filename* = string : filename to save figure
        """
        # TO DO:
        # - Colorbar in contourplots
        colorbar = kwds.get('colorbar', True)
        cmap = kwds.get('cmap', 'jet')
        alpha = kwds.get('alpha', 1)
        rescale = kwds.get('rescale', False)
        ve = kwds.get('ve', 1.)
        figsize = kwds.get('figsize', (8,4))
        geomod_coord =  kwds.get('geomod_coord', False)
        contour = kwds.get('contour', False)
        linewidth = kwds.get("linewidth", 1)
        levels = kwds.get("plot_layer", None)

        if not kwds.has_key('ax'):
            colorbar = kwds.get('colorbar', True)
            # create new axis for plot
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        else:
            colorbar = False
            ax = kwds['ax']

        if direction == 'x':
            if type(cell_pos) == str:
                # decipher cell position
                if cell_pos == 'center' or cell_pos == 'centre':
                    pos = self.nx / 2
                elif cell_pos == 'min':
                    pos = 0
                elif cell_pos == 'max':
                    pos = self.nx
            else:
                pos = cell_pos
            grid_slice = self.grid[pos,:,:]
            grid_slice = grid_slice.transpose()
            #print grid_slice
            aspect = self.extent_z/self.extent_x * ve
            if geomod_coord:

                ax.set_xticks(np.linspace(0,self.ny-1,6, endpoint = True, dtype = int))
                ax.set_yticks(np.linspace(0,self.nz-1,6, endpoint = True, dtype = int))
                ax.set_xticklabels(np.linspace(self.ymin,self.ymax,6,dtype = int, endpoint = True ))
                ax.set_yticklabels(np.linspace(self.zmin,self.zmax,6,dtype = int, endpoint = True ))


                ax.set_ylabel("z[m]")
                ax.set_xlabel("y[m]")
            else:
                ax.set_ylabel("z[voxels]")
                ax.set_xlabel("y[voxels]")

            if contour:
                ry = np.arange(self.nz)
                rx = np.arange(self.ny)

        elif direction == 'y':
            if type(cell_pos) == str:
                # decipher cell position
                if cell_pos == 'center' or cell_pos == 'centre':
                    pos = self.ny / 2
                elif cell_pos == 'min':
                    pos = 0
                elif cell_pos == 'max':
                    pos = self.ny
            else:
                pos = cell_pos
            grid_slice = self.grid[:,pos,:]
            grid_slice = grid_slice.transpose()
            aspect = self.extent_z/self.extent_y * ve
            if geomod_coord:
                #print np.linspace(0,self.extent_x,11), np.linspace(0,self.extent_x,11, endpoint = True)



                ax.set_xticks(np.linspace(0,self.nx-1,6, endpoint = True, dtype = int))
                ax.set_yticks(np.linspace(0,self.nz-1,6, endpoint = True, dtype = int))

                ax.set_xticklabels(np.linspace(self.xmin,self.xmax,6, endpoint = True, dtype = int))
                ax.set_yticklabels(np.linspace(self.zmin,self.zmax,6, dtype = int,endpoint = True ))
                #ax.invert_yaxis
                ax.set_ylabel("z[m]")
                ax.set_xlabel("x[m]")
            else:
                ax.set_ylabel("z[voxels]")
                ax.set_xlabel("x[voxels]")

            if contour:
                ry = np.arange(self.nz)
                rx = np.arange(self.nx)

        elif direction == 'z' :

            if type(cell_pos) == str:
                # decipher cell position
                if cell_pos == 'center' or cell_pos == 'centre':
                    pos = self.nz / 2
                elif cell_pos == 'min':
                    pos = 0
                elif cell_pos == 'max':
                    pos = self.nz
            else:
                pos = cell_pos
            grid_slice = self.grid[:,:,pos].transpose()
            aspect = 1.
            # setting labels
            if geomod_coord:

            #    print self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax
            #    print np.linspace(self.xmin,self.xmax,6, endpoint = False, dtype = int)
            #    print np.linspace(0,self.extent_y,6, endpoint = False, dtype = int)
                ax.set_xticks(np.linspace(0,self.nx-1,6,dtype = int, endpoint = True ))
                ax.set_yticks(np.linspace(0,self.ny-1,6, endpoint = True, dtype = int))
                ax.set_xticklabels(np.linspace(self.xmin,self.xmax,6,dtype = int, endpoint = True ))
                ax.set_yticklabels(np.linspace(self.ymin,self.ymax,6,dtype = int, endpoint = True ))
                ax.set_ylabel("y[m]")
                ax.set_xlabel("x[m]")
            else:
                ax.set_ylabel("y[voxels]")
                ax.set_xlabel("x[voxels]")

            if contour:
                ry = np.arange(self.ny)
                rx = np.arange(self.nx)

        if not hasattr(self, 'unit_ids'):
            self.determine_geology_ids()
        if rescale:
            vmin = np.min(grid_slice)
            vmax = np.max(grid_slice)
        else: # use global range for better comparison
            vmin = min(self.unit_ids)
            vmax = max(self.unit_ids)


        if contour:
            Rx, Ry = np.meshgrid(rx, ry)
            #print np.amax(grid_slice)
            im = ax.contour(Rx, Ry, grid_slice, int(np.amax(grid_slice)+1),
                interpolation='nearest', cmap = cmap, alpha = alpha, linewidths = linewidth,
                 antialiased = True, levels = levels)


        else:
            im = ax.imshow(grid_slice, interpolation='nearest',
                      cmap = cmap,
                       origin='lower_left',
                       vmin = vmin,
                      vmax = vmax,
                      aspect = aspect)
        if colorbar:
       #            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
       #            cax = divider.append_axes("bottom", size="5%", pad=0.2)
            cbar1 = fig.colorbar(im, orientation="horizontal")
            ticks = np.arange(vmin, vmax+0.1, int(np.log2(vmax-vmin)/1.2), dtype='int')
            cbar1.set_ticks(ticks)
         #         cbar1.set_ticks(self.unit_ids[::int(np.log2(len(self.unit_ids)/2))])
            cbar1.set_label("Geology ID")
        #         cax.xaxis.set_major_formatter(FormatStrFormatter("%d"))

        if kwds.has_key("ax"):
            # return image and do not show
            return im

        if kwds.has_key('savefig') and kwds['savefig']:
            # save to file
            filename = kwds.get("fig_filename", "grid_section_direction_%s_pos_%d.png" %
                                (direction, cell_pos))
            plt.savefig(filename, transparent = True)

        else:

            plt.show()






    def export_to_vtk(self, vtk_filename="geo_grid", real_coords = True, **kwds):
        """Export grid to VTK for visualisation

        **Arguments**:
            - *vtk_filename* = string : vtk filename (obviously...)
            - *real_coords* = bool : model extent in "real world" coordinates

        **Optional Keywords**:
            - *grid* = numpy grid : grid to save to vtk (default: self.grid)
            - *var_name* = string : name of variable to plot (default: Geology)

        Note: requires pyevtk, available at: https://bitbucket.org/pauloh/pyevtk
        """
        grid = kwds.get("grid", self.grid)
        var_name = kwds.get("var_name", "Geology")
        #from evtk.hl import gridToVTK
        import pyevtk
        from pyevtk.hl import gridToVTK
        # define coordinates
        x = np.zeros(self.nx + 1)
        y = np.zeros(self.ny + 1)
        z = np.zeros(self.nz + 1)
        x[1:] = np.cumsum(self.delx)
        y[1:] = np.cumsum(self.dely)
        z[1:] = np.cumsum(self.delz)



        # plot in coordinates
        if real_coords:
            x += self.xmin
            y += self.ymin
            z += self.zmin


        gridToVTK(vtk_filename, x, y, z,
                  cellData = {var_name: grid})

    def export_to_csv(self, filename = "geo_grid.csv"):
        """Export grid to x,y,z,value pairs in a csv file

        Ordering is x-dominant (first increase in x, then y, then z)

        **Arguments**:
            - *filename* = string : filename of csv file (default: geo_grid.csv)
        """
        f = open(filename, 'w')
        for zz in self.delz:
            for yy in self.dely:
                for xx in self.delx:
                    f.write("%.1f,%.1f,%.1f,%.d" % (xx,yy,zz,self.grid[xx,yy,zz]))
        f.close()


    def determine_geology_ids(self):
        """Determine all ids assigned to cells in the grid"""
        self.unit_ids = np.unique(self.grid)

    def get_name_mapping_from_file(self, filename):
        """Get the mapping between unit_ids in the model and real geological names
        from a csv file (e.g. the SHEMAT property file)

        **Arguments**:
            - *filename* = string : filename of csv file with id, name entries
        """
        self.unit_name = {}
        filelines = open(filename, 'r').readlines()[1:]
        for line in filelines:
            l = line.split(",")
            self.unit_name[int(l[1])] = l[0]

    def get_name_mapping_from_dict(self, unit_name_dict):
        """Get the name mapping directly from a dictionary

        **Arguments**:
            - *unit_name_dict* = dict with "name" : unit_id (int) pairs
        """
        self.unit_name = unit_name_dict


    def remap_ids(self, mapping_dictionary):
        """Remap geological unit ids to new ids as defined in mapping dictionary

        **Arguments**:
            - *mapping_dictionary* = dict : {1 : 1, 2 : 3, ...} : e.g.: retain
            id 1, but map id 2 to 3 (note: if id not specified, it will be retained)
        """
        # first step: create a single mesh for each id to avoid accidential
        # overwriting below (there might be a better solution...)
        if not hasattr(self, 'unit_ids'):
            self.determine_geology_ids()
        geol_grid_ind = {}
        for k,v in mapping_dictionary.items():
            geol_grid_ind[k] = self.grid == k
            print("Remap id %d -> %d" % (k,v))
        # now reassign values in actual grid
        for k,v in mapping_dictionary.items():
            print("Reassign id %d to grid" % v)
            self.grid[geol_grid_ind[k]] = v
        # update global geology ids
        self.determine_geology_ids()

    def determine_cell_volumes(self):
        """Determine cell volumes for each cell (e.g. for total formation volume calculation)"""
        self.cell_volume = np.ndarray(np.shape(self.grid))
        for k,dz in enumerate(self.delz):
            for j,dy in enumerate(self.dely):
                for i,dx in enumerate(self.delx):
                    self.cell_volume[i,j,k] = dx * dy * dz


    def determine_indicator_grids(self):
        """Determine indicator grids for all geological units"""
        self.indicator_grids = {}
        if not hasattr(self, 'unit_ids'):
            self.determine_geology_ids()
        grid_ones = np.ones(np.shape(self.grid))
        for unit_id in self.unit_ids:
            self.indicator_grids[unit_id] = grid_ones * (self.grid == unit_id)

    def determine_id_volumes(self):
        """Determine the total volume of each unit id in the grid

        (for example for cell discretisation studies, etc."""
        if not hasattr(self, 'cell_volume'):
            self.determine_cell_volumes()
        if not hasattr(self, 'indicator_grids'):
            self.determine_indicator_grids()
        self.id_volumes = {}
        for unit_id in self.unit_ids:
            self.id_volumes[unit_id] = np.sum(self.indicator_grids[unit_id] * self.cell_volume)

    def print_unit_names_volumes(self):
        """Formatted output to STDOUT of unit names (or ids, if names are note
        defined) and calculated volumes
        """
        if not hasattr(self, 'id_vikumes'):
            self.determine_id_volumes()

        if hasattr(self, "unit_name"):
            # print with real geological names
            print("Total volumes of modelled geological units:\n")
            for unit_id in self.unit_ids:
                print("%26s : %.2f km^3" % (self.unit_name[unit_id],
                                            self.id_volumes[unit_id]/1E9))
        else:
            # print with unit ids only
            print("Total volumes of modelled geological units:\n")
            for unit_id in self.unit_ids:
                print("%3d : %.2f km^3" % (unit_id,
                                            self.id_volumes[unit_id]/1E9))


    def extract_subgrid(self, subrange, **kwds):
        """Extract a subgrid model from existing grid

        **Arguments**:
            - *subrange* = (x_from, x_to, y_from, y_to, z_from, z_to) : range for submodel in either cell or world coords

        **Optional keywords**:
            - *range_type* = 'cell', 'world' : define if subrange in cell ids (default) or real-world coordinates
        """
        range_type = kwds.get('range_type', 'cell')

        if not hasattr(self, 'boundaries_x'):
            self.determine_cell_boundaries()

        if range_type == 'world':
            # determine cells
            subrange[0] = np.argwhere(self.boundaries_x > subrange[0])[0][0]
            subrange[1] = np.argwhere(self.boundaries_x < subrange[1])[-1][0]
            subrange[2] = np.argwhere(self.boundaries_y > subrange[2])[0][0]
            subrange[3] = np.argwhere(self.boundaries_y < subrange[3])[-1][0]
            subrange[4] = np.argwhere(self.boundaries_z > subrange[4])[0][0]
            subrange[5] = np.argwhere(self.boundaries_z < subrange[5])[-1][0]

        # create a copy of the original grid
        import copy
        subgrid = copy.deepcopy(self)

        # extract grid
        subgrid.grid = self.grid[subrange[0]:subrange[1],
                                 subrange[2]:subrange[3],
                                 subrange[4]:subrange[5]]

        subgrid.nx = subrange[1] - subrange[0]
        subgrid.ny = subrange[3] - subrange[2]
        subgrid.nz = subrange[5] - subrange[4]

        # update extent
        subgrid.xmin = self.boundaries_x[subrange[0]]
        subgrid.xmax = self.boundaries_x[subrange[1]]
        subgrid.ymin = self.boundaries_y[subrange[2]]
        subgrid.ymax = self.boundaries_y[subrange[3]]
        subgrid.zmin = self.boundaries_z[subrange[4]]
        subgrid.zmax = self.boundaries_z[subrange[5]]

        subgrid.extent_x = subgrid.xmax - subgrid.xmin
        subgrid.extent_y = subgrid.ymax - subgrid.ymin
        subgrid.extent_z = subgrid.zmax - subgrid.zmin

        # update cell spacings
        subgrid.delx = self.delx[subrange[0]:subrange[1]]
        subgrid.dely = self.dely[subrange[2]:subrange[3]]
        subgrid.delz = self.delz[subrange[4]:subrange[5]]

        # now: update other attributes:
        subgrid.determine_cell_centers()
        subgrid.determine_cell_boundaries()
        subgrid.determine_cell_volumes()
        subgrid.determine_geology_ids()

        # finally: return subgrid
        return subgrid



# ******************************************************************************
#  Some additional helper functions
# ******************************************************************************

def combine_grids(G1, G2, direction, merge_type = 'keep_first', **kwds):
    """Combine two grids along one axis

    ..Note: this implementation assumes (for now) that the overlap is perfectly matching,
    i.e. grid cell sizes identical and at equal positions, or that they are perfectly adjacent!

    **Arguments**:
        - G1, G2 = GeoGrid : grids to be combined
        - direction = 'x', 'y', 'z': direction in which grids are combined
        - merge_type = method to combine grid:
            'keep_first' : keep elements of first grid (default)
            'keep_second' : keep elements of second grid
            'random' : randomly choose an element to retain
        ..Note: all other dimensions must be matching perfectly!!

    **Optional keywords**:
        - *overlap_analysis* = bool : perform a detailed analysis of the overlapping area, including
        mismatch. Also returns a second item, a GeoGrid with information on mismatch!

    **Returns**:
        - *G_comb* = GeoGrid with combined grid
        - *G_overlap* = Geogrid with analysis of overlap (of overlap_analysis=True)

    """
    overlap_analysis = kwds.get("overlap_analysis", False)
    # first step: determine overlap
    if direction == 'x':
        if G2.xmax > G1.xmax:
            overlap_min = G2.xmin
            overlap_max = G1.xmax
            # identifier alias for grids with higher/ lower values
            G_high = G2
            G_low = G1
        else:
            overlap_min = G1.xmin
            overlap_max = G2.xmax
            # identifier alias for grids with higher/ lower values
            G_high = G1
            G_low = G2

        # check if all other dimensions are perfectly matching
        if (G1.ymin != G2.ymin) or (G1.zmin != G2.zmin) or \
            (G1.ymax != G2.ymax) or (G1.zmax != G2.zmax):
            raise ValueError("Other dimensions (apart from %s) not perfectly matching! Check and try again!" % direction)

    elif direction == 'y':
        if G2.ymax > G1.ymax:
            overlap_min = G2.ymin
            overlap_max = G1.ymax
            # identifier alias for grids with higher/ lower values
            G_high = G2
            G_low = G1
        else:
            overlap_min = G1.ymin
            overlap_max = G2.ymax
            # identifier alias for grids with higher/ lower values
            G_high = G1
            G_low = G2

        # check if all other dimensions are perfectly matching
        if (G1.xmin != G2.xmin) or (G1.zmin != G2.zmin) or \
            (G1.xmax != G2.xmax) or (G1.zmax != G2.zmax):
            raise ValueError("Other dimensions (apart from %s) not perfectly matching! Check and try again!" % direction)

    elif direction == 'z':
        if G2.zmax > G1.zmax:
            overlap_min = G2.zmin
            overlap_max = G1.zmax
            # identifier alias for grids with higher/ lower values
            G_high = G2
            G_low = G1
        else:
            overlap_min = G1.zmin
            overlap_max = G2.zmax
            # identifier alias for grids with higher/ lower values
            G_high = G1
            G_low = G2

        # check if all other dimensions are perfectly matching
        if (G1.ymin != G2.ymin) or (G1.xmin != G2.xmin) or \
            (G1.ymax != G2.ymax) or (G1.xmax != G2.xmax):
            raise ValueError("Other dimensions (apart from %s) not perfectly matching! Check and try again!" % direction)

    overlap = overlap_max - overlap_min

    if overlap == 0:
        print("Grids perfectly adjacent")
    elif overlap < 0:
        raise ValueError("No overlap between grids! Check and try again!")
    else:
        print("Positive overlap in %s direction of %f meters" % (direction, overlap))

    # determine cell centers
    G1.determine_cell_centers()
    G2.determine_cell_centers()

    # intialise new grid
    G_comb = GeoGrid()
    # initialise overlap grid, if analyis performed
    if overlap_analysis:
        G_overlap = GeoGrid()


    if direction == 'x':
        pass
    elif direction == 'y':
        #=======================================================================
        # Perform overlap analysis
        #=======================================================================

        # initialise overlap grid with dimensions of overlap
        G_overlap.set_dimensions(dim = (G1.xmin, G1.xmax, overlap_min, overlap_max, G1.zmin, G1.zmax))
        G_low_ids = np.where(G_low.cell_centers_y > overlap_min)[0]
        G_high_ids = np.where(G_high.cell_centers_y < overlap_max)[0]
        delx = G1.delx
        dely = G_low.dely[G_low_ids]
        delz = G1.delz
        G_overlap.set_delxyz((delx, dely, delz))
        # check if overlap region is identical
        if not (len(G_low_ids) == len(G_high_ids)):
            raise ValueError("Overlap length not identical, please check and try again!")
        # now: determine overlap mismatch
        G_overlap.grid = G_low.grid[:,G_low_ids,:] - G_high.grid[:,G_high_ids,:]
        # for some very strange reason, this next step is necessary to enable the VTK
        # export with pyevtk - looks like a bug in pyevtk...
        G_overlap.grid = G_overlap.grid + np.zeros(G_overlap.grid.shape)
        #

        #=======================================================================
        # Set up combined grid
        #=======================================================================

        G_comb.set_dimensions(dim = (G1.xmin, G1.xmax, G_low.ymin, G_high.ymax, G1.zmin, G1.zmax))
        # combine dely arrays
        dely = np.hstack((G_low.dely[:G_low_ids[0]], G_high.dely))
        G_comb.set_delxyz((delx, dely, delz))

        #=======================================================================
        # Now merge grids
        #=======================================================================
        if merge_type == 'keep_first':
            if G1.ymax > G2.ymax:
                G_comb.grid = np.concatenate((G2.grid[:,:G_low_ids[0],:], G1.grid), axis=1)
            else:
                G_comb.grid = np.concatenate((G1.grid, G2.grid[:,:G_low_ids[0],:]), axis=1)

        elif merge_type == 'keep_second':
            pass
        elif merge_type == 'random':
            pass
        else:
            raise ValueError("Merge type %s not recognised! Please check and try again!" % merge_type)

    elif direction == 'z':
        pass




    # Return combined grid and results of overlap analysis, if determined
    if overlap_analysis:
        return G_comb, G_overlap
    else:
        return G_comb

def optimial_cell_increase(starting_cell_width, n_cells, width):
    """Determine an array with optimal cell width for a defined starting cell width,
    total number of cells, and total width

    Basically, this function optimised a factor between two cells to obtain a total
    width

    **Arguments**:
    	- *starting_cell_width* = float : width of starting/ inner cell
	- *n_cells* = int : total number of cells
	- *total_width* = float : total width (sum over all elements in array)

    **Returns**:
        del_array : numpy.ndarray with cell discretisations

    Note: optmisation with scipy.optimize - better (analytical?) methods might exist but
    I can't think of them at the moment
    """
    import scipy.optimize
    # define some helper functions

    def width_sum(inc_factor, inner_cell, n_cells, total_width):
        return sum(del_array(inc_factor, inner_cell, n_cells)) - total_width

    def del_array(inc_factor, inner_cell, n_cells):
        return np.array([inner_cell * inc_factor**i for i in range(n_cells)])

    # now the actual optimisation step:
    opti_factor = scipy.optimize.fsolve(width_sum, 1.1, (starting_cell_width, n_cells, width))

    # return the discretisation array
    return del_array(opti_factor, starting_cell_width, n_cells).flatten()



if __name__ == '__main__':
    pass
