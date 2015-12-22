#include <iostream>
#include "geomodeller-c-api.h" //import geomodeller api
#include "geomodeller_compute.h" //sequence for clean computation with several gmapi commands
#include <stdio.h>
#include <cstdlib>
#include <omp.h>

//#include <unistd.h>


#define DLLEXPORT extern "C" __declspec(dllexport)

//******************************************************************
//***********************Get Model Boundaries***********************
//******************************************************************
DLLEXPORT int* get_model_bounds(char* filename)
{

  // parameter declaration
  double xmin, ymin, zmin, xmax, ymax, zmax;


  // set up licensing - querying models is allowed for anyone, write access needs a valid licence
  printf("\n\t***** Initialize Read API *****\n\n");
  gmapi_initialiseReadAPI();
  printf("\n\t***** Initialize Write API *****\n\n");
  gmapi_initialiseWriteAPINoArgs();


  // open the project without displaying it
  printf("\n\t***** Open Project *****\n\n");
  printf("Project to open is: %s\n", filename);
  if(gmapi_openProjectWithoutDisplay(filename) != 0 )
  {
    for (int i = 1; i < 80; i++) {printf("*");}
    printf("\n\n\tProblem with return value!\n\n\n");
    //return = -1;
  }


  // get the model bounds and print to stdout
  printf("\n\t***** Read Model Boundaries *****\n\n");
  gmapi_getModelBounds(xmin, ymin, zmin, xmax, ymax, zmax);
  int* boundaries = new int[6];
  boundaries[0] = xmin;
  boundaries[1] = xmax;
  boundaries[2] = ymin;
  boundaries[3] = ymax;
  boundaries[4] = zmin;
  boundaries[5] = zmax;

  return boundaries;
}


//******************************************************************
//*****************Compute and Export Regular Grid******************
//******************************************************************
DLLEXPORT int* compute_regular_grid(char* filename, int disc[])
{
  // parameter declaration
  double xmin, ymin, zmin, xmax, ymax, zmax;


  // set up licensing - querying models is allowed for anyone, write access needs a valid licence
  printf("\n\t***** Initialize Read API *****\n\n");
  gmapi_initialiseReadAPI();
  printf("\n\t***** Initialize Write API *****\n\n");
  gmapi_initialiseWriteAPINoArgs();


  // open the project without displaying it
  printf("\n\t***** Open Project *****\n\n");
  printf("Project to open is: %s\n", filename);
  if(gmapi_openProjectWithoutDisplay(filename) != 0 )
  {
    for (int i = 1; i < 80; i++) {printf("*");}
    printf("\n\n\tProblem with return value!\n\n\n");
    //return = -1;
  }


  // get the model bounds and print to stdout
  printf("\n\t***** Read Model Boundaries *****\n\n");
  gmapi_getModelBounds(xmin, ymin, zmin, xmax, ymax, zmax);
  printf("ModelBounds MIN xmin=%lf ymin=%lf zmin=%lf\n", xmin, ymin, zmin);
  printf("ModelBounds MAX xmax=%lf ymax=%lf zmax=%lf\n\n", xmax, ymax, zmax);


  // Determine increment size in each direction
  double dX;
  double dY;
  double dZ;
  dX = (xmax - xmin) / float(disc[0]);
  dY = (ymax - ymin) / float(disc[1]);
  dZ = abs(zmin - zmax) / float(disc[2]);
  // printf("%lf\n",atof(argv[4]));
  int nodes_x = int(disc[0]);
  int nodes_y = int(disc[1]);
  int nodes_z = int(disc[2]);
  printf("Number of cells\t: nx = %d, ny = %d, nz = %d\n", nodes_x, nodes_y, nodes_z);
  printf("Cell sizes\t: dX = %.2lf, dY = %.2lf, dZ = %.2lf\n", dX, dY, dZ);
  printf("\n\n");


  // compute the model, see below
  // function can't be defined in functions, therefore compute() is loaded from compute.h
  printf("\n\t***** Compute model *****\n\n");
  if(compute() != 0)
  {
    for (int i = 1; i < 80; i++) {printf("*");}
    printf("\n\n\tProblem with return value!\n\n\n");
    //return = -1;
  }
  else
  {
    printf("Model computation successful!\n\n");
  }


  // Open file to export lithological grid
  printf("\n\t***** Export lithology on regular grid *****\n\n");
  FILE *fp;
  fp = fopen("grid_export//exported_regular_grid.txt","w");

  int n_total = nodes_x * nodes_y * nodes_z;
  int j = 1;
  int k = 0;
  int p = 0;
  int q = 0;
  int r = 0;
  double X = xmin;
  double Y = ymin;
  double Z = zmin;
  int litho=-1;
  int iErrorLitho = -1;
  int* formation = new int[n_total];


  while(Z < zmax){
    while(Y < ymax){
        while(X < xmax){
            // printf("%f, %f, %f\n", X, Y, Z);
            litho=-1;
            iErrorLitho = -1;
            iErrorLitho = gmapi_getComputedLithologyXYZ(X+dX/2.0,Y+dY/2.0,Z+dZ/2.0,litho);
            // iErrorLitho = GetComputedLithologyXYZ(X+dX/2.0,Y+dY/2.0,Z+dZ/2.0,litho);
            // printf("at Z=%f, lithology=%d and error is %d\n",Z,litho,iErrorLitho);
            // printf("%d,",litho);
            fprintf(fp, "%d,",litho);
            //formation[j-1] = litho;
            formation[r*nodes_x*nodes_y + p*nodes_x + q] = litho;
            printf("Exporting cell %8d of %8d (%7.3f %%)\r",j,n_total,double(j)/n_total*100);
            X += dX;
            j += 1;
            p += 1;
        }
        X = xmin;
        Y += dY;
        q += 1;
        p = 0;
    }
    // printf("\n");
    fprintf(fp, "\n");
    Z += dZ;
    k ++;
    Y = ymin;
    r += 1;
    q = 0;
  }
  printf("\n");

  fclose(fp);

  // open file to save delx, dely, delz
  FILE *fp2;
  fp2 = fopen("grid_export//delxyz.txt","w");
  fprintf(fp2, "%d*%.2f\n",nodes_x, dX);
  fprintf(fp2, "%d*%.2f\n",nodes_y, dY);
  for (int iz = 0; iz < nodes_z; iz++){
    fprintf(fp2, "%.2f,",dZ);
  }
  fclose(fp2);

  // open file to save project dimensions (because they are not stored in SHEMAT...)
  FILE *fp3;
  fp3 = fopen("grid_export//project_dimensions.txt","w");
  fprintf(fp3, "xmin, xmax, ymin, ymax, zmin, zmax\n");
  fprintf(fp3, "%.1f, %.1f, %.1f, %.1f, %.1f, %.1f\n",xmin,xmax,ymin,ymax,zmin,zmax);
  fclose(fp3);


  return formation;
}


//******************************************************************
//****************Compute and Return Irregular Grid*****************
//******************************************************************
DLLEXPORT int* compute_irregular_grid(char* filename, double coord[], int coord_len, int cores)
{
  // parameter declaration
  double xmin, ymin, zmin, xmax, ymax, zmax;


  // set up licensing - querying models is allowed for anyone, write access needs a valid licence
  printf("\n\t***** Initialize Read API *****\n\n");
  gmapi_initialiseReadAPI();
  printf("\n\t***** Initialize Write API *****\n\n");
  gmapi_initialiseWriteAPINoArgs();


  // open the project without displaying it
  printf("\n\t***** Open Project *****\n\n");
  printf("Project to open is: %s\n", filename);
  if(gmapi_openProjectWithoutDisplay(filename) != 0 )
  {
    for (int i = 1; i < 80; i++) {printf("*");}
    printf("\n\n\tProblem with return value!\n\n\n");
    //return = -1;
  }


  // compute the model, see below
  // function can't be defined in functions, therefore compute() is loaded from compute.h
  printf("\n\t***** Compute model *****\n\n");
  if(compute() != 0)
  {
    for (int i = 1; i < 80; i++) {printf("*");}
    printf("\n\n\tProblem with return value!\n\n\n");
    //return = -1;
  }
  else
  {
    printf("Model computation successful!\n\n");
  }


  // Save model (disable because it makes problems in mpi mode)
  /*char *path=NULL;
  size_t size;
  path=getcwd(path,size);
  printf("\n\n\tCurrent path is: %s",path);
  printf("\n\t***** Save model *****\n\n");
  gmapi_saveProjectAs(filename);
  */

  // Open file to export lithological grid with given coordinates
  printf("\n\t***** Export lithology on irregular grid *****\n\n");
  int* formation = new int[coord_len/3];
  int form_value;
  int litho = -1;

  omp_set_num_threads(cores);

#pragma omp parallel for
  printf("Number of cores are: %s\n",cores);

  for(int m=0;m<coord_len/3;m++)
  {
    form_value = gmapi_getComputedLithologyXYZ(coord[0+3*m],coord[1+3*m],coord[2+3*m],litho);
    formation[m] = litho;
  }
  printf("The lenght of litho is: %s\n", litho);

  return formation;
}

notepad.cc



























DLLEXPORT int* init_api(char* filename)
{
  // set up licensing - querying models is allowed for anyone, write access needs a valid licence
  printf("\n\t***** Initialize Read API *****\n\n");
  gmapi_initialiseReadAPI();
  printf("\n\t***** Initialize Write API *****\n\n");
  gmapi_initialiseWriteAPINoArgs();


  // open the project without displaying it
  printf("\n\t***** Open Project *****\n\n");
  printf("Project to open is: %s\n", filename);
  if(gmapi_openProjectWithoutDisplay(filename) != 0 )
  {
    for (int i = 1; i < 80; i++) {printf("*");}
    printf("\n\n\tProblem with return value!\n\n\n");
    //return = -1;
  }


  // compute the model, see below
  // function can't be defined in functions, therefore compute() is loaded from compute.h
  printf("\n\t***** Compute model *****\n\n");
  if(compute() != 0)
  {
    for (int i = 1; i < 80; i++) {printf("*");}
    printf("\n\n\tProblem with return value!\n\n\n");
    //return = -1;
  }
  else
  {
    printf("Model computation successful!\n\n");
  }

 }

DLLEXPORT int* compute_irregular_grid(int coord_x, int coord_y, int coord_z)
{

  printf("\n\t***** Export lithology on irregular grid *****\n\n");
  int form_value;
  int litho = -1;
  form_value = gmapi_getComputedLithologyXYZ( coord_x, coord_y, coord_z,litho);
  printf("The lenght of litho is: %s\n", litho);

  return litho;
}





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


        # call cpp function
                    #Detection of operative system:
        if platform.system() == "Linux":
            lib = ctypes.CDLL('./libgeomod.so') #linux
        elif platform.system() == "Windows":
            lib = ctypes.windll.LoadLibrary(os.path.dirname(os.path.abspath(__file__)) + os.path.sep +"libgeomodwin1.dll") #windows
        else:
            print("Your operative system is not supported")

        lib.compute_irregular_grid.restype = ndpointer(dtype=ctypes.c_int, shape=(coord_len/3,))
        # This is the function wich needs GPU!!!
        lib.init_api(filename_ctypes)

        formations_raw = []
        # check if cell centers are defined - if not, do so!
        if not hasattr(self, 'cell_centers_x'):
            self.determine_cell_centers()
            for k in range(self.nz):
                for j in range(self.ny):
                    for i in range(self.nx):
                        ids.append((i,j,k))
                        formations_raw.append(lib.compute_irregular_grid(self.cell_centers_x[i], self.cell_centers_y[j], self.cell_centers_z[k]))


        # re-sort formations into array
        for i in range(len(formations_raw)):
            self.grid[ids[i][0],ids[i][1],ids[i][2]] = formations_raw[i]
