//sequence for clean computation with several gmapi commands

#include <stdio.h>

int compute(void) 
{
  // clear the list of sections/series/faults that are used to compute the geological model
  gmapi_clearThisModel();

  // Define the sections we want to include in the computation of the model
  // The keyword "all" includes all sections present in the project. 
  if(gmapi_addSectionToModel("all") !=0 )
  {
    printf("Error while adding section...\n");
	for (int i = 1; i < 80; i++) {printf("*");}
	printf("\n\n\tProblem with return value!\n\n\n");
	//return = -1;
  }

  // Define the series we want to include in the computation of the model
  // The keyword "all" includes all series present in the project. 
  if(gmapi_addSeriesToModel("all") !=0 )
  {
    printf("Error while adding series...\n");
	for (int i = 1; i < 80; i++) {printf("*");}
	printf("\n\n\tProblem with return value!\n\n\n");
	//return = -1;
  }

  // Add faults to model
  if(gmapi_addFaultToModel("all") !=0 )
  {
    printf("Error while adding faults...\n");
	for (int i = 1; i < 80; i++) {printf("*");}
	printf("\n\n\tProblem with return value!\n\n\n");
	//return = -1;
  }

  // compute the geological model with the data defined above
  if(gmapi_computeThisModel() !=0 )
  {
    printf("Error computing model!\n");
	for (int i = 1; i < 80; i++) {printf("*");}
	printf("\n\n\tProblem with return value!\n\n\n");
	//return = -1;
  }

  return 0;
}
