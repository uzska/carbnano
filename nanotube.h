#ifndef NANOTUBE_H_
#define NANOTUBE_H_

/*
 * N x 4 matrix. The four different types of elements are 
 * pos of face1, pos of face2, radius, bin #
 */
gsl_matrix *generate_All_Nanotubes(FILE *f);

#endif
