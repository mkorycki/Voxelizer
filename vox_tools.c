//
//  vox_tools.c
//  Voxelizer
//
//  Created by Michał Korycki.
//  Copyright © 2017 Michał Korycki. All rights reserved.
//

#include <math.h>


void vox_cross_prod(double uX, double uY, double uZ,
                double vX, double vY, double vZ,
                double *crossX, double *crossY, double *crossZ){
    
    *crossX = uY * vZ - uZ * vY;
    *crossY = uZ * vX - uX * vZ;
    *crossZ = uX * vY - uY * vX;
    
}




void vox_find_minmax(int npoints, double *inX, double *inY, double *inZ,
                 double *MinX, double *MaxX,
                 double *MinY, double *MaxY,
                 double *MinZ, double *MaxZ){
    int ipoint;
    
    *MinX = inX[0];
    *MaxX = inX[0];
    *MinY = inY[0];
    *MaxY = inY[0];
    *MinZ = inZ[0];
    *MaxZ = inZ[0];
    
    for (ipoint=1; ipoint<npoints; ipoint++) {
        if (inX[ipoint]>*MaxX) {
            *MaxX = inX[ipoint];
        }
        if (inX[ipoint]<*MinX) {
            *MinX = inX[ipoint];
        }
        
        if (inY[ipoint]>*MaxY) {
            *MaxY = inY[ipoint];
        }
        if (inY[ipoint]<*MinY) {
            *MinY = inY[ipoint];
        }
        
        if (inZ[ipoint]>*MaxZ) {
            *MaxZ = inZ[ipoint];
        }
        if (inZ[ipoint]<*MinZ) {
            *MinZ = inZ[ipoint];
        }
    }
    
}


void vox_rotate(double angle, double iX, double iY, double iZ, double *oX, double *oY, double *oZ){
    *oX = iX * cos(angle) - iY * sin(angle);
    *oY = iX * sin(angle) + iY * cos(angle);
    *oZ = iZ;
}




double vox_magnitude(double uX, double uY, double uZ){
    return ( sqrt(uX*uX + uY*uY + uZ*uZ) );
}
