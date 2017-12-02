//
//  vox_tools.h
//  Voxelizer
//
//  Created by Michał Korycki.
//  Copyright © 2017 Michał Korycki. All rights reserved.
//

#ifndef VOX_TOOLS
#define VOX_TOOLS

void vox_cross_prod(double uX, double uY, double uZ,
                double vX, double vY, double vZ,
                double *crossX, double *crossY, double *crossZ);

void vox_find_minmax(int npoints, double *inX, double *inY, double *inZ,
                 double *MinX, double *MaxX,
                 double *MinY, double *MaxY,
                 double *MinZ, double *MaxZ);

void vox_rotate(double angle, double iX, double iY, double iZ, double *oX, double *oY, double *oZ);

double vox_magnitude(double uX, double uY, double uZ);

#endif
