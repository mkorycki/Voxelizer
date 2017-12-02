//
//  vox_grid_params.h
//  Voxelizer
//
//  Created by Michał Korycki.
//  Copyright © 2017 Michał Korycki. All rights reserved.
//

#ifndef VOX_GRID_PARAMS
#define VOX_GRID_PARAMS

typedef struct{
   double startX;
   double startY;
   double startZ;
   double dx;
   double dy;
   double dz;
   int nx;
   int ny;
   int nz;
}vox_grid_params;

#endif
