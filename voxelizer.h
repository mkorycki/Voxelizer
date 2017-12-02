//
//  voxelizer.h
//  Voxelizer
//
//  Created by Michał Korycki.
//  Copyright © 2017 Michał Korycki. All rights reserved.
//

#ifndef VOXELIZER
#define VOXELIZER

#include "vox_grid_params.h"
#include "vox_polygons.h"
#include "vox_rays.h"

#include <stdio.h>
#include <stdlib.h>

void vox_poly_prepare(vox_polygons polygons);

int voxelize (vox_polygons polygons, vox_grid_params gridparams,
              vox_rays rays, int *grids);

void vox_postproc(int nrays, vox_grid_params gridparams, int *grids, int *voxels);

#endif
