//
//  vox_postproc.c
//  Voxelizer
//
//  Created by Michał Korycki.
//  Copyright © 2017 Michał Korycki. All rights reserved.
//

#include "vox_grid_params.h"

void vox_postproc(int nrays, vox_grid_params gridparams, int *grids, int *voxels){

    int ix, iy, iz, iray;
    int nodds;
    float floatnodds;
    float floatnrays;
    
    int gincrement;
    int vincrement;

    floatnrays = (float) nrays;
        
    for (ix = 0; ix<gridparams.nx; ix++){
        for (iy = 0; iy<gridparams.ny; iy++){
            for (iz = 0; iz<gridparams.nz; iz++){
                vincrement = ix + gridparams.nx*iy + gridparams.nx*gridparams.ny*iz;
                nodds = 0;
                for (iray = 0; iray<nrays; iray++){
                    gincrement = ix + gridparams.nx*iy + gridparams.nx*gridparams.ny*iz +
                                 gridparams.nx*gridparams.ny*gridparams.nz*iray;
                    
                    nodds += grids[gincrement]%2;
                }
                floatnodds = (float) nodds;
                if (floatnodds/floatnrays>=0.5) {
                    voxels[vincrement] = 1;
                }
                else{
                    voxels[vincrement] = 0;
                }
            }
        }
    }
}
