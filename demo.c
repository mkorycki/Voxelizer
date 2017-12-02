//
//  demo.c
//  Voxelizer
//
//  Created by Michał Korycki.
//  Copyright © 2017 Michał Korycki. All rights reserved.
//

#include "voxelizer.h"

void print_voxels(vox_grid_params gridparams, int *voxels);

int main()
{
    vox_grid_params grid_params;
    
    grid_params.startX = -10.0;
    grid_params.startY = -10.0;
    grid_params.startZ = 0.0;
    
    grid_params.dx = 1.0;
    grid_params.dy = 1.0;
    grid_params.dz = 1.0;
    
    grid_params.nx = 20;
    grid_params.ny = 20;
    grid_params.nz = 5;
    
    
    
    vox_polygons poly;
    poly.npolygons = 2;
    poly.polygon = malloc(sizeof(vox_polygon)*poly.npolygons);
    
    poly.polygon[0].nparts = 1;
    poly.polygon[0].parts = malloc(sizeof(vox_polygon_part)*poly.polygon[0].nparts);
    
    poly.polygon[1].nparts = 1;
    poly.polygon[1].parts = malloc(sizeof(vox_polygon_part)*poly.polygon[1].nparts);
    
    
    poly.polygon[0].parts[0].npoints = 5;
    poly.polygon[0].parts[0].pointX = malloc(sizeof(double)*poly.polygon[0].parts[0].npoints);
    poly.polygon[0].parts[0].pointY = malloc(sizeof(double)*poly.polygon[0].parts[0].npoints);
    poly.polygon[0].parts[0].pointZ = malloc(sizeof(double)*poly.polygon[0].parts[0].npoints);
    
    poly.polygon[1].parts[0].npoints = 5;
    poly.polygon[1].parts[0].pointX = malloc(sizeof(double)*poly.polygon[1].parts[0].npoints);
    poly.polygon[1].parts[0].pointY = malloc(sizeof(double)*poly.polygon[1].parts[0].npoints);
    poly.polygon[1].parts[0].pointZ = malloc(sizeof(double)*poly.polygon[1].parts[0].npoints);
    
    
    poly.polygon[0].parts[0].pointX[0] = 5.0;
    poly.polygon[0].parts[0].pointY[0] = 5.0;
    poly.polygon[0].parts[0].pointZ[0] = 3.0;
    
    poly.polygon[0].parts[0].pointX[1] = 5.0;
    poly.polygon[0].parts[0].pointY[1] =-5.0;
    poly.polygon[0].parts[0].pointZ[1] = 3.0;
    
    poly.polygon[0].parts[0].pointX[2] =-5.0;
    poly.polygon[0].parts[0].pointY[2] =-5.0;
    poly.polygon[0].parts[0].pointZ[2] = 3.0;
    
    poly.polygon[0].parts[0].pointX[3] =-5.0;
    poly.polygon[0].parts[0].pointY[3] = 5.0;
    poly.polygon[0].parts[0].pointZ[3] = 3.0;
    
    poly.polygon[0].parts[0].pointX[4] = 5.0;
    poly.polygon[0].parts[0].pointY[4] = 5.0;
    poly.polygon[0].parts[0].pointZ[4] = 3.0;
    
    
    poly.polygon[1].parts[0].pointX[0] = 5.0;
    poly.polygon[1].parts[0].pointY[0] = 5.0;
    poly.polygon[1].parts[0].pointZ[0] = 1.0;
    
    poly.polygon[1].parts[0].pointX[1] =-5.0;
    poly.polygon[1].parts[0].pointY[1] = 5.0;
    poly.polygon[1].parts[0].pointZ[1] = 1.0;
    
    poly.polygon[1].parts[0].pointX[2] =-5.0;
    poly.polygon[1].parts[0].pointY[2] =-5.0;
    poly.polygon[1].parts[0].pointZ[2] = 1.0;
    
    poly.polygon[1].parts[0].pointX[3] = 5.0;
    poly.polygon[1].parts[0].pointY[3] =-5.0;
    poly.polygon[1].parts[0].pointZ[3] = 1.0;
    
    poly.polygon[1].parts[0].pointX[4] = 5.0;
    poly.polygon[1].parts[0].pointY[4] = 5.0;
    poly.polygon[1].parts[0].pointZ[4] = 1.0;
    
    
    vox_poly_prepare(poly);
    
    
    vox_rays rays;
    rays.nrays = 2;
    rays.vX = malloc(sizeof(double)*rays.nrays);
    rays.vY = malloc(sizeof(double)*rays.nrays);
    rays.vZ = malloc(sizeof(double)*rays.nrays);
    
    rays.vX[0] = 0.0;
    rays.vY[0] = 0.0;
    rays.vZ[0] = 1.0;
    
    rays.vX[1] = 0.0;
    rays.vY[1] = 0.0;
    rays.vZ[1] =-1.0;
    
    int *grids = malloc(sizeof(int)*grid_params.nx
                                   *grid_params.ny
                                   *grid_params.nz
                                   *rays.nrays);
    
    
    voxelize (poly, grid_params, rays, grids);
    
    
    int *voxels = malloc(sizeof(int)*grid_params.nx
                                    *grid_params.ny
                                    *grid_params.nz);
    
    vox_postproc(rays.nrays, grid_params, grids, voxels);
    
    print_voxels(grid_params, voxels);
    
    return 0;
}

void print_voxels(vox_grid_params gridparams, int *voxels)
{
    
    int i, ix, iy, iz;
    
    for (iz = 0; iz < gridparams.nz; iz++) {

        printf("\nlayer: %d\n",iz);

        for (iy = gridparams.ny - 1; iy >= 0; iy--) {
            for (ix=0; ix < gridparams.nx; ix++) {
                i = iz * gridparams.ny * gridparams.nx +
                    iy * gridparams.nx +
                    ix;

                printf("%d ",voxels[i]);

            }

            printf("\n");

        }
    }
    
}
