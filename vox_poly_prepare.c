//
//  vox_poly_prepare.c
//  Voxelizer
//
//  Created by Michał Korycki.
//  Copyright © 2017 Michał Korycki. All rights reserved.
//

#include "vox_polygons.h"
#include "vox_grid_params.h"
#include "vox_tools.h"
#include <stdlib.h>

void vox_poly_prepare(vox_polygons polygons){
    int ipolygon;
    int ipart;
    int ipoint;
    double MinX, MaxX;
    double MinY, MaxY;
    double MinZ, MaxZ;
    char *pnts;
    char Maxpnts;
    
    int ipointp, ipointm;
    double uX, uY, uZ, vX, vY, vZ;
    double crossX, crossY, crossZ;
    double mag;
    
    for (ipolygon=0; ipolygon<polygons.npolygons; ipolygon++) {
        for (ipart=0; ipart<polygons.polygon[ipolygon].nparts; ipart++) {
            pnts = malloc(sizeof(char)*polygons.polygon[ipolygon].parts[ipart].npoints);
            for (ipoint=0; ipoint<polygons.polygon[ipolygon].parts[ipart].npoints; ipoint++){
                pnts[ipoint] = 0;
            }
            
            
            vox_find_minmax(polygons.polygon[ipolygon].parts[ipart].npoints,
                            polygons.polygon[ipolygon].parts[ipart].pointX,
                            polygons.polygon[ipolygon].parts[ipart].pointY,
                            polygons.polygon[ipolygon].parts[ipart].pointZ,
                            &MinX, &MaxX, &MinY, &MaxY, &MinZ, &MaxZ);
            
            
            if (ipart==0) {
                polygons.polygon[ipolygon].MinX = MinX;
                polygons.polygon[ipolygon].MaxX = MaxX;
                polygons.polygon[ipolygon].MinY = MinY;
                polygons.polygon[ipolygon].MaxY = MaxY;
                polygons.polygon[ipolygon].MinZ = MinZ;
                polygons.polygon[ipolygon].MaxZ = MaxZ;
            }
            
            Maxpnts = 0;
            for (ipoint=0; ipoint<polygons.polygon[ipolygon].parts[ipart].npoints; ipoint++){
                if (polygons.polygon[ipolygon].parts[ipart].pointX[ipoint]==MinX) {
                    pnts[ipoint]++;
                }
                if (polygons.polygon[ipolygon].parts[ipart].pointX[ipoint]==MaxX) {
                    pnts[ipoint]++;
                }
                if (polygons.polygon[ipolygon].parts[ipart].pointY[ipoint]==MinY) {
                    pnts[ipoint]++;
                }
                if (polygons.polygon[ipolygon].parts[ipart].pointY[ipoint]==MaxY) {
                    pnts[ipoint]++;
                }
                if (polygons.polygon[ipolygon].parts[ipart].pointZ[ipoint]==MinZ) {
                    pnts[ipoint]++;
                }
                if (polygons.polygon[ipolygon].parts[ipart].pointZ[ipoint]==MaxZ) {
                    pnts[ipoint]++;
                }
                if (pnts[ipoint]>Maxpnts) {
                    Maxpnts = pnts[ipoint];
                }
            }
            
            for (ipoint=0; ipoint<polygons.polygon[ipolygon].parts[ipart].npoints; ipoint++){
                if (pnts[ipoint]==Maxpnts) {
                    ipointm = ipoint-1;
                    ipointp = ipoint+1;
                    if (ipoint==0 || ipoint==polygons.polygon[ipolygon].parts[ipart].npoints-1) {
                        ipointp = 1;
                        ipointm = polygons.polygon[ipolygon].parts[ipart].npoints-2;
                    }
                    
                    uX = polygons.polygon[ipolygon].parts[ipart].pointX[ipointm] -
                    polygons.polygon[ipolygon].parts[ipart].pointX[ipoint];
                    uY = polygons.polygon[ipolygon].parts[ipart].pointY[ipointm] -
                    polygons.polygon[ipolygon].parts[ipart].pointY[ipoint];
                    uZ = polygons.polygon[ipolygon].parts[ipart].pointZ[ipointm] -
                    polygons.polygon[ipolygon].parts[ipart].pointZ[ipoint];
                    
                    vX = polygons.polygon[ipolygon].parts[ipart].pointX[ipointp] -
                    polygons.polygon[ipolygon].parts[ipart].pointX[ipoint];
                    vY = polygons.polygon[ipolygon].parts[ipart].pointY[ipointp] -
                    polygons.polygon[ipolygon].parts[ipart].pointY[ipoint];
                    vZ = polygons.polygon[ipolygon].parts[ipart].pointZ[ipointp] -
                    polygons.polygon[ipolygon].parts[ipart].pointZ[ipoint];
                    
                    vox_cross_prod(uX, uY, uZ, vX, vY, vZ, &crossX, &crossY, &crossZ);
                    
                    mag = vox_magnitude(crossX, crossY, crossZ);
                    
                    if (mag==0.0) continue;
                    
                    if (ipart==0) {
                        polygons.polygon[ipolygon].A = crossX/mag;
                        polygons.polygon[ipolygon].B = crossY/mag;
                        polygons.polygon[ipolygon].C = crossZ/mag;
                        polygons.polygon[ipolygon].parts[ipart].isfull = 1;
                    }
                    
                    if (ipart!=0 &&
                        (polygons.polygon[ipolygon].A==-crossX/mag ||
                         polygons.polygon[ipolygon].B==-crossY/mag ||
                         polygons.polygon[ipolygon].C==-crossZ/mag )) {
                        polygons.polygon[ipolygon].parts[ipart].isfull = 0;
                    }
                    else{
                        polygons.polygon[ipolygon].parts[ipart].isfull = 1;
                    }
                }
            }
            
            
            free(pnts);
            
            polygons.polygon[ipolygon].D =
            polygons.polygon[ipolygon].parts[0].pointX[0]*polygons.polygon[ipolygon].A +
            polygons.polygon[ipolygon].parts[0].pointY[0]*polygons.polygon[ipolygon].B +
            polygons.polygon[ipolygon].parts[0].pointZ[0]*polygons.polygon[ipolygon].C;
            
        }
    }
    
}
