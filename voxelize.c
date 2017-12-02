//
//  voxelize.c
//  Voxelizer
//
//  Created by Michał Korycki.
//  Copyright © 2017 Michał Korycki. All rights reserved.
//

#include "vox_polygons.h"
#include "vox_grid_params.h"
#include "vox_rays.h"
#include "vox_tools.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>





typedef struct{
    double *pointK;
    double *pointL;
}projection_2d_part;

typedef struct{
    projection_2d_part *parts;
    double Kx, Ky, Kz;
    double Lx, Ly, Lz;
}projection_2d_polygon;





#define EPSILON 0.0000001


int between(double b1, double b2, double p){
    if (b1>b2) {
        if (p<=b1&&p>b2) {
            return(1);
        }
    }
    if (b1<b2) {
        if (p>=b1&&p<b2) {
            return(1);
        }
    }
    return(0);
}

                 
int intersection(double  x, double  y, double  z,
                 double vx, double vy, double vz,
                 vox_polygon *polygon, projection_2d_polygon *polygon2d){
    
    int isinside;
    int ipart;
    int ipoint;
    double t; //distance from point to polygon
    double Px, Py, Pz; //point on plane
    double Pk, Pl;
    double dkdl, dldk;
    double dk, dl;
    double tmpK, tmpL;
    int interKp, interKm, interLp, interLm;
    int interPoints;
    
    double vect_norm_angle;
    
    vect_norm_angle = acos(vx*polygon->A+vy*polygon->B+vz*polygon->C);
    //Ray vectors and plane normal are always 1 unit lenght
    if (vect_norm_angle<=M_PI_2+EPSILON && vect_norm_angle>=M_PI_2-EPSILON) {
        return(0);
    }
    
    if (vect_norm_angle<=M_PI/2.0) {
        t = -(- x*polygon->A -  y*polygon->B -  z*polygon->C + polygon->D)/
             (-vx*polygon->A - vy*polygon->B - vz*polygon->C);
    }
    else{
        t = -( x*polygon->A +  y*polygon->B +  z*polygon->C - polygon->D)/
             (vx*polygon->A + vy*polygon->B + vz*polygon->C);
    }
    
    if (t<0.0) return (0);


    Px = x + t*vx;
    Py = y + t*vy;
    Pz = z + t*vz;
    
    if (Px<polygon->MinX-EPSILON || Px>polygon->MaxX+EPSILON ||
        Py<polygon->MinY-EPSILON || Py>polygon->MaxY+EPSILON ||
        Pz<polygon->MinZ-EPSILON || Pz>polygon->MaxZ+EPSILON) {
        return(0);
    }
    
    Pk=
    (Px-polygon->parts[0].pointX[0])*polygon2d->Kx+
    (Py-polygon->parts[0].pointY[0])*polygon2d->Ky+
    (Pz-polygon->parts[0].pointZ[0])*polygon2d->Kz;
    
    Pl=
    (Px-polygon->parts[0].pointX[0])*polygon2d->Lx+
    (Py-polygon->parts[0].pointY[0])*polygon2d->Ly+
    (Pz-polygon->parts[0].pointZ[0])*polygon2d->Lz;
    

    
    isinside = 0;
    for (ipart = 0; ipart<polygon->nparts; ipart++) {
        
        interKp = 0;
        interKm = 0;
        interLp = 0;
        interLm = 0;
        interPoints = 0;
        
        for (ipoint = 0; ipoint<polygon->parts[ipart].npoints - 1; ipoint++) {
            if (between(polygon2d->parts[ipart].pointK[ipoint],polygon2d->parts[ipart].pointK[ipoint+1],Pk)) {
                dk = polygon2d->parts[ipart].pointK[ipoint+1] - polygon2d->parts[ipart].pointK[ipoint];
                dl = polygon2d->parts[ipart].pointL[ipoint+1] - polygon2d->parts[ipart].pointL[ipoint];
                dldk = dl/dk;
                
                tmpL = polygon2d->parts[ipart].pointL[ipoint] + (Pk-polygon2d->parts[ipart].pointK[ipoint]) * dldk;
                if (tmpL>Pl) {
                    interLp++;
                }
                if (tmpL<Pl) {
                    interLm++;
                }
            }
            if (between(polygon2d->parts[ipart].pointL[ipoint],polygon2d->parts[ipart].pointL[ipoint+1],Pl)) {
                dk = polygon2d->parts[ipart].pointK[ipoint+1] - polygon2d->parts[ipart].pointK[ipoint];
                dl = polygon2d->parts[ipart].pointL[ipoint+1] - polygon2d->parts[ipart].pointL[ipoint];
                dkdl = dk/dl;
                
                tmpK = polygon2d->parts[ipart].pointK[ipoint] + (Pl-polygon2d->parts[ipart].pointL[ipoint]) * dkdl;
                if (tmpK>Pk) {
                    interKp++;
                }
                if (tmpK<Pk) {
                    interKm++;
                }
            }
        }
        if (interLp%2==1) {
            interPoints++;
        }
        if (interLm%2==1) {
            interPoints++;
        }
        if (interKp%2==1) {
            interPoints++;
        }
        if (interKm%2==1) {
            interPoints++;
        }
        if (polygon->parts[ipart].isfull==1) {
            if (interPoints>0) {
                isinside=1;
            }
        }
        if (polygon->parts[ipart].isfull==0) {
            if (interPoints>0) {
                isinside=0;
            }
        }
    }
    
    return(isinside);

}


void make_2d_projections(vox_polygons polygons, projection_2d_polygon *polygon2d){
    
    int ipolygon;
    int ipart;
    int ipoint;
    double mag, magtmp;
    double tmpKx, tmpKy, tmpKz;
    double zeroX, zeroY, zeroZ;
    
    
    
    
    for (ipolygon = 0; ipolygon<polygons.npolygons; ipolygon++){
        mag = 0;
        for (ipoint=1; ipoint<polygons.polygon[ipolygon].parts[0].npoints-1; ipoint++) {
            tmpKx = polygons.polygon[ipolygon].parts[0].pointX[ipoint] -
                                     polygons.polygon[ipolygon].parts[0].pointX[0];
            tmpKy = polygons.polygon[ipolygon].parts[0].pointY[ipoint] -
                                     polygons.polygon[ipolygon].parts[0].pointY[0];
            tmpKz = polygons.polygon[ipolygon].parts[0].pointZ[ipoint] -
                                     polygons.polygon[ipolygon].parts[0].pointZ[0];
            
            magtmp = vox_magnitude(tmpKx, tmpKy, tmpKz);

            if (magtmp>mag) {
                mag = magtmp;
                polygon2d[ipolygon].Kx = tmpKx;
                polygon2d[ipolygon].Ky = tmpKy;
                polygon2d[ipolygon].Kz = tmpKz;
            }
        }
        polygon2d[ipolygon].Kx /= mag;
        polygon2d[ipolygon].Ky /= mag;
        polygon2d[ipolygon].Kz /= mag;
        
        vox_cross_prod(polygons.polygon[ipolygon].A,polygons.polygon[ipolygon].B,polygons.polygon[ipolygon].C,
                       polygon2d[ipolygon].Kx, polygon2d[ipolygon].Ky, polygon2d[ipolygon].Kz,
                       &(polygon2d[ipolygon].Lx), &(polygon2d[ipolygon].Ly), &(polygon2d[ipolygon].Lz));

        
        
        zeroX = polygons.polygon[ipolygon].parts[0].pointX[0];
        zeroY = polygons.polygon[ipolygon].parts[0].pointY[0];
        zeroZ = polygons.polygon[ipolygon].parts[0].pointZ[0];
        for (ipart=0; ipart<polygons.polygon[ipolygon].nparts; ipart++) {
            for (ipoint=0; ipoint<polygons.polygon[ipolygon].parts[ipart].npoints; ipoint++) {
                polygon2d[ipolygon].parts[ipart].pointK[ipoint] =
                (polygons.polygon[ipolygon].parts[ipart].pointX[ipoint]-zeroX)*polygon2d[ipolygon].Kx+
                (polygons.polygon[ipolygon].parts[ipart].pointY[ipoint]-zeroY)*polygon2d[ipolygon].Ky+
                (polygons.polygon[ipolygon].parts[ipart].pointZ[ipoint]-zeroZ)*polygon2d[ipolygon].Kz;
                polygon2d[ipolygon].parts[ipart].pointL[ipoint] =
                (polygons.polygon[ipolygon].parts[ipart].pointX[ipoint]-zeroX)*polygon2d[ipolygon].Lx+
                (polygons.polygon[ipolygon].parts[ipart].pointY[ipoint]-zeroY)*polygon2d[ipolygon].Ly+
                (polygons.polygon[ipolygon].parts[ipart].pointZ[ipoint]-zeroZ)*polygon2d[ipolygon].Lz;
            }
        }
        
    }
    
}

int voxelize (vox_polygons polygons, vox_grid_params gridparams,
              vox_rays rays, int *grids){

    projection_2d_polygon *polygon2d;
    
    int ipolygon;
    int ipart;
    int iray;

    int ix, iy, iz;
    int increment;

    double x, y, z;
    
    for (iray = 0; iray<rays.nrays; iray++){
        for (ix = 0; ix<gridparams.nx; ix++){
            for (iy = 0; iy<gridparams.ny; iy++){
                for (iz = 0; iz<gridparams.nz; iz++){
                    increment = ix + gridparams.nx*iy + gridparams.nx*gridparams.ny*iz +
                    gridparams.nx*gridparams.ny*gridparams.nz*iray;
                    grids[increment] = 0;
                }
            }
        }
    }
    
    
    
    polygon2d = malloc(sizeof(projection_2d_polygon)*polygons.npolygons);
    for (ipolygon = 0; ipolygon<polygons.npolygons; ipolygon++){
        polygon2d[ipolygon].parts = malloc(sizeof(projection_2d_part)*polygons.polygon[ipolygon].nparts);
        for (ipart=0; ipart<polygons.polygon[ipolygon].nparts; ipart++) {
            polygon2d[ipolygon].parts[ipart].pointK =
            malloc(sizeof(double)*polygons.polygon[ipolygon].parts[ipart].npoints);
            polygon2d[ipolygon].parts[ipart].pointL =
            malloc(sizeof(double)*polygons.polygon[ipolygon].parts[ipart].npoints);
        }
    }
    
    make_2d_projections(polygons, polygon2d);

    #pragma omp parallel private(x,y,z,ix,iy,iz,iray,ipolygon,increment)
    {

    #pragma omp for schedule(dynamic) collapse(3)
    for (ipolygon = 0; ipolygon<polygons.npolygons; ipolygon++)
    {
       for (iray = 0; iray<rays.nrays; iray++)
       {
          for (ix = 0; ix<gridparams.nx; ix++){
          x = gridparams.startX + ix * gridparams.dx + gridparams.dx/2.0;
          for (iy = 0; iy<gridparams.ny; iy++){
          y = gridparams.startY + iy * gridparams.dy + gridparams.dy/2.0;
          for (iz = 0; iz<gridparams.nz; iz++){
          z = gridparams.startZ + iz * gridparams.dz + gridparams.dz/2.0;
              increment = ix + gridparams.nx*iy + gridparams.nx*gridparams.ny*iz +
              gridparams.nx*gridparams.ny*gridparams.nz*iray;
            
              grids[increment] += intersection(x, y, z,
                                  rays.vX[iray], rays.vY[iray], rays.vZ[iray],
                                  &(polygons.polygon[ipolygon]),&(polygon2d[ipolygon]));
          }
          }
          }
       }
    }
    }
    return(0);
}
