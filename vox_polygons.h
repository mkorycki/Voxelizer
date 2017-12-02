//
//  vox_polygons.h
//  Voxelizer
//
//  Created by Michał Korycki.
//  Copyright © 2017 Michał Korycki. All rights reserved.
//

#ifndef VOX_POLYGONS
#define VOX_POLYGONS

typedef struct{
    int npoints;
    short isfull;
    double *pointX;
    double *pointY;
    double *pointZ;
}vox_polygon_part;

typedef struct{
    int nparts;
    vox_polygon_part *parts;
    double A;
    double B;
    double C;
    double D;
    double MinX, MaxX;
    double MinY, MaxY;
    double MinZ, MaxZ;
}vox_polygon;

typedef struct{
    int npolygons;
    vox_polygon *polygon;
}vox_polygons;

#endif
