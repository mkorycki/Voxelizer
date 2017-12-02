# Voxelizer
A simple but accurate tool for voxelization.

## Features
- OpenMP parallelization
- freedom of ray tracing configuration
- no graphic card required

## Basic data types
To take advantage of Voxelizer You need to get to know basic data types. It is your responsibility to correctly set the data.

### Grid Parameters
The grid for voxelization process is orthogonal and regular. To set grid parameters use a single instance of `vox_grid_params` structure.
```
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
```
The lower left corner of the grid is defined by `startX`, `startY` and `startZ` coordinates.
Spatial grid increment in three dimensions is defined by `dx`, `dy`, `dz`. Every increment must be greater than 0.
A number of grid boxes in every dimension is defined by `nx`, `ny`, `nz`. Every number must be greater than 0.

### Polygons
Polygons are the most complex structure set in Voxelizer. Fortunately, You don't have to set all values in them.
You just allocate a suitable amount of memory, set all polygons points and let the `vox_poly_prepare` function do the rest.

The most important thing to understand is how polygons are organized.  Polygon is made of one or more parts of type `vox_polygon_part`.
If a polygon is to be represented by a single set of points where every point is connected to two other points, only one part is suitable.
If a polygon is to contain single or more 'holes' in its surface, additional parts are needed. For example, if You want to create w polygon describing a wheel, You
would create it from two parts: outer circle and inner circle. Inner circle would describe a 'hole' in wheel surface. In Voxelizer all polygons describing outer
point set in clockwise order and 'holes' are organized counterclockwise or simply in opposite order to the first part.
```
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
```
Your responsibility is to set `npolygons` to a total number of polygons and allocate `polygon`. Later set `nparts` in every polygon and allocate `parts`.
Finally set `npoints` and allocate `pointX`, `pointY`, `pointZ`. You will be all set to put every point's coordinates in their's arrays remembering that
always first and last point in every part must be the same. Every other variable is calculated by `vox_poly_prepare` function.

### Rays for ray tracing
One of the biggest advantages of Voxelizer is the freedom of configuring own rays. It is userfull because sometimes your data may have view polygons missing.
In that case a ray pointed in some arbitrary chosen direction would trace an incorrect number of interseciotns causing errors. If Your data has no missing polygons
You do not have to worry about it and create any possible ray tracing direction. However if there are some missing you would like to point rays in different directions
or make more of them.
```
typedef struct{
    int nrays;
    double *vX;
    double *vY;
    double *vZ;
}vox_rays;
```
The only thing to doo is to set `nrays` to a number of rays to operate on and set them in vector form in `vX`, `vY` and `vZ` arrays. Normalized vectors are recommended.
If you do not have any missing polygons in your data You could create only one ray in any direction.

## Performing voxelization
After exporting your data to `vox_polygons`, creating `vox_rays` and `vox_grid_params` you can start performing calculations. The only thing left to do
is to fill all missing variables in polygons structures. Fortunately, this is achieved easily by:
```
void vox_poly_prepare(vox_polygons polygons);
```
After preparing polygons, allocate `grids` and  perform voxelization:
```
int *grids = malloc(sizeof(int)*
                 grid_params.nx*
                 grid_params.ny*
                 grid_params.nz*
                 rays.nrays);
```
and call this function:
```
int voxelize (vox_polygons polygons, vox_grid_params gridparams, vox_rays rays, int *grids);
```
Depending on the size of the grid, the number of polygons and rays this may take some time. `voxelize` will produce numbers of intersections created in process
of ray tracing. That's why voxels are not ready yet. Voxels are created by using function `postproc`:
```
int *voxels = malloc(sizeof(int)*
                  grid_params.nx*
                  grid_params.ny*
                  grid_params.nz);
```
and call function:
```
void vox_postproc(int nrays, vox_grid_params gridparams, int *grids, int *voxels);
```

## Demo
This repository contains a `demo.c` file. It is a simple program that creates two parallel polygons and performs voxelization. These two polygons are left this way intentionally not to form a solid.
If You use rays that could intersect these polygons then voxelization will be successful. That's why
```
rays.vX[0] = 0.0;
rays.vY[0] = 0.0;
rays.vZ[0] = 1.0;

rays.vX[1] = 0.0;
rays.vY[1] = 0.0;
rays.vZ[1] =-1.0;
```
configuration will be successful.
However, if you use for example:
```
rays.vX[0] = 0.0;
rays.vY[0] = 1.0;
rays.vZ[0] = 0.0;

rays.vX[1] = 0.0;
rays.vY[1] =-1.0;
rays.vZ[1] = 0.0;
```
will not bring any voxels.
Polygons defined in the demo are both parallel to XY plane. If a ray tracing vectors are parallel to those polygons no intersections will be found.
