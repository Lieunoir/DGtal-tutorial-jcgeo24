#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>


using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;

int main()
{
  polyscope::init();

  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();

  auto h=0.5; //gridstep
  params( "polynomial", "goursat" )( "gridstep", h );
  // Build an implicit surface from a given formula
  auto implicit_shape  = SH3::makeImplicitShape3D  ( params );
  // Add some digitization infos (gridstep, bounding box size)
  auto digitized_shape = SH3::makeDigitizedImplicitShape3D( implicit_shape, params );
  // Keep these digitization infos
  auto K               = SH3::getKSpace( params );
  // Build the actual digitization as a "3d image"
  auto binary_image    = SH3::makeBinaryImage( digitized_shape, params );
  // Recover the boundary of the voxels as a surface graph
  auto surface         = SH3::makeLightDigitalSurface( binary_image, K, params );
  SH3::Cell2Index c2i;
  // Build a set of vertices and faces from the graph
  auto primalSurface   = SH3::makePrimalPolygonalSurface(c2i, surface);

  //Need to convert the faces
  std::vector<std::vector<std::size_t>> faces;
  for(auto &face: primalSurface->allFaces())
    faces.push_back(primalSurface->verticesAroundFace( face ));

  auto digsurf = polyscope::registerSurfaceMesh("Primal surface", primalSurface->positions(), faces);
  digsurf->setEdgeWidth(1.0);
  digsurf->setEdgeColor({1.,1.,1.});

  polyscope::show();
  return 0;
}
