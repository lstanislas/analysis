#include <cmath>

#include <TGeoManager.h>

#include "DataFormatsMCH/Cluster.h"
#include "MathUtils/Cartesian.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHSimulation/Response.h"

using o2::mch::Cluster;

//_________________________________________________________________________________________________
Cluster MakeCluster(uint32_t uid, float x, float y, float ex = 0.2, float ey = 0.2)
{
  /// create a cluster in the global coordinate system from a local position on a DE
  /// the DE ID is extracted from the uid, which must be defined as in Cluster::buildUniqueId(...)
  /// the geometry must be loaded (e.g. from the CCDB) prior to use this function

  static auto transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);

  auto de = Cluster::getDEId(uid);
  o2::math_utils::Point3D<float> local{x, y, 0.};
  auto global = transformation(de)(local);

  return Cluster{global.x(), global.y(), global.z(), ex, ey, uid, 0, 0};
}

//_________________________________________________________________________________________________
o2::math_utils::Point3D<float> GlobalToLocal(int de, float x, float y, float z, bool run2 = false)
{
  /// return the cluster (or track) position in the local coordinate system

  static o2::mch::geo::TransformationCreator transformation;
  if (!transformation) {
    if (run2) {
      std::ifstream geoFile("AlignedGeometry.json");
      if (!geoFile.is_open()) {
        std::cout << "cannot open geometry file AlignedGeometry.json" << std::endl;
        exit(-1);
      }
      transformation = o2::mch::geo::transformationFromJSON(geoFile);
    } else {
      transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
    }
  }

  o2::math_utils::Point3D<float> global{x, y, z};
  return transformation(de) ^ global;
}
//_________________________________________________________________________________________________
o2::math_utils::Point3D<float> LocalToGlobal(int de, float x, float y, float z, bool run2 = false)
{
  /// return the cluster (or track) position in the global coordinate system

  static o2::mch::geo::TransformationCreator transformation;
  if (!transformation) {
    if (run2) {
      std::ifstream geoFile("AlignedGeometry.json");
      if (!geoFile.is_open()) {
        std::cout << "cannot open geometry file AlignedGeometry.json" << std::endl;
        exit(-1);
      }
      transformation = o2::mch::geo::transformationFromJSON(geoFile);
    } else {
      transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
    }
  }

  o2::math_utils::Point3D<float> local{x, y, z};
  return transformation(de)(local);
}
//_________________________________________________________________________________________________
float DistanceToClosestWire(int de, float x)
{
  /// return the distance (x) to the closest wire in the local coordinate system
  /// from the cluster position in the local coordinate system

  static const o2::mch::Response response[] = {{o2::mch::Station::Type1}, {o2::mch::Station::Type2345}};

  int iSt = (de < 300) ? 0 : 1;
  float wire = response[iSt].getAnod(x);

  return x - wire;
};

//_________________________________________________________________________________________________
float DistanceToClosestWire(int de, float x, float y, float z, bool run2 = false)
{
  /// return the distance (x) to the closest wire in the local coordinate system
  /// from the cluster (or track) position in the global coordinate system

  auto local = GlobalToLocal(de, x, y, z, run2);

  return DistanceToClosestWire(de, local.x());
};
