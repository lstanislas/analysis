#include <map>

#include <TGeoManager.h>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "MCHTracking/TrackExtrap.h"

using o2::mch::TrackExtrap;

//_________________________________________________________________________________________________
void InitFromCCDB(int run, bool local, bool loadGeom, bool loadField)
{
  /// load necessary objects from CCDB and prepare track extrapolation if needed

  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();

  if (local) {

    // setup CCDB from local snapshot
    ccdb.setURL("file://ccdb");

    // in this case we cannot retrieve the timestamp of the run so we must know it
    static const std::map<int, long> timeStamps{{529691, 1669611720310},
                                                {544490, 1697060413764},
                                                {549586, 1712392003826},
                                                {559410, 1730513833932},
                                                {562849, 1746885436504},
                                                {562967, 1747291152329},
                                                {562968, 1747292471339},
                                                {562969, 1747293553416},
                                                {562970, 1747294630226},
                                                {562971, 1747295593716},
                                                {562972, 1747296663426},
                                                {562975, 1747299985774},
                                                {562976, 1747301213235},
                                                {562979, 1747304198913}};
    if (auto ts = timeStamps.find(run); ts != timeStamps.end()) {
      ccdb.setTimestamp(ts->second);
    } else {
      LOG(error) << "unknown time stamp for run " << run;
      exit(-1);
    }

  } else {

    // retrieve the timestamp of the run from the ALICE CCDB server
    auto [tStart, tEnd] = ccdb.getRunDuration(run);
    ccdb.setTimestamp(tEnd);
  }

  // load geometry if requested
  if (loadGeom) {
    auto geom = ccdb.get<TGeoManager>("GLO/Config/GeometryAligned");
  }

  // load magnetic field and prepare track extrapolation if requested
  if (loadField) {
    auto grp = ccdb.get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
    o2::base::Propagator::initFieldFromGRP(grp);
    TrackExtrap::setField();
  }
}
