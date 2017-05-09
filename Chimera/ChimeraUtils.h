#ifndef CHIMERAUTILS_H
#define CHIMERAUTILS_H

#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"

inline double X2Tick(double x, size_t plane){
  auto ts = larutil::TimeService::GetME();
  auto larp = larutil::LArProperties::GetME();

  // (X [cm] / Drift Velocity [cm/us] - TPC waveform tick 0 offset) ... [us]
  double tick = (x / larp->DriftVelocity() - ts->TriggerOffsetTPC());
  // 1st plane's tick
  if(plane==0) return tick * 2;// 2 ticks/us
  // 2nd plane needs more drift time for plane0=>plane1 gap (0.3cm) difference
  tick -= 0.3 / larp->DriftVelocity(larp->Efield(1));
  // 2nd plane's tick
  if(plane==1) return tick * 2;
  // 3rd plane needs more drift time for plane1=>plane2 gap (0.3cm) difference
  tick -= 0.3 / larp->DriftVelocity(larp->Efield(2));
  return tick * 2;
}

#endif
