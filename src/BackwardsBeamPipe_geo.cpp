//==========================================================================
//
//      <detector id="Pipe_in" name ="DetName" type="BackwardsBeamPipe" >
//      <layer id="#(int)" inner_r="#(double)" outer_z="#(double)" >
//      <slice material="string" thickness="#(double)" >
//      </layer>
//      </detector>
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */)
{

  using namespace ROOT::Math;
  xml_det_t  x_det    = e;
  string     det_name = x_det.nameStr();
  xml_dim_t  rot      = x_det.rotation();
  DetElement sdet(det_name, x_det.id());
  Assembly   assembly(det_name + "_assembly");
  Material   m_Al     = det.material("Aluminum");
  Material   m_Vacuum = det.material("Vacuum");
  string     vis_name = dd4hep::getAttrOrDefault(x_det, _Unicode(vis), "GrayVis");

  xml::Component IP_pipe_c = x_det.child(_Unicode(Pipe));

  // IP
  double thickness = IP_pipe_c.attr<double>(_Unicode(wall_thickness));
  double innerD1   = IP_pipe_c.hasAttr(_Unicode(innerD1)) ? IP_pipe_c.attr<double>(_Unicode(innerD1))
                                                          : IP_pipe_c.attr<double>(_Unicode(outerD1)) - 2 * thickness;
  double innerD2   = IP_pipe_c.hasAttr(_Unicode(innerD2)) ? IP_pipe_c.attr<double>(_Unicode(innerD2))
                                                          : IP_pipe_c.attr<double>(_Unicode(outerD2)) - 2 * thickness;
  double xoff      = dd4hep::getAttrOrDefault(IP_pipe_c, _Unicode(xoff), 0.0 );
  double end1      = IP_pipe_c.attr<double>(_Unicode(end1));
  double end2      = IP_pipe_c.attr<double>(_Unicode(end2));

  double length = abs(end2 - end1);

  // -----------------------------
  // IP beampipe
  ConeSegment tube_vacuum(length / 2.0, 0.0, innerD1 / 2.0, 0.0, innerD2 / 2.0);
  ConeSegment tube_tube(length / 2.0, innerD1 / 2.0, innerD1 / 2.0 + thickness, innerD2 / 2.0,
                        innerD2 / 2.0 + thickness);

  Volume v_vacuum("v_vacuum", tube_vacuum, m_Vacuum);
  Volume v_tube("v_tube", tube_tube, m_Al);

  sdet.setAttributes(det, v_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

  assembly.placeVolume(v_vacuum, Position(0, 0, -length / 2.0));
  assembly.placeVolume(v_tube, Position(0, 0, -length / 2.0));

  // -----------------------------
  // final placement
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly, Transform3D(RotationZYX(rot.x(),rot.y(),rot.z()),Position(xoff, 0.0, end1)));
  pv_assembly.addPhysVolID("system", sdet.id()).addPhysVolID("barrel", 1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox();
  return sdet;
}

DECLARE_DETELEMENT(BackwardsBeamPipe, create_detector)
