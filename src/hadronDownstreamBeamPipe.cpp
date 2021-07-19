//==========================================================================
//
//      <detector name ="DetName" type="Beampipe" >
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

/** \addtogroup beamline Beamline Instrumentation 
 */

/** \addtogroup IRChamber Interaction Region Vacuum Chamber.
 * \brief Type: **IRChamber**.
 * \ingroup beamline
 *
 *
 * \code
 *   <detector>
 *   </detector>
 * \endcode
 *
 */
static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector sens)  {

  using namespace ROOT::Math;
  xml_det_t  x_det     = e;
  string     det_name  = x_det.nameStr();
  Material   air       = det.air();
  DetElement sdet        (det_name,x_det.id());
  Assembly   assembly    (det_name+"_assembly");
  Material   m_Cu    = det.material("Copper");
  Material   m_Al    = det.material("Aluminum");
  Material   m_Be    = det.material("Beryllium");
  Material   m_SS    = det.material("StainlessSteel");
  string     vis_name  = x_det.visStr();

  xml::Component pos   = x_det.position();
  xml::Component rot   = x_det.rotation();
    
  /// hard-code defintion here, then refine and make more general
    
  double b0_hadron_tube_inner_r = 30.0; // mm
  double b0_hadron_tube_outer_r = 32.0; //mm
  double b0_hadron_tube_length  = 1200.0; //mm
    
    
  Tube b0_hadron_tube(b0_hadron_tube_inner_r, b0_hadron_tube_outer_r, b0_hadron_tube_length/2.0);
  
  Volume v_b0_hadron_tube("b0_hadron_tube", b0_hadron_tube, m_Be);

  //v_upstream_IP_tube.setVisAttributes(det,"GrayVis");
  //v_downstream_IP_tube.setVisAttributes(det,"RedVis");
  sdet.setAttributes(det, v_b0_hadron_tube  , x_det.regionStr(), x_det.limitsStr(), vis_name);
  //sdet.setAttributes(det, v_downstream_IP_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

  auto pv_b0_hadron_tube = assembly.placeVolume( b0_hadron_tube, Position(pos.x(), pos.y(), pos.z()));

 
  //Cone upstream_conic_section(upstream_conic_length / 2.0,
  //                            IP_beampipe_ID / 2.0, IP_beampipe_OD / 2.0,
  //                            IP_beampipe_ID / 2.0 + upstream_delta_r,
  //                            IP_beampipe_OD / 2.0+ upstream_delta_r  + IP_beampipe_thickness);
  
  
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly);
  pv_assembly.addPhysVolID("system",sdet.id()).addPhysVolID("barrel",1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox() ;
  return sdet;
}

DECLARE_DETELEMENT(handronDownstreamBeamPipe,create_detector)
