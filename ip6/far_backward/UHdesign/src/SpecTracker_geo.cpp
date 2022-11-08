#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  sens.setType("tracker");

  xml_det_t     x_det           = e;
  xml_comp_t    x_dim           = x_det.dimensions();
  xml_comp_t    x_pos           = x_det.position();
  xml_comp_t    x_rot           = x_det.rotation();
  //
  string        det_name        = x_det.nameStr();
  int           det_ID          = x_det.id();
  //
  double        sizeX           = x_dim.x();
  double        sizeY           = x_dim.y();
  double        sizeZ           = x_dim.z();
  double        posX            = x_pos.x();
  double        posY            = x_pos.y();
  double        posZ            = x_pos.z();
  double        rotX            = x_rot.x();
  double        rotY            = x_rot.y();
  double        rotZ            = x_rot.z();

  // Create main detector to be returned 
  DetElement det( det_name, det_ID );
  
  Volume motherVol = description.pickMotherVolume( det );
  
  // Option to pick a specific mother volume.  
  // However, once detector is placed within, mother volume is invisible...
  // std::string   motherName = dd4hep::getAttrOrDefault<std::string>( x_det, _Unicode(place_into), "" );
  // VolumeManager man        = VolumeManager::getVolumeManager( description );
  // DetElement    mdet       = man.detector().child( motherName );
  // Volume        motherVol  = mdet.volume();

  // Build detector
  Box box( sizeX, sizeY, sizeZ );
  Volume vol( det_name + "_vol", box, description.material( "Silicon" ) );
  vol.setVisAttributes( description.visAttributes( x_det.visStr() ) );
  vol.setSensitiveDetector( sens );

  // Transform
  Transform3D tr( Translation3D(posX, posY, posZ) * RotationZYX(rotZ, rotY, rotX) );

  // Place detector
  PlacedVolume detPV = motherVol.placeVolume( vol, tr );
  
  // Connect readout 
  detPV.addPhysVolID( "system", det_ID );

  det.setPlacement( detPV );

  return det;
}

DECLARE_DETELEMENT(LumiSpecTracker, create_detector)
