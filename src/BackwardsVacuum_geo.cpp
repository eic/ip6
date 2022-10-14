#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>

//////////////////////////////////////////////////
// Vacuum drift volume in far backwards region
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();

  xml::Component  dim = x_det.child(_Unicode(dimensions));
  double    WidthL    = dim.attr<double>(_Unicode(xL));
  double    WidthR    = dim.attr<double>(_Unicode(xL));

  double    Width     = (WidthL+WidthR)/2;
  double    Height    = dim.y();
  double    Thickness = dim.z();

  xml::Component pos = x_det.child(_Unicode(focus));
  xml_dim_t rot  = x_det.rotation();
  double    off  = pos.z();
  double    wall = dd4hep::getAttrOrDefault(x_det, _Unicode(wall), 1*mm);

  //Make bounding box to make IntersectionSolid with other components
  xml::Component BB = x_det.child(_Unicode(bounding));
  double    BB_MinX = BB.xmin();
  double    BB_MinY = BB.ymin();
  double    BB_MinZ = BB.zmin();
  double    BB_MaxX = BB.xmax();
  double    BB_MaxY = BB.ymax();
  double    BB_MaxZ = BB.zmax();

  double BB_X = abs(BB_MaxX-BB_MinX);
  double BB_Y = abs(BB_MaxY-BB_MinY);
  double BB_Z = abs(BB_MaxZ-BB_MinZ);

  Box Far_Backwards_Box(BB_X,BB_Y,BB_Z);

  //Central pipe box
  Box Extended_Beam_Box(Width+wall,Height+wall,Thickness);

  //Central vacuum box 
  Box Extended_Vacuum_Box(Width,Height,Thickness);

  //Entry box
  xml::Component EB = x_det.child(_Unicode(exitdim));
  double ED_X = EB.x();
  double ED_Y = EB.y();
  double ED_Z = off-EB.attr<double>(_Unicode(lumiZ));
  double Lumi_R = EB.attr<double>(_Unicode(lumiR));

  Box  Entry_Beam_Box(ED_X+wall,ED_Y+wall,ED_Z);
  Box  Entry_Vacuum_Box(ED_X,ED_Y,ED_Z-wall);
  Tube Lumi_Exit(0,Lumi_R,ED_Z);

  //Materials
  Material Vacuum  = desc.material("Vacuum");
  Material Air     = desc.material("Air");
  Material Silicon = desc.material("Silicon");
  Material Steel   = desc.material("StainlessSteel");
  
  Solid Wall_Box   = Extended_Beam_Box;
  Solid Vacuum_Box = Extended_Vacuum_Box;

  //Add entry boxes to main beamline
  Wall_Box   = UnionSolid(Wall_Box,   Entry_Beam_Box,   Transform3D(RotationY(-rot.theta())));
  Vacuum_Box = UnionSolid(Vacuum_Box, Entry_Vacuum_Box, Transform3D(RotationY(-rot.theta())));
  Vacuum_Box = UnionSolid(Vacuum_Box, Lumi_Exit,        Transform3D(RotationY(-rot.theta())));



  Assembly DetAssembly(detName + "_assembly");

  //Construct full vacuum chamber
  for (xml_coll_t mod(x_det, _Unicode(module)); mod; ++mod ) {

    int    moduleID  = dd4hep::getAttrOrDefault(mod, _Unicode(id), 0);
    double tagoff    = dd4hep::getAttrOrDefault(mod, _Unicode(offset_min), 50.0*mm    );
    double thetamin  = dd4hep::getAttrOrDefault(mod, _Unicode(theta_min),  0.030*rad)-rot.theta();
    double thetamax  = dd4hep::getAttrOrDefault(mod, _Unicode(theta_max),  0.030*rad)-rot.theta();
    bool   max_align = dd4hep::getAttrOrDefault(mod, _Unicode(max_align),  false    );
    double w         = dd4hep::getAttrOrDefault(mod, _Unicode(width),      200.0*mm   )/2;
    double h         = dd4hep::getAttrOrDefault(mod, _Unicode(height),     100.0*mm   )/2;
    double tagboxL   = dd4hep::getAttrOrDefault(mod, _Unicode(taglength),  0.0*mm   );
    
    auto tag_w = w;
    auto tag_h = h;
    auto vac_w = w+wall/2;
    auto vac_h = h;
//     auto vac_w = w-wall/2;
//     auto vac_h = h-wall;
    w+=wall;
    h+=wall;
     
    auto theta      = thetamin+0.0001*rad;
    auto thetaA     = thetamin; //HACK
    auto offsetx    = -w*(cos(theta));
    auto offsetz    = w*(sin(theta));
    auto vacoffsetx = -vac_w*(cos(theta));
    auto vacoffsetz = vac_w*(sin(theta));
    auto l          = (tagoff)/(sin(theta));

    auto tagoffsetx = vacoffsetx-(l+tagboxL/2)*sin(thetaA);
    auto tagoffsetz = vacoffsetz-(l+tagboxL/2)*cos(thetaA);

    if(max_align){
      theta  = thetamax+0.0001*rad;
      thetaA = thetamax;
      offsetx = -offsetx;
      offsetz    = -offsetz;
      vacoffsetx = (vac_w+wall)*(cos(theta));
      vacoffsetz = -(vac_w+wall)*(sin(theta));
      l      = (2*offsetx+tagoff)/sin(theta);
      tagoffsetx = -wall+vacoffsetx-(l+tagboxL/2)*sin(thetaA);
      tagoffsetz = vacoffsetz-(l+tagboxL/2)*cos(thetaA);
    }


    std::cout << thetamin << " "  << thetamax << std::endl; 
    std::cout << max_align << " " << w << " "  << h << " "  << l << " "  << theta << std::endl; 

    Box TagWallBox(w,h,l+tagboxL);
    Box TagVacBox(vac_w,vac_h,l+tagboxL);

    RotationZYX rotate(0,thetaA,0);
    Wall_Box   = UnionSolid(Wall_Box,   TagWallBox, Transform3D(rotate, Position(offsetx,0,offsetz)));
    Vacuum_Box = UnionSolid(Vacuum_Box, TagVacBox,  Transform3D(rotate, Position(vacoffsetx,0,vacoffsetz)));

    Box DetectorBox(vac_w,vac_h,tagboxL/2);
    Volume tagVol("Tracker_Box", DetectorBox, Vacuum);
    tagVol.setVisAttributes(desc.visAttributes("YellowVis"));

    Volume Tagger_Air;
    
    double airThickness = 0;
    double vacuumThickness = tagboxL;
    
    // Add window layer and air-vacuum boxes
    for (xml_coll_t lay(mod, _Unicode(windowLayer)); lay; ++lay) {
      
      string layerType      = dd4hep::getAttrOrDefault(lay, _Unicode(type), "window");
      string layerVis       = dd4hep::getAttrOrDefault(lay, _Unicode(vis), "RedVis");
      double layerZ         = dd4hep::getAttrOrDefault(lay, _Unicode(z), 0*mm);
      double layerThickness = dd4hep::getAttrOrDefault(lay, _Unicode(sensor_thickness), 1 * mm);
      string layerMaterial  = dd4hep::getAttrOrDefault(lay, _Unicode(material), "Copper");
      
      Material WindowMaterial  = desc.material(layerMaterial);
      
      airThickness    = tagboxL-layerZ;
      vacuumThickness = tagboxL-airThickness;
      
      Box Box_Air(vac_w,vac_h, airThickness/2);
      Tagger_Air = Volume("AirVolume",Box_Air,Air);
      Tagger_Air.setVisAttributes(desc.visAttributes("BlueVis"));
      tagVol.placeVolume(Tagger_Air, Position(0, 0, -tagboxL/2 + airThickness/2));
      
      Box    Window_Box(tag_w,tag_h, layerThickness / 2);
      Volume layVol("WindowVolume", Window_Box, WindowMaterial);
      layVol.setVisAttributes(desc.visAttributes(layerVis));
      
      Tagger_Air.placeVolume(layVol, Position(-wall/2, 0, airThickness/2-layerThickness/2));

      Box    Tag_Wall_Box(wall/2,tag_h, airThickness/2);
      Volume wallVol("TagWallVolume", Tag_Wall_Box, Steel);
      wallVol.setVisAttributes("WhiteVis");
      Tagger_Air.placeVolume(wallVol, Position(vac_w-wall/2, 0, 0));
      

    }
    
    // Add Hodoscope layers
    int N_layers = 0;
    for (xml_coll_t lay(mod, _Unicode(trackLayer)); lay; ++lay, ++N_layers) {
      
      int    layerID        = dd4hep::getAttrOrDefault(lay, _Unicode(id), 0);
      string layerType      = dd4hep::getAttrOrDefault(lay, _Unicode(type), "timepix");
      string layerVis       = dd4hep::getAttrOrDefault(lay, _Unicode(vis), "GreenVis");
      double layerZ         = dd4hep::getAttrOrDefault(lay, _Unicode(z), 0*mm);
      double layerThickness = dd4hep::getAttrOrDefault(lay, _Unicode(sensor_thickness), 200 * um);
      
      Volume mother = tagVol;
      double MotherThickness = tagboxL/2;
      
      if(layerZ>vacuumThickness){
	mother = Tagger_Air;
	layerZ -= vacuumThickness;
	MotherThickness = airThickness/2;
      }
      
      Box    Layer_Box(tag_w,tag_h, layerThickness / 2);
      Volume layVol("TrackerVolume", Layer_Box, Silicon);
      layVol.setSensitiveDetector(sens);
      layVol.setVisAttributes(desc.visAttributes(layerVis));
      
      // string module_name = detName + _toString(N_layers,"_TrackLayer_%d");
      
      PlacedVolume pv_mod = mother.placeVolume(layVol, Position(-wall/2, 0, MotherThickness - layerZ + layerThickness/2));
      pv_mod.addPhysVolID("layer", layerID);
    }

    // Add Calorimeter layers
    for (xml_coll_t lay(mod, _Unicode(calorimeter)); lay; ++lay, ++N_layers) {
      
      int    layerID        = dd4hep::getAttrOrDefault(lay, _Unicode(id),   0);
      string layerType      = dd4hep::getAttrOrDefault(lay, _Unicode(type), "TaggerCalPbWO4");
      string layerVis       = dd4hep::getAttrOrDefault(lay, _Unicode(vis),  "RedVis");
      double layerThickness = dd4hep::getAttrOrDefault(lay, _Unicode(thickness), 180 * mm);
      
      Box    Layer_Box(tag_w,tag_h, layerThickness / 2);
      Volume layVol("CalorimeterVolume", Layer_Box, desc.material("PbWO4"));
      layVol.setSensitiveDetector(sens);
      layVol.setVisAttributes(desc.visAttributes(layerVis));
      
      // string module_name = detName + _toString(N_layers,"_TrackLayer_%d");
      
      PlacedVolume pv_mod = Tagger_Air.placeVolume(layVol, Position(-wall/2, 0, -airThickness/2 + layerThickness/2));
      pv_mod.addPhysVolID("layer", layerID+100);
    }

    
    PlacedVolume pv_mod2 = DetAssembly.placeVolume(tagVol, Transform3D(rotate, Position(tagoffsetx,-wall*(((moduleID-1)*2)-1),tagoffsetz)));//Very strange y offset needs correcting for...
    pv_mod2.addPhysVolID("module", moduleID);

    
 
  }

  // Cut pipe solids so they are only in the far backwards box 
  RotationY   rotate2(-rot.theta());
  //Position    position(0,0,(off-BB_Z)/cos(rot.theta()));
  Position    position(0,0,(BB_MinZ-off-BB_Z)/cos(rot.theta()));

  IntersectionSolid Wall_Box_Sub(Wall_Box,Far_Backwards_Box,Transform3D(rotate2, position));
  IntersectionSolid Vacuum_Box_Sub(Vacuum_Box,Far_Backwards_Box,Transform3D(rotate2, position));
  
  Volume vacVol("Vacuum_Box", Vacuum_Box_Sub, Vacuum);
  vacVol.setVisAttributes(desc.visAttributes("YellowVis"));
  vacVol.placeVolume(DetAssembly); 
  Volume wallVol("Tagger_Box", Wall_Box_Sub, Steel);
  wallVol.placeVolume(vacVol);

  DetElement det(x_det.nameStr(), detID);

  // placement in mother volume
  Transform3D  tr(RotationY(rot.theta()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = desc.pickMotherVolume(det).placeVolume(wallVol, tr);
  detPV.addPhysVolID("system", detID);
  det.setPlacement(detPV);
  


  //sens.setType("tracker"); //PUT ME BACK!!!!!!!!

//   // Create Global Volume
//   Box Tagger_Box(Width, Height, Thickness);
 
//   Volume Tagger_Air;

//   double airThickness = 0;
//   double vacuumThickness = Thickness*2;

//   Volume detVol("Tagger_Box", Tagger_Box, Vacuum);
//   detVol.setVisAttributes(desc.visAttributes(x_det.visStr()));

//   // Add window layer and air-vacuum boxes
//   for (xml_coll_t lay(x_det, _Unicode(windowLayer)); lay; ++lay ) {

//     string layerType      = dd4hep::getAttrOrDefault(lay, _Unicode(type), "window");
//     string layerVis       = dd4hep::getAttrOrDefault(lay, _Unicode(vis), "RedVis");
//     double layerZ         = dd4hep::getAttrOrDefault(lay, _Unicode(z), 0*mm);
//     double layerThickness = dd4hep::getAttrOrDefault(lay, _Unicode(sensor_thickness), 1 * mm);
//     string layerMaterial  = dd4hep::getAttrOrDefault(lay, _Unicode(material), "Copper");
    
//     Material WindowMaterial  = desc.material(layerMaterial);

//     airThickness    = Thickness*2-layerZ;
//     vacuumThickness = Thickness*2-airThickness;

//     Box Box_Air(Width, Height, airThickness/2);
//     Tagger_Air = Volume("AirVolume",Box_Air,Air);
//     Tagger_Air.setVisAttributes(desc.visAttributes("BlueVis"));
//     detVol.placeVolume(Tagger_Air, Position(0, 0, -Thickness + airThickness/2));

//     Box    Window_Box(Width, Height, layerThickness / 2);
//     Volume layVol("WindowVolume", Window_Box, WindowMaterial);
//     layVol.setVisAttributes(desc.visAttributes(layerVis));

//     Tagger_Air.placeVolume(layVol, Position(0, 0, airThickness/2-layerThickness/2));
//   }

//   // Add Hodoscope layers
//   int N_layers = 0;
//   for (xml_coll_t lay(x_det, _Unicode(trackLayer)); lay; ++lay, ++N_layers) {

//     string layerType      = dd4hep::getAttrOrDefault(lay, _Unicode(type), "timepix");
//     string layerVis       = dd4hep::getAttrOrDefault(lay, _Unicode(vis), "GreenVis");
//     double layerZ         = dd4hep::getAttrOrDefault(lay, _Unicode(z), 0*mm);
//     double layerThickness = dd4hep::getAttrOrDefault(lay, _Unicode(sensor_thickness), 200 * um);

//     Volume mother = detVol;
//     double MotherThickness = Thickness;

//     if(layerZ>vacuumThickness){
//       mother = Tagger_Air;
//       layerZ -= vacuumThickness;
//       MotherThickness = airThickness/2;
//     }

//     Box    Layer_Box(Width, Height, layerThickness / 2);
//     Volume layVol("LayerVolume", Layer_Box, Silicon);
//     layVol.setSensitiveDetector(sens);
//     layVol.setVisAttributes(desc.visAttributes(layerVis));

//     // string module_name = detName + _toString(N_layers,"_TrackLayer_%d");

//     PlacedVolume pv_mod = mother.placeVolume(layVol, Position(0, 0, MotherThickness - layerZ + layerThickness/2));
//     pv_mod.addPhysVolID("layer", N_layers);
//   }



//   // mother volume for the tracker
//   std::string   mother_nam = dd4hep::getAttrOrDefault(x_det, _Unicode(place_into), "");
//   VolumeManager man        = VolumeManager::getVolumeManager(desc);
//   DetElement    mdet       = man.detector().child(mother_nam);

//   // placement in mother volume
//   Transform3D  tr(RotationZYX(0, 0, 0), Position(pos.x(), pos.y(), pos.z()));
//   PlacedVolume detPV = mdet.volume().placeVolume(detVol, tr);
//   detPV.addPhysVolID("system", detID);
//   DetElement det(detName, detID);
//   det.setPlacement(detPV);

  return det;
}

DECLARE_DETELEMENT(BackwardsVacuum, createDetector)
