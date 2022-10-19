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

static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector /*sens*/)
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

  //Materials
  Material Vacuum  = desc.material("Vacuum");
  Material Steel   = desc.material("StainlessSteel");

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

  //Entry box description
  xml::Component EB = x_det.child(_Unicode(exitdim));
  double ED_X = EB.x();
  double ED_Y = EB.y();
  double ED_Z = off-EB.attr<double>(_Unicode(lumiZ));
  double Lumi_R = EB.attr<double>(_Unicode(lumiR));
  double exitTheta = EB.attr<double>(_Unicode(maxTheta));

  // Generic box for making intersection solid with
  double xbox = 10*m;
  double ybox = 10*m;
  double zbox = 50*m;
  Box Cut_Box(xbox,ybox,zbox);

  //Central pipe box
  Box Extended_Beam_Box(Width+wall,Height+wall,Thickness);

  //Central vacuum box 
  Box Extended_Vacuum_Box(Width,Height,Thickness);

  Solid Wall_Box   = Extended_Beam_Box;
  Solid Vacuum_Box = Extended_Vacuum_Box;

  Assembly DetAssembly(detName + "_assembly");
  Assembly DetAssemblyAir(detName + "_assembly_air");

  DetElement det(detName, detID);

  //-----------------------------------------------------------------
  //Add Tagger box containers and vacuum box extension for modules
  //-----------------------------------------------------------------
  for (xml_coll_t mod(x_det, _Unicode(module)); mod; ++mod ) {

    int    moduleID  = dd4hep::getAttrOrDefault(mod, _Unicode(id), 0);
    string moduleName= dd4hep::getAttrOrDefault(mod, _Unicode(modname), "Tagger0");
    double tagoff    = dd4hep::getAttrOrDefault(mod, _Unicode(offset_min), 50.0*mm    );
    double thetamin  = dd4hep::getAttrOrDefault(mod, _Unicode(theta_min),  0.030*rad)-rot.theta();
    double thetamax  = dd4hep::getAttrOrDefault(mod, _Unicode(theta_max),  0.030*rad)-rot.theta();
    bool   max_align = dd4hep::getAttrOrDefault(mod, _Unicode(max_align),  false    );
    bool   extend_vacuum = dd4hep::getAttrOrDefault(mod, _Unicode(extend_vacuum),  true    );

    xml_dim_t moddim = mod.child(_Unicode(dimensions));
    double w         = moddim.x();
    double h         = moddim.y();
    double tagboxL   = moddim.z();
    
//     auto tag_w = w;
//     auto tag_h = h;
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

    Volume   mother    = DetAssemblyAir;
    
    if(extend_vacuum){
      Wall_Box   = UnionSolid(Wall_Box,   TagWallBox, Transform3D(rotate, Position(offsetx,0,offsetz)));
      Vacuum_Box = UnionSolid(Vacuum_Box, TagVacBox,  Transform3D(rotate, Position(vacoffsetx,0,vacoffsetz)));
      mother     = DetAssembly;
    }
    
    Assembly TaggerAssembly("tagAssembly");
    
    PlacedVolume pv_mod2 = mother.placeVolume(TaggerAssembly, Transform3D(rotate, Position(tagoffsetx,0,tagoffsetz)));//Very strange y offset needs correcting for...
    pv_mod2.addPhysVolID("module", moduleID);
    DetElement moddet(moduleName, moduleID);
    moddet.setPlacement(pv_mod2);
    det.add(moddet);
    
  }

  Wall_Box   = IntersectionSolid(Wall_Box,   Cut_Box,   Position(-xbox+Width+wall,0,0));
  Vacuum_Box = IntersectionSolid(Vacuum_Box, Cut_Box,   Position(-xbox+Width,0,0));

  // Luminosity connecting box 
  bool addLumi = dd4hep::getAttrOrDefault(x_det, _Unicode(lumi),  true    );

  if(addLumi){
    
    Box  Entry_Beam_Box(ED_X+wall,ED_Y+wall,ED_Z);
    Box  Entry_Vacuum_Box(ED_X,ED_Y,ED_Z-wall);
    Tube Lumi_Exit(0,Lumi_R,ED_Z);
 
    //Add entry boxes to main beamline volume
    Wall_Box   = UnionSolid(Wall_Box,   Entry_Beam_Box,   Transform3D(RotationY(-rot.theta())));
    Vacuum_Box = UnionSolid(Vacuum_Box, Entry_Vacuum_Box, Transform3D(RotationY(-rot.theta())));
    Vacuum_Box = UnionSolid(Vacuum_Box, Lumi_Exit,        Transform3D(RotationY(-rot.theta())));

  }
  
  double exitDist = BB_MinZ-off;
  // Restrict tagger boxes into region right, connected to beampipe magnet exit
  double cutX = (ED_X-exitDist*tan(-rot.theta()))*cos(rot.theta());
  double cutZ = (ED_X-exitDist*tan(-rot.theta()))*sin(rot.theta())+exitDist*cos(rot.theta());
    

  Wall_Box   = IntersectionSolid(Wall_Box,   Cut_Box,   Transform3D(RotationY(exitTheta),Position(xbox-cutX,0,cutZ)));
  Vacuum_Box = IntersectionSolid(Vacuum_Box, Cut_Box,   Transform3D(RotationY(exitTheta),Position(xbox-cutX,0,cutZ)));

  // Cut pipe solids so they are only in the far backwards box 
  RotationY   rotate2(-rot.theta());
  //Position    position(0,0,(off-BB_Z)/cos(rot.theta()));
  Position    position(0,0,(exitDist-BB_Z)/cos(rot.theta()));

  IntersectionSolid Wall_Box_Sub  (Wall_Box,   Far_Backwards_Box, Transform3D(rotate2, position));
  IntersectionSolid Vacuum_Box_Sub(Vacuum_Box, Far_Backwards_Box, Transform3D(rotate2, position));
  

  Volume vacVol("Vacuum_Box", Vacuum_Box_Sub, Vacuum);
  vacVol.setVisAttributes(desc.visAttributes("YellowVis"));
  vacVol.placeVolume(DetAssembly); 
  Volume wallVol("Tagger_Box", Wall_Box_Sub, Steel);
  wallVol.placeVolume(vacVol);

  Assembly backAssembly("assembly");
  backAssembly.placeVolume(wallVol);
  backAssembly.placeVolume(DetAssemblyAir);


  // placement in mother volume
  Transform3D  tr(RotationY(rot.theta()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = desc.pickMotherVolume(det).placeVolume(backAssembly, tr);
  //PlacedVolume detPV = desc.pickMotherVolume(det).placeVolume(wallVol, tr);
  detPV.addPhysVolID("system", detID);
  det.setPlacement(detPV);
 
  for (xml_coll_t mod(x_det, _Unicode(module)); mod; ++mod ) {

    string moduleName= dd4hep::getAttrOrDefault(mod, _Unicode(modname), "Tagger0");
    desc.declareParent(moduleName, det);
  }

  return det;
}

DECLARE_DETELEMENT(BackwardsVacuum, create_detector)
