//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================

// Framework include files
#include "DD4hep/Detector.h"
#include "DD4hep/Memory.h"
#include "DD4hep/DD4hepUI.h"
#include "DD4hep/Factories.h"
#include "DD4hep/Printout.h"
#include "DD4hep/DetectorTools.h"
#include "DD4hep/DetFactoryHelper.h"
#include "XML/DocumentHandler.h"
#include "XML/Utilities.h"

// ROOT includes
#include "TInterpreter.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TGDMLParse.h"
#include "TGDMLWrite.h"
#include "TUri.h"

using namespace std;
using namespace dd4hep;

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,13,0)

/// Factory to import subdetectors from GDML fragment
static Ref_t create_detector(Detector& description, xml_h e, Ref_t /* sens_det */)  {
  using namespace dd4hep::detail;
  xml_det_t   x_det = e;
  int         id    = x_det.hasAttr(_U(id)) ? x_det.id() : 0;
  xml_dim_t   x_pos  (x_det.child(_U(position),false));
  xml_dim_t   x_rot  (x_det.child(_U(rotation),false));
  xml_dim_t   x_gdml (x_det.child(_U(gdmlFile)));
  xml_dim_t   x_par  (x_det.child(_U(parent)));
  string      name    = x_det.nameStr();
  string      par_nam = x_par.nameStr();
  string      gdml    = x_gdml.attr<string>(_U(ref));
  string      gdml_name = x_gdml.attr<string>(_U(name));
  DetElement  det_parent = description.detector(par_nam);
  TGDMLParse parser;
  if ( !gdml.empty() && gdml[0] == '/' )  {
    TUri uri(gdml.c_str());
    gdml = uri.GetRelativePart();
  }
  else {
    string path = xml::DocumentHandler::system_path(e, gdml);
    TUri uri(path.c_str());
    gdml = uri.GetRelativePart();
  }
  if ( !det_parent.isValid() )  {
    except(name,"+++ Cannot access detector parent: %s",par_nam.c_str());
  }
  DetElement  sdet(name, id);
  Volume volume = parser.GDMLReadFile(gdml.c_str());
  if ( !volume.isValid() )   {
    except("ROOTGDMLParse","+++ Failed to parse GDML file:%s",gdml.c_str());
  }
  volume.import(); // We require the extensions in dd4hep.
  printout(INFO,"ROOTGDMLParse","+++ Attach GDML volume %s", volume.name());
  Volume mother = det_parent.volume();
  PlacedVolume pv;

  if ( !gdml_name.empty() ) {
    PlacedVolume node = volume->FindNode(gdml_name.c_str());
    if ( !node.isValid() ) {
      volume->PrintNodes();
      except("ROOTGDMLParse","+++ Failed to parse GDML file:%s for node:%s",
        gdml.c_str(), gdml_name.c_str());
    }
    volume = node.volume();
  }

  if ( x_pos && x_rot )   {
    Rotation3D rot(RotationZYX(x_rot.z(),x_rot.y(),x_rot.x()));
    Transform3D transform(rot,Position(x_pos.x(),x_pos.y(),x_pos.z()));
    pv = mother.placeVolume(volume,transform);
  }
  else if ( x_rot )  {
    Rotation3D rot(RotationZYX(x_rot.z(),x_rot.y(),x_rot.x()));
    Transform3D transform(rot,Position(0,0,0));
    pv = mother.placeVolume(volume,transform);
  }
  else if ( x_pos )   {
    pv = mother.placeVolume(volume,Position(x_pos.x(),x_pos.y(),x_pos.z()));
  }
  else  {
    pv = mother.placeVolume(volume);
  }
  volume.setVisAttributes(description, x_det.visStr());
  volume.setLimitSet(description, x_det.limitsStr());
  volume.setRegion(description, x_det.regionStr());
  if ( id != 0 )  {
    pv.addPhysVolID("system", id);
  }
  sdet.setPlacement(pv);
  return sdet;
}

// first argument is the type from the xml file
DECLARE_DETELEMENT(DD4hep_GdmlDetector,create_detector)

#endif
