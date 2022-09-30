#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
//////////////////////////////////////////////////
// Far Forward B0 Electromagnetic Calorimeter
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

static std::tuple<int, int> add_disk(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                     int id);

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();
  DetElement   det(detName, detID);
  sens.setType("calorimeter");
  
  // assembly
  Assembly detVol(detName);

  xml_dim_t pos = x_det.position();
  xml_dim_t rot = x_det.rotation();
  
  // module placement
  xml_comp_t plm = x_det.child(_Unicode(placements));
  std::map<int, int> sectorModuleNumbers;
  auto addModuleNumbers = [&sectorModuleNumbers](int sector, int nmod) {
    auto it = sectorModuleNumbers.find(sector);
    if (it != sectorModuleNumbers.end()) {
      it->second += nmod;
    } else {
      sectorModuleNumbers[sector] = nmod;
    }
  };
  
  int sector_id = 1;
  
  for (xml::Collection_t mod(plm, _Unicode(individuals)); mod; ++mod) {
    auto [sector, nmod] = add_individuals(desc, assembly, mod, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }
  
  for (xml::Collection_t disk(plm, _Unicode(disk)); disk; ++disk) {
    auto [sector, nmod] = add_disk(desc, assembly, disk, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }
 
  // position and rotation of parent volume 
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  det.setPlacement(detPV);
  return det;
}

// helper function to build module with or w/o wrapper
std::tuple<Volume, Position> build_module(Detector& desc, xml::Collection_t& plm, SensitiveDetector& sens)
{
  auto   mod = plm.child(_Unicode(module));
  auto   sx  = mod.attr<double>(_Unicode(sizex));
  auto   sy  = mod.attr<double>(_Unicode(sizey));
  auto   sz  = mod.attr<double>(_Unicode(sizez));
  Box    modShape(sx / 2., sy / 2., sz / 2.);
  auto   modMat = desc.material(mod.attr<std::string>(_Unicode(material)));
  Volume modVol("module_vol", modShape, modMat);
  modVol.setSensitiveDetector(sens);
  modVol.setVisAttributes(desc.visAttributes(mod.attr<std::string>(_Unicode(vis))));

  // no wrapper
  if (!plm.hasChild(_Unicode(wrapper))) {
    return std::make_tuple(modVol, Position{sx, sy, sz});
    // build wrapper
  } else {
    auto wrp       = plm.child(_Unicode(wrapper));
    auto thickness = wrp.attr<double>(_Unicode(thickness));
    if (thickness < 1e-12 * mm) {
      return std::make_tuple(modVol, Position{sx, sy, sz});
    }
    auto   wrpMat = desc.material(wrp.attr<std::string>(_Unicode(material)));
    Box    wrpShape((sx + thickness) / 2., (sy + thickness) / 2., sz / 2.);
    Volume wrpVol("wrapper_vol", wrpShape, wrpMat);
    wrpVol.placeVolume(modVol, Position(0., 0., 0.));
    wrpVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis))));
    return std::make_tuple(wrpVol, Position{sx + thickness, sy + thickness, sz});
  }
}

// place modules, id must be provided
static std::tuple<int, int> add_individuals(Detector& desc, Assembly& env, xml::Collection_t& plm,
                                            SensitiveDetector& sens, int sid)
{
  auto [modVol, modSize] = build_module(desc, plm, sens);
  int sector_id          = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  int nmodules           = 0;
  for (xml::Collection_t pl(plm, _Unicode(placement)); pl; ++pl) {
    Position    pos(dd4hep::getAttrOrDefault<double>(pl, _Unicode(x), 0.),
                    dd4hep::getAttrOrDefault<double>(pl, _Unicode(y), 0.),
                    dd4hep::getAttrOrDefault<double>(pl, _Unicode(z), 0.));
    Position    rot(dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotx), 0.),
                    dd4hep::getAttrOrDefault<double>(pl, _Unicode(roty), 0.),
                    dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotz), 0.));
    auto        mid   = pl.attr<int>(_Unicode(id));
    Transform3D tr    = Translation3D(pos.x(), pos.y(), pos.z()) * RotationZYX(rot.z(), rot.y(), rot.x());
    auto        modPV = env.placeVolume(modVol, tr);
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", mid);
    nmodules++;
  }

  return {sector_id, nmodules};
}

// place disk of modules
static std::tuple<int, int> add_disk(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                     int sid)
{
  auto [modVol, modSize] = build_module(desc, plm, sens);
  int    sector_id       = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  int    id_begin        = dd4hep::getAttrOrDefault<int>(plm, _Unicode(id_begin), 1);
  double rmin            = plm.attr<double>(_Unicode(rmin));
  double rmax            = plm.attr<double>(_Unicode(rmax));
  double phimin          = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimin), 0.);
  double phimax          = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimax), 2. * M_PI);

  // placement inside mother
  auto pos = get_xml_xyz(plm, _Unicode(position));
  auto rot = get_xml_xyz(plm, _Unicode(rotation));

  // optional envelope volume
  bool        has_envelope = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(envelope), false);
  Material    material     = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "Air"));
  Tube        solid(rmin, rmax, modSize.z() / 2.0, phimin, phimax);
  Volume      env_vol(std::string(env.name()) + "_envelope", solid, material);
  Transform3D tr_global = RotationZYX(rot.z(), rot.y(), rot.x()) * Translation3D(pos.x(), pos.y(), pos.z());
  if (has_envelope) {
    env.placeVolume(env_vol, tr_global);
  }

  // local placement of modules
  int  mid    = 0;
  auto points = epic::geo::fillRectangles({0., 0.}, modSize.x(), modSize.y(), rmin, rmax, phimin, phimax);
  for (auto& p : points) {
    Transform3D tr_local = RotationZYX(0.0, 0.0, 0.0) * Translation3D(p.x(), p.y(), 0.0);
    auto modPV = (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", id_begin + mid++);
  }
  return {sector_id, mid};
}

DECLARE_DETELEMENT(B0_ECAL, createDetector)
