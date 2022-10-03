#include "Math/Point2D.h"
#include <vector>
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>

//////////////////////////////////////////////////
// Far Forward B0 Electromagnetic Calorimeter
//////////////////////////////////////////////////

using std::string;
using std::tuple;
using std::make_tuple;
using std::vector;
using std::map;
using namespace dd4hep;

static tuple<int, int> add_individuals(Detector& desc, Assembly& env, xml::Collection_t& plm,
                                            SensitiveDetector& sens, int id);
static tuple<int, int> add_disk(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                     int id);
using Point = ROOT::Math::XYPoint;
vector<Point> fillRectangles(Point ref, double sx, double sy, double rmin, double rintermediate, double rmax, double phmin = -M_PI,
                                    double phmax = M_PI);

// helper function to get x, y, z if defined in a xml component
template <class XmlComp>
Position get_xml_xyz(XmlComp& comp, dd4hep::xml::Strng_t name)
{
  Position pos(0., 0., 0.);
  if (comp.hasChild(name)) {
    auto child = comp.child(name);
    pos.SetX(dd4hep::getAttrOrDefault<double>(child, _Unicode(x), 0.));
    pos.SetY(dd4hep::getAttrOrDefault<double>(child, _Unicode(y), 0.));
    pos.SetZ(dd4hep::getAttrOrDefault<double>(child, _Unicode(z), 0.));
  }
  return pos;
}

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
  map<int, int> sectorModuleNumbers;
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
    auto [sector, nmod] = add_individuals(desc, detVol, mod, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }
  
  for (xml::Collection_t disk(plm, _Unicode(disk)); disk; ++disk) {
    auto [sector, nmod] = add_disk(desc, detVol, disk, sens, sector_id++);
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
tuple<Volume, Position> build_module(Detector& desc, xml::Collection_t& plm, SensitiveDetector& sens)
{
  auto   mod = plm.child(_Unicode(module));
  auto   sx  = mod.attr<double>(_Unicode(sizex));
  auto   sy  = mod.attr<double>(_Unicode(sizey));
  auto   sz  = mod.attr<double>(_Unicode(sizez));
  Box    modShape(sx / 2., sy / 2., sz / 2.);
  auto   modMat = desc.material(mod.attr<string>(_Unicode(material)));
  Volume modVol("module_vol", modShape, modMat);
  modVol.setSensitiveDetector(sens);
  modVol.setVisAttributes(desc.visAttributes(mod.attr<string>(_Unicode(vis))));

  // no wrapper
  if (!plm.hasChild(_Unicode(wrapper))) {
    return make_tuple(modVol, Position{sx, sy, sz});
    // build wrapper
  } else {
    auto wrp       = plm.child(_Unicode(wrapper));
    auto thickness = wrp.attr<double>(_Unicode(thickness));
    if (thickness < 1e-12 * mm) {
      return make_tuple(modVol, Position{sx, sy, sz});
    }
    auto   wrpMat = desc.material(wrp.attr<string>(_Unicode(material)));
    Box    wrpShape((sx + thickness) / 2., (sy + thickness) / 2., sz / 2.);
    Volume wrpVol("wrapper_vol", wrpShape, wrpMat);
    wrpVol.placeVolume(modVol, Position(0., 0., 0.));
    wrpVol.setVisAttributes(desc.visAttributes(wrp.attr<string>(_Unicode(vis))));
    return make_tuple(wrpVol, Position{sx + thickness, sy + thickness, sz});
  }
}

// place modules, id must be provided
static tuple<int, int> add_individuals(Detector& desc, Assembly& env, xml::Collection_t& plm,
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
static tuple<int, int> add_disk(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                     int sid)
{
  auto [modVol, modSize] = build_module(desc, plm, sens);
  int    sector_id       = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  int    id_begin        = dd4hep::getAttrOrDefault<int>(plm, _Unicode(id_begin), 1);
  double rmin            = plm.attr<double>(_Unicode(rmin));
  double rintermediate   = plm.attr<double>(_Unicode(rintermediate));
  double envelopeclearance = dd4hep::getAttrOrDefault<int>(plm, _Unicode(envelopeclearance), 0);
  double rmax            = plm.attr<double>(_Unicode(rmax));
  double phimin          = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimin), 0.);
  double phimax          = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimax), 2. * M_PI);

  // placement inside mother
  auto pos = get_xml_xyz(plm, _Unicode(position));
  auto rot = get_xml_xyz(plm, _Unicode(rotation));

  // optional envelope volume
  bool        has_envelope = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(envelope), false);
  Material    material     = desc.material(getAttrOrDefault<string>(plm, _U(material), "Air"));
  Tube        inner_solid(rmin-envelopeclearance, rintermediate+envelopeclearance, modSize.z() / 2.0, 0, 2. * M_PI);
  Tube        outer_solid(rintermediate, rmax+envelopeclearance, modSize.z() / 2.0, phimin, phimax);
  UnionSolid  solid(inner_solid, outer_solid);
  Volume      env_vol(string(env.name()) + "_envelope", solid, material);
  Transform3D tr_global = RotationZYX(rot.z(), rot.y(), rot.x()) * Translation3D(pos.x(), pos.y(), pos.z());
  if (has_envelope) {
    env.placeVolume(env_vol, tr_global);
  }

  // local placement of modules
  int  mid    = 0;
  auto points = fillRectangles({0., 0.}, modSize.x(), modSize.y(), rmin, rintermediate, rmax, phimin, phimax);
  for (auto& p : points) {
    Transform3D tr_local = RotationZYX(0.0, 0.0, 0.0) * Translation3D(p.x(), p.y(), 0.0);
    auto modPV = (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", id_begin + mid++);
  }
  return {sector_id, mid};
}

// check if a 2d point is already in the container
bool already_placed(const Point& p, const vector<Point>& vec, double xs = 1.0, double ys = 1.0,
                    double tol = 1e-6)
{
  for (auto& pt : vec) {
    if ((std::abs(pt.x() - p.x()) / xs < tol) && std::abs(pt.y() - p.y()) / ys < tol) {
      return true;
    }
  }
  return false;
}

// check if a point is in a ring
inline bool rec_in_ring(const Point& pt, double sx, double sy, double rmin, double rintermediate, double rmax, double phmin, double phmax)
{
  double rmax_pacman = (pt.phi() < phmin || pt.phi() > phmax) ? rintermediate : rmax;
  if (pt.r() > rmax_pacman || pt.r() < rmin) {
    return false;
  }

  // check four corners
  vector<Point> pts{
                           Point(pt.x() - sx / 2., pt.y() - sy / 2.),
                           Point(pt.x() - sx / 2., pt.y() + sy / 2.),
                           Point(pt.x() + sx / 2., pt.y() - sy / 2.),
                           Point(pt.x() + sx / 2., pt.y() + sy / 2.),
                        };
  bool inside = false;
  for (auto& p : pts) {
    rmax_pacman = (p.phi() < phmin || p.phi() > phmax) ? rintermediate : rmax;
    inside += (p.r() <= rmax_pacman && p.r() >= rmin);
  }
  return inside;
}

// a helper function to recursively fill square in a ring
void add_rectangle(Point p, vector<Point>& res, double sx, double sy, double rmin, double rintermediate, double rmax, double phmin,
                   double phmax, int max_depth = 20, int depth = 0)
{
  // exceeds the maximum depth in searching or already placed
  if ((depth > max_depth) || (already_placed(p, res, sx, sy))) {
    return;
  }

  bool in_ring = rec_in_ring(p, sx, sy, rmin, rintermediate, rmax, phmin, phmax);
  if (in_ring) {
    res.emplace_back(p);
  }
  // continue search for a good placement or if no placement found yet
  if (in_ring || res.empty()) {
    // check adjacent squares
    add_rectangle(Point(p.x() + sx, p.y()), res, sx, sy, rmin, rintermediate, rmax, phmin, phmax, max_depth, depth + 1);
    add_rectangle(Point(p.x() - sx, p.y()), res, sx, sy, rmin, rintermediate, rmax, phmin, phmax, max_depth, depth + 1);
    add_rectangle(Point(p.x(), p.y() + sy), res, sx, sy, rmin, rintermediate, rmax, phmin, phmax, max_depth, depth + 1);
    add_rectangle(Point(p.x(), p.y() - sy), res, sx, sy, rmin, rintermediate, rmax, phmin, phmax, max_depth, depth + 1);
  }
}

// fill squares
vector<Point> fillRectangles(Point ref, double sx, double sy, double rmin, double rintermediate, double rmax, double phmin,
                                  double phmax)
{
  // convert (0, 2pi) to (-pi, pi)
  if (phmax > M_PI) {
    phmin -= M_PI;
    phmax -= M_PI;
  }
  // start with a seed square and find one in the ring
  // move to center
  ref = ref - Point(int(ref.x() / sx) * sx, int(ref.y() / sy) * sy);
  vector<Point> res;
  add_rectangle(ref, res, sx, sy, rmin, rintermediate, rmax, phmin, phmax, (int(rmax / sx) + 1) * (int(rmax / sy) + 1) * 2);
  return res;
}

DECLARE_DETELEMENT(B0_ECAL, createDetector)
