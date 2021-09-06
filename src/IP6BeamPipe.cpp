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

#include <Acts/Plugins/DD4hep/ActsExtension.hpp>

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
static Ref_t create_beampipe_central(Detector& det, xml_h e, SensitiveDetector /* sens */)  {

  using namespace ROOT::Math;
  xml_det_t  x_det     = e;
  string     det_name  = x_det.nameStr();
  DetElement sdet        (det_name,x_det.id());
  Material   m_Vacuum  = det.material("Vacuum");
  string     vis_name  = x_det.visStr();

  // -----------------------------
  // IP beampipe:
  //
  // This needs to be placed directly, not as part of the assembly, for ACTS reasons

  xml::Component beampipe_c = x_det.child(_Unicode(beampipe));

  double upstream_straight_length   = beampipe_c.attr<double>(_Unicode(upstream_straight_length));
  double downstream_straight_length = beampipe_c.attr<double>(_Unicode(downstream_straight_length));
  double straight_dz = (upstream_straight_length + downstream_straight_length) / 2.0;
  double straight_z0 = (upstream_straight_length - downstream_straight_length) / 2.0;
  double rmax = beampipe_c.attr<double>(_Unicode(OD)) / 2.0;

  Tube     envelope(0.0, rmax, straight_dz);
  Volume v_envelope(det_name + "_envelope", envelope, m_Vacuum);

  size_t i_layer = 0;
  for (xml_coll_t x_layer_i(beampipe_c, _Unicode(layer)); x_layer_i; ++x_layer_i) {
    xml_comp_t x_layer = x_layer_i;
    double thickness = x_layer.attr<double>(_Unicode(thickness));
    Material material = det.material(x_layer.materialStr());
    Tube     layer(rmax - thickness, rmax, straight_dz);
    Volume v_layer(det_name + i_layer, layer, material);
    sdet.setAttributes(det, v_layer, x_det.regionStr(), x_det.limitsStr(), vis_name);
    v_envelope.placeVolume(v_layer);
    rmax -= thickness;
    i_layer++;
  }
  
  auto pv_envelope = det.pickMotherVolume(sdet).placeVolume(v_envelope, Position(0, 0, straight_z0));
  pv_envelope.addPhysVolID("system",sdet.id()).addPhysVolID("barrel",1);
  sdet.setPlacement(pv_envelope);

  Acts::ActsExtension* beamPipeExtension = new Acts::ActsExtension();
  beamPipeExtension->addType("beampipe", "layer");
  sdet.addExtension<Acts::ActsExtension>(beamPipeExtension);

  return sdet;
}

static Ref_t create_beampipe_assembly(Detector& det, xml_h e, SensitiveDetector /* sens */)  {

  using namespace ROOT::Math;
  xml_det_t  x_det     = e;
  string     det_name  = x_det.nameStr();
  DetElement sdet        (det_name,x_det.id());
  Material   m_Al      = det.material("Aluminum");
  Material   m_Vacuum  = det.material("Vacuum");
  string     vis_name  = x_det.visStr();

  // ---------------------------------
  // Upstream and downstream beampipes
  //
  // Helper function to create polycone pairs (shell and vacuum)
  auto zplane_to_polycones = [](xml::Component& x_pipe) {
    std::vector<double> zero, rmax, rmin, z;
    for (xml_coll_t x_zplane_i(x_pipe, _Unicode(zplane)); x_zplane_i; ++x_zplane_i) {
      xml_comp_t x_zplane = x_zplane_i;
      auto thickness = getAttrOrDefault(x_zplane, _Unicode(thickness), x_pipe.thickness());
      thickness += getAttrOrDefault(x_zplane, _Unicode(extra_thickness), 0.0);
      zero.push_back(0);
      rmax.push_back(x_zplane.attr<double>(_Unicode(OD)) / 2.0);
      rmin.push_back(x_zplane.attr<double>(_Unicode(OD)) / 2.0 - thickness);
      z.push_back(x_zplane.attr<double>(_Unicode(z)));
    }
    return std::make_pair<Polycone,Polycone>(
      {0, 2.0 * M_PI, rmin, rmax, z},
      {0, 2.0 * M_PI, zero, rmin, z}
    );
  };

  auto create_volumes = [&](const std::string& name,
                            xml::Component& x_pipe1,
                            xml::Component& x_pipe2,
                            xml_coll_t& x_additional_subtraction_i,
                            bool subtract_vacuum_from_matter = true,
                            bool subtract_matter_from_vacuum = false) {
    auto pipe1_polycones = zplane_to_polycones(x_pipe1);
    auto pipe2_polycones = zplane_to_polycones(x_pipe2);

    auto crossing_angle = getAttrOrDefault(x_pipe2, _Unicode(crossing_angle), 0.0);
    auto axis_intersection = getAttrOrDefault(x_pipe2, _Unicode(axis_intersection), 0.0);

    auto tf = Transform3D(Position(0,0,axis_intersection)) *
              Transform3D(RotationY(crossing_angle)) *
              Transform3D(Position(0,0,-axis_intersection));

    // union of all matter and vacuum
    UnionSolid matter_union(pipe1_polycones.first,
                            pipe2_polycones.first,
                            tf);
    UnionSolid vacuum_union(pipe1_polycones.second,
                            pipe2_polycones.second,
                            tf);

    // subtract vacuum from matter
    BooleanSolid matter;
    if (subtract_vacuum_from_matter) {
      matter = SubtractionSolid(matter_union, vacuum_union);
    } else {
      matter = matter_union;
    }
    // subtract matter from vacuum
    BooleanSolid vacuum;
    if (subtract_matter_from_vacuum) {
      vacuum = SubtractionSolid(vacuum_union, matter_union);
    } else {
      vacuum = vacuum_union;
    }

    // subtract additional vacuum from matter
    for (; x_additional_subtraction_i; ++x_additional_subtraction_i) {
      xml_comp_t x_additional_subtraction = x_additional_subtraction_i;
      auto additional_polycones = zplane_to_polycones(x_additional_subtraction);
      auto additional_crossing_angle = getAttrOrDefault(x_additional_subtraction, _Unicode(crossing_angle), 0.0);
      auto additional_axis_intersection = getAttrOrDefault(x_additional_subtraction, _Unicode(axis_intersection), 0.0);
      auto additional_tf = Transform3D(Position(0,0,additional_axis_intersection)) *
                           Transform3D(RotationY(additional_crossing_angle)) *
                           Transform3D(Position(0,0,-additional_axis_intersection));
      matter = SubtractionSolid(matter, additional_polycones.second, additional_tf);
    }

    return std::make_pair<Volume,Volume>(
      {"v_" + name + "_matter", matter, m_Al},
      {"v_" + name + "_vacuum", vacuum, m_Vacuum}
    );
  };

  // -----------------------------
  // Upstream:
  // - incoming hadron tube: straight section, tapered section, straight section
  // - outgoing electron tube: tapered section, straight section
  // -----------------------------
  // Downstream:
  // - incoming electron tube: tube with tube cut out
  // - outgoing hadron tube: cone centered at scattering angle
  // (incoming electron tube internal touching to outgoing hadron tube)

  Assembly assembly(det_name + "_assembly");

  xml::Component beampipe_c = x_det.child(_Unicode(beampipe));
  xml::Component pipe1_c = beampipe_c.child(_Unicode(pipe1));
  xml::Component pipe2_c = beampipe_c.child(_Unicode(pipe2));
  xml_coll_t additional_subtractions(beampipe_c, _Unicode(additional_subtraction));
  bool subtract_vacuum = getAttrOrDefault<bool>(beampipe_c, _Unicode(subtract_vacuum), true);
  bool subtract_matter = getAttrOrDefault<bool>(beampipe_c, _Unicode(subtract_matter), true);
  auto volumes = create_volumes(beampipe_c.nameStr(),
                                pipe1_c,
                                pipe2_c,
                                additional_subtractions,
                                subtract_vacuum,
                                subtract_matter
                                );

  auto tf = Transform3D(RotationZYX(0, 0, 0));
  if (getAttrOrDefault<bool>(beampipe_c, _Unicode(reflect), true)) {
    tf = Transform3D(RotationZYX(0, M_PI, 0));
  }
  assembly.placeVolume(volumes.first, tf);
  if (getAttrOrDefault<bool>(beampipe_c, _Unicode(place_vacuum), true)) {
    assembly.placeVolume(volumes.second, tf);
  }
  assembly->GetShape()->ComputeBBox() ;

  // -----------------------------
  // Final assembly placement
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly);
  pv_assembly.addPhysVolID("system",sdet.id()).addPhysVolID("barrel",1);
  sdet.setPlacement(pv_assembly);
  return sdet;
}

DECLARE_DETELEMENT(IP6BeamPipeCentral,create_beampipe_central)
DECLARE_DETELEMENT(IP6BeamPipeAssembly,create_beampipe_assembly)
