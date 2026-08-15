// Reconstruct renderable geometry for each component from the solution + mass properties.
// Ported from the original Three.js viewer (reconstructGeometry). Dimensions come from the soln;
// centers come from each component's mass-properties CG (so dragged parts move with their CG).

import type { MassProperties } from "./types";

export type Geom =
  | { type: "box"; center: [number, number, number]; Lx: number; Ly: number; Lz: number; mass: number }
  | { type: "cylinder"; center: [number, number, number]; radius: number; length: number; axis: "x"; mass: number }
  | { type: "point"; center: [number, number, number]; radius: number; mass: number };

export function reconstructGeometry(
  soln: Record<string, any>,
  mp: MassProperties,
  fuselageRadius: number,
): Record<string, Geom> {
  const mainWing = soln["Main Wing"] || {};
  const geometry = soln.Geometry || {};
  const hstab = soln.HStab || {};
  const vstab = soln["V Stab"] || {};
  const power = soln.Power || {};
  const propulsion = soln.Propulsion || {};
  const airfoils = soln.airfoils || {};

  const b_w = mainWing.wingspan || mainWing.b_w || 0;
  const c_root = mainWing.chordlen || 0;
  const boom_len = geometry.boom_length || 0;
  const boom_y = geometry.boom_y || 0;
  const vstab_span = vstab.vstab_span || 0;
  const vstab_root_chord = vstab.vstab_root_chord || 0;
  const hstab_chord = hstab.hstab_chordlen || 0;
  const hstab_span = hstab.hstab_span || 0;

  const t_wing = (airfoils.wing_t_over_c || 0.12) * c_root;
  const t_tail = (airfoils.tail_t_over_c || 0.1) * Math.max(hstab_chord, 1e-6);
  const boom_radius = geometry.boom_radius || 0.01;
  const solar_side = power.solar_panel_side_length || 0.125;
  const fuselage_length = 0.5;
  const batt_box = [0.2, 0.1, 0.03];
  const solar_area = (power.solar_panels_n || 0) * solar_side * solar_side;
  const solar_mean_chord = solar_area / Math.max(b_w, 1e-6);
  const prop_radius = 0.5 * (propulsion.propeller_diameter || 0.4);

  const components = mp.components || {};
  const masses = Object.values(components).map((c) => c.mass || 0);
  const maxMass = Math.max(1e-9, ...(masses.length ? masses : [1e-9]));

  const out: Record<string, Geom> = {};
  for (const [name, comp] of Object.entries(components)) {
    const center: [number, number, number] = [comp.x_cg || 0, comp.y_cg || 0, comp.z_cg || 0];
    const mass = comp.mass || 0;

    if (name === "Main wing") {
      out[name] = { type: "box", center, Lx: 0.75 * c_root, Ly: b_w, Lz: Math.max(t_wing, 0.005), mass };
    } else if (name === "Horizontal stabilizer") {
      out[name] = { type: "box", center, Lx: hstab_chord, Ly: hstab_span, Lz: Math.max(t_tail, 0.003), mass };
    } else if (name.includes("Vertical stabilizer")) {
      out[name] = {
        type: "box",
        center,
        Lx: Math.max(0.5 * (vstab_root_chord + hstab_chord), 1e-6),
        Ly: Math.max(t_tail, 0.003),
        Lz: Math.max(vstab_span, 1e-6),
        mass,
      };
    } else if (name === "Boom L" || name === "Boom R") {
      out[name] = { type: "cylinder", center, radius: boom_radius, length: boom_len, axis: "x", mass };
    } else if (name === "Fuselage L" || name === "Fuselage R") {
      out[name] = { type: "cylinder", center, radius: fuselageRadius, length: fuselage_length, axis: "x", mass };
    } else if (name === "Solar cells") {
      out[name] = { type: "box", center, Lx: Math.max(solar_mean_chord, 0.05), Ly: b_w, Lz: 0.002, mass };
    } else if (name === "Batteries") {
      out[name] = { type: "box", center, Lx: batt_box[0], Ly: 2 * boom_y, Lz: batt_box[2], mass };
    } else if (name === "Prop L" || name === "Prop R") {
      out[name] = { type: "cylinder", center, radius: Math.max(prop_radius, 0.01), length: 0.001, axis: "x", mass };
    } else {
      out[name] = { type: "point", center, radius: 0.05 * (mass / maxMass) + 0.02, mass };
    }
  }
  return out;
}

export const CATEGORY_COLORS: Record<string, number> = {
  structure: 0x4a90e2,
  power: 0xe74c3c,
  avionics: 0x27ae60,
  propulsion: 0xe67e22,
};

export function colorFor(name: string, category?: string): number {
  if (category && CATEGORY_COLORS[category]) return CATEGORY_COLORS[category];
  if (name.includes("Boom") || name.includes("Fuselage")) return 0x8b7355;
  if (name.includes("Solar")) return 0xf1c40f;
  return 0x9b59b6;
}
