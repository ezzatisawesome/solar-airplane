// Client-side port of the mass-balance math (process_mass.py).
//
// Each exported component MassProperties stores mass, its CG position (x_cg,y_cg,z_cg) and its
// inertia (Ixx..Ixz) ABOUT ITS OWN CG. Dragging a point mass is a rigid translation: the
// component's own-CG inertia is unchanged, only its CG moves. The total is then re-summed with
// the parallel-axis theorem using AeroSandbox's exact convention
// (get_inertia_tensor_about_point): R = point - cg, J = I + m(R·R·δ - R⊗R), products J_ab = I_ab - m·R_a·R_b.

import type { ComponentMP, MassProperties, PositionInfo } from "./types";

export function combineTotal(components: Record<string, ComponentMP>): ComponentMP {
  let m = 0, mx = 0, my = 0, mz = 0;
  for (const c of Object.values(components)) {
    m += c.mass;
    mx += c.mass * c.x_cg;
    my += c.mass * c.y_cg;
    mz += c.mass * c.z_cg;
  }
  const x = m ? mx / m : 0;
  const y = m ? my / m : 0;
  const z = m ? mz / m : 0;

  let Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Iyz = 0, Ixz = 0;
  for (const c of Object.values(components)) {
    const Rx = x - c.x_cg;
    const Ry = y - c.y_cg;
    const Rz = z - c.z_cg;
    const RR = Rx * Rx + Ry * Ry + Rz * Rz;
    Ixx += c.Ixx + c.mass * (RR - Rx * Rx);
    Iyy += c.Iyy + c.mass * (RR - Ry * Ry);
    Izz += c.Izz + c.mass * (RR - Rz * Rz);
    Ixy += c.Ixy - c.mass * Rx * Ry;
    Iyz += c.Iyz - c.mass * Ry * Rz;
    Ixz += c.Ixz - c.mass * Rz * Rx;
  }
  return { mass: m, x_cg: x, y_cg: y, z_cg: z, Ixx, Iyy, Izz, Ixy, Iyz, Ixz };
}

export type Layout = Record<string, [number, number, number]>;

/** Recompute mass properties with a set of dragged point-mass positions applied. */
export function applyLayout(mp: MassProperties, layout: Layout): MassProperties {
  const components: Record<string, ComponentMP> = {};
  for (const [name, c] of Object.entries(mp.components)) {
    const ov = layout[name];
    const draggable = mp.positions[name]?.draggable;
    components[name] = ov && draggable ? { ...c, x_cg: ov[0], y_cg: ov[1], z_cg: ov[2] } : { ...c };
  }

  const positions: Record<string, PositionInfo> = {};
  for (const [name, p] of Object.entries(mp.positions)) {
    const ov = layout[name];
    positions[name] = ov && p.draggable ? { ...p, xyz: ov } : p;
  }

  return { ...mp, components, positions, total: combineTotal(components), applied_overrides: layout };
}
