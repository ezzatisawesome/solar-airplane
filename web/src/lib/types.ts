// Shapes of the static JSON exported by `python export_site.py`.

export interface ComponentMP {
  mass: number;
  x_cg: number;
  y_cg: number;
  z_cg: number;
  Ixx: number;
  Iyy: number;
  Izz: number;
  Ixy: number;
  Iyz: number;
  Ixz: number;
}

export interface PositionInfo {
  xyz: [number, number, number];
  kind: "point" | "box" | "cylinder" | "disk";
  draggable: boolean;
  category: string;
}

export interface MassProperties {
  components: Record<string, ComponentMP>;
  total: ComponentMP;
  positions: Record<string, PositionInfo>;
  applied_overrides?: Record<string, [number, number, number]>;
  total_mass_verification?: Record<string, number>;
}

export interface Build {
  name: string;
  masses: Record<string, number>;
  meta: Record<string, unknown>;
}

export interface SpecSection {
  title: string;
  rows: [string, string][];
}

export interface SpecTab {
  id: string;
  label: string;
  sections: SpecSection[];
}

export interface RunState {
  run: string;
  soln: Record<string, any>;
  layout: Record<string, [number, number, number]>;
  massProperties: MassProperties;
  massChart: any;
  energyBalance: any | null;
  xflrLoads: any | null;
  hasSensitivity?: boolean;
  specs: SpecTab[];
  builds: Build[];
  constants: { fuselage_radius: number };
  airfoils: Airfoil[];
}

export interface Airfoil {
  surface: string;
  name: string;
  coords: [number, number][];
}

export interface RunMetrics {
  total_mass: number | null;
  x_cg: number | null;
  static_margin: number | null;
  L_over_D: number | null;
  wingspan: number | null;
  airspeed: number | null;
}

export interface ManifestEntry {
  id: string;
  label: string;
  date: string | null;
  metrics: RunMetrics;
  builds: string[];
}

export interface Manifest {
  runs: ManifestEntry[];
}
