// Sun geometry + energy-march helpers derived from the solution.
//
// Everything here is computed client-side from data the exporter already emits
// (Mission lat/day, Power.battery_states + capacity, solar-array geometry). It powers the
// Power tab's SOC march, the night-flight metrics, and the Overview mission summary.

const DEG = Math.PI / 180;

/** Day-of-year (1-365, non-leap) -> "Apr 10" style label. */
export function dayOfYearLabel(n: number): string {
  const days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
  const months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
  let d = Math.max(1, Math.round(n));
  for (let m = 0; m < 12; m++) {
    if (d <= days[m]) return `${months[m]} ${d}`;
    d -= days[m];
  }
  return `Day ${Math.round(n)}`;
}

export interface SunGeometry {
  dayOfYear: number;
  dateLabel: string;
  latDeg: number;
  declDeg: number; // solar declination
  daylightHours: number;
  nightHours: number;
  sunriseHod: number; // solar-time hour-of-day
  sunsetHod: number;
  noonElevationDeg: number; // sun elevation at solar noon
  polarDay: boolean; // sun never sets
  polarNight: boolean; // sun never rises
}

/** Standard sunrise-equation solar geometry for a given latitude + day-of-year. */
export function sunGeometry(latDeg: number, dayOfYear: number): SunGeometry {
  const decl = 23.45 * Math.sin(2 * Math.PI * (284 + dayOfYear) / 365); // degrees
  const cosH0 = -Math.tan(latDeg * DEG) * Math.tan(decl * DEG);
  const polarDay = cosH0 < -1;
  const polarNight = cosH0 > 1;
  const H0 = Math.acos(Math.min(1, Math.max(-1, cosH0))); // radians
  const daylight = polarDay ? 24 : polarNight ? 0 : (2 * H0) / DEG / 15;
  return {
    dayOfYear,
    dateLabel: dayOfYearLabel(dayOfYear),
    latDeg,
    declDeg: decl,
    daylightHours: daylight,
    nightHours: 24 - daylight,
    sunriseHod: 12 - daylight / 2,
    sunsetHod: 12 + daylight / 2,
    noonElevationDeg: 90 - Math.abs(latDeg - decl),
    polarDay,
    polarNight,
  };
}

export interface EnergyMarch {
  n: number;
  hours: number[]; // mission time axis, 0..24 h
  wh: number[]; // battery energy (Wh)
  soc: number[]; // 0..1 state of charge
  capacityWh: number;
  socMin: number;
  socMax: number;
  socMinHour: number; // mission time of the SOC minimum ("dawn")
  depthOfDischarge: number; // socMax - socMin
  nightEnergyWh: number; // energy drawn between the SOC peak and the following minimum
  floor: number; // 20% reserve line
  hasData: boolean;
}

const FLOOR = 0.2;

/** Battery state-of-charge march reconstructed from Power.battery_states / battery_capacity. */
export function energyMarch(soln: Record<string, any>): EnergyMarch {
  const P = soln.Power || {};
  const wh: number[] = P.battery_states ?? [];
  const capacity: number = P.battery_capacity ?? Math.max(1e-9, ...(wh.length ? wh : [1]));
  const n = wh.length;
  const hours = wh.map((_, i) => (n > 1 ? (i / (n - 1)) * 24 : 0));
  const soc = wh.map((v) => (capacity > 0 ? v / capacity : 0));

  let iMin = 0;
  let iMax = 0;
  soc.forEach((s, i) => {
    if (s < soc[iMin]) iMin = i;
    if (s > soc[iMax]) iMax = i;
  });
  const socMin = soc.length ? soc[iMin] : 0;
  const socMax = soc.length ? soc[iMax] : 0;
  // Night draw: from the peak charge down to the following minimum (fall back to full swing).
  const peakBeforeMin = soc.slice(0, iMin + 1).reduce((m, s) => Math.max(m, s), soc[0] ?? 0);
  const nightEnergyWh = Math.max(0, (peakBeforeMin - socMin) * capacity);

  return {
    n,
    hours,
    wh,
    soc,
    capacityWh: capacity,
    socMin,
    socMax,
    socMinHour: hours[iMin] ?? 0,
    depthOfDischarge: Math.max(0, socMax - socMin),
    nightEnergyWh,
    floor: FLOOR,
    hasData: n > 1,
  };
}

export interface SolarArray {
  nCells: number;
  cellSide: number; // m
  areaM2: number;
  onWing: number;
  onHstab: number;
  cellEff: number; // STC cell efficiency
  systemEff: number; // STC system efficiency
  peakPowerW: number; // area * 1000 W/m^2 * system efficiency
  nMppt: number;
  nPacks: number;
}

/** Solar-array summary (counts, area, peak power) from the Power block. */
export function solarArray(soln: Record<string, any>): SolarArray {
  const P = soln.Power || {};
  const side = P.solar_panel_side_length ?? 0.125;
  const nCells = P.solar_panel_n ?? 0;
  const area = nCells * side * side;
  const systemEff = P.solar_system_eff_stc ?? 0;
  return {
    nCells,
    cellSide: side,
    areaM2: area,
    onWing: P.panels_on_wing ?? nCells,
    onHstab: P.panels_on_hstab ?? 0,
    cellEff: P.cell_efficiency_stc ?? 0,
    systemEff,
    peakPowerW: area * 1000 * systemEff, // 1 kW/m^2 STC irradiance
    nMppt: P.n_mppt ?? 0,
    nPacks: P.num_packs ?? 0,
  };
}

/** Load proxy: total cruise electrical power (all motors) — the night-flight bus draw. */
export function cruiseLoadW(soln: Record<string, any>): number {
  const perf = soln.Performance || {};
  return perf["power_cruise (all motors)"] ?? 0;
}
