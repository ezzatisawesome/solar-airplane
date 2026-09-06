// Flat, searchable index of every label→value datum in a run, so an engineer can jump
// straight to a number without knowing which tab it lives on.

import type { RunState } from "./types";
import { parseGeometry } from "./geometry";
import { energyMarch, sunGeometry } from "./solar";

export interface SearchEntry {
  tab: string; // tab id to jump to
  tabLabel: string;
  section: string;
  label: string;
  value: string;
}

export function buildIndex(run: RunState): SearchEntry[] {
  const out: SearchEntry[] = [];

  // Every spec row across every tab.
  for (const t of run.specs ?? []) {
    for (const s of t.sections) {
      for (const [label, value] of s.rows) {
        out.push({ tab: t.id, tabLabel: t.label, section: s.title, label, value });
      }
    }
  }

  // Derived sun / night / energy figures (they live on Overview + Power but aren't spec rows).
  const lat = run.soln.Mission?.operating_lat ?? 0;
  const day = run.soln.Mission?.mission_date ?? 172;
  const sun = sunGeometry(lat, day);
  const e = energyMarch(run.soln);
  const derived: SearchEntry[] = [
    { tab: "overview", tabLabel: "Overview", section: "Sun & night", label: "Design day", value: `${sun.dateLabel} (day ${Math.round(day)})` },
    { tab: "overview", tabLabel: "Overview", section: "Sun & night", label: "Daylight hours", value: `${sun.daylightHours.toFixed(2)} h` },
    { tab: "overview", tabLabel: "Overview", section: "Sun & night", label: "Night hours", value: `${sun.nightHours.toFixed(2)} h` },
    { tab: "overview", tabLabel: "Overview", section: "Sun & night", label: "Solar noon elevation", value: `${sun.noonElevationDeg.toFixed(1)}°` },
  ];
  if (e.hasData) {
    derived.push(
      { tab: "power", tabLabel: "Power", section: "Energy budget", label: "Min SOC (dawn)", value: `${(100 * e.socMin).toFixed(1)}%` },
      { tab: "power", tabLabel: "Power", section: "Energy budget", label: "Depth of discharge", value: `${(100 * e.depthOfDischarge).toFixed(0)}%` },
      { tab: "power", tabLabel: "Power", section: "Energy budget", label: "Overnight discharge", value: `${e.nightEnergyWh.toFixed(0)} Wh` },
      { tab: "power", tabLabel: "Power", section: "Energy budget", label: "Night closes", value: e.socMin >= e.floor ? "yes" : "no" },
    );
  }
  out.push(...derived);

  // Geometry dimensions (also surfaced on Airfoils & Planform).
  const g = parseGeometry(run.soln);
  const dims: [string, number][] = [
    ["Wingspan", g.b], ["Wing chord", g.cRoot], ["Wing tip chord", g.cTip],
    ["H-stab span", g.hstabSpan], ["H-stab chord", g.hstabChord],
    ["V-stab height", g.vstabSpan], ["V-stab root chord", g.vstabRootChord],
    ["Boom length", g.boomLen], ["Boom Y", g.boomY], ["Fuselage length", g.fuseLen],
  ];
  for (const [label, v] of dims) {
    out.push({ tab: "shapes", tabLabel: "Airfoils & Planform", section: "Dimensions", label, value: `${v.toFixed(3)} m` });
  }

  return out;
}

/** Case-insensitive token match on label + value + section. */
export function searchIndex(index: SearchEntry[], query: string): SearchEntry[] {
  const q = query.trim().toLowerCase();
  if (!q) return [];
  const tokens = q.split(/\s+/);
  return index
    .filter((e) => {
      const hay = `${e.label} ${e.value} ${e.section} ${e.tabLabel}`.toLowerCase();
      return tokens.every((t) => hay.includes(t));
    })
    .slice(0, 40);
}
