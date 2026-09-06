import type { Airfoil, RunState } from "../lib/types";
import { LABEL } from "../lib/ui";

const PITCH = 0.13; // solar-cell pitch (m), matches solarCells()
const HL = "#3d6d99"; // solar-cell blue

/** Which lifting surface a foil belongs to, for solar coverage. */
function surfaceKind(surface: string): "wing" | "hstab" | null {
  const s = surface.toLowerCase();
  if (s.includes("wing")) return "wing";
  if (s.includes("hstab") || s.includes("horizontal") || s.includes("tail")) return "hstab";
  return null;
}

/**
 * 2D airfoil section outlines (Selig coordinates) for each lifting surface, with the quarter-chord
 * marked and the solar-cell chordwise coverage band drawn on the upper surface.
 */
export default function Airfoils({ run }: { run: RunState }) {
  const foils = run.airfoils ?? [];
  if (!foils.length) return <div className="p-4 text-[var(--muted)]">No airfoil data for this run.</div>;

  const cRoot = run.soln["Main Wing"]?.chordlen ?? 0;
  const hChord = run.soln.HStab?.hstab_chordlen ?? 0;

  // Chordwise fraction covered by cells, per surface (matches the packing in solarCells()).
  const cover = (kind: "wing" | "hstab" | null): [number, number] | null => {
    if (kind === "wing" && cRoot > 0) {
      const rows = Math.floor(cRoot / PITCH);
      return [0, Math.min(1, (rows * PITCH) / cRoot)]; // LE back to last full row
    }
    if (kind === "hstab" && hChord > 0) return [0.1, 0.7]; // 10% LE margin, 60% band
    return null;
  };

  return (
    <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
      {foils.map((f) => (
        <AirfoilPlot key={f.surface} foil={f} coverage={cover(surfaceKind(f.surface))} />
      ))}
    </div>
  );
}

function AirfoilPlot({ foil, coverage }: { foil: Airfoil; coverage: [number, number] | null }) {
  const pts = foil.coords;
  const xs = pts.map((p) => p[0]);
  const ys = pts.map((p) => p[1]);
  const minX = Math.min(...xs);
  const maxX = Math.max(...xs);
  const minY = Math.min(...ys);
  const maxY = Math.max(...ys);
  const W = 640;
  const H = 200;
  const pad = 24;
  const s = Math.min((W - 2 * pad) / (maxX - minX || 1), (H - 2 * pad) / (maxY - minY || 1));
  const tx = (x: number) => pad + (x - minX) * s;
  const ty = (y: number) => H / 2 - (y - (minY + maxY) / 2) * s;
  const d = pts.map((p, i) => `${i ? "L" : "M"}${tx(p[0]).toFixed(1)},${ty(p[1]).toFixed(1)}`).join(" ") + " Z";
  const chord = maxX - minX || 1;
  const c4x = minX + 0.25 * chord;

  // Upper-surface points inside the solar coverage band (normalized chord fraction).
  let solarPath = "";
  if (coverage) {
    const [f0, f1] = coverage;
    const upper = pts
      .filter((p) => p[1] >= 0)
      .map((p) => [p, (p[0] - minX) / chord] as const)
      .filter(([, xn]) => xn >= f0 && xn <= f1)
      .sort((a, b) => a[1] - b[1]);
    if (upper.length >= 2) {
      solarPath = upper.map(([p], i) => `${i ? "L" : "M"}${tx(p[0]).toFixed(1)},${(ty(p[1]) - 2).toFixed(1)}`).join(" ");
    }
  }

  return (
    <div>
      <div className={`${LABEL} mb-2`}>
        {foil.surface} — <span className="text-[var(--ink)]">{foil.name}</span>
      </div>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full">
        <line x1={tx(minX)} y1={ty(0)} x2={tx(maxX)} y2={ty(0)} stroke="var(--border)" strokeDasharray="4 4" />
        <line x1={tx(c4x)} y1={ty(minY)} x2={tx(c4x)} y2={ty(maxY)} stroke="var(--accent)" strokeDasharray="3 3" opacity={0.7} />
        <text x={tx(c4x) + 4} y={pad - 6} fill="var(--accent)" fontSize={11}>c/4</text>
        <path d={d} fill="var(--card-2)" stroke="var(--ink)" strokeWidth={1.5} />
        {solarPath && <path d={solarPath} fill="none" stroke={HL} strokeWidth={4} strokeLinecap="round" />}
      </svg>
      {coverage && <div className="text-[11px] text-[var(--muted)] mt-1">Solar cells cover the highlighted upper-surface band.</div>}
    </div>
  );
}
