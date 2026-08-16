import type { Airfoil, RunState } from "../lib/types";
import { LABEL } from "../lib/ui";

/** 2D airfoil section outlines (Selig coordinates) for each lifting surface. */
export default function Airfoils({ run }: { run: RunState }) {
  const foils = run.airfoils ?? [];
  if (!foils.length) return <div className="p-4 text-[var(--muted)]">No airfoil data for this run.</div>;
  return (
    <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
      {foils.map((f) => (
        <AirfoilPlot key={f.surface} foil={f} />
      ))}
    </div>
  );
}

function AirfoilPlot({ foil }: { foil: Airfoil }) {
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
  return (
    <div>
      <div className={`${LABEL} mb-2`}>
        {foil.surface} — <span className="text-[var(--ink)]">{foil.name}</span>
      </div>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full">
        <line x1={tx(minX)} y1={ty(0)} x2={tx(maxX)} y2={ty(0)} stroke="var(--border)" strokeDasharray="4 4" />
        <path d={d} fill="var(--card-2)" stroke="var(--ink)" strokeWidth={1.5} />
      </svg>
    </div>
  );
}
