import type { RunState } from "../lib/types";
import Field from "./Field";

type Rect = { x: number; y: number; w: number; h: number; fill?: string };
type Line = { x1: number; y1: number; x2: number; y2: number; dash?: boolean };
type Label = { x: number; y: number; text: string; anchor?: "start" | "middle" | "end" };

/**
 * 2D orthographic views of the airframe (top + side) reconstructed from the solution geometry,
 * annotated with the main dimensions of the wing, stabilizers, boom and fuselage.
 */
export default function Planform({ run }: { run: RunState }) {
  const s = run.soln;
  const mw = s["Main Wing"] || {};
  const g = s.Geometry || {};
  const h = s.HStab || {};
  const v = s["V Stab"] || {};

  const b = mw.wingspan || 0;
  const c = mw.chordlen || 0;
  const bl = g.boom_length || 0;
  const by = g.boom_y || 0;
  const hs = h.hstab_span || 0;
  const hc = h.hstab_chordlen || 0;
  const vspan = v.vstab_span || 0;
  const vrc = v.vstab_root_chord || 0;
  const fuseLen = 0.5; // matches the mass model
  const fuseR = run.constants?.fuselage_radius ?? 0.04;

  if (!b || !c) return <div className="p-4 text-[var(--muted)]">No geometry for this run.</div>;

  const hx = bl + vrc / 4; // hstab streamwise station

  // ---- TOP VIEW: horizontal = span (y), vertical = streamwise (x, nose at top) ----
  const topRects: Rect[] = [
    { x: -b / 2, y: 0, w: b, h: c }, // main wing
    { x: -hs / 2, y: hx, w: hs, h: hc }, // hstab
  ];
  const topLines: Line[] = [
    { x1: 0, y1: 0, x2: 0, y2: Math.max(c, hx + hc), dash: true }, // centerline
    { x1: by, y1: 0.25 * c, x2: by, y2: bl }, // boom R
    { x1: -by, y1: 0.25 * c, x2: -by, y2: bl }, // boom L
  ];
  const topLabels: Label[] = [
    { x: 0, y: -0.06 * c, text: `b = ${b.toFixed(2)} m`, anchor: "middle" },
    { x: -b / 2 - 0.02, y: c / 2, text: `c = ${c.toFixed(2)} m`, anchor: "end" },
    { x: by + 0.02, y: (0.25 * c + bl) / 2, text: `boom ${bl.toFixed(2)} m`, anchor: "start" },
    { x: 0, y: hx + hc + 0.08, text: `hstab ${hs.toFixed(2)}×${hc.toFixed(2)} m`, anchor: "middle" },
  ];

  // ---- SIDE VIEW: horizontal = streamwise (x), vertical = height (z, up) ----
  const sideRects: Rect[] = [
    { x: -fuseLen, y: -fuseR, w: fuseLen, h: 2 * fuseR }, // fuselage
    { x: 0, y: -0.01, w: c, h: 0.02 }, // wing (edge-on)
    { x: bl, y: 0, w: vrc, h: vspan }, // vstab
    { x: hx, y: vspan - 0.01, w: hc, h: 0.02 }, // hstab at top of vstab
  ];
  const sideLines: Line[] = [
    { x1: 0.25 * c, y1: -0.02, x2: bl, y2: -0.02 }, // boom
    { x1: -fuseLen, y1: 0, x2: bl + vrc, y2: 0, dash: true }, // waterline
  ];
  const sideLabels: Label[] = [
    { x: -fuseLen / 2, y: -fuseR - 0.04, text: `fuse ${fuseLen.toFixed(2)} m`, anchor: "middle" },
    { x: bl + vrc + 0.03, y: vspan / 2, text: `vstab h ${vspan.toFixed(2)} m`, anchor: "start" },
    { x: bl + vrc / 2, y: -0.06, text: `vstab c ${vrc.toFixed(2)} m`, anchor: "middle" },
  ];

  // Shared scale so the top and side views are drawn to the same scale (a metre is the same size
  // in both). Driven by the widest horizontal extent across the two views.
  const W = 900;
  const margin = 70;
  const spanH = (r: Rect[], l: Line[]) => {
    const xs = [...r.flatMap((x) => [x.x, x.x + x.w]), ...l.flatMap((x) => [x.x1, x.x2])];
    return Math.max(...xs) - Math.min(...xs);
  };
  const scale = (W - 2 * margin) / Math.max(spanH(topRects, topLines), spanH(sideRects, sideLines) || 1);

  return (
    <div className="space-y-10">
      <View title="Top view" rects={topRects} lines={topLines} labels={topLabels} flipV W={W} margin={margin} scale={scale} />
      <View title="Side view" rects={sideRects} lines={sideLines} labels={sideLabels} W={W} margin={margin} scale={scale} />
      <DimsTable
        dims={[
          ["Wingspan", b],
          ["Wing chord", c],
          ["H-stab span", hs],
          ["H-stab chord", hc],
          ["V-stab height", vspan],
          ["V-stab root chord", vrc],
          ["Boom length", bl],
          ["Boom Y", by],
          ["Fuselage length", fuseLen],
          ["Fuselage radius", fuseR],
        ]}
      />
    </div>
  );
}

/** Generic auto-fitting 2D view. `flipV` renders the vertical axis increasing downward (top view). */
function View({
  title,
  rects,
  lines,
  labels,
  flipV = false,
  W,
  margin,
  scale,
}: {
  title: string;
  rects: Rect[];
  lines: Line[];
  labels: Label[];
  flipV?: boolean;
  W: number;
  margin: number;
  scale: number;
}) {
  const xs: number[] = [];
  const ys: number[] = [];
  for (const r of rects) {
    xs.push(r.x, r.x + r.w);
    ys.push(r.y, r.y + r.h);
  }
  for (const l of lines) {
    xs.push(l.x1, l.x2);
    ys.push(l.y1, l.y2);
  }
  const xmin = Math.min(...xs);
  const xmax = Math.max(...xs);
  const ymin = Math.min(...ys);
  const ymax = Math.max(...ys);

  const m = margin;
  const contentW = (xmax - xmin) * scale;
  const offX = Math.max(m, (W - contentW) / 2); // center horizontally so both views align to scale
  const H = (ymax - ymin) * scale + 2 * m;
  const X = (x: number) => offX + (x - xmin) * scale;
  const Y = (y: number) => (flipV ? m + (y - ymin) * scale : H - m - (y - ymin) * scale);

  return (
    <div>
      <div className="label mb-3">{title}</div>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full max-h-[55vh]">
        {lines.map((l, i) => (
          <line key={i} x1={X(l.x1)} y1={Y(l.y1)} x2={X(l.x2)} y2={Y(l.y2)} stroke="var(--muted)" strokeWidth={1.5} strokeDasharray={l.dash ? "6 6" : undefined} />
        ))}
        {rects.map((r, i) => (
          <rect key={i} x={X(r.x)} y={Y(flipV ? r.y : r.y + r.h)} width={r.w * scale} height={r.h * scale} fill={r.fill ?? "var(--card-2)"} stroke="var(--ink)" strokeWidth={1.5} />
        ))}
        {labels.map((t, i) => (
          <text key={i} x={X(t.x)} y={Y(t.y)} fill="var(--muted)" fontSize="12" textAnchor={t.anchor ?? "middle"} className="mono">
            {t.text}
          </text>
        ))}
      </svg>
    </div>
  );
}

function DimsTable({ dims }: { dims: [string, number][] }) {
  return (
    <dl className="grid grid-cols-2 sm:grid-cols-3 gap-x-8 gap-y-1 text-sm max-w-3xl">
      {dims.map(([k, val]) => (
        <Field key={k} label={k} value={`${val.toFixed(3)} m`} />
      ))}
    </dl>
  );
}
