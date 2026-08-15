import { useState } from "react";
import type { RunState } from "../lib/types";
import Field from "./Field";

type Rect = { x: number; y: number; w: number; h: number; fill?: string };
type Line = { x1: number; y1: number; x2: number; y2: number; dash?: boolean };
// A dimension: the measurement line + its label, and the aircraft edge it refers to (highlighted on hover).
type Dim = { measure: Line; label: { x: number; y: number; text: string; anchor?: "start" | "middle" | "end" }; edge: Line };

const HL = "#f59e0b"; // hover-highlight color

/**
 * 2D orthographic views (top + side) reconstructed from the solution geometry. Hovering a
 * dimension highlights the aircraft edge it measures. Top and side views share one scale.
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
  const vspan = v.vstab_span || 0;
  const vrc = v.vstab_root_chord || 0;
  const fuseLen = 0.5;
  const fuseR = run.constants?.fuselage_radius ?? 0.04;

  if (!b || !c) return <div className="p-4 text-[var(--muted)]">No geometry for this run.</div>;
  const hx = bl + vrc / 4;

  // TOP VIEW: horizontal = span (y), vertical = streamwise (x, nose at top).
  const topRects: Rect[] = [
    { x: -b / 2, y: 0, w: b, h: c },
    { x: -hs / 2, y: hx, w: hs, h: 0.12 },
  ];
  const topLines: Line[] = [
    { x1: 0, y1: 0, x2: Math.max(c, hx + 0.12), y2: 0, dash: true },
    { x1: by, y1: 0.25 * c, x2: by, y2: bl },
    { x1: -by, y1: 0.25 * c, x2: -by, y2: bl },
  ];
  const topDims: Dim[] = [
    { measure: { x1: -b / 2, y1: -0.1, x2: b / 2, y2: -0.1 }, label: { x: 0, y: -0.16, text: `b = ${b.toFixed(2)} m` }, edge: { x1: -b / 2, y1: 0, x2: b / 2, y2: 0 } },
    { measure: { x1: -b / 2 - 0.06, y1: 0, x2: -b / 2 - 0.06, y2: c }, label: { x: -b / 2 - 0.1, y: c / 2, text: `c = ${c.toFixed(2)} m`, anchor: "end" }, edge: { x1: -b / 2, y1: 0, x2: -b / 2, y2: c } },
    { measure: { x1: by + 0.06, y1: 0.25 * c, x2: by + 0.06, y2: bl }, label: { x: by + 0.1, y: (0.25 * c + bl) / 2, text: `boom ${bl.toFixed(2)} m`, anchor: "start" }, edge: { x1: by, y1: 0.25 * c, x2: by, y2: bl } },
    { measure: { x1: -hs / 2, y1: hx + 0.18, x2: hs / 2, y2: hx + 0.18 }, label: { x: 0, y: hx + 0.26, text: `hstab ${hs.toFixed(2)} m` }, edge: { x1: -hs / 2, y1: hx, x2: hs / 2, y2: hx } },
  ];

  // SIDE VIEW: horizontal = streamwise (x), vertical = height (z, up).
  const sideRects: Rect[] = [
    { x: -fuseLen, y: -fuseR, w: fuseLen, h: 2 * fuseR },
    { x: 0, y: -0.01, w: c, h: 0.02 },
    { x: bl, y: 0, w: vrc, h: vspan },
    { x: hx, y: vspan - 0.01, w: 0.12, h: 0.02 },
  ];
  const sideLines: Line[] = [{ x1: 0.25 * c, y1: -0.02, x2: bl, y2: -0.02 }];
  const sideDims: Dim[] = [
    { measure: { x1: -fuseLen, y1: -fuseR - 0.05, x2: 0, y2: -fuseR - 0.05 }, label: { x: -fuseLen / 2, y: -fuseR - 0.09, text: `fuse ${fuseLen.toFixed(2)} m` }, edge: { x1: -fuseLen, y1: fuseR, x2: 0, y2: fuseR } },
    { measure: { x1: bl + vrc + 0.05, y1: 0, x2: bl + vrc + 0.05, y2: vspan }, label: { x: bl + vrc + 0.09, y: vspan / 2, text: `vstab h ${vspan.toFixed(2)} m`, anchor: "start" }, edge: { x1: bl, y1: 0, x2: bl, y2: vspan } },
    { measure: { x1: bl, y1: -0.06, x2: bl + vrc, y2: -0.06 }, label: { x: bl + vrc / 2, y: -0.1, text: `vstab c ${vrc.toFixed(2)} m` }, edge: { x1: bl, y1: 0, x2: bl + vrc, y2: 0 } },
  ];

  const W = 900;
  const margin = 80;
  const spanH = (r: Rect[], l: Line[]) => {
    const xs = [...r.flatMap((x) => [x.x, x.x + x.w]), ...l.flatMap((x) => [x.x1, x.x2])];
    return Math.max(...xs) - Math.min(...xs);
  };
  const scale = (W - 2 * margin) / Math.max(spanH(topRects, topLines), spanH(sideRects, sideLines) || 1);

  return (
    <div className="space-y-10">
      <View title="Top view" rects={topRects} lines={topLines} dims={topDims} flipV W={W} margin={margin} scale={scale} />
      <View title="Side view" rects={sideRects} lines={sideLines} dims={sideDims} W={W} margin={margin} scale={scale} />
      <DimsTable
        dims={[
          ["Wingspan", b],
          ["Wing chord", c],
          ["H-stab span", hs],
          ["H-stab chord", h.hstab_chordlen || 0],
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

function View({
  title,
  rects,
  lines,
  dims,
  flipV = false,
  W,
  margin,
  scale,
}: {
  title: string;
  rects: Rect[];
  lines: Line[];
  dims: Dim[];
  flipV?: boolean;
  W: number;
  margin: number;
  scale: number;
}) {
  const [hover, setHover] = useState<number | null>(null);

  const xs: number[] = [];
  const ys: number[] = [];
  const push = (x: number, y: number) => {
    xs.push(x);
    ys.push(y);
  };
  for (const r of rects) {
    push(r.x, r.y);
    push(r.x + r.w, r.y + r.h);
  }
  for (const l of lines) {
    push(l.x1, l.y1);
    push(l.x2, l.y2);
  }
  for (const d of dims) {
    push(d.measure.x1, d.measure.y1);
    push(d.measure.x2, d.measure.y2);
    push(d.label.x, d.label.y);
  }
  const xmin = Math.min(...xs);
  const xmax = Math.max(...xs);
  const ymin = Math.min(...ys);
  const ymax = Math.max(...ys);

  const m = margin;
  const contentW = (xmax - xmin) * scale;
  const offX = Math.max(m, (W - contentW) / 2);
  const H = (ymax - ymin) * scale + 2 * m;
  const X = (x: number) => offX + (x - xmin) * scale;
  const Y = (y: number) => (flipV ? m + (y - ymin) * scale : H - m - (y - ymin) * scale);
  const L = (l: Line, extra: Record<string, unknown>) => <line x1={X(l.x1)} y1={Y(l.y1)} x2={X(l.x2)} y2={Y(l.y2)} {...extra} />;

  return (
    <div>
      <div className="label mb-3">{title}</div>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full max-h-[55vh]">
        {lines.map((l, i) => (
          <g key={i}>{L(l, { stroke: "var(--muted)", strokeWidth: 1.5, strokeDasharray: l.dash ? "6 6" : undefined })}</g>
        ))}
        {rects.map((r, i) => (
          <rect key={i} x={X(r.x)} y={Y(flipV ? r.y : r.y + r.h)} width={r.w * scale} height={r.h * scale} fill={r.fill ?? "var(--card-2)"} stroke="var(--ink)" strokeWidth={1.5} />
        ))}
        {dims.map((d, i) => {
          const on = hover === i;
          return (
            <g key={i} onMouseEnter={() => setHover(i)} onMouseLeave={() => setHover(null)} style={{ cursor: "default" }}>
              {/* highlighted edge on the aircraft */}
              {on && L(d.edge, { stroke: HL, strokeWidth: 4 })}
              {/* dimension line + label */}
              {L(d.measure, { stroke: on ? HL : "var(--muted)", strokeWidth: on ? 2 : 1 })}
              <text x={X(d.label.x)} y={Y(d.label.y)} fill={on ? HL : "var(--muted)"} fontSize="12" textAnchor={d.label.anchor ?? "middle"} className="mono">
                {d.label.text}
              </text>
              {/* invisible fat hit-line for easy hovering */}
              {L(d.measure, { stroke: "transparent", strokeWidth: 16 })}
            </g>
          );
        })}
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
