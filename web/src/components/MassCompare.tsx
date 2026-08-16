import { useEffect, useMemo, useState } from "react";
import type { ChartConfiguration } from "chart.js";
import type { RunState } from "../lib/types";
import { MASS_COMPONENTS } from "../lib/compare";
import ChartCanvas from "./ChartCanvas";
import { BTN, BTN_ACCENT, LABEL, MONO } from "../lib/ui";

// One color per component "floor" of the tower.
const PALETTE = [
  "#4a90e2", "#e74c3c", "#27ae60", "#e67e22", "#f1c40f", "#9b59b6",
  "#1abc9c", "#e91e63", "#8b7355", "#3498db", "#2ecc71", "#95a5a6",
  "#d35400", "#16a085", "#c0392b", "#7f8c8d", "#f39c12",
];
const DEV = import.meta.env.DEV;

interface EditBuild {
  name: string;
  masses: Record<string, number>;
}

/**
 * Expected (design) vs as-built masses. Horizontal grouped bar chart + a comparison table whose
 * as-built columns are editable (dev only) and save straight to output/<run>/builds.json.
 */
export default function MassCompare({ run }: { run: RunState }) {
  const design = (run.soln.Masses ?? {}) as Record<string, number>;
  const [builds, setBuilds] = useState<EditBuild[]>(() =>
    (run.builds ?? []).map((b) => ({ name: b.name, masses: { ...b.masses } })),
  );
  const [saving, setSaving] = useState(false);
  const [msg, setMsg] = useState<string | null>(null);

  // Reset editable state when switching runs.
  useEffect(() => {
    setBuilds((run.builds ?? []).map((b) => ({ name: b.name, masses: { ...b.masses } })));
    setMsg(null);
  }, [run.run]);

  const rows = useMemo(
    () =>
      MASS_COMPONENTS.map((c) => ({
        key: c.key,
        label: c.label,
        design: Number(design[c.key]) || 0,
        builds: builds.map((b) => Number(b.masses[c.key]) || 0),
      }))
        .filter((r) => r.design > 0 || r.builds.some((v) => v > 0))
        .sort((a, b) => Math.max(b.design, ...b.builds) - Math.max(a.design, ...a.builds)),
    [design, builds],
  );

  // One tower per entity (Expected + each build); each component is a stacked "floor" and the
  // total tower height is the total mass.
  const config = useMemo<ChartConfiguration>(
    () => ({
      type: "bar",
      data: {
        labels: ["Expected", ...builds.map((b) => b.name)],
        datasets: rows.map((r, idx) => ({
          label: r.label,
          data: [r.design, ...r.builds],
          backgroundColor: PALETTE[idx % PALETTE.length],
          categoryPercentage: 0.55,
          barPercentage: 0.9,
        })),
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          x: { stacked: true, ticks: { color: "#ccc", font: { size: 13 } }, grid: { display: false } },
          y: { stacked: true, beginAtZero: true, title: { display: true, text: "Total mass (kg)", color: "#888" }, ticks: { color: "#888" }, grid: { color: "rgba(128,128,128,0.15)" } },
        },
        plugins: {
          legend: { position: "bottom", labels: { color: "#bbb", boxWidth: 12, font: { size: 10 } } },
          tooltip: { callbacks: { label: (c) => `${c.dataset.label}: ${(c.raw as number).toFixed(3)} kg` } },
        },
      },
    }),
    [rows, builds],
  );

  const setMass = (bi: number, key: string, value: string) =>
    setBuilds((prev) => prev.map((b, i) => (i === bi ? { ...b, masses: { ...b.masses, [key]: Number(value) } } : b)));

  const addBuild = () =>
    setBuilds((prev) => [
      ...prev,
      { name: `Build ${prev.length + 1}`, masses: Object.fromEntries(MASS_COMPONENTS.map((c) => [c.key, Number(design[c.key]) || 0])) },
    ]);

  const save = async () => {
    setSaving(true);
    setMsg(null);
    try {
      const res = await fetch(`/api/builds/${run.run}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ builds }),
      });
      if (!res.ok) throw new Error(await res.text());
      setMsg("Saved. Reloading…");
      setTimeout(() => window.location.reload(), 500);
    } catch (e) {
      setMsg(`Save failed: ${e}`);
      setSaving(false);
    }
  };

  if (!rows.length) return <div className="text-sm text-[var(--muted)]">No mass data.</div>;

  return (
    <div className="grid grid-cols-1 md:grid-cols-2 gap-x-10 gap-y-6 items-start">
      <div className="min-w-0">
        <ChartCanvas title="Total mass — expected vs as-built" config={config} heightClass="h-[460px]" />
      </div>

      <div className="min-w-0">
        <div className="flex items-center gap-3 mb-2">
          <h3 className={`${LABEL}`}>Comparison {DEV && "(edit as-built values, then save)"}</h3>
          {DEV && (
            <>
              <button className={BTN} onClick={addBuild}>+ Add build</button>
              <button className={BTN_ACCENT} onClick={save} disabled={saving}>Save to builds.json</button>
              {msg && <span className="text-xs text-[var(--muted)]">{msg}</span>}
            </>
          )}
        </div>
        <table className="text-sm w-full">
          <thead>
            <tr className="text-[var(--muted)] border-b border-[var(--border)]">
              <th className="text-left font-medium py-1">Component</th>
              <th className="text-right font-medium py-1">Expected</th>
              {builds.map((b, i) => (
                <th key={i} className="text-right font-medium py-1">{b.name}</th>
              ))}
              {builds.length === 1 && <th className="text-right font-medium py-1">Δ / %</th>}
            </tr>
          </thead>
          <tbody>
            {rows.map((r) => (
              <tr key={r.key} className="border-b border-[var(--border)]/40">
                <td className="py-1 text-[var(--muted)]">{r.label}</td>
                <td className={`py-1 text-right ${MONO}`}>{r.design.toFixed(3)}</td>
                {builds.map((_, i) => (
                  <td key={i} className="py-1 text-right">
                    {DEV ? (
                      <input
                        type="number"
                        step="0.001"
                        value={r.builds[i]}
                        onChange={(e) => setMass(i, r.key, e.target.value)}
                        className={`bg-[var(--card-2)] border border-[var(--border)] rounded px-2 py-0.5 text-right w-24 ${MONO}`}
                      />
                    ) : (
                      <span className={`${MONO}`}>{r.builds[i].toFixed(3)}</span>
                    )}
                  </td>
                ))}
                {builds.length === 1 && <DeltaCell design={r.design} build={r.builds[0]} />}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function DeltaCell({ design, build }: { design: number; build: number }) {
  const d = build - design;
  const pct = design ? (d / design) * 100 : build > 0 ? Infinity : 0;
  const cls = Math.abs(d) < 1e-6 ? "text-[var(--muted)]" : d > 0 ? "text-[var(--err)]" : "text-[var(--ok)]";
  return (
    <td className={`py-1 text-right ${MONO} ${cls}`}>
      {d > 0 ? "+" : ""}{d.toFixed(3)} {Number.isFinite(pct) ? `(${pct > 0 ? "+" : ""}${pct.toFixed(1)}%)` : ""}
    </td>
  );
}
