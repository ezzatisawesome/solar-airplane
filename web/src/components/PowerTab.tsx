import { useMemo } from "react";
import type { ChartConfiguration, Plugin } from "chart.js";
import type { RunState, SpecTab } from "../lib/types";
import { energyMarch, sunGeometry } from "../lib/solar";
import ChartCanvas from "./ChartCanvas";
import Field from "./Field";
import { LABEL } from "../lib/ui";

function cssVar(name: string, fallback: string): string {
  if (typeof window === "undefined") return fallback;
  const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return v || fallback;
}

/**
 * The Power tab: the battery state-of-charge march (night shaded, dawn low + 20% reserve marked),
 * followed by the detailed power / solar / propulsion spec tables.
 */
export default function PowerTab({ run, specTab }: { run: RunState; specTab?: SpecTab }) {
  const e = useMemo(() => energyMarch(run.soln), [run]);
  const lat = run.soln.Mission?.operating_lat ?? 0;
  const day = run.soln.Mission?.mission_date ?? 172;
  const sun = useMemo(() => sunGeometry(lat, day), [lat, day]);

  const socConfig = useMemo<ChartConfiguration | null>(() => {
    if (!e.hasData) return null;
    const accent = cssVar("--accent", "#8ab4f8");
    const bad = cssVar("--err", "#f28b82");
    const warn = cssVar("--warn", "#fdd663");
    const grid = "rgba(128,128,128,0.14)";
    const muted = cssVar("--muted", "#9aa0a6");

    const nightPlugin: Plugin = {
      id: "nightBands",
      beforeDatasetsDraw(chart) {
        const { ctx, chartArea } = chart;
        const x = chart.scales.x;
        if (!x || !chartArea) return;
        ctx.save();
        ctx.fillStyle = "rgba(90,110,150,0.14)";
        for (const [a, b] of [[0, sun.sunriseHod], [sun.sunsetHod, 24]] as [number, number][]) {
          if (b <= a) continue;
          const xa = x.getPixelForValue(a);
          const xb = x.getPixelForValue(b);
          ctx.fillRect(xa, chartArea.top, xb - xa, chartArea.bottom - chartArea.top);
        }
        ctx.restore();
      },
    };

    return {
      type: "line",
      data: {
        datasets: [
          {
            label: "State of charge",
            data: e.hours.map((h, i) => ({ x: h, y: 100 * e.soc[i] })),
            borderColor: accent,
            backgroundColor: "rgba(128,128,128,0.15)",
            fill: true,
            pointRadius: 0,
            borderWidth: 2,
            tension: 0.2,
          },
          {
            label: "20% reserve",
            data: [{ x: 0, y: 100 * e.floor }, { x: 24, y: 100 * e.floor }],
            borderColor: bad,
            borderDash: [6, 4],
            pointRadius: 0,
            borderWidth: 1,
          },
          {
            label: "Dawn low",
            data: [{ x: e.socMinHour, y: 100 * e.socMin }],
            borderColor: warn,
            backgroundColor: warn,
            pointRadius: 5,
            showLine: false,
          },
        ],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        interaction: { mode: "index", intersect: false },
        scales: {
          x: {
            type: "linear",
            min: 0,
            max: 24,
            title: { display: true, text: "Local time (h) — night shaded", color: muted },
            ticks: { color: muted, stepSize: 3, callback: (v) => `${v}:00` },
            grid: { color: grid },
          },
          y: {
            min: 0,
            max: 100,
            title: { display: true, text: "SOC %", color: muted },
            ticks: { color: muted, callback: (v) => `${v}%` },
            grid: { color: grid },
          },
        },
        plugins: { legend: { labels: { color: cssVar("--ink", "#e8eaed"), boxWidth: 12 } } },
      },
      plugins: [nightPlugin],
    };
  }, [e, sun]);

  return (
    <div className="space-y-10">
      {socConfig ? (
        <ChartCanvas title="Battery state of charge" config={socConfig} heightClass="h-80" />
      ) : (
        <div className="text-[var(--muted)] text-sm">No battery time-series for this run.</div>
      )}

      {specTab && (
        <div className="grid grid-cols-1 md:grid-cols-2 gap-x-10 gap-y-6">
          {specTab.sections.filter((s) => s.rows.length > 0).map((s) => (
            <div key={s.title}>
              <h3 className={`${LABEL} mb-2`}>{s.title}</h3>
              <dl className="text-sm">
                {s.rows.map(([label, value]) => (
                  <Field key={label} label={label} value={value} />
                ))}
              </dl>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
