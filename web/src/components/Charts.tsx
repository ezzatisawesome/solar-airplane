import { useMemo } from "react";
import type { ChartConfiguration } from "chart.js";
import type { RunState } from "../lib/types";
import ChartCanvas from "./ChartCanvas";

export function BatteryChart({ run }: { run: RunState }) {
  const battery: number[] = run.soln.Power?.battery_states ?? [];
  const capacity: number = run.soln.Power?.battery_capacity ?? 0;

  const config = useMemo<ChartConfiguration | null>(() => {
    if (!battery.length) return null;
    const hours = battery.map((_, i) => ((i / (battery.length - 1)) * 24).toFixed(1));
    return {
      type: "line",
      data: {
        labels: hours,
        datasets: [
          {
            label: "Battery energy (Wh)",
            data: battery,
            borderColor: "#e8e8e8",
            backgroundColor: "rgba(232,232,232,0.08)",
            fill: true,
            pointRadius: 0,
            borderWidth: 2,
            tension: 0.25,
          },
        ],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          x: { title: { display: true, text: "Mission time (h)", color: "#888" }, ticks: { color: "#888", maxTicksLimit: 13 }, grid: { color: "rgba(128,128,128,0.12)" } },
          y: { beginAtZero: true, title: { display: true, text: "Battery (Wh)", color: "#888" }, ticks: { color: "#888" }, grid: { color: "rgba(128,128,128,0.12)" }, suggestedMax: capacity || undefined },
        },
        plugins: { legend: { labels: { color: "#bbb" } } },
      },
    };
  }, [battery, capacity]);

  if (!config) return <div className="text-[var(--muted)] text-sm">No battery time-series for this run.</div>;
  return <ChartCanvas title="Battery energy throughout flight" config={config} />;
}
