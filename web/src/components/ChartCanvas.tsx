import { useEffect, useRef } from "react";
import {
  BarController, BarElement, CategoryScale, Chart, type ChartConfiguration,
  Filler, Legend, LinearScale, LineController, LineElement, PointElement, Tooltip,
} from "chart.js";
import { LABEL } from "../lib/ui";

// Register every controller/element used anywhere in the app, once.
Chart.register(
  LineController, LineElement, PointElement,
  BarController, BarElement,
  LinearScale, CategoryScale, Filler, Tooltip, Legend,
);

/** Shared Chart.js canvas — owns the chart lifecycle (create on config change, destroy on unmount). */
export default function ChartCanvas({
  config,
  title,
  heightClass = "h-72",
}: {
  config: ChartConfiguration;
  title?: string;
  heightClass?: string;
}) {
  const ref = useRef<HTMLCanvasElement>(null);
  useEffect(() => {
    if (!ref.current) return;
    const chart = new Chart(ref.current, config);
    return () => chart.destroy();
  }, [config]);
  return (
    <div>
      {title && <h3 className={`${LABEL} mb-3`}>{title}</h3>}
      <div className={heightClass}>
        <canvas ref={ref} />
      </div>
    </div>
  );
}
