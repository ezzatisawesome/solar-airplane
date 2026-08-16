import { useEffect, useState } from "react";
import { Link } from "react-router-dom";
import { fetchManifest } from "../lib/data";
import type { Manifest, ManifestEntry } from "../lib/types";
import { MONO } from "../lib/ui";

function fmt(v: number | null | undefined, digits = 2, suffix = "") {
  return v == null ? "—" : v.toFixed(digits) + suffix;
}

export default function Gallery() {
  const [manifest, setManifest] = useState<Manifest | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetchManifest().then(setManifest).catch((e) => setError(String(e)));
  }, []);

  if (error) return <div className="p-6 text-[var(--err)]">{error}</div>;
  if (!manifest) return <div className="p-6 text-[var(--muted)]">Loading runs…</div>;

  return (
    <div className="p-6">
      <h1 className="text-base font-semibold mb-4">Runs &amp; Airframes</h1>
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
        {manifest.runs.map((r) => (
          <RunCard key={r.id} run={r} />
        ))}
      </div>
    </div>
  );
}

function RunCard({ run }: { run: ManifestEntry }) {
  const m = run.metrics;
  return (
    <Link
      to={`/run/${run.id}`}
      className="block rounded-md border border-[var(--border)] p-4 bg-[var(--card)] hover:border-[var(--muted)] transition-colors"
    >
      <div className="flex items-start justify-between">
        <span className="font-medium">{run.label}</span>
        {run.date && <span className="text-xs text-[var(--muted)]">{run.date.slice(0, 10)}</span>}
      </div>
      <dl className="grid grid-cols-2 gap-x-3 gap-y-1 text-sm mt-2">
        <Metric label="Mass" value={fmt(m.total_mass, 2, " kg")} />
        <Metric label="CG x" value={fmt(m.x_cg, 3, " m")} />
        <Metric label="L/D" value={fmt(m.L_over_D, 1)} />
        <Metric label="Span" value={fmt(m.wingspan, 2, " m")} />
      </dl>
      {run.builds.length > 0 && (
        <div className="mt-3 flex flex-wrap gap-1">
          {run.builds.map((b) => (
            <span key={b} className="text-xs px-2 py-0.5 rounded-full bg-[var(--card-2)] text-[var(--muted)] border border-[var(--border)]">
              {b}
            </span>
          ))}
        </div>
      )}
    </Link>
  );
}

function Metric({ label, value }: { label: string; value: string }) {
  return (
    <>
      <dt className="text-[var(--muted)]">{label}</dt>
      <dd className={`text-right ${MONO}`}>{value}</dd>
    </>
  );
}
