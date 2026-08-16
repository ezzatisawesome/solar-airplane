import { useEffect, useState } from "react";
import { useMatch, useNavigate } from "react-router-dom";
import { fetchManifest } from "../lib/data";
import type { ManifestEntry } from "../lib/types";

/** Compact dropdown in the top bar to pick which run/aircraft to view. */
export default function RunSwitcher() {
  const [runs, setRuns] = useState<ManifestEntry[]>([]);
  const navigate = useNavigate();
  const match = useMatch("/run/:id");
  const currentId = match?.params.id ?? "";

  useEffect(() => {
    fetchManifest()
      .then((m) => setRuns(m.runs))
      .catch(() => setRuns([]));
  }, []);

  if (!runs.length) return null;

  return (
    <label className="flex items-center gap-2 text-sm">
      <select
        value={currentId}
        onChange={(e) => navigate(`/run/${e.target.value}`)}
        className="bg-[var(--card-2)] border border-[var(--border)] rounded-md px-2 py-1 text-[var(--ink)] max-w-[9rem] sm:max-w-none"
      >
        {runs.map((r) => (
          <option key={r.id} value={r.id}>
            {r.label}
            {r.builds.length ? ` · ${r.builds.length} built` : ""}
          </option>
        ))}
      </select>
    </label>
  );
}
