import type { Manifest, RunState } from "./types";

// All data is static under public/ (copied to the bundle root at build). Use a relative base so
// the site works from any host or subpath.
const base = import.meta.env.BASE_URL || "/";

async function getJSON<T>(path: string): Promise<T> {
  const res = await fetch(base + path);
  if (!res.ok) throw new Error(`Failed to load ${path}: ${res.status}`);
  return res.json();
}

export const fetchManifest = () => getJSON<Manifest>("manifest.json");
export const fetchRun = (id: string) => getJSON<RunState>(`runs/${id}.json`);
