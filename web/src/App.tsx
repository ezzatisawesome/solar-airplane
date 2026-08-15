import { useEffect, useState } from "react";
import { Link, Navigate, Route, Routes, useMatch } from "react-router-dom";
import Gallery from "./views/Gallery";
import RunDetail from "./views/RunDetail";
import RunSwitcher from "./components/RunSwitcher";
import { fetchManifest } from "./lib/data";
import { useRunState } from "./lib/store";
import { downloadJSON } from "./lib/download";

// Where the AircraftSim UI runs (override with VITE_AIRCRAFTSIM_URL at build time).
const AIRCRAFTSIM_URL = import.meta.env.VITE_AIRCRAFTSIM_URL || "http://localhost:3001";

function PlaneIcon() {
  return (
    <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
      <path d="M21 16v-2l-8-5V3.5A1.5 1.5 0 0 0 11.5 2 1.5 1.5 0 0 0 10 3.5V9l-8 5v2l8-2.5V19l-2 1.5V22l3.5-1 3.5 1v-1.5L13 19v-5.5l8 2.5Z" />
    </svg>
  );
}

// Simulate the currently-viewed aircraft in AircraftSim (only shown on a run page).
function SimulateButton() {
  const match = useMatch("/run/:id");
  const id = match?.params.id;
  if (!id) return null;
  return (
    <button
      className="btn btn-accent h-7 px-2.5 text-xs flex items-center gap-1.5"
      onClick={() => window.open(`${AIRCRAFTSIM_URL}/?aircraft=${encodeURIComponent(id)}`, "_blank", "noopener")}
    >
      <PlaneIcon />
      Simulate
    </button>
  );
}

function ThemeToggle() {
  const [light, setLight] = useState(() => localStorage.getItem("theme") === "light");
  useEffect(() => {
    document.documentElement.classList.toggle("light", light);
    localStorage.setItem("theme", light ? "light" : "dark");
  }, [light]);
  return (
    <button className="btn h-7 w-7 px-0 flex items-center justify-center" onClick={() => setLight((v) => !v)} title="Toggle theme">
      {light ? "☾" : "☀"}
    </button>
  );
}

// Current run name, centered in the header.
function HeaderTitle() {
  const match = useMatch("/run/:id");
  if (!match?.params.id) return null;
  return (
    <div className="absolute left-1/2 -translate-x-1/2 text-[15px] font-semibold pointer-events-none">
      {match.params.id}
    </div>
  );
}

function EyeIcon() {
  return (
    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="M2 12s3.5-7 10-7 10 7 10 7-3.5 7-10 7-10-7-10-7Z" />
      <circle cx="12" cy="12" r="3" />
    </svg>
  );
}

// Download the current run's state (including any live dragged layout) from the header.
function DownloadButton() {
  const { run, layout, total } = useRunState();
  if (!run) return null;
  return (
    <button
      className="btn h-7 px-2.5 text-xs"
      onClick={() =>
        downloadJSON(`${run.run}.state.json`, {
          run: run.run,
          layout,
          total: total ?? run.massProperties.total,
        })
      }
    >
      Download
    </button>
  );
}

// Landing: jump straight to the latest run's 3D viewer (no selection page first).
function HomeRedirect() {
  const [latest, setLatest] = useState<string | null>(null);
  const [error, setError] = useState(false);
  useEffect(() => {
    fetchManifest()
      .then((m) => {
        const runs = m.runs;
        setLatest(runs.length ? runs[runs.length - 1].id : null);
      })
      .catch(() => setError(true));
  }, []);
  if (error) return <div className="p-6 text-red-400">Could not load runs. Run the exporter.</div>;
  if (!latest) return <div className="p-6 text-neutral-400">Loading…</div>;
  return <Navigate to={`/run/${latest}`} replace />;
}

export default function App() {
  return (
    <div className="h-screen flex flex-col overflow-hidden">
      <header className="relative flex items-center gap-3 px-5 py-2 border-b border-[var(--border)] bg-[var(--bg)]">
        <Link to="/" className="flex items-center gap-2 font-semibold">
          <EyeIcon />
          AircraftView
        </Link>
        <img src="power.png" alt="power" className="h-6 w-auto" />
        <Link to="/runs" className="text-sm text-[var(--muted)] hover:text-[var(--ink)]">
          All runs
        </Link>
        <RunSwitcher />
        <HeaderTitle />
        <div className="ml-auto flex items-center gap-2">
          <DownloadButton />
          <SimulateButton />
          <ThemeToggle />
        </div>
      </header>
      <main className="flex-1 min-h-0">
        <Routes>
          <Route path="/" element={<HomeRedirect />} />
          <Route path="/run/:id" element={<RunDetail />} />
          <Route path="/runs" element={<Gallery />} />
        </Routes>
      </main>
    </div>
  );
}
