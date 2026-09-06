import { useCallback, useEffect, useMemo, useState } from "react";
import { useParams } from "react-router-dom";
import { fetchRun } from "../lib/data";
import type { ComponentMP, RunState, SpecTab } from "../lib/types";
import type { Layout } from "../lib/mass";
import AircraftView from "../components/AircraftView";
import MassCompare from "../components/MassCompare";
import Airfoils from "../components/Airfoils";
import Planform from "../components/Planform";
import OverviewTab from "../components/OverviewTab";
import PowerTab from "../components/PowerTab";
import Field from "../components/Field";
import { setRunState } from "../lib/store";
import { buildIndex, searchIndex, type SearchEntry } from "../lib/search";
import { LABEL, MONO, TAB, TAB_ACTIVE } from "../lib/ui";

export default function RunDetail() {
  const { id = "" } = useParams();
  const [run, setRun] = useState<RunState | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [tab, setTab] = useState("view");
  const [query, setQuery] = useState("");

  useEffect(() => {
    let alive = true;
    setRun(null);
    setError(null);
    setTab("view");
    setQuery("");
    document.title = `AircraftView - ${id}`;
    fetchRun(id)
      .then((r) => {
        if (!alive) return;
        setRun(r);
        setRunState({ run: r, layout: {}, total: r.massProperties.total });
      })
      .catch((e) => alive && setError(String(e)));
    return () => {
      alive = false;
    };
  }, [id]);

  const onLayoutChange = useCallback((layout: Layout, total: ComponentMP) => {
    setRunState({ layout, total });
  }, []);

  const tabs = useMemo(() => {
    const spec = (run?.specs ?? []).map((s) => ({ id: s.id, label: s.label }));
    return [{ id: "view", label: "AircraftView" }, ...spec, { id: "shapes", label: "Airfoils & Planform" }];
  }, [run]);

  const index = useMemo(() => (run ? buildIndex(run) : []), [run]);
  const results = useMemo(() => searchIndex(index, query), [index, query]);

  if (error) return <div className="p-6 text-[var(--err)]">{error}</div>;
  if (!run) return <div className="p-6 text-[var(--muted)]">Loading {id}…</div>;

  const specTab = run.specs.find((s) => s.id === tab);
  const jump = (e: SearchEntry) => {
    setTab(e.tab);
    setQuery("");
  };

  return (
    <div className="h-full flex flex-col min-h-0 px-4 pb-4 pt-3 gap-3">
      <div className="flex flex-wrap items-center gap-3 border-b border-[var(--border)] pb-2">
        <nav className="flex flex-wrap gap-1">
          {tabs.map((t) => (
            <button key={t.id} onClick={() => setTab(t.id)} className={`${TAB} ${tab === t.id ? TAB_ACTIVE : ""}`}>
              {t.label}
            </button>
          ))}
        </nav>
        <div className="relative ml-auto">
          <SearchIcon />
          <input
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="Search all data…"
            className="h-8 w-52 sm:w-64 pl-7 pr-2 rounded-md bg-[var(--card-2)] border border-[var(--border)] text-[13px] text-[var(--ink)] placeholder:text-[var(--faint)] focus:outline-none focus:border-[var(--accent)]"
          />
          {query.trim() && (
            <div className="absolute right-0 z-30 mt-1 w-80 max-h-96 overflow-y-auto rounded-md border border-[var(--border)] bg-[var(--card)] shadow-lg">
              {results.length === 0 ? (
                <div className="px-3 py-2 text-xs text-[var(--muted)]">No matches for “{query}”.</div>
              ) : (
                results.map((e, i) => (
                  <button
                    key={i}
                    onClick={() => jump(e)}
                    className="w-full text-left px-3 py-1.5 hover:bg-[var(--card-2)] border-b border-[var(--border)]/40 last:border-0"
                  >
                    <div className="flex justify-between gap-3 text-[13px]">
                      <span className="text-[var(--ink)] truncate">{e.label}</span>
                      <span className={`${MONO} text-[var(--ink)] shrink-0`}>{e.value}</span>
                    </div>
                    <div className="text-[10px] text-[var(--faint)] mt-0.5">
                      {e.tabLabel} · {e.section}
                    </div>
                  </button>
                ))
              )}
            </div>
          )}
        </div>
      </div>

      {/* Fixed-height content area so switching tabs never reflows the page. */}
      <div className="flex-1 min-h-0">
        {tab === "view" && <AircraftView run={run} onLayoutChange={onLayoutChange} />}
        {tab === "overview" && (
          <div className="h-full overflow-y-auto">
            <OverviewTab run={run} specTab={specTab} />
          </div>
        )}
        {tab === "power" && (
          <div className="h-full overflow-y-auto">
            <PowerTab run={run} specTab={specTab} />
          </div>
        )}
        {tab === "shapes" && (
          <div className="h-full overflow-y-auto space-y-10">
            <Planform run={run} />
            <Airfoils run={run} />
          </div>
        )}
        {specTab && tab !== "overview" && tab !== "power" && (
          <div className="h-full overflow-y-auto space-y-10">
            <SpecSections tab={specTab} columns={tab === "geometry" ? 2 : 3} />
            {tab === "geometry" && <MassCompare run={run} />}
          </div>
        )}
      </div>
    </div>
  );
}

function SearchIcon() {
  return (
    <svg
      className="absolute left-2 top-1/2 -translate-y-1/2 text-[var(--faint)] pointer-events-none"
      width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"
    >
      <circle cx="11" cy="11" r="7" />
      <path d="m21 21-4.3-4.3" />
    </svg>
  );
}

function SpecSections({ tab, columns = 3 }: { tab?: SpecTab; columns?: 1 | 2 | 3 }) {
  if (!tab) return null;
  const gridCls =
    columns === 1
      ? "grid-cols-1"
      : columns === 2
        ? "grid-cols-1 md:grid-cols-2"
        : "grid-cols-1 md:grid-cols-2 xl:grid-cols-3";
  return (
    <div className={`grid ${gridCls} gap-x-10 gap-y-6`}>
      {tab.sections.filter((s) => s.rows.length > 0).map((s) => (
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
  );
}
