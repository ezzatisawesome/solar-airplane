import { useCallback, useEffect, useMemo, useState } from "react";
import { useParams } from "react-router-dom";
import { fetchRun } from "../lib/data";
import type { ComponentMP, RunState, SpecTab } from "../lib/types";
import type { Layout } from "../lib/mass";
import AircraftView from "../components/AircraftView";
import { BatteryChart } from "../components/Charts";
import MassCompare from "../components/MassCompare";
import Airfoils from "../components/Airfoils";
import Planform from "../components/Planform";
import Field from "../components/Field";
import { setRunState } from "../lib/store";
import { LABEL, TAB, TAB_ACTIVE } from "../lib/ui";

export default function RunDetail() {
  const { id = "" } = useParams();
  const [run, setRun] = useState<RunState | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [tab, setTab] = useState("view");

  useEffect(() => {
    let alive = true;
    setRun(null);
    setError(null);
    setTab("view");
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

  if (error) return <div className="p-6 text-[var(--err)]">{error}</div>;
  if (!run) return <div className="p-6 text-[var(--muted)]">Loading {id}…</div>;

  const specTab = run.specs.find((s) => s.id === tab);

  return (
    <div className="h-full flex flex-col min-h-0 px-4 pb-4 pt-3 gap-3">
      <nav className="flex flex-wrap gap-1 border-b border-[var(--border)] pb-2">
        {tabs.map((t) => (
          <button key={t.id} onClick={() => setTab(t.id)} className={`${TAB} ${tab === t.id ? TAB_ACTIVE : ""}`}>
            {t.label}
          </button>
        ))}
      </nav>

      {/* Fixed-height content area so switching tabs never reflows the page. */}
      <div className="flex-1 min-h-0">
        {tab === "view" && <AircraftView run={run} onLayoutChange={onLayoutChange} />}
        {tab === "shapes" && (
          <div className="h-full overflow-y-auto space-y-10">
            <Planform run={run} />
            <Airfoils run={run} />
          </div>
        )}
        {specTab && (
          <div className="h-full overflow-y-auto space-y-10">
            <SpecSections tab={specTab} columns={tab === "geometry" ? 2 : 3} />
            {tab === "geometry" && <MassCompare run={run} />}
            {tab === "power" && <BatteryChart run={run} />}
          </div>
        )}
      </div>
    </div>
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
