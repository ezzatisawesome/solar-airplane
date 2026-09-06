import { useMemo } from "react";
import type { RunState, SpecTab } from "../lib/types";
import { energyMarch, sunGeometry } from "../lib/solar";
import Field from "./Field";
import { LABEL } from "../lib/ui";

const hod = (h: number) => {
  const hh = Math.floor(h);
  const mm = Math.round((h - hh) * 60);
  return `${String(hh).padStart(2, "0")}:${String(mm % 60).padStart(2, "0")}`;
};

/**
 * Overview tab: a mission-and-sun summary (design day, where it flies, how long the night is),
 * then the detailed spec sections.
 */
export default function OverviewTab({ run, specTab }: { run: RunState; specTab?: SpecTab }) {
  const e = useMemo(() => energyMarch(run.soln), [run]);
  const lat = run.soln.Mission?.operating_lat ?? 0;
  const day = run.soln.Mission?.mission_date ?? 172;
  const alt = run.soln.Mission?.operating_altitude ?? 0;
  const sun = useMemo(() => sunGeometry(lat, day), [lat, day]);

  return (
    <div className="space-y-10">
      <div className="grid grid-cols-1 md:grid-cols-2 gap-x-10 gap-y-6">
        <div>
          <h3 className={`${LABEL} mb-2`}>Mission & design day</h3>
          <dl className="text-sm">
            <Field label="Run ID" value={run.soln.Meta?.run_id ?? run.run} />
            <Field label="Design day" value={`${sun.dateLabel} (day ${Math.round(sun.dayOfYear)})`} />
            <Field label="Operating latitude" value={`${lat.toFixed(3)}°`} />
            <Field label="Operating altitude (MSL)" value={`${alt.toFixed(0)} m`} />
            <Field label="Solar declination" value={`${sun.declDeg.toFixed(1)}°`} />
          </dl>
        </div>
        <div>
          <h3 className={`${LABEL} mb-2`}>Sun & night</h3>
          <dl className="text-sm">
            <Field label="Daylight" value={`${sun.daylightHours.toFixed(2)} h`} />
            <Field label="Night" value={`${sun.nightHours.toFixed(2)} h`} />
            <Field label="Sunrise → sunset" value={sun.polarDay ? "midnight sun" : sun.polarNight ? "polar night" : `${hod(sun.sunriseHod)} → ${hod(sun.sunsetHod)}`} />
            <Field label="Solar noon elevation" value={`${sun.noonElevationDeg.toFixed(1)}°`} />
            {e.hasData && <Field label="Overnight discharge" value={`${e.nightEnergyWh.toFixed(0)} Wh (${(100 * e.depthOfDischarge).toFixed(0)}%)`} />}
          </dl>
        </div>
      </div>

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
