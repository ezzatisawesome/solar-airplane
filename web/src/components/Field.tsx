/** Shared label→value row used by spec sections, dimension lists, etc. */
import { MONO } from "../lib/ui";
export default function Field({ label, value }: { label: string; value: string }) {
  return (
    <div className="flex justify-between border-b border-[var(--border)]/50 py-1.5">
      <span className="text-[var(--muted)]">{label}</span>
      <span className={`${MONO} text-right`}>{value}</span>
    </div>
  );
}
