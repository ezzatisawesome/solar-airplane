// Shared header/button recipes — plain Tailwind utilities (no .btn CSS class),
// matching AircraftSim so the two apps render identical controls. Colors read
// AircraftView's :root vars via arbitrary values; sizing is 28px tall / 6px
// radius / 12px text, same as AircraftSim's BTN.
const BASE =
  "inline-flex items-center justify-center gap-1.5 h-8 rounded-md text-[13px] transition-colors cursor-pointer disabled:opacity-40 disabled:cursor-not-allowed";
const NEUTRAL =
  "font-medium bg-[var(--card-2)] text-[var(--ink)] border border-[var(--border)] hover:brightness-110";
const ACCENT =
  "font-semibold bg-[var(--accent)] text-[var(--accent-fg)] border border-transparent hover:brightness-90";

export const BTN = `${BASE} px-3 ${NEUTRAL}`;
export const BTN_ACCENT = `${BASE} px-3 ${ACCENT}`;
export const BTN_ICON = `${BASE} w-8 ${NEUTRAL}`; // square icon-only button, no side padding

// Former .label / .mono / .tab CSS classes, now plain Tailwind utilities.
export const LABEL = "font-mono text-[10px] uppercase tracking-[0.12em] text-[var(--muted)]";
export const MONO = "font-mono tabular-nums";
export const TAB =
  "px-3.5 py-1.5 rounded-md text-[13px] text-[var(--muted)] cursor-pointer transition-colors hover:text-[var(--ink)] hover:bg-[var(--card-2)]";
export const TAB_ACTIVE = "text-[var(--ink)] bg-[var(--card-2)]";
