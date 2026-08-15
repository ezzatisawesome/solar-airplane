import { useSyncExternalStore } from "react";
import type { ComponentMP, RunState } from "./types";
import type { Layout } from "./mass";

// Tiny shared store so the global header (Download/Simulate) can see the currently-open run and
// its live dragged layout, without prop-drilling through the router.
interface AppState {
  run: RunState | null;
  layout: Layout;
  total: ComponentMP | null;
}

let state: AppState = { run: null, layout: {}, total: null };
const subs = new Set<() => void>();

export function setRunState(patch: Partial<AppState>) {
  state = { ...state, ...patch };
  subs.forEach((f) => f());
}

export function useRunState(): AppState {
  return useSyncExternalStore(
    (cb) => {
      subs.add(cb);
      return () => subs.delete(cb);
    },
    () => state,
  );
}
