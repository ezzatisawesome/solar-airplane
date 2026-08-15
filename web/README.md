# Solar Airplane — Design Explorer (web)

A fully static Vite + React + TypeScript + Tailwind site that displays and compares aircraft runs
and as-built airframes/flights. No backend at runtime: run data is pre-exported to JSON and the
mass balance recomputes in the browser.

## Develop

```bash
cd web
npm install
npm run dev      # http://localhost:3002
```

`npm run dev` (and `npm run build`) **auto-run `export_site.py`** first, and the dev server
re-exports + hot-reloads whenever a run's `soln.json` / `state.json` / `builds.json` changes under
`../output/`. You never run the exporter by hand.

Manual refresh (rarely needed): `npm run export`.

## Build & publish (static, host anywhere)

```bash
npm run build    # -> web/dist/ (auto-exports data first)
npm run preview  # optional local preview of the production build
```

`web/dist/` is a self-contained static bundle — no Python/API at runtime (verified: it makes zero
`/api` calls). Deploy it to any static host:

- **Netlify / Vercel**: build command `npm run build`, publish directory `dist`.
- **GitHub Pages / S3 / nginx**: upload the contents of `dist/`.

Routing uses `HashRouter` and `base: "./"`, so it works from any subpath with no server rewrites.

## Data model

`export_site.py` (repo root) writes into `public/`:

- `manifest.json` — index of runs (id, label, date, headline metrics, build names).
- `runs/<id>.json` — full state: soln, mass properties (+ per-component `positions`), mass chart,
  energy/xflr data, specs HTML, `builds[]` (as-built masses), and `constants`.
- `runs/<id>.sensitivity.json` — large sensitivity data, loaded on demand.

## Structure

- `src/lib/mass.ts` — client-side mass/CG/inertia recompute (parity-tested vs Python to 1e-9).
- `src/lib/geometry.ts` — component geometry for the 3D view.
- `src/lib/compare.ts` — run/build/upload entity model for comparisons.
- `src/components/AircraftView.tsx` — Three.js 3D view + drag-to-retune point masses.
- `src/components/Charts.tsx` — Chart.js visualizations.
- `src/views/` — `Gallery`, `RunDetail`, `Compare`.
