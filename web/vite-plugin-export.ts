import { execFileSync } from "node:child_process";
import { existsSync, writeFileSync } from "node:fs";
import path from "node:path";
import type { Plugin } from "vite";

// Vite runs with cwd = web/, so the repo root is one level up.
const ROOT = path.resolve(process.cwd(), "..");
const PY = existsSync(path.join(ROOT, ".venv/bin/python"))
  ? path.join(ROOT, ".venv/bin/python")
  : "python3";

function runExport(): void {
  try {
    execFileSync(PY, ["export_site.py"], { cwd: ROOT, stdio: "inherit" });
  } catch {
    console.warn("\n[export] export_site.py failed — run it manually to refresh run data.\n");
  }
}

/**
 * Runs `export_site.py` automatically so run data is always fresh:
 *  - once when the dev server or a build starts, and
 *  - again (debounced) whenever a run's soln/state/builds JSON changes under output/,
 *    then triggers a full page reload.
 * This removes the manual "export, then serve" two-step.
 */
export function exportData(): Plugin {
  let timer: ReturnType<typeof setTimeout> | undefined;
  return {
    name: "solar-export-data",
    buildStart() {
      runExport();
    },
    configureServer(server) {
      // Dev-only endpoint: persist as-built builds to output/<run>/builds.json, then re-export.
      server.middlewares.use((req, res, next) => {
        if (req.method !== "POST" || !req.url?.startsWith("/api/builds/")) return next();
        const runId = decodeURIComponent(req.url.slice("/api/builds/".length).split("?")[0]);
        if (!/^run_[a-z0-9]+$/i.test(runId) || !existsSync(path.join(ROOT, "output", runId))) {
          res.statusCode = 400;
          return res.end("invalid run");
        }
        let body = "";
        req.on("data", (c) => (body += c));
        req.on("end", () => {
          try {
            const data = JSON.parse(body); // {builds:[{name,masses}]}
            writeFileSync(path.join(ROOT, "output", runId, "builds.json"), JSON.stringify(data.builds ?? data, null, 2));
            runExport();
            res.setHeader("Content-Type", "application/json");
            res.end(JSON.stringify({ ok: true }));
          } catch (e) {
            res.statusCode = 500;
            res.end(String(e));
          }
        });
      });

      server.watcher.add(path.join(ROOT, "output"));
      const onChange = (file: string) => {
        if (!/[\\/]output[\\/].*(soln|state|builds)\.json$/.test(file)) return;
        clearTimeout(timer);
        timer = setTimeout(() => {
          runExport();
          server.ws.send({ type: "full-reload" });
        }, 300);
      };
      server.watcher.on("add", onChange);
      server.watcher.on("change", onChange);
      server.watcher.on("unlink", onChange);       // a run's json deleted -> re-export prunes it
      server.watcher.on("unlinkDir", (dir: string) => {
        if (/[\\/]output[\\/]run_[a-z0-9]+$/i.test(dir)) onChange(path.join(dir, "soln.json"));
      });
    },
  };
}
