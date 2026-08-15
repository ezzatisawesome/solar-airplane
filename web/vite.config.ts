import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tailwindcss from "@tailwindcss/vite";
import { exportData } from "./vite-plugin-export";

// Fully static site: run data is auto-exported (export_site.py) into public/ on dev/build start
// and re-exported when a run changes. No API/backend in production. Relative base so the bundle
// works from any static host / subpath.
export default defineConfig({
  base: "./",
  plugins: [exportData(), react(), tailwindcss()],
  server: { port: 3002, strictPort: true },
  build: { outDir: "dist" },
});
