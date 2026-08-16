import { useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { TransformControls } from "three/addons/controls/TransformControls.js";
import type { ComponentMP, RunState } from "../lib/types";
import { applyLayout, type Layout } from "../lib/mass";
import { colorFor, reconstructGeometry, type Geom } from "../lib/geometry";
import { BTN, LABEL, MONO } from "../lib/ui";

interface SceneApi {
  select: (name: string | null) => void;
  reset: () => void;
  setOpacity: (v: number) => void;
  setDimsVisible: (v: boolean) => void;
}

/**
 * The "AircraftView" tab: a 3D scene on the left and a selectable element list on the right.
 * Selecting a movable element attaches translate-arrows (a gizmo) so it can be nudged along each
 * axis — no free-dragging. Mass balance recomputes in-browser on every move.
 */
export default function AircraftView({
  run,
  onLayoutChange,
}: {
  run: RunState;
  onLayoutChange?: (layout: Layout, total: ComponentMP) => void;
}) {
  const mountRef = useRef<HTMLDivElement>(null);
  const apiRef = useRef<SceneApi | null>(null);
  const [total, setTotal] = useState<ComponentMP>(run.massProperties.total);
  const [dirty, setDirty] = useState(false);
  const [selected, setSelected] = useState<string | null>(null);
  const [opacity, setOpacity] = useState(0.28);

  const cbRef = useRef(onLayoutChange);
  cbRef.current = onLayoutChange;
  const opacityRef = useRef(opacity);
  opacityRef.current = opacity;

  // Overlay panels (view controls + mass balance) start collapsed on small screens so they
  // don't cover the 3D canvas; open by default on wider viewports.
  const [panelsOpen, setPanelsOpen] = useState(() => window.innerWidth >= 768);

  const [showDims, setShowDims] = useState(true);
  const dimsRef = useRef(showDims);
  dimsRef.current = showDims;

  // Element list grouped by category, sorted, with mass + draggability.
  const elements = useMemo(() => {
    const mp = run.massProperties;
    return Object.entries(mp.positions)
      .map(([name, p]) => ({
        name,
        category: p.category,
        draggable: p.draggable,
        mass: mp.components[name]?.mass ?? 0,
        color: `#${colorFor(name, p.category).toString(16).padStart(6, "0")}`,
      }))
      .sort((a, b) => a.category.localeCompare(b.category) || b.mass - a.mass);
  }, [run]);

  useEffect(() => {
    const mount = mountRef.current;
    if (!mount) return;
    setDirty(false);
    setSelected(null);
    setTotal(run.massProperties.total);

    const mp = run.massProperties;
    const geom = reconstructGeometry(run.soln, mp, run.constants?.fuselage_radius ?? 0.015);
    const layout: Layout = {};

    const scene = new THREE.Scene();
    // Transparent clear so the themed container background shows through (light/dark aware).
    // Rotate the aircraft frame (x fore, y span, z up) into a nice view; gizmo lives outside it.
    const group = new THREE.Group();
    group.rotation.z = Math.PI;
    group.rotation.x = -Math.PI / 2;
    scene.add(group);

    const { clientWidth: w, clientHeight: h } = mount;
    const camera = new THREE.PerspectiveCamera(50, w / Math.max(h, 1), 0.01, 100);
    camera.position.set(3, 2.5, 3);

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setClearColor(0x000000, 0);
    renderer.setSize(w, h);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    mount.appendChild(renderer.domElement);

    const orbit = new OrbitControls(camera, renderer.domElement);
    orbit.enableDamping = true;

    scene.add(new THREE.AmbientLight(0xffffff, 0.75));
    const dir = new THREE.DirectionalLight(0xffffff, 0.85);
    dir.position.set(5, 6, 4);
    scene.add(dir);
    const grid = new THREE.GridHelper(6, 24, 0x888888, 0x555555);
    grid.rotation.x = -Math.PI / 2;
    group.add(grid);

    // Meshes
    const meshes: Record<string, THREE.Mesh> = {};
    const structuralMats: THREE.MeshPhongMaterial[] = [];
    const makeMesh = (g: Geom, color: number): THREE.Mesh => {
      let geometry: THREE.BufferGeometry;
      if (g.type === "box") geometry = new THREE.BoxGeometry(g.Lx, g.Ly, g.Lz);
      else if (g.type === "cylinder") {
        geometry = new THREE.CylinderGeometry(g.radius, g.radius, g.length, 20);
        geometry.rotateZ(Math.PI / 2);
      } else geometry = new THREE.SphereGeometry(g.radius, 16, 12);
      const mesh = new THREE.Mesh(geometry, new THREE.MeshPhongMaterial({ color }));
      mesh.position.set(g.center[0], g.center[1], g.center[2]);
      return mesh;
    };

    const setEmissive = (mesh: THREE.Mesh, hex: number, intensity: number) => {
      const m = mesh.material as THREE.MeshPhongMaterial;
      m.emissive.setHex(hex);
      m.emissiveIntensity = intensity;
    };
    const baseEmissive = (mesh: THREE.Mesh) =>
      setEmissive(mesh, mesh.userData.draggable ? 0x2ecc71 : 0x000000, mesh.userData.draggable ? 0.16 : 0);

    for (const [name, g] of Object.entries(geom)) {
      const info = mp.positions[name];
      const mesh = makeMesh(g, colorFor(name, info?.category));
      mesh.userData = { name, draggable: info?.draggable === true };
      // Make the airframe structure translucent so the internal masses are visible through it.
      if (!mesh.userData.draggable) {
        const m = mesh.material as THREE.MeshPhongMaterial;
        m.transparent = true;
        m.opacity = opacityRef.current;
        m.depthWrite = opacityRef.current >= 1;
        structuralMats.push(m);
      }
      baseEmissive(mesh);
      group.add(mesh);
      meshes[name] = mesh;
    }

    const cg = new THREE.Mesh(
      new THREE.SphereGeometry(0.05, 16, 12),
      new THREE.MeshPhongMaterial({ color: 0xff3b30, emissive: 0xff3b30, emissiveIntensity: 0.5 }),
    );
    const t0 = mp.total;
    cg.position.set(t0.x_cg, t0.y_cg, t0.z_cg);
    group.add(cg);

    // ---- In-scene dimension annotations (move/rotate with the aircraft) ----
    const dimsGroup = new THREE.Group();
    dimsGroup.visible = dimsRef.current;
    group.add(dimsGroup);

    const so = run.soln;
    const mw = so["Main Wing"] || {};
    const gg = so.Geometry || {};
    const hh = so.HStab || {};
    const vv = so["V Stab"] || {};
    const b = mw.wingspan || 0;
    const c = mw.chordlen || 0;
    const bl = gg.boom_length || 0;
    const boomY = gg.boom_y || 0;
    const hs = hh.hstab_span || 0;
    const vspan = vv.vstab_span || 0;
    const vrc = vv.vstab_root_chord || 0;
    const wingX = mp.positions["Main wing"]?.xyz[0] ?? 0;
    const hx = bl + vrc / 4;

    // Each dimension bundles a line + number label, linked to the component(s) it measures.
    const dimObjects: { line: THREE.Line; label: THREE.Sprite; targets: string[] }[] = [];
    const addDim = (
      a: [number, number, number],
      bb: [number, number, number],
      text: string,
      labelPos: [number, number, number],
      targets: string[],
    ) => {
      const geo = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(...a), new THREE.Vector3(...bb)]);
      const line = new THREE.Line(geo, new THREE.LineBasicMaterial({ color: 0x9aa0a6 }));
      dimsGroup.add(line);

      const cv = document.createElement("canvas");
      cv.width = 256;
      cv.height = 64;
      const ctx = cv.getContext("2d")!;
      ctx.font = "30px monospace";
      ctx.fillStyle = "#ffffff";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(text, 128, 34);
      const label = new THREE.Sprite(
        new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(cv), depthTest: false, transparent: true }),
      );
      label.position.set(...labelPos);
      label.scale.set(0.6, 0.15, 1);
      dimsGroup.add(label);

      dimObjects.push({ line, label, targets });
    };

    addDim([wingX, -b / 2, 0.14], [wingX, b / 2, 0.14], `${b.toFixed(2)} m`, [wingX, 0, 0.24], ["Main wing"]);
    addDim([wingX - c / 2, b / 2 + 0.12, 0], [wingX + c / 2, b / 2 + 0.12, 0], `${c.toFixed(2)} m`, [wingX, b / 2 + 0.22, 0], ["Main wing"]);
    addDim([0, boomY, 0.06], [bl, boomY, 0.06], `${bl.toFixed(2)} m`, [bl / 2, boomY, 0.14], ["Boom R", "Boom L"]);
    addDim([bl, boomY, 0], [bl, boomY, vspan], `${vspan.toFixed(2)} m`, [bl, boomY, vspan + 0.1], ["Vertical stabilizer R", "Vertical stabilizer L"]);
    addDim([hx, -hs / 2, vspan], [hx, hs / 2, vspan], `${hs.toFixed(2)} m`, [hx, 0, vspan + 0.14], ["Horizontal stabilizer"]);
    addDim([-0.5, -boomY, -0.1], [0, -boomY, -0.1], "0.50 m", [-0.25, -boomY, -0.2], ["Fuselage L", "Fuselage R"]);

    // Translate gizmo (arrows). Lives in the (unrotated) scene; local space aligns to aircraft axes.
    const gizmo = new TransformControls(camera, renderer.domElement);
    gizmo.setMode("translate");
    gizmo.setSpace("local");
    gizmo.setSize(0.8);
    gizmo.addEventListener("dragging-changed", (e) => {
      orbit.enabled = !(e as unknown as { value: boolean }).value;
    });
    gizmo.addEventListener("objectChange", () => {
      const obj = gizmo.object as THREE.Mesh | undefined;
      if (!obj) return;
      const name = obj.userData.name as string;
      layout[name] = [+obj.position.x.toFixed(4), +obj.position.y.toFixed(4), +obj.position.z.toFixed(4)];
      const res = applyLayout(mp, layout);
      cg.position.set(res.total.x_cg, res.total.y_cg, res.total.z_cg);
      setTotal(res.total);
      setDirty(Object.keys(layout).length > 0);
      cbRef.current?.({ ...layout }, res.total);
    });
    scene.add(gizmo);

    // Hover a dimension line → highlight it and the component edge(s) it measures.
    const dimRay = new THREE.Raycaster();
    dimRay.params.Line.threshold = 0.05;
    const dimPointer = new THREE.Vector2();
    let hoveredDim: (typeof dimObjects)[number] | null = null;
    const clearDimHL = () => {
      if (!hoveredDim) return;
      (hoveredDim.line.material as THREE.LineBasicMaterial).color.setHex(0x9aa0a6);
      (hoveredDim.label.material as THREE.SpriteMaterial).color.setHex(0xffffff);
      hoveredDim.targets.forEach((n) => meshes[n] && baseEmissive(meshes[n]));
      hoveredDim = null;
    };
    const onDimMove = (e: PointerEvent) => {
      if (!dimsGroup.visible) return;
      const r = renderer.domElement.getBoundingClientRect();
      dimPointer.x = ((e.clientX - r.left) / r.width) * 2 - 1;
      dimPointer.y = -((e.clientY - r.top) / r.height) * 2 + 1;
      dimRay.setFromCamera(dimPointer, camera);
      // Hover either the line or its number label.
      const hit = dimRay.intersectObjects(dimObjects.flatMap((d) => [d.line, d.label]))[0];
      const dim = hit ? dimObjects.find((d) => d.line === hit.object || d.label === hit.object) ?? null : null;
      if (dim === hoveredDim) return;
      clearDimHL();
      if (dim) {
        hoveredDim = dim;
        (dim.line.material as THREE.LineBasicMaterial).color.setHex(0xf59e0b);
        (dim.label.material as THREE.SpriteMaterial).color.setHex(0xf59e0b);
        dim.targets.forEach((n) => meshes[n] && setEmissive(meshes[n], 0xf59e0b, 0.6));
      }
    };
    renderer.domElement.addEventListener("pointermove", onDimMove);

    let current: THREE.Mesh | null = null;
    apiRef.current = {
      select: (name) => {
        if (current) baseEmissive(current);
        gizmo.detach();
        current = name ? meshes[name] ?? null : null;
        if (!current) return;
        setEmissive(current, 0xffffff, 0.6);
        if (current.userData.draggable) gizmo.attach(current);
      },
      setOpacity: (v) => {
        for (const m of structuralMats) {
          m.opacity = v;
          m.transparent = v < 1;
          m.depthWrite = v >= 1;
          m.needsUpdate = true;
        }
      },
      setDimsVisible: (v) => {
        dimsGroup.visible = v;
      },
      reset: () => {
        for (const k of Object.keys(layout)) delete layout[k];
        const res = applyLayout(mp, layout);
        for (const [name, p] of Object.entries(res.positions)) {
          if (meshes[name]) meshes[name].position.set(p.xyz[0], p.xyz[1], p.xyz[2]);
        }
        cg.position.set(res.total.x_cg, res.total.y_cg, res.total.z_cg);
        setTotal(res.total);
        setDirty(false);
        cbRef.current?.({}, res.total);
      },
    };

    // Dev-only handle for tooling/tests.
    if (import.meta.env.DEV) {
      (window as unknown as { __av?: unknown }).__av = { scene, camera, meshes, cg, gizmo, canvas: renderer.domElement };
    }

    let raf = 0;
    const animate = () => {
      raf = requestAnimationFrame(animate);
      orbit.update();
      renderer.render(scene, camera);
    };
    animate();

    const onResize = () => {
      const cw = mount.clientWidth;
      const ch = mount.clientHeight;
      camera.aspect = cw / Math.max(ch, 1);
      camera.updateProjectionMatrix();
      renderer.setSize(cw, ch);
    };
    const ro = new ResizeObserver(onResize);
    ro.observe(mount);

    return () => {
      cancelAnimationFrame(raf);
      ro.disconnect();
      apiRef.current = null;
      if (import.meta.env.DEV) delete (window as unknown as { __av?: unknown }).__av;
      gizmo.detach();
      gizmo.dispose();
      orbit.dispose();
      scene.traverse((o) => {
        const m = o as THREE.Mesh;
        if (m.geometry) m.geometry.dispose();
        const mat = m.material as (THREE.Material & { map?: THREE.Texture }) | THREE.Material[] | undefined;
        if (Array.isArray(mat)) mat.forEach((x) => x.dispose());
        else if (mat) {
          mat.map?.dispose(); // sprite label textures
          mat.dispose();
        }
      });
      grid.dispose();
      renderer.domElement.removeEventListener("pointermove", onDimMove);
      renderer.forceContextLoss();
      renderer.dispose();
      mount.removeChild(renderer.domElement);
    };
  }, [run]);

  // Drive selection from the element list.
  useEffect(() => {
    apiRef.current?.select(selected);
  }, [selected]);

  useEffect(() => {
    apiRef.current?.setDimsVisible(showDims);
  }, [showDims]);

  return (
    <div className="flex flex-col md:flex-row h-full min-h-0">
      <aside className="order-2 md:order-1 w-full md:w-64 flex flex-col min-h-0 flex-1 md:flex-none border-t md:border-t-0 md:border-r border-[var(--border)]">
        <div className="flex items-center justify-between px-3 py-2">
          <span className={`${LABEL}`}>Elements</span>
          <button className={BTN} disabled={!dirty} onClick={() => apiRef.current?.reset()}>
            Reset
          </button>
        </div>
        <ul className="flex-1 overflow-y-auto">
          {elements.map((el) => (
            <li key={el.name}>
              <button
                onClick={() => setSelected(el.name === selected ? null : el.name)}
                className={`w-full flex items-center gap-2 px-3 py-1.5 text-left text-sm transition-colors ${
                  selected === el.name ? "bg-[var(--card-2)] text-[var(--ink)]" : "hover:bg-[var(--card-2)]"
                }`}
              >
                <span className="w-2.5 h-2.5 rounded-sm shrink-0" style={{ background: el.color }} />
                <span className="truncate flex-1">{el.name}</span>
                {el.draggable && <span className="text-[var(--faint)] text-xs" title="movable">↔</span>}
                <span className={`${MONO} text-xs text-[var(--muted)]`}>{el.mass.toFixed(2)}</span>
              </button>
            </li>
          ))}
        </ul>
        <div className="px-3 py-2 border-t border-[var(--border)] text-xs text-[var(--muted)]">
          {selected
            ? elements.find((e) => e.name === selected)?.draggable
              ? "Drag the arrows to move; CG updates live."
              : "This element is fixed by geometry."
            : "Select an element to move it."}
        </div>
      </aside>

      <div className="order-1 md:order-2 relative min-w-0 bg-[var(--bg)] h-[55vh] shrink-0 md:h-auto md:flex-1">
        <div ref={mountRef} className="absolute inset-0 overflow-hidden" />

        {/* Toggle for the overlay panels — keeps the canvas usable on small screens. */}
        <button
          className={`${BTN} absolute top-2 right-2 z-10`}
          onClick={() => setPanelsOpen((v) => !v)}
        >
          {panelsOpen ? "Hide panels" : "Panels"}
        </button>

        <div
          className={`absolute top-11 right-2 w-44 sm:w-56 max-h-[calc(100%-3.5rem)] overflow-y-auto flex-col gap-2 ${panelsOpen ? "flex" : "hidden"}`}
        >
          {/* View controls, stacked above the mass balance panel */}
          <div className="bg-[var(--card)] border border-[var(--border)] rounded-md px-3 py-2">
            <label className="flex items-center justify-between cursor-pointer mb-2">
              <span className={`${LABEL}`}>Dimensions</span>
              <input type="checkbox" checked={showDims} onChange={(e) => setShowDims(e.target.checked)} className="accent-white" />
            </label>
            <div className="flex items-center justify-between mb-1">
              <span className={`${LABEL}`}>Structure opacity</span>
              <span className={`${MONO} text-xs text-[var(--muted)]`}>{Math.round(opacity * 100)}%</span>
            </div>
            <input
              type="range"
              min={0.05}
              max={1}
              step={0.01}
              value={opacity}
              onChange={(e) => {
                const v = +e.target.value;
                setOpacity(v);
                apiRef.current?.setOpacity(v);
              }}
              className="w-full accent-white"
            />
          </div>

          <div className="bg-[var(--card)] border border-[var(--border)] rounded-md px-3 py-2">
            <div className={`${LABEL} mb-1`}>Mass Balance</div>
          <Row label="Total" value={`${total.mass.toFixed(3)} kg`} />
          <Row label="CG x" value={`${total.x_cg.toFixed(4)} m`} />
          <Row label="CG y" value={`${total.y_cg.toFixed(4)} m`} />
          <Row label="CG z" value={`${total.z_cg.toFixed(4)} m`} />
          <div className={`${LABEL} mt-2 mb-1`}>Moments of Inertia (kg·m²)</div>
          <Row label="Ixx" value={total.Ixx.toFixed(3)} />
          <Row label="Iyy" value={total.Iyy.toFixed(3)} />
          <Row label="Izz" value={total.Izz.toFixed(3)} />
          <Row label="Ixy" value={total.Ixy.toFixed(3)} />
          <Row label="Ixz" value={total.Ixz.toFixed(3)} />
          <Row label="Iyz" value={total.Iyz.toFixed(3)} />
          </div>
        </div>
      </div>
    </div>
  );
}

function Row({ label, value }: { label: string; value: string }) {
  return (
    <div className="flex justify-between gap-6 py-0.5 text-xs">
      <span className="text-[var(--muted)]">{label}</span>
      <b className={`${MONO}`}>{value}</b>
    </div>
  );
}
