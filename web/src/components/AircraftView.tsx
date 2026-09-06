import { useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { TransformControls } from "three/addons/controls/TransformControls.js";
import type { ComponentMP, RunState } from "../lib/types";
import { applyLayout, type Layout } from "../lib/mass";
import { colorFor, parseGeometry, reconstructGeometry, solarCells, type Geom } from "../lib/geometry";
import { BTN, LABEL, MONO } from "../lib/ui";

interface SceneApi {
  select: (name: string | null) => void;
  reset: () => void;
  setOpacity: (v: number) => void;
  setDimsVisible: (v: boolean) => void;
  setView: (name: "iso" | "top" | "side" | "front") => void;
  setLayer: (name: "structure" | "cells" | "cg", v: boolean) => void;
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
  const [showStructure, setShowStructure] = useState(true);
  const [showCells, setShowCells] = useState(true);
  const [showCG, setShowCG] = useState(true);

  const cbRef = useRef(onLayoutChange);
  cbRef.current = onLayoutChange;
  const opacityRef = useRef(opacity);
  opacityRef.current = opacity;

  // Lets the in-scene click handler drive React selection (toggles off if re-clicked).
  const selectRef = useRef<(name: string) => void>(() => {});
  selectRef.current = (name: string) => setSelected((p) => (p === name ? null : name));

  // Quarter-chord (aerodynamic reference) x, geometrically consistent with the drawn wing:
  // wing mesh LE sits at center - 0.4*rootChord, so c/4 = center - 0.4c + 0.25c = center - 0.15c.
  const x_c4 = useMemo(() => {
    const wingX = run.massProperties.components["Main wing"]?.x_cg
      ?? run.massProperties.positions["Main wing"]?.xyz[0] ?? 0;
    return wingX - 0.15 * parseGeometry(run.soln).cRoot;
  }, [run]);

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

  const maxMass = useMemo(() => Math.max(0, ...elements.map((e) => e.mass)), [elements]);

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
    const structureMeshes: THREE.Mesh[] = [];
    const makeMesh = (g: Geom, color: number): THREE.Mesh => {
      let geometry: THREE.BufferGeometry;
      if (g.type === "box") geometry = new THREE.BoxGeometry(g.Lx, g.Ly, g.Lz);
      else if (g.type === "cylinder") {
        geometry = new THREE.CylinderGeometry(g.radius, g.radius, g.length, 20);
        geometry.rotateZ(Math.PI / 2);
      } else if (g.type === "wing" && g.sections && g.sections.length >= 2) {
        // Lofted surface from spanwise sections -> shows taper AND polyhedral (z rise).
        // Mirror the half-span sections across the centerline, flat LE at x=0.
        const half = g.sections;
        const full = [...half.slice(1).reverse().map((s) => ({ ...s, y: -s.y })), ...half];
        const x0 = -0.4 * g.rootChord; // roughly center chord on the CG
        const pos: number[] = [];
        const idx: number[] = [];
        for (const s of full) {
          pos.push(x0, s.y, s.z);                 // LE
          pos.push(x0 + s.chord, s.y, s.z);       // TE
        }
        for (let i = 0; i < full.length - 1; i++) {
          const a = 2 * i, b = 2 * i + 1, c = 2 * (i + 1), d = 2 * (i + 1) + 1;
          idx.push(a, b, d, a, d, c); // two triangles per panel (LE_i,TE_i,TE_i+1 / LE_i,TE_i+1,LE_i+1)
        }
        geometry = new THREE.BufferGeometry();
        geometry.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
        geometry.setIndex(idx);
        geometry.computeVertexNormals();
      } else if (g.type === "wing") {
        // Tapered planform: flat LE at x=0, chord tapers from root to tip. Built in the
        // (x=chord, span) plane and extruded by thickness; fins are rotated span -> z.
        const { rootChord: cr, tipChord: ct, span: sp, thickness: th } = g;
        const shape = new THREE.Shape();
        if (g.axis === "y") {
          // horizontal surface: flat LE at x=0, symmetric double taper (tip chord at both ends)
          shape.moveTo(0, -sp / 2);
          shape.lineTo(0, sp / 2);
          shape.lineTo(ct, sp / 2);
          shape.lineTo(cr, 0);
          shape.lineTo(ct, -sp / 2);
        } else {
          // vertical fin: root (cr) at bottom span=-sp/2, tip (ct = hstab chord) at top.
          // Vertical trailing edge at x=cr; raked leading edge (not a square front).
          shape.moveTo(0, -sp / 2);       // LE root (bottom, front)
          shape.lineTo(cr, -sp / 2);      // TE root (bottom, back)
          shape.lineTo(cr, sp / 2);       // TE tip  (top, back)  -> vertical TE
          shape.lineTo(cr - ct, sp / 2);  // LE tip  (top, front) -> raked LE
        }
        shape.closePath();
        geometry = new THREE.ExtrudeGeometry(shape, { depth: th, bevelEnabled: false });
        geometry.translate(-0.4 * cr, 0, -th / 2); // roughly center chord + thickness on CG
        if (g.axis === "z") geometry.rotateX(Math.PI / 2); // span y -> z, root stays at bottom
      } else geometry = new THREE.SphereGeometry(g.radius, 16, 12);
      const mesh = new THREE.Mesh(geometry, new THREE.MeshPhongMaterial({ color, side: THREE.DoubleSide }));
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
      if (name === "Solar cells") continue; // drawn as an actual cell grid below
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
        structureMeshes.push(mesh);
      }
      baseEmissive(mesh);
      group.add(mesh);
      meshes[name] = mesh;
    }

    // Solar cells: draw the actual row-aware cell grid on the wing top surface (ported from
    // solar-uav-design's _cell_quads). Cells are positioned in the wing-local frame, so we
    // anchor them at the "Main wing" component center.
    const wingCenter = geom["Main wing"]?.center ?? [0, 0, 0];
    const hstabCenter = geom["Horizontal stabilizer"]?.center ?? [0, 0, 0];
    const { wing: wingCells, hstab: hstabCells, side: cellSide } = solarCells(run.soln);
    // One instanced mesh covers both surfaces; each cell is anchored at its surface's center.
    const placed: [number, number, number][] = [
      ...wingCells.map((c) => [wingCenter[0] + c[0], wingCenter[1] + c[1], wingCenter[2] + c[2]] as [number, number, number]),
      ...hstabCells.map((c) => [hstabCenter[0] + c[0], hstabCenter[1] + c[1], hstabCenter[2] + c[2]] as [number, number, number]),
    ];
    if (placed.length) {
      const cellGeo = new THREE.BoxGeometry(cellSide, cellSide, 0.003);
      const cellMat = new THREE.MeshPhongMaterial({ color: 0x1a3a6b, emissive: 0x0a1a3a, emissiveIntensity: 0.35 });
      const cells = new THREE.InstancedMesh(cellGeo, cellMat, placed.length);
      const m4 = new THREE.Matrix4();
      placed.forEach((p, i) => {
        m4.setPosition(p[0], p[1], p[2]);
        cells.setMatrixAt(i, m4);
      });
      cells.instanceMatrix.needsUpdate = true;
      cells.userData = { name: "Solar cells", draggable: false };
      group.add(cells);
      meshes["Solar cells"] = cells as unknown as THREE.Mesh;
    }

    const cg = new THREE.Mesh(
      new THREE.SphereGeometry(0.05, 16, 12),
      new THREE.MeshPhongMaterial({ color: 0xff3b30, emissive: 0xff3b30, emissiveIntensity: 0.5 }),
    );
    const t0 = mp.total;
    cg.position.set(t0.x_cg, t0.y_cg, t0.z_cg);
    group.add(cg);

    // Wing quarter-chord marker (blue) — the aerodynamic reference the CG is balanced against.
    const wingCenterC4 = geom["Main wing"]?.center ?? [0, 0, 0];
    const c4 = new THREE.Mesh(
      new THREE.SphereGeometry(0.03, 12, 10),
      new THREE.MeshPhongMaterial({ color: 0x8ab4f8, emissive: 0x8ab4f8, emissiveIntensity: 0.4 }),
    );
    c4.position.set(x_c4, 0, wingCenterC4[2]);
    group.add(c4);

    // ---- In-scene dimension annotations (move/rotate with the aircraft) ----
    const dimsGroup = new THREE.Group();
    dimsGroup.visible = dimsRef.current;
    group.add(dimsGroup);

    const {
      b, cRoot: c, boomLen: bl, boomY, fuseLen,
      hstabSpan: hs, hstabChord: htc, vstabSpan: vspan, vstabRootChord: vrc,
    } = parseGeometry(run.soln);
    const wingX = mp.positions["Main wing"]?.xyz[0] ?? 0;
    const hx = bl + vrc / 4 + 0.25 * htc; // hstab quarter-chord (matches the CG-anchored mesh)

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
    addDim([-fuseLen, -boomY, -0.1], [0, -boomY, -0.1], `${fuseLen.toFixed(2)} m`, [-fuseLen / 2, -boomY, -0.2], ["Fuselage L", "Fuselage R"]);

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

    // Click a part directly in the scene to select it (a small drag = orbit, so ignore those).
    const selRay = new THREE.Raycaster();
    let downAt: { x: number; y: number } | null = null;
    const onDown = (e: PointerEvent) => (downAt = { x: e.clientX, y: e.clientY });
    const onUp = (e: PointerEvent) => {
      if (!downAt || Math.hypot(e.clientX - downAt.x, e.clientY - downAt.y) > 6) return;
      if ((gizmo as unknown as { dragging: boolean }).dragging) return;
      const r = renderer.domElement.getBoundingClientRect();
      const pt = new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
      selRay.setFromCamera(pt, camera);
      const hits = selRay.intersectObjects(Object.values(meshes), false);
      if (hits[0]) selectRef.current((hits[0].object as THREE.Mesh).userData.name as string);
    };
    renderer.domElement.addEventListener("pointerdown", onDown);
    renderer.domElement.addEventListener("pointerup", onUp);

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
      setView: (name) => {
        const R = 3.4;
        const q = group.quaternion;
        const dir = (x: number, y: number, z: number) => new THREE.Vector3(x, y, z).applyQuaternion(q).normalize();
        const fore = dir(1, 0, 0), span = dir(0, 1, 0), up = dir(0, 0, 1);
        group.updateMatrixWorld();
        const target = group.localToWorld(cg.position.clone());
        let offset: THREE.Vector3;
        if (name === "top") offset = up.multiplyScalar(R);
        else if (name === "side") offset = span.multiplyScalar(R);
        else if (name === "front") offset = fore.multiplyScalar(R);
        else offset = fore.multiplyScalar(1.1).add(span.multiplyScalar(0.9)).add(up.multiplyScalar(0.9)).normalize().multiplyScalar(R);
        camera.position.copy(target.clone().add(offset));
        orbit.target.copy(target);
        orbit.update();
      },
      setLayer: (name, v) => {
        if (name === "structure") structureMeshes.forEach((m) => (m.visible = v));
        else if (name === "cells") { if (meshes["Solar cells"]) meshes["Solar cells"].visible = v; }
        else if (name === "cg") { cg.visible = v; c4.visible = v; }
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
      renderer.domElement.removeEventListener("pointerdown", onDown);
      renderer.domElement.removeEventListener("pointerup", onUp);
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

  useEffect(() => {
    apiRef.current?.setLayer("structure", showStructure);
  }, [showStructure]);
  useEffect(() => {
    apiRef.current?.setLayer("cells", showCells);
  }, [showCells]);
  useEffect(() => {
    apiRef.current?.setLayer("cg", showCG);
  }, [showCG]);

  return (
    <div className="flex h-full min-h-0">
      <div className="relative min-w-0 flex-1 bg-[var(--bg)]">
        <div ref={mountRef} className="absolute inset-0 overflow-hidden" />

        <div className="absolute top-2 left-2 w-44 sm:w-56 max-h-[calc(100%-1rem)] overflow-y-auto flex flex-col gap-2">
          {/* View controls */}
          <Panel title="View">
            <div className="flex gap-1 mb-2">
              {(["iso", "top", "side", "front"] as const).map((v) => (
                <button
                  key={v}
                  onClick={() => apiRef.current?.setView(v)}
                  className="flex-1 px-1 py-0.5 text-[11px] capitalize rounded-sm bg-[var(--card-2)] text-[var(--muted)] hover:text-[var(--ink)] hover:brightness-110 transition-colors"
                >
                  {v}
                </button>
              ))}
            </div>
            <Toggle label="Dimensions" checked={showDims} onChange={setShowDims} />
            <Toggle label="Structure" checked={showStructure} onChange={setShowStructure} />
            <Toggle label="Solar cells" checked={showCells} onChange={setShowCells} />
            <Toggle label="CG / c/4" checked={showCG} onChange={setShowCG} />
            <div className="flex gap-3 my-1.5 text-[10px] text-[var(--muted)]">
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full" style={{ background: "#ff3b30" }} />CG</span>
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full" style={{ background: "#8ab4f8" }} />c/4</span>
            </div>
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
          </Panel>

          {/* Per-component mass, drawn as a clickable horizontal bar table.
              Click an element to select it; movable ones (↔) attach the move gizmo. */}
          <Panel title="Components & Mass">
            <div className="flex items-center justify-between">
              <Row label="Total" value={`${total.mass.toFixed(3)} kg`} />
              <button className={BTN} disabled={!dirty} onClick={() => apiRef.current?.reset()}>
                Reset
              </button>
            </div>
            <div className="mt-1.5 flex flex-col gap-1">
              {elements.map((el) => (
                <button
                  key={el.name}
                  onClick={() => setSelected(el.name === selected ? null : el.name)}
                  className={`w-full text-left px-1 py-0.5 rounded-sm transition-colors ${
                    selected === el.name ? "bg-[var(--card-2)]" : "hover:bg-[var(--card-2)]"
                  }`}
                >
                  <div className="flex items-center justify-between gap-2 text-[11px]">
                    <span className="flex items-center gap-1 truncate text-[var(--muted)]">
                      <span className="w-2 h-2 rounded-sm shrink-0" style={{ background: el.color }} />
                      <span className="truncate">{el.name}</span>
                      {el.draggable && <span className="text-[var(--faint)] shrink-0" title="movable">↔</span>}
                    </span>
                    <span className={`${MONO} shrink-0`}>
                      {el.mass.toFixed(3)}
                      <span className="text-[var(--faint)]">
                        {" "}
                        {total.mass > 0 ? ((el.mass / total.mass) * 100).toFixed(1) : "0.0"}%
                      </span>
                    </span>
                  </div>
                  <div className="mt-0.5 h-1.5 rounded-sm bg-[var(--card-2)] overflow-hidden">
                    <div
                      className="h-full rounded-sm"
                      style={{ width: `${maxMass > 0 ? (el.mass / maxMass) * 100 : 0}%`, background: el.color }}
                    />
                  </div>
                </button>
              ))}
            </div>
            <div className="mt-2 text-[11px] text-[var(--muted)]">
              {selected
                ? elements.find((e) => e.name === selected)?.draggable
                  ? "Drag the arrows to move; CG updates live."
                  : "This element is fixed by geometry."
                : "Select an element to move it."}
            </div>
          </Panel>

          <Panel title="Mass Balance">
            <Row label="Total" value={`${total.mass.toFixed(3)} kg`} />
            <Row label="CG x" value={`${total.x_cg.toFixed(4)} m`} />
            <Row label="CG y" value={`${total.y_cg.toFixed(4)} m`} />
            <Row label="CG z" value={`${total.z_cg.toFixed(4)} m`} />
            <Row label="c/4 x" value={`${x_c4.toFixed(4)} m`} />
            <Row label="CG − c/4" value={`${(total.x_cg - x_c4).toFixed(4)} m`} />
            <div className={`${LABEL} mt-2 mb-1`}>Moments of Inertia (kg·m²)</div>
            <Row label="Ixx" value={total.Ixx.toFixed(3)} />
            <Row label="Iyy" value={total.Iyy.toFixed(3)} />
            <Row label="Izz" value={total.Izz.toFixed(3)} />
            <Row label="Ixy" value={total.Ixy.toFixed(3)} />
            <Row label="Ixz" value={total.Ixz.toFixed(3)} />
            <Row label="Iyz" value={total.Iyz.toFixed(3)} />
          </Panel>
        </div>
      </div>
    </div>
  );
}

/**
 * A collapsible overlay panel. The header doubles as the toggle; a rotating chevron
 * indicates open/closed state (replaces the old global "Hide panels" button).
 */
function Panel({ title, children }: { title: string; children: ReactNode }) {
  const [open, setOpen] = useState(true);
  return (
    <div className="bg-[var(--card)] border border-[var(--border)] rounded-md">
      <button
        onClick={() => setOpen((v) => !v)}
        className="w-full flex items-center justify-between px-3 py-2"
      >
        <span className={`${LABEL}`}>{title}</span>
        <span
          className={`text-[var(--muted)] text-xs transition-transform ${open ? "rotate-180" : ""}`}
        >
          ▾
        </span>
      </button>
      {open && <div className="px-3 pb-2">{children}</div>}
    </div>
  );
}

function Toggle({ label, checked, onChange }: { label: string; checked: boolean; onChange: (v: boolean) => void }) {
  return (
    <label className="flex items-center justify-between cursor-pointer mb-1.5">
      <span className={`${LABEL}`}>{label}</span>
      <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} className="accent-white" />
    </label>
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
