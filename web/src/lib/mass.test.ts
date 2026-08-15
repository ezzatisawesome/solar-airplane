import { describe, expect, it } from "vitest";
import { applyLayout } from "./mass";
import type { MassProperties } from "./types";
import fixture from "./__fixtures__/run_31sn.json";

// Parity against the Python ground truth (process_mass.compute_mass_properties). The fixture holds
// the base mass properties plus, for several dragged layouts, Python's resulting total.
describe("mass.ts parity with Python", () => {
  const base = fixture.massProperties as unknown as MassProperties;
  const keys = ["mass", "x_cg", "y_cg", "z_cg", "Ixx", "Iyy", "Izz", "Ixy", "Iyz", "Ixz"] as const;

  fixture.cases.forEach((c: any, i: number) => {
    it(`case ${i} matches Python within 1e-9`, () => {
      const total = applyLayout(base, c.layout).total as any;
      for (const k of keys) {
        expect(Math.abs(total[k] - c.total[k])).toBeLessThan(1e-9);
      }
    });
  });
});
