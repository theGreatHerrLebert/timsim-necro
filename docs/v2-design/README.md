# timsim v2 — design record

The documents that designed v2, kept as written. They are a **historical record, not a status page**:
each states the state of the world when it was written (several say "S0 complete, v1 untouched"), and
some decisions were later overruled by the code. For what exists *today*, read the repo `README.md`
and the top-level planning docs (`PORTING_ROADMAP.md`, `REALISM_PLAN.md`, `HYE_QUANT.md`,
`PHOSPHO_FLR.md`, `RAMP.md`).

Read in this order:

| document | what it covers |
|---|---|
| **`TIMSIM_V2_SPEC.md`** | The design spec (Revision 4). The axes (structure / quantity / design / measurement), the artifact + schema contract, and the reasoning behind the split. Where an earlier revision and the code disagreed, the code won — and the revision says so explicitly. |
| **`TIMSIM_V2_PLAN.md`** | The implementation plan (Revision 4) — the strangler staging, what shipped in the first stage, what it cost, and what building taught us that planning did not. |
| **`TIMSIM_V2_RENDER.md`** | The design for the **render**: the last and heaviest measurement stage, turning the instrument-independent feature space into a real Bruker `.d` in bounded memory (target: peak RSS < 1 GB), so cost scales with the elution window rather than the run length. |
| **`TIMSIM_V2_MS2.md`** | Scoping MS2: DIA-PASEF fragment frames gated by the diagonal quadrupole transmission. Grounded in a full read of the v1 render — the conclusion being that most of the physics already existed in Rust, and the new work was a scan-resolved sweep emitter plus wiring. |

## Reviews

Each review is an independent (codex) critique of the document it is named after, kept alongside it
because the corrections are the interesting part.

| review | the correction it forced |
|---|---|
| `TIMSIM_V2_RENDER.codex-review.md` | Format and instrument assumptions checked against primary Bruker TDF documentation — the decisive risks being the on-disk layout and simulator realism. |
| `TIMSIM_V2_RENDER.codex-review2.md` | The scaling experiment establishes independence from *run duration*, not a universal few-MB bound. Working set is a function of local precursor arrival rate, elution-width distribution, mobility width, payload and acquisition geometry — not of run length. |
| `TIMSIM_V2_MS2.codex-review.md` | A precursor need not fall in exactly one window group (overlapping diagonals, soft transmission edges); co-isolation is expected rather than a defect; and the first version must be labelled explicitly as a b/y-only, complete-fragmentation idealization. |
