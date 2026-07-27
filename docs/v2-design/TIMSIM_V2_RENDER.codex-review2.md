§7 is a useful prototype result, but it validates a narrower claim than the text currently implies: with a stationary synthetic arrival process, fixed elution-window distribution, and fixed per-precursor payload, the sweep’s active set and one-frame accumulation buffer do not grow with run duration. That is exactly the expected sweep-line property, and the 30k/3k-frame vs 60k/6k-frame comparison is a sound first scaling experiment.

It is not circular, but it is conditional. “Fixed elution density” deliberately holds the quantity that determines active-set size constant; that establishes independence from total duration, not a universal few-MB bound. The plan should say this plainly:

`working set = f(local precursor arrival rate, elution-width distribution, mobility width, ion/fragment payload, acquisition geometry, and peak density)`, not `f(run length)`.

The missing breaking regimes are important:

- A concentrated RT/enrichment burst changes local arrival rate and can make active-set and frame-buffer peaks much larger.
- Broad or long-tailed elution windows increase overlap even at identical average density.
- Gradient mapping can compress many `rt_index` values into a short time region, especially at gradient edges.
- “Density grows with run length” is a different workload model: if a longer run is also loaded with proportionally more material per second, flatness is neither expected nor required.
- MS2 can be materially worse than MS1 because fragment multiplicity and co-isolation, rather than precursor count alone, determine `(scan,tof)` occupancy.

The stated stress fixture is therefore not optional follow-up; it is the test that turns the result into a capacity claim. Report p50/p95/p99 and maximum local arrival rate, active precursors, unique bins, and RSS—not only averages or one peak. Run a matrix varying local density, elution width/tails, mobility width, isotope depth, fragments per precursor, and cap behavior.

The conservation result is partly meaningful, but the prose overstates it. Comparing emitted mass to an expectation constructed from the same truncated/clamped integration bounds verifies that the implementation consistently accounts for those bounds; it can catch some active-set lifetime or accumulation mistakes. It does not independently establish “exactly once” if the expected calculation shares the same frame/scan partitioning, clamping, index mapping, or loop structure. A duplicated frame interval in both paths can pass.

A stronger MS1 oracle should be independent in both decomposition and implementation:

- For a tiny fixture, compute the expected total analytically from the continuous RT/mobility/isotope distributions and independently derived event boundaries, then compare with emitted bins.
- Independently enumerate each precursor/event/scan contribution into a reference tensor, using a deliberately simple representation and no production active-set, buffer, or frame-index code. Compare every bin, not only totals.
- Add metamorphic tests: split one event exposure into two adjacent events and require identical integrated signal; permute precursor input order; render overlapping chunks and compare their union to a single render; duplicate one precursor and require exactly 2× only its expected bins.
- Track per-precursor, per-charge, per-event emitted mass with stable identity tags in test builds. This directly catches enter/leave off-by-one, one-frame overlap, and leakage between scans/frames.

For MS2, equivalent micro-fixtures should independently calculate diagonal `Q(scan,mz)`, including soft-edge points, and exercise all three transmission modes. A total-intensity check is insufficient.

The proposed fixes are directionally right. Column projection plus a sequential, `frame_start`-sorted artifact is the right remedy for the observed RSS problem; avoid an in-memory RT map entirely. Co-sort or denormalise fragment payloads with the render artifact (or partition them identically) is likely necessary to avoid replacing RAM pressure with random I/O stalls.

A dense per-scan TOF buffer is also the right hot-path structure, provided it is not naively cleared over the full TOF range every frame. Use touched-TOF lists/epochs and deterministic scan-major, TOF-sorted emission. Its fixed allocation should be included in RSS reporting. This aligns naturally with TDF’s per-scan sorting/delta coding. It does not conflict with diagonal transmission—the latter determines which contributions are added—nor DDA, whose event stream merely changes events online. It does make the buffer-cap policy more consequential: dense storage bounds allocation, but not emitted peak count or compression time.

The throughput finding changes sequencing. “Prove the win, then port” was right when the main uncertainty was memory architecture. Memory has now been provisionally proven; throughput is the new load-bearing uncertainty. Before committing heavily to the rustdf writer, build a writer-free production-shaped benchmark that includes:

- projected, sequential pre-sorted input;
- dense/touched-bin accumulation;
- realistic MS2 fragments and diagonal transmission;
- sorting/dedup and the actual TDF block compression encoder, writing to a sink or temp file.

That benchmark should establish both a throughput target and an end-to-end RSS target. Otherwise the TDF port risks validating append-only format correctness while concealing an hours-long renderer.

Finally, §7 presently tests only simplified MS1 isotope emission. It should explicitly test same-precursor contributions across multiple active frames, which is where lifetime bugs live; multiple co-eluting precursors colliding into identical `(scan,tof)` bins; MS2 fragment fan-out; edge transmission; reference-frame lockstep; bounded compression backpressure; and deterministic behavior under cap/spill. Those are the remaining tests most likely to invalidate the practical—not merely asymptotic—claim.
