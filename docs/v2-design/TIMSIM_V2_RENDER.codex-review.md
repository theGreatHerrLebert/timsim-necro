I’ll assess the format and instrument assumptions against current primary documentation, since the TDF layout and simulator realism are the decisive risks.
web search: 
web search: Bruker TDF analysis.tdf_bin binary layout Frames TimsId ScanNumBegin ScanNumEnd ...
exec
/bin/bash -lc "rg -n \"TimsId|tdf_bin|CompressionType|NumScans|NumPeaks\" . -g '*.{rs,py,sql}' | head -160" in /scratch/timsim-demo/timsim-necro
 succeeded in 0ms:
./rustims/scripts/repair_timsim_midia_frames.py:19:Id / Time / MsMsType to its rightful slot. The binary tdf_bin is not
./rustims/scripts/repair_timsim_midia_frames.py:20:touched — only the metadata is wrong; TimsId still points to a valid
./rustims/scripts/repair_timsim_midia_frames.py:70:        "SELECT TimsId FROM Frames WHERE MsMsType = -1 ORDER BY TimsId"
./rustims/scripts/repair_timsim_midia_frames.py:113:    # 3. Pair each broken row (sorted by TimsId = chronological binary
./rustims/scripts/repair_timsim_midia_frames.py:126:    #    update by TimsId, not by Id.
./rustims/scripts/repair_timsim_midia_frames.py:145:        print(f"  TimsId={u[3]:>12}  ->  Id={u[0]:>6}  Time={u[1]:.4f}  MsMsType={u[2]}")
./rustims/scripts/repair_timsim_midia_frames.py:151:    # 6. Apply in a single transaction; we update by TimsId since it's the
./rustims/scripts/repair_timsim_midia_frames.py:154:        "UPDATE Frames SET Id = ?, Time = ?, MsMsType = ? WHERE TimsId = ?",
./rustims/scripts/parity_tof_writer.py:2:Bruker tdf_bin realdata for the same raw (scan, tof, intensity) frame input.
./rustims/rustdf/src/data/utility.rs:46:/// The Bruker `tdf_bin` layout requires this ordering: [`modify_tofs`] delta-
./rustims/rustdf/src/data/meta.rs:283:            "TimsCompressionType" => {
./rustims/rustdf/src/data/meta.rs:286:            "MaxNumPeaksPerScan" => {
./rustims/rustdf/src/data/meta.rs:326:        "TimsId",
./rustims/rustdf/src/data/meta.rs:329:        "NumScans",
./rustims/rustdf/src/data/meta.rs:330:        "NumPeaks",
./rustims/rustdf/src/data/handle.rs:587:        file_path.push("analysis.tdf_bin");
./rustims/rustdf/src/data/handle.rs:738:                panic!("TimsCompressionType is not 1 or 2.")
./rustims/rustdf/src/data/handle.rs:762:        file_path.push("analysis.tdf_bin");
./rustims/rustdf/src/data/handle.rs:812:                panic!("TimsCompressionType is not 1 or 2.")
./rustims/rustdf/src/data/handle.rs:1061:        file_path.push("analysis.tdf_bin");
./rustims/rustdf/src/data/handle.rs:1152:        file_path.push("analysis.tdf_bin");
./rustims/rustdf/src/data/handle.rs:1190:        file_path.push("analysis.tdf_bin");
./rustims/rustdf/src/data/handle.rs:1261:        file_path.push("analysis.tdf_bin");
./rustims/tims-viewer/src/data/meta.rs:4://! and a total-point estimate before touching `analysis.tdf_bin`, so the camera and
./rustims/rustdf/examples/compare_calibration.rs:61:            "SELECT T1, T2, NumScans, MzCalibration, TimsCalibration FROM Frames WHERE Id = ?1",
./rustims/rustdf/examples/compare_all_converters.rs:48:        "SELECT T1,T2,NumScans,MzCalibration,TimsCalibration FROM Frames WHERE Id=?1",
./rustims/packages/imspy-simulation/src/imspy_simulation/tdf.py:56:        self.binary_file = self.full_path / "analysis.tdf_bin"
./rustims/packages/imspy-simulation/src/imspy_simulation/tdf.py:62:        # When True, the per-frame tdf_bin realdata is produced by the Rust
./rustims/packages/imspy-simulation/src/imspy_simulation/tdf.py:107:            print(f"Setting up TDF file meta data, created: {self.full_path}/analysis.tdf and {self.full_path}/analysis.tdf_bin")
./rustims/packages/imspy-simulation/src/imspy_simulation/tdf.py:210:        r.TimsId = frame_start_pos
./rustims/packages/imspy-simulation/src/imspy_simulation/tdf.py:213:        r.NumScans = self.helper_handle.num_scans
./rustims/packages/imspy-simulation/src/imspy_simulation/tdf.py:214:        r.NumPeaks = len(intensity)

exec
/bin/bash -lc "sed -n '1,170p' rustims/rustdf/src/data/utility.rs && sed -n '180,245p' rustims/packages/imspy-simulation/src/imspy_simulation/tdf.py && sed -n '1040,1215p' rustims/rustdf/src/data/handle.rs" in /scratch/timsim-demo/timsim-necro
 succeeded in 0ms:
use byteorder::{ByteOrder, LittleEndian};
use mscore::timstof::frame::TimsFrame;
use rayon::iter::IntoParallelRefIterator;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::io;
use std::io::{Read, Write};

/// Decompresses a ZSTD compressed byte array
///
/// # Arguments
///
/// * `compressed_data` - A byte slice that holds the compressed data
///
/// # Returns
///
/// * `decompressed_data` - A vector of u8 that holds the decompressed data
///
pub fn zstd_decompress(compressed_data: &[u8]) -> io::Result<Vec<u8>> {
    let mut decoder = zstd::Decoder::new(compressed_data)?;
    let mut decompressed_data = Vec::new();
    decoder.read_to_end(&mut decompressed_data)?;
    Ok(decompressed_data)
}

/// Compresses a byte array using ZSTD
///
/// # Arguments
///
/// * `decompressed_data` - A byte slice that holds the decompressed data
///
/// # Returns
///
/// * `compressed_data` - A vector of u8 that holds the compressed data
///
pub fn zstd_compress(decompressed_data: &[u8], compression_level: i32) -> io::Result<Vec<u8>> {
    let mut encoder = zstd::Encoder::new(Vec::new(), compression_level)?;
    encoder.write_all(decompressed_data)?;
    let compressed_data = encoder.finish()?;
    Ok(compressed_data)
}

/// Deduplicate `(scan, tof)` pairs (summing their intensities) and return the
/// arrays sorted ascending by `(scan, tof)`.
///
/// The Bruker `tdf_bin` layout requires this ordering: [`modify_tofs`] delta-
/// encodes TOF within each scan (so TOFs must ascend within a scan) and
/// [`get_peak_cnts`] walks scans assuming they ascend. The Python writer
/// enforces the same invariant via an `np.unique` dedup + `np.lexsort((tof,
/// scan))` before encoding; the Rust write path previously fed raw, unsorted
/// frame data straight into the encoder, producing negative/garbage TOF deltas
/// that vendor readers (e.g. DiaNN) reject. Mirroring the Python preprocessing
/// here keeps the two writers byte-for-byte identical.
fn sort_dedup_scan_tof(
    scans: &[u32],
    tofs: &[u32],
    intensities: &[u32],
) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    use std::collections::HashMap;
    let mut acc: HashMap<(u32, u32), u64> = HashMap::with_capacity(scans.len());
    for i in 0..scans.len() {
        *acc.entry((scans[i], tofs[i])).or_insert(0) += intensities[i] as u64;
    }
    let mut pairs: Vec<((u32, u32), u64)> = acc.into_iter().collect();
    // Sort by scan, then tof — matches numpy's lexsort((tof, scan)).
    pairs.sort_unstable_by_key(|&((s, t), _)| (s, t));

    let n = pairs.len();
    let mut out_scan = Vec::with_capacity(n);
    let mut out_tof = Vec::with_capacity(n);
    let mut out_int = Vec::with_capacity(n);
    for ((s, t), inten) in pairs {
        out_scan.push(s);
        out_tof.push(t);
        out_int.push(inten.min(u32::MAX as u64) as u32);
    }
    (out_scan, out_tof, out_int)
}

pub fn reconstruct_compressed_data(
    scans: Vec<u32>,
    tofs: Vec<u32>,
    intensities: Vec<u32>,
    total_scans: u32,
    compression_level: i32,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    // Ensuring all vectors have the same length
    assert_eq!(scans.len(), tofs.len());
    assert_eq!(scans.len(), intensities.len());

    // Dedup + sort by (scan, tof) so TOF delta-encoding stays monotonic.
    let (scans, mut tofs, intensities) = sort_dedup_scan_tof(&scans, &tofs, &intensities);

    // Modify TOFs based on scans
    modify_tofs(&mut tofs, &scans);

    // Get peak counts from total scans and scans
    let peak_cnts = get_peak_cnts(total_scans, &scans);

    // Interleave TOFs and intensities
    let mut interleaved = Vec::new();
    for (&tof, &intensity) in tofs.iter().zip(intensities.iter()) {
        interleaved.push(tof);
        interleaved.push(intensity);
    }

    // Get real data using the custom loop logic
    let real_data = get_realdata(&peak_cnts, &interleaved);

    // Compress real_data using zstd_compress
    let compressed_data = zstd_compress(&real_data, compression_level)?;

    // Final data preparation with compressed data
    let mut final_data = Vec::new();

    // Include the length of the compressed data as a header (4 bytes)
    final_data.extend_from_slice(&(compressed_data.len() as u32 + 8).to_le_bytes());

    // Include total_scans as part of the header
    final_data.extend_from_slice(&total_scans.to_le_bytes());

    // Include the compressed data itself
    final_data.extend_from_slice(&compressed_data);

    Ok(final_data)
}

pub fn compress_collection(
    frames: Vec<TimsFrame>,
    max_scan_count: u32,
    compression_level: i32,
    num_threads: usize,
) -> Vec<Vec<u8>> {
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let result = pool.install(|| {
        frames
            .par_iter()
            .map(|frame| {
                let compressed_data = reconstruct_compressed_data(
                    frame.scan.iter().map(|&x| x as u32).collect(),
                    frame.tof.iter().map(|&x| x as u32).collect(),
                    frame
                        .ims_frame
                        .intensity
                        .iter()
                        .map(|&x| x as u32)
                        .collect(),
                    max_scan_count,
                    compression_level,
                )
                .unwrap();
                compressed_data
            })
            .collect()
    });
    result
}

/// Parses the decompressed bruker binary data
///
/// # Arguments
///
/// * `decompressed_bytes` - A byte slice that holds the decompressed data
///
/// # Returns
///
            frame: TimsFrame,
            scan_mode: int,
            frame_start_pos: int,
            only_frame_one: bool = False
    ):
        """Build a row for the frame meta data table from a TimsFrame object.
            Arguments:
                intensity: NDArray
                frame: TimsFrame object
                scan_mode: int
                frame_start_pos: int
                only_frame_one: bool
        """
        try:
            max_index = self.helper_handle.meta_data.Id.max()
        except AttributeError as e:
            max_index = self.helper_handle.meta_data.frame_id.max()

        r = self.helper_handle.meta_data.iloc[0, :].copy()
        if not only_frame_one:
            # check for index out of bounds since ref data handle might not hold same number of frames
            if frame.frame_id > max_index:
                r = self.helper_handle.meta_data.iloc[max_index - 1, :].copy()
            else:
                r = self.helper_handle.meta_data.iloc[frame.frame_id - 1, :].copy()

        r.Id = frame.frame_id
        r.Time = frame.retention_time
        r.ScanMode = scan_mode
        r.MsMsType = frame.ms_type
        r.TimsId = frame_start_pos
        r.MaxIntensity = int(np.max(intensity)) if len(intensity) > 0 else 0
        r.SummedIntensities = int(np.sum(intensity)) if len(intensity) > 0 else 0
        r.NumScans = self.helper_handle.num_scans
        r.NumPeaks = len(intensity)

        return r

    def compress_frame(self, frame: TimsFrame, only_frame_one: bool = False) -> (NDArray, bytes):
        """Compress a single frame using zstd.
            Arguments:
                frame: TimsFrame object
                only_frame_one: bool

            Returns:
                bytes: intensities, compressed data
        """
        # either use frame 1 or the ref handle frame for writing of calibration data and call to conversion function
        i = 1 if only_frame_one else frame.frame_id

        try:
            max_index = self.helper_handle.meta_data.Id.max()
        except AttributeError as e:
            max_index = self.helper_handle.meta_data.frame_id.max()

        if frame.frame_id > max_index and not only_frame_one:
            i = max_index

        # transform mz and mobility to tof and scan
        tof = self.mz_to_tof(i, frame.mz).astype(np.uint32)
        scan = self.inv_mobility_to_scan(i, frame.mobility).astype(np.uint32)
        intensity = frame.intensity.astype(np.uint32)

        # Since, mz -> tof is not bijective, we need to check for duplicates
        # stack scan and tof to form a 2D array for unique grouping
        scan_tof = np.stack((scan, tof), axis=1)
        scan_max_index: u32,
        im_lower: f64,
        im_upper: f64,
        tof_max_index: u32,
        mz_lower: f64,
        mz_upper: f64,
    ) -> Self {
        let raw_data_layout = TimsRawDataLayout::new(data_path);
        let index_converter = build_index_converter(
            bruker_lib_path,
            data_path,
            use_bruker_sdk,
            scan_max_index,
            im_lower,
            im_upper,
            tof_max_index,
            mz_lower,
            mz_upper,
        );

        let mut file_path = PathBuf::from(data_path);
        file_path.push("analysis.tdf_bin");
        let mut infile = File::open(file_path).unwrap();
        let mut data = Vec::new();
        infile.read_to_end(&mut data).unwrap();

        TimsDataLoader::InMemory(TimsInMemoryLoader {
            raw_data_layout,
            index_converter,
            compressed_data: data,
        })
    }

    /// Create a lazy loader with pre-computed ion mobility calibration lookup table.
    ///
    /// This method enables accurate ion mobility calibration with fast parallel extraction.
    /// The im_lookup table should be pre-computed using the Bruker SDK.
    ///
    /// # Arguments
    /// * `data_path` - Path to the .d folder
    /// * `bruker_lib_path` - Path to the Bruker SDK shared library; used to
    ///   derive an accurate m/z calibration. Pass "NO_SDK" (or an empty
    ///   string) to skip and use the 2-point boundary m/z model.
    /// * `tof_max_index` - Maximum TOF index (from GlobalMetaData)
    /// * `mz_lower` - Minimum m/z value (from GlobalMetaData)
    /// * `mz_upper` - Maximum m/z value (from GlobalMetaData)
    /// * `im_lookup` - Pre-computed scan→1/K0 lookup table
    ///
    /// # Returns
    /// A new TimsDataLoader with LookupIndexConverter
    pub fn new_lazy_with_calibration(
        data_path: &str,
        bruker_lib_path: &str,
        tof_max_index: u32,
        mz_lower: f64,
        mz_upper: f64,
        im_lookup: Vec<f64>,
    ) -> Self {
        let raw_data_layout = TimsRawDataLayout::new(data_path);

        let index_converter = TimsIndexConverter::Lookup(build_lookup_converter(
            bruker_lib_path,
            data_path,
            tof_max_index,
            mz_lower,
            mz_upper,
            im_lookup,
        ));

        TimsDataLoader::Lazy(TimsLazyLoder {
            raw_data_layout,
            index_converter,
        })
    }

    /// Create an in-memory loader with pre-computed ion mobility calibration lookup table.
    ///
    /// This method enables accurate ion mobility calibration with fast parallel extraction.
    /// The im_lookup table should be pre-computed using the Bruker SDK.
    ///
    /// # Arguments
    /// * `data_path` - Path to the .d folder
    /// * `bruker_lib_path` - Path to the Bruker SDK shared library; used to
    ///   derive an accurate m/z calibration. Pass "NO_SDK" (or an empty
    ///   string) to skip and use the 2-point boundary m/z model.
    /// * `tof_max_index` - Maximum TOF index (from GlobalMetaData)
    /// * `mz_lower` - Minimum m/z value (from GlobalMetaData)
    /// * `mz_upper` - Maximum m/z value (from GlobalMetaData)
    /// * `im_lookup` - Pre-computed scan→1/K0 lookup table
    ///
    /// # Returns
    /// A new TimsDataLoader with LookupIndexConverter
    pub fn new_in_memory_with_calibration(
        data_path: &str,
        bruker_lib_path: &str,
        tof_max_index: u32,
        mz_lower: f64,
        mz_upper: f64,
        im_lookup: Vec<f64>,
    ) -> Self {
        let raw_data_layout = TimsRawDataLayout::new(data_path);

        let index_converter = TimsIndexConverter::Lookup(build_lookup_converter(
            bruker_lib_path,
            data_path,
            tof_max_index,
            mz_lower,
            mz_upper,
            im_lookup,
        ));

        let mut file_path = PathBuf::from(data_path);
        file_path.push("analysis.tdf_bin");
        let mut infile = File::open(file_path).unwrap();
        let mut data = Vec::new();
        infile.read_to_end(&mut data).unwrap();

        TimsDataLoader::InMemory(TimsInMemoryLoader {
            raw_data_layout,
            index_converter,
            compressed_data: data,
        })
    }

    /// Create a lazy loader using the exact SDK-free Bruker calibration formulas.
    ///
    /// Builds a [`BrukerFormulaConverter`] from the `MzCalibration` /
    /// `TimsCalibration` tables (frame `calibration_frame_id`, default 1 — the
    /// coefficients are near-constant per run). Needs no Bruker SDK at build or
    /// runtime; 1/K0 is machine-exact and m/z is bit-exact for MzCalibration
    /// ModelType 1 (few ppm for ModelType 2).
    pub fn new_lazy_with_bruker_formula(data_path: &str, calibration_frame_id: u32) -> Self {
        let raw_data_layout = TimsRawDataLayout::new(data_path);
        let index_converter = TimsIndexConverter::BrukerFormula(
            BrukerFormulaConverter::from_d_folder(data_path, calibration_frame_id).unwrap(),
        );
        TimsDataLoader::Lazy(TimsLazyLoder {
            raw_data_layout,
            index_converter,
        })
    }

    /// In-memory counterpart of [`Self::new_lazy_with_bruker_formula`].
    pub fn new_in_memory_with_bruker_formula(data_path: &str, calibration_frame_id: u32) -> Self {
        let raw_data_layout = TimsRawDataLayout::new(data_path);
        let index_converter = TimsIndexConverter::BrukerFormula(
            BrukerFormulaConverter::from_d_folder(data_path, calibration_frame_id).unwrap(),
        );

        let mut file_path = PathBuf::from(data_path);
        file_path.push("analysis.tdf_bin");
        let mut infile = File::open(file_path).unwrap();
        let mut data = Vec::new();
        infile.read_to_end(&mut data).unwrap();

        TimsDataLoader::InMemory(TimsInMemoryLoader {
            raw_data_layout,
            index_converter,
            compressed_data: data,
        })
    }

    /// Create a lazy loader with full calibration (both m/z and IM).
    ///
    /// This method uses regression-derived m/z calibration coefficients instead of
    /// the simple boundary model, providing more accurate m/z conversion.
    ///
    /// # Arguments
    /// * `data_path` - Path to the .d folder
    /// * `tof_intercept` - Intercept for sqrt(mz) = intercept + slope * tof
    /// * `tof_slope` - Slope for sqrt(mz) = intercept + slope * tof
    /// * `im_min` - Minimum 1/K0 value
    /// * `im_max` - Maximum 1/K0 value
    /// * `scan_max_index` - Maximum scan index
    pub fn new_lazy_with_mz_calibration(
        data_path: &str,

The sweep-line is the right core, but the document currently understates two important dimensions: input ordering and diaPASEF’s *scan-dependent* isolation.

1. **Streaming/memory.** An active set keyed by `frame_end` is the correct decomposition. Its bound is genuinely proportional to local chromatographic density × elution support, not run duration—provided precursors arrive in nondecreasing `frame_start`. That proviso is material: ordinary Parquet scans are not automatically RT-sorted, so either write a start-sorted artifact or perform an external/partitioned sort; do not accidentally materialize all precursors to sort them.

   The active heap is not the only memory pressure. A dense DIA frame can contain many co-isolated precursors × fragments × mobility bins; “a few thousand” sparse bins is optimistic. Hash-map overhead, duplicate aggregation, reference-frame data, and parallel compression queues can dominate. Bound all of those: bounded producer/compressor channels, per-frame peak cap/overflow policy, peak-RSS telemetry, and a stress fixture with concentrated RT and broad elution tails. “Cap/spill” is acceptable only if the spill preserves deterministic reduction order and does not silently alter intensities.

2. **TDF writer: append-only is plausible and likely the lowest-risk result, not an assumption.** `analysis.tdf_bin` is a concatenation of self-delimiting per-frame blocks; `Frames.TimsId` is the byte offset to each frame. The repository’s current encoder already emits a 4-byte block length followed by 4-byte scan count and compressed data, and writes `TimsId` as the current byte position. That supports append/write-frame → insert metadata row in one transaction, with no binary global index or final offset pass. AlphaTims independently reads a frame by seeking to `Frames.TimsId` then reading its length header. [AlphaTims binary reader](https://alphatims.readthedocs.io/en/latest/_modules/alphatims/bruker.html)

   Still, prove this against the target Bruker reader: compression type must agree with `GlobalMetadata.TimsCompressionType`; scan count, `NumPeaks`, sorted/deduplicated `(scan,tof)`, and calibration/metadata tables must be valid. The local writer notes that unsorted TOFs yield invalid delta coding. SQLite aggregates such as run/frame count can be finalized with end-of-run SQL updates; they do not require retaining frame payloads. Add a crash/partial-output policy—an interrupted `.d` must be explicitly invalid or recoverable, never superficially complete.

3. **§3.4 is too idealized for benchmark-realistic DIA.** Factorization is a good base model and conservation oracle, but `MS2 if m/z ∈ W` is wrong for diaPASEF. The quadrupole position can change during the TIMS ramp, so transmission is a function of scan/mobility and m/z—effectively a 2-D, often diagonal acquisition mask, not one fixed frame window. Standard diaPASEF deliberately positions isolation as a function of TIMS elution time. [diaPASEF description](https://pmc.ncbi.nlm.nih.gov/articles/PMC9465115/)

   Replace the binary predicate with an instrument/method-specific transmission term:
   `I = A × E(frame) × M(scan) × Q(mz, scan, frame) × CE(mz,charge,scan) × fragment_payload`.
   `Q` should have measured/parameterized edge roll-off and the actual stepping/diagonal window geometry; `CE` should control charge- and m/z-dependent fragmentation yield/pattern. Otherwise edge precursors are unrealistically fully transmitted, chimericity is wrong, and the benchmark can overstate search performance. Also model detector saturation/dynamic range and decide explicitly whether peptide abundance is apportioned once across charge states (`charge_fraction`) before isotope emission—otherwise conservation can double-count.

4. **DIA-first is right; DDA need not imply a rewrite if one seam is designed now.** Make the sweep consume an `AcquisitionEvent` stream rather than a precomputed frame role/window. DIA supplies it statically. DDA’s online controller consumes rendered/centroided MS1 evidence, exclusion state, and timing, then emits PASEF MS2 events. The active precursor store, RT sweep, and frame writer remain valid; what changes is scheduling and selection, not feature streaming. A two-pass DDA implementation would be less faithful because selection depends on the realized noisy MS1 signal.

5. **Validation is too weak.** TIC correlation and search recovery are aggregate, forgiving tests; a search engine can identify the same peptides while isolation leakage, mobility placement, cofragmentation, intensity rank, and calibration are wrong. Keep them, but add:

   - Golden, noise-free microfixtures with exact per-bin `(frame, scan, tof, intensity)` comparisons after quantization.
   - Per-precursor 2-D RT×mobility moments/widths, isotope ratios, fragment rank/intensity distributions, and isolation-edge transmission curves.
   - Per-frame peak count, base-peak, saturation/clipping, m/z/mobility marginal distributions, and MS2 interference/chimericity metrics.
   - Reader interoperability: vendor SDK/AlphaTims/DIA-NN open and round-trip the emitted TDF, plus metadata/offset consistency checks.
   - Search accuracy stratified by abundance, charge, mobility, window edge, and coelution—not just total IDs.

Finally, clarify whether “frame” means one TIMS ramp or a full MS1+MS2 duty cycle; the 36k-at-100-ms arithmetic and schedule wording currently blur those. That affects RT integration, event timing, and DIA window assignment.
tokens used
21,538
The sweep-line is the right core, but the document currently understates two important dimensions: input ordering and diaPASEF’s *scan-dependent* isolation.

1. **Streaming/memory.** An active set keyed by `frame_end` is the correct decomposition. Its bound is genuinely proportional to local chromatographic density × elution support, not run duration—provided precursors arrive in nondecreasing `frame_start`. That proviso is material: ordinary Parquet scans are not automatically RT-sorted, so either write a start-sorted artifact or perform an external/partitioned sort; do not accidentally materialize all precursors to sort them.

   The active heap is not the only memory pressure. A dense DIA frame can contain many co-isolated precursors × fragments × mobility bins; “a few thousand” sparse bins is optimistic. Hash-map overhead, duplicate aggregation, reference-frame data, and parallel compression queues can dominate. Bound all of those: bounded producer/compressor channels, per-frame peak cap/overflow policy, peak-RSS telemetry, and a stress fixture with concentrated RT and broad elution tails. “Cap/spill” is acceptable only if the spill preserves deterministic reduction order and does not silently alter intensities.

2. **TDF writer: append-only is plausible and likely the lowest-risk result, not an assumption.** `analysis.tdf_bin` is a concatenation of self-delimiting per-frame blocks; `Frames.TimsId` is the byte offset to each frame. The repository’s current encoder already emits a 4-byte block length followed by 4-byte scan count and compressed data, and writes `TimsId` as the current byte position. That supports append/write-frame → insert metadata row in one transaction, with no binary global index or final offset pass. AlphaTims independently reads a frame by seeking to `Frames.TimsId` then reading its length header. [AlphaTims binary reader](https://alphatims.readthedocs.io/en/latest/_modules/alphatims/bruker.html)

   Still, prove this against the target Bruker reader: compression type must agree with `GlobalMetadata.TimsCompressionType`; scan count, `NumPeaks`, sorted/deduplicated `(scan,tof)`, and calibration/metadata tables must be valid. The local writer notes that unsorted TOFs yield invalid delta coding. SQLite aggregates such as run/frame count can be finalized with end-of-run SQL updates; they do not require retaining frame payloads. Add a crash/partial-output policy—an interrupted `.d` must be explicitly invalid or recoverable, never superficially complete.

3. **§3.4 is too idealized for benchmark-realistic DIA.** Factorization is a good base model and conservation oracle, but `MS2 if m/z ∈ W` is wrong for diaPASEF. The quadrupole position can change during the TIMS ramp, so transmission is a function of scan/mobility and m/z—effectively a 2-D, often diagonal acquisition mask, not one fixed frame window. Standard diaPASEF deliberately positions isolation as a function of TIMS elution time. [diaPASEF description](https://pmc.ncbi.nlm.nih.gov/articles/PMC9465115/)

   Replace the binary predicate with an instrument/method-specific transmission term:
   `I = A × E(frame) × M(scan) × Q(mz, scan, frame) × CE(mz,charge,scan) × fragment_payload`.
   `Q` should have measured/parameterized edge roll-off and the actual stepping/diagonal window geometry; `CE` should control charge- and m/z-dependent fragmentation yield/pattern. Otherwise edge precursors are unrealistically fully transmitted, chimericity is wrong, and the benchmark can overstate search performance. Also model detector saturation/dynamic range and decide explicitly whether peptide abundance is apportioned once across charge states (`charge_fraction`) before isotope emission—otherwise conservation can double-count.

4. **DIA-first is right; DDA need not imply a rewrite if one seam is designed now.** Make the sweep consume an `AcquisitionEvent` stream rather than a precomputed frame role/window. DIA supplies it statically. DDA’s online controller consumes rendered/centroided MS1 evidence, exclusion state, and timing, then emits PASEF MS2 events. The active precursor store, RT sweep, and frame writer remain valid; what changes is scheduling and selection, not feature streaming. A two-pass DDA implementation would be less faithful because selection depends on the realized noisy MS1 signal.

5. **Validation is too weak.** TIC correlation and search recovery are aggregate, forgiving tests; a search engine can identify the same peptides while isolation leakage, mobility placement, cofragmentation, intensity rank, and calibration are wrong. Keep them, but add:

   - Golden, noise-free microfixtures with exact per-bin `(frame, scan, tof, intensity)` comparisons after quantization.
   - Per-precursor 2-D RT×mobility moments/widths, isotope ratios, fragment rank/intensity distributions, and isolation-edge transmission curves.
   - Per-frame peak count, base-peak, saturation/clipping, m/z/mobility marginal distributions, and MS2 interference/chimericity metrics.
   - Reader interoperability: vendor SDK/AlphaTims/DIA-NN open and round-trip the emitted TDF, plus metadata/offset consistency checks.
   - Search accuracy stratified by abundance, charge, mobility, window edge, and coelution—not just total IDs.

Finally, clarify whether “frame” means one TIMS ramp or a full MS1+MS2 duty cycle; the 36k-at-100-ms arithmetic and schedule wording currently blur those. That affects RT integration, event timing, and DIA window assignment.
