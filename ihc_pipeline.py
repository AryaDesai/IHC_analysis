"""IHC gradient profiles along the PSM: masks, band measurements and phase plots.

The notebook IHC_analysis.ipynb runs the analysis through this module. The
image files are read without being changed, and the masks, ROIs and results are
written to Config.output_root. Each setting is defined in Config together with
its description, which parameter_table() shows in the notebook.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import asdict, dataclass, field, fields, replace
from functools import lru_cache
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import warnings

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.signal import savgol_filter
from skimage.filters import threshold_isodata

# Timestamped progress messages and progress bars. Setting ihc.PROGRESS = False
# in the notebook turns them off.
PROGRESS = True


def _log(message):
    """Print a timestamped progress message."""
    if PROGRESS:
        from tqdm.auto import tqdm

        # tqdm.write clears any progress bar, prints the line and redraws the bar,
        # so the message does not land in the middle of a half-drawn bar. The bars
        # draw on stderr, so the message must go there too.
        tqdm.write(f"[{time.strftime('%H:%M:%S')}] {message}", file=sys.stderr)


def _bar(iterable=None, desc="", total=None, disable=False, **kwargs):
    """Progress bar around an iterable, hidden when PROGRESS is False."""
    from tqdm.auto import tqdm

    return tqdm(
        iterable,
        desc=desc,
        total=total,
        disable=disable or not PROGRESS,
        leave=False,
        dynamic_ncols=True,
        **kwargs,
    )


def _since(start):
    """Seconds since start, as text."""
    return f"{time.perf_counter() - start:.1f} s"


def setting(default, group, description, choices=None):
    """A Config field whose default, group and description feed parameter_table()."""
    metadata = dict(group=group, description=description, choices=choices)
    if isinstance(default, dict):  # each Config gets its own copy of a dict default
        return field(default_factory=lambda: dict(default), metadata=metadata)
    return field(default=default, metadata=metadata)


@dataclass
class Config:
    experimental_data: str | Path = setting(
        "", "Files",
        "Folder that contains the replicate folders (N_1, N_2, ...). The folder names can "
        "have more text after the number, as in N_1_oct_6.",
    )
    output_root: str | Path = setting(
        "", "Files",
        "Folder where masks, ROIs, measurements and figures are saved. It has to be "
        "separate from the image folders, and a run stops if it lies inside one.",
    )
    replicates: tuple = setting(
        (), "Files",
        "Replicates to analyse, for example ('N_1', 'N_2'). Leave empty to use every replicate "
        "found.",
    )
    image_glob: str = setting(
        "*.tiff", "Files",
        "Filename pattern of the exported multichannel TIFFs. TIFFs with fewer channels "
        "than the settings need are skipped.",
    )
    include_nd2: bool = setting(
        True, "Files",
        "Also look for ND2 files. When a position has also been exported as a TIFF, the "
        "TIFF is used for that position.",
    )
    denoised_tokens: tuple = setting(
        (" - Denoised", "_denoised"), "Files",
        "Text that marks a filename as an NIS-denoised export. The original image has the "
        "same name with this text removed.",
    )
    nis_preference: str = setting(
        "nd2", "Files",
        "The denoised export that supplies the DAPI when a sample has both an ND2 and a "
        "TIFF version. If the preferred export is missing, the other one is used.",
        ("nd2", "tiff"),
    )
    phase_files: dict = setting(
        {}, "Phases",
        "Phase table for each replicate, for example {'N_1': '/path/phases.csv'}. A replicate "
        "not listed here uses phases.csv in its folder, if there is one.",
    )
    phase_overrides: dict = setting(
        {}, "Phases",
        "Phases set by hand, for example {'N_2/triple_ihc_014_P1': '2'}, which take the "
        "place of the entries in the phase table.",
    )
    dapi_channel: int = setting(
        1, "Channels",
        "DAPI channel in the original image, counting from 1. An NIS export that contains "
        "only DAPI is read from its single channel.",
    )
    protein_channels: dict = setting(
        {"Her7-Venus": 2, "pMEK": 3, "ppERK": 4}, "Channels",
        "Protein names and their channels in the original image, counting from 1. The "
        "protein signal is read from the original image, and NIS exports are used only for "
        "the DAPI.",
    )
    pixel_size_um: float | None = setting(
        None, "Calibration",
        "Pixel size in µm. None reads it from the image file.",
    )
    z_step_um: float | None = setting(
        None, "Calibration",
        "Distance between Z-planes in µm. None reads it from the image file.",
    )
    use_ai_denoise: bool = setting(
        True, "DAPI",
        "If True, the DAPI is taken from the sample's NIS-denoised export, so every sample "
        "needs one. If False, the raw DAPI is median filtered, background subtracted and "
        "contrast enhanced with CLAHE before segmentation.",
    )
    use_cellpose: bool = setting(
        True, "Mask",
        "If True, find nuclei with Cellpose. If False, threshold each Z-plane with the IsoData "
        "method.",
    )
    mask_source: str = setting(
        "compute", "Mask",
        "Where the nuclear mask comes from. 'compute' makes it here, 'existing' reads the "
        "Mask.tif in the image's own folder, and 'edited' reads a mask corrected in Fiji and "
        "saved as edited_masks/<replicate>/<sample>/Mask.tif in the output folder.",
        ("compute", "existing", "edited"),
    )
    dapi_median_size: int = setting(
        3, "DAPI filtering",
        "Width in pixels of the median filter applied to the raw DAPI when use_ai_denoise is "
        "False. Use an odd number, or 0 to skip the filter.",
    )
    background_size_px: int = setting(
        101, "DAPI filtering",
        "Width in pixels of the square used to estimate the smooth background under the "
        "nuclei, which is then subtracted. The square should be wider than a nucleus, and "
        "0 skips this step.",
    )
    clahe_window_px: int = setting(
        31, "Contrast",
        "Approximate width in pixels of the tiles that CLAHE equalizes separately.",
    )
    clahe_clip: float = setting(
        3.0, "Contrast",
        "CLAHE clip limit, as a multiple of the average histogram bin (the 'slope' in Fiji). "
        "Higher values increase local contrast more.",
    )
    cellpose_model: str = setting(
        "cpsam_v2", "Cellpose",
        "Cellpose 4 model name, or the path to a weights file.",
    )
    cellpose_diameter_px: float | None = setting(
        None, "Cellpose",
        "Typical nuclear diameter in pixels. When it is set, Cellpose rescales the image "
        "so the nuclei match the size it was trained on. None leaves the image at its own "
        "scale.",
    )
    cellpose_percentiles: tuple = setting(
        (1.0, 99.0), "Cellpose",
        "DAPI percentiles that Cellpose maps to 0 and 1 before segmenting.",
    )
    cellprob_threshold: float = setting(
        0.0, "Cellpose",
        "Cellpose's threshold on its nucleus probability score. Lower values make masks larger.",
    )
    flow_threshold: float = setting(
        0.4, "Cellpose",
        "How closely a nucleus's shape must agree with Cellpose's predicted flows to be kept. "
        "Lower values are stricter.",
    )
    min_size_px: int = setting(
        15, "Cellpose",
        "Smallest nucleus kept, as an area in pixels within one Z-plane.",
    )
    max_size_fraction: float = setting(
        0.4, "Cellpose",
        "Largest nucleus kept, as a fraction of the image area.",
    )
    flow_iterations: int | None = setting(
        None, "Cellpose",
        "Number of steps Cellpose takes when following its flows to group pixels into nuclei. "
        "None uses Cellpose's default.",
    )
    tile_overlap: float = setting(
        0.1, "Cellpose",
        "Fraction by which Cellpose's 256-pixel tiles overlap.",
    )
    augment: bool = setting(
        False, "Cellpose",
        "Also segment flipped copies of each tile and average the results, which takes "
        "longer.",
    )
    batch_size: int = setting(
        32, "Compute",
        "Number of image tiles Cellpose sends to the GPU together. Larger batches are faster "
        "until GPU memory runs out.",
    )
    device: str = setting(
        "auto", "Compute",
        "Where Torch runs. 'auto' uses an NVIDIA GPU (CUDA) if there is one, then an Apple GPU "
        "(MPS), then the CPU.",
        ("auto", "cpu", "mps", "cuda"),
    )
    image_workers: int = setting(
        1, "Compute",
        "Number of samples measured in parallel, each in its own process. The masks are "
        "made beforehand, one sample at a time, so that only one copy of the Cellpose "
        "model is on the GPU.",
    )
    threads_per_worker: int = setting(
        8, "Compute",
        "Threads each process uses for filtering Z-planes and for Torch.",
    )
    reuse_cache: bool = setting(
        True, "Compute",
        "Load a sample's saved mask when one exists, instead of making it again. The "
        "settings and code that made the saved mask are not checked, so set this to False "
        "after changing either. Measurements are made again on every run, and RELOAD_SAVED "
        "in the notebook loads the saved ones instead.",
    )
    export_csv: bool = setting(
        False, "Files",
        "Also save each sample's measurements as CSV, next to the Parquet file.",
    )
    width_override_px: int | None = setting(
        None, "ROI",
        "Band width in pixels to use for every side. None uses the width saved with each ROI.",
    )
    interpolation: str = setting(
        "linear", "ROI",
        "How the band reads the image and the mask between pixel centres.",
        ("linear", "nearest"),
    )
    boundary_mode: str = setting(
        "zero", "ROI",
        "What to do where a band runs off the image. 'zero' counts those points as 0, 'exclude'"
        " leaves them out of the average, and 'error' stops the run.",
        ("zero", "exclude", "error"),
    )
    methods: tuple = setting(
        ("block_mip", "per_slice", "pooled", "full_mip"), "Measurement",
        "Ways of handling Z to measure. 'per_slice' measures each plane, 'pooled' adds the "
        "planes together, 'block_mip' takes maximum projections of blocks of planes, and "
        "'full_mip' projects the whole stack.",
    )
    z_start: int = setting(1, "Measurement", "First Z-plane to use, counting from 1.")
    z_stop: int | None = setting(
        None, "Measurement",
        "Last Z-plane to use, counting from 1. None uses the last plane.",
    )
    mip_block_size: int = setting(
        3, "Measurement",
        "Number of Z-planes in each block maximum projection.",
    )
    mip_start: int = setting(
        2, "Measurement",
        "First Z-plane of the first block, counting from 1. With 31 planes, starting at 2 with "
        "blocks of 3 gives 2–4, 5–7, ..., 29–31.",
    )
    incomplete_block: str = setting(
        "drop", "Measurement",
        "What to do with a last block that has too few planes: drop it, include it, or stop "
        "with an error.",
        ("drop", "include", "error"),
    )
    min_mask_fraction: float = setting(
        0.05, "Measurement",
        "Smallest fraction of the band a compartment must cover for its mean to be "
        "reported. Where the coverage is lower, the mean is left empty (NaN).",
    )
    protein_median_size: int = setting(
        3, "Protein filtering",
        "Width in pixels of the median filter applied to the protein channels. Use an odd "
        "number, or 0 to skip the filter.",
    )
    protein_outliers: bool = setting(
        True, "Protein filtering",
        "Replace pixels that are much brighter than their surroundings with the local median.",
    )
    outlier_radius_px: float = setting(
        2.0, "Protein filtering",
        "Radius in pixels of the neighbourhood used to find bright outliers.",
    )
    outlier_threshold: float = setting(
        50.0, "Protein filtering",
        "How far above the local median, in intensity units, a pixel must be to be replaced.",
    )
    mask_cap: float | None = setting(
        4095.0, "Protein filtering",
        "Protein intensities above this value are clipped before the nuclear and "
        "cytoplasmic measurements. 4095 is the maximum of the 12-bit camera, and None "
        "turns clipping off. The whole-band measurement uses the unclipped values.",
    )
    exclusions_csv: str = setting(
        "", "Exclusions",
        "Optional table of distance ranges to leave out of the analysis, with columns "
        "replicate, sample_id, side, start_um, stop_um, channels and reason.",
    )
    analysis_method: str = setting(
        "block_mip", "Analysis",
        "Z method used for the profiles and phase plots, which has to be one of those "
        "listed in methods.",
    )
    z_average: str = setting(
        "pool", "Analysis",
        "How the Z-planes or blocks of a side are combined. 'pool' weights each by its mask "
        "coverage, 'mean' weights them equally, and 'median' takes the middle value.",
        ("pool", "mean", "median"),
    )
    combine_sides: bool = setting(
        False, "Analysis",
        "If True, combine an embryo's sides into one profile. If False, each side is its own "
        "profile.",
    )
    side_average: str = setting(
        "mean", "Analysis",
        "How sides are combined when combine_sides is True. 'mean' weights sides equally and "
        "'pool' weights them by mask coverage.",
        ("mean", "pool"),
    )
    align_profiles: int = setting(
        0, "Analysis",
        "End of the band that profiles are aligned to and distances are measured from. 0 "
        "aligns the posterior ends (distance runs posterior to anterior) and 1 aligns the "
        "anterior ends, the last point of each band (distance runs anterior to posterior). "
        "Exclusions and day_reference_um are always measured from the posterior end.",
        (0, 1),
    )
    background_subtract: bool = setting(
        True, "Normalization",
        "Subtract a baseline from each embryo's profiles. The baseline is a low percentile "
        "of the profiles themselves rather than a separate measurement of the optical "
        "background.",
    )
    background_percentile: float = setting(
        5.0, "Normalization",
        "Percentile of the smoothed side profiles that sets the embryo's baseline.",
    )
    normalization_window_um: float = setting(
        40.0, "Normalization",
        "Width in µm of the rolling median that smooths the profiles before the baseline "
        "percentile is taken. The baseline is subtracted from the unsmoothed profile.",
    )
    day_normalize: bool = setting(
        False, "Normalization",
        "After background subtraction, divide each imaging day's profiles by that day's median "
        "value, separately for each channel and compartment. The imaging day is the ND2 "
        "acquisition date, or the replicate when there is no date.",
    )
    day_reference_um: tuple = setting(
        (0.0, 400.0), "Normalization",
        "Distance range in µm, measured from the posterior end, whose values set each day's "
        "median.",
    )
    average: str = setting(
        "mean", "Statistics",
        "Average used for the summary across all profiles.",
        ("mean", "median"),
    )
    error: str = setting(
        "sd", "Statistics",
        "Spread shown around the overall average: standard deviation, standard error, or "
        "interquartile range.",
        ("sd", "sem", "iqr"),
    )
    minimum_lines: int = setting(
        5, "Statistics",
        "Fewest profiles needed at a distance for the overall average to be drawn there.",
    )
    phase_average: str = setting(
        "mean", "Phases",
        "Average used for each phase.",
        ("mean", "median"),
    )
    phase_error: str = setting(
        "sem", "Phases",
        "Spread shown around each phase average: standard deviation, standard error, or "
        "interquartile range.",
        ("sd", "sem", "iqr"),
    )
    phase_minimum_lines: int = setting(
        2, "Phases",
        "Fewest profiles needed at a distance for a phase average to be drawn there.",
    )
    phase_labels: tuple = setting(
        ("1", "2", "3"), "Phases",
        "Phase labels to summarise and plot, in order. Samples with other labels or no "
        "label are shown with the individual profiles and left out of the phase averages.",
    )
    phase_compartments: dict = setting(
        {"Her7-Venus": "nuclear", "ppERK": "nuclear", "pMEK": "cytoplasmic"}, "Phases",
        "Compartment shown for each channel in the phase plots. Every compartment is "
        "measured, so this setting only changes what is plotted.",
    )
    smooth_average: bool = setting(
        True, "Plot",
        "Smooth the plotted average with a Savitzky–Golay filter. The saved statistics and "
        "the error bands are calculated from the unsmoothed profiles.",
    )
    savgol_window_um: float = setting(60.0, "Plot", "Width in µm of the smoothing window.")
    savgol_order: int = setting(1, "Plot", "Degree of the smoothing polynomial.")
    plot_channels: tuple = setting(
        ("Her7-Venus", "ppERK", "pMEK"), "Plot",
        "Channels to plot. Each must be a name in protein_channels.",
    )
    plot_compartments: tuple = setting(
        ("nuclear", "cytoplasmic", "all"), "Plot",
        "Compartments to plot.",
    )
    phase_colormap: str = setting(
        "viridis", "Plot",
        "Matplotlib colormap for the phases, which are given evenly spaced colours from "
        "it.",
    )
    phase_colormap_range: tuple = setting(
        (0.0, 0.85), "Plot",
        "Part of the colormap to use, between 0 and 1. The default skips the pale yellow end of"
        " viridis, which is hard to see on a light background.",
    )
    show_individuals: bool = setting(
        True, "Plot",
        "Draw each profile behind the phase average in the per-phase panels.",
    )
    dpi: int = setting(150, "Plot", "Resolution of saved figures, in dots per inch.")
    qc_percentiles: tuple = setting(
        (1.0, 99.7), "Plot",
        "Intensity percentiles that set the display contrast of the QC figure.",
    )

    def validate(self):
        """Check every setting before any image is read, and return self."""
        import matplotlib

        def fail_unless(ok, message):
            if not ok:
                raise ValueError(message)

        fail_unless(
            str(self.experimental_data) and Path(self.experimental_data).is_dir(),
            "experimental_data must be an existing folder",
        )
        fail_unless(str(self.output_root), "Set output_root to a folder outside the image folders")
        for f in fields(self):
            choices = f.metadata.get("choices")
            fail_unless(
                not choices or getattr(self, f.name) in choices, f"{f.name}: choose one of {choices}"
            )
        for name in (
            "batch_size", "image_workers", "threads_per_worker", "mip_block_size", "z_start",
            "mip_start", "minimum_lines", "phase_minimum_lines", "dpi", "dapi_channel",
            "clahe_window_px",
        ):  # fmt: skip
            x = getattr(self, name)
            fail_unless(
                isinstance(x, int) and not isinstance(x, bool) and x >= 1,
                f"{name} must be a positive integer",
            )
        for name in ("protein_median_size", "dapi_median_size"):
            x = getattr(self, name)
            fail_unless(x == 0 or (x > 0 and x % 2 == 1), f"{name} must be 0 or an odd number")
        for name in ("pixel_size_um", "z_step_um", "cellpose_diameter_px", "width_override_px", "mask_cap"):
            x = getattr(self, name)
            fail_unless(x is None or (np.isfinite(x) and x > 0), f"{name} must be positive or None")
        for name in (
            "outlier_radius_px", "normalization_window_um",
            "savgol_window_um", "clahe_clip",
        ):  # fmt: skip
            x = getattr(self, name)
            fail_unless(np.isfinite(x) and x > 0, f"{name} must be positive")
        for name in (
            "min_size_px", "flow_threshold", "outlier_threshold", "savgol_order", "background_size_px",
        ):  # fmt: skip
            fail_unless(getattr(self, name) >= 0, f"{name} cannot be negative")
        fail_unless(
            self.width_override_px is None or int(self.width_override_px) == self.width_override_px,
            "width_override_px must be a whole number",
        )
        fail_unless(0 < self.min_mask_fraction <= 1, "min_mask_fraction must be in (0, 1]")
        fail_unless(0 <= self.tile_overlap < 1, "tile_overlap must be in [0, 1)")
        fail_unless(0 < self.max_size_fraction <= 1, "max_size_fraction must be in (0, 1]")
        fail_unless(0 <= self.background_percentile <= 100, "background_percentile must be in [0, 100]")
        fail_unless(
            self.flow_iterations is None or self.flow_iterations >= 1,
            "flow_iterations must be positive or None",
        )
        for name in ("cellpose_percentiles", "qc_percentiles"):
            x = getattr(self, name)
            fail_unless(len(x) == 2 and 0 <= x[0] < x[1] <= 100, f"{name} must be (low, high) within [0, 100]")
        fail_unless(
            0 <= self.day_reference_um[0] < self.day_reference_um[1],
            "day_reference_um must be an increasing (start, stop) range from 0",
        )
        fail_unless(
            self.methods and set(self.methods) <= set(METHODS), f"methods must come from {METHODS}"
        )
        fail_unless(self.analysis_method in self.methods, "analysis_method must be one of methods")
        fail_unless(
            all(v >= 1 and int(v) == v for v in self.protein_channels.values()),
            "protein_channels must map to channel numbers counting from 1",
        )
        fail_unless(set(self.plot_channels) <= self.protein_channels.keys(), "Unknown plot channel")
        fail_unless(set(self.plot_compartments) <= set(COMPARTMENTS), "Unknown plot compartment")
        fail_unless(
            set(self.phase_compartments) <= self.protein_channels.keys()
            and set(self.phase_compartments.values()) <= set(COMPARTMENTS),
            "phase_compartments must map protein channels to compartments",
        )
        fail_unless(
            self.phase_labels and len(set(self.phase_labels)) == len(self.phase_labels),
            "phase_labels must be a non-empty list without repeats",
        )
        fail_unless(
            self.phase_colormap in matplotlib.colormaps,
            f"phase_colormap: unknown Matplotlib colormap {self.phase_colormap!r}",
        )
        fail_unless(
            len(self.phase_colormap_range) == 2 and all(0 <= v <= 1 for v in self.phase_colormap_range),
            "phase_colormap_range must be two values within [0, 1]",
        )
        return self


COMPARTMENTS = ("nuclear", "cytoplasmic", "all")
METHODS = ("per_slice", "pooled", "block_mip", "full_mip")


def parameter_table(config=None):
    """Table of every setting with its group, value, allowed choices and description."""
    c = config or Config()
    return pd.DataFrame(
        [
            dict(
                group=f.metadata.get("group"),
                parameter=f.name,
                value=getattr(c, f.name),
                choices=f.metadata.get("choices"),
                explanation=f.metadata.get("description"),
            )
            for f in fields(c)
        ]
    )


def jsonable(obj):
    """Convert the paths and NumPy values json.dumps cannot handle by itself."""
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(type(obj).__name__)


def write_json(path, data):
    """Save data as indented JSON, creating the folder if needed."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, default=jsonable))


def file_identity(path):
    """Path, size and modification time of a file, recorded as provenance.

    Reading the size and modification time is much faster than hashing a
    multi-GB image.
    """
    p = Path(path).resolve()
    info = p.stat()
    return dict(path=str(p), bytes=info.st_size, modified_ns=info.st_mtime_ns)


def git_state(folder=None):
    """The git commit the code is at, and the tracked files changed since that commit.

    Both are None when git is not installed or the folder is not a git clone.
    """
    folder = Path(folder or Path(__file__).parent)

    def git(*args):
        return subprocess.run(
            ["git", "-C", str(folder), *args], capture_output=True, text=True, check=True, timeout=10
        ).stdout

    try:
        commit = git("rev-parse", "HEAD").strip()
        changed = [line[3:] for line in git("status", "--porcelain", "--untracked-files=no").splitlines()]
    except (OSError, subprocess.SubprocessError):
        return dict(git_commit=None, git_uncommitted=None)
    return dict(git_commit=commit, git_uncommitted=changed)


def version_info():
    """Versions of the main packages, the git commit of the code and a hash of this file."""
    result = {}
    for name in ("numpy", "pandas", "scipy", "scikit-image", "tifffile", "roifile", "torch", "cellpose"):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = None
    result.update(git_state())
    result["pipeline_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result


def _zcyx(a, axes):
    """Reorder an image to Z, C, Y, X, adding Z or C if absent and dropping length-1 axes."""
    for axis in list(axes):
        if axis not in "ZCYX":
            i = axes.index(axis)
            if a.shape[i] != 1:
                raise ValueError(f"Unsupported non-singleton axis: {axes}")
            a = np.squeeze(a, axis=i)
            axes = axes[:i] + axes[i + 1 :]
    for axis in "ZC":
        if axis not in axes:
            a = np.expand_dims(a, 0)
            axes = axis + axes
    if set(axes) != set("ZCYX"):
        raise ValueError(f"Expected ZCYX axes, got {axes}")
    return np.asarray(a).transpose(*(axes.index(x) for x in "ZCYX"))


def read_image(path, position=1, config=None):
    """Read one image (one position of an ND2) as a ZCYX array, plus its calibration."""
    p = Path(path)
    start = time.perf_counter()
    if p.suffix.lower() == ".nd2":
        import nd2

        with nd2.ND2File(p) as f:
            sizes = dict(f.sizes)
            axes = "".join(sizes)
            _log(f"Reading {p.name}, position {position} of {sizes.get('P', 1)} ({sizes})")
            a = f.to_dask()
            if "P" in axes:
                if not 1 <= position <= sizes["P"]:
                    raise ValueError(f"{p.name} has no position {position}")
                index = [slice(None)] * len(axes)
                index[axes.index("P")] = position - 1
                a = a[tuple(index)]
                axes = axes.replace("P", "")
            elif position != 1:
                raise ValueError(f"{p.name} has only one position")
            if PROGRESS:
                from tqdm.dask import TqdmCallback

                with TqdmCallback(desc=f"Reading {p.name}", leave=False, unit="plane"):
                    a = a.compute()
            else:
                a = a.compute()
            a = _zcyx(a, axes)
            voxel = f.voxel_size()
            xy = float(voxel.x)
            z = float(voxel.z)
    else:
        if position != 1:
            raise ValueError("A TIFF holds one position, so position must be 1")
        _log(f"Reading {p.name} ({p.stat().st_size / 1024**3:.2f} GB)")
        with tifffile.TiffFile(p) as f:
            a = _zcyx(f.asarray(), f.series[0].axes)
            meta = f.imagej_metadata or {}
            xy = None
            z = meta.get("spacing")
            if f.ome_metadata:
                from xml.etree import ElementTree as ET

                pix = ET.fromstring(f.ome_metadata).find(".//{*}Pixels")
                if pix is not None:
                    unit = pix.get("PhysicalSizeXUnit", "µm")
                    if pix.get("PhysicalSizeX") and unit in ("µm", "um", "micrometer"):
                        xy = float(pix.get("PhysicalSizeX"))
                        z = float(pix.get("PhysicalSizeZ", z or 1))
            unit = meta.get("unit", "").lower().replace("μ", "µ")
            res = f.pages[0].tags.get("XResolution")
            if xy is None and res and unit in ("um", "µm", "micron", "microns", "micrometer"):
                xy = res.value[1] / res.value[0]
    if config:
        xy = config.pixel_size_um if config.pixel_size_um is not None else xy
        z = config.z_step_um if config.z_step_um is not None else z
    if xy is None or xy <= 0:
        raise ValueError(f"{p} has no pixel size in its metadata. Set pixel_size_um.")
    _log(f"  read {p.name}: ZCYX {tuple(a.shape)} in {_since(start)}")
    return a, dict(
        pixel_size_um=float(xy), z_step_um=float(z) if z else None, axes="ZCYX", shape=list(a.shape)
    )


def _canonical_stem(stem, config):
    """Filename stem with the denoised_tokens removed, i.e. the original image's name."""
    for token in config.denoised_tokens:
        stem = re.sub(re.escape(token), "", stem, flags=re.I)
    return stem


def _is_denoised(path, config):
    """Whether a filename marks an NIS-denoised export."""
    return _canonical_stem(Path(path).stem, config) != Path(path).stem


def sample_dir(config, sample, kind):
    """Output folder output_root/<kind>/<replicate>/<sample> for one sample."""
    return Path(config.output_root) / kind / sample["replicate"] / sample["sample_id"]


def _roi_sources(sample, config):
    """ROI files of a sample, in order of preference.

    The converted CSV under output_root comes first. Otherwise a rois.csv, or
    else PSM_line_side*.zip files, in the folder named after the image file.
    """
    canonical = sample_dir(config, sample, "rois") / "rois.csv"
    if canonical.exists():
        return [canonical]
    folders = [Path(sample["raw_path"]).with_suffix("")]
    if sample.get("nis_path"):
        folders.append(Path(sample["nis_path"]).with_suffix(""))
    for folder in folders:
        if (folder / "rois.csv").exists():
            return [folder / "rois.csv"]
    for folder in folders:
        paths = sorted(folder.glob("PSM_line_side*.zip"))
        if paths:
            return paths
    return []


def _acquisition_day(nd2_file):
    """The ND2's acquisition date as YYYY-MM-DD, or "" if the file has none."""
    # NIS-Elements writes e.g. "10/6/2026  5:04:28 PM" (month first, US locale).
    info = getattr(nd2_file, "text_info", None) or {}
    stamp = pd.to_datetime(info.get("date", ""), errors="coerce")
    return "" if pd.isna(stamp) else stamp.strftime("%Y-%m-%d")


def discover_samples(config):
    """Find every sample in the N_* replicate folders and return one row per sample.

    A sample is one embryo image: an exported TIFF, or one position of an ND2
    file. When the same position exists as both, the TIFF is used. Each row also
    records the matching NIS-denoised file, if any, the ROI files, the imaging
    day and the phase. The notebook shows this table before anything runs.
    """
    c = config.validate()
    root = Path(c.experimental_data).resolve()
    output = Path(c.output_root).resolve()
    discovery_start = time.perf_counter()
    _log(f"Searching {root} for N_* replicate folders (slow on a network drive)...")
    folders = [
        p
        for p in [root, *root.rglob("N_*")]
        if p.is_dir()
        and re.match(r"^N_\d+(?:_|$)", p.name)
        and p.resolve() != output
        and output not in p.resolve().parents
    ]
    # Folders inside a replicate can also be named N_2 (for example rois/N_2),
    # so only the outermost matching folder counts as the replicate.
    replicate_paths = {p.resolve() for p in folders}
    folders = [
        p for p in folders if not any(parent in replicate_paths for parent in p.resolve().parents)
    ]
    by_rep = {}
    for folder in folders:
        rep = "N_" + re.match(r"^N_(\d+)", folder.name)[1]
        if c.replicates and rep not in c.replicates:
            continue
        if rep in by_rep:
            raise ValueError(
                f"Found two {rep} folders under {root}: {by_rep[rep]} and {folder}. "
                "Point experimental_data at a single experiment."
            )
        by_rep[rep] = folder
    _log(f"Found {len(by_rep)} replicate(s): {', '.join(sorted(by_rep)) or 'none'}")
    result = []
    for rep, folder in sorted(by_rep.items(), key=lambda kv: int(kv[0][2:])):
        _log(f"{rep}: listing {c.image_glob} files under {folder.name}...")
        all_tiffs = sorted(
            p
            for p in folder.rglob(c.image_glob)
            if Path(c.output_root).resolve() not in p.resolve().parents
        )
        records = {}
        nis_tiffs = {}
        for p in _bar(all_tiffs, desc=f"{rep}: inspecting TIFFs", unit="file"):
            if _is_denoised(p, c):
                key = _canonical_stem(p.stem, c)
                if key in nis_tiffs:
                    raise ValueError(f"Ambiguous NIS counterpart for {rep}/{key}")
                nis_tiffs[key] = p
                continue
            if p.stem in ("Mask", "Cyto_Mask", "labels", "processed_dapi"):
                continue
            with tifffile.TiffFile(p) as f:
                s = f.series[0]
                channels = s.shape[s.axes.index("C")] if "C" in s.axes else 1
            if channels < max(c.protein_channels.values()):
                continue
            if p.stem in records:
                raise ValueError(f"Duplicate original TIFF for {rep}/{p.stem}")
            records[p.stem] = dict(raw_path=str(p), raw_position=1, nis_path="", nis_position=1)
        if c.include_nd2:
            import nd2

            nds = sorted(folder.glob("*.nd2"))
            denoised = {_canonical_stem(p.stem, c): p for p in nds if _is_denoised(p, c)}
            for p in _bar(nds, desc=f"{rep}: inspecting ND2 files", unit="file"):
                if _is_denoised(p, c):
                    continue
                with nd2.ND2File(p) as f:
                    npos = f.sizes.get("P", 1)
                    day = _acquisition_day(f)
                for position in range(1, npos + 1):
                    key = f"{p.stem}_P{position}"
                    record = records.setdefault(
                        key,
                        dict(raw_path=str(p), raw_position=position, nis_path="", nis_position=1),
                    )
                    record["imaging_day"] = day  # same acquisition as a matching TIFF
                    # The DAPI source is chosen separately from the protein source, so a
                    # TIFF original can still take its DAPI from the denoised ND2.
                    if p.stem in denoised and c.nis_preference == "nd2":
                        record.update(nis_path=str(denoised[p.stem]), nis_position=position)
                    elif key in nis_tiffs:
                        record["nis_path"] = str(nis_tiffs[key])
                    elif p.stem in denoised:
                        record.update(nis_path=str(denoised[p.stem]), nis_position=position)
        undated = []
        for key, record in records.items():
            if key in nis_tiffs and not record["nis_path"]:
                record.update(nis_path=str(nis_tiffs[key]), nis_position=1)
            # Without an acquisition date, the replicate stands in for the imaging day.
            if not record.get("imaging_day"):
                record["imaging_day"] = rep
                undated.append(key)
            s = dict(
                replicate=rep,
                replicate_path=str(folder),
                sample_id=key,
                embryo_id=f"{rep}/{key}",
                **record,
            )
            sources = _roi_sources(s, c)
            s.update(roi_sources=[str(p) for p in sources], has_rois=bool(sources))
            result.append(s)
        found = [r for r in result if r["replicate"] == rep]
        _log(
            f"{rep}: {len(found)} image(s), {sum(r['has_rois'] for r in found)} with ROIs, "
            f"{sum(bool(r['nis_path']) for r in found)} with a denoised NIS file; "
            f"imaging day(s) {', '.join(sorted({r['imaging_day'] for r in found}))}"
        )
        if undated:
            _log(f"{rep}: WARNING {undated} have no acquisition date, so their imaging day is set to '{rep}'")
    if not result:
        raise ValueError("No original multichannel images found in N_* folders")
    _log(f"Sample discovery finished in {_since(discovery_start)}")
    return attach_phases(pd.DataFrame(result), c)


def attach_phases(samples, config):
    """Add each sample's phase from its replicate's phase table.

    Each replicate uses one table: phase_files[replicate] if given, otherwise
    phases.csv in the replicate folder. phase_overrides are applied last.
    """
    table = samples.copy()
    table["phase"] = ""
    table["phase_source"] = ""
    for rep, group in table.groupby("replicate"):
        folder = Path(group.iloc[0]["replicate_path"])
        explicit = config.phase_files.get(rep)
        path = Path(explicit) if explicit else folder / "phases.csv"
        alternatives = sorted(folder.glob("*phase*.csv"))
        if not path.exists():
            if alternatives and not explicit:
                warnings.warn(
                    f"{rep} has no phases.csv but has {[p.name for p in alternatives]}. Choose one with phase_files."
                )
            if explicit:
                raise FileNotFoundError(path)
            continue
        _log(f"{rep}: phase annotations from {path.name}")
        annotations = pd.read_csv(path, dtype=str, keep_default_na=False)
        if "phase" not in annotations:
            raise ValueError(f"No phase column: {path}")
        if "sample_id" in annotations:
            if "replicate" in annotations:
                annotations = annotations[annotations.replicate == rep]
            keys = annotations.sample_id
        elif "file" in annotations:

            def key(x):
                match = re.fullmatch(r"(.+)\.nd2(?:\s+(P\d+))?", x.strip(), re.I)
                return f"{match[1]}_{match[2] or 'P1'}" if match else Path(x).stem

            keys = annotations.file.map(key)
        else:
            raise ValueError(f"Phase CSV needs sample_id or file column: {path}")
        if keys.duplicated().any():
            raise ValueError(f"Duplicate phase assignments in {path}")
        mapping = dict(zip(keys, annotations.phase))
        for i, row in group.iterrows():
            table.loc[i, "phase"] = mapping.get(row.sample_id, "")
            table.loc[i, "phase_source"] = str(path)
        unknown = set(mapping) - set(group.sample_id)
        if unknown:
            warnings.warn(f"{rep}: {len(unknown)} rows of the phase table match no sample: {sorted(unknown)}")
    for i, row in table.iterrows():
        key = f"{row.replicate}/{row.sample_id}"
        if key in config.phase_overrides:
            table.loc[i, "phase"] = str(config.phase_overrides[key])
            table.loc[i, "phase_source"] = "config.phase_overrides"
    return table


@dataclass
class ROI:
    roi_id: str
    xy: np.ndarray
    width_px: float

    def validate(self):
        self.xy = np.asarray(self.xy, dtype=float)
        if (
            self.xy.ndim != 2
            or self.xy.shape[1] != 2
            or len(self.xy) < 2
            or not np.isfinite(self.xy).all()
        ):
            raise ValueError("ROI vertices must be finite XY pairs")
        keep = np.r_[True, np.linalg.norm(np.diff(self.xy, axis=0), axis=1) > 0]
        self.xy = self.xy[keep]
        if len(self.xy) < 2 or self.width_px < 1 or not np.isfinite(self.width_px):
            raise ValueError("ROI needs nonzero length and width >=1 pixel")
        if not re.fullmatch(r"[A-Za-z0-9_-]+", self.roi_id):
            raise ValueError("ROI IDs use letters, digits, underscores and hyphens")
        return self


def load_rois(paths):
    """Read ROIs from a vertex CSV or from ImageJ line-ROI ZIP files."""
    if isinstance(paths, (str, Path)):
        paths = [paths]
    result = []
    for path in paths:
        p = Path(path)
        if p.suffix.lower() == ".csv":
            t = pd.read_csv(p)
            if not {"roi_id", "vertex_index", "x_px", "y_px", "width_px"} <= set(t):
                raise ValueError(f"{p} needs the columns roi_id, vertex_index, x_px, y_px and width_px")
            for name, g in t.groupby("roi_id", sort=False):
                g = g.sort_values("vertex_index")
                if g.vertex_index.duplicated().any() or not np.array_equal(
                    g.vertex_index.to_numpy(), np.arange(len(g))
                ):
                    raise ValueError(
                        "ROI vertex indices must be unique, contiguous, and start at 0"
                    )
                if g.width_px.nunique() != 1:
                    raise ValueError("Width must be constant within each ROI")
                result.append(
                    ROI(
                        str(name), g[["x_px", "y_px"]].to_numpy(), float(g.width_px.iloc[0])
                    ).validate()
                )
        else:
            import roifile

            items = roifile.roiread(p)
            items = items if isinstance(items, list) else [items]
            for j, r in enumerate(items):
                if r.roitype not in (roifile.ROI_TYPE.POLYLINE, roifile.ROI_TYPE.LINE):
                    raise ValueError("Only line/polyline ROIs are supported")
                match = re.search(r"side\d+", p.stem)
                name = match[0] if match else (r.name or f"side{j}")
                result.append(
                    ROI(
                        name, r.coordinates(), r.float_stroke_width or r.stroke_width or 1
                    ).validate()
                )
    if len({r.roi_id for r in result}) != len(result):
        raise ValueError("Duplicate ROI IDs")
    return result


def save_rois(path, rois):
    """Save ROIs as a vertex CSV.

    The CSV is written under a temporary name and then renamed, so the file on
    disk holds either the old or the new ROIs even if writing is interrupted.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for r in rois:
        r.validate()
        records.extend(
            dict(roi_id=r.roi_id, vertex_index=i, x_px=x, y_px=y, width_px=r.width_px)
            for i, (x, y) in enumerate(r.xy)
        )
    if not records:
        raise ValueError("No ROI vertices to save")
    if len({r.roi_id for r in rois}) != len(rois):
        raise ValueError("Duplicate ROI IDs")
    with tempfile.NamedTemporaryFile(
        mode="w", dir=p.parent, suffix=".csv", delete=False, newline=""
    ) as f:
        temporary = Path(f.name)
        pd.DataFrame(records).to_csv(f, index=False, float_format="%.9g")
    temporary.replace(p)


def convert_rois(samples, config, overwrite=False):
    """Copy each sample's ROIs to output_root/rois/<replicate>/<sample>/rois.csv.

    A CSV that already exists is kept, so that edits made in Fiji are preserved,
    unless overwrite=True.
    """
    records = []
    for s in _bar(samples.to_dict("records"), desc="Converting ROIs", unit="image"):
        path = sample_dir(config, s, "rois") / "rois.csv"
        sources = _roi_sources(s, config)
        if path.exists() and not overwrite:
            status = "already present"
        elif sources:
            save_rois(path, load_rois(sources))
            status = "converted"
        else:
            status = "needs drawing"
        records.append(
            dict(
                replicate=s["replicate"],
                sample_id=s["sample_id"],
                roi_path=str(path),
                status=status,
            )
        )
    counts = pd.Series([r["status"] for r in records]).value_counts().to_dict()
    _log(f"ROI conversion: {counts}")
    return pd.DataFrame(records)


def centerline(roi):
    """Smooth centerline through the ROI's vertices, one point per pixel of length.

    This is the spline ImageJ fits to a polyline (a natural cubic spline whose
    knots are spaced by the square root of the distance between vertices), so the
    measured band follows the same curve as the line drawn in Fiji.
    """
    xy = roi.validate().xy.astype(np.float32)
    origin = np.floor(xy.min(axis=0))
    points = xy - origin
    delta = np.diff(points, axis=0)
    chord = np.sqrt(np.sum(delta.astype(float) ** 2, axis=1))
    increments = np.sqrt(np.sqrt(np.sum(delta**2, axis=1).astype(float))).astype(np.float32)
    nodes = np.r_[
        np.float32(0), np.cumsum(np.maximum(increments, np.float32(0.001)), dtype=np.float32)
    ]
    n = len(nodes)
    second = np.zeros((n, 2))
    upper = np.zeros(n)
    right = float(nodes[1] - nodes[0])
    slope = (points[1] - points[0]).astype(float) / right
    for j in range(1, n - 1):
        left, old = right, slope
        right = float(nodes[j + 1]) - float(nodes[j])
        slope = (points[j + 1].astype(float) - points[j].astype(float)) / right
        inv = 1 / (2 * (right + left) - left * upper[j - 1])
        second[j] = (6 * (slope - old) - left * second[j - 1]) * inv
        upper[j] = right * inv
    for j in range(n - 2, -1, -1):
        second[j] -= upper[j] * second[j + 1]
    count = int(chord.sum() * 2) + 1
    q = np.arange(count) * (float(nodes[-1]) / (count - 1))
    hi = np.clip(np.searchsorted(nodes, q, side="left"), 1, n - 1)
    lo = hi - 1
    span = (nodes[hi] - nodes[lo]).astype(float)[:, None]
    a = (q - nodes[lo])[:, None]
    b = (nodes[hi] - q)[:, None]
    fine = (
        (
            a * points[hi]
            + b * points[lo]
            + ((a * a - span * span) * a * second[hi] + (b * b - span * span) * b * second[lo]) / 6
        )
        / span
    ).astype(np.float32)
    length = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(fine.astype(float), axis=0), axis=1))]
    target = np.linspace(0, length[-1], int(np.floor(length[-1] + 0.5)) + 1)
    return (
        np.column_stack([np.interp(target, length, fine[:, d]) for d in (0, 1)]).astype(np.float32)
        + origin
    )


def band_grid(roi, width_override=None):
    """Sampling points of the band around an ROI, and the distance along it in pixels.

    The grid has one row per pixel of band width and one column per pixel of
    length. Each row is the centerline shifted sideways along its normal.
    """
    xy = centerline(roi).astype(float)
    delta = np.vstack([xy[1] - xy[0], np.diff(xy, axis=0)])
    normal = (
        np.column_stack([-delta[:, 1], delta[:, 0]])
        / np.linalg.norm(delta, axis=1).astype(np.float32).astype(float)[:, None]
    )
    width = int(round(width_override if width_override is not None else roi.width_px))
    offsets = np.arange(width) - (width - 1) / 2
    grid = xy[None, :, :] + offsets[:, None, None] * normal[None, :, :]
    distance = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
    return grid, distance


def write_mask(path, mask, pixel_size_um, z_step_um=None):
    """Save a mask as a 16-bit ZYX ImageJ TIFF (0 outside, 4095 inside) that Fiji opens."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(axes="ZYX", unit="um")
    if z_step_um:
        metadata["spacing"] = z_step_um
    a = (np.asarray(mask) > 0).astype(np.uint16) * 4095
    tifffile.imwrite(
        p, a, imagej=True, resolution=(1 / pixel_size_um, 1 / pixel_size_um), metadata=metadata
    )


def _parallel(fn, a, workers, desc=None):
    """Apply fn to each Z-plane in parallel threads and stack the results."""
    with ThreadPoolExecutor(workers) as pool:
        results = pool.map(fn, a)
        if desc:
            results = _bar(results, desc=desc, total=len(a), unit="plane")
        return np.stack(list(results))


def choose_device(requested):
    """Torch device name for the device setting; 'auto' picks the fastest one present."""
    import torch

    available = {"cuda": torch.cuda.is_available(), "mps": torch.backends.mps.is_available()}
    if requested == "auto":
        return next((name for name, ok in available.items() if ok), "cpu")
    if not available.get(requested, True):
        raise RuntimeError(f"{requested} is not available in this Python environment")
    return requested


def _contrast(a, c):
    """CLAHE on each Z-plane, returning values from 0 to 1."""
    import cv2  # installed with cellpose

    def one(x):
        low, high = np.percentile(x, [1, 99])
        if high <= low:
            return np.zeros_like(x, dtype=np.float32)
        # OpenCV's CLAHE needs integers, so the 1st-99th percentile range is stretched
        # to 0-255. With 8 bits each tile's histogram has 256 bins, which fill well
        # enough for the clip limit to work. With 16 bits the bins are nearly empty,
        # hardly anything is clipped, and the threshold then over-segments.
        x = np.round(np.clip((x - low) / (high - low), 0, 1) * 255).astype(np.uint8)
        # OpenCV divides the plane into a grid of tiles. The grid is sized so each
        # tile is about clahe_window_px wide.
        h, w = x.shape
        grid = (max(1, round(w / c.clahe_window_px)), max(1, round(h / c.clahe_window_px)))
        clahe = cv2.createCLAHE(clipLimit=c.clahe_clip, tileGridSize=grid)
        return clahe.apply(x).astype(np.float32) / 255

    return _parallel(one, a, c.threads_per_worker, desc="DAPI: CLAHE contrast")


def _median(a, size, c, desc):
    """Square median filter per Z-plane, edges padded by repeating the border pixel."""
    import cv2  # installed with cellpose

    def one(x):
        # OpenCV gives the same result as SciPy about 100 times faster, but it only
        # supports sizes 3 and 5 on 16-bit or float32 images.
        if size in (3, 5) and x.dtype in (np.uint16, np.float32):
            return cv2.medianBlur(x, size)
        return ndi.median_filter(x, size, mode="nearest")

    return _parallel(one, a, c.threads_per_worker, desc=desc)


def _filter_dapi(dapi, c):
    """Prepare raw DAPI for segmentation: median filter, background subtraction, CLAHE."""
    import cv2  # installed with cellpose

    _log("DAPI: median filter, background subtraction and CLAHE")
    a = dapi
    if c.dapi_median_size:
        a = _median(a, c.dapi_median_size, c, "DAPI: median filter")
    if c.background_size_px:
        # An opening (erosion, then dilation) with a square wider than a nucleus
        # removes the nuclei and keeps the smooth background beneath them.
        # Subtracting it leaves the nuclei on a flat background. OpenCV gives the
        # same result as scipy's grey_opening about 6 times faster.
        size = (c.background_size_px, c.background_size_px)
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT, size)
        a = _parallel(
            lambda x: np.maximum(x - cv2.morphologyEx(x, cv2.MORPH_OPEN, kernel), 0),
            a,
            c.threads_per_worker,
            desc="DAPI: background subtraction",
        )
    return _contrast(a, c)


@lru_cache(maxsize=2)
def _cellpose_model(model_name, device, models_dir):
    """Load a Cellpose 4 model once per session (the first use downloads its weights)."""
    import torch

    os.environ["CELLPOSE_LOCAL_MODELS_PATH"] = models_dir
    os.environ.setdefault("MPLCONFIGDIR", str(Path(models_dir).parent / "runtime" / "matplotlib"))
    from cellpose import models

    models.MODEL_DIR = Path(models_dir)
    models.MODEL_DIR.mkdir(parents=True, exist_ok=True)
    if not Path(model_name).is_file() and model_name not in models.MODEL_NAMES:
        raise ValueError(f"Unknown Cellpose model, or missing weights file: {model_name}")
    start = time.perf_counter()
    _log(f"Loading Cellpose model '{model_name}' on {device}")
    # bfloat16 is Cellpose 4's default. It roughly halves the GPU memory the model
    # needs and is faster on recent GPUs.
    model = models.CellposeModel(
        pretrained_model=model_name,
        device=torch.device(device),
        gpu=device != "cpu",
        use_bfloat16=True,
    )
    _log(f"  model ready in {_since(start)}")
    return model


def _safe_output(sample, c):
    """Stop if output_root is inside the sample's image folder, then create it."""
    root = Path(c.output_root).resolve()
    raw = Path(sample["raw_path"]).resolve()
    if root == raw.parent or raw.parent in root.parents:
        raise ValueError("output_root must be outside the folder that holds the images")
    root.mkdir(parents=True, exist_ok=True)


def mask_kind(c):
    """Name of the folder for computed masks, e.g. 'masks_cellpose_denoised'.

    Each combination of segmentation method and DAPI source has its own folder,
    so masks made with different settings are stored separately.
    """
    method = "cellpose" if c.use_cellpose else "threshold"
    source = "denoised" if c.use_ai_denoise else "raw"
    return f"masks_{method}_{source}"


def _read_mask(path):
    """Read a mask TIFF as a boolean ZYX array; any value above 0 is inside."""
    mask = tifffile.imread(path) > 0
    return mask[None] if mask.ndim == 2 else mask


# Settings that change a computed mask. They are saved with it in mask.json.
MASK_SETTINGS = (
    "use_ai_denoise", "use_cellpose", "dapi_channel", "dapi_median_size", "background_size_px", "clahe_window_px", "clahe_clip",
    "cellpose_model", "cellpose_diameter_px", "cellpose_percentiles", "cellprob_threshold",
    "flow_threshold", "min_size_px", "max_size_fraction", "flow_iterations", "tile_overlap",
    "augment",
)  # fmt: skip


def segment_sample(sample, config):
    """Return one sample's nuclear mask (boolean ZYX), its calibration and its provenance.

    With mask_source='compute', the mask is made from the sample's DAPI and saved
    in <mask_kind>/<replicate>/<sample>/, or read from there when reuse_cache is
    True and a saved mask exists. The other mask sources read a mask made
    outside the pipeline or corrected in Fiji.
    """
    c = config.validate()
    s = dict(sample)
    _safe_output(s, c)
    _log(f"Mask: {s['replicate']}/{s['sample_id']} (mask_source='{c.mask_source}')")
    folder = sample_dir(c, s, mask_kind(c))
    meta_path, mask_path = folder / "mask.json", folder / "Mask.tif"
    if c.mask_source == "compute" and c.reuse_cache and meta_path.exists() and mask_path.exists():
        _log(f"Mask: using the saved {mask_path}")
        saved = json.loads(meta_path.read_text())
        return _read_mask(mask_path), saved["provenance"]["calibration"], dict(saved, cached=True)
    if c.mask_source == "compute" and c.use_ai_denoise and not s.get("nis_path"):
        raise FileNotFoundError(
            f"{s['replicate']}/{s['sample_id']} has no NIS-denoised file. Export one from "
            "NIS-Elements, or set use_ai_denoise=False to filter the raw DAPI instead."
        )
    raw, cal = read_image(s["raw_path"], int(s.get("raw_position", 1)), c)
    if c.dapi_channel > raw.shape[1]:
        raise ValueError(f"The image has no channel {c.dapi_channel} for DAPI")
    if c.mask_source != "compute":
        if c.mask_source == "edited":
            paths = [sample_dir(c, s, "edited_masks") / "Mask.tif"]
        else:
            paths = [Path(s["raw_path"]).with_suffix("") / "Mask.tif"]
            if s.get("nis_path"):
                paths.append(Path(s["nis_path"]).with_suffix("") / "Mask.tif")
        matches = [p for p in paths if p.exists()]
        if not matches:
            raise FileNotFoundError(f"No {c.mask_source} mask for {s['sample_id']}")
        _log(f"Mask: reading {matches[0]}")
        mask = _read_mask(matches[0])
        if mask.shape != (raw.shape[0], *raw.shape[-2:]):
            raise ValueError(f"{matches[0]} does not have the image's Z, Y, X shape")
        return mask, cal, dict(source=c.mask_source, file=file_identity(matches[0]))

    device = choose_device(c.device)
    import torch

    torch.set_num_threads(c.threads_per_worker)
    folder.mkdir(parents=True, exist_ok=True)
    stamp = dict(
        original=file_identity(s["raw_path"]),
        position=int(s.get("raw_position", 1)),
        settings={k: getattr(c, k) for k in MASK_SETTINGS},
        calibration=cal,
        versions=version_info(),
        device=device,
    )
    if c.use_ai_denoise:
        stamp.update(
            source="NIS-denoised DAPI",
            nis=file_identity(s["nis_path"]),
            nis_position=int(s.get("nis_position", 1)),
        )
    else:
        stamp["source"] = "filtered raw DAPI"
    method = "Cellpose" if c.use_cellpose else "IsoData threshold"
    _log(f"Mask: making it from {stamp['source']} with {method} on {device}")
    start = time.perf_counter()
    if c.use_ai_denoise:
        a, nis_cal = read_image(s["nis_path"], stamp["nis_position"], c)
        ch = 0 if a.shape[1] == 1 else c.dapi_channel - 1
        if ch >= a.shape[1] or a.shape[0] != raw.shape[0] or a.shape[-2:] != raw.shape[-2:]:
            raise ValueError("The NIS export's size or channels do not match the original image")
        if not np.isclose(nis_cal["pixel_size_um"], cal["pixel_size_um"], rtol=1e-4):
            raise ValueError("The NIS export's pixel size differs from the original image's")
        processed = a[:, ch].astype(np.float32)
    else:
        processed = _filter_dapi(raw[:, c.dapi_channel - 1], c)
    if c.use_cellpose:
        model = _cellpose_model(
            c.cellpose_model, device, str((Path(c.output_root) / "models").resolve())
        )
        stamp["segmentation_weights"] = file_identity(model.pretrained_model)
        labels = []
        for a in _bar(processed, desc="Mask: Cellpose segmentation", unit="plane"):
            lab, _, _ = model.eval(
                a,
                batch_size=c.batch_size,
                diameter=c.cellpose_diameter_px,
                normalize={"normalize": True, "percentile": list(c.cellpose_percentiles)},
                flow_threshold=c.flow_threshold,
                cellprob_threshold=c.cellprob_threshold,
                min_size=c.min_size_px,
                max_size_fraction=c.max_size_fraction,
                niter=c.flow_iterations,
                augment=c.augment,
                tile_overlap=c.tile_overlap,
                resample=True,
                do_3D=False,
            )
            labels.append(lab.astype(np.uint32))
        labels = np.stack(labels)
        mask = labels > 0
        tifffile.imwrite(folder / "labels.tif", labels, metadata={"axes": "ZYX"})
    else:
        mask = _parallel(
            lambda x: x > threshold_isodata(x) if np.ptp(x) > 0 else np.zeros(x.shape, bool),
            processed,
            c.threads_per_worker,
            desc="Mask: IsoData threshold",
        )
    _log(f"Mask: made in {_since(start)}; saving Mask.tif, Cyto_Mask.tif and processed_dapi.tif")
    write_mask(mask_path, mask, cal["pixel_size_um"], cal["z_step_um"])
    write_mask(folder / "Cyto_Mask.tif", ~mask, cal["pixel_size_um"], cal["z_step_um"])
    tifffile.imwrite(folder / "processed_dapi.tif", processed, metadata={"axes": "ZYX"})
    saved = dict(provenance=stamp, seconds=time.perf_counter() - start, cached=False)
    write_json(meta_path, saved)
    return mask, cal, saved


def _protein_filter(a, c, name="protein"):
    """Median filter a protein channel, then replace bright outlier pixels."""
    if c.protein_median_size:
        a = _median(a, c.protein_median_size, c, f"{name}: median filter")
    if c.protein_outliers:
        r = int(np.ceil(c.outlier_radius_px))
        y, x = np.mgrid[-r : r + 1, -r : r + 1]
        footprint = x * x + y * y <= c.outlier_radius_px**2 + 1

        def clean(x):
            m = ndi.median_filter(x, footprint=footprint, mode="nearest")
            return np.where(x.astype(float) - m > c.outlier_threshold, m, x)

        a = _parallel(clean, a, c.threads_per_worker, desc=f"{name}: outlier removal")
    return a.astype(np.float32)


def _windows(nz, c, method):
    """(first, stop) Z-plane ranges, counting from 0, that a Z method measures."""
    first = c.z_start - 1
    last = c.z_stop if c.z_stop is not None else nz
    if not 0 <= first < last <= nz:
        raise ValueError("z_start and z_stop lie outside the image's Z-planes")
    if method == "per_slice":
        return [(z, z + 1) for z in range(first, last)]
    if method in ("pooled", "full_mip"):
        return [(first, last)]
    start = max(first, c.mip_start - 1)
    out = []
    for z in range(start, last, c.mip_block_size):
        stop = min(z + c.mip_block_size, last)
        if stop - z < c.mip_block_size:
            if c.incomplete_block == "error":
                raise ValueError("Incomplete final MIP block")
            if c.incomplete_block == "drop":
                continue
        out.append((z, stop))
    if not out:
        raise ValueError("No complete block MIP windows")
    return out


def _band_sampler(grid, h, w, order):
    """Pixel indices and weights that sample an image along the band.

    Linear interpolation of a point is a weighted sum of its 4 surrounding
    pixels, so sampling is (image.ravel()[idx] * wts).sum(-1). Computing idx and
    wts once lets every slice, channel and mask reuse them. This reproduces
    ndi.map_coordinates(order, mode="constant", cval=0): points outside the
    image get weight 0. order=0 takes the nearest pixel instead.
    """
    r, c = grid[..., 1].ravel(), grid[..., 0].ravel()
    inside = (r >= 0) & (r <= h - 1) & (c >= 0) & (c <= w - 1)
    if order == 0:
        rn = np.clip(np.floor(r + 0.5), 0, h - 1).astype(np.int64)
        cn = np.clip(np.floor(c + 0.5), 0, w - 1).astype(np.int64)
        return (rn * w + cn)[:, None], inside[:, None].astype(float)
    # Top-left neighbour. Clipping keeps the +1 neighbours inside the image. At the
    # edge their weight is 0, so the clipped index does not change the result.
    r0 = np.clip(np.floor(r), 0, max(h - 2, 0)).astype(np.int64)
    c0 = np.clip(np.floor(c), 0, max(w - 2, 0)).astype(np.int64)
    fr, fc = r - r0, c - c0
    idx = np.stack([r0 * w + c0, r0 * w + c0 + 1, (r0 + 1) * w + c0, (r0 + 1) * w + c0 + 1], -1)
    wts = np.stack([(1 - fr) * (1 - fc), (1 - fr) * fc, fr * (1 - fc), fr * fc], -1)
    wts[~inside] = 0
    return idx, wts


def measure_arrays(channels, mask, rois, cal, config):
    """Measure each protein along each ROI band and return one row per position.

    At each position along a band, the masked image and the mask itself are
    sampled at the same points and summed across the band's width. Their ratio
    is the mean intensity inside the compartment (masked_mean). Dividing by the
    band width instead gives band_mean. The pooled method adds the sums of all
    Z-planes before dividing. Maximum projections are taken after masking, so
    each compartment is projected separately.
    """
    c = config
    names = list(c.protein_channels)
    if channels.ndim != 4 or channels.shape[0] != len(names) or mask.shape != channels.shape[1:]:
        raise ValueError("Expected CZYX channels and matching ZYX mask")
    _, nz, h, w = channels.shape
    order = 1 if c.interpolation == "linear" else 0
    windows = {method: _windows(nz, c, method) for method in c.methods}

    # Band geometry for each side.
    sides = []
    for roi in rois:
        grid, distance = band_grid(roi, c.width_override_px)
        inside = (
            (grid[..., 0] >= 0)
            & (grid[..., 0] <= w - 1)
            & (grid[..., 1] >= 0)
            & (grid[..., 1] <= h - 1)
        )
        if c.boundary_mode == "error" and not inside.all():
            raise ValueError(f"{roi.roi_id}: band extends outside image")
        # valid marks which band points count (shape: band width x length).
        valid = inside.astype(float) if c.boundary_mode == "exclude" else np.ones(inside.shape)
        sides.append((roi, distance, inside, valid, *_band_sampler(grid, h, w, order)))

    _log(
        f"Measure: {len(rois)} side(s) x methods {list(c.methods)}, "
        f"each x {len(COMPARTMENTS)} compartments x {len(names)} channels"
    )
    records = []
    bar = _bar(desc="Measuring bands", total=len(COMPARTMENTS) * len(rois), unit="side")
    for comp in COMPARTMENTS:
        # Mask the whole stack once for all sides and methods.
        m = mask if comp == "nuclear" else ~mask if comp == "cytoplasmic" else np.ones_like(mask)
        values = channels
        if comp != "all" and c.mask_cap is not None:
            values = np.minimum(channels, c.mask_cap)
        masked = values * m  # (channel, z, y, x)
        # Each frame is projected once, and every side then samples the projection.
        projections = {
            (method, z0): (masked[:, z0:z1].max(axis=1), m[z0:z1].max(axis=0))
            for method in c.methods
            if method in ("block_mip", "full_mip")
            for z0, z1 in windows[method]
        }
        for roi, distance, inside, valid, idx, wts in sides:

            def band_sum(images):
                """Sample (..., y, x) images on the band, sum across its width -> (..., length)."""
                sampled = (images.reshape(*images.shape[:-2], -1)[..., idx] * wts).sum(-1)
                return (sampled.reshape(*sampled.shape[:-1], *valid.shape) * valid).sum(-2)

            # Sums for each plane, shaped (channel, z, length) and (z, length). The
            # per_slice and pooled methods both use them.
            slice_num, slice_weight = band_sum(masked), band_sum(m)
            one_band = valid.sum(axis=0)
            for method in c.methods:
                for z0, z1 in windows[method]:
                    if (method, z0) in projections:
                        image, mask_image = projections[method, z0]
                        numerator, weight, denom = band_sum(image), band_sum(mask_image), one_band
                    else:
                        # A pooled band covers one band area per plane.
                        numerator = slice_num[:, z0:z1].sum(axis=1)
                        weight = slice_weight[z0:z1].sum(axis=0)
                        denom = one_band * (z1 - z0 if method == "pooled" else 1)
                    fraction = np.divide(weight, denom, out=np.zeros_like(weight), where=denom > 0)
                    for ci, name in enumerate(names):
                        mean = np.divide(
                            numerator[ci],
                            weight,
                            out=np.full_like(weight, np.nan),
                            where=(weight > 0) & (fraction >= c.min_mask_fraction),
                        )
                        band = np.divide(
                            numerator[ci], denom, out=np.full_like(weight, np.nan), where=denom > 0
                        )
                        records.append(
                            pd.DataFrame(
                                dict(
                                    side=roi.roi_id,
                                    channel=name,
                                    compartment=comp,
                                    method=method,
                                    z_start=z0 + 1,
                                    z_stop=z1,
                                    position_index=np.arange(len(distance)),
                                    distance_um=np.arange(len(distance)) * cal["pixel_size_um"],
                                    geometry_distance_um=distance * cal["pixel_size_um"],
                                    masked_sum=numerator[ci],
                                    mask_weight=weight,
                                    band_weight=denom,
                                    mask_fraction=fraction,
                                    masked_mean=mean,
                                    band_mean=band,
                                    outside_fraction=1 - inside.mean(axis=0),
                                    width_px=int(round(c.width_override_px or roi.width_px)),
                                    saved_width_px=roi.width_px,
                                )
                            )
                        )
            bar.update()
    bar.close()
    return pd.concat(records, ignore_index=True)


# Settings that change a measurement. They are saved with it in measurement.json.
MEASURE_SETTINGS = (
    "protein_channels", "protein_median_size", "protein_outliers", "outlier_radius_px",
    "outlier_threshold", "mask_cap", "width_override_px", "interpolation", "boundary_mode",
    "methods", "z_start", "z_stop", "mip_block_size", "mip_start", "incomplete_block",
    "min_mask_fraction",
)  # fmt: skip


def run_sample(sample, config):
    """Measure one sample and save its table to measurements/<replicate>/<sample>/."""
    c = config.validate()
    s = dict(sample)
    start = time.perf_counter()
    _safe_output(s, c)
    sources = _roi_sources(s, c)
    if not sources:
        raise ValueError(f"{s['sample_id']} has no ROIs. Convert them with convert_rois or draw them in Fiji")
    rois = load_rois(sources)
    _log(
        f"=== {s['replicate']}/{s['sample_id']}: {len(rois)} side(s) "
        f"[{', '.join(r.roi_id for r in rois)}] from {Path(sources[0]).name} ==="
    )
    folder = sample_dir(c, s, "measurements")
    table_path = folder / "profiles.parquet"
    mask, cal, mask_meta = segment_sample(s, c)
    _log(f"Measure: mask ready after {_since(start)}")
    folder.mkdir(parents=True, exist_ok=True)
    provenance = dict(
        raw=file_identity(s["raw_path"]),
        raw_position=int(s.get("raw_position", 1)),
        rois=[file_identity(p) for p in sources],
        calibration=cal,
        settings={k: getattr(c, k) for k in MEASURE_SETTINGS},
        versions=version_info(),
    )
    _log("Measure: reading the protein channels from the original image")
    raw, _ = read_image(s["raw_path"], int(s.get("raw_position", 1)), c)
    if max(c.protein_channels.values()) > raw.shape[1]:
        raise ValueError("protein_channels names a channel the image does not have")
    step = time.perf_counter()
    filtered = np.stack(
        [
            _protein_filter(raw[:, index - 1], c, name)
            for name, index in c.protein_channels.items()
        ]
    )
    _log(f"Measure: protein filtering done in {_since(step)}")
    step = time.perf_counter()
    t = measure_arrays(filtered, mask, rois, cal, c)
    _log(f"Measure: {len(t):,} rows measured in {_since(step)}; saving {table_path}")
    for key in ("replicate", "imaging_day", "sample_id", "embryo_id", "phase", "phase_source"):
        t[key] = s.get(key, "")
    t.to_parquet(table_path, index=False, compression="zstd")
    if c.export_csv:
        _log("Measure: also writing profiles.csv (export_csv=True)")
        t.to_csv(folder / "profiles.csv", index=False)
    write_json(
        folder / "measurement.json",
        dict(provenance=provenance, mask=mask_meta, seconds=time.perf_counter() - start),
    )
    _log(f"=== {s['replicate']}/{s['sample_id']} done in {_since(start)} ===")
    return t


def run(samples, config):
    """Make masks for the selected samples, then measure them.

    Masks are made first, one sample at a time in this process, so the GPU holds
    one copy of the model. The samples are then measured by image_workers
    processes, which load the saved masks and do not use the GPU.
    """
    c = config.validate()
    records = samples.to_dict("records")
    if not records:
        raise ValueError("Select at least one sample")
    start = time.perf_counter()
    if c.mask_source == "compute":
        _log(f"Run: stage 1, masks for {len(records)} sample(s) on {choose_device(c.device)}")
        for s in _bar(records, desc="Masks", unit="sample"):
            segment_sample(s, c)
    measure_config = replace(c, reuse_cache=True)
    _log(f"Run: stage 2, measuring with image_workers={c.image_workers}")
    if c.image_workers == 1:
        tables = [
            run_sample(s, measure_config) for s in _bar(records, desc="Samples", unit="sample")
        ]
    else:
        import multiprocessing as mp

        _log(
            "Run: each worker process prints its steps to the kernel or terminal log, "
            "not to this cell. The bar below advances as samples finish, in order."
        )
        with ProcessPoolExecutor(c.image_workers, mp_context=mp.get_context("spawn")) as pool:
            tables = list(
                _bar(
                    pool.map(run_sample, records, [measure_config] * len(records)),
                    desc="Samples",
                    total=len(records),
                    unit="sample",
                )
            )
    _log(f"Run: all {len(records)} sample(s) finished in {_since(start)}; writing run metadata")
    result = _categorize(pd.concat(tables, ignore_index=True))
    folder = Path(c.output_root) / "runs" / time.strftime("%Y%m%d_%H%M%S")
    folder.mkdir(parents=True, exist_ok=True)
    write_json(folder / "config.json", asdict(c))
    write_json(folder / "samples.json", records)
    write_json(folder / "versions.json", version_info())
    samples.to_csv(folder / "samples.csv", index=False)
    result.attrs["run_directory"] = str(folder)
    return result


def load_measurements(samples, config):
    """Load the saved measurement tables of the given samples.

    This is enough when only analysis or plot settings changed. After changing
    masks, ROIs, filtering or Z settings, measure again with run().
    """
    tables = []
    for s in _bar(samples.to_dict("records"), desc="Loading saved measurements", unit="sample"):
        t = pd.read_parquet(sample_dir(config, s, "measurements") / "profiles.parquet")
        t["phase"] = str(s.get("phase", ""))
        t["imaging_day"] = s.get("imaging_day", s["replicate"])
        t["phase_source"] = s.get("phase_source", "")
        tables.append(t)
    result = _categorize(pd.concat(tables, ignore_index=True))
    _log(f"Loaded {len(result):,} saved measurement rows for {len(tables)} sample(s)")
    return result


LABELS = ("replicate", "imaging_day", "sample_id", "embryo_id", "side", "channel", "compartment", "method",
          "phase", "phase_source")


def _categorize(t):
    """Store the label columns as categories.

    Millions of rows share a few dozen labels, so storing each label once and
    small integer codes per row makes filtering and grouping much faster.
    """
    for name in LABELS:
        if name in t:
            t[name] = t[name].astype("category")
    return t


def apply_exclusions(table, config):
    """Mark the rows that fall in an exclusions_csv range (in a copy of the table)."""
    t = table.copy()
    t["excluded"] = False
    t["exclusion_reason"] = ""
    if not config.exclusions_csv:
        return t
    annotations = pd.read_csv(config.exclusions_csv, keep_default_na=False)
    required = {"replicate", "sample_id", "side", "start_um", "stop_um", "channels", "reason"}
    if not required <= set(annotations):
        raise ValueError(f"Exclusions need {sorted(required)}")
    for _, a in annotations.iterrows():
        side = str(a.side)
        side = side if side.startswith("side") else "side" + side
        if float(a.start_um) >= float(a.stop_um):
            raise ValueError("Exclusion start must precede stop")
        chosen = (
            (t.replicate == a.replicate)
            & (t.sample_id == a.sample_id)
            & (t.side == side)
            & t.distance_um.between(float(a.start_um), float(a.stop_um))
        )
        if a.channels not in ("all", "*"):
            chosen &= t.channel.isin(str(a.channels).split(";"))
        t.loc[chosen, "excluded"] = True
        t.loc[chosen, "exclusion_reason"] = a.reason
    return t


def _rolling(profile, window_um, step_um, kind, order=1):
    """Rolling median or Savitzky–Golay smoothing of a profile with gaps.

    Gaps (NaN) are filled by linear interpolation before filtering and put back
    afterwards, so the smoothing does not fill in values where there was no data.
    """
    a = np.asarray(profile, float)
    known = np.isfinite(a)
    if not known.any():
        return np.full(len(a), np.nan)
    width = 2 * int(window_um / step_um / 2) + 1
    filled = np.interp(np.arange(len(a)), np.flatnonzero(known), a[known])
    if kind == "median":
        filtered = ndi.median_filter(filled, size=width, mode="nearest")
    else:
        if width <= order or known.sum() <= width:
            return np.full(len(a), np.nan)
        filtered = savgol_filter(filled, width, order)
    return np.where(known, filtered, np.nan)


def side_profiles(measurements, config):
    """Turn measurements into one profile per side (or per embryo).

    Excluded ranges are removed, then the Z-planes or blocks of each side are
    combined. A baseline shared by an embryo's sides is subtracted, and with
    day_normalize each imaging day is scaled to a common level. The baseline is
    found on rolling-median-smoothed profiles but subtracted from the unsmoothed
    ones. Returns the profiles and a table of the baselines and day factors.
    """
    c = config
    start = time.perf_counter()
    _log(
        f"Profiles: applying exclusions, then reducing Z ({c.analysis_method}, {c.z_average}) "
        f"for {len(measurements):,} measured rows"
    )
    # Select the method first so the exclusions work on fewer rows.
    view = apply_exclusions(measurements[measurements.method == c.analysis_method], c)
    if view.empty:
        raise ValueError(f"No profiles for {c.analysis_method}")
    if "imaging_day" not in view:
        view["imaging_day"] = view.replicate
    group = [
        "replicate",
        "imaging_day",
        "sample_id",
        "embryo_id",
        "side",
        "channel",
        "compartment",
        "position_index",
    ]
    for name in ("masked_sum", "mask_weight", "band_weight"):
        view[name] = view[name].where(~view.excluded, 0.0)
    view["masked_mean"] = view.masked_mean.where(~view.excluded)
    t = view.groupby(group, sort=False, as_index=False).agg(
        distance_um=("distance_um", "first"),
        masked_sum=("masked_sum", "sum"),
        mask_weight=("mask_weight", "sum"),
        band_weight=("band_weight", "sum"),
        z_mean=("masked_mean", "mean"),
        z_median=("masked_mean", "median"),
        phase=("phase", "first"),
        phase_source=("phase_source", "first"),
        excluded=("excluded", "all"),
    )
    # Back to plain text. This table is much smaller, and text joins like
    # embryo_id + "/" + side need it.
    t = t.astype({k: str for k in LABELS if k in t and k != "phase"})
    t["phase"] = t.phase.astype(object).fillna("").astype(str)  # a missing phase becomes ""
    t["mask_fraction"] = t.mask_weight.div(t.band_weight.where(t.band_weight > 0)).fillna(0)
    t["raw_value"] = (
        t.masked_sum.div(t.mask_weight.where(t.mask_weight > 0)).where(
            t.mask_fraction >= c.min_mask_fraction
        )
        if c.z_average == "pool"
        else t["z_mean" if c.z_average == "mean" else "z_median"]
    )
    t["posterior_distance_um"] = t.distance_um
    if c.align_profiles == 1:
        # Measure each side from its anterior end. Every position of a band is kept
        # (excluded ones as gaps), so the last position is the band's anterior end.
        # This is done before sides are combined so that they line up at that end.
        band = ["replicate", "sample_id", "embryo_id", "side", "channel", "compartment"]
        g = t.groupby(band, sort=False)
        t["position_index"] = g.position_index.transform("max") - t.position_index
        t["distance_um"] = g.distance_um.transform("max") - t.distance_um
    if c.combine_sides:
        for name in ("masked_sum", "mask_weight"):
            t[name] = t[name].where(np.isfinite(t.raw_value), 0.0)
        t = t.groupby([k for k in group if k != "side"], sort=False, as_index=False).agg(
            distance_um=("distance_um", "first"),
            posterior_distance_um=("posterior_distance_um", "mean"),
            raw_value=("raw_value", "mean"),
            masked_sum=("masked_sum", "sum"),
            mask_weight=("mask_weight", "sum"),
            mask_fraction=("mask_fraction", "mean"),
            phase=("phase", "first"),
            phase_source=("phase_source", "first"),
            excluded=("excluded", "all"),
        )
        t["side"] = "combined"
        if c.side_average == "pool":
            t["raw_value"] = t.masked_sum.div(t.mask_weight.where(t.mask_weight > 0))
    normalization = []
    keys = ["replicate", "imaging_day", "embryo_id", "channel", "compartment"]
    groups = t.groupby(keys, sort=False)
    _log(
        f"Profiles: baseline for {groups.ngroups} embryo x channel x compartment group(s)"
    )
    for values, g in _bar(groups, desc="Baselines", total=groups.ngroups, unit="group"):
        medians = []
        for _, side in g.groupby("side", sort=False):
            side = side.sort_values("position_index")
            d = side.distance_um.to_numpy()
            step = float(np.median(np.diff(d))) if len(d) > 1 else 1.0
            m = _rolling(side.raw_value.to_numpy(), c.normalization_window_um, step, "median")
            medians.extend(m[np.isfinite(m)])
        background = (
            float(np.percentile(medians, c.background_percentile))
            if c.background_subtract and medians
            else 0.0
        )
        status = "ok" if medians or not c.background_subtract else "no baseline: profile empty"
        t.loc[g.index, "baseline"] = background
        t.loc[g.index, "value"] = g.raw_value - background
        t.loc[g.index, "normalization_status"] = status
        normalization.append(
            dict(
                zip(keys, values),
                baseline=background,
                status=status,
            )
        )
    normalization = pd.DataFrame(normalization)
    if c.day_normalize:
        # Each imaging day, channel and compartment gets one factor, the median of
        # the background-subtracted values within day_reference_um over all of
        # that day's embryos. Dividing by it puts the days on one scale, while
        # embryos imaged on the same day keep their relative brightness. Using a
        # fixed distance range keeps the factor independent of how far each ROI
        # was traced.
        day = ["imaging_day", "channel", "compartment"]
        lo, hi = c.day_reference_um
        factor = t[t.posterior_distance_um.between(lo, hi)].groupby(day).value.median()
        factor = factor.where(factor > 0).rename("day_factor")  # a factor of 0 or less cannot be used for scaling
        t = t.join(factor, on=day)
        t["value"] = t.value / t.day_factor
        t.loc[t.day_factor.isna(), "normalization_status"] = "no day factor: reference median <= 0"
        normalization = normalization.join(factor, on=day)
        _log(f"Profiles: day factors {factor.round(2).to_dict()}")
    t["profile_unit"] = t.embryo_id if c.combine_sides else t.embryo_id + "/" + t.side
    _log(f"Profiles: {t.profile_unit.nunique()} profile unit(s) ready in {_since(start)}")
    return t, normalization


def summarize(profiles, config, by_phase=True):
    """Average the profiles at each distance, per phase or over all of them.

    Profiles of different lengths are averaged where they have data, so the
    count n falls toward the anterior end. n counts sides, or embryos when
    combine_sides is True. n_embryos and n_replicates are reported alongside.
    All profiles must share one pixel size, since they are averaged position by
    position.
    """
    c = config
    _log(f"Summary: {'per phase' if by_phase else 'overall'}")
    t = profiles.copy()
    if by_phase:
        t = t[t.phase.isin(c.phase_labels)]
    group = (["phase"] if by_phase else []) + ["channel", "compartment", "position_index"]
    average = c.phase_average if by_phase else c.average
    error = c.phase_error if by_phase else c.error
    minimum = c.phase_minimum_lines if by_phase else c.minimum_lines
    columns = group + [
        "distance_um",
        "center",
        "low",
        "high",
        "n",
        "n_embryos",
        "n_replicates",
        "unit",
    ]
    if t.empty:
        return pd.DataFrame(columns=columns)
    t["valid_embryo"] = t.embryo_id.where(np.isfinite(t.value))
    t["valid_replicate"] = t.replicate.where(np.isfinite(t.value))
    g = t.groupby(group, sort=False)
    stats = g.agg(
        distance_um=("distance_um", "first"),
        dmin=("distance_um", "min"),
        dmax=("distance_um", "max"),
        center=("value", average),
        sd=("value", "std"),
        n=("value", "count"),
        n_embryos=("valid_embryo", "nunique"),
        n_replicates=("valid_replicate", "nunique"),
    )
    if ((stats.dmax - stats.dmin) > np.maximum(1e-4, abs(stats.distance_um) * 1e-4)).any():
        raise ValueError(
            "Profiles have different pixel sizes, so their positions do not line up "
            "for averaging."
        )
    if error == "iqr":
        quantiles = g.value.quantile([0.25, 0.75]).unstack()
        stats["low"] = quantiles[0.25]
        stats["high"] = quantiles[0.75]
    else:
        spread = stats.sd / np.sqrt(stats.n) if error == "sem" else stats.sd
        stats["low"] = stats.center - spread
        stats["high"] = stats.center + spread
    stats.loc[stats.n < minimum, ["center", "low", "high"]] = np.nan
    stats["unit"] = "embryo" if c.combine_sides else "side"
    return stats.reset_index()[columns]


# Colours and fonts of the phase figures.
SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
PHASE_FIGURE_STYLE = {
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": MUTED,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 18,
}


def _plotted_center(g, c):
    """Distances and average of one summary curve, smoothed when smooth_average is True."""
    d, center = g.distance_um.to_numpy(), g.center.to_numpy()
    if c.smooth_average and len(d) > 1:
        step = float(np.median(np.diff(d)))
        center = _rolling(center, c.savgol_window_um, step, "savgol", c.savgol_order)
    return d, center


def _phase_colors(c):
    """{phase label: RGBA}, evenly spaced along phase_colormap_range."""
    import matplotlib

    cmap = matplotlib.colormaps[c.phase_colormap]
    positions = np.linspace(*c.phase_colormap_range, len(c.phase_labels))
    if len(c.phase_labels) == 1:
        positions = [np.mean(c.phase_colormap_range)]
    return {str(phase): cmap(p) for phase, p in zip(c.phase_labels, positions)}


def _phase_y_label(c):
    """Y-axis label that says whether profiles are background subtracted and normalised."""
    steps = [s for s, on in (("background subtracted", c.background_subtract),
                             ("normalised", c.day_normalize)) if on]
    return f"Intensity ({', '.join(steps) or 'raw'}) (AU)"


def _distance_label(c):
    """X-axis label naming the end that profiles are aligned to."""
    if c.align_profiles == 1:
        return "Anterior to posterior distance (µm)"
    return "Posterior to anterior distance (µm)"


def _method_label(c):
    """Readable name of the Z method used for the profiles."""
    return {
        "block_mip": f"{c.mip_block_size}-slice MIP",
        "per_slice": "Per slice",
        "pooled": "Pooled",
        "full_mip": "Full MIP",
    }.get(c.analysis_method, c.analysis_method)


def _phase_title(c, channel, compartment):
    """Figure title naming the channel, compartment and how profiles were averaged."""
    unit = "embryo" if c.combine_sides else "side"
    return (
        f"{channel}, {compartment} ({_method_label(c)}): {c.z_average} per {unit}, "
        f"{c.phase_average} ± {c.phase_error} across {unit}s"
    )


def _phase_counts(profiles, c, channel, compartment, phase=None):
    """Text such as 'n = 8 sides from 4 embryos with N = 2', counting profiles with data.

    With phase=None, profiles of every phase are counted.
    """
    keep = (
        (profiles.channel == channel)
        & (profiles.compartment == compartment)
        & np.isfinite(profiles.value)
    )
    if phase is not None:
        keep &= profiles.phase == phase
    selected = profiles[keep]
    n, embryos = selected.profile_unit.nunique(), selected.embryo_id.nunique()
    replicates = selected.replicate.nunique()
    if c.combine_sides:
        return f"n = {n} embryo{'s' * (n != 1)} with N = {replicates}"
    return f"n = {n} side{'s' * (n != 1)} from {embryos} embryo{'s' * (embryos != 1)} with N = {replicates}"


def plot_phases(profiles, config, summary=None):
    """One figure per channel with each phase's average profile and spread."""
    import matplotlib.pyplot as plt

    c = config
    stats = summarize(profiles, c, True) if summary is None else summary
    colors = _phase_colors(c)
    figures = {}
    for channel, compartment in c.phase_compartments.items():
        with plt.rc_context(PHASE_FIGURE_STYLE):
            fig, ax = plt.subplots(figsize=(14, 7.6), constrained_layout=True, dpi=c.dpi)
            subset = (
                stats[(stats.channel == channel) & (stats.compartment == compartment)]
                if not stats.empty
                else stats
            )
            for phase, g in subset.groupby("phase", sort=False):
                g = g.sort_values("position_index")
                d, center = _plotted_center(g, c)
                color = colors.get(str(phase))
                label = f"Phase {phase} ({_phase_counts(profiles, c, channel, compartment, phase)})"
                ax.fill_between(d, g.low.to_numpy(), g.high.to_numpy(), alpha=0.2, color=color, lw=0)
                ax.plot(d, center, label=label, color=color, lw=3)
            ax.set_xlabel(_distance_label(c))
            ax.set_ylabel(_phase_y_label(c))
            if not subset.empty:
                ax.legend(frameon=False, labelcolor=MUTED, fontsize=14)
            fig.suptitle(_phase_title(c, channel, compartment), x=0.01, ha="left", color=INK, fontsize=22)
        figures[channel] = fig
    return figures


def plot_individual_profiles(profiles, config):
    """Grid of every profile in grey with the overall average and spread on top.

    There is one panel per channel and compartment. Phases are not needed, and
    the average, error and minimum_lines settings apply.
    """
    import matplotlib
    import matplotlib.pyplot as plt

    c = config
    stats = summarize(profiles, c, by_phase=False)
    color = matplotlib.colormaps[c.phase_colormap](c.phase_colormap_range[0])
    unit = "embryo" if c.combine_sides else "side"
    with plt.rc_context(PHASE_FIGURE_STYLE):
        fig, axes = plt.subplots(
            len(c.plot_channels),
            len(c.plot_compartments),
            figsize=(8 * len(c.plot_compartments), 4 * len(c.plot_channels)),
            squeeze=False,
            sharex=True,
            constrained_layout=True,
            dpi=c.dpi,
        )
        for row, channel in enumerate(c.plot_channels):
            for col, compartment in enumerate(c.plot_compartments):
                axis = axes[row, col]
                selected = profiles[
                    (profiles.channel == channel) & (profiles.compartment == compartment)
                ]
                n_units = max(selected.profile_unit.nunique(), 1)
                for _, g in selected.groupby("profile_unit", sort=False):
                    g = g.sort_values("position_index")
                    axis.plot(g.distance_um, g.value, color=MUTED, lw=0.3, alpha=min(0.5, 40 / n_units))
                g = (
                    stats[(stats.channel == channel) & (stats.compartment == compartment)]
                    if not stats.empty
                    else stats
                ).sort_values("position_index")
                if not g.empty:
                    d, center = _plotted_center(g, c)
                    axis.fill_between(d, g.low.to_numpy(), g.high.to_numpy(), color=color, alpha=0.3, lw=0)
                    axis.plot(d, center, color=color, lw=2)
                values = selected.value.to_numpy(dtype=float)
                if np.isfinite(values).any():
                    axis.set_ylim(*np.nanpercentile(values, [0.2, 99.8]))
                axis.set_title(
                    f"{channel}, {compartment}\n{_phase_counts(profiles, c, channel, compartment)}",
                    loc="left",
                    color=INK,
                )
                if row == len(c.plot_channels) - 1:
                    axis.set_xlabel(_distance_label(c))
                if col == 0:
                    axis.set_ylabel(_phase_y_label(c))
        fig.suptitle(
            f"{_method_label(c)}: all {unit} profiles (grey), {c.average} ± {c.error} across {unit}s; "
            f"average shown where n ≥ {c.minimum_lines}",
            x=0.01,
            ha="left",
            color=INK,
        )
    return fig


def plot_phase_panels(profiles, config, summary=None):
    """One figure per channel with a panel per phase, showing its profiles and average."""
    import matplotlib.pyplot as plt

    c = config
    stats = summarize(profiles, c, True) if summary is None else summary
    figures = {}
    for channel, compartment in c.phase_compartments.items():
        with plt.rc_context(PHASE_FIGURE_STYLE):
            fig, axes = plt.subplots(
                1,
                len(c.phase_labels),
                figsize=(8 * len(c.phase_labels), 7),
                sharex=True,
                sharey=True,
                squeeze=False,
                constrained_layout=True,
                dpi=c.dpi,
            )
            colors = _phase_colors(c)
            for axis, phase in zip(axes[0], c.phase_labels):
                selected = profiles[
                    (profiles.channel == channel)
                    & (profiles.compartment == compartment)
                    & (profiles.phase == phase)
                ]
                if c.show_individuals:
                    for _, g in selected.groupby("profile_unit"):
                        g = g.sort_values("position_index")
                        axis.plot(g.distance_um, g.value, color=MUTED, lw=0.5, alpha=0.5)
                if not stats.empty:
                    g = stats[
                        (stats.channel == channel)
                        & (stats.compartment == compartment)
                        & (stats.phase == phase)
                    ].sort_values("position_index")
                    d, center = _plotted_center(g, c)
                    axis.fill_between(
                        d,
                        g.low.to_numpy(),
                        g.high.to_numpy(),
                        color=colors.get(str(phase)),
                        alpha=0.2,
                        lw=0,
                    )
                    axis.plot(d, center, color=colors.get(str(phase)), lw=2)
                axis.set_title(
                    f"Phase {phase}\n{_phase_counts(profiles, c, channel, compartment, phase)}",
                    loc="left",
                    color=INK,
                )
                axis.set_xlabel(_distance_label(c))
            axes[0, 0].set_ylabel(_phase_y_label(c))
            fig.suptitle(_phase_title(c, channel, compartment), x=0.01, ha="left", color=INK)
        figures[channel] = fig
    return figures


def export_analysis(profiles, normalization, summary, config, folder=None):
    """Save the profiles, baselines, phase summary and settings as CSV and JSON files."""
    folder = Path(folder or Path(config.output_root) / "analysis" / time.strftime("%Y%m%d_%H%M%S"))
    folder.mkdir(parents=True, exist_ok=True)
    _log(f"Saving analysis tables to {folder}")
    profiles.to_csv(folder / "side_or_embryo_profiles.csv", index=False)
    normalization.to_csv(folder / "normalization.csv", index=False)
    summary.to_csv(folder / "phase_summary.csv", index=False)
    write_json(folder / "config.json", asdict(config))
    write_json(folder / "versions.json", version_info())
    return folder



def qc_figure(sample, config, z=None):
    """Raw DAPI at one Z-plane, beside the same plane with the mask and ROI bands drawn on it."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    s = dict(sample)
    _log(f"QC figure: {s['replicate']}/{s['sample_id']}")
    a, cal = read_image(s["raw_path"], int(s.get("raw_position", 1)), config)
    mask, _, _ = segment_sample(s, config)
    z = len(a) // 2 if z is None else z
    _log(f"QC figure: drawing Z={z + 1}")
    fig, ax = plt.subplots(1, 2, figsize=(18, 12), constrained_layout=True, dpi=400)
    dapi = a[z, config.dapi_channel - 1]
    lo, hi = np.percentile(dapi, config.qc_percentiles)
    for axis in ax:
        axis.imshow(dapi, cmap="gray", vmin=lo, vmax=hi)
    rgba = np.zeros((*dapi.shape, 4))
    rgba[mask[z]] = [0, 1, 0, 0.3]
    ax[1].imshow(rgba)
    for r in load_rois(_roi_sources(s, config)):
        grid, _ = band_grid(r, config.width_override_px)
        poly = np.vstack([grid[0], grid[-1, ::-1]])
        for axis in ax:
            axis.add_patch(Polygon(poly, fill=False, edgecolor="skyblue", lw=2))
            axis.plot(r.xy[0, 0], r.xy[0, 1], "ro", ms=5)
    ax[0].set_title(f"{s['replicate']}/{s['sample_id']}, Z={z+1}: DAPI (raw)")
    ax[1].set_title("Nuclear mask + ROI")
    return fig

