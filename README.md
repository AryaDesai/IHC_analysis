# IHC gradient analysis

This pipeline measures how Her7-Venus, pMEK and ppERK immunostaining changes
along the presomitic mesoderm (PSM) of zebrafish embryos, and compares the
profiles between oscillation phases.

For each embryo image it makes a nuclear mask from the DAPI channel. It then
measures each protein inside and outside the nuclei along hand-drawn bands
that run from posterior to anterior. The profiles are background subtracted,
averaged and plotted per phase.

The files are:

- `IHC_analysis.ipynb`: the notebook you work in.
- `ihc_pipeline.py`: all of the calculations.
- `ihc_roi_editor.py`: a napari window for drawing bands and correcting masks.
- `IHC_ROI_only_CSV.ijm`: a Fiji macro that draws bands in the same file format.
- `PARAMETERS_IHC.md`: a description of every setting.

## Installation

The project uses [uv](https://docs.astral.sh/uv/). `pyproject.toml` lists the
dependencies, and `uv.lock` pins their exact versions so every machine gets the
same environment.

On a Mac with Apple Silicon:

```sh
uv sync --locked
```

On a Windows or Linux machine with an NVIDIA GPU, add the CUDA build of PyTorch:

```sh
uv sync --locked --extra cuda
uv run --locked --extra cuda python -c "import torch; print(torch.cuda.is_available())"
```

On that machine, keep `--extra cuda` in every uv command. Otherwise uv
reinstalls the CPU-only PyTorch. The NVIDIA driver itself is installed with
the operating system, not by uv.

To move the project to another computer, copy the project files and run
`uv sync` there. Do not copy `.venv`, because it only works on the machine
that created it. The images can live anywhere; point `EXPERIMENTAL_DATA` and
`OUTPUT_ROOT` in the notebook at them.

In VS Code, select the project's `.venv` as the notebook kernel. To use
JupyterLab instead:

```sh
uv run --locked --group notebook jupyter lab IHC_analysis.ipynb
```

## Workflow

1. In the notebook's settings cell, set the folder that holds the replicate
   folders (`N_1...`, `N_2...`) and a separate output folder.
2. Run the inventory cell. It lists every sample it found, the NIS-denoised file
   paired with it, its ROIs and its phase. A sample is one embryo image: an
   exported TIFF, or one position of an ND2 file. When a position exists as both,
   the TIFF is used.
3. Convert existing ImageJ ROI ZIPs to CSV, or draw missing bands in napari or Fiji.
4. Measure. Start with one sample and look at the QC figure before running all of them.
5. Run the profile and phase cells.

After the first run, set `RELOAD_SAVED = True` and `RUN_MEASUREMENTS = False` to
change analysis or plot settings without measuring again. Changing masks, ROIs,
protein filtering or Z settings requires measuring again.

## Nuclear masks

The DAPI that is segmented comes from one of two places:

- With `use_ai_denoise=True`, the sample's NIS-denoised export is used as it is.
  The pairing is by filename: `triple_ihc_013 - Denoised.nd2` belongs to
  `triple_ihc_013.nd2`. If both a denoised ND2 and a denoised TIFF exist,
  `nis_preference` picks one. A sample without an NIS-denoised export stops the
  run with an error.
- With `use_ai_denoise=False`, the raw DAPI is median filtered, its smooth
  background is subtracted, and CLAHE evens out the contrast.

`use_cellpose=True` segments the DAPI with Cellpose 4. With `False`, each
Z-plane is thresholded with the IsoData method instead.

Masks are saved under `masks/` in the output folder and reused on later runs
while `reuse_cache` is True. The saved mask is not compared with the current
settings, so set `reuse_cache=False` after changing a mask setting.

`mask_source='existing'` reads a `Mask.tif` from the folder named after the image
(for example one made in Fiji), and `mask_source='edited'` reads a mask corrected
in the napari editor.

The cytoplasmic mask is everything that is not nucleus. Masks are saved as 16-bit
ImageJ TIFFs (0 outside, 4095 inside) with their pixel size, so Fiji opens them
directly.

## Bands (ROIs)

Each sample's bands are stored in one CSV file:

```csv
roi_id,vertex_index,x_px,y_px,width_px
side0,0,1967,390,100
side0,1,1682,337,100
```

Each row is one clicked point, in image pixels counted from 0. The first point of
a side is its posterior end. A smooth spline through the points (the one ImageJ
uses for polylines) is the centre of the band, and `width_px` is its width.
Fiji and napari read and write this format, and both draw the band from the same
spline. The converted files are kept in `rois/` in the output folder; the
original ZIPs are left as they are.

Before saving, the napari editor checks that the file has not changed on disk
since it was loaded, so two open editors cannot overwrite each other's work.

## Measurements

At each position along a band, the pipeline samples the protein image and the
mask at the same points across the band's width and adds them up:

```text
masked_sum  = sum of (protein × mask)
mask_weight = sum of mask
masked_mean = masked_sum / mask_weight    mean intensity inside the compartment
band_mean   = masked_sum / band_weight    signal of the compartment per band area
```

When a compartment covers less than `min_mask_fraction` of the band at a
position, its `masked_mean` is left blank there.

Before measuring, each protein channel is median filtered, pixels much brighter
than their surroundings are replaced by the local median, and values inside the masks are capped at `mask_cap`
(4095, the 12-bit camera maximum). Each step can be turned off. Proteins are
always measured from the original image, not from an NIS export.

Z is handled in four ways, all measured in the same run:

- `per_slice` measures each Z-plane.
- `pooled` adds the sums of all planes before dividing.
- `block_mip` takes maximum projections of blocks of planes (by default 2–4,
  5–7, ..., 29–31).
- `full_mip` projects the whole stack.

Projections are taken after masking, so nuclear and cytoplasmic signal are
projected separately. `analysis_method` chooses which of these the profiles
use.

## Profiles and phases

The planes or blocks of each side are combined (`z_average`). A baseline is
then subtracted: the profiles are smoothed with a rolling median, and the 5th
percentile over all of an embryo's sides becomes that embryo's baseline. With
`day_normalize=True`, each imaging day is also divided by its median value over
the first 400 µm. The imaging day is read from the ND2 acquisition date, so
one replicate imaged on several days gets a factor for each day.

Phases come from each replicate's `phases.csv`. This file has either the
columns `file,phase` (with entries such as `triple_ihc_013.nd2 P2,3`) or
`sample_id,phase`. A different file can be chosen per replicate with
`phase_files`. Samples without a phase still appear in the individual profiles.

The phase plots show the average of each phase with its spread. By default each
side is its own profile, so `n` counts sides; the plots also give the number of
embryos and replicates. With `combine_sides=True`, an embryo's sides become one
profile. Profiles end where their band ends and are not stretched to a common
length, so `n` falls toward the anterior end. The plotted average can be
smoothed; the saved statistics are not.

Ranges that should be left out, such as a bright clump, go in an exclusions CSV
(`IHC_exclusions.csv`) with the columns
`replicate,sample_id,side,start_um,stop_um,channels,reason`.

## Output folder

- `masks/`: computed masks, the DAPI that was segmented, Cellpose labels, and
  `mask.json` with the settings and input files used.
- `edited_masks/`: masks corrected in napari.
- `rois/`: band CSVs.
- `measurements/`: one `profiles.parquet` per sample, with `measurement.json`.
- `runs/`: the settings and sample list of each run.
- `analysis/`: profile and summary CSVs, normalization factors and figures.

## Committing the notebook

The repository stores the notebook without its outputs, which hold figures and
local file paths. Git removes them when a commit is made, through a filter named
in `.gitattributes`. After cloning, set the filter up once in the clone:

```sh
git config filter.strip-notebook.clean "python -c \"import json,sys; nb=json.load(sys.stdin); [c.update(outputs=[], execution_count=None) for c in nb['cells'] if c['cell_type']=='code']; nb['metadata'].pop('widgets', None); sys.stdout.write(json.dumps(nb, indent=1, ensure_ascii=False) + '\\n')\""
git config filter.strip-notebook.smudge cat
git config filter.strip-notebook.required true
```

The outputs stay in your local notebook; only the committed copy is cleared.

## Tests

```sh
uv run --locked -m pytest -q -p no:napari test_ihc_pipeline.py
```

The tests use small synthetic images and do not need a GPU.

## References

- [Cellpose segmentation settings](https://cellpose.readthedocs.io/en/latest/settings.html)
- [ImageJ macro functions](https://imagej.net/ij/developer/macro/functions.html)
