# IHC gradient analysis

This pipeline measures how Her7-Venus, pMEK and ppERK immunostaining changes
along the presomitic mesoderm (PSM) of zebrafish embryos and compares the
profiles between oscillation phases.

For each embryo image, the pipeline segments the nuclei in the DAPI channel and
uses that mask to separate nuclear from cytoplasmic signal. The signal of each
protein is then measured along hand-drawn bands that run from the posterior to
the anterior end of the PSM. The resulting profiles are background subtracted,
averaged within each phase and plotted.

The repository contains:

- `IHC_analysis.ipynb`, the notebook used to run the analysis.
- `ihc_pipeline.py`, which contains the image processing and analysis code.
- `ihc_roi_editor.py`, a napari window for drawing bands and correcting masks.
- `IHC_ROI_only_CSV.ijm`, a Fiji macro for drawing bands in the same file format.
- `PARAMETERS_IHC.md`, which describes every setting.
- `sync.bat` and `sync.sh`, which copy data and results between a computer and
  the lab server.

## Installation

The Python environment is managed with [uv](https://docs.astral.sh/uv/).
`pyproject.toml` lists the dependencies and `uv.lock` records their exact
versions, so every computer ends up with the same environment.

On a Mac with Apple Silicon, install the environment with:

```sh
uv sync --locked
```

Windows and Linux computers with an NVIDIA GPU need the CUDA build of PyTorch,
which is selected with the `cuda` extra. The second command checks that
PyTorch can see the GPU:

```sh
uv sync --locked --extra cuda
uv run --locked --extra cuda python -c "import torch; print(torch.cuda.is_available())"
```

On these computers, `--extra cuda` has to be included in every uv command,
because uv otherwise switches back to the CPU build of PyTorch. The NVIDIA
driver comes from the operating system's driver installation.

To set up a new computer, clone the repository and run `uv sync`. Each computer
builds its own `.venv`, since a virtual environment is tied to the machine it
was created on. The section on working across several computers explains how
the data and results are copied over.

By default, the notebook reads images from `Exp_1__mek_h7v_pperk`, the folder
next to `analysis_scripts`, and writes results to `analysis_scripts/ihc_results`.
Other locations can be set with `EXPERIMENTAL_DATA` and `OUTPUT_ROOT` in the
notebook, or with the environment variables `IHC_EXPERIMENTAL_DATA` and
`IHC_OUTPUT_ROOT`.

To run the notebook in VS Code, select the project's `.venv` as the kernel.
JupyterLab can be started with:

```sh
uv run --locked --group notebook jupyter lab IHC_analysis.ipynb
```

## Workflow

1. In the notebook's settings cell, check the folder that holds the replicate
   folders (`N_1...`, `N_2...`) and the output folder.
2. Run the inventory cell. It lists each sample with its NIS-denoised file, its
   ROIs and its phase. A sample is a single embryo image, which is either an
   exported TIFF or one position in an ND2 file. When a position has been
   exported as a TIFF, the pipeline uses the TIFF.
3. Convert any ImageJ ROI ZIPs to CSV, and draw bands in napari or Fiji for
   samples that have none.
4. Measure a single sample and check its QC figure, then measure the remaining
   samples.
5. Run the profile and phase cells.

Once the measurements exist, setting `RELOAD_SAVED = True` and
`RUN_MEASUREMENTS = False` loads them from disk. This is useful when only the
analysis or plot settings change. Changes to the masks, ROIs, protein filtering
or Z settings need a new measurement.

## Nuclear masks

The DAPI image used for segmentation depends on `use_ai_denoise`. When it is
True, the pipeline uses the sample's NIS-denoised export. Exports are matched to
raw files by name, so `triple_ihc_013 - Denoised.nd2` is paired with
`triple_ihc_013.nd2`. If a sample has both a denoised ND2 and a denoised TIFF,
`nis_preference` decides which one is used, and a sample with no denoised export
stops the run with an error. When `use_ai_denoise` is False, the pipeline
prepares the raw DAPI itself. It applies a median filter, subtracts the slowly
varying background and evens out the contrast with CLAHE.

Segmentation uses Cellpose 4 when `use_cellpose` is True. When it is False,
each Z-plane is thresholded with the IsoData method.

Computed masks are written to a folder whose name records how they were made:
`masks_cellpose_denoised/`, `masks_cellpose_raw/`, `masks_threshold_denoised/`
or `masks_threshold_raw/`. Masks made with different combinations of
`use_cellpose` and `use_ai_denoise` are therefore kept in separate folders. When
`reuse_cache` is True, a mask saved by an earlier run is loaded instead of being
computed again. The pipeline does not check whether the other Cellpose or
filtering settings, or the code itself, have changed since the mask was saved,
so after such a change `reuse_cache=False` makes new masks.

Masks can also come from other sources. `mask_source='existing'` loads a
`Mask.tif` from the folder named after the image, such as a mask made in Fiji,
and `mask_source='edited'` loads a mask that was corrected in the napari editor.

The cytoplasmic mask covers all pixels outside the nuclei. Masks are saved as
16-bit ImageJ TIFFs, with 0 outside and 4095 inside, and include the pixel size
so that Fiji opens them at the correct scale.

## Bands (ROIs)

The bands of each sample are stored in one CSV file:

```csv
roi_id,vertex_index,x_px,y_px,width_px
side0,0,1967,390,100
side0,1,1682,337,100
```

Each row is a point that was clicked when the band was drawn, given in image
pixels counted from 0. The first point of each side is at the posterior end.
The centre line of the band is a smooth spline through these points, the same
spline ImageJ uses for polylines, and `width_px` sets the width of the band.
Fiji and napari both read and write this format and draw the band from the same
spline, so a band looks the same in either program. ImageJ ROI ZIPs that are
converted to this format are saved in `rois/` in the output folder, and the ZIP
files themselves are kept unchanged.

Before saving, the napari editor checks whether the ROI file has changed on disk
since it was opened. This stops one editor window from overwriting changes that
were saved from another.

## Measurements

At each position along a band, the pipeline samples the protein image and the
mask at a row of points across the band's width and sums them:

```text
masked_sum  = sum of (protein × mask)
mask_weight = sum of mask
masked_mean = masked_sum / mask_weight
band_mean   = masked_sum / band_weight
```

`masked_mean` is the mean intensity within the compartment, and `band_mean` is
the compartment's signal divided by the area of the band. Where a compartment
covers less than `min_mask_fraction` of the band, its `masked_mean` is left
empty.

Before the measurement, each protein channel goes through three steps, each of
which can be switched off. A median filter reduces noise, pixels that are much
brighter than their surroundings are replaced by the local median, and values
inside the masks are capped at `mask_cap`, which is 4095, the maximum of the
12-bit camera. The protein signal is read from the original image, and the
NIS-denoised export is used only for the DAPI.

Each run measures the Z-stack in four ways:

- `per_slice` measures every Z-plane separately.
- `pooled` adds up the sums from all planes before dividing.
- `block_mip` measures maximum projections of blocks of planes, by default
  planes 2–4, 5–7 and so on up to 29–31.
- `full_mip` measures a maximum projection of the whole stack.

The projections are made after masking, so nuclear and cytoplasmic signal are
projected separately. `analysis_method` selects which of the four the profiles
are built from.

## Profiles and phases

For each side, the planes or blocks are first combined according to
`z_average`, and then a baseline is subtracted. To find the baseline, the
profiles are smoothed with a rolling median, and the 5th percentile across all
sides of an embryo is used as that embryo's baseline. With `day_normalize=True`,
the profiles from each imaging day are also divided by that day's median over
the first 400 µm. The imaging day comes from the ND2 acquisition date, so a
replicate imaged over several days gets a separate factor for each day.

Phases are read from the `phases.csv` file of each replicate. The file has
either the columns `file,phase`, with entries such as
`triple_ihc_013.nd2 P2,3`, or the columns `sample_id,phase`. `phase_files` can
point a replicate to a different file. Samples that have no phase are included
in the individual profile plots.

The phase plots show the mean profile of each phase and its spread. By default
every side counts as a separate profile, so `n` is the number of sides, and the
plots also report how many embryos and replicates these come from. With
`combine_sides=True`, the sides of an embryo are merged into a single profile.
Each profile ends where its band ends, so `n` decreases toward the anterior end.
The plotted mean can be smoothed for display, while the saved statistics are
calculated from the unsmoothed profiles.

Regions to leave out of the analysis, such as a bright clump, are listed in an
exclusions CSV (`IHC_exclusions.csv`) with the columns
`replicate,sample_id,side,start_um,stop_um,channels,reason`.

## Output folder

- `masks_<method>_<source>/` holds the computed masks, the DAPI image that was
  segmented, the Cellpose labels and `mask.json`, which records the settings
  and input files.
- `edited_masks/` holds masks corrected in napari.
- `rois/` holds the band CSVs.
- `measurements/` holds a `profiles.parquet` and a `measurement.json` for each
  sample.
- `runs/` records the settings and sample list of each run.
- `analysis/` holds the profile and summary tables, normalization factors and
  figures.

`mask.json`, `measurement.json` and each run's `versions.json` also record the
git commit the code was at, under `git_commit`, and the tracked files that had
uncommitted changes, under `git_uncommitted`. An empty list means the code
matched that commit exactly.

## Working on several computers

The raw data and results are stored on the lab server and the code is stored on
GitHub. Each computer keeps a local copy of the project in the same layout as
the server:

```text
MEK_oscillation_project/
  Exp_1__mek_h7v_pperk/     raw data, copied from the server
  analysis_scripts/         git clone, with its own .venv
    ihc_results/            results
```

The scripts `sync.bat` (Windows) and `sync.sh` (Mac) copy files between the
local copy and the server. They are run from `analysis_scripts`:

```sh
sync pull        # on a Mac: ./sync.sh pull
sync push        # on a Mac: ./sync.sh push
```

`pull` updates the code with `git pull` and then copies new raw data and
results from the server. `push` copies the results from this computer to the
server and lists any code changes that have not been committed. Adding
`--dry-run` lists what would be copied without copying anything.

The raw data is mirrored, which means the local copy is made to match the
server, and local raw files that are not on the server are deleted. Results are
copied in both directions without deleting anything, and a file is replaced
only when the incoming copy is newer. The Cellpose models, the
matplotlib cache and the napari settings belong to each computer and are not
copied.

A session usually starts with `sync pull` and ends with `sync push`. ROIs and
masks should be edited on one computer at a time. If the same file is edited on
two computers, the newer copy replaces the older one when they are synced, and
the changes in the older copy are lost.

## Committing the notebook

The notebook's outputs contain figures and local file paths, so the repository
stores the notebook with its outputs removed. Git removes them during each
commit, using a filter that is named in `.gitattributes`. The filter has to be
set up once in each clone:

```sh
git config filter.strip-notebook.clean "python -c \"import json,sys; nb=json.load(sys.stdin); [c.update(outputs=[], execution_count=None) for c in nb['cells'] if c['cell_type']=='code']; nb['metadata'].pop('widgets', None); sys.stdout.write(json.dumps(nb, indent=1, ensure_ascii=False) + '\\n')\""
git config filter.strip-notebook.smudge cat
git config filter.strip-notebook.required true
```

The notebook on your computer keeps its outputs, and only the copy stored in
the repository is cleared.

## Tests

The tests use small synthetic images and run on computers without a GPU:

```sh
uv run --locked -m pytest -q -p no:napari test_ihc_pipeline.py
```

## References

- [Cellpose segmentation settings](https://cellpose.readthedocs.io/en/latest/settings.html)
- [ImageJ macro functions](https://imagej.net/ij/developer/macro/functions.html)
