# Settings

Every setting of the pipeline, grouped as in the notebook's parameter table. This file is
generated from `Config` in `ihc_pipeline.py`. In the notebook, `ihc.parameter_table(cfg)`
shows the values currently in use.

## Files

| Setting | Default | Description |
|---|---|---|
| `experimental_data` | `''` | Folder that contains the replicate folders (N_1, N_2, ...). The folder names can have more text after the number, as in N_1_oct_6. |
| `output_root` | `''` | Folder where masks, ROIs, measurements and figures are saved. It has to be separate from the image folders, and a run stops if it lies inside one. |
| `replicates` | `()` | Replicates to analyse, for example ('N_1', 'N_2'). Leave empty to use every replicate found. |
| `image_glob` | `'*.tiff'` | Filename pattern of the exported multichannel TIFFs. TIFFs with fewer channels than the settings need are skipped. |
| `include_nd2` | `True` | Also look for ND2 files. When a position has also been exported as a TIFF, the TIFF is used for that position. |
| `denoised_tokens` | `(' - Denoised', '_denoised')` | Text that marks a filename as an NIS-denoised export. The original image has the same name with this text removed. |
| `nis_preference` | `'nd2'` | The denoised export that supplies the DAPI when a sample has both an ND2 and a TIFF version. If the preferred export is missing, the other one is used. Choices: `nd2`, `tiff`. |
| `export_csv` | `False` | Also save each sample's measurements as CSV, next to the Parquet file. |

## Phases

| Setting | Default | Description |
|---|---|---|
| `phase_files` | `{}` | Phase table for each replicate, for example {'N_1': '/path/phases.csv'}. A replicate not listed here uses phases.csv in its folder, if there is one. |
| `phase_overrides` | `{}` | Phases set by hand, for example {'N_2/triple_ihc_014_P1': '2'}, which take the place of the entries in the phase table. |
| `phase_average` | `'mean'` | Average used for each phase. Choices: `mean`, `median`. |
| `phase_error` | `'sem'` | Spread shown around each phase average: standard deviation, standard error, or interquartile range. Choices: `sd`, `sem`, `iqr`. |
| `phase_minimum_lines` | `2` | Fewest profiles needed at a distance for a phase average to be drawn there. |
| `phase_labels` | `('1', '2', '3')` | Phase labels to summarise and plot, in order. Samples with other labels or no label are shown with the individual profiles and left out of the phase averages. |
| `phase_compartments` | `{'Her7-Venus': 'nuclear', 'ppERK': 'nuclear', 'pMEK': 'cytoplasmic'}` | Compartment shown for each channel in the phase plots. Every compartment is measured, so this setting only changes what is plotted. |

## Channels

| Setting | Default | Description |
|---|---|---|
| `dapi_channel` | `1` | DAPI channel in the original image, counting from 1. An NIS export that contains only DAPI is read from its single channel. |
| `protein_channels` | `{'Her7-Venus': 2, 'pMEK': 3, 'ppERK': 4}` | Protein names and their channels in the original image, counting from 1. The protein signal is read from the original image, and NIS exports are used only for the DAPI. |

## Calibration

| Setting | Default | Description |
|---|---|---|
| `pixel_size_um` | `None` | Pixel size in µm. None reads it from the image file. |
| `z_step_um` | `None` | Distance between Z-planes in µm. None reads it from the image file. |

## DAPI

| Setting | Default | Description |
|---|---|---|
| `use_ai_denoise` | `True` | If True, the DAPI is taken from the sample's NIS-denoised export, so every sample needs one. If False, the raw DAPI is median filtered, background subtracted and contrast enhanced with CLAHE before segmentation. |

## Mask

| Setting | Default | Description |
|---|---|---|
| `use_cellpose` | `True` | If True, find nuclei with Cellpose. If False, threshold each Z-plane with the IsoData method. |
| `mask_source` | `'compute'` | Where the nuclear mask comes from. 'compute' makes it here, 'existing' reads the Mask.tif in the image's own folder, and 'edited' reads a mask corrected in Fiji and saved as edited_masks/<replicate>/<sample>/Mask.tif in the output folder. Choices: `compute`, `existing`, `edited`. |

## DAPI filtering

| Setting | Default | Description |
|---|---|---|
| `dapi_median_size` | `3` | Width in pixels of the median filter applied to the raw DAPI when use_ai_denoise is False. Use an odd number, or 0 to skip the filter. |
| `background_size_px` | `101` | Width in pixels of the square used to estimate the smooth background under the nuclei, which is then subtracted. The square should be wider than a nucleus, and 0 skips this step. |

## Contrast

| Setting | Default | Description |
|---|---|---|
| `clahe_window_px` | `31` | Approximate width in pixels of the tiles that CLAHE equalizes separately. |
| `clahe_clip` | `3.0` | CLAHE clip limit, as a multiple of the average histogram bin (the 'slope' in Fiji). Higher values increase local contrast more. |

## Cellpose

| Setting | Default | Description |
|---|---|---|
| `cellpose_model` | `'cpsam_v2'` | Cellpose 4 model name, or the path to a weights file. |
| `cellpose_diameter_px` | `None` | Typical nuclear diameter in pixels. When it is set, Cellpose rescales the image so the nuclei match the size it was trained on. None leaves the image at its own scale. |
| `cellpose_percentiles` | `(1.0, 99.0)` | DAPI percentiles that Cellpose maps to 0 and 1 before segmenting. |
| `cellprob_threshold` | `0.0` | Cellpose's threshold on its nucleus probability score. Lower values make masks larger. |
| `flow_threshold` | `0.4` | How closely a nucleus's shape must agree with Cellpose's predicted flows to be kept. Lower values are stricter. |
| `min_size_px` | `15` | Smallest nucleus kept, as an area in pixels within one Z-plane. |
| `max_size_fraction` | `0.4` | Largest nucleus kept, as a fraction of the image area. |
| `flow_iterations` | `None` | Number of steps Cellpose takes when following its flows to group pixels into nuclei. None uses Cellpose's default. |
| `tile_overlap` | `0.1` | Fraction by which Cellpose's 256-pixel tiles overlap. |
| `augment` | `False` | Also segment flipped copies of each tile and average the results, which takes longer. |

## Compute

| Setting | Default | Description |
|---|---|---|
| `batch_size` | `32` | Number of image tiles Cellpose sends to the GPU together. Larger batches are faster until GPU memory runs out. |
| `device` | `'auto'` | Where Torch runs. 'auto' uses an NVIDIA GPU (CUDA) if there is one, then an Apple GPU (MPS), then the CPU. Choices: `auto`, `cpu`, `mps`, `cuda`. |
| `image_workers` | `1` | Number of samples measured in parallel, each in its own process. The masks are made beforehand, one sample at a time, so that only one copy of the Cellpose model is on the GPU. |
| `threads_per_worker` | `8` | Threads each process uses for filtering Z-planes and for Torch. |
| `reuse_cache` | `True` | Load a sample's saved mask when one exists, instead of making it again. The settings and code that made the saved mask are not checked, so set this to False after changing either. Measurements are made again on every run, and RELOAD_SAVED in the notebook loads the saved ones instead. |

## ROI

| Setting | Default | Description |
|---|---|---|
| `width_override_px` | `None` | Band width in pixels to use for every side. None uses the width saved with each ROI. |
| `interpolation` | `'linear'` | How the band reads the image and the mask between pixel centres. Choices: `linear`, `nearest`. |
| `boundary_mode` | `'zero'` | What to do where a band runs off the image. 'zero' counts those points as 0, 'exclude' leaves them out of the average, and 'error' stops the run. Choices: `zero`, `exclude`, `error`. |

## Measurement

| Setting | Default | Description |
|---|---|---|
| `methods` | `('block_mip', 'per_slice', 'pooled', 'full_mip')` | Ways of handling Z to measure. 'per_slice' measures each plane, 'pooled' adds the planes together, 'block_mip' takes maximum projections of blocks of planes, and 'full_mip' projects the whole stack. |
| `z_start` | `1` | First Z-plane to use, counting from 1. |
| `z_stop` | `None` | Last Z-plane to use, counting from 1. None uses the last plane. |
| `mip_block_size` | `3` | Number of Z-planes in each block maximum projection. |
| `mip_start` | `2` | First Z-plane of the first block, counting from 1. With 31 planes, starting at 2 with blocks of 3 gives 2–4, 5–7, ..., 29–31. |
| `incomplete_block` | `'drop'` | What to do with a last block that has too few planes: drop it, include it, or stop with an error. Choices: `drop`, `include`, `error`. |
| `min_mask_fraction` | `0.05` | Smallest fraction of the band a compartment must cover for its mean to be reported. Where the coverage is lower, the mean is left empty (NaN). |

## Protein filtering

| Setting | Default | Description |
|---|---|---|
| `protein_median_size` | `3` | Width in pixels of the median filter applied to the protein channels. Use an odd number, or 0 to skip the filter. |
| `protein_outliers` | `True` | Replace pixels that are much brighter than their surroundings with the local median. |
| `outlier_radius_px` | `2.0` | Radius in pixels of the neighbourhood used to find bright outliers. |
| `outlier_threshold` | `50.0` | How far above the local median, in intensity units, a pixel must be to be replaced. |
| `mask_cap` | `4095.0` | Protein intensities above this value are clipped before the nuclear and cytoplasmic measurements. 4095 is the maximum of the 12-bit camera, and None turns clipping off. The whole-band measurement uses the unclipped values. |

## Exclusions

| Setting | Default | Description |
|---|---|---|
| `exclusions_csv` | `''` | Optional table of distance ranges to leave out of the analysis, with columns replicate, sample_id, side, start_um, stop_um, channels and reason. |

## Analysis

| Setting | Default | Description |
|---|---|---|
| `analysis_method` | `'block_mip'` | Z method used for the profiles and phase plots, which has to be one of those listed in methods. |
| `z_average` | `'pool'` | How the Z-planes or blocks of a side are combined. 'pool' weights each by its mask coverage, 'mean' weights them equally, and 'median' takes the middle value. Choices: `pool`, `mean`, `median`. |
| `combine_sides` | `False` | If True, combine an embryo's sides into one profile. If False, each side is its own profile. |
| `side_average` | `'mean'` | How sides are combined when combine_sides is True. 'mean' weights sides equally and 'pool' weights them by mask coverage. Choices: `mean`, `pool`. |

## Normalization

| Setting | Default | Description |
|---|---|---|
| `background_subtract` | `True` | Subtract a baseline from each embryo's profiles. The baseline is a low percentile of the profiles themselves rather than a separate measurement of the optical background. |
| `background_percentile` | `5.0` | Percentile of the smoothed side profiles that sets the embryo's baseline. |
| `normalization_window_um` | `40.0` | Width in µm of the rolling median that smooths the profiles before the baseline percentile is taken. The baseline is subtracted from the unsmoothed profile. |
| `day_normalize` | `False` | After background subtraction, divide each imaging day's profiles by that day's median value, separately for each channel and compartment. The imaging day is the ND2 acquisition date, or the replicate when there is no date. |
| `day_reference_um` | `(0.0, 400.0)` | Distance range in µm, measured from the posterior end, whose values set each day's median. |

## Statistics

| Setting | Default | Description |
|---|---|---|
| `average` | `'mean'` | Average used for the summary across all profiles. Choices: `mean`, `median`. |
| `error` | `'sd'` | Spread shown around the overall average: standard deviation, standard error, or interquartile range. Choices: `sd`, `sem`, `iqr`. |
| `minimum_lines` | `5` | Fewest profiles needed at a distance for the overall average to be drawn there. |

## Plot

| Setting | Default | Description |
|---|---|---|
| `smooth_average` | `True` | Smooth the plotted average with a Savitzky–Golay filter. The saved statistics and the error bands are calculated from the unsmoothed profiles. |
| `savgol_window_um` | `60.0` | Width in µm of the smoothing window. |
| `savgol_order` | `1` | Degree of the smoothing polynomial. |
| `plot_channels` | `('Her7-Venus', 'ppERK', 'pMEK')` | Channels to plot. Each must be a name in protein_channels. |
| `plot_compartments` | `('nuclear', 'cytoplasmic', 'all')` | Compartments to plot. |
| `phase_colormap` | `'viridis'` | Matplotlib colormap for the phases, which are given evenly spaced colours from it. |
| `phase_colormap_range` | `(0.0, 0.85)` | Part of the colormap to use, between 0 and 1. The default skips the pale yellow end of viridis, which is hard to see on a light background. |
| `show_individuals` | `True` | Draw each profile behind the phase average in the per-phase panels. |
| `dpi` | `150` | Resolution of saved figures, in dots per inch. |
| `qc_percentiles` | `(1.0, 99.7)` | Intensity percentiles that set the display contrast of the QC figure. |
