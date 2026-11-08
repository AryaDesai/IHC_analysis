"""Tests of the measurements, file formats and analysis steps, which run without a GPU."""

from dataclasses import replace
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import tifffile
import ihc_pipeline as ihc


def cfg(tmp_path, **changes):
    base = dict(
        experimental_data=tmp_path,
        output_root=tmp_path / "results",
        include_nd2=False,
        protein_channels={"p": 2},
        plot_channels=("p",),
        phase_compartments={"p": "nuclear"},
        protein_median_size=0,
        protein_outliers=False,
        mask_cap=None,
        methods=("per_slice", "pooled"),
        analysis_method="per_slice",
        min_mask_fraction=0.05,
        background_subtract=False,
        smooth_average=False,
        phase_minimum_lines=1,
    )
    base.update(changes)
    return ihc.Config(**base)


def test_roi_roundtrip(tmp_path):
    r = ihc.ROI("side0", np.array([[2.25, 3.5], [8.2, 4.1], [12.6, 9.0]]), 12)
    path = tmp_path / "rois.csv"
    ihc.save_rois(path, [r])
    loaded = ihc.load_rois(path)[0]
    assert np.allclose(loaded.xy, r.xy) and loaded.width_px == 12
    assert np.allclose(ihc.band_grid(loaded)[0], ihc.band_grid(r)[0])


def test_concentration_is_not_occupancy(tmp_path):
    c = cfg(tmp_path)
    roi = ihc.ROI("side0", np.array([[5.0, 10.0], [20.0, 10.0]]), 10)
    image = np.full((1, 4, 24, 30), 100.0, np.float32)
    results = []
    for cutoff in [8, 10]:
        mask = np.zeros((4, 24, 30), bool)
        mask[:, :cutoff] = True
        t = ihc.measure_arrays(image, mask, [roi], {"pixel_size_um": 1.0}, c)
        r = t[(t.method == "pooled") & (t.compartment == "nuclear")]
        assert np.allclose(r.masked_mean, 100)
        results.append(r.band_mean.to_numpy())
    assert np.all(results[0] < results[1])


def test_voxel_partition_and_pooling(tmp_path):
    c = cfg(tmp_path)
    roi = ihc.ROI("side0", np.array([[4.0, 8.0], [10.0, 8.0]]), 5)
    rng = np.random.default_rng(1)
    image = rng.uniform(0, 5000, (1, 4, 18, 18)).astype(np.float32)
    mask = rng.random((4, 18, 18)) > 0.4
    t = ihc.measure_arrays(image, mask, [roi], {"pixel_size_um": 1.0}, c)
    for method in c.methods:
        p = t[t.method == method].pivot(
            index=["z_start", "position_index"], columns="compartment", values="masked_sum"
        )
        assert np.allclose(p.nuclear + p.cytoplasmic, p["all"], rtol=1e-6)
    single = (
        t[(t.method == "per_slice") & (t.compartment == "nuclear")]
        .groupby("position_index")[["masked_sum", "mask_weight"]]
        .sum()
    )
    pooled = t[(t.method == "pooled") & (t.compartment == "nuclear")].set_index("position_index")
    assert np.allclose(single.masked_sum, pooled.masked_sum)
    assert np.allclose(single.mask_weight, pooled.mask_weight)


def test_mip_masks_before_projection(tmp_path):
    c = cfg(
        tmp_path, methods=("block_mip",), analysis_method="block_mip", mip_start=1, mip_block_size=2
    )
    roi = ihc.ROI("side0", np.array([[3.0, 5.0], [6.0, 5.0]]), 1)
    image = np.full((1, 2, 12, 12), 200.0, np.float32)
    image[:, 0] = 10
    mask = np.zeros((2, 12, 12), bool)
    mask[0] = True
    t = ihc.measure_arrays(image, mask, [roi], {"pixel_size_um": 1.0}, c)
    assert np.allclose(t[t.compartment == "nuclear"].masked_mean, 10)
    assert np.allclose(t[t.compartment == "cytoplasmic"].masked_mean, 200)


def test_missing_mask_support_and_bounds(tmp_path):
    c = cfg(tmp_path, boundary_mode="exclude")
    roi = ihc.ROI("side0", np.array([[3.0, 0.0], [6.0, 0.0]]), 5)
    image = np.full((1, 4, 12, 12), 30.0, np.float32)
    mask = np.zeros((4, 12, 12), bool)
    t = ihc.measure_arrays(image, mask, [roi], {"pixel_size_um": 1.0}, c)
    assert t[t.compartment == "nuclear"].masked_mean.isna().all()
    assert np.allclose(t[t.compartment == "all"].masked_mean, 30)
    with pytest.raises(ValueError):
        ihc.measure_arrays(
            image, mask, [roi], {"pixel_size_um": 1.0}, replace(c, boundary_mode="error")
        )


def test_fiji_mask_interchange(tmp_path):
    a = np.zeros((4, 12, 14), bool)
    a[:, 2:5, 3:7] = True
    p = tmp_path / "Mask.tif"
    ihc.write_mask(p, a, 0.35825, 0.423)
    with tifffile.TiffFile(p) as f:
        assert f.series[0].axes == "ZYX" and f.imagej_metadata["unit"] == "um"
        b = f.asarray()
        assert set(np.unique(b)) == {0, 4095}
    assert np.array_equal(b > 0, a)


def test_discovery_and_phase_aliases(tmp_path):
    for rep in ["N_1_date", "N_10_date"]:
        p = tmp_path / rep / "tiffs"
        p.mkdir(parents=True)
        a = np.ones((4, 2, 12, 14), np.uint16)
        for name in ["example_P1.tiff", "example - Denoised_P1.tiff"]:
            tifffile.imwrite(
                p / name,
                a,
                imagej=True,
                resolution=(2, 2),
                metadata={"axes": "ZCYX", "unit": "um", "spacing": 1},
            )
        (p.parent / "phases.csv").write_text("file,phase\nexample.nd2,2\n")
    c = cfg(tmp_path)
    s = ihc.discover_samples(c)
    assert list(s.replicate) == ["N_1", "N_10"] and len(s) == 2
    assert s.nis_path.str.contains("Denoised").all() and (s.phase == "2").all()
    assert s.embryo_id.nunique() == 2
    # N_ folders inside the output root are not replicates.
    (Path(c.output_root) / "masks" / "N_1").mkdir(parents=True)
    assert len(ihc.discover_samples(c)) == 2
    # Nor are N_ folders inside a replicate.
    (tmp_path / "N_1_date" / "tiffs_N_1" / "rois" / "N_1").mkdir(parents=True)
    assert len(ihc.discover_samples(c)) == 2
    # Two experiment folders with the same replicate number are an error.
    (tmp_path / "N_1_other_date").mkdir()
    with pytest.raises(ValueError, match="two N_1 folders"):
        ihc.discover_samples(c)


@pytest.mark.parametrize("preference", ["nd2", "tiff"])
def test_nis_preference_keeps_original_and_position(tmp_path, monkeypatch, preference):
    import nd2

    folder = tmp_path / "N_1_date"
    tiffs = folder / "tiffs"
    tiffs.mkdir(parents=True)
    raw_nd2 = folder / "example.nd2"
    nis_nd2 = folder / "example - Denoised.nd2"
    raw_nd2.touch()
    nis_nd2.touch()
    raw_tiff = tiffs / "example_P2.tiff"
    nis_tiff = tiffs / "example - Denoised_P2.tiff"
    for path in [raw_tiff, nis_tiff]:
        tifffile.imwrite(path, np.ones((3, 2, 12, 14), np.uint16),
                         imagej=True, metadata={"axes": "ZCYX"})

    class TwoPositionND2:
        def __init__(self, path):
            self.sizes = {"P": 2}

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

    monkeypatch.setattr(nd2, "ND2File", TwoPositionND2)
    c = cfg(tmp_path, include_nd2=True, nis_preference=preference)
    samples = ihc.discover_samples(c).set_index("sample_id")
    assert len(samples) == 2
    assert samples.loc["example_P1", "raw_path"] == str(raw_nd2)
    assert samples.loc["example_P1", "nis_path"] == str(nis_nd2)
    p2 = samples.loc["example_P2"]
    assert p2.raw_path == str(raw_tiff)  # the proteins are read from the original image
    assert p2.nis_path == str(nis_nd2 if preference == "nd2" else nis_tiff)
    assert p2.nis_position == (2 if preference == "nd2" else 1)

    # Without a denoised ND2, both settings use the denoised TIFF.
    nis_nd2.unlink()
    fallback = ihc.discover_samples(c).set_index("sample_id").loc["example_P2"]
    assert fallback.nis_path == str(nis_tiff) and fallback.nis_position == 1


def test_side_choice_and_phase_counts(tmp_path):
    c = cfg(tmp_path)
    records = []
    for rep, value in [("N_1", 10.0), ("N_2", 30.0)]:
        for side in ["side0", "side1"]:
            for z in [1, 2]:
                for x in range(3):
                    records.append(
                        dict(
                            replicate=rep,
                            sample_id="a",
                            embryo_id=rep + "/a",
                            side=side,
                            channel="p",
                            compartment="nuclear",
                            method="per_slice",
                            z_start=z,
                            z_stop=z,
                            position_index=x,
                            distance_um=float(x),
                            masked_sum=value * 5,
                            mask_weight=5.0,
                            band_weight=10.0,
                            masked_mean=value,
                            mask_fraction=0.5,
                            phase="1",
                            phase_source="test",
                        )
                    )
    table = pd.DataFrame(records)
    for combined, n in [(False, 4), (True, 2)]:
        current = replace(c, combine_sides=combined)
        profiles, _ = ihc.side_profiles(table, current)
        summary = ihc.summarize(profiles, current)
        assert (
            (summary.n == n).all()
            and (summary.n_embryos == 2).all()
            and (summary.n_replicates == 2).all()
        )
        assert np.allclose(summary.center, 20)


def test_day_normalization_removes_day_gain(tmp_path):
    # One replicate imaged on two days with the same profile shape. The second day
    # has twice the gain and a different floor, and each day gets its own factor.
    shape = np.array([0.0, 1, 3, 6, 4, 2, 1, 0])
    records = []
    for day, gain, floor in [("2026-10-06", 1.0, 50.0), ("2026-10-08", 2.0, 60.0)]:
        for x, s in enumerate(shape):
            value = floor + gain * s
            records.append(
                dict(replicate="N_1", imaging_day=day, sample_id=day, embryo_id="N_1/" + day, side="side0",
                     channel="p", compartment="nuclear", method="per_slice", z_start=1,
                     z_stop=1, position_index=x, distance_um=float(x), masked_sum=value * 5,
                     mask_weight=5.0, band_weight=10.0, masked_mean=value, mask_fraction=0.5,
                     phase="1", phase_source="test")
            )
    c = cfg(tmp_path, background_subtract=True, background_percentile=0.0,
            normalization_window_um=1.0, day_normalize=True, day_reference_um=(0.0, 100.0))
    profiles, normalization = ihc.side_profiles(pd.DataFrame(records), c)
    days = profiles.pivot(index="position_index", columns="imaging_day", values="value")
    assert np.allclose(days["2026-10-06"], days["2026-10-08"])  # floor and gain both removed
    assert np.allclose(normalization.set_index("imaging_day").day_factor, [1.5, 3.0])  # median of shape x gain


def test_unknown_or_invalid_parameters(tmp_path):
    with pytest.raises(TypeError):
        ihc.Config(nonexistent=1)
    with pytest.raises(ValueError):
        cfg(tmp_path, methods=("unknown",)).validate()
    with pytest.raises(ValueError):
        cfg(tmp_path, width_override_px=0).validate()


def test_mask_cap_applies_to_compartments_only(tmp_path):
    roi = ihc.ROI("side0", np.array([[4.0, 8.0], [10.0, 8.0]]), 3)
    values = np.full((1, 4, 18, 18), 6000.0, np.float32)
    mask = np.ones((4, 18, 18), bool)
    c = cfg(tmp_path, mask_cap=4095.0)
    t = ihc.measure_arrays(values, mask, [roi], {"pixel_size_um": 1.0}, c)
    assert np.allclose(t[t.compartment == "nuclear"].masked_mean, 4095.0)
    assert np.allclose(t[t.compartment == "all"].masked_mean, 6000.0)
    t = ihc.measure_arrays(
        values, mask, [roi], {"pixel_size_um": 1.0}, replace(c, mask_cap=None)
    )
    assert np.allclose(t[t.compartment == "nuclear"].masked_mean, 6000.0)


def test_measurement_records_its_settings(tmp_path):
    rawdir = tmp_path / "N_1_test" / "tiffs"
    rawdir.mkdir(parents=True)
    path = rawdir / "sample_P1.tiff"
    values = np.full((4, 2, 24, 30), 100, np.uint16)
    tifffile.imwrite(
        path,
        values,
        imagej=True,
        resolution=(2, 2),
        metadata={"axes": "ZCYX", "unit": "um", "spacing": 1.0},
    )
    mask = np.zeros((4, 24, 30), bool)
    mask[:, :12] = True
    path.with_suffix("").mkdir()
    ihc.write_mask(path.with_suffix("") / "Mask.tif", mask, 0.5, 1.0)
    ihc.save_rois(
        path.with_suffix("") / "rois.csv",
        [ihc.ROI("side0", np.array([[5.0, 12.0], [20.0, 12.0]]), 6)],
    )
    c = cfg(tmp_path, mask_source="existing", use_cellpose=False, device="cpu")
    s = ihc.discover_samples(c).iloc[0]
    metadata = ihc.sample_dir(c, s, "measurements") / "measurement.json"

    ihc.run_sample(s, replace(c, width_override_px=12))
    settings = json.loads(metadata.read_text())["provenance"]["settings"]
    assert settings["width_override_px"] == 12
    assert set(settings) == set(ihc.MEASURE_SETTINGS)  # the analysis and plot settings are left out


def test_band_sampler_matches_map_coordinates():
    """Check the band sampler with both interpolation orders, using a band that runs off the image."""
    from scipy import ndimage as ndi

    image = np.random.default_rng(0).random((20, 30)).astype(np.float32)
    x, y = np.meshgrid(np.linspace(-3, 33, 50), np.linspace(-2.5, 21.7, 7))
    grid = np.stack([x, y], -1)  # (band width, length, xy), partly outside the image
    for order in (0, 1):
        idx, wts = ihc._band_sampler(grid, *image.shape, order)
        mine = (image.ravel()[idx] * wts).sum(-1).reshape(x.shape)
        ref = ndi.map_coordinates(image, [y, x], order=order, mode="constant", cval=0, prefilter=False)
        assert np.allclose(mine, ref, atol=1e-6)


def test_run_parallel_workers_match_single(tmp_path):
    rawdir = tmp_path / "N_1_test" / "tiffs"
    rawdir.mkdir(parents=True)
    rng = np.random.default_rng(1)
    for name in ("a_P1", "b_P1"):
        path = rawdir / f"{name}.tiff"
        tifffile.imwrite(
            path,
            rng.integers(0, 4000, (4, 2, 24, 30), dtype=np.uint16),
            imagej=True,
            resolution=(2, 2),
            metadata={"axes": "ZCYX", "unit": "um", "spacing": 1.0},
        )
        path.with_suffix("").mkdir()
        ihc.write_mask(path.with_suffix("") / "Mask.tif", rng.random((4, 24, 30)) > 0.5, 0.5, 1.0)
        ihc.save_rois(
            path.with_suffix("") / "rois.csv",
            [ihc.ROI("side0", np.array([[5.0, 12.0], [20.0, 12.0]]), 6)],
        )
    c = cfg(tmp_path, mask_source="existing", use_cellpose=False, device="cpu")
    samples = ihc.discover_samples(c)
    one = ihc.run(samples, c)
    two = ihc.run(samples, replace(c, image_workers=2))
    pd.testing.assert_frame_equal(one, two)


def test_filtered_dapi_threshold_finds_nuclei_on_a_slope(tmp_path):
    # Two bright "nuclei" on a background that rises across the image. Background
    # subtraction should flatten the slope so the threshold finds just the nuclei.
    # Within half a background width of the image edge the opening cannot follow
    # the slope, so the background check leaves that margin out.
    y, x = np.mgrid[:120, :160]
    plane = 200 + 3.0 * x
    plane[30:40, 20:30] += 800
    plane[70:80, 120:130] += 800
    dapi = np.stack([plane, plane]).astype(np.uint16)
    c = cfg(tmp_path, background_size_px=31, clahe_window_px=40, threads_per_worker=1)
    filtered = ihc._filter_dapi(dapi, c)
    from skimage.filters import threshold_isodata

    mask = filtered[0] > threshold_isodata(filtered[0])
    assert mask[30:40, 20:30].mean() > 0.9 and mask[70:80, 120:130].mean() > 0.9
    assert mask[15:-15, 15:-15].mean() < 0.02


def test_ai_denoise_needs_an_nis_export(tmp_path):
    rawdir = tmp_path / "N_1_test" / "tiffs"
    rawdir.mkdir(parents=True)
    tifffile.imwrite(rawdir / "sample_P1.tiff", np.ones((2, 2, 12, 14), np.uint16), imagej=True,
                     resolution=(2, 2), metadata={"axes": "ZCYX", "unit": "um"})
    c = cfg(tmp_path, use_ai_denoise=True, reuse_cache=False, device="cpu")
    s = ihc.discover_samples(c).iloc[0]
    with pytest.raises(FileNotFoundError, match="no NIS-denoised file"):
        ihc.segment_sample(s, c)


def test_mask_folders_are_separate_for_each_method(tmp_path):
    rawdir = tmp_path / "N_1_test" / "tiffs"
    rawdir.mkdir(parents=True)
    plane = np.zeros((2, 2, 40, 40), np.uint16)
    plane[:, 0, 10:20, 10:20] = 1000
    tifffile.imwrite(rawdir / "sample_P1.tiff", plane, imagej=True, resolution=(2, 2),
                     metadata={"axes": "ZCYX", "unit": "um"})
    c = cfg(tmp_path, use_ai_denoise=False, use_cellpose=False, background_size_px=0,
            clahe_window_px=20, device="cpu")
    s = ihc.discover_samples(c).iloc[0]
    ihc.segment_sample(s, c)
    assert (ihc.sample_dir(c, s, "masks_threshold_raw") / "Mask.tif").exists()
    assert ihc.mask_kind(replace(c, use_cellpose=True, use_ai_denoise=True)) == "masks_cellpose_denoised"


def test_git_state_records_commit_and_changed_files(tmp_path):
    import subprocess
    def git(*args):
        subprocess.run(["git", "-C", str(tmp_path), *args], check=True, capture_output=True)
    git("init", "-q")
    (tmp_path / "code.py").write_text("a = 1\n")
    git("add", "code.py")
    git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "start")
    assert ihc.git_state(tmp_path)["git_uncommitted"] == []
    (tmp_path / "code.py").write_text("a = 2\n")
    state = ihc.git_state(tmp_path)
    assert len(state["git_commit"]) == 40 and state["git_uncommitted"] == ["code.py"]
    assert ihc.git_state(tmp_path / "not_a_repo") == dict(git_commit=None, git_uncommitted=None)
