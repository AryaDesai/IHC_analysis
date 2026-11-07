"""napari editor for drawing ROIs and correcting nuclear masks.

Open it from the notebook with launch_editor, or run
    python ihc_roi_editor.py --config /path/config.json
Each side of the embryo is a separate Shapes layer with its own band width. The
orange bands are computed by the pipeline's own band_grid, so they show exactly
the pixels that will be measured.
"""

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from ihc_pipeline import (
    Config,
    ROI,
    _read_mask,
    _roi_sources,
    band_grid,
    discover_samples,
    file_identity,
    load_rois,
    read_image,
    sample_dir,
    save_rois,
    write_mask,
)


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--config", required=True)
    p.add_argument("--samples", help="Optional resolved sample JSON written by launch_editor")
    p.add_argument(
        "--smoke",
        action="store_true",
        help="Load the first sample, save its ROIs and mask, and exit without showing a window",
    )
    args = p.parse_args()
    config = Config(**json.loads(Path(args.config).read_text())).validate()
    # napari saves its preferences under output_root rather than the home folder.
    for kind in ("CACHE", "CONFIG", "DATA", "STATE"):
        os.environ[f"XDG_{kind}_HOME"] = str(
            Path(config.output_root) / "viewer_settings" / kind.lower()
        )
    import napari
    from qtpy.QtWidgets import (
        QWidget,
        QVBoxLayout,
        QLabel,
        QPushButton,
        QComboBox,
        QSpinBox,
        QMessageBox,
    )

    samples = (
        pd.DataFrame(json.loads(Path(args.samples).read_text()))
        if args.samples
        else discover_samples(config)
    )

    class Editor(QWidget):
        def __init__(self, viewer):
            super().__init__()
            self.viewer = viewer
            self.layers = []
            self.dirty = False
            self.loading = False
            layout = QVBoxLayout(self)
            self.sample_box = QComboBox()
            self.sample_box.addItems(
                [f"{s.replicate}/{s.sample_id}" for _, s in samples.iterrows()]
            )
            layout.addWidget(self.sample_box)
            self.status = QLabel()
            self.status.setWordWrap(True)
            layout.addWidget(self.status)
            self.source = QComboBox()
            self.source.addItems(["Original channels", "NIS DAPI (if available)"])
            layout.addWidget(self.source)
            self.source.currentIndexChanged.connect(self.change_source)
            self.channel = QComboBox()
            layout.addWidget(self.channel)
            self.channel.currentIndexChanged.connect(self.change_channel)
            self.width = QSpinBox()
            self.width.setRange(1, 10000)
            layout.addWidget(QLabel("Selected ROI width (pixels)"))
            layout.addWidget(self.width)
            self.width.valueChanged.connect(self.change_width)
            for label, action in [
                ("Load selected image / reload ROIs", self.load),
                ("New side [N]", self.new_side),
                ("Edit selected side [E]", self.edit),
                ("Delete selected side", self.delete_side),
                ("Preview actual measurement bands", self.preview),
                ("Save ROIs [S]", self.save),
                ("Save ROIs and next image", self.next_image),
                ("Reload mask from disk", self.load_mask),
                ("Save edited nuclear mask", self.save_mask),
            ]:
                button = QPushButton(label)
                button.clicked.connect(action)
                layout.addWidget(button)
            help_text = (
                "Draw each side from posterior to anterior and double-click to finish.\n"
                "Each side has its own layer. Select it before changing its width.\n"
                "An ROI applies to every Z-plane.\n"
                "To correct the mask, select the 'Nuclear mask' layer and paint 1 for "
                "nucleus or 0 for background. Edited masks are saved separately and are "
                "measured when mask_source='edited'."
            )
            if config.width_override_px is not None:
                help_text += (
                    f"\nwidth_override_px is {config.width_override_px}, so the preview shows "
                    "that width. The width control still changes the width saved with the ROI."
                )
            layout.addWidget(QLabel(help_text))
            viewer.bind_key("n", lambda v: self.new_side(), overwrite=True)
            viewer.bind_key("e", lambda v: self.edit(), overwrite=True)
            viewer.bind_key("s", lambda v: self.save(), overwrite=True)
            viewer.layers.selection.events.active.connect(self.selection_changed)
            self.load()

        def guard_unsaved(self):
            if not self.dirty:
                return True
            answer = QMessageBox.question(
                self,
                "Unsaved ROI edits",
                "Save ROI edits before loading another image?",
                QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
            )
            if answer == QMessageBox.Cancel:
                return False
            if answer == QMessageBox.Save:
                return self.save()
            return True

        def load(self):
            if hasattr(self, "sample") and not self.guard_unsaved():
                return
            self.loading = True
            self.sample = samples.iloc[self.sample_box.currentIndex()].to_dict()
            self.raw, self.cal = read_image(
                self.sample["raw_path"], int(self.sample.get("raw_position", 1)), config
            )
            self.view_data = self.raw
            self.viewer.layers.clear()
            self.layers = []
            self.image = self.viewer.add_image(
                self.raw[:, 0], name="Original image", colormap="gray"
            )
            self.bands = self.viewer.add_shapes(
                ndim=2,
                name="Measurement bands",
                shape_type="polygon",
                face_color="orange",
                edge_color="orange",
                opacity=config.viewer_opacity,
            )
            self.bands.editable = False
            self.channel.blockSignals(True)
            self.channel.clear()
            self.channel.addItems([f"Channel {i+1}" for i in range(self.raw.shape[1])])
            self.channel.setCurrentIndex(min(config.viewer_channel - 1, self.raw.shape[1] - 1))
            self.channel.blockSignals(False)
            self.change_channel(self.channel.currentIndex())
            self.viewer.dims.set_current_step(0, len(self.raw) // 2)
            self.roi_path = sample_dir(config, self.sample, "rois") / "rois.csv"
            self.roi_stamp = file_identity(self.roi_path) if self.roi_path.exists() else {}
            for r in load_rois(_roi_sources(self.sample, config)):
                self.add_roi(r)
            self.mask_layer = None
            self.load_mask()
            self.preview()
            if self.layers:
                self.viewer.layers.selection.active = self.layers[0]
            self.dirty = False
            self.loading = False
            self.status.setText(
                f"{self.sample['replicate']}/{self.sample['sample_id']}\nROIs: {self.roi_path}"
            )

        def change_channel(self, index):
            if index < 0 or not hasattr(self, "raw"):
                return
            a = self.view_data[:, min(index, self.view_data.shape[1] - 1)]
            self.image.data = a
            lo, hi = np.percentile(a[len(a) // 2], config.viewer_percentiles)
            self.image.contrast_limits = (float(lo), float(max(hi, lo + 1)))

        def change_source(self, index):
            if not hasattr(self, "sample"):
                return
            if index == 1 and self.sample.get("nis_path"):
                a, _ = read_image(
                    self.sample["nis_path"], int(self.sample.get("nis_position", 1)), config
                )
                ch = 0 if a.shape[1] == 1 else config.dapi_channel - 1
                self.view_data = a[:, ch : ch + 1]
            else:
                self.view_data = self.raw
            self.image.name = (
                "NIS DAPI" if index == 1 and self.sample.get("nis_path") else "Original image"
            )
            self.change_channel(self.channel.currentIndex())

        def mark_dirty(self, event=None):
            if not self.loading:
                self.dirty = True

        def add_roi(self, roi=None):
            used = {l.metadata["roi_id"] for l in self.layers}
            i = 0
            while f"side{i}" in used:
                i += 1
            name = roi.roi_id if roi else f"side{i}"
            width = roi.width_px if roi else config.new_roi_width_px
            layer = self.viewer.add_shapes(
                [roi.xy[:, ::-1]] if roi else [],
                ndim=2,
                name=name,
                shape_type="path",
                edge_color="yellow",
                edge_width=2,
                face_color="transparent",
            )
            layer.metadata.update(roi_id=name, width_px=width)
            layer.events.data.connect(self.mark_dirty)
            self.layers.append(layer)
            self.viewer.layers.selection.active = layer
            layer.mode = "direct" if roi else "add_polyline"
            return layer

        def new_side(self):
            self.add_roi()
            self.dirty = True

        def edit(self):
            layer = self.viewer.layers.selection.active
            if layer in self.layers:
                layer.mode = "direct"

        def delete_side(self):
            layer = self.viewer.layers.selection.active
            if layer in self.layers:
                self.layers.remove(layer)
                self.viewer.layers.remove(layer)
                self.dirty = True
                self.preview()

        def selection_changed(self, event=None):
            layer = self.viewer.layers.selection.active
            if layer in self.layers:
                self.width.blockSignals(True)
                self.width.setValue(int(round(layer.metadata["width_px"])))
                self.width.blockSignals(False)

        def change_width(self, value):
            layer = self.viewer.layers.selection.active
            if layer in self.layers:
                layer.metadata["width_px"] = value
                self.dirty = True
                self.preview()

        def rois(self):
            rois = []
            for layer in self.layers:
                if len(layer.data) != 1:
                    raise ValueError(
                        f"{layer.name} must hold a single line. Use New side to draw another side."
                    )
                rois.append(
                    ROI(
                        layer.metadata["roi_id"],
                        np.asarray(layer.data[0])[:, ::-1],
                        layer.metadata["width_px"],
                    ).validate()
                )
            return rois

        def preview(self):
            polygons = []
            for layer in self.layers:
                if len(layer.data) != 1 or len(layer.data[0]) < 2:
                    continue
                r = ROI(
                    layer.metadata["roi_id"],
                    np.asarray(layer.data[0])[:, ::-1],
                    layer.metadata["width_px"],
                )
                g, _ = band_grid(r, config.width_override_px)
                polygons.append(np.vstack([g[0], g[-1, ::-1]])[:, ::-1])
            self.bands.data = polygons
            if polygons:
                self.bands.shape_type = ["polygon"] * len(polygons)

        def save(self):
            try:
                save_rois(self.roi_path, self.rois(), expected=self.roi_stamp)
                self.roi_stamp = file_identity(self.roi_path)
                self.dirty = False
                self.preview()
                self.status.setText(f"Saved {self.roi_path}")
                return True
            except Exception as error:
                self.status.setText(str(error))
                return False

        def next_image(self):
            if self.save() and self.sample_box.currentIndex() + 1 < len(samples):
                self.sample_box.setCurrentIndex(self.sample_box.currentIndex() + 1)
                self.load()

        def load_mask(self):
            if not hasattr(self, "raw"):
                return
            paths = [
                sample_dir(config, self.sample, "edited_masks") / "Mask.tif",
                sample_dir(config, self.sample, "masks") / "Mask.tif",
                Path(self.sample["raw_path"]).with_suffix("") / "Mask.tif",
            ]
            if self.sample.get("nis_path"):
                paths.append(Path(self.sample["nis_path"]).with_suffix("") / "Mask.tif")
            source = next((p for p in paths if p.exists()), None)
            shape = (len(self.raw), *self.raw.shape[-2:])
            a = _read_mask(source) if source else np.zeros(shape, bool)
            if a.shape != shape:
                raise ValueError(f"{source} does not have the image's Z, Y, X shape")
            if self.mask_layer is not None:
                self.viewer.layers.remove(self.mask_layer)
            self.mask_layer = self.viewer.add_labels(
                a.astype(np.uint8), name="Nuclear mask", opacity=0.3
            )
            self.mask_path = sample_dir(config, self.sample, "edited_masks") / "Mask.tif"
            self.mask_stamp = file_identity(self.mask_path) if self.mask_path.exists() else {}

        def save_mask(self):
            current = file_identity(self.mask_path) if self.mask_path.exists() else {}
            if current != self.mask_stamp:
                self.status.setText("The edited mask changed on disk. Reload it before saving.")
                return
            mask = np.asarray(self.mask_layer.data) > 0
            write_mask(self.mask_path, mask, self.cal["pixel_size_um"], self.cal["z_step_um"])
            write_mask(
                self.mask_path.with_name("Cyto_Mask.tif"),
                ~mask,
                self.cal["pixel_size_um"],
                self.cal["z_step_um"],
            )
            self.mask_stamp = file_identity(self.mask_path)
            self.status.setText(f"Saved {self.mask_path}\nSet mask_source='edited' to measure it.")

    viewer = napari.Viewer(show=not args.smoke)
    editor = Editor(viewer)
    viewer.window.add_dock_widget(editor, name="IHC ROIs and masks", area="right")
    if args.smoke:
        assert editor.layers and editor.rois()
        editor.preview()
        editor.change_channel(min(1, editor.raw.shape[1] - 1))
        assert editor.save()
        editor.save_mask()
        viewer.close()
        print("Editor check passed: loaded ROIs, drew bands, switched channel, saved ROIs and mask")
    else:
        napari.run()


if __name__ == "__main__":
    main()
