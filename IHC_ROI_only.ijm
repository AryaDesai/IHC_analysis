// IHC_ROI_only.ijm draws the PSM lines of each image in Fiji. The masks and
// measurements are made afterwards in the notebook.
//
// The lines of each image are saved as rois/<replicate>/<sample>/rois.csv in the
// notebook's output folder, which is where the notebook reads them. Each row of the
// CSV is one clicked point (roi_id, vertex_index, x_px, y_px, width_px) in image
// pixels counted from 0, starting at the posterior end of the side.
//
// The images can be ND2 files or exported TIFFs. Each position of an ND2 file is
// opened separately and named <file>_P<position>, the name the notebook gives that
// sample. NIS-denoised exports are skipped, and the images themselves are not saved.

requires("1.54");
input = getDirectory("Folder with the images (for ND2 files, the replicate folder)");
output = getDirectory("Output folder of the notebook (ihc_results)");

// The replicate is the N_<number> part of a folder name in the input path.
replicate = "N_1";
parts = split(replace(input, "\\", "/"), "/");
for (p = 0; p < parts.length; p++) {
	if (matches(parts[p], "N_[0-9]+.*")) {
		fields = split(parts[p], "_");
		replicate = "N_" + fields[1];
	}
}

Dialog.create("IHC ROI settings");
Dialog.addChoice("Image files:", newArray("ND2", "TIFF"));
Dialog.addString("Replicate:", replicate);
Dialog.addNumber("Line width for new lines (pixels):", 100);
Dialog.addCheckbox("Redo images that already have ROIs", false);
Dialog.show();
type = Dialog.getChoice();
replicate = Dialog.getString();
lineWidth = Dialog.getNumber();
redo = Dialog.getCheckbox();
if (!matches(replicate, "N_[0-9]+")) exit("The replicate should look like N_1 or N_2.");
if (lineWidth < 1) exit("The line width should be at least 1 pixel.");

// List the images. An ND2 file adds one entry for each of its positions.
files = getFileList(input);
Array.sort(files);
paths = newArray(0);
series = newArray(0);
samples = newArray(0);
if (type == "ND2") run("Bio-Formats Macro Extensions");
for (i = 0; i < files.length; i++) {
	name = files[i];
	lower = toLowerCase(name);
	if (indexOf(lower, "denoised") >= 0) continue;
	if (type == "ND2") {
		if (!endsWith(lower, ".nd2")) continue;
		Ext.setId(input + name);
		Ext.getSeriesCount(positions);
		Ext.close();
		stem = substring(name, 0, lengthOf(name) - 4);
		for (pos = 1; pos <= positions; pos++) {
			paths = Array.concat(paths, input + name);
			series = Array.concat(series, pos);
			samples = Array.concat(samples, stem + "_P" + pos);
		}
	} else {
		if (!endsWith(lower, ".tiff")) continue;
		paths = Array.concat(paths, input + name);
		series = Array.concat(series, 0);
		samples = Array.concat(samples, substring(name, 0, lengthOf(name) - 5));
	}
}
total = paths.length;
if (total == 0) exit("No " + type + " images were found in\n" + input);

File.makeDirectory(output + "rois/");
File.makeDirectory(output + "rois/" + replicate + "/");
luts = newArray("Blue", "Green", "Red", "Magenta"); // display colours for channels 1 to 4

for (im = 0; im < total; im++) {
	sample = samples[im];
	folder = output + "rois/" + replicate + "/" + sample + "/";
	csvPath = folder + "rois.csv";
	if (File.exists(csvPath) && !redo) continue;
	original = "";
	if (File.exists(csvPath)) original = File.openAsString(csvPath);

	if (series[im] > 0) {
		run("Bio-Formats Importer", "open=[" + paths[im] + "] color_mode=Composite view=Hyperstack stack_order=XYCZT series_" + series[im]);
	} else {
		open(paths[im]);
	}
	traceID = getImageID();
	Stack.getDimensions(imageWidth, imageHeight, channels, slices, frames);
	if (channels > 1) Stack.setDisplayMode("composite");
	Stack.setSlice(floor(slices / 2) + 1);
	for (ch = 1; ch <= channels; ch++) {
		Stack.setChannel(ch);
		run(luts[(ch - 1) % 4]);
		run("Enhance Contrast", "saturated=0.35"); // changes the display, not the pixel values
	}
	Stack.setChannel(1);

	Dialog.create("How is this sample mounted?");
	Dialog.addMessage("Image " + (im + 1) + " of " + total + ": " + sample);
	Dialog.addChoice("Mounting:", newArray("Lateral", "Flat", "Skip"));
	beep();
	Dialog.show();
	mount = Dialog.getChoice();
	sides = 0;
	if (mount == "Flat") sides = 2;
	if (mount == "Lateral") sides = 1;
	if (sides == 0) {
		selectImage(traceID);
		close();
		continue;
	}

	result = "roi_id,vertex_index,x_px,y_px,width_px\n";
	for (side = 0; side < sides; side++) {
		id = "side" + side;
		selectImage(traceID);
		run("Select None");
		loaded = loadCSVLine(original, id);
		if (!loaded) run("Line Width...", "line=" + lineWidth);
		setTool("polyline");
		message = "Image " + (im + 1) + " of " + total + ": " + sample + "\nLine " + (side + 1) + " of " + sides + "\n \n";
		if (loaded) message = message + "The saved line is shown and can be edited.\n";
		message = message + "Draw the line from posterior to anterior, then click OK.\nThe line applies to every channel and Z-plane.";
		beep();
		waitForUser(message);
		selectImage(traceID);
		while (selectionType() != 5 && selectionType() != 6) {
			waitForUser("No line was found on " + sample + ".\nDraw a line from posterior to anterior, then click OK.");
			selectImage(traceID);
		}
		getSelectionCoordinates(xs, ys);
		width = Roi.getStrokeWidth;
		if (width < 1) width = lineWidth;
		for (vertex = 0; vertex < xs.length; vertex++) {
			result = result + id + "," + vertex + "," + d2s(xs[vertex], 9) + "," + d2s(ys[vertex], 9) + "," + d2s(width, 9) + "\n";
		}
	}
	selectImage(traceID);
	close();

	// Stop if the CSV was changed by another program while this image was open, so
	// that those changes are kept.
	current = "";
	if (File.exists(csvPath)) current = File.openAsString(csvPath);
	if (current != original) exit("The ROIs of " + sample + " were changed elsewhere while the image was open.\nNothing was saved for it: " + csvPath);
	// The CSV is written to a temporary file and then renamed, so the file on disk
	// holds either the old or the new ROIs even if saving is interrupted.
	File.makeDirectory(folder);
	temp = folder + "rois.csv.tmp";
	File.saveString(result, temp);
	if (File.exists(csvPath)) {
		backup = folder + "rois.csv.previous";
		if (File.exists(backup)) deleted = File.delete(backup);
		if (!File.rename(csvPath, backup)) exit("Could not keep the previous CSV of " + sample + ": " + csvPath);
	}
	if (!File.rename(temp, csvPath)) exit("Could not save the CSV of " + sample + ". The previous version is in rois.csv.previous");
}
showMessage("ROI selection", "Finished. The ROI CSVs are saved under:\n" + output + "rois/" + replicate);

// Draw the saved line with this id from a rois.csv text as the current selection.
// Returns false when the CSV has no such line.
function loadCSVLine(text, id) {
	if (text == "") return false;
	lines = split(text, "\n");
	if (replace(lines[0], "\r", "") != "roi_id,vertex_index,x_px,y_px,width_px") exit("The ROI CSV does not have the expected header.");
	count = 0;
	for (j = 1; j < lines.length; j++) {
		row = split(replace(lines[j], "\r", ""), ",");
		if (row.length >= 5 && row[0] == id) count++;
	}
	if (count < 2) return false;
	x = newArray(count);
	y = newArray(count);
	lineWidthSaved = 0;
	found = 0;
	for (j = 1; j < lines.length; j++) {
		row = split(replace(lines[j], "\r", ""), ",");
		if (row.length >= 5 && row[0] == id) {
			index = parseInt(row[1]);
			if (index != found || index >= count) exit("The vertices of " + id + " should be numbered 0, 1, 2, ... in order.");
			x[index] = parseFloat(row[2]);
			y[index] = parseFloat(row[3]);
			lineWidthSaved = parseFloat(row[4]);
			found++;
		}
	}
	makeSelection("polyline", x, y);
	Roi.setStrokeWidth(lineWidthSaved);
	return true;
}
