// correct_and_inspect_rois.ijm opens chosen samples in Fiji with their saved PSM
// lines drawn on them, so that the lines can be checked and corrected.
//
// The samples are named as in IHC_ROI_only.ijm: <file>_P<position> for ND2 files and
// the file name without .tiff for TIFFs. They are picked from a checklist of every
// sample in the folder.
//
// The lines are read from and saved to rois/<replicate>/<sample>/rois.csv in the
// notebook's output folder, in the same format as IHC_ROI_only.ijm. A sample whose
// lines are accepted without changes is not saved again, so its rois.csv.previous is
// kept. A sample without saved lines can be drawn from scratch.
//
// Choosing 'Skip' in the mounting dialog leaves a sample as it was. Skipping a side
// that has a saved line removes that line when the sample is saved.

requires("1.54");
input = getDirectory("Folder with the images (for ND2 files, the replicate folder)");
output = defaultOutputRoot(input);

// The replicate is the N_<number> part of a folder name in the input path.
replicate = "N_1";
parts = split(replace(input, "\\", "/"), "/");
for (p = 0; p < parts.length; p++) {
	if (matches(parts[p], "N_[0-9]+.*")) replicate = replace(parts[p], "^(N_[0-9]+).*", "$1");
}

Dialog.create("Correct and inspect ROIs");
Dialog.addChoice("Image files:", newArray("ND2", "TIFF"));
Dialog.addString("Replicate:", replicate);
Dialog.addNumber("Line width for new lines (pixels):", 200);
Dialog.addCheckbox("Use the polling loop (space accepts, Alt skips)", true);
Dialog.show();
type = Dialog.getChoice();
replicate = Dialog.getString();
lineWidth = Dialog.getNumber();
use_polling_loop = Dialog.getCheckbox();
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
if (samples.length == 0) exit("No " + type + " images were found in\n" + input);
roiRoot = output + "rois/" + replicate + "/";

// Choose the samples from a checklist.
chosen = newArray(samples.length);
labels = newArray(samples.length);
for (s = 0; s < samples.length; s++) {
	labels[s] = samples[s];
	if (!File.exists(roiRoot + samples[s] + "/rois.csv")) labels[s] = labels[s] + " (no ROIs)";
}
columns = 3;
rows = Math.ceil(samples.length / columns);
Dialog.create("Choose the samples");
Dialog.addMessage("Tick the samples to open.");
Dialog.addCheckboxGroup(rows, columns, labels, newArray(samples.length));
Dialog.show();
for (s = 0; s < samples.length; s++) chosen[s] = Dialog.getCheckbox();
order = newArray(0);
for (s = 0; s < samples.length; s++) {
	if (chosen[s]) order = Array.concat(order, s);
}
total = order.length;
if (total == 0) exit("No samples were chosen.");

File.makeDirectory(output);
File.makeDirectory(output + "rois/");
File.makeDirectory(roiRoot);
luts = newArray("Blue", "Green", "Red", "Magenta"); // display colours for channels 1 to 4
saved = 0;
unchanged = 0;
skipped = 0;

for (k = 0; k < total; k++) {
	im = order[k];
	sample = samples[im];
	folder = roiRoot + sample + "/";
	csvPath = folder + "rois.csv";
	original = "";
	if (File.exists(csvPath)) original = File.openAsString(csvPath);

	if (series[im] > 0) {
		run("Bio-Formats Importer", "open=[" + paths[im] + "] color_mode=Colorized view=Hyperstack stack_order=XYCZT series_" + series[im]);
	} else {
		open(paths[im]);
	}
	traceID = getImageID();
	Stack.getDimensions(imageWidth, imageHeight, channels, slices, frames);
	if (channels > 1) Stack.setDisplayMode("color"); // one channel at a time, scrolled with the C slider
	Stack.setSlice(floor(slices / 2) + 1);
	for (ch = 1; ch <= channels; ch++) {
		Stack.setChannel(ch);
		run(luts[(ch - 1) % 4]);
		run("Enhance Contrast", "saturated=0.35"); // changes the display, not the pixel values
	}
	Stack.setChannel(1);

	// Show the saved lines together before asking, so the sample can be inspected
	// as a whole. A saved side1 means the sample was drawn as flat.
	showSavedLines(original, traceID);
	mountDefault = "Lateral";
	if (original != "" && csvHasSide(original, "side1")) mountDefault = "Flat";
	status = "no saved ROIs";
	if (original != "") status = "saved ROIs are shown";
	Dialog.createNonBlocking("How is this sample mounted?");
	Dialog.addMessage("Sample " + (k + 1) + " of " + total + ": " + sample + " (" + status + ")\n'Skip' leaves this sample as it was.");
	Dialog.addChoice("Mounting:", newArray("Flat", "Lateral", "Skip"), mountDefault);
	beep();
	Dialog.show();
	mount = Dialog.getChoice();
	selectImage(traceID);
	Overlay.remove;
	sides = 0;
	if (mount == "Flat") sides = 2;
	if (mount == "Lateral") sides = 1;
	if (sides == 0) {
		selectImage(traceID);
		close();
		skipped++;
		continue;
	}

	result = "roi_id,vertex_index,x_px,y_px,width_px\n";
	before = "roi_id,vertex_index,x_px,y_px,width_px\n"; // the saved lines in the same format, to detect changes
	drawn = 0;
	for (side = 0; side < sides; side++) {
		id = "side" + side;
		selectImage(traceID);
		run("Select None");
		// Set the drawing width before restoring a line. This also gives restored
		// lines Fiji's translucent wide-line display and sets the width of a redraw.
		run("Line Width...", "line=" + lineWidth);
		loaded = loadCSVLine(original, id);
		if (loaded) before = before + formatLine(id, lineWidth);
		setTool("polyline");
		message = "Sample " + (k + 1) + " of " + total + ": " + sample + "\nLine " + (side + 1) + " of " + sides + "\n \n";
		if (loaded) message = message + "The saved line is shown and can be edited.\n";
		else message = message + "This side has no saved line.\n";
		beep();
		if (use_polling_loop) {
			message = message + "Draw or correct the line from posterior to anterior.\nThe line applies to every channel and Z-plane.\n \nPress space to accept or press Alt to skip this side.";
			skip = pollSide(message, traceID, sample);
		} else {
			message = message + "Draw or correct the line from posterior to anterior, then click OK.\nThe line applies to every channel and Z-plane.\n \nTick 'Skip this side' to leave this side out.";
			skip = askSide(message);
			selectImage(traceID);
			while (!skip && selectionType() != 5 && selectionType() != 6) {
				skip = askSide("No line was found on " + sample + ".\nDraw a line from posterior to anterior, then click OK,\nor tick 'Skip this side'.");
				selectImage(traceID);
			}
		}
		selectImage(traceID);
		if (skip) continue;
		drawn++;
		result = result + formatLine(id, lineWidth);
	}
	// The saved lines of sides that were not opened (side1 of a sample now marked
	// lateral) count as part of the saved state, so dropping them is a change.
	for (side = sides; side < 2; side++) {
		if (original != "" && csvHasSide(original, "side" + side)) before = before + "dropped side" + side + "\n";
	}
	selectImage(traceID);
	close();
	// With every side skipped there is nothing to save, and a saved CSV is left as it was.
	if (drawn == 0) {
		skipped++;
		continue;
	}
	if (result == before) {
		unchanged++;
		continue;
	}

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
	saved++;
}
showMessage("Correct and inspect ROIs", "Finished.\nSaved: " + saved + "\nUnchanged: " + unchanged + "\nSkipped: " + skipped + "\n \nThe ROI CSVs are under:\n" + roiRoot);

// Write the current line selection as CSV rows with this id.
function formatLine(id, fallbackWidth) {
	getSelectionCoordinates(xs, ys);
	width = Roi.getStrokeWidth;
	if (width < 1) width = fallbackWidth;
	text = "";
	for (vertex = 0; vertex < xs.length; vertex++) {
		text = text + id + "," + vertex + "," + d2s(xs[vertex], 9) + "," + d2s(ys[vertex], 9) + "," + d2s(width, 9) + "\n";
	}
	return text;
}

// Add every saved line to the overlay of the image so that all sides are seen at once.
function showSavedLines(text, imageID) {
	if (text == "") return;
	selectImage(imageID);
	Overlay.remove;
	for (side = 0; side < 2; side++) {
		if (loadCSVLine(text, "side" + side)) {
			Roi.setStrokeColor("yellow");
			Overlay.addSelection;
		}
	}
	run("Select None");
}

// Return true when the rois.csv text has rows for this id.
function csvHasSide(text, id) {
	lines = split(text, "\n");
	for (j = 1; j < lines.length; j++) {
		row = split(replace(lines[j], "\r", ""), ",");
		if (row.length >= 5 && row[0] == id) return true;
	}
	return false;
}

// Show the message in a dialog that leaves the image usable for drawing, and
// return true when 'Skip this side' was ticked.
function askSide(message) {
	Dialog.createNonBlocking("Draw the line");
	Dialog.addMessage(message);
	Dialog.addCheckbox("Skip this side", false);
	Dialog.show();
	return Dialog.getCheckbox();
}

// Show the message in the Log window beside the image and keep the length of the
// line up to date on the last line of the Log. The image stays usable for drawing.
// Returns false when space is pressed with a line drawn and true when Alt is pressed.
function pollSide(message, imageID, sample) {
	selectImage(imageID);
	getPixelSize(unit, pixelWidth, pixelHeight);
	getLocationAndSize(windowX, windowY, windowWidth, windowHeight);
	print("\\Clear");
	print(message);
	if (unit == "pixel" || unit == "pixels") print("This image has no pixel size. The length is given in pixels.");
	print(" ");
	print("No line has been drawn yet.");
	// The Log window is placed to the right of the image window and kept on the
	// screen. The image is then made active again because only the active window
	// receives the key presses.
	selectWindow("Log");
	setLocation(Math.min(windowX + windowWidth, screenWidth - 400), windowY);
	selectImage(imageID);
	shown = "No line has been drawn yet.";
	while (true) {
		if (!isOpen(imageID)) exit("The image of " + sample + " was closed while its line was being drawn.\nNothing was saved for it.");
		// Selection types 5 and 6 are a straight line and a segmented line.
		isLine = selectionType() == 5 || selectionType() == 6;
		if (isKeyDown("alt")) {
			skip = true;
			break;
		}
		if (isKeyDown("space") && isLine) {
			skip = false;
			break;
		}
		text = "No line has been drawn yet.";
		// The length is measured along the drawn line in the calibrated unit of the image.
		if (isLine) text = "The line is " + d2s(getValue("Length"), 1) + " " + unit + " long.";
		// The Log is rewritten only when the text changes. "\\Update:" replaces its last line.
		if (text != shown) {
			print("\\Update:" + text);
			shown = text;
		}
	}
	// Wait until the key is released so that one press does not also end the next side.
	while (isKeyDown("space") || isKeyDown("alt")) wait(50);
	return skip;
}

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

// Find the notebook's default output beside the experimental-data folder.
// Walking up also handles TIFFs selected from a subfolder of a replicate.
function defaultOutputRoot(imageFolder) {
	folder = replace(imageFolder, "\\", "/");
	while (endsWith(folder, "/")) folder = substring(folder, 0, lengthOf(folder) - 1);
	while (folder != "") {
		if (File.exists(folder + "/ihc_pipeline.py")) return folder + "/ihc_results/";
		analysis = folder + "/analysis_scripts";
		if (File.exists(analysis + "/ihc_pipeline.py")) return analysis + "/ihc_results/";
		parent = File.getParent(folder);
		if (parent == folder) break;
		folder = parent;
	}
	exit("Could not find analysis_scripts/ihc_pipeline.py above the image folder:\n" + imageFolder);
}
