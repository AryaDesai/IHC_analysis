// PSM_phase_sorting.ijm - manual clock-phase sorting.
//
// Shows every position of every .nd2 file as an average intensity projection of all z slices, starting on the
// clock channel (use the C slider to see the other channels), and asks for its phase.
// Next to the image, a plot shows the clock-channel profile along the PSM lines drawn with IHC_ROI_only.ijm
// (read from <ihc_results>/rois/<replicate>/<file>_P<position>/rois.csv). The profile is the mean across the
// line width, measured on the projection. Distances are counted from the anterior end by default, as with
// align_profiles = 1 in ihc_pipeline.py: only the distances are renumbered, the values are not changed.
// Only the mean of the sides is drawn: at each distance, the mean of whichever sides reach that distance.
// The dialog can also denoise the image (Remove Outliers + Despeckle) before you decide; the profile is
// then measured again on the denoised image.
// Answers are saved to phases.csv in the image folder, with two columns: file, phase (1, 2, 3 or unclear).
// Positions already in phases.csv are skipped, so you can stop (Cancel) and continue later.
// Denoising and contrast only change the display; the .nd2 files and ROI files are not modified.

var clockChannel = 2;       // channel shown first (2 = 488 clock)
var profileChannel = 2;     // channel measured for the profile plot
var alignTo = "anterior";   // end the profile distances are measured from: "anterior" or "posterior"
var saturated = 3;          // % of pixels shown fully saturated; higher = brighter display
var outlierRadius = 2;      // Remove Outliers: radius (pixels)
var outlierThreshold = 50;  // Remove Outliers: how much brighter than its surroundings a pixel must be to count as an outlier
var plotWidth = 800;        // plot frame size in screen pixels; larger = more detail
var plotHeight = 450;
var plotFontSize = 16;
var plotLineWidth = 2;

if (alignTo != "anterior" && alignTo != "posterior") exit("alignTo should be \"anterior\" or \"posterior\".");

dir = getDirectory("Choose the folder with the .nd2 files");

// Read ROI lines from the same default output folder as IHC_ROI_only.ijm.
roiRoot = defaultOutputRoot(dir);
replicate = "N_1";
parts = split(replace(dir, "\\", "/"), "/");
for (p = 0; p < parts.length; p++) {
	if (matches(parts[p], "N_[0-9]+.*")) replicate = replace(parts[p], "^(N_[0-9]+).*", "$1");
}
Dialog.create("Phase sorting settings");
Dialog.addString("Replicate:", replicate);
Dialog.addMessage("Profiles are aligned to the " + alignTo + " end (set alignTo at the top of the macro).");
Dialog.show();
replicate = Dialog.getString();

// phases.csv: create it, or read what is already sorted
csvPath = dir + "phases.csv";
if (!File.exists(csvPath)) File.saveString("file,phase\n", csvPath);
sorted = File.openAsString(csvPath);

run("Close All");
setBatchMode(true);   // each image stays hidden until it is ready, so it only appears once
run("Bio-Formats Macro Extensions");
list = getFileList(dir);
Array.sort(list);
for (i = 0; i < list.length; i++) {
	if (!endsWith(toLowerCase(list[i]), ".nd2")) continue;

	// An .nd2 file can hold several positions
	Ext.setId(dir + list[i]);
	Ext.getSeriesCount(nSeries);
	Ext.close();

	for (s = 1; s <= nSeries; s++) {
		name = list[i];
		if (nSeries > 1) name = list[i] + " P" + s;
		if (indexOf(sorted, "\n" + name + ",") >= 0) continue;   // already sorted
		// Sample name used by IHC_ROI_only.ijm and the notebook. A "<file> - Denoised.nd2" export
		// uses the ROIs drawn on "<file>.nd2", since the ROI macro skips denoised files.
		stem = substring(list[i], 0, lengthOf(list[i]) - 4);
		stem = replace(stem, "(?i)\\s*-\\s*denoised$", "");
		roiCsv = roiRoot + "rois" + File.separator + replicate + File.separator + stem + "_P" + s + File.separator + "rois.csv";

		// Average intensity projection of all z slices, shown one channel at a time
		run("Bio-Formats Importer", "open=[" + dir + list[i] + "] color_mode=Colorized view=Hyperstack stack_order=XYCZT use_virtual_stack series_" + s);
		stackID = getImageID();
		run("Z Project...", "projection=[Average Intensity]");
		rename(name);
		selectImage(stackID);
		close();
		selectWindow(name);
		if (nChannels() > 1) Stack.setDisplayMode("color");
		brighten(name);
		setBatchMode("show");
		showProfile(name, roiCsv);

		// Ask for the phase; ticking "denoise" cleans the image and asks again
		phase = "";
		while (phase == "") {
			Dialog.createNonBlocking("Phase: " + name);
			Dialog.addMessage(name + "\n \n" +
				"Phase 1: large posterior clock domain + separate stripe in the anterior PSM\n" +
				"Phase 2: posterior domain reaches about mid-PSM + thin band at the anterior end\n" +
				"Phase 3: only a thin band at the tail-bud tip + broad stripe in mid-to-anterior PSM\n \n" +
				"Cancel stops sorting; answers so far are saved.");
			Dialog.addRadioButtonGroup("Phase", newArray("1", "2", "3", "unclear"), 1, 4, "1");
			Dialog.addCheckbox("Denoise this image first (Remove Outliers + Despeckle), then ask again", false);
			Dialog.show();
			choice = Dialog.getRadioButton();
			if (Dialog.getCheckbox()) {
				denoise(name);
				showProfile(name, roiCsv);
			}
			else phase = choice;
		}
		File.append(name + "," + phase, csvPath);

		closeWindow("Profile: " + name);
		closeWindow(name);
	}
}
setBatchMode(false);
showMessage("Phase sorting", "Done. Phases saved in\n" + csvPath);


// Contrast-enhances every channel, then shows the clock channel.
function brighten(title) {
	selectWindow(title);
	for (c = 1; c <= nChannels(); c++) {
		Stack.setChannel(c);
		run("Enhance Contrast", "saturated=" + saturated);
	}
	Stack.setChannel(minOf(clockChannel, nChannels()));
}

// Removes bright single-pixel noise and speckle from all channels, then re-brightens.
function denoise(title) {
	selectWindow(title);
	run("Remove Outliers...", "radius=" + outlierRadius + " threshold=" + outlierThreshold + " which=Bright stack");
	run("Despeckle", "stack");
	brighten(title);
	updateDisplay();
}

function nChannels() {
	getDimensions(w, h, nc, nz, nt);
	return nc;
}

function closeWindow(title) {
	if (isOpen(title)) {
		selectWindow(title);
		close();
	}
}

// Plots the profile of profileChannel along each ROI line of the image, next to the image.
// Profiles run posterior to anterior, the order the lines were drawn in. With alignTo = "anterior"
// the point at index k of a side with n points is plotted at distance (n - 1 - k) * pixel size,
// the same renumbering as ihc_pipeline.py; the values themselves are not changed.
function showProfile(title, roiPath) {
	plotTitle = "Profile: " + title;
	closeWindow(plotTitle);
	selectWindow(title);
	getPixelSize(unit, pw, ph);
	Stack.getPosition(shownChannel, sl, fr);
	channel = minOf(profileChannel, nChannels());
	Stack.setChannel(channel);

	// Read the lines from rois.csv (roi_id, vertex_index, x_px, y_px, width_px)
	ids = newArray(0);
	values = newArray(0);    // every side's profile, one after the other
	lengths = newArray(0);
	note = "";
	if (roiRoot == "" || !File.exists(roiPath)) note = "No ROIs found for this position";
	else {
		rows = split(File.openAsString(roiPath), "\n");
		if (replace(rows[0], "\r", "") != "roi_id,vertex_index,x_px,y_px,width_px") note = "Unexpected header in rois.csv";
		else {
			for (j = 1; j < rows.length; j++) {
				row = split(replace(rows[j], "\r", ""), ",");
				if (row.length >= 5 && !contains(ids, row[0])) ids = Array.concat(ids, row[0]);
			}
			for (k = 0; k < ids.length; k++) {
				count = 0;
				for (j = 1; j < rows.length; j++) {
					row = split(replace(rows[j], "\r", ""), ",");
					if (row.length >= 5 && row[0] == ids[k]) count++;
				}
				xs = newArray(count);
				ys = newArray(count);
				width = 1;
				for (j = 1; j < rows.length; j++) {
					row = split(replace(rows[j], "\r", ""), ",");
					if (row.length >= 5 && row[0] == ids[k]) {
						v = parseInt(row[1]);
						xs[v] = parseFloat(row[2]);
						ys[v] = parseFloat(row[3]);
						width = parseFloat(row[4]);
					}
				}
				makeSelection("polyline", xs, ys);
				Roi.setStrokeWidth(width);
				profile = getProfile();   // mean across the line width, posterior to anterior
				n = profile.length;
				values = Array.concat(values, profile);
				lengths = Array.concat(lengths, n);
			}
			run("Select None");
		}
	}
	Stack.setChannel(shownChannel);
	getLocationAndSize(imageX, imageY, imageW, imageH);

	Plot.create(plotTitle, "Distance from " + alignTo + " end (" + unit + ")", "Mean intensity, channel " + channel);
	Plot.setFrameSize(plotWidth, plotHeight);
	Plot.setFontSize(plotFontSize);
	Plot.setAxisLabelSize(plotFontSize);
	if (lengths.length == 0) {
		Plot.setLimits(0, 1, 0, 1);
		Plot.addText(note, 0.1, 0.5);
	} else {
		longest = 0;
		for (k = 0; k < lengths.length; k++) longest = maxOf(longest, lengths[k]);
		sum = newArray(longest);
		count = newArray(longest);
		start = 0;
		for (k = 0; k < lengths.length; k++) {
			n = lengths[k];
			for (j = 0; j < n; j++) {
				position = j;   // posterior alignment: index from the posterior end
				if (alignTo == "anterior") position = n - 1 - j;
				sum[position] = sum[position] + values[start + j];
				count[position] = count[position] + 1;
			}
			start = start + n;
		}
		// Mean of the sides that reach each distance, as when the pipeline combines sides
		// (with one side this is that side's profile)
		mx = newArray(longest);
		my = newArray(longest);
		for (j = 0; j < longest; j++) {
			mx[j] = j * pw;
			my[j] = sum[j] / count[j];
		}
		Plot.setColor("black");
		Plot.setLineWidth(plotLineWidth);
		Plot.add("line", mx, my);
		label = "1 side";
		if (lengths.length > 1) label = "mean of " + lengths.length + " sides";

	}
	Plot.show();
	setBatchMode("show");
	// Place the plot to the right of the image, or below it when there is no room
	getLocationAndSize(px, py, plotW, plotH);
	if (imageX + imageW + plotW <= screenWidth) setLocation(imageX + imageW, imageY);
	else setLocation(imageX, minOf(imageY + imageH, screenHeight - plotH));
	selectWindow(title);
}

function contains(array, value) {
	for (a = 0; a < array.length; a++) {
		if (array[a] == value) return true;
	}
	return false;
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
