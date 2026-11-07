// Draw or edit the measurement bands of one image in Fiji.
// The bands are saved to output_root/rois/<replicate>/<sample>/rois.csv, the same
// file the notebook and the napari editor use, so choose the notebook's output_root.
// Each CSV row is one clicked point: roi_id,vertex_index,x_px,y_px,width_px, in
// image pixels counted from 0, starting at the posterior end of the side.
// A CSV that already exists is loaded so its points can be edited. The image
// itself is never saved.

requires("1.54");
input = getDirectory("Folder containing the original exported TIFF images");
output = getDirectory("Python output_root (contains or will contain rois/)");
replicate = "N_1";
parts = split(replace(input,"\\","/"),"/");
for (p=0;p<parts.length;p++) {
    if (matches(parts[p],"N_[0-9]+.*")) {
        fields=split(parts[p],"_"); replicate="N_"+fields[1];
    }
}
Dialog.create("IHC ROI settings");
Dialog.addString("Replicate ID",replicate);
Dialog.addNumber("Width for NEW sides (pixels)",100);
Dialog.addCheckbox("Edit images with existing CSV ROIs",false);
Dialog.addCheckbox("Display matched NIS TIFF when available",true);
Dialog.show();
replicate=Dialog.getString();newWidth=Dialog.getNumber();redo=Dialog.getCheckbox();showNIS=Dialog.getCheckbox();
if (!matches(replicate,"N_[0-9]+")) exit("Replicate must look like N_1 or N_2.");
if (newWidth<1) exit("Width must be positive.");
File.makeDirectory(output+"rois/");File.makeDirectory(output+"rois/"+replicate+"/");
files=getFileList(input);Array.sort(files);
for (im=0;im<files.length;im++) {
    name=files[im];lower=toLowerCase(name);
    if (!(endsWith(lower,".tiff") || endsWith(lower,".tif"))) continue;
    if (indexOf(lower,"denoised")>=0 || name=="Mask.tif" || name=="Cyto_Mask.tif") continue;
    dot=lastIndexOf(name,".");sample=substring(name,0,dot);
    folder=output+"rois/"+replicate+"/"+sample+"/";csvPath=folder+"rois.csv";
    if (File.exists(csvPath) && !redo) continue;
    original="";if (File.exists(csvPath)) original=File.openAsString(csvPath);
    displayPath=input+name;
    if (showNIS) {
        pos=lastIndexOf(sample,"_P");
        if (pos>=0) paired=substring(sample,0,pos)+" - Denoised"+substring(sample,pos)+substring(name,dot);
        else paired=sample+" - Denoised"+substring(name,dot);
        if (File.exists(input+paired)) displayPath=input+paired;
    }
    open(displayPath);imageID=getImageID();
    Stack.getDimensions(w,h,channels,slices,frames);
    if (channels>1) Stack.setDisplayMode("color"); // show one channel at a time
    Stack.setChannel(1);Stack.setSlice(floor(slices/2)+1);
    run("Enhance Contrast", "saturated=0.35"); // changes the display, not the pixel values
    Dialog.create(sample+": number of sides");
    Dialog.addChoice("Mounting",newArray("Flat: two sides","Lateral: one side","Skip"));
    Dialog.show();mount=Dialog.getChoice();
    sides=2;if (mount=="Lateral: one side") sides=1;if (mount=="Skip") sides=0;
    if (sides==0) {close();continue;}
    result="roi_id,vertex_index,x_px,y_px,width_px\n";
    for (side=0;side<sides;side++) {
        selectImage(imageID);run("Select None");id="side"+side;
        loaded=loadCSVLine(original,id);
        if (!loaded) {run("Line Width...","line="+newWidth);setTool("polyline");}
        else setTool("polyline");
        waitForUser(sample+" — "+id,
            "Draw or edit ONE polyline, posterior first.\n"+
            "Existing sides retain their saved width. Use Edit > Options > Line Width to change it.\n"+
            "Change channel or Z freely; the line applies to the whole stack.\n"+
            "Click OK when finished.");
        selectImage(imageID);
        if (selectionType()!=6 && selectionType()!=5) exit("Expected a polyline or straight line. No new CSV saved for "+sample);
        getSelectionCoordinates(xs,ys);width=Roi.getStrokeWidth;
        if (width<1) width=newWidth;
        for (vertex=0;vertex<xs.length;vertex++) {
            result=result+id+","+vertex+","+d2s(xs[vertex],9)+","+d2s(ys[vertex],9)+","+d2s(width,9)+"\n";
        }
    }
    // Stop if the CSV changed since it was loaded, so edits made elsewhere are not lost.
    current="";if (File.exists(csvPath)) current=File.openAsString(csvPath);
    if (current!=original) exit("ROIs changed in another application. Reload them before saving: "+csvPath);
    File.makeDirectory(folder);
    // Write to a temporary file and then rename it, so a half-written CSV is never left behind.
    temp=folder+"rois.csv.tmp";File.saveString(result,temp);
    if (File.exists(csvPath)) {
        backup=folder+"rois.csv.previous";
        if (File.exists(backup)) File.delete(backup);
        if (!File.rename(csvPath,backup)) exit("Could not back up existing CSV: "+csvPath);
    }
    if (!File.rename(temp,csvPath)) exit("Could not publish CSV; previous version remains in rois.csv.previous");
    close();
}
showMessage("ROI selection","Finished. Shared CSV ROIs are under:\n"+output+"rois/"+replicate);

function loadCSVLine(text,id) {
    if (text=="") return false;
    lines=split(text,"\n");count=0;
    if (replace(lines[0],"\r","")!="roi_id,vertex_index,x_px,y_px,width_px") exit("Unexpected ROI CSV header");
    for (j=1;j<lines.length;j++) {
        row=split(replace(lines[j],"\r",""),",");
        if (row.length>=5 && row[0]==id) count++;
    }
    if (count<2) return false;
    x=newArray(count);y=newArray(count);width=0;found=0;
    for (j=1;j<lines.length;j++) {
        row=split(replace(lines[j],"\r",""),",");
        if (row.length>=5 && row[0]==id) {
            index=parseInt(row[1]);
            if (index!=found || index>=count) exit("ROI vertices must be ordered 0,1,2,...");
            x[index]=parseFloat(row[2]);y[index]=parseFloat(row[3]);
            width=parseFloat(row[4]);found++;
        }
    }
    makeSelection("polyline",x,y);Roi.setStrokeWidth(width);
    return true;
}
