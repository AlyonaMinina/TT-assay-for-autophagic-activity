//This macro processes a folder with .czi images to quantify intensity of fluorescence (in arbitrary units) within user-selected ROIs. For each image, ROIs are saved as.zip files and quantifications are saved as .csv file

//Step by step:
	//1. Create a folder with .czi images selected for analysis (do not export them as other file formats!)
	//2. Drag and drop macro file into ImageJ to have access to the code
	//3. NB! make sure that at least one chennel is called "GFP". If needed edit the default number of ROIs and their dimensions in the dialog window
	//4. The macro will open one image at a time and wait for the user to adjust ROI position size.
	//5. After all ROIs are adjusted  -> click ok. The macro will process all ROIs present in the ROI Manager. The ROI.zip file will be saved for each image file individually, while quantification data will be compiled into a single file.csv containing information about image name, ROI number, ROI area in um2, fluorescence intensity in each channel (Integrated Density), mean intensity for the second channel (IntDen/area) and ratio of fluorescence intensity Ch2/Ch1
   //NB! You can rerun the macro on the same folder, in this case the script will automatically load the saved ROIs for corresponding image and offer you to adjust them. To enable this do not change content or names of the folders and subfolders you analyzed



//Clear the log window if it was open
if (isOpen("Log")){
	selectWindow("Log");
	run("Close");
}
	
//Print the unnecessary greeting
print(" ");
print("Welcome to the macro for TT assay!");
print(" ");
print("Please select the folder with images for analysis");
print(" ");

//Find the original directory and create a new one for quantification results
original_dir = getDirectory("Select a directory");
original_folder_name = File.getName(original_dir);
output_dir = original_dir +"Results" + File.separator;
File.makeDirectory(output_dir);

// Get a list of all the files in the directory
file_list = getFileList(original_dir);

//Create a shorter list containing .czi files only
czi_list = newArray(0);
for(z = 0; z < file_list.length; z++) {
	if(endsWith(file_list[z], ".czi")) {
		czi_list = Array.concat(czi_list, file_list[z]);
		}
	}
	
//inform the user about how many images will be analyzed in the selected folder
print(czi_list.length + " images were detected for analysis");
print("");

//Obtain coordinates to draw ROIs in the center of the image
		

/// Request info from the user about the number and dimensions of the ROIs they wish to analyze
Channel_1 = "GFP";	  
Channel_2 = "RFP";
number_of_ROIs = 1;
ROI_height = 159;
ROI_width = 159;

Dialog.create("Please provide ROIs parameters for your images");
Dialog.addString("Channel 1 name:", Channel_1);
Dialog.addToSameRow();
Dialog.addString("Channel 2 name:", Channel_2);
Dialog.addNumber("Number of ROIs to be analyzed on each image:", number_of_ROIs);
Dialog.addNumber("Dimensions of ROIs. ROI height in um:", ROI_height);
Dialog.addNumber("ROI width in um:", ROI_width);
Dialog.show();

Channel_1 = Dialog.getString();
Channel_2 = Dialog.getString();

number_of_ROIs = Dialog.getNumber();
ROI_height = Dialog.getNumber();
ROI_width = Dialog.getNumber();	

	
//Create the table for all results
Table.create("Image Results");
	
//Loop analysis through the list of . czi files
for (i = 0; i < czi_list.length; i++){
	path = original_dir + czi_list[i];
	run("Bio-Formats Windowless Importer",  "open=path");
		      
//Get the image file title and remove the extension from it    
title = getTitle();
a = lengthOf(title);
b = a-4;
short_name = substring(title, 0, b);
		
//Print for the user what image is being processed
print ("Processing image " + i+1 + " out of " + czi_list.length + ":");
print(title);
print("");
					
//Adjust the ROIs for each micrograph
run("ROI Manager...");
		
// Make sure ROI Manager is clean of any previously used ROIs
roiManager("reset");
	
// Obtain coordinates to draw ROIs in the center of the image
x = getWidth()/2;
toScaled(x);
//x_coordinate =  parseInt(x);
x_coordinate =  parseInt(0);

y = getHeight()/2;
toScaled(y);
//y_coordinate =  parseInt(y);
y_coordinate =  parseInt(0);
	
//Draw ROIs of the user-provided number and dimensions. Automatically load in already existing ROIs for the image (if not desired, comment out lines 107-114 and the line 124.
ROIset = output_dir + short_name + "_ROIs.zip";
f = File.exists(ROIset);
	if(f>0){ 
	roiManager("Open", ROIset);
	roiManager("Show All");
	roiManager("Show All with labels");
	}
		
		else {
			for (no_roi = 0; no_roi < number_of_ROIs; no_roi++) {
				makeRectangle(x_coordinate, y_coordinate, ROI_height, ROI_width);
			    run("Specify...", "width=ROI_width height=ROI_height x=x_coordinate y=y_coordinate slice=1 scaled");
		        roiManager("Add");
			    roiManager("Select", no_roi);
		        roiManager("Rename", no_roi + 1);
		        roiManager("Show All");
				roiManager("Show All with labels");
				}
			}
//Wait for the user to adjust the ROIs size and position
waitForUser("Adjust each ROI, then hit OK"); 

//Relable ROIs again in case Ayona wanted to delete redraw all of them
for (no_roi = 0; no_roi < number_of_ROIs; no_roi++) {
		roiManager("Select", no_roi);
        roiManager("Rename", no_roi + 1);
        roiManager("Show All");
		roiManager("Show All with labels");
	}
		
						
//Perform Fluorescence intensity for each ROI and save the results into a custom table
run("ROI Manager...");
ROI_number = roiManager("count");
for ( r=0; r<ROI_number; r++ ) {
	roiManager("Select", r);
	current_last_row = Table.size("Image Results");
	Table.set("File name", current_last_row, short_name, "Image Results");
	Table.set("ROI number", current_last_row, r+1, "Image Results");
	
	//just in case if Results were open from anoter analysis
	run("Clear Results");
	run("Set Measurements...", "area redirect=None decimal=3");
	run("Measure");
	area = getResult("Area", 0);
	Table.set("Area in um2", current_last_row, area, "Image Results");
	run("Clear Results");
	
//Quantify fluorescence on the first channel of the image.
	setSlice(1);
	run("Set Measurements...", "integrated redirect=None decimal=3");
	run("Measure");
	IntDenCh1 = getResult("IntDen",  0);
	Column_1 = "IntDen " + Channel_1;
	Table.set(Column_1, current_last_row, IntDenCh1, "Image Results");
	RawIntDenCh1 = getResult("RawIntDen",  0);
	run("Clear Results");
			
			
//Quantify fluorescence on the second channel of the image.
	setSlice(2);
	run("Set Measurements...", "integrated redirect=None decimal=3");
	run("Measure");
	IntDenCh2 = getResult("IntDen",  0);
	Column_2 = "IntDen " + Channel_2;
	Table.set(Column_2, current_last_row, IntDenCh2, "Image Results");
	run("Clear Results");
			
			
//Calculate the ratio of fluorescence Ch2/Ch1 based on the IntDen values
	Ch1 = Table.get(Column_1, current_last_row,"Image Results");
	//ParseInt is needed, because Mac OS retrieves values from the table as strings
	Ch1_Int = parseInt(Ch1);
	Ch2 = Table.get(Column_2, current_last_row,"Image Results");
	Ch2_Int = parseInt(Ch2);

//Detect which of the channels contains GFP signal to determine how to calculate the ratio
	if(Channel_1 == "GFP") {
		ratio = Ch2_Int/Ch1_Int;
	Column_3 = "Fluorescence ratio of "+ Channel_2 + " to " + Channel_1;
	Table.set(Column_3, current_last_row, ratio, "Image Results");	
		}
	else {
		ratio = Ch1_Int/Ch2_Int;
	Column_3 = "Fluorescence ratio of "+ Channel_1 + " to " + Channel_2;
	Table.set(Column_3, current_last_row, ratio, "Image Results");	
		}
	}

//Save maxima quantification as .csv file and ROIs as a .zip file
	roiManager("Save", output_dir + short_name +"_ROIs.zip");
	run("Close All");
	roiManager("reset");
	run("Clear Results");
	}		

//Save the quantification results into a .csv table file
Table.save(output_dir + "Fluorescence quantification for " + original_folder_name + ".csv");
 
//A feeble attempt to close those pesky ImageJ windows		
run("Close All");
roiManager("reset");
run("Clear Results");
selectWindow("Image Results");
run("Close");
selectWindow("Results");
run("Close");
selectWindow("ROI Manager");
run("Close");

//Print the final message
print(" ");
print("All Done!");
print("Your quantification results are saved in the folder " + output_dir);
print(" "); 
print(" ");
print("Alyona Minina. 2024");
