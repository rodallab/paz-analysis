/*
 * To do:
 * Documentation
 * Port to python (long term)
 */


BATCHMODE = true;
DEBUG=false; // Log progress on 'TS_Log' to help catch cryptic bugs
HEATMAP=false; // Whether to generate heatmap

// *** Basic info
// For CHANS, use all caps and standardized names: BRP, CLC, DAP160, NWK, DYN, FASII
CHANS=newArray("NWK", "FASII", "BRP"); // List of channels
AZ=3; //Channel containing active zone marker for analysis PUT 0 if none
// FOR VERSION, USE STANDARDIZED CODE: VERSIONNUMBER_EXPERIMENTLABEL
VERSION="218_AMS020";
IMAGEEXT="tif";

// *** Whole NMJ segmentation parameters
PREBGRB=50; // Rolling ball radius to apply BG subtraction to whole image prior to any seg or analysis. 0=none.
SIGMA=1; //This is the amount to blur your image to make the mask
SEGMENT=0; //This is the channel to use to make a mask. Put 0 to combine channels (see below)
COMBINE=newArray(1,2,3); // Which channels to add together to make mask IF segment==0.
ERODE=3; //This is the number of times to erode the NMJ mask
METHOD="Li";

// *** Mesh segmentation parameters
tolerances=newArray(2560, 3840, 8000);   // Per channel tolerance for minima detection
comp_tolerance = 1280; // Tolerance for composite PAZ image for minima detection
IG_RADIUS = 4; // internal gradient sigma for edge sharpening for single PAZ channel
COMP_IG = 5; // internal gradient sigma for edge sharpening for composite PAZ channel
POST_IG = 1;  // sigma for blur after gradient filter
TRACEWIDTH = 2; 
logr = 2; // sigma to log filer AZ prior to segmentation
AZ_INTENSITY_FILTER = 1;

// *** Assorted analysis parameters
EDM_FILTER_MIN = -1; // EDM filter to prune edge meshes. Put -1 for none.
JACCARDWIDTH = 2;  // with of mesh trace for jaccard index
Z = 2.5; // Std Dev above mean to clip values for presenting AZ-PAZ heatmap

HEADERS=newArray("Entropy");  // Columns to keep for haralick's
NEWHEADERS=newArray("Entr"); // Fix to remove spaces and other weirdness
GLCMdmax=8; // Largest distance to measure haralick's
dSTEP=2; // amount to increment haralick d

// **** RADIAL PROFILE PARAMETERS ****
SIZE=11; // Size to set all radial profile arrays to
minAREA=0.25; // Minimum area size to analyze for radial profile
minCIRC=0.6; // Minimum circularity to analyze for radial profile


temp = split(VERSION, "_");
EXPERIMENT=temp[1];

run("Colors...", "foreground=white background=black selection=yellow");
run("Options...", "iterations=1 count=1 black");
if(roiManager("count")>0)
	roiManager("reset");
	
if(selectionType()>-1)
	run("Select None");
	
close("*");
if(DEBUG==true)
	showText("TS_Log", "Start Log\n");
// Get path to images, assemble image list
// imgdone is an image counter to serve as 
// a row label in case there are non-image
// files in the image folder.

path=getDirectory("Select a folder to analyze ...");
//path="Z:\\Current members\\DelSignore\\PAZ Organization Analysis Project\\PAZ\\SD570a-Clc-BRP-Dyn_Processed-LateMaxIP\\test";
savepath = path+"\\Analysis_V"+VERSION;
centroidpath = savepath+"\\"+EXPERIMENT+"_CentroidData";
if(File.exists(savepath) != true)
	File.makeDirectory(savepath);
if(File.exists(centroidpath) != true)
	File.makeDirectory(centroidpath);
	
count=getImgCount(path);
files=newArray(count);
fl=getFileList(path);
imgdone=0;
setBatchMode(BATCHMODE);

// Table NMJ_Data_AllImages will hold a compiled list for all images of all NMJ level data
// eg: NMJ intensities, coloc, and mesh properties
eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('NMJ_Data_AllImages');");
eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('Comp-ALLTHECELLS');");
eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('Comp-ALLTHEPROFILES');");
eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('Comp-NMJAverages');");
eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('Comp-Profiles');");
eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('AZ-Profiles');");
eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('AZ-ALLTHEPROFILES');");
for(ch=0; ch<CHANS.length; ch++)
	{
	tablename=CHANS[ch]+"-ALLTHECELLS";
	eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('"+tablename+"');");
	totalcellrows=newArray(CHANS.length);		
	tablename=CHANS[ch]+"-NMJAverages";
	eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('"+tablename+"');");
	}

for(i=0; i<fl.length; i++)
	{
	if(endsWith(fl[i], IMAGEEXT)==true)
		{
		if(DEBUG==true)	
			print("[TS_Log]", fl[i]+ " line 94\n");	
		open(path+"\\"+fl[i]);
		getDimensions(width, height, channels, slices, frames);
		getPixelSize(unit, pw, ph, pd);
		SCALE=1/pw;
		
		imgpath = savepath+"\\"+stripX(fl[i]);
		if(File.exists(imgpath) != true)
			File.makeDirectory(imgpath);
		
		// Create image-specific tables
		eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('CentroidData');");
		for(ch=0; ch<channels; ch++)
			{
			totalcellrows[ch]=Table.size(CHANS[ch]+"-ALLTHECELLS");
			Table.set("File", imgdone, fl[i], CHANS[ch]+"-NMJAverages");
			tablename=CHANS[ch]+"_NetData";
			eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('"+tablename+"');");
			}
			
		total_comp_rows=Table.size("Comp-ALLTHECELLS");
		eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('Comp_NetData');");
		Table.set("File", imgdone, fl[i], "Comp-NMJAverages");
	
		// *******************************************************************
		// *****   Segmentation and pre-processing of NMJ and AZ-PAZ   *******
		// *******************************************************************
		
		if(PREBGRB>0)
			run("Subtract Background...", "rolling="+PREBGRB+" stack");

		// Generate tight whole NMJ mask ("mask") and distance map ("EDM")
		// for downstream analysis and segmentations
		segment(fl[i]);

		// Normalize each channel by mean for downstream analysis. Called mean_norm
		meanNorm(fl[i], "mask");

		// Normalize each channel by min-max, clipping high values for downstream seg. Called min_max
		minmaxNorm(fl[i], "mask", 3);
		
		if(DEBUG==true)
			print("[TS_Log]", fl[i]+ " line 139\n");	
		// Create LoG filtered image for additional CoV measurements
		selectImage("mean_norm");
		run("Duplicate...", "title=log_filtered duplicate");
		for(c=1; c<channels+1; c++)
			{
			Stack.setChannel(c);
			run("Mexican Hat Filter", "radius="+logr+" slice");
			}

		// Make single channel mesh images for each individual channel (called CHANS[x]_mesh)
		// Save centroid locations to image-specific table "CentroidData"
		// If a channel is AZ, also generate AZ object segmentation "AZ_objects"
		// params = img, mask, channel, gaussian sigma, tolerance, LoG sigma
		for(c=1; c<CHANS.length+1; c++)
			{
			if(DEBUG==true)
				print("[TS_Log]", fl[i]+ " line 156\n");
			centroids = makeSingleChannelNetwork("min_max", "mask", c, 1, tolerances[c-1], logr);
			if(DEBUG==true)
				print("[TS_Log]", fl[i]+ " line 158\n");
			cx=Array.slice(centroids, 0, centroids.length/2); 
			cy=Array.slice(centroids, centroids.length/2);	
			for(row=0; row<cx.length; row++)
				{	
				Table.set(""+CHANS[c-1]+"_X", row, cx[row], "CentroidData");
				Table.set(""+CHANS[c-1]+"_Y", row, cy[row], "CentroidData");
				}

			// Make mesh mask for each individual channel
			selectImage(""+CHANS[c-1]+"_mesh");
			run("Select All");
			run("Analyze Particles...", "clear add");  //1.53m bug
			//run("Analyze Particles...", "clear add exclude");
			makeMeshMask(fl[i]);
			rename("C"+c+"_MeshMask");
			if(roiManager("count")>0)
				roiManager("reset");
			}
		if(DEBUG==true)
			print("[TS_Log]", fl[i]+ " line 177\n");
		// Make mesh binary map for composite of PAZ and AZ channels
		// Called 'Composite_mesh'
		// Save centroids for composite PAZ
		centroids = makeCompositeNetwork(fl[i], "mask", 1, 1000);
		cx=Array.slice(centroids, 0, centroids.length/2); 
		cy=Array.slice(centroids, centroids.length/2);	
		for(row=0; row<cx.length; row++)
			{	
			Table.set("Comp_X", row, cx[row], "CentroidData");
			Table.set("Comp_Y", row, cy[row], "CentroidData");
			}
		if(DEBUG==true)
			print("[TS_Log]", fl[i]+ " line 191\n");
		// Save centroid data
		Table.save(centroidpath+"\\"+stripX(fl[i])+"_centroids.csv", "CentroidData");
		// Make mesh mask for each individual channel
		selectImage("Composite_mesh");
		run("Select All");
		run("Analyze Particles...", "clear add");  //1.53m bug
		//run("Analyze Particles...", "clear add exclude");
		makeMeshMask(fl[i]);
		rename("Composite_MeshMask");
		roiManager("reset");
		if(DEBUG==true)
			print("[TS_Log]", fl[i]+ " line 202\n");
		
		// *******************************************************************
		// *********************   Analyses of Stuffs    *********************
		// *******************************************************************
		
		// ****************************
		// Analysis of whole NMJ values
		// ****************************
		
		selectImage("mask");
		run("Select All");
		run("Create Selection");
		roiManager("Add");
		if(DEBUG==true)
			print("[TS_Log]", fl[i]+ " line 217\n");
		// Get area and std_dev from mean normalized images
		for(c=0; c<CHANS.length; c++)
			{
			selectImage("mean_norm");
			roiManager("Select", 0);
			Stack.setChannel(c+1);
			nmj_area = getValue("Area");
			Table.set("NMJ_Area", imgdone, nmj_area, "NMJ_Data_AllImages");
			nmj_std = getValue("StdDev");
			Table.set("StdDev_"+CHANS[c], imgdone, nmj_std, "NMJ_Data_AllImages");
			selectImage("log_filtered");
			roiManager("Select", 0);
			Stack.setChannel(c+1);
			pix=getPixValues();
			if(DEBUG==true)
				print("[TS_Log]", fl[i]+ " line 233\n");
			rms = RMS(pix);
			Table.set("logrms_"+CHANS[c], imgdone, rms, "NMJ_Data_AllImages");
			}
		
		//Measure PCCS across whole NMJ
		for(c1=1; c1<channels+1; c1++)
			{
			if(channels>1)
				Stack.setChannel(c1);
			c1pix=getPixValues();
			if(DEBUG==true)
				print("[TS_Log]", fl[i]+ " line 245\n");
			for(c2=c1+1; c2<channels+1; c2++)
				{
				Stack.setChannel(c2);
				c2pix=getPixValues();
				if(DEBUG==true)
					print("[TS_Log]", fl[i]+ " before pcc\n");
				pcc=PCC(c1pix, c2pix);
				if(DEBUG==true)
					print("[TS_Log]", fl[i]+ " after pcc\n");
				Table.set("PCC_"+CHANS[c1-1]+"-"+CHANS[c2-1], imgdone, pcc, "NMJ_Data_AllImages");
				}
			}
		if(roiManager("count")>0)
			roiManager("reset");

		
		// ********************************
		// **  Analysis of architecture  **
		// ********************************

		// Get centroids for composite mesh
		compx=Table.getColumn("Comp_X", "CentroidData");
		compy=Table.getColumn("Comp_Y", "CentroidData");

		// Get and compare all meshes
		for(c1=0; c1<channels; c1++)
			{
			c1x=Table.getColumn(""+CHANS[c1] + "_X", "CentroidData");
			c1y=Table.getColumn(""+CHANS[c1] + "_Y", "CentroidData");
			
			for(c2=c1+1; c2<channels; c2++)
				{
				c2x=Table.getColumn(""+CHANS[c2]+"_X", "CentroidData");
				c2y=Table.getColumn(""+CHANS[c2]+"_Y", "CentroidData");
				
				//Calculate Centroid deltas & randomized centroid deltas
				if(c1x.length>1 && c2x.length>1)
					{
					centroidd12=nnd(c1x,c1y,c2x,c2y);
					Table.save(imgpath+"\\"+CHANS[c1]+"--"+CHANS[c2]+"_nnd.csv");
					} 
					else {centroidd12=newArray(NaN, NaN);}
				if(DEBUG==true)
					print("[TS_Log]", fl[i]+ " line 285\n");

				Table.set("HexE_"+CHANS[c1]+"-"+CHANS[c2], imgdone, centroidd12[2], "NMJ_Data_AllImages");	
				
				// Calculate jaccard index for mesh-mesh pairs of binary mesh masks
				ji=calculateJaccard("C"+c1+1+"_MeshMask", "C"+c2+1+"_MeshMask");
				Table.set("Jaccard_"+CHANS[c1]+"-"+CHANS[c2], imgdone, ji, "NMJ_Data_AllImages");
			
				}

			// Compare C1 to composite
			if(c1x.length>1 && compx.length>1)
				{
				centroiddxc=nnd(c1x,c1y,compx,compy);
				Table.save(imgpath+"\\"+CHANS[c1]+"--Composite_nnd.csv");
				}
				else {centroidxc=newArray(NaN, NaN);}
				
			
			Table.set("HexE_"+CHANS[c1]+"-Composite", imgdone, centroiddxc[2], "NMJ_Data_AllImages");	
			//Table.set("HexE-Random_"+CHANS[c1]+"-"+CHANS[c2], imgdone, centroid1rand, "NMJ_Data_AllImages");
				
			ji=calculateJaccard("C"+c1+1+"_MeshMask", "Composite_MeshMask");
			Table.set("Jaccard_"+CHANS[c1]+"-Composite", imgdone, ji, "NMJ_Data_AllImages");
			}		

		// *****************************************
		// *******  Analyze AZ properties   ********
		// *****************************************
		
		// This is largely convenience for finding and recording AZ later
		selectImage("AZ_objects");
		run("Analyze Particles...", "size=3-infinity pixel clear add");  //1.53m bug
		//run("Analyze Particles...", "size=3-infinity pixel clear add exclude");
		n_AZ = roiManager("count");
		if(DEBUG==true)
			print("[TS_Log]", fl[i]+ " line 158\n");
		// Create arrays for AZ location and intensity
		azx = newArray(n_AZ);
		azy = newArray(n_AZ);
		az_mean = newArray(n_AZ);
		az_area = newArray(n_AZ);
		selectImage("mean_norm");
		for(az=0; az<n_AZ; az++)
			{
			roiManager("select", az);
			x = getValue("X");
			y = getValue("Y");
			toUnscaled(x,y);
			azx[az] = x; azy[az] = y;

			az_mean[az] = getValue("Mean");
			az_area[az] = getValue("Area");
			}
		roiManager("reset");
		
		// *****************************************
		// **** Measure mesh/cell pixel values  ****
		// *****************************************

		// Analysis of composite mesh, including AZ-PAZ analysis and radial profiling
		// **************************************************************************

		// Get mesh rois for composite mesh. 
		// Filter anything touching edges and anything smaller than X 
		selectImage("Composite_mesh");
		run("Select All");
		run("Analyze Particles...", "size=12-infinity pixel clear add"); //1.53m bug
		//run("Analyze Particles...", "size=12-infinity pixel clear add exclude");
		
		// Exclude from analysis anything at the edge of the NMJ
		nROI = roiManager("count");
		selectImage("EDM");
		for(roi = 0; roi<nROI; roi++)
			{
			roiManager("select", roi);
			EDM_min = getValue("Min");
			if(EDM_min <= EDM_FILTER_MIN)
				{
				roiManager("Delete");
				nROI -= 1;
				roi -=1;
				}
			}
		
		if(nROI>0)
			{
			saveROIs(imgpath, "composite_filtered.roi");
			az_ids = newArray(nROI); // array holding total AZ int dens for each comp. unit
			az_areas = newArray(nROI); // array holding total AZ areas for each comp. unit
			az_n = newArray(nROI); // array to keep track of AZ number in each comp. unit
			
			// Analyze AZ properties in PAZ mesh units
			for(r=0; r<roiManager("count"); r++)
				{				
				// Intialize tables with 0 values. 
				// This ensures consistent ordering of critical count and ID columns
				// This also ensures meshes with 0 are recorded
				roiManager("Select", r);
				az_count = 0;
				az_meantotal = 0;
				az_areatotal = 0;
	
				// Record image specific data
				Table.set("AZcount", r, az_count, "Comp_NetData");
				Table.set("AZidtotal", r, az_meantotal*az_areatotal, "Comp_NetData");
				Table.set("AZareatotal", r, az_areatotal, "Comp_NetData");
				//Table.set("AZ"+az_count+"_xy", r, "("+azx[az]+","+azy[az]+")", "Comp_NetData");
				//Table.set("AZ"+az_count+"_id", r, az_id[az], "Comp_NetData");
	
				// Record data on compiled table
				Table.set("Image", total_comp_rows+r, fl[i], "Comp-ALLTHECELLS");
				Table.set("AZcount", total_comp_rows + r, az_count, "Comp-ALLTHECELLS");
				Table.set("AZidtotal", total_comp_rows+r, az_meantotal*az_areatotal, "Comp-ALLTHECELLS");
				Table.set("AZareatotal", total_comp_rows+r, az_areatotal, "Comp-ALLTHECELLS");
				
				// Check whether mesh unit contains AZ(s) and record
				for(az=0; az<n_AZ; az++)
					{
					if(Roi.contains(azx[az], azy[az]))
						{
						az_count += 1;
						az_meantotal += az_mean[az];
						az_areatotal += az_area[az];
						Table.set("AZcount", r, az_count, "Comp_NetData");
						Table.set("AZidtotal", r, az_meantotal*az_areatotal, "Comp_NetData");
						Table.set("AZareatotal", r, az_areatotal, "Comp_NetData");
						Table.set("AZ"+az_count+"_xy", r, "("+azx[az]+","+azy[az]+")", "Comp_NetData");
						Table.set("AZ"+az_count+"_id", r, az_mean[az]*az_area[az], "Comp_NetData");
	
						Table.set("Image", total_comp_rows+r, fl[i], "Comp-ALLTHECELLS");
						Table.set("AZcount", total_comp_rows + r, az_count, "Comp-ALLTHECELLS");
						Table.set("AZidtotal", total_comp_rows+r, az_meantotal*az_areatotal, "Comp-ALLTHECELLS");
						Table.set("AZareatotal", total_comp_rows+r, az_areatotal, "Comp-ALLTHECELLS");
						}
					}
				az_ids[r] = az_meantotal*az_areatotal;
				az_areas[r] = az_areatotal;
				az_n[r] = az_count;
				}
	
			// ******** Analyze cell and mesh data for composite PAZ
			for(c=1; c<=CHANS.length; c++)
				analyzeMesh(0, c, true);
	
			copyToCompiled("Comp");
				
			// ******** Do radial profiling for each channel 
			for(measChan=1; measChan<channels+1; measChan++)
				{
				if(DEBUG==true)
					print("[TS_Log]", fl[i]+" Line 391\n");
				
				for(roi=0; roi<roiManager("count"); roi++)
					{
					profile = radialProfile("mean_norm", roi, measChan);
					Table.set("Image", total_comp_rows + roi, fl[i], "Comp-ALLTHEPROFILES");
					Table.set("AZcount", total_comp_rows + roi, az_n[roi], "Comp-ALLTHEPROFILES");
					Table.set("AZidtotal", total_comp_rows + roi, az_ids[roi], "Comp-ALLTHEPROFILES");
					selectImage("EDM");
					roiManager("select", roi);
					EDM_min=getValue("Min");
					EDM_mean = getValue("Mean");
					Table.set("EDMmin", total_comp_rows + roi, EDM_min, "Comp-ALLTHEPROFILES");
					Table.set("EDMmean", total_comp_rows + roi, EDM_mean, "Comp-ALLTHEPROFILES");
					for(p=0; p<profile.length; p++)
						Table.set(CHANS[measChan-1]+":"+p, total_comp_rows + roi, profile[p], "Comp-ALLTHEPROFILES");
					}
				// Get average radial profile for the channel measChan using the composite mesh
				profile=AvgRadialProfile("mean_norm", measChan); // Note this creates a table called Profiles storing per-mesh profiles. Return is the avg profile

				if(DEBUG==true)
					print("[TS_Log]", fl[i]+" Line 404\n");
				Table.save(imgpath+"\\"+CHANS[measChan-1]+"-RadialProfs.csv", "Profiles");	
				close("Profiles");
	
				// Store the average radial profile for the whole (filtered) mesh in the compiled radial profile table for this mesh.
				Table.setColumn(CHANS[measChan-1]+"_"+stripX(fl[i]), profile, "Comp-Profiles");
				}
			roiManager("reset");
			}
		
		// Analysis of individual meshes
		// *****************************
		if(DEBUG==true)
			print("[TS_Log]", fl[i]+" Line 415\n");
		for(seg_c=1; seg_c<=CHANS.length; seg_c++)
			{
			selectImage(CHANS[seg_c-1]+"_mesh");
			run("Select All");
			run("Analyze Particles...", "size=12-infinity pixel clear add"); //1.53
			//run("Analyze Particles...", "size=12-infinity pixel clear add exclude"); 

			// Exclude from analysis anything at the edge of the NMJ
			nROI = roiManager("count");
			selectImage("EDM");
			for(roi = 0; roi<nROI; roi++)
				{
				roiManager("select", roi);
				EDM_min = getValue("Min");
				if(EDM_min <= EDM_FILTER_MIN)
					{
					roiManager("Delete");
					nROI -= 1;
					roi -=1;
					}
				}
			if(nROI>0)
				{
				if(DEBUG==true)
					print("[TS_Log]", fl[i]+" Line 480\n");
				saveROIs(imgpath, ""+CHANS[seg_c-1]+"_filtered.roi");
				az_n = newArray(nROI);
				az_ids = newArray(nROI);
				// If analyzing AZ mesh, measure all channels, AZ info, and do radial profiling
				if(seg_c==AZ)
					{
					// Analyze AZ properties in PAZ mesh units
					for(r=0; r<roiManager("count"); r++)
						{				
						// Intialize tables with 0 values. 
						// This ensures consistent ordering of critical count and ID columns
						// This also ensures meshes with 0 are recorded
						roiManager("Select", r);
						az_count = 0;
						az_idtotal = 0;
						az_areatotal = 0;
			
						// Record image specific data
						Table.set("AZcount", r, az_count, CHANS[seg_c-1]+"_NetData");
						Table.set("AZidtotal", r, az_idtotal, CHANS[seg_c-1]+"_NetData");
						Table.set("AZareatotal", r, az_idtotal, CHANS[seg_c-1]+"_NetData");
						//Table.set("AZ"+az_count+"_xy", r, "("+azx[az]+","+azy[az]+")", "Comp_NetData");
						//Table.set("AZ"+az_count+"_id", r, az_id[az], "Comp_NetData");
			
						// Record data on compiled table
						Table.set("Image", totalcellrows[AZ-1]+r, fl[i], CHANS[seg_c-1]+"-ALLTHECELLS");
						Table.set("AZcount", totalcellrows[AZ-1] + r, az_count, CHANS[seg_c-1]+"-ALLTHECELLS");
						Table.set("AZidtotal", totalcellrows[AZ-1]+r, az_idtotal, CHANS[seg_c-1]+"-ALLTHECELLS");
						Table.set("AZareatotal", totalcellrows[AZ-1]+r, az_idtotal, CHANS[seg_c-1]+"-ALLTHECELLS");
						// Check whether mesh unit contains AZ(s) and record
						for(az=0; az<n_AZ; az++)
							{
							if(Roi.contains(azx[az], azy[az]))
								{
								az_count += 1;
								az_idtotal += az_mean[az]*az_area[az];
								az_areatotal += az_area[az];
								Table.set("AZcount", r, az_count, CHANS[seg_c-1]+"_NetData");
								Table.set("AZidtotal", r, az_idtotal, CHANS[seg_c-1]+"_NetData");
								Table.set("AZareatotal", r, az_areatotal, CHANS[seg_c-1]+"_NetData");
								Table.set("AZ"+az_count+"_xy", r, "("+azx[az]+","+azy[az]+")", CHANS[seg_c-1]+"_NetData");
								Table.set("AZ"+az_count+"_id", r, az_mean[az]*az_area[az], CHANS[seg_c-1]+"_NetData");
			
								Table.set("Image", totalcellrows[AZ-1]+r, fl[i], CHANS[seg_c-1]+"-ALLTHECELLS");
								Table.set("AZcount", totalcellrows[AZ-1] + r, az_count, CHANS[seg_c-1]+"-ALLTHECELLS");
								Table.set("AZidtotal", totalcellrows[AZ-1]+r, az_idtotal, CHANS[seg_c-1]+"-ALLTHECELLS");
								Table.set("AZareatotal", totalcellrows[AZ-1]+r, az_areatotal, CHANS[seg_c-1]+"-ALLTHECELLS");
								}
							}
						az_n[r] = az_count;
						az_ids[r] = az_idtotal;
						}

					if(DEBUG==true)
						print("[TS_Log]", fl[i]+" Line 528\n");
						
					for(meas_c=1; meas_c<=CHANS.length; meas_c++)
						analyzeMesh(seg_c, meas_c, true);
					
					// ******** Do radial profiling for each channel 
					for(measChan=1; measChan<channels+1; measChan++)
						{
						if(DEBUG==true)
							print("[TS_Log]", fl[i]+" Line 531\n");
						
						for(roi=0; roi<roiManager("count"); roi++)
							{
							profile = radialProfile("mean_norm", roi, measChan);
							Table.set("Image", totalcellrows[AZ-1]+roi, fl[i], "AZ-ALLTHEPROFILES");
							Table.set("AZcount", totalcellrows[AZ-1]+roi, az_n[roi], "AZ-ALLTHEPROFILES");
							Table.set("AZidtotal", totalcellrows[AZ-1]+roi, az_ids[roi], "AZ-ALLTHEPROFILES");
							selectImage("EDM");
							roiManager("select", roi);
							EDM_min = getValue("Min");
							EDM_mean = getValue("Mean");
							Table.set("EDMmin", totalcellrows[AZ-1]+roi, EDM_min, "AZ-ALLTHEPROFILES");
							Table.set("EDMmean", totalcellrows[AZ-1]+roi, EDM_mean, "AZ-ALLTHEPROFILES");
							for(p=0; p<profile.length; p++)
								Table.set(CHANS[measChan-1]+":"+p, totalcellrows[AZ-1]+roi, profile[p], "AZ-ALLTHEPROFILES");
							}
						// Get average radial profile for the channel measChan using the composite mesh
						profile=AvgRadialProfile("mean_norm", measChan); // Note this creates a table called Profiles storing per-mesh profiles. Return is the avg profile
			
						if(DEBUG==true)
							print("[TS_Log]", fl[i]+" Line 404\n");
						Table.save(imgpath+"\\"+CHANS[measChan-1]+"-AZ_RadialProfs.csv", "Profiles");	
						close("Profiles");
			
						// Store the average radial profile for the whole (filtered) mesh in the compiled radial profile table for this mesh.
						Table.setColumn(CHANS[measChan-1]+"_"+stripX(fl[i]), profile, "AZ-Profiles");
						}
					
					}
				else 
					analyzeMesh(seg_c, seg_c, false);
				// For individual channel meshes, only analyze the same channel
			
				copyToCompiled(CHANS[seg_c-1]);
				roiManager("reset");
				}
			}
		if(DEBUG==true)
			print("[TS_Log]", fl[i]+" Line 428\n");
		
		// **** Save image-specific tables
		for(c=0; c<CHANS.length; c++)
			Table.save(imgpath+"\\"+CHANS[c]+"_NetData.csv", CHANS[c]+"_NetData");
		Table.save(imgpath+"\\Comp_NetData.csv", "Comp_NetData");
		close("*");
		imgdone+=1;
		}
	}

// *** Analyze AZ-PAZ comparison
if(HEATMAP==true)
	{
	AZ_PAZ_heatmap("Comp");	
	AZ_PAZ_heatmap(CHANS[AZ-1]);
	}

// *** Save compiled tables
for(c=0; c<CHANS.length; c++)
	{
	Table.save(savepath+"\\"+CHANS[c]+"_ALLTHECELLS.csv", CHANS[c]+"-ALLTHECELLS");
	Table.save(savepath+"\\"+CHANS[c]+"_NMJAverages.csv", CHANS[c]+"-NMJAverages");
	close(CHANS[c] + "-ALLTHECELLS");
	close(CHANS[c] + "-NMJAverages");
	close(CHANS[c] + "_NetData");
	}
Table.save(savepath+"\\Comp_ALLTHECELLS.csv", "Comp-ALLTHECELLS");
Table.save(savepath+"\\Comp_ALLTHEPROFILES.csv", "Comp-ALLTHEPROFILES");
Table.save(savepath+"\\AZ_ALLTHEPROFILES.csv", "AZ-ALLTHEPROFILES");
Table.save(savepath+"\\Comp_NMJAverages.csv", "Comp-NMJAverages");
Table.save(savepath+"\\NMJ_Data.csv", "NMJ_Data_AllImages");
Table.save(savepath+"\\Comp_Radial_Profiles.csv", "Comp-Profiles");
Table.save(savepath+"\\AZ_Radial_Profiles.csv", "AZ-Profiles");
close("Comp-ALLTHECELLS");
close("Comp-ALLTHEPROFILES");
close("AZ-ALLTHEPROFILES");
close("Comp-NMJAverages");
close("Comp_NetData");
close("Comp-Profiles");
close("D");
close("Profiles");
close("Profiles");
close("CentroidData");

close("NMJ_Data_AllImages");
close("*");

















// ***********************************************
// *******    Segmentation Functions    **********
// ***********************************************

function segment(imp)
	{
	total_rois = roiManager("count");
	// Make whole NMJ mask
	selectImage(imp);
	getDimensions(width, height, channels, slices, frames);
	if(SEGMENT==0)
		{
		run("Select All");
		run("Duplicate...", "title=norm duplicate channels=1-"+channels); //mask for whole NMJ
		getDimensions(width, height, channels, slices, frames);
		run("32-bit");
		newImage("mask", "32-bit black", width, height, slices);
		
		// normalize each channel to mean to combine evenly and add to mask
		for(c=0; c<COMBINE.length; c++)
			{
			selectImage("norm");
			run("Duplicate...", "title="+COMBINE[c]+" duplicate channels="+""+COMBINE[c]); //mask for whole cell
			run("Select All");
			overallmean=getValue("Mean");
			run("Divide...", "value="+overallmean);
			imageCalculator("Add", "mask", ""+COMBINE[c]);
			rename("mask");
			close(COMBINE[c]);
			}
		close("norm");
		selectImage("mask");
		run("16-bit");
		}
		else run("Duplicate...", "title=mask duplicate channels="+SEGMENT); //mask for whole cell
	
	// Blur, threshold, and erode if specified
	run("Gaussian Blur...", "sigma="+SIGMA+" stack");
	setAutoThreshold(METHOD+" dark stack");
	run("Convert to Mask", "method="+METHOD+" background=Dark black stack");
	run("Fill Holes", "stack");
	
	// Keep a copy of the un-eroded mask to clear background for AZ marker ID
	rename("full_mask");
	run("Create Selection");
	roiManager("add");
	run("Select All");
	run("Duplicate...", "title=mask");
	for(e=0; e<ERODE; e++)
		run("Erode");

	run("Create Selection");
	roiManager("add");
	roiManager("select", newArray(roiManager("count")-2, roiManager("count")-1));
	roiManager("save", imgpath+"\\NMJ_masks.zip");
	roiManager("select", newArray(roiManager("count")-2, roiManager("count")-1));
	roiManager("delete");
	run("Select All");
	
	// Duplicate mask and make into a distance map
	// For downstream spatial analysis and filtering out of edge mesh units
	run("Duplicate...", "title=EDM duplicate"); 
	run("Distance Map");
	}

function makeSingleChannelNetwork(img, mask, chan, sigma, tolerance, LoGr)
	{
	// Generate an image containing a mesh-like representation of a single input channel
	// Return an array of (x0, ..., xn, y0, ..., yn) describing the x,y centroid 
	// locations of n mesh units. 
	
	// Copy image to segment (to be processed, but not gradient mapped)
	total_rois = roiManager("count");
	selectImage(img);
	run("Select All");
	run("Duplicate...", "title=toSegment duplicate channels="+chan);
	setMinAndMax(0, 1);


	if(chan == AZ)
		{
		run("16-bit");			
	   	run("Mexican Hat Filter", "radius="+LoGr);
	   	
	   	// Clear outside of NMJ area for segmentation
	   	selectImage(mask);
	   	selectImage("full_mask");
		run("Create Selection");
		roiManager("add");
		selectImage("toSegment");
		roiManager("select", total_rois);
		run("Clear Outside");
		run("Select All");
		run("Invert");
		
		if(sigma>0)
    		run("Gaussian Blur...", "sigma="+sigma+" stack");
    		
		roiManager("select", total_rois);
		roiManager("delete");
		run("Select All");
		run("Duplicate...", "title=forMinima duplicate");
		}


	else {
		run("16-bit");
		
		// Clear outside of NMJ area for segmentation
		selectImage(mask);
		run("Select All");
		run("Create Selection");
		roiManager("add");
		selectImage("toSegment");
		roiManager("select", total_rois);
		max=getValue("Max");
		run("Enlarge...", "enlarge=-1 pixel");
		run("Make Inverse");
		setColor(max);
		fill();	
		roiManager("select", total_rois);
		roiManager("delete");
		
		run("Morphological Filters", "operation=[Internal Gradient] element=Octagon radius="+IG_RADIUS);
		rename("forMinima");
		if(POST_IG>0)
    		run("Gaussian Blur...", "sigma="+POST_IG+" stack");
		}

	// Find local minima to serve as seeds for watershedding
	selectImage("forMinima");
	run("Find Maxima...", "prominence="+tolerance+" exclude light output=[Single Points]");
	rename("BinMark");

	// For AZ channels, generate a mask for all AZ objects 
	// Then refine seeds after EDM-watershed of objects
	if(chan == AZ)
		{
		graymarker=binaryToGrayMarker("BinMark");

		run("Seeded Region Growing ...", "image=toSegment seeds="+graymarker+" stack=[Current slice only]");
		rename("AZ_objects");

		// do edm watershed on object to separate touching AZ
		// Note this will update in place "BinMark"
		watershed_az("AZ_objects");
		selectWindow("AZ_objects");
		run("Create Selection");
		if(selectionType()<0)
			run("Select All");
		roiManager("add");
		roiManager("select", roiManager("count")-1);
		roiManager("save", imgpath+"\\AZ_objects.roi");
		roiManager("select", roiManager("count")-1);
		roiManager("delete");
		close(graymarker);
		}
		
	// Create mesh-like rois from seeds, using the original channel (not the gradient) as the input
	run("Marker-controlled Watershed", "input=toSegment marker=BinMark mask="+mask+" binary calculate use");
	rename("SeededRegions");
	run("16-bit");
	
	// Create binary mesh representation and save centroid coordinates
	setThreshold(1, 66000);
	run("Convert to Mask");
	run("Analyze Particles...", "clear add");
	//run("Analyze Particles...", "clear add exclude");
	x1 = newArray(roiManager("count"));
	y1 = newArray(roiManager("count"));
	for(r=0; r<roiManager("count"); r++)
		{
		roiManager("Select", r);
		x = getValue("X"); 
		y = getValue("Y");
		toUnscaled(x,y);
		x1[r] = x;
		y1[r] = y;
		}
	
	saveROIs(imgpath, ""+CHANS[chan-1]+"_Mesh.roi");
	rename(""+CHANS[chan-1]+"_mesh");
	roiManager("reset");
	
	close("BinMark");
	close("forMinima");
	close("toSegment");
	roiManager("reset");
	centroids=Array.concat(x1,y1);
	return centroids;
	}
	
function makeCompositeNetwork(img, mask, sigma, tolerance)
	{
	// Generate an image containing a mesh-like representation of a composite of all channels
	// Return an array of (x0, ..., xn, y0, ..., yn) describing the x,y centroid 
	// locations of n mesh units. 
	total_rois = roiManager("count");
	selectImage(img);
	getDimensions(width, height, channels, slices, frames);
	run("Duplicate...", "title=Normed duplicate channels=1-"+channels);
	run("Gaussian Blur...", "sigma=1 stack");
	newImage("toSegment", "32-bit black", width, height, slices);
	
	// Add PAZ channels and subtract AZ channels together to create
	// composite image of all labels called 'toSegment'
	for(c=1; c<channels+1; c++)
		{
		selectImage("Normed");
		run("Select All");
		run("Duplicate...", "title=Ch"+c+" duplicate channels="+c);	
		if(c!=AZ)
			{
			toAdd=getTitle();
			imageCalculator("Add", "toSegment", toAdd);
			rename("toSegment");
			close(toAdd);
			}
		else {
			selectImage("Ch"+c);
			imageCalculator("Subtract", "toSegment", "Ch"+c);
			rename("toSegment");
			}
		close("Ch"+c);
		}
	close("Normed");
	
	// Apply 'internal gradient' filter with octagon element 
	// This sharpens PAZ edges for detection of local minima. Called 'forMinima'
	selectImage("toSegment");
	run("16-bit");
	run("Morphological Filters", "operation=[Internal Gradient] element=Octagon radius="+IG_RADIUS+1);
	rename("forMinima");

    if(POST_IG>0)
    	run("Gaussian Blur...", "sigma="+POST_IG+" stack");
	
	// Find local minima in composite forMinima image
	// Generate mesh representation within masked region using marker-controlled watershed
	run("Find Maxima...", "prominence="+tolerance+" exclude light output=[Single Points]");
	rename("BinMark");
	run("Marker-controlled Watershed", "input=toSegment marker=BinMark mask="+mask+" binary calculate use");
	rename("Comp_Mesh_Regions");
	
	// Create binary mesh representation and save centroid coordinates
	setThreshold(1, 66000);
	run("Convert to Mask");
	run("Analyze Particles...", "clear add");
	//run("Analyze Particles...", "clear add exclude");
	x1 = newArray(roiManager("count"));
	y1 = newArray(roiManager("count"));
	for(r=0; r<roiManager("count"); r++)
		{
		roiManager("select", r);
		x = getValue("X"); 
		y = getValue("Y");
		toUnscaled(x,y);
		x1[r] = x;
		y1[r] = y;
		}
		
	rename("Composite_mesh");
	close("BinMark");
	close("forMinima");
	close("toSegment");
	saveROIs(imgpath, "Composite_Mesh.roi");
	roiManager("reset");
	centroids=Array.concat(x1,y1);
	return centroids;
	}

function makeMeshMask(img)
	{
	// Make a binary mask that labels mesh elements 1 and everything else 0.
	// Used for jaccard index. 
	getDimensions(width, height, channels, slices, frames);
	newImage("meshmask", "8-bit black", width, height, 1);
	for(r=0; r<roiManager("count"); r++)
		{
		roiManager("select", r);
		run("Enlarge...", "enlarge="+JACCARDWIDTH/2);
		checkarea=getValue("Area");
		setColor(255);
		run("Fill", "slice");
		run("Enlarge...", "enlarge=-"+JACCARDWIDTH);
		setColor(0);
		if(getValue("Area")<checkarea)
			run("Fill", "slice");
		}
	return "meshmask";
	}


function AZ_By_Intensity(img)
	{
	// Generate AZ ROIs from a grayscale (instance) labeled image
	// Filter out background ROIs by rejecting objects with below-average AZ channel intensity
	
	selectImage(img);
	depth = bitDepth();
	if((depth != 16))
		run("16-bit");

	Stack.getDimensions(width, height, channels, slices, frames);
	run("Select All");
	
	maxval=0;
	if(depth==8)
		maxval=256;
	if(depth==16);
		maxval=65536;

	getHistogram(values, counts, maxval);
	
	for(v=0; v<counts.length; v++)
		{
		if(counts[v]>0)
			{
			selectImage(img);
			setThreshold(v, v);
			run("Create Selection");
			if(getValue("Area")>3)
				{
				roiManager("add");
				selectImage("mean_norm");
				roiManager("select", roiManager("count")-1);
				if(CHANS.length>1)
					Stack.setChannel(AZ);
				brp_mean = getValue("Mean");
				if(brp_mean<=AZ_INTENSITY_FILTER)
					roiManager("delete");
				}
			}
		}
	}


function minmaxNorm(img, mask, clip)
	{
	// Normalize image 0-1 min-max. 
	// Optionally, clip values by StdDev (pass 0 for no clipping)
	// Optionally, measure values only within mask (pass empty string to measure whole image)
	roiCount=roiManager("count");
	selectWindow(img);
	getDimensions(width, height, channels, slices, frames);
	run("Duplicate...", "title=min_max duplicate channels=1-"+channels);
	run("32-bit");
	
	if(mask!="")
		{
		selectWindow(mask);
		run("Create Selection");
		roiManager("add");
		}
	
	selectImage("min_max");	
	for(c=1; c<channels+1; c++)
		{
		run("Select All");
		if(mask!="")
			roiManager("select", roiCount);
		Stack.setChannel(c);
		getStatistics(area, mean, min, max, std, histogram);
		run("Select All");
		
		// Calculate the clipping threshold and assign as max
		if(clip>0)
			max=mean+std*clip;
			
		// Normalize all pixels 0-1 (if unclipped)
		run("Subtract...", "value="+min+" slice");
		run("Divide...", "value="+(max-min)+" slice");

		// If we are clipping, the replace all values over 1 (the max/clipping threshold) to 1
		if(clip>0)
			changeValues(1, 66000, 1);
		}
	roiManager("reset");
	}

function meanNorm(img, mask)
	{
	// Normalize img by the mean intensiy of image.
	// Optionally, measure values only within mask (pass empty string to measure whole image)
	roiCount=roiManager("count");
	selectWindow(img);
	getDimensions(width, height, channels, slices, frames);
	run("Duplicate...", "title=mean_norm duplicate channels=1-"+channels);
	run("32-bit");
	if(mask!="")
		{
		selectWindow(mask);
		run("Create Selection");
		roiManager("add");
		}
	
	selectImage("mean_norm");	
	for(c=1; c<channels+1; c++)
		{
		run("Select All");
		if(mask!="")
			roiManager("select", roiCount);
		Stack.setChannel(c);
		mean=getValue("Mean");
		run("Select All");
		run("Divide...", "value="+mean+" slice");
		}
	roiManager("reset");
	}

function watershed_az(img)
	{
	// Take as input a grayscale object map output
	// from seeded region growing. Convert objects to binary
	// Then do distance transform watershed. 
	// Produce new object and centroid map.
	start_rois = roiManager("count");
	AZ_By_Intensity(img);

	selectImage(img);
	run("Select All");
	run("Duplicate...", "title=AZ_binary");
	run("Select All");
	run("Clear");
	run("8-bit");
	setColor(255);
	for(roi=start_rois; roi<roiManager("count"); roi++)
		{
		selectImage("AZ_binary");
		roiManager("select", roi);
		fill();
		}
	run("Distance Transform Watershed", "distances=[Weights (5,7)] output=[16 bits] dynamic=0 connectivity=8");
	dtw = getTitle();
	close(img);
	close("AZ_binary");
	selectImage(dtw);
	rename(img);
	setThreshold(1, 66000);
	run("Make Binary");
	if(roiManager("count")>start_rois)
		{
		to_delete = getSequence(start_rois, roiManager("count")-1, 1);
		roiManager("select", to_delete);
		roiManager("delete");
		}

	// Update binary centroid map "BinMark" 
	selectImage(img);
	run("Select All");
	run("Analyze Particles...", "add");
	selectImage("BinMark");
	run("Select All");
	run("Clear");
	for(roi=start_rois; roi<roiManager("count"); roi++)
		{
		roiManager("select", roi);
		x=getValue("X");
		y=getValue("Y");
		toUnscaled(x,y);
		setPixel(x, y, 255);
		}
	
	if(roiManager("count")>start_rois)
		{
		to_delete = getSequence(start_rois, roiManager("count")-1, 1);
		roiManager("select", to_delete);
		roiManager("delete");
		}
	}
	
// ***********************************************
// ********** Analysis Functions *****************
// ***********************************************

function analyzeMesh(mesh, c_anal, pearson)
	{
	// analyze core and mesh regions for intensity etc
	// assumes an roi manager containing rois for all mesh units to analyze
	// mesh is the channel number (1,2,...) used to make the mesh
	// for composite mesh, mesh = 0
	// c_anal is the channel to analyze (1,2,...)
	// pearson = [true, false] whether to calculate pearsons between all channels

	if(mesh==0)
		total_rows = total_comp_rows;
	else total_rows = totalcellrows[mesh-1];
	
	if(mesh>0)
		c_label = CHANS[mesh-1];
		else c_label = "Comp";
		
	nROI = roiManager("count");
	for(r = 0; r<nROI; r++)
		{
		// get initial area estimate to check for when we shrink later
		// Need get on mean_norm image bc it is scaled (EDM not)
		selectImage("mean_norm");	
		roiManager("select", r);
		circ=getValue("Circ.");
		getStatistics(carea, xmean, xmin, xmax, xstd, xhistogram);

		// get EDM scores
		selectImage("EDM");
		roiManager("select", r);
		getStatistics(xarea, EDMmean, EDMmin, EDMmax, EDMstd);
		
		Table.set("Area", r, carea, c_label+"_NetData");
		Table.set("Circ.", r, circ, c_label+"_NetData");		
		Table.set("EDMmin", r, EDMmin, c_label+"_NetData");
		Table.set("EDMmean", r, EDMmean, c_label+"_NetData");
		
		Table.set("Image", total_rows + r, fl[i], c_label+"-ALLTHECELLS");
		Table.set("Area", total_rows + r, carea, c_label+"-ALLTHECELLS");
		Table.set("Circ.", total_rows + r, circ, c_label+"-ALLTHECELLS");
		Table.set("EDMmin", total_rows + r, EDMmin, c_label+"-ALLTHECELLS");
		Table.set("EDMmean", total_rows + r, EDMmean, c_label+"-ALLTHECELLS");
		
		// Shrink mesh ROI by half tracewidth and get pixel intensity data for core region
		selectImage("mean_norm");	
		roiManager("select", r);
		shrink=TRACEWIDTH/2;
		run("Enlarge...", "enlarge=-"+shrink+" pixel");
		if(CHANS.length>1)
			Stack.setChannel(c_anal);
		getStatistics(carea2, cmean, cmin, cmax, cstd);

		// For very small meshes, ROI may not shrink. If so, try smaller shrink factor.
		if(carea2==carea)
			shrink=TRACEWIDTH/4;	
		run("Enlarge...", "enlarge=-"+shrink+" pixel");
			
		getStatistics(carea2, cmean, cmin, cmax, cstd);
		c1pix=getPixValues();
		
		// Get LoG Filtered data for CoV	
		selectImage("log_filtered");	
		roiManager("select", r);
		run("Enlarge...", "enlarge=-"+shrink+" pixel");
		if(CHANS.length>1)
			Stack.setChannel(c_anal);		
		vals=getPixValues();
		clog_rms=RMS(vals);
		Table.set("CoreArea", r, carea2, c_label+"_NetData");					
		Table.set("CoreMean_"+CHANS[c_anal-1], r, cmean, c_label+"_NetData");
		Table.set("CoreCoV_"+CHANS[c_anal-1], r, cstd/cmean, c_label+"_NetData");
		Table.set("CoreLoGrms_"+CHANS[c_anal-1], r, clog_rms, c_label+"_NetData");

		Table.set("CoreArea", total_rows + r, carea2, c_label+"-ALLTHECELLS");
		Table.set("CoreMean_"+CHANS[c_anal-1], total_rows + r, cmean, c_label+"-ALLTHECELLS");
		Table.set("CoreCoV_"+CHANS[c_anal-1], total_rows + r, cstd/cmean, c_label+"-ALLTHECELLS");
		Table.set("CoreLoGrms_"+CHANS[c_anal-1], total_rows + r, clog_rms, c_label+"-ALLTHECELLS");
		
		//Get Pearson Comps to this channel if selected
		if(pearson == true)
			{
			selectImage("mean_norm");
			for(c2=c_anal+1; c2<=channels; c2++)
				{
				Stack.setChannel(c2);
				c2pix=getPixValues();
				pcc=PCC(c1pix, c2pix);
				Table.set("CorePCC:"+CHANS[c_anal-1]+"-"+CHANS[c2-1], r, pcc, c_label+"_NetData");
				Table.set("CorePCC:"+CHANS[c_anal-1]+"-"+CHANS[c2-1], total_rows + r, pcc, c_label+"-ALLTHECELLS");
				}
			}
		
		// analyze meshy part
		// Measure pixel values (Int, CoV, PCC) in each reticulum

		// create a mesh ROI as xor of ROI shrunken and expanded by tracewidth/2
		roiManager("select", r);
		run("Enlarge...", "enlarge=-"+shrink+" pixel");
		roiManager("add");
		roiManager("select", r);
		run("Enlarge...", "enlarge="+TRACEWIDTH/2+" pixel");
		roiManager("add");
		roiManager("select", newArray(nROI, nROI+1));
		roiManager("XOR");
		roiManager("add");
		roiManager("delete");
		// Now we have original rois intact but with new 'mesh' roi at position nROI. 
		// We will remove this after measurements for this roi are done.

		// Measure basic intensity values
		selectImage("mean_norm");
		roiManager("select", nROI);
		if(CHANS.length>1)
			Stack.setChannel(c_anal);
		getStatistics(retarea, retmean, retmin, retmax, retstd);
		
		selectImage("log_filtered");
		roiManager("select", nROI);
		if(CHANS.length>1)
			Stack.setChannel(c_anal);								
		c1pix=getPixValues();
		log_rms=RMS(c1pix);
		Table.set("MeshArea", r, retarea, c_label+"_NetData");
		Table.set("MeshMean_"+CHANS[c_anal-1], r, retmean, c_label+"_NetData");
		Table.set("MeshCoV_"+CHANS[c_anal-1], r, retstd/retmean, c_label+"_NetData");
		Table.set("MeshLoGrms_"+CHANS[c_anal-1], r, log_rms, c_label+"_NetData");
		Table.set("MeshRatio_"+CHANS[c_anal-1], r, retmean/cmean, c_label+"_NetData");

		Table.set("MeshArea", total_rows + r, retarea, c_label+"-ALLTHECELLS");
		Table.set("MeshMean_"+CHANS[c_anal-1], total_rows + r, retmean, c_label+"-ALLTHECELLS");
		Table.set("MeshCoV_"+CHANS[c_anal-1], total_rows + r, retstd/retmean, c_label+"-ALLTHECELLS");
		Table.set("MeshLoGrms_"+CHANS[c_anal-1], total_rows + r, log_rms, c_label+"-ALLTHECELLS");
		Table.set("MeshRatio_"+CHANS[c_anal-1], total_rows + r, retmean/cmean, c_label+"-ALLTHECELLS");

		//Get Pearson Comps to this channel
		if(pearson == true)
			{
			selectImage("mean_norm");
			for(c2=c_anal+1; c2<=channels; c2++)
				{
				Stack.setChannel(c2);
				c2pix=getPixValues();
				pcc=PCC(c1pix, c2pix);
				Table.set("MeshPCC:"+CHANS[c_anal-1]+"-"+CHANS[c2-1], r, pcc, c_label+"_NetData");
				Table.set("MeshPCC:"+CHANS[c_anal-1]+"-"+CHANS[c2-1], total_rows + r, pcc, c_label+"-ALLTHECELLS");						
				}
			}

		// delete added mesh roi to restore original ROI list
		roiManager("select", nROI);
		roiManager("delete");
		}

	// Analyze Haralick's textures for each channel per ROI
	GLCM("mean_norm", c_anal);	
	
	// Copy columns from GLCMtemp to appropriate data tables
	glcmcols=split(Table.headings("GLCMtemp"));
	for(col=0; col<glcmcols.length; col++)
		{
		tempcol=Table.getColumn(glcmcols[col], "GLCMtemp");
		Table.setColumn(glcmcols[col], tempcol, c_label+"_NetData");
		for(row=0; row<tempcol.length; row++)
			Table.set(glcmcols[col], total_rows+row, tempcol[row], c_label+"-ALLTHECELLS");
		}
	close("GLCMtemp");
	}


function AZ_PAZ_heatmap(c_label)
	{
	/* 
	 *  Do core AZ-PAZ comparison:
	 *  Generate a heatmap of all mesh metrics
	 *  Create a new table AZ-PAZ_PCC to measure PCC
	 *  between PAZ metrics and AZ integrated density
	 *  
	 *  While this could in principle be applied to any mesh,
	 *  currently only recording AZ data for the composite mesh.
	 *  
	 *  Note - ROIs are already filtered to remove edge ROIs (EDM min = 1)
	 */

	 
	//get array of column names for compiled data (cols[0] will contain file names, cols[1] has AZ (x,y) coordinates)
	cols=split(Table.headings(c_label+"-ALLTHECELLS"));  

	// Remove file and coordinate columns to get list of column headings to be compared to AZ-ID
	compcols=Array.slice(cols, 2, cols.length);
	azids=Table.getColumn("AZidtotal", c_label+"-ALLTHECELLS");

	
	// Remove NaNs from azids
	azidsFilt=Array.deleteValue(azids, NaN);
	
	// Sort azidsfilt in ascending value for heatmapping:
	azidsSorted=reSort(azidsFilt, azidsFilt);
	azidsSortedNorm=normArray(azidsSorted, 65535);
	//Calculate Pearson R for each column compared to AZ-ID
	pccs=newArray(compcols.length);

	// Make table with columns grouped thematically (in same order analyzed)
	// 40 pixels * however many columns in width by 5 pixels * however many cells there are to display
	newImage("AZ-PAZ_Heatmap", "16-bit black", 40*(compcols.length), 5*azidsFilt.length, 1); 
	
	for(col=0; col<compcols.length; col++)
		{
		// Retrieve column to be compared to azid
		tempcol=Table.getColumn(compcols[col], c_label + "-ALLTHECELLS");
		// First filter the column to be compared to remove any rows in which azid==NaN:
		compFilt=filterNaN(azids, tempcol);
		
		// Calculate PCC between NaN filtered azid and comp column arrays, and append to pccs
		// This will be used to order the sorted heatmap:
		pccs[col]=PCC(azidsFilt, compFilt);
		
		// Sort the compared column by the value of azid:
		tempSorted=reSort(azidsFilt, compFilt);
		tempNorm=normArray(tempSorted, 65535);
		for(n=0; n<tempSorted.length; n++)
			{
			makeRectangle(40*(col), n*5, 40, 5);
			run("Add...", "value="+tempNorm[n]);
			}
		}

	run("Rotate 90 Degrees Left");
	run("ICA3");
	saveAs("tiff", savepath+"\\"+c_label+"_AZ-PAZ_Heatmap.tif");


	Table.create("AZ-PAZ_PCC");
	Table.set("Type", 0, "PCC_to_az-id");
	for(col=0; col<compcols.length; col++)
		Table.set(compcols[col], 0, pccs[col]);

	Table.save(savepath+"\\"+c_label+"_AZ-PAZ_PCCs.csv", "AZ-PAZ_PCC");
	// For sorted table
	
	// Sort column headings by the magnitude of their PCC to azid:
	colsSorted = reSort(pccs, compcols);
	
	newImage("AZ-PAZ_HeatmapSorted", "16-bit black", 40*(compcols.length), 5*azidsFilt.length, 1); // 40 pixels * however many columns in width by 5 pixels * however many cells there are to display
	for(col=0; col<compcols.length; col++)
		{
		tempcol=Table.getColumn(colsSorted[col], c_label + "-ALLTHECELLS");
		
		// First filter the column to be compared to remove any rows in which azid==NaN:
		compFilt=filterNaN(azids, tempcol);
		
		// Sort the compared column by the value of azid:
		tempSorted=reSort(azidsFilt, compFilt);
		tempNorm=normArray(tempSorted, 65535);
		for(n=0; n<tempSorted.length; n++)
			{
			makeRectangle(40*(col), n*5, 40, 5);
			run("Add...", "value="+tempNorm[n]);
			}
		}
	
	run("Select All");
	run("Rotate 90 Degrees Left");
	run("ICA3");
	saveAs("tiff", savepath+"\\"+c_label+"_AZ-PAZ_HeatmapSorted.tif");

	print("AZ-PAZ_Heatmap.tif is a heatmap of " + azids.length + " subcell values compiled from " + imgdone + " NMJs");
	print("Values are normalized from -2.2 to +2.2 StdDev from mean calculated from pooled subcell data for each measurement");
	print("Shown are " + compcols.length + " different measurements arranged by column as follows:");
	Array.print(compcols);

	print("");
	print("AZ-PAZ_HeatmapSorted.tif is a heatmap of " + azids.length + " subcell values compiled from " + imgdone + " NMJs");
	print("Values are normalized from -2.2 to +2.2 StdDev from mean calculated from pooled subcell data for each measurement");
	print("Shown are " + colsSorted.length + " different measurements arranged by column as follows:");
	
	Array.print(colsSorted);

	selectWindow("Log");
	saveAs("Text", savepath+"\\"+c_label+"_Heatmap labels and info.txt");
	close("Log");

	}

function getStatsNaN(array)
	{
	// Function to get stats from an array ignoring NaNs
	filtarray=Array.deleteValue(array, NaN);
	Array.getStatistics(filtarray, min, max, mean, std);
	return newArray(min, max, mean, std);
	}

function GLCM(img, chan)
	{
	/*  For image img, iterate through all ROIs in ROI manager
	 *  Create a straightened (rectangular) representation of each
	 *  mesh ROI, then analyze haralick's texture measures
	 *  over range of distances from 1 to GLCMmax for the indicated channel.
	 *  
	 */
	
	selectImage(img);
	getDimensions(width, height, channels, slices, frames);
	run("Select All");
	run("Duplicate...", "title=temp duplicate channels="+chan);
	selectImage("temp");
	run("Select All");
	max = getValue("Max");
	min = getValue("Min");

	nroi=roiManager("Count");
	Table.create("GLCMtemp");
	for(r=0; r<nroi; r++)
		{
		roiManager("Select", r);
		run("Interpolate", "interval=1");	
		getSelectionCoordinates(xs, ys);
		newx=Array.deleteIndex(xs,0); newy=Array.deleteIndex(ys,0);
		makeSelection("freeline", newx, newy);
		run("Straighten...", "title=straightened line="+TRACEWIDTH+" process");
		setMinAndMax(min, max);
		run("8-bit");

		run("Select All");
		for(d=1; d<GLCMdmax+1; d++)
			{
			d+=dSTEP-1;
			run("GLCM Texture", "enter="+d+" select=[0 degrees] angular contrast correlation inverse entropy");
			for(col=0; col<HEADERS.length; col++)
				{
				tempval=Table.get(HEADERS[col], 0, "Results");
				Table.set(NEWHEADERS[col]+"D"+d+"_"+CHANS[chan-1], r, parseFloat(tempval), "GLCMtemp");
				}
			}
		close("straightened");
		}
	close("temp");
	}
	
	
	
function copyToCompiled(c_label)
	{
	// Calculate average values for all data and add to composite data table
	cols=split(Table.headings(c_label+"_NetData"));
	for(col=0; col<cols.length; col++)
		{
		tempvals=Table.getColumn(cols[col], c_label + "_NetData");
		// Get stats for tempvals, ignoring NaN. return is (Min, max, mean, std)
		stats=getStatsNaN(tempvals);
		Table.set(cols[col], imgdone, stats[2], c_label + "-NMJAverages");
		Table.set(cols[col]+"_Stdv", imgdone, stats[3], c_label + "-NMJAverages");
		}
	}
	
function nnd(x1,y1,x2,y2)
	{
	// Calculate the density weighted dimensionless 
	// nearest neighbor distance between two sets of x,y points
	selectImage("full_mask");
	run("Create Selection");
	nmjarea=getValue("Area");
	nmjarea=nmjarea/SCALE/SCALE;
	x1=Array.deleteValue(x1, NaN);
	x2=Array.deleteValue(x2, NaN);
	y1=Array.deleteValue(y1, NaN);
	y2=Array.deleteValue(y2, NaN);
	
	// Calculate nnd1-2 and nnd2-1 for each set of points 1 and 2
	eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('D');");
	//Table.create("D");
	d12=newArray(x1.length);
	Array.fill(d12, 1000000000000000000000000000000000000000000000000000000000);
	for(i1=0; i1<x1.length; i1++)
		{		
		for(i2=0; i2<x2.length; i2++)
			{
			di=sqrt((pow((x2[i2]-x1[i1]), 2)+pow((y2[i2]-y1[i1]),2)));
			if(di<d12[i1])
				{
				d12[i1]=di;		
				keep1i=i2;
				}
			}
		
		d12[i1]=d12[i1]/SCALE; //convert pixel length to microns

		Table.set("d12-Ch1x", i1, x1[i1], "D");
		Table.set("d12-Ch1y", i1, y1[i1], "D");
		Table.set("d12-Ch2x", i1, x2[keep1i]);
		Table.set("d12-Ch2y", i1, y2[keep1i], "D");
		Table.set("1-2delta", i1, d12[i1], "D");
		}
		
	Array.getStatistics(d12, min12, max, mean12, stdDev);
	
	d21=newArray(x2.length);
	Array.fill(d21, 1000000000000000000000000000000000000000000000000000000000);
	
	for(i1=0; i1<x2.length; i1++)
		{		
		for(i2=0; i2<x1.length; i2++)
			{
			di=sqrt((pow((x1[i2]-x2[i1]), 2)+pow((y1[i2]-y2[i1]),2)));
			if(di<d21[i1])
				{
				d21[i1]=di;		
				keep2i=i2;
				}
			}
		
		d21[i1]=d21[i1]/SCALE; //convert pixel length to microns
	
		Table.set("d21-Ch2x", i1, x2[i1], "D");
		Table.set("d21-Ch2y", i1, y2[i1], "D");
		Table.set("d21-Ch1x", i1, x1[keep2i], "D");
		Table.set("d21-Ch1y", i1, y1[keep2i], "D");
		Table.set("2-1delta", i1, d21[i1], "D");
		}
	Table.update;
	Array.getStatistics(d21, min21, max, mean21, stdDev);
	dens1=x1.length/nmjarea;
	dens2=x2.length/nmjarea;
	wavg=(x1.length*mean12+x2.length*mean21)/(x1.length+x2.length);
	
	dlesswavg=wavg/ (( (dens1/sqrt(dens2)) + (dens2/sqrt(dens1)) )  /(dens1+dens2));

	poissonError=(0.5-dlesswavg)/0.5*100;
	hexError=(0.4-dlesswavg)/0.4*100;
	
	return newArray(dlesswavg,poissonError,hexError);
	//return newArray(hexTA-wavg, spotTA-wavg);
	}

function calculateJaccard(img1, img2)
	{
	// Calculate the Jaccard overlap index for two binary masks
	// Jaccard index is the fraction of overlapping foreground pixels (AND) compared to total foreground pixels (OR)
	
	imageCalculator("OR create", img1, img2);
	getHistogram(values, OR_counts, 256);
	imageCalculator("AND create", img1, img2);
	getHistogram(values, AND_counts, 256);
	jaccard=AND_counts[255]/OR_counts[255];
	return jaccard;
	}

function PCC(X,Y)
	{
	// Calculate the Pearson Correlation coefficient between two arrays
	if(X.length != Y.length)
		exit("array lengths must match in PCC");
	nan=false;
	// First check whether either array contains nan and remove these elements from both
	// do check first bc this is slow. The built in imagej 'filter' is useless for this :(
	for(i=0; i<X.length; i++)
		if(X[i]!=X[i])
			nan=true;
	if(nan==true)
		{
		y_filt = filterNaN(X, Y);
		x_filt = filterNaN(X, X);
		}
		else
		{
		x_filt = X;
		y_filt = Y;
		}
		
	nan=false;
	for(i=0; i<y_filt.length; i++)
		if(y_filt[i]!=y_filt[i])
			nan=true;
	if(nan==true)
		{
		x_filt = filterNaN(y_filt, x_filt);
		y_filt = filterNaN(y_filt, y_filt);	
		}
	Array.getStatistics(x_filt, Xmin, Xmax, Xmean, XstdDev);
	Array.getStatistics(y_filt, Ymin, Ymax, Ymean, YstdDev);

	xy=0; xsq=0; ysq=0;
	for(n=0; n<x_filt.length; n++)
		{
		x=x_filt[n]-Xmean;
		y=y_filt[n]-Ymean;
		xy+=(x*y);
		xsq+=(x*x);
		ysq+=(y*y);
		}

	pearson=(xy/sqrt(xsq*ysq));

	return pearson;
	}

function radialProfile(img, roi, chan)
	{
	// Calculate average radial profile for a single mesh unit.
	// Return the average profile as an array.
	selectImage(img);
	getDimensions(width, height, channels, slices, frames);
	roiManager("select", roi);
	area=getValue("Area");
	circ=getValue("Circ.");
	
	if((area>minAREA)&&(circ>minCIRC))
		{
		run("Interpolate", "interval=1");
		getSelectionCoordinates(xs, ys);
		centX=getValue("XM")*SCALE;
		centY=getValue("YM")*SCALE;
		// For each point in ROI perimeter, make a line to centroid and get profile for each channel
		selectImage(img);
		Stack.setChannel(chan);
		tempprofile=newArray(SIZE);
		count=0;
		for(p=0; p<xs.length; p++)
			{
			makeLine(xs[p], ys[p], centX, centY);
			profile=getProfile();
			if(profile.length>2)
				{
				profile=Array.resample(profile, SIZE);
				tempprofile=addArrays(tempprofile, profile);
				count+=1;
				}
			}
			
		avgprofile=newArray(SIZE);
		for(i=0; i<SIZE; i++)
			avgprofile[i]=tempprofile[i]/count;
		} else
			{
			avgprofile=newArray(SIZE);
			for(i=0; i<SIZE; i++)
				avgprofile[i]=NaN;
			}
	return avgprofile;
	}	
	
function AvgRadialProfile(img, chan)
	{
	// Calculate the average radial profile
	// For an image containing many mesh units.
	// Store each mesh profile in table Profiles
	// Return average profile.
	
	Table.create("Profiles");
	nRoi=roiManager("count");
	cumulativeProfile=newArray(SIZE);
	count=0;
	for(r=0; r<nRoi; r++)
		{
		prof=radialProfile(img, r, chan);
		Table.setColumn(r, prof, "Profiles");
		if(prof[0]>0)
			{
			nprof=normMinMax(prof);
			cumulativeProfile=addArrays(cumulativeProfile, nprof);
			count+=1;
			}
		}
	for(i=0; i<cumulativeProfile.length; i++)
		cumulativeProfile[i]=cumulativeProfile[i]/count;
	return cumulativeProfile;
	}

// ***********************************************
// ********** Utility Functions ******************
// ***********************************************
function addArrays(a1, a2)
	{
	// Add two arrays of the same length element-wise
	temparray=newArray(a1.length);
	if(a1.length!=a2.length)
		exit("Trying to add arrays; must be same length");
	for(i=0; i<a1.length; i++)
		temparray[i]=a1[i]+a2[i];
	return temparray;
	}
	
function binaryToGrayMarker(binMarker)
	{
	// Convert an image containing binary semantic labels to grayscale instance labels
	// Return title of gray marker image "GrayMark"
	start_rois = roiManager("count");
	selectImage(binMarker);
	run("Select All");
	getDimensions(width, height, channels, slices, frames);
	nRois=roiManager("count");
	run("Analyze Particles...", "add");
	newImage("GrayMark", "8-bit black", width, height, 1);
	for(r=start_rois; r<roiManager("count"); r++)
		{
		roiManager("select", r);
		x=getValue("X"); y=getValue("Y");
		setPixel(x, y, r-start_rois+1);
		}
	if(roiManager("count")>start_rois)
		{
		toRemove = getSequence(start_rois, roiManager("count")-1, 1);
		roiManager("select", toRemove);
		roiManager("delete");
		}
	return getTitle();
	}

function filterNaN(array1, array2)
	{
	// Remove array2 positions where array1 == NaN
	if(array1.length != array2.length)
		exit("exception in filterNaN: array1 and array2 must be same length");
	filt1=newArray(0);
	filt2=newArray(0);
	for(i=0; i<array1.length; i++)
		{
		if(isNaN(array1[i])==false)
			{
			val1=newArray(1); val1[0]=array1[i];
			val2=newArray(1); val2[0]=array2[i];
			filt1=Array.concat(filt1, val1);
			filt2=Array.concat(filt2, val2);
			}
		}
	return filt2;
	}
	
function getImgCount(path)
	{
	// Sift an array of filenames and return count of images matching IMAGEEXT
	count=0;
	fl=getFileList(path);
	for(n=0; n<fl.length; n++)
		{
		if(endsWith(fl[n], IMAGEEXT)==true)
			count=count+1;
		}
	return count;
	}

function getSequence(start, stop, interval)
	{
	// Return a sequence whose first value is start
	// whose last value is less than or equal to stop
	// and which increments by interval
	size = floor((stop-start+interval)/interval);
	ar = newArray(size);
	for(ix=0; ix<size; ix++)
		{
		val = start+ix*interval;
		ar[ix] = val;
		}
	return ar;
	}

function getPixValues()
	{
	// Get pixel intensity values for all points within an ROI and return as an array
	Roi.getContainedPoints(xpoints, ypoints);
	vals=newArray(xpoints.length);
	for(n=0; n<xpoints.length; n++)
		vals[n]=getPixel(xpoints[n], ypoints[n]);
	return vals;
	}

function normArray(a, scale)
	{
	// Normalize array from 0 - scale, clipping values at zhi and zlo.
	Array.getStatistics(a, amin, amax, amean, astdDev);
	zhi=amean+Z*astdDev;
	zlo=amean-Z*astdDev;
	new=newArray(a.length);
	for(n=0; n<a.length; n++)
		{
		new[n]=(a[n]-zlo)/(zhi-zlo)*scale;
		if((a[n]-zlo)/(zhi-zlo)>1)
			{new[n]=scale;}
		if((a[n]-zlo)/(zhi-zlo)<0)
			{new[n]=0;}
		}
	return new;	
	}
	
function normMinMax(a)
	{
	// Normalize an array from 0-1 min-max
	Array.getStatistics(a, amin, amax, amean, astdDev);
	new=newArray(a.length);
	for(n=0; n<a.length; n++)
		new[n]=(a[n]-amin)/(amax-amin);
	return new;
	}

function reSort(guide, a)
	{
	// Sort positions of array 'a' by the rank order of array 'guide'
	rank=Array.rankPositions(guide);
	new=newArray(a.length);
	for(z=0; z<a.length; z++)
		{
		new[z]=a[rank[z]];
		}
	return new;
	}
	
function RMS(array)
	{
	// Calculate the root mean squared of array
	squared=newArray(array.length);
	for(i=0; i<array.length; i++)
		{
		squared[i]=array[i]*array[i];
		}
	Array.getStatistics(squared, min, max, mean, stdDev);
	rt=sqrt(mean);
	return rt;
	}
	
function saveROIs(path, imp)
	{
	// Save all ROIs in Roi list
	rois=newArray(roiManager("count"));
	for (n=0; n<rois.length; n++)
		{
		rois[n]=n;
		}
	roiManager("select", rois);
	roiManager("Save", path+"\\"+stripX(imp)+"ROIs.zip");
	}
	
function stripX(string)
	{
	// This is because Macro language doesn't have a general use name without extension
	return substring(string, 0, lastIndexOf(string, "."));
	}

// VERSIONS
/*
 *  211 - Bug fix where tracewidth wasn't being implemented
 *  212 - Parameterized the az intensity filter used in AZ_By_Intensity
 *  213 - Bug fix - EDM filtering == 1 instead of <= 1, so not excluding EDM 0. 
 *        Parameterized EDM filter range.
 *  214 - Bug fix - using unscaled nmj area for centroid error measurements
 *  	- Radial profiles now also exported in ALLTHEPROFILES format, with EDM and AZ data 
 *  	- cleaned up how columns and filenames are saved
 *  216 - Fixed bug where jaccard and centroid delta weren't being calculated for composite vs channel 3 pairing
 *  
 */