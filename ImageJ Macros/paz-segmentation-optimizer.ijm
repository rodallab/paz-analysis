// To do:
// note outlier bad segmentations?
// test different post-gradient blurs
// 
DEBUG=false;
BATCHMODE = true;
RANK=true;
VERSION="206_AMS020_comp-ig5_post-ig1";
CHANS=newArray("Nwk", "FasII", "BRP"); // List of channels
AZ=3; //Channel containing active zone marker for analysis PUT 0 if none
AREABONUS = true; // If getting undersegmentation set to false
IG_RADIUS = 4;
POST_IG = 1;  // sigma for blur after gradient filter
// Set range of tolerances to test here
// Note tolerances for AZ channels will be multiplied by 10.
//tolerances = newArray(5, 10, 20, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25000);
tolerances = newArray(25, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800);
//tolerances = newArray(600, 1800);


// *** Basic info
FILTER=true;
IMAGEEXT="tif";

ITERATIONS = 2;
SPACE = 0.2; // how big of a step to generate new tolerance values (fraction)
// Note: SPACE will scale as SPACE/iteration on each iteration.
RANGE = 6; // how many values to spread for new tolerance values

// *** Whole NMJ segmentation parameters
PREBGRB=50; // Rolling ball radius to apply BG subtraction to whole image prior to any seg or analysis. 0=none.
SIGMA=1; //This is the amount to blur your image to make the mask
SEGMENT=0; //This is the channel to use to make a mask. Put 0 to combine channels (see below)
COMBINE=newArray(1,2,3); // Which channels to add together to make mask IF segment==0.
ERODE=3; //This is the number of times to erode the NMJ mask
METHOD="Li";

// *** Mesh segmentation parameters
bgrbs = newArray(0, 0, 0); // radius for additional background subtraction while mesh segmenting. 0 for none.
comp_bgrb = 0;
TRACEWIDTH = 4;
logr = 2; // sigma to log filer AZ prior to segmentation
seg_sigma = 1;

// Initialize path and tables
path=getDirectory("Select a folder to analyze ...");
//path="Z:\\Current members\\DelSignore\\PAZ Organization Analysis Project\\SD596_CompositePAZPlanning\\test\\";
savepath=path+"\\OptimizeV-"+VERSION;
startTime=getTime();
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
showText("Experiment_Log", "Start: "+hour+":"+minute);
if(DEBUG==true)
	showText("DEBUG_LOG", "Start: "+hour+":"+minute+"\n");
File.makeDirectory(savepath);
comp_tolerances = tolerances;
if(roiManager("count")>0)
	roiManager("reset");
close("*");
run("Colors...", "foreground=white background=black selection=yellow");
run("Options...", "iterations=1 count=1 black");
if(roiManager("count")>0)
	roiManager("reset");
setBatchMode(BATCHMODE);

for(c=0; c<CHANS.length; c++)
	{
	tablename=CHANS[c]+"-Ranks";
	eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('"+tablename+"');");
	}

tablename="Comp-Ranks";
eval("js", "rt=new ResultsTable(); rt.setNaNEmptyCells(true); rt.show('"+tablename+"');");
count=getImgCount(path);
files=newArray(count);

for(iteration = 1; iteration<=ITERATIONS; iteration++)
	{
	print("[Experiment_Log]", "\n");
	fl=getFileList(path);
	imgdone=0;
	for(i=0; i<fl.length; i++)
		{
		if(endsWith(fl[i], IMAGEEXT)==true)
			{
			imgpath=savepath+"\\"+stripX(fl[i]);
			if(File.exists(imgpath))
				File.delete(imgpath);
			File.makeDirectory(imgpath);
			open(path+""+fl[i]);
			print("[Experiment_Log]", "\n File "+imgdone+1+": "+fl[i]);
			
			for(c=0; c<CHANS.length; c++)
				Table.set("File", imgdone, fl[i], CHANS[c]+"-Ranks");
			Table.set("File", imgdone, fl[i], "Comp-Ranks");
			Stack.getDimensions(width,height,channels,slices,frames);
			nmjmeans=newArray(channels);
			
			if(channels!=CHANS.length)
				exit("CHANS (channel labels) does not contain the correct number of labels");
			
			run("Select All");
			if(PREBGRB>0)
				run("Subtract Background...", "rolling="+PREBGRB+" stack");
			
			// Generate tight whole NMJ mask ("mask") and distance map ("EDM")
			// for downstream analysis and segmentations
			segment(fl[i]);
			selectImage("mask");
			run("Create Selection");
			roiManager("add");
			roiManager("select", roiManager("count")-1);
			roiManager("save", imgpath+"\\nmj_mask.roi");
	
			// Normalize each channel by mean for downstream analysis. Called mean_norm
			meanNorm(fl[i], "mask");
	
			// Normalize each channel by min-max, clipping high values for downstream seg. Called min_max
			minmaxNorm(fl[i], "mask", 3);
	
			// ********************************************************
			// ********  Do optimizing of tolerance value  ************
			// ********************************************************
			selectImage("mask");
			run("Create Selection");
			nmj_area = getValue("Area");
	
			for(c=1; c<=CHANS.length; c++)
				{
				// Update/refine tolerances based on prior iteration
				if(iteration > 1)
					tolerances = spread(c_bests[c-1], RANGE, SPACE/(iteration-1));
				col_labels = newArray(tolerances.length);
				for(col_i=0; col_i<tolerances.length; col_i++)
					col_labels[col_i] = ""+iteration+"_"+toString(tolerances[col_i]);					

				// Do the optimization for given channel and tolerances
				ranks = optimizeSC(c, tolerances);
				
				setRow(""+CHANS[c-1]+"-Ranks", col_labels, ranks, imgdone);
				Table.update(""+CHANS[c-1]+"-Ranks");
				}

			if(iteration > 1)
				comp_tolerances = spread(comp_best, RANGE, SPACE/(iteration-1));

			comp_labels = newArray(comp_tolerances.length);
			for(col_i=0; col_i<comp_tolerances.length; col_i++)
				comp_labels[col_i] = ""+iteration+"_"+toString(comp_tolerances[col_i]);

			// Do the optimization for the composite PAZ
			ranks = optimizeComp(comp_tolerances);
			setRow("Comp-Ranks", comp_labels, ranks, imgdone);
			close("*");
			imgdone+=1;
			}
		}
	print("[Experiment_Log]", "\n");
	print("[Experiment_Log]", "\nIteration "+iteration+" results:");
	
	// Average preliminary ranks for tolerances across images
	// For single channel images
	c_bests = newArray(CHANS.length);
	for(c=1; c<=CHANS.length; c++)
		{
		// Get column headings (file then all tolerances)
		colnames=split(Table.headings(CHANS[c-1]+"-Ranks"));

		// Only look at the current iteration
		colnames = Array.slice(colnames, colnames.length - tolerances.length, colnames.length);
		tolerances = newArray(colnames.length);
		avg_score = newArray(tolerances.length);
		for(col_i = 0; col_i<tolerances.length; col_i++)
			{
			tolerances[col_i] = substring(colnames[col_i], 2);
			temp_col = Table.getColumn(colnames[col_i], CHANS[c-1]+"-Ranks");
			Array.getStatistics(temp_col, min, max, mean, stdDev);
			avg_score[col_i] = mean;	
			}
		temp_ranks = Array.rankPositions(avg_score);
		
		if(RANK==true)
			best_index = 0;
		else best_index = temp_ranks.length-1;
		
		index_of_max = temp_ranks[best_index];
		c_bests[c-1] = tolerances[index_of_max];
		FACTOR = 1;
		if(c==AZ)
			FACTOR=10;
		print("[Experiment_Log]", "\n"+CHANS[c-1]+": tolerance: "+parseInt(tolerances[index_of_max])*FACTOR); 
		}

	// Average preliminary ranks for tolerances across images
	// For composite images

	// Get column headings (file then all tolerances)
	colnames=split(Table.headings("Comp-Ranks"));

	// Only look at the current iteration
	colnames = Array.slice(colnames, colnames.length - comp_tolerances.length, colnames.length);
	avg_score = newArray(comp_tolerances.length);
	for(col_i = 0; col_i<comp_tolerances.length; col_i++)
		{
		temp_col = Table.getColumn(colnames[col_i], "Comp-Ranks");
		Array.getStatistics(temp_col, min, max, mean, stdDev);
		avg_score[col_i] = mean;	
		}
	temp_ranks = Array.rankPositions(avg_score);
	
	if(RANK==true)
		best_index = 0;
	else best_index = temp_ranks.length-1;
	
	index_of_max = temp_ranks[best_index];
	comp_best = comp_tolerances[index_of_max];
	
	print("[Experiment_Log]", "\nComp tolerance: "+comp_tolerances[index_of_max]); 
	
	
	// Organize and save the current results for the whole range of tolerances
	// ***********************************************************************
	Table.create("Optimized_Tolerances");
	// Get column headings (file then all tolerances)
	colnames=split(Table.headings("Comp-Ranks"));
	colnames = Array.slice(colnames, 1, colnames.length);
	avg_score = newArray(colnames.length);
	for(col_i = 0; col_i<colnames.length; col_i++)
		{
		temp_col = Table.getColumn(colnames[col_i], "Comp-Ranks");
		Array.getStatistics(temp_col, min, max, mean, stdDev);
		avg_score[col_i] = mean;	
		}

	colnames = reSort(avg_score, colnames);
	avg_score = reSort(avg_score, avg_score);

	Table.setColumn("Comp_tol", colnames, "Optimized_Tolerances");
	Table.setColumn("Comp_ranks", avg_score, "Optimized_Tolerances");

	for(c=0; c<CHANS.length; c++)
		{
		// Get column headings (file then all tolerances)
		colnames=split(Table.headings(""+CHANS[c]+"-Ranks"));
		colnames = Array.slice(colnames, 1, colnames.length);
		avg_score = newArray(colnames.length);
		for(col_i = 0; col_i<colnames.length; col_i++)
			{
			temp_col = Table.getColumn(colnames[col_i], ""+CHANS[c]+"-Ranks");
			Array.getStatistics(temp_col, min, max, mean, stdDev);
			avg_score[col_i] = mean;	
			}
		colnames = reSort(avg_score, colnames);
		avg_score = reSort(avg_score, avg_score);
		
		Table.setColumn(CHANS[c]+"_tol", colnames, "Optimized_Tolerances");
		Table.setColumn(CHANS[c]+"_ranks", avg_score, "Optimized_Tolerances");
		}
	Table.save(savepath+"\\Optimized_Tolerances.csv", "Optimized_Tolerances");
	}
endTime=getTime();
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
print("[Experiment_Log]", "\nEnd: "+ hour+":"+minute+ "\nRuntime:"+round((endTime-startTime)/60000)+" min");
selectWindow("Experiment_Log");
saveAs("Text", savepath+"\\Experiment_Log.txt");
close("Experiment_Log.txt");
close("Experiment_Log");
for(c=0; c<CHANS.length; c++)
	close(""+CHANS[c]+"-Ranks");
close("Comp-Ranks");
close("log");
close("*");


	
function optimizeSC(c, tolerances)
	{
	
	scores = newArray(tolerances.length);
	preprocessSCmesh("min_max", "mask", c, seg_sigma, logr);
	
	if(DEBUG == true)
		print("[DEBUG_LOG]", "Channel: "+CHANS[c-1]+"\n"); 
	for(tol_i=0; tol_i<tolerances.length; tol_i++)
		{
		if(DEBUG == true)
			print("[DEBUG_LOG]", "Starting tolerance: "+tolerances[tol_i]+"   "); 
		pre_rois = roiManager("count");
		// Make single channel mesh images for each individual channel (called CHANS[x]_mesh)
		// params = img, mask, channel, gaussian sigma, tolerance, LoG sigma
		factor=1;
		if(c==AZ)
			factor=10;
		if(DEBUG == true)
			print("[DEBUG_LOG]", "Line291\n"); 
		makeSCmesh("mask", c, tolerances[tol_i]*factor);
	
		if(DEBUG==true)
			print("[DEBUG_LOG]", "Line293\n");
		selectImage(CHANS[c-1]+"_mesh");
		run("Select All");
		run("Analyze Particles...", "size=12-infinity pixel add");
	
		nROI = roiManager("count") - pre_rois;
		if(nROI>0)
			{
			if(nROI>1)
				{
				to_save = getSequence(pre_rois, roiManager("count")-1, 1);
				roiManager("select", to_save);
				roiManager("OR");
				roiManager("add");
				roiManager("select", roiManager("count")-1);
				}
			roiManager("select", roiManager("count")-1);
			roiManager("rename", factor*tolerances[tol_i]+"_all");
			}
		if(DEBUG==true)
			print("[DEBUG_LOG]", "Line312\n");
		// Note - here we've added a composite ROI for all unfiltered mesh units
		// AFTER the individual mesh units. We need to keep track of this while
		// Still analyzing individual rois.
			
		// Exclude from analysis anything at the edge of the NMJ
		if(FILTER==true)
			{
			selectImage("EDM");
			for(roi = pre_rois; roi<roiManager("count")-1; roi++)
				{
				roiManager("select", roi);
				EDM_min = getValue("Min");
				if(EDM_min <= 1)
					{
					roiManager("Delete");
					nROI -= 1;
					roi -=1;
					}
				}
			}

		mesh_scores = newArray(nROI);
		areas = newArray(nROI);
		for(roi=pre_rois; roi<roiManager("count")-1; roi++)
			{
			target = "mean_norm";
			measures = measureROI(target, roi, c);
			areas[roi-pre_rois] = measures[0]; cmean = measures[1]; retmean = measures[2];
			if(c==AZ)
				mesh_scores[roi-pre_rois] = cmean;
			else mesh_scores[roi-pre_rois] = retmean/cmean;

			}
		if(DEBUG==true)
			print("[DEBUG_LOG]", "Line346\n");
		Array.getStatistics(mesh_scores, mr_min, mr_max, mr_mean, mr_stdDev);
		Array.getStatistics(areas, amin, amax, amean, astdDev);
		
		// If segmentation doesn't result in a valid segmentation
		// store penalty value
		if(isNaN(mr_mean))
			mr_mean = -1;

		// penalize segmentation for omitting larger areas of NMJ
		// Give bonus for larger mesh areas
		// This is to address systematic oversegmentation
		// This greatly improves PAZ optimization but not AZ optimization

		root_area = sqrt(amean*areas.length);
		toUnscaled(root_area);
		unscaled_area = pow(root_area, 2);
		area_penalty = 1/(nmj_area - unscaled_area);

		// Don't do area penalty for AZ segmentation
		if(c==AZ)
			amean = 1;

		if(AREABONUS==false)
			amean = 1;
		
		score = pow(mr_mean*area_penalty*area_penalty*amean, 1/4);
		scores[tol_i]=score;			
		close(CHANS[c-1]+"_mesh");

		if(nROI>0)
			{
			// if only 1 roi, need to select second to last
			roiManager("select", roiManager("count")-2);
			
			// If we have ROIs, merge them into one and keep
			if(nROI>1)
				{
				rois=getSequence(pre_rois, roiManager("count")-2, 1);
				roiManager("select", rois);
				roiManager("OR");
				roiManager("add");
				roiManager("delete");
				roiManager("select", roiManager("count")-1);	
				}

			roiManager("rename", tolerances[tol_i]*factor);
			}
		else {
			selectImage(fl[i]);
			run("Select All");
			roiManager("add");
			roiManager("select", roiManager("count")-1);
			roiManager("rename", tolerances[tol_i]*factor);
			}		
		
		}
		
	
	close("forMinima");
	close("toSegment");
	saveROIs(imgpath, ""+iteration+"_"+CHANS[c-1]+"-all_.zip");
	if(roiManager("count")>0)
		roiManager("reset");

	if(DEBUG == true)
		print("[DEBUG_LOG]", "Ending tolerances\n");
	// Replace NaN with a 'bad' score (-1)
	scores = replaceNaN(scores, -1);
	score_ranks = rankArray(scores, false);
	//max = score_ranks[score_ranks.length-1];
	if(RANK==true)
		return score_ranks;
	else return scores;
	}

function optimizeComp(tolerances)
	{
	scores = newArray(tolerances.length);
	makeCompositeImage("min_max", "mask", seg_sigma);
	for(tol_i=0; tol_i<tolerances.length; tol_i++)
		{
		// Make single channel mesh images for each individual channel (called CHANS[x]_mesh)
		// params = img, mask, channel, gaussian sigma, tolerance, LoG sigma
		makeCompositeMesh("min_max", "mask", tolerances[tol_i]);
		pre_rois = roiManager("count"); // number of pre-stored Rois
		selectImage("Composite_mesh");
		run("Select All");
		run("Analyze Particles...", "size=12-infinity pixel add");
		nROI = roiManager("count") - pre_rois;
		
		if(nROI>0)
			{
			if(nROI>1)
				{
				to_save = getSequence(pre_rois, roiManager("count")-1, 1);
				roiManager("select", to_save);
				roiManager("OR");
				roiManager("add");
				}
			roiManager("select", roiManager("count")-1);
			roiManager("rename", tolerances[tol_i]+"_all");
			}
		
		// Note - here we've added a composite ROI for all unfiltered mesh units
		// AFTER the individual mesh units. We need to keep track of this while
		// Still analyzing individual rois.
		
		// Exclude from analysis anything at the edge of the NMJ
		if(FILTER==true)
			{
			selectImage("EDM");
			for(roi = pre_rois; roi<roiManager("count")-1; roi++)  
				{
				roiManager("select", roi);
				EDM_min = getValue("Min");
				if(EDM_min <= 1)
					{
					roiManager("Delete");
					nROI -= 1;
					roi -= 1;
					}
				}
			}
		mesh_scores = newArray(nROI);
		areas = newArray(nROI);

		for(roi=pre_rois; roi<roiManager("count")-1; roi++)
			{
			cmean = 0;
			retmean = 0;
			for(meas_c=1; meas_c<=CHANS.length; meas_c++)
				{
				if(meas_c != AZ)
					{
					measures = measureROI("mean_norm", roi, meas_c);
					areas[roi-pre_rois] += measures[0]; cmean += measures[1]; retmean += measures[2];
					}
				}
			mesh_scores[roi-pre_rois] = retmean/cmean;
			}
		Array.getStatistics(mesh_scores, mr_min, mr_max, mr_mean, mr_stdDev);
		Array.getStatistics(areas, amin, amax, amean, astdDev);
		// If segmentation doesn't result in a valid segmentation
		// store penalty value
		if(isNaN(mr_mean))
			mr_mean = -1;

		// penalize segmentation for omitting larger areas of NMJ
		// Give bonus for larger mesh areas
		// This is to address systematic oversegmentation of PAZ.
		// This greatly improves PAZ optimization but hurts AZ optimization
		root_area = sqrt(amean*areas.length);
		toUnscaled(root_area);
		unscaled_area = pow(root_area, 2);
		spatial_penalty = 1/(nmj_area - unscaled_area);

		if(AREABONUS==false)
			amean = 1;
		
		score = pow(mr_mean*amean*spatial_penalty*spatial_penalty, 1/4);
		scores[tol_i]=score;
		close("Composite_mesh");

		if(nROI>0)
			{
			// If we have ROIs, merge them into one and keep
			if(nROI>1)
				{
				rois=getSequence(pre_rois, roiManager("count")-2, 1);
				roiManager("select", rois);
				roiManager("OR");
				roiManager("add");
				roiManager("delete");
				}
			roiManager("select", roiManager("count")-1);
			roiManager("rename", tolerances[tol_i]);
			}
		else {
			selectImage(fl[i]);
			run("Select All");
			roiManager("add");
			roiManager("select", roiManager("count")-1);
			roiManager("rename", tolerances[tol_i]);
			}
		}
	
	close("forMinima");
	close("toSegment");	
	
	saveROIs(imgpath, ""+iteration+"_Comp-all_.zip");
	if(roiManager("count")>0)
		roiManager("reset");

	// Replace NaN with 'bad' score so we can calculate average and rank
	scores = replaceNaN(scores, -1);
	score_ranks = rankArray(scores, false);
	//max = score_ranks[score_ranks.length-1];
	if(RANK==true)
		return score_ranks;
	else return scores;
	}
	
function measureROI(img, roi, chan)
	{
	// *** Get core and mesh intensities for each channel
	// Shrink mesh ROI by half tracewidth and get pixel intensity data for core region
	total_rois = roiManager("count");
	selectImage(img);	
	roiManager("select", roi);
	getStatistics(carea, cmean, cmin, cmax, cstd);
	run("Enlarge...", "enlarge=-"+TRACEWIDTH/2+" pixel");
	if(CHANS.length>1)
		Stack.setChannel(chan);
	getStatistics(carea2, cmean, cmin, cmax, cstd);
	
	// For very small meshes, ROI may not shrink. If so, try smaller shrink factor.
	if(carea2==carea)
		run("Enlarge...", "enlarge=-"+TRACEWIDTH/4+" pixel");
		
	getStatistics(carea2, cmean, cmin, cmax, cstd);

	// create a mesh ROI as xor of ROI shrunken and expanded by tracewidth/2
	roiManager("select", roi);
	run("Enlarge...", "enlarge=-2 pixel");
	roiManager("add");
	run("Enlarge...", "enlarge=2 pixel");
	roiManager("add");
	roiManager("select", newArray(total_rois, total_rois+1));
	roiManager("XOR");
	roiManager("add");
	roiManager("delete");
	// Now we have original rois intact but with new 'mesh' roi at position total_rois. 
	// We will remove this after measurements for this roi are done.

	// Measure basic intensity values
	selectImage(img);
	roiManager("select", total_rois);
	if(CHANS.length>1)
		Stack.setChannel(chan);
	getStatistics(retarea, retmean, retmin, retmax, retstd);

	// delete added mesh roi to restore original ROI list
	roiManager("select", total_rois);
	roiManager("delete");

	if(img != "mean_norm")
		close(img);
	return newArray(carea, cmean, retmean);
	}
	
function segment(imp)
	{
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
	rename("full_mask");
	run("Duplicate...", "title=mask");
	for(e=0; e<ERODE; e++)
		run("Erode");

	// Duplicate mask and make into a distance map
	// For downstream spatial analysis and filtering out of edge mesh units
	run("Duplicate...", "title=EDM duplicate"); 
	run("Distance Map");
	}


function preprocessSCmesh(img, mask, chan, sigma, LoGr)
	{
	// Do preprocessing of images for subsequent segmentation 
	// create image for segmentation 'toSegment'
	// and image for detection of minima 'forMinima'
	
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
	   	// testing whether full mask better for finding AZ
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
	}

function makeSCmesh(mask, chan, tolerance)
	{
	// To do: AZ object segmentation refinement
	
	// Generate an image containing a mesh-like representation of a single input channel
	// Return an array of (x0, ..., xn, y0, ..., yn) describing the x,y centroid 
	// locations of n mesh units. 
	total_rois = roiManager("count");
	selectImage("forMinima");
	run("Find Maxima...", "prominence="+tolerance+" exclude light output=[Single Points]");
	rename("BinMark");
		// If channel is AZ, generate a mask for all AZ objects
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
		roiManager("save", imgpath+"\\AZ_objects_"+tolerance+".roi");
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
	run("Analyze Particles...", "add");
	x1 = newArray(roiManager("count"));
	y1 = newArray(roiManager("count"));
	for(r=total_rois; r<roiManager("count"); r++)
		{
		roiManager("Select", r);
		x = getValue("X"); 
		y = getValue("Y");
		toUnscaled(x,y);
		x1[r] = x;
		y1[r] = y;
		}
	
	rename(""+CHANS[chan-1]+"_mesh");
	centroids=Array.concat(x1,y1);
	if(roiManager("count")>total_rois)
		{
		rois = getSequence(total_rois, roiManager("count")-1, 1);
		roiManager("select", rois);
		roiManager("delete");
		}
	close("BinMark");
	close("AZ_objects");
	return centroids;
	}
	
function makeCompositeImage(img, mask, sigma)
	{
	// Generate an image containing a mesh-like representation of a composite of all channels
	// Return an array of (x0, ..., xn, y0, ..., yn) describing the x,y centroid 
	// locations of n mesh units. 
	total_rois = roiManager("count");
	selectImage(img);
	getDimensions(width, height, channels, slices, frames);
	run("Duplicate...", "title=Normed duplicate channels=1-"+channels);
	run("Gaussian Blur...", "sigma="+seg_sigma+" stack");
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
	
	// Apply 'internal gradient' filter with octagon element radius 4
	// This sharpens PAZ edges for detection of local minima. Called 'forMinima'
	selectImage("toSegment");
	run("16-bit");
	run("Morphological Filters", "operation=[Internal Gradient] element=Octagon radius="+IG_RADIUS+1);
	rename("forMinima");

    if(POST_IG>0)
    	run("Gaussian Blur...", "sigma="+POST_IG+" stack");

	}

function makeCompositeMesh(img, mask,  tolerance)
	{
	// Find local minima in composite forMinima image
	// Generate mesh representation within masked region using marker-controlled watershed
	total_rois = roiManager("count");
	selectImage("forMinima");
	run("Find Maxima...", "prominence="+tolerance+" exclude light output=[Single Points]");
	rename("BinMark");
	run("Marker-controlled Watershed", "input=toSegment marker=BinMark mask="+mask+" binary calculate use");
	rename("Comp_Mesh_Regions");
	
	// Create binary mesh representation and save centroid coordinates
	setThreshold(1, 66000);
	run("Convert to Mask");
	run("Analyze Particles...", "add");
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
	if(roiManager("count")>total_rois)
		{
		rois = getSequence(total_rois, roiManager("count")-1, 1);
		roiManager("select", rois);
		roiManager("delete");
		}
	centroids=Array.concat(x1,y1);
	return centroids;
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

function rankArray(array, lowhi)
	{
	if(lowhi==false)
		value = -1;
	else value = 100000000000000;	
	to_rank = replaceNaN(array, value);
	ranks=Array.rankPositions(to_rank);
	if(lowhi==false)
		ranks=Array.reverse(ranks);
	rankOrder=newArray(ranks.length);
	for(i=0; i<ranks.length; i++)
		{
		r=0;
		while(ranks[r]!=i)
			{r+=1;}
		rankOrder[i]=r;
		}
	for(i=0; i<array.length; i++)
		{
		// handle equal values/assign tied ranks
		found=ArrayFind(array, array[i]);
		if(found.length>1)
			{
			ranktotal=0;
			for(r=0; r<found.length; r++)
				ranktotal+=rankOrder[found[r]];
			rankavg=ranktotal/found.length;
			for(r=0; r<found.length; r++)
				rankOrder[found[r]]=rankavg;
			}
		}
	return rankOrder;
	}

function ArrayFind(array, search)
	{
	found=newArray(0);
	for(i=0; i<array.length; i++)
		{
		if(array[i]==search)
			{
			a=newArray(1); 
			a[0]=i;
			found=Array.concat(found, a);
			}
		}
	return found;
	}
	
function replaceNaN(array, value)
	{
	// Replace all NaN in array with value
	new=newArray(array.length);
	for(i=0; i<array.length; i++)
		{
		if(isNaN(array[i])==true)
			new[i] = value;
		else new[i] = array[i];
		}
	return new;
	}

function AZByIntensity(img)
	{
	// Generate ROIs from a grayscale (instance) labeled image
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
				if(brp_mean<1)
					roiManager("delete");
				}
			}
		}
	}
	
function reSort(guide, a)
	{
	// Sorts a by rank of guide
	rank=Array.rankPositions(guide);
	new=newArray(a.length);
	for(z=0; z<a.length; z++)
		{
		new[z]=a[rank[z]];
		}
	return new;
	}
	
function watershed_az(img)
	{
	// Take as input a grayscale object map output
	// from seeded region growing. Convert objects to binary
	// Then do distance transform watershed. 
	// Produce new object and centroid map.
	start_rois = roiManager("count");
	AZByIntensity(img);

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
	
function saveROIs(path, imp)
	{
	// Save all ROIs in Roi list
	rois=newArray(roiManager("count"));
	for (n=0; n<rois.length; n++)
		{
		rois[n]=n;
		}
	roiManager("select", rois);
	savename = path+"\\"+stripX(imp)+"ROIs.zip";
	roiManager("Save", savename);
	}

function setRow(table, colnames, values, row)
	{
	for(col=0; col<colnames.length; col++)
		{
		Table.set(colnames[col], row, values[col], table);
		}
	}
	
function spread(num, len, space)
	{
	spread_array = newArray(len);
	start_factor = 1 - floor((len-1)/2)*space;
	for(i=0; i<len; i++)
		spread_array[i] = (start_factor + i*space) * num;
	return spread_array;
	}

function test()
	{
	tolerances = newArray(100, 500, 1000, 2000, 3000, 4000, 5000, 10000, 15000);
	c_scores = newArray(1.705, 1.696, 1.704, 1.664, 1.633, 1.633, 1.612, 1.494, 1.475);
	Fit.doFit("3rd Degree Polynomial", tolerances, c_scores);
	Fit.plot();
	deriv = newArray(Fit.nParams()-1);
	for(p=1; p<Fit.nParams(); p++)
		deriv[p-1]=Fit.p(p);
	Fit.plot();
	}

function stripX(string)
	{
	// This is because Macro language doesn't have a general use name without extension
	return substring(string, 0, lastIndexOf(string, "."));
	}

