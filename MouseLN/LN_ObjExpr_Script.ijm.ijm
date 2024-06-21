open("C:/Users/schromec/Desktop/SPACE/Figures/Fig2_OrigImage_SegMask.tif");
open("C:/Users/schromec/Desktop/SPACE/Figures/Fig2_OrigImage.tif");
run("Stack to Images");
for(i=1; i<43; i++) {
	if(i<10) {
		arg = "input=Fig2_OrigImage-000" + i + " labels=Fig2_OrigImage_SegMask.tif mean";
		file = "C:/Users/schromec/Desktop/SPACE/Figures/Fig2_ObjExpr/Ch0" + i + ".csv";
	} else {
		arg = "input=Fig2_OrigImage-00"  + i + " labels=Fig2_OrigImage_SegMask.tif mean";
		file = "C:/Users/schromec/Desktop/SPACE/Figures/Fig2_ObjExpr/Ch"  + i + ".csv";
	}
	run("Intensity Measurements 2D/3D", arg);
	saveAs("Results", file);
}
selectImage("Fig2_OrigImage_SegMask.tif");
run("Analyze Regions", "centroid");
saveAs("Results", "C:/Users/schromec/Desktop/SPACE/Figures/Fig2_ObjExpr/Centroids.csv");
run("Close All");