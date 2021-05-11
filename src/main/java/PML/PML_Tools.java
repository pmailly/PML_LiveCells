package PML;


import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Concatenator;
import ij.plugin.GaussianBlur3D;
import ij.plugin.ZProjector;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.FastFilters3D;
import mcib3d.image3d.regionGrowing.Watershed3D;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import loci.formats.meta.IMetadata;
import mpicbg.ij.integral.RemoveOutliers;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;
        
 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose dots_Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
public class PML_Tools {
   
    // min nucleus volume in mic^3
    public double minNuc = 500;
    // max nucleus volume in micron^3
    public double maxNuc = 5000;
    
    // min volume in pixels^3 for dots
    public double minPML = 0.05;
    // max volume in pixels^3 for dots
    public double maxPML = 60;
    private Calibration cal = new Calibration(); 
    
    // dots dog paramaters
    private int sig1 = 2;
    private int sig2 = 4;
    // dots threshold method
    private String thMet = "Moments";
    // pml dot dilatation factor for diffuse mask
    private float dilate = 0.5f;
    
    // Trackmate dialog parameters
    public double radius = 0.75;
    public double threshold = 35;
    
    public CLIJ2 clij2 = CLIJ2.getInstance();
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
        
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal = new Calibration();  
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return(cal);
    }
    
    
    /**
     *
     * @param img
     */
    public void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    public Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
   
    
    
     /**  
     * median 3D box filter
     * Using CLIJ2
     * @param imgCL
     * @param sizeX
     * @param sizeY
     * @param sizeZ
     * @return imgOut
     */ 
    public ClearCLBuffer medianFilter(ClearCLBuffer imgCL, double sizeX, double sizeY, double sizeZ) {
        ClearCLBuffer imgIn = clij2.push(imgCL);
        ClearCLBuffer imgOut = clij2.create(imgIn);
        clij2.median3DBox(imgIn, imgOut, sizeX, sizeY, sizeZ);
        clij2.release(imgCL);
        return(imgOut);
    }
    
    /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param imgCL
     * @param size1
     * @param size2
     * @return imgGauss
     */ 
    public ClearCLBuffer DOG(ClearCLBuffer imgCL, double size1, double size2) {
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, size1, size1, size1, size2, size2, size2);
        clij2.release(imgCL);
        return(imgCLDOG);
    }
    
    /**
     * Fill hole
     * USING CLIJ2
     */
    private void fillHole(ClearCLBuffer imgCL) {
        long[] dims = clij2.getSize(imgCL);
        ClearCLBuffer slice = clij2.create(dims[0], dims[1]);
        ClearCLBuffer slice_filled = clij2.create(slice);
        for (int z = 0; z < dims[2]; z++) {
            clij2.copySlice(imgCL, slice, z);
            clij2.binaryFillHoles(slice, slice_filled);
            clij2.copySlice(slice_filled, imgCL, z);
        }
        clij2.release(slice);
        clij2.release(slice_filled);
    }
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param imgCL
     * @param thMed
     * @param fill 
     */
    public ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed, boolean fill) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        if (fill)
            fillHole(imgCLBin);
        return(imgCLBin);
    }
        
    /**
     * Dialog 
     * 
     * @return 
     */
    public String dialog() {
        
        String dir = "";
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addDirectoryField("Choose Directory Containing Image Files : ", "");
        gd.addNumericField("Min PML size (µm3) : ", minPML, 3);
        gd.addNumericField("Max PML size (µm3) : ", maxPML, 3);
        gd.addMessage("Diffuse analyze", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("PML dilatation factor (µm) :", dilate, 3);
        gd.addMessage("Trackmate parameters", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("PML dots radius (µm) :", radius, 3);
        gd.addNumericField("PML threshold        :", threshold, 3);
        gd.addImage(icon);
        gd.showDialog();
        dir = gd.getNextString()+File.separator;
        minPML = gd.getNextNumber();
        maxPML = gd.getNextNumber();
        dilate = (float)gd.getNextNumber();
        radius = gd.getNextNumber();
        threshold = gd.getNextNumber();
        return(dir);
    }
    
    /**
     * Remove Outliers
     * 
     * @param img
     * @param radX
     * @param radY
     * @param factor
     * @return img
     */
    public ImagePlus removeOutliers(ImagePlus img, int radX, int radY, float factor) {
        
        for (int i = 0; i < img.getNSlices(); i++) {
            img.setSlice(i);
            ImageProcessor ip = img.getProcessor();
            RemoveOutliers removeOut = new RemoveOutliers(ip.convertToFloatProcessor());
            removeOut.removeOutliers(radX, radY, factor);
        }
        return(img);
    }
    
    
/**
     * Nucleus segmentation
     * @param imgNuc
     * @return 
     */
    public Object3D findnucleus(ImagePlus imgNuc) {
        removeOutliers(imgNuc, 20, 20, 1);
        ImageStack stack = new ImageStack(imgNuc.getWidth(), imgNuc.getHeight());
        for (int i = 1; i <= imgNuc.getStackSize(); i++) {
            IJ.showStatus("Finding nucleus section "+i+" / "+imgNuc.getStackSize());
            imgNuc.setZ(i);
            imgNuc.updateAndDraw();
            IJ.run(imgNuc, "Nuclei Outline", "blur=20 blur2=30 threshold_method=Li outlier_radius=0 outlier_threshold=0 max_nucleus_size=500"
                    + " min_nucleus_size=10 erosion=0 expansion_inner=0 expansion=0 results_overlay");
            imgNuc.setZ(1);
            imgNuc.updateAndDraw();
            ImagePlus mask = new ImagePlus("mask", imgNuc.createRoiMask().getBufferedImage());
            ImageProcessor ip =  mask.getProcessor();
            ip.invertLut();
            stack.addSlice(ip);
        }
        ImagePlus imgStack = new ImagePlus("Nucleus", stack);
        imgStack.setCalibration(imgNuc.getCalibration());
//        imgStack.show();
//        new WaitForUserDialog("test").show();
        Objects3DPopulation nucPop = new Objects3DPopulation(getPopFromImage(imgStack).getObjectsWithinVolume(minNuc, maxNuc, true));
        nucPop.removeObjectsTouchingBorders(imgStack, false);
        //nucPop.updateNamesAndValues();
        Object3D nucObj = nucPop.getObject(0);
        closeImages(imgStack);
        return(nucObj);
    }
    
    
    // Threshold images and fill holes
    public void threshold(ImagePlus img, AutoThresholder.Method thMed, boolean fill) {
        //  Threshold and binarize
        img.setZ(img.getNSlices()/2);
        img.updateAndDraw();
        IJ.setAutoThreshold(img, thMed.toString()+" dark");
        Prefs.blackBackground = false;
        IJ.run(img, "Convert to Mask", "method="+thMed.toString()+" background=Dark");
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    
    /**
     *Z project
     * @param img
     * @return 
     */
    public ImagePlus stackProj(ImagePlus img) {
       ZProjector proj = new ZProjector(img);
       ImagePlus imgProj = proj.run(img, "max all");
       return(imgProj); 
    }
    
    /** 
     * Find dots
     * @param img channel
     * @return dots population
     */
    public Objects3DPopulation findDots(ImagePlus img, Object3D nucObj) {
        ClearCLBuffer imgCL = clij2.push(img);
        //ClearCLBuffer imgCLMed = clij2.create(imgCL);
        //clij2.mean3DBox(imgCL, imgCLMed, 1, 1, 1);
        //clij2.copy(imgCLMed, imgCL);
        //clij2.release(imgCLMed);
        //clij2.release(imgCL);
        ClearCLBuffer imgCLDOG = DOG(imgCL, sig1, sig2);
        clij2.release(imgCL);
        ImagePlus imgBin = clij2.pull(threshold(imgCLDOG, thMet, false));
        clij2.release(imgCLDOG);
        imgBin.setCalibration(img.getCalibration());
        //imgBin.show();
        //new WaitForUserDialog("test").show();
        Objects3DPopulation pmlPop = new Objects3DPopulation(getPopFromImage(imgBin).getObjectsWithinVolume(minPML, maxPML, true));
        for ( int i = 0; i < pmlPop.getNbObjects(); i++)
        {
            Object3D obj = pmlPop.getObject(i);
            if ( !obj.hasOneVoxelColoc(nucObj) )
            {
                pmlPop.removeObject(i);
                i--;
            }
            
        }
        closeImages(imgBin);
      
        return(pmlPop);
    } 
    
    
    public Objects3DPopulation getPopFromImage(ImagePlus img, Calibration cal) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
    /**
     * Read diffus PML intensity
     * fill PML voxel with zero in PML channel
     * Add to nucObj comment mean intensity pre-processing(integrated = false)
     * or integrated intensity (integrated = true)
     * 
     * @param pmlPop
     * @param nucObj
     * @param imgPML
     * @param mean
     * @return 
     */
    public double pmlDiffus(Objects3DPopulation pmlPop, Object3D nucObj, ImagePlus imgPML) {
        ImageHandler imhDotsDiffuse = ImageHandler.wrap(imgPML.duplicate());
        double pmlIntDiffuse ;
        int pmlDots = pmlPop.getNbObjects();
        double volPMLDilated = 0;
        float dilate = 0.5f;
        for (int p = 0; p < pmlDots; p++) {
            Object3D pmlObj = pmlPop.getObject(p);
            // dilate 
            Object3DVoxels pmlDilatedObj = pmlObj.getDilatedObject(dilate, dilate, dilate);
            pmlDilatedObj.draw(imhDotsDiffuse, 0);
            volPMLDilated += pmlDilatedObj.getVolumeUnit();
        }
        double nucVolume = nucObj.getVolumeUnit();
        pmlIntDiffuse = nucObj.getIntegratedDensity(imhDotsDiffuse); 
        imhDotsDiffuse.closeImagePlus();
        return(pmlIntDiffuse);
    }
    
    
    /**
     * Save image objects
     */
    public ImagePlus saveImageObjects(ArrayList<Objects3DPopulation> pmlPopList, Objects3DPopulation nucPop, ImagePlus[] imgArray, String pathName) {
        ImagePlus[] hyperBin = new ImagePlus[imgArray.length];
       for (int i=0; i<imgArray.length; i++)
       {
           ImagePlus img = imgArray[i];
           
           // image PML at time point i
           ImageHandler imhObjects = ImageHandler.wrap(img).createSameDimensions();
            nucPop.getObject(i).draw(imhObjects, 64);
            Objects3DPopulation pmlPop = pmlPopList.get(i);
            for (int o = 0; o < pmlPop.getNbObjects(); o++) 
            {
                Object3D pmlObj = pmlPop.getObject(o);
                pmlObj.draw(imhObjects, 255);
            // labelsObject(pmlObj, imhObjects.getImagePlus(), o, 255);
            }
            imhObjects.getImagePlus().setCalibration(cal);
            imhObjects.getImagePlus().setSlice(imhObjects.getImagePlus().getNSlices()/2);
            IJ.resetMinAndMax(imhObjects.getImagePlus());
            hyperBin[i] = imhObjects.getImagePlus();
        }
       
       ImagePlus hyperRes = new Concatenator().concatenate(hyperBin, false);
       //hyperRes.show();
        // save image for objects population
        FileSaver ImgObjectsFile = new FileSaver(hyperRes);
        ImgObjectsFile.saveAsTiff(pathName);
        return hyperRes;
    }
    
     /**
     * Save image objects
     */
    public ImagePlus saveImagePMLs(Objects3DPopulation nucPop, ImagePlus[] imgArray, String pathName) {
       int[] bb =  {10000,0,10000,0,10000,0};
       for (int i=0; i<imgArray.length; i++)
       {
           Object3D nuc = nucPop.getObject(i);
           int[] boundingBox = nuc.getBoundingBox();
           for ( int j = 0; j < 6; j+=2) 
           {
               if ( boundingBox[j] < bb[j]) bb[j] = boundingBox[j];
           }
           for ( int j = 1; j < 6; j+=2) 
           {
               if ( boundingBox[j] > bb[j]) bb[j] = boundingBox[j];
           }
        }
        ImagePlus hyperPML = new Concatenator().concatenate(imgArray, false);
        Roi roi = new Roi(bb[0], bb[2], bb[1]-bb[0], bb[3]-bb[2]);
        hyperPML.setRoi(roi);       
       //hyperRes.show();
        // save image for objects population
        FileSaver ImgObjectsFile = new FileSaver(hyperPML);
        ImgObjectsFile.saveAsTiff(pathName);
        return hyperPML;
    }
    
 

    /**
     * Create diffuse PML image
     * Fill zero in pml dots
     * @param pmlPopList
     * @param imgArray
     * @param pathName

     */
    public void saveDiffusImage(ArrayList<Objects3DPopulation> pmlPopList, ImagePlus[] imgArray, String pathName) {
        ImagePlus[] hyperDifuse = new ImagePlus[imgArray.length];
        for (int i = 0; i < imgArray.length; i++) {
           ImageHandler imh = ImageHandler.wrap(imgArray[i]);
           Objects3DPopulation pmlPop = pmlPopList.get(i);
           for (int j = 0; j < pmlPop.getNbObjects(); j++) {
               Object3D pmlObj = pmlPop.getObject(j);
               // dilate 
                Object3DVoxels pmlDilatedObj = pmlObj.getDilatedObject(dilate, dilate, dilate);
                pmlDilatedObj.draw(imh, 0);
            }
           hyperDifuse[i] = imh.getImagePlus(); 
        }
        ImagePlus hyperDifuseTime = new Concatenator().concatenate(hyperDifuse, true);
        // Save diffus
        FileSaver imgDiffus = new FileSaver(hyperDifuseTime);
        imgDiffus.saveAsTiff(pathName);
        closeImages(hyperDifuseTime);
    }
           
    
    public ImagePlus WatershedSplit(ImagePlus binaryMask, float rad) {
        float resXY = 1;
        float resZ = 1;
        float radXY = rad;
        float radZ = rad;
        Calibration cal = binaryMask.getCalibration();
        if (cal != null) {
            resXY = (float) cal.pixelWidth;
            resZ = (float) cal.pixelDepth;
            radZ = radXY * (resXY / resZ);
        }
        ImageInt imgMask = ImageInt.wrap(binaryMask);
        ImageFloat edt = EDT.run(imgMask, 0, resXY, resZ, false, 0);
        ImageHandler edt16 = edt.convertToShort(true);
        ImagePlus edt16Plus = edt16.getImagePlus();
        GaussianBlur3D.blur(edt16Plus, 2.0, 2.0, 2.0);
        edt16 = ImageInt.wrap(edt16Plus);
        edt16.intersectMask(imgMask);
        // seeds
        ImageHandler seedsImg = FastFilters3D.filterImage(edt16, FastFilters3D.MAXLOCAL, radXY, radXY, radZ, 0, false);
        Watershed3D water = new Watershed3D(edt16, seedsImg, 0, 0);
        water.setLabelSeeds(true);
        return(water.getWatershedImage3D().getImagePlus());
    }
    
    
    /**
     * Get total pml volume
     * @param pmlPop
     * @return 
     */
    public DescriptiveStatistics getPMLVolume(Objects3DPopulation pmlPop) {
        DescriptiveStatistics pmlVolume = new DescriptiveStatistics();
        for (int i = 0; i < pmlPop.getNbObjects(); i++) {
            Object3D pmlObj = pmlPop.getObject(i);
            pmlVolume.addValue(pmlObj.getVolumeUnit());
        }
        return(pmlVolume);
    }
    

    /**
     * Get pml intensity
     * @param pmlPop
     * @param img
     * @return 
     */
    public DescriptiveStatistics getPMLIntensity(Objects3DPopulation pmlPop, ImagePlus img) {
        DescriptiveStatistics pmlInt = new DescriptiveStatistics();
        ImageHandler imh = ImageHandler.wrap(img);
        for (int i = 0; i < pmlPop.getNbObjects(); i++) {
            Object3D pmlObj = pmlPop.getObject(i);
            pmlInt.addValue(pmlObj.getIntegratedDensity(imh));
        }
        imh.closeImagePlus();
        return(pmlInt);
    }
}