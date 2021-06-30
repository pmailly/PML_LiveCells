package PML;


import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.GaussianBlur3D;
import ij.plugin.ImageCalculator;
import ij.plugin.SubHyperstackMaker;
import ij.plugin.Thresholder;
import ij.plugin.ZProjector;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
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
import loci.plugins.util.ImageProcessorReader;
import mpicbg.ij.integral.RemoveOutliers;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

import StardistPML.StarDist2D;
import ij.plugin.RoiEnlarger;
import ij.plugin.RoiScaler;
import ij.plugin.frame.RoiManager;
import java.util.HashMap;
import net.imagej.Dataset;
import net.imagej.ImageJ;

        
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
    public double minPML = 0.02;
    // max volume in pixels^3 for dots
    public double maxPML = 60;
    private Calibration cal = new Calibration(); 
    

    // dots threshold method
    private String thMet = "Otsu";
    // pml dot dilatation factor for diffuse mask
    private float dilate = 0.5f;
    
    // Trackmate dialog parameters
    public double radius = 0.75;
    public double threshold = 35;
    public double merging_dist = 0.8;
    public String trackMate_Detector_Method = "LoG";
    
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
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("uk.ac.sussex.gdsc.utils.DifferenceOfGaussians_PlugIn");
        } catch (ClassNotFoundException e) {
            IJ.log("GDSC Suite not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("TurboReg_");
        } catch (ClassNotFoundException e) {
            IJ.log("TurboReg not installed, please install from http://bigwww.epfl.ch/thevenaz/turboreg/");
            return false;
        }
        return true;
    }
    
    /**
     * Find images in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
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
     * Find channels name
     * @param imageName
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public static String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader, boolean bioformat) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                        if (!bioformat) {
                            channels[n] = channels[n].replace("_", "-");
                            channels[n] = "w"+(n+1)+channels[n];
                        }
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            default :
                for (int n = 0; n < chs; n++)
                    channels[0] = Integer.toString(n);
        }
        return(channels);         
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
    
    public Calibration getCalib()
    {
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return cal;
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
    
   
    
    
    /* Median filter 
     * Using CLIJ2
     * @param ClearCLBuffer
     * @param sizeXY
     * @param sizeZ
     */ 
    public ClearCLBuffer median_filter(ClearCLBuffer  imgCL, double sizeXY, double sizeZ) {
        ClearCLBuffer imgCLMed = clij2.create(imgCL);
        clij2.mean3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
        clij2.release(imgCL);
        return(imgCLMed);
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
     * Find image type
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                    return fileExt;
                case "czi" :
                  return fileExt;
                case "lif"  :
                    return fileExt;
                case "isc2" :
                   return fileExt;
                default :
                   ext = fileExt;
                   break; 
            }
        }
        return(ext);
    }
    
    
    /**
     * Dialog 
     * 
     * @param channels
     * @return 
     */
    public int[] dialog(String[] channels) {
        String[] thMethods = new Thresholder().methods;
        String[] TrackMate_Detector = {"DoG", "LoG"};
        String[] chNames = {"Nucleus", "PML"};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String chName : channels) {
            gd.addChoice(chNames[index]+" : ", channels, channels[index]);
            index++;
        }
        gd.addMessage("PML parameters", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min PML size (µm3) : ", minPML, 3);
        gd.addNumericField("Max PML size (µm3) : ", maxPML, 3);
        gd.addMessage("Diffuse analyze", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("Dots threshold method : ", thMethods, thMet);
        gd.addNumericField("PML dilatation factor (µm) :", dilate, 3);
        gd.addMessage("Trackmate parameters", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("Dots detector method :", TrackMate_Detector, TrackMate_Detector[1]);
        gd.addNumericField("Merging / spliting max distance : ", merging_dist,3);
        gd.addNumericField("PML dots radius (µm) :", radius, 3);
        gd.addNumericField("PML threshold        :", threshold, 3);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Calibration xy (µm)  :", cal.pixelWidth, 3);
        if ( cal.pixelDepth == 1) cal.pixelDepth = 0.5;
        gd.addNumericField("Calibration z (µm)  :", cal.pixelDepth, 3);
        gd.showDialog();
        int[] chChoices = new int[channels.length];
        for (int n = 0; n < chChoices.length; n++) {
            chChoices[n] = ArrayUtils.indexOf(channels, gd.getNextChoice());
        }
        minPML = gd.getNextNumber();
        maxPML = gd.getNextNumber();
        thMet = gd.getNextChoice();
        dilate = (float)gd.getNextNumber();
        trackMate_Detector_Method = gd.getNextChoice();
        merging_dist = gd.getNextNumber();
        radius = gd.getNextNumber();
        threshold = gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelHeight = cal.pixelWidth;
        cal.pixelDepth = gd.getNextNumber();
        if (gd.wasCanceled())
                chChoices = null;
        return(chChoices);
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
        ImagePlus imgRemo = new Duplicator().run(img);
        for (int i = 1; i <= imgRemo.getNSlices(); i++) {
            imgRemo.setSlice(i);
            ImageProcessor ip = imgRemo.getProcessor();
            RemoveOutliers removeOut = new RemoveOutliers(ip.convertToFloatProcessor());
            removeOut.removeOutliers(radX, radY, factor);
        }
        return(imgRemo);
    }
  /**
     * Nucleus segmentation
     * @param imgNuc
     * @return 
     */
    public Object3D findnucleus(ImagePlus imgNuc, int t) {
        
        IJ.run(imgNuc, "Median 3D...", "x=2 y=2 z=2");
        IJ.run(imgNuc, "Remove Outliers", "block_radius_x=40 block_radius_y=40 standard_deviations=1 stack");
        // try to get rid of extreme slices without nuclei
        ImagePlus globalBin = new Duplicator().run(imgNuc);
 
        IJ.setAutoThreshold(globalBin, "Otsu dark stack");
        Prefs.blackBackground = false;
        IJ.run(globalBin, "Convert to Mask", "method=Otsu background=Dark stack");
        ImageStack stack = new ImageStack(imgNuc.getWidth(), imgNuc.getHeight());
        for (int i = 1; i <= imgNuc.getStackSize(); i++) {
            IJ.showStatus("Finding nucleus section "+i+" / "+imgNuc.getStackSize());
            globalBin.setSlice(i);
            ImageStatistics stat = globalBin.getStatistics();
            // contains part of the nucleus
            if (stat.mean > 10)
            {
                imgNuc.setZ(i);
                imgNuc.updateAndDraw();
                IJ.run(imgNuc, "Nuclei Outline", "blur=2 blur2=300 threshold_method=Li outlier_radius=0 outlier_threshold=0 max_nucleus_size=500"
                        + " min_nucleus_size=10 erosion=0 expansion_inner=0 expansion=0 results_overlay");
                imgNuc.setZ(1);
                imgNuc.updateAndDraw();
                ImagePlus mask = new ImagePlus("mask", imgNuc.createRoiMask().getBufferedImage());
                ImageProcessor ip =  mask.getProcessor();
                stack.addSlice(ip);
            }
            // empty slice
            else {
                ImagePlus mask = IJ.createImage("mask", "8-bit black", imgNuc.getWidth(), imgNuc.getHeight(), 1);
                ImageProcessor ip =  mask.getProcessor();
                stack.addSlice(ip);  
            }
        }
        closeImages(globalBin);
        ImagePlus imgStack = new ImagePlus("Nucleus", stack);
        imgStack.setCalibration(imgNuc.getCalibration());
        // Apply transformation and rebinarize
        /**if (t > 0){
            trans.get(t - 1).doTransformation(imgStack);
            IJ.setAutoThreshold(imgStack, "Default dark stack");
            Prefs.blackBackground = false;
            IJ.run(imgStack, "Convert to Mask", "method=Default background=Dark stack");
         }*/
          
        Objects3DPopulation nucPop = new Objects3DPopulation(getPopFromImage(imgStack).getObjectsWithinVolume(minNuc, maxNuc, true));
        // Find bigger object
        int ind = 0;    
        if (nucPop.getNbObjects()>1){
            double maxv = 0;
            for (int k=0; k<nucPop.getNbObjects(); k++) {
                Object3D obj = nucPop.getObject(k);
                double vol = obj.getVolumePixels();
                if (vol>maxv){
                    ind = k;
                    maxv = vol;
                }
            }
            
        }
        if ( nucPop.getNbObjects()==0 ) {
            imgStack.show();
            ind = 0;
            //IJ.Log("No nucleus found, t "+t);
            return null;
        }
        Object3D nucObj = nucPop.getObject(ind);
        closeImages(imgStack);
        return(nucObj);
    } 
    
    
    // Threshold images and fill holes
    public void threshold(ImagePlus img, String thMed, boolean fill) {
        //  Threshold and binarize
        img.setZ(img.getNSlices()/2);
        img.updateAndDraw();
        IJ.setAutoThreshold(img, thMed+" dark");
        Prefs.blackBackground = false;
        IJ.run(img, "Convert to Mask", "method="+thMed+" background=Dark");
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
     * Find dots with LoG
     * @param img channel
     * @param nucObj
     * @return dots population
     */
    public Objects3DPopulation findDotsLoG(ImagePlus img, Object3D nucObj) {
        ImagePlus imgLaplacien = new Duplicator().run(img); 
        double sigma = radius/(3*img.getCalibration().pixelWidth);
        IJ.run(imgLaplacien, "Laplacian of Gaussian", "sigma="+sigma+" negate stack");
        imgLaplacien.setSlice(imgLaplacien.getNSlices()/2);
        IJ.setAutoThreshold(imgLaplacien, thMet+" dark");
        Prefs.blackBackground = false;
        IJ.run(imgLaplacien, "Convert to Mask", "method="+thMet+" background=Default");
        
        Objects3DPopulation pmlPop = new Objects3DPopulation(getPopFromImage(imgLaplacien).getObjectsWithinVolume(minPML, maxPML, true));
        for ( int i = 0; i < pmlPop.getNbObjects(); i++)
        {
            Object3D obj = pmlPop.getObject(i);
            if ( !obj.hasOneVoxelColoc(nucObj) )
            {
                pmlPop.removeObject(i);
                i--;
            }
            
        }
        closeImages(imgLaplacien);
      
        return(pmlPop);
    } 
    
    /** 
     * Find dots with DOG method
     * @param img channel
     * @return dots population
     */
    public Objects3DPopulation findDotsDoG(ImagePlus img, Object3D nucObj) {
        ImagePlus imgDOG = new Duplicator().run(img);
        double sig1 = radius;
        double sig2 = radius/3;
        IJ.run(imgDOG,"Difference of Gaussians", "  sigma1="+sig1+" sigma2="+sig2+" scaled stack");
        imgDOG.setSlice(imgDOG.getNSlices()/2);
        IJ.setAutoThreshold(imgDOG, thMet+" dark");
        Prefs.blackBackground = false;
        IJ.run(imgDOG, "Convert to Mask", "method="+thMet+" background=Default");

        Objects3DPopulation pmlPop = new Objects3DPopulation(getPopFromImage(imgDOG).getObjectsWithinVolume(minPML, maxPML, true));
        for ( int i = 0; i < pmlPop.getNbObjects(); i++)
        {
            Object3D obj = pmlPop.getObject(i);
            if ( !obj.hasOneVoxelColoc(nucObj) )
            {
                pmlPop.removeObject(i);
                i--;
            }
            
        }
        closeImages(imgDOG);
      
        return(pmlPop);
    } 
    
    /** 
     * Find dots with DOG CLIJ Method
     * @param img channel
     * @return dots population
     */
    public Objects3DPopulation findDotsDoGCLIJ(ImagePlus img, Object3D nucObj) {
        ClearCLBuffer imgCL = clij2.push(img);
        double sig1 = radius/(3*img.getCalibration().pixelWidth);
        double sig2 = radius/img.getCalibration().pixelWidth;
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
    public ImagePlus saveImageObjects(ArrayList<Objects3DPopulation> pmlPopList, Objects3DPopulation nucPop, ImagePlus[] imgArray, String pathName, ArrayList<Transformer> trans) {
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
            // labelsObject(pmlObj, imhObjects.getImagePlus(), o, 255);saveI
            }
            imhObjects.getImagePlus().setCalibration(cal);
            imhObjects.getImagePlus().setSlice(imhObjects.getImagePlus().getNSlices()/2);
            IJ.resetMinAndMax(imhObjects.getImagePlus());
            hyperBin[i] = imhObjects.getImagePlus();
             if (i>0) (trans.get(i-1)).doTransformation(hyperBin[i]);
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
    public ImagePlus saveImagePMLs(ArrayList<Objects3DPopulation> pmlPopList, ImagePlus[] imgArray, String pathName) {
       ImagePlus[] hyperPMLOrg = new ImagePlus[imgArray.length];
       ImagePlus[] hyperPML = new ImagePlus[imgArray.length];
        for (int i = 0; i < imgArray.length; i++) {
           ImageHandler imh = ImageHandler.wrap(imgArray[i]).createSameDimensions();
           Objects3DPopulation pmlPop = pmlPopList.get(i);
           pmlPop.draw(imh, 255);
           hyperPML[i] = imh.getImagePlus(); 
           hyperPMLOrg[i] = imgArray[i];
          }
        ImagePlus hyperPMLTime = new Concatenator().concatenate(hyperPML, true);
        ImagePlus hyperPMLOrgTime = new Concatenator().concatenate(hyperPMLOrg, true);
        // Save diffus
        FileSaver imgPMLOrg = new FileSaver(hyperPMLOrgTime);
        imgPMLOrg.saveAsTiff(pathName);
        closeImages(hyperPMLOrgTime);
        return(hyperPMLTime);
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
        ImagePlus hyperDifuseTime = new Concatenator().concatenate(hyperDifuse, false);
        // Save diffus
        FileSaver imgDiffus = new FileSaver(hyperDifuseTime);
        imgDiffus.saveAsTiff(pathName);
        closeImages(hyperDifuseTime);
    }
    
     
     
    public ImagePlus drawNucleus(Objects3DPopulation pop, ImagePlus[] imArray, boolean label, int nnuc) {
        ImagePlus[] hyperBin = new ImagePlus[imArray.length];
        // Draw at each time
        for (int i=0; i<pop.getNbObjects(); i++)
        {
            ImageHandler imhObjects = ImageHandler.wrap(imArray[i]).createSameDimensions();
            pop.getObject(i).draw(imhObjects, 255); 
            imhObjects.getImagePlus().setCalibration(cal);
            hyperBin[i] = imhObjects.getImagePlus();
            if ( label ){
                Object3D ob = pop.getObject(i);
                //double meanrad = Math.pow(ob.getVolumePixels()*3.0/(4.0*Math.PI), 1.0/3.0); // add mean radius to be outside object to write label
                double meanrad = 0; // write in the center
                int z = (int) ob.getCenterZ();
                int x = (int) (ob.getCenterX() + meanrad*1.3);
                int y = (int) (ob.getCenterY() + meanrad*1.3);
                hyperBin[i].setSlice(z);
                ImageProcessor proc = hyperBin[i].getProcessor();
                proc.setColor(255);
                proc.setFontSize(30);
                proc.drawString(""+nnuc, x, y);
                hyperBin[i].updateAndDraw();
            }
        }
        // if no objects on some images at the end: Change by putting pop in a map to associate with the slice ?
        if (pop.getNbObjects()<imArray.length){
            for (int j=pop.getNbObjects(); j<imArray.length; j++){
                ImageHandler imhObjects = ImageHandler.wrap(imArray[j]).createSameDimensions();
                imhObjects.getImagePlus().setCalibration(cal);
                hyperBin[j] = imhObjects.getImagePlus();
            }
        }
        return new Concatenator().concatenate(hyperBin, false);
    }
    
     public ImagePlus drawPMLs(ArrayList<Objects3DPopulation> pmlPopList, ImagePlus[] imgArray) {
       ImagePlus[] hyperPML = new ImagePlus[imgArray.length];
        for (int i = 0; i < pmlPopList.size(); i++) {
           ImageHandler imh = ImageHandler.wrap(imgArray[i]).createSameDimensions();
           Objects3DPopulation pmlPop = pmlPopList.get(i);
           pmlPop.draw(imh, 255);
           hyperPML[i] = imh.getImagePlus(); 
        }
        // if no objects on some images at the end: Change by putting pop in a map to associate with the slice ?
        if (pmlPopList.size()<imgArray.length){
        for (int j=pmlPopList.size(); j<imgArray.length; j++){
            ImageHandler imhObjects = ImageHandler.wrap(imgArray[j]).createSameDimensions();
            imhObjects.getImagePlus().setCalibration(cal);
            hyperPML[j] = imhObjects.getImagePlus();
        }
        }
        return new Concatenator().concatenate(hyperPML, false);
    }
    
     /** align */
     public ImagePlus alignStack(ArrayList<Transformer> transformers, ImagePlus imp) {
         ImagePlus[] hyper = new ImagePlus[transformers.size()+1];
         SubHyperstackMaker sub = new SubHyperstackMaker();
         hyper[0] = sub.makeSubhyperstack(imp, "1-1", "1-"+imp.getNSlices(), "1-1");
         for ( int i = 1; i <= transformers.size(); i++) {
            Transformer trans = transformers.get(i-1);
            ImagePlus cur = sub.makeSubhyperstack(imp, "1-1", "1-"+imp.getNSlices(), (i+1)+"-"+(i+1));
            trans.doTransformation(cur);
            hyper[i] = cur;
        }
        return new Concatenator().concatenate(hyper, false);
     }
     
     public void saveWholeImage(ImagePlus[] imgs, String fileName) {
        
         ImagePlus img = new Concatenator().concatenate(imgs, false);
        double factor = 750.0/img.getWidth();  // diminue la taille pour pas enregistrer une image trop big
        ImagePlus resized = img.resize((int)(img.getWidth()*factor), (int)(img.getHeight()*factor), "bilinear");
        FileSaver fileImg = new FileSaver(img);
        fileImg.saveAsTiff(fileName);
        closeImages(img);
     }
     
     /**
     * Save image objects
     */
    public ImagePlus alignAndSave(ArrayList<Objects3DPopulation> pmlPopList, Objects3DPopulation nucPop, ImagePlus[] imgArray, String root, int index, boolean saveObj) {

        ImagePlus unnucl = drawNucleus(nucPop, imgArray, false, 0);
        // Get transformations to do to align stack
        ArrayList<Transformer> trans = new StackReg_Plus().stackRegister(stackProj(unnucl));
        
        ImagePlus nucleus = null;
        if (saveObj) {
        // align nucleus stack
        nucleus = alignStack(trans, unnucl);
        // rebinarize
         IJ.setAutoThreshold(nucleus, "Default dark stack");
         Prefs.blackBackground = false;
         IJ.run(nucleus, "Convert to Mask", "method=Default background=Dark stack");
         IJ.run(nucleus, "Multiply...", "value=0.5 stack"); 
        }
        closeImages(unnucl);

         // draw and align PMLs
         ImagePlus unpml = drawPMLs(pmlPopList, imgArray);
         // align pml stack
         ImagePlus pml = alignStack(trans, unpml);
         closeImages(unpml);
         
        // rebinarize
         IJ.setAutoThreshold(pml, "Default dark stack");
         Prefs.blackBackground = false;
         IJ.run(pml, "Convert to Mask", "method=Default background=Dark stack");
         
         if (saveObj){
         // Save images objects 
        ImageCalculator ic = new ImageCalculator();
        ic.run("Add stack", nucleus, pml);
        FileSaver fileObjects = new FileSaver(nucleus);
        fileObjects.saveAsTiff(root+"_Objects-"+index+".tif");
        closeImages(nucleus);
        }
          
        FileSaver filePML = new FileSaver(pml);
        filePML.saveAsTiff(root+"_PMLs-"+index+".tif");
        //closeImages(pml);
        return pml;      
    }
    
     /**
     * Save image objects
     */
    public void drawOnWholeImage(ArrayList<Objects3DPopulation> pmlPopList, Objects3DPopulation nucPop, ImagePlus[] imgArray, Roi roi, int index) {
        translateRoiBack(nucPop, roi);
        ImagePlus imgnucl = drawNucleus(nucPop, imgArray, true, index);
        // draw and align PMLs
        translateRoiBack(pmlPopList, roi);
        ImagePlus imgpml = drawPMLs(pmlPopList, imgArray);
        // Save images objects 
        ImageCalculator ic = new ImageCalculator();
        for (int i=0; i<imgArray.length; i++){
            ImagePlus framenuc = new Duplicator().run(imgnucl, 1, 1, 1, imgnucl.getNSlices(), i+1, i+1);  
            ImagePlus framepml = new Duplicator().run(imgpml, 1, 1, 1, imgpml.getNSlices(), i+1, i+1);
            IJ.run(framenuc, "Multiply...", "value=0.5 stack");        
            ic.run("Add stack", framenuc, framepml);
            closeImages(framepml);
        
            ic.run("Add stack", imgArray[i], framenuc);
            closeImages(framenuc);
        }     
        closeImages(imgnucl);
        closeImages(imgpml);
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
        //ImageHandler seedsImg = ImageInt.wrap(binaryMask);
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
    
     /** \brief keep only where enough signal and draw binary mask */
	public void makeNucleiMask(ImagePlus ip, ImagePlus bin)
	{
                RoiManager rm = RoiManager.getInstance();
		Roi[] rois = rm.getRoisAsArray();
		rm.reset();
		IJ.setForegroundColor(255, 255, 255);
                Roi roi = new Roi(0,0, ip.getWidth(), ip.getHeight());
                ip.setRoi(roi);
                IJ.run(ip, "Clear", "stack");
		for ( int i=0; i < rois.length; i++ )
		{
			Roi cur = rois[i];
                        bin.setSlice(cur.getPosition());
       
                        ImageStatistics stat = bin.getStatistics();
                        // contains part of the nucleus
                        if (stat.mean > 10)
                        {
                            ip.setSlice(cur.getPosition());
                            ip.setRoi(cur);
                            IJ.run(ip, "Fill", "slice");
                        }
		}
		rm.runCommand(ip,"Deselect");
	}
        
         /** \brief Change the size of the Roi to scale it to new image size, keep only where enough signal and draw binary mask */
	public void makeNucleiResizeMask(ImagePlus ip, double factxy, ImagePlus bin)
	{
                RoiManager rm = RoiManager.getInstance();
		Roi[] rois = rm.getRoisAsArray();
		rm.reset();
		RoiScaler scaler = new RoiScaler();
                IJ.setForegroundColor(255, 255, 255);
                Roi roi = new Roi(0,0, ip.getWidth(), ip.getHeight());
                ip.setRoi(roi);
                IJ.run(ip, "Clear", "stack");
		for ( int i=0; i < rois.length; i++ )
		{
			Roi cur = rois[i];
                        bin.setSlice(cur.getPosition());
       
                        ImageStatistics stat = bin.getStatistics();
                        // contains part of the nucleus
                        if (stat.mean > 10)
                        {
                            RoiEnlarger re = new RoiEnlarger();
                            Roi shrink = re.enlarge(cur, -2);  // very small shrink or larger and then compensate ??
                            Roi scaled = scaler.scale(shrink, factxy, factxy, false);   
                            ip.setT(cur.getPosition());
                            scaled.setPosition(cur.getPosition());
                            ip.setRoi(scaled);
                            //rm.addRoi(scaled);
                            IJ.run(ip, "Fill", "slice");
                        }
		}
		rm.runCommand(ip,"Deselect");
	}
    
    
    /** Look for all nuclei
     Do z slice by slice stardist 
     */
    public ImagePlus stardistNuclei(ImagePlus imgNuc){
        double factor = 300.0/imgNuc.getWidth();
        ImagePlus resized = imgNuc.resize((int)(imgNuc.getWidth()*factor), (int)(imgNuc.getHeight()*factor), "bilinear");
        StarDist2D star = new StarDist2D();
        star.loadInput(resized);
        star.run();
    
        // try to get rid of extreme slices without nuclei
        ImagePlus globalBin = new Duplicator().run(imgNuc);
        IJ.setAutoThreshold(globalBin, "Otsu dark stack");
        Prefs.blackBackground = false;
        IJ.run(globalBin, "Convert to Mask", "method=Otsu background=Dark stack");
    
        //ImagePlus back = resized.resize(imgNuc.getWidth(), imgNuc.getHeight(), "bilinear");
        //closeImages(resized);
        closeImages(imgNuc);
        
        /**makeNucleiResizeMask(back, 1.0/factor, globalBin);
        back.setDimensions(1, back.getNSlices(), 1);
        back.setCalibration(cal);
        // back.show();
       // new WaitForUserDialog("test").show(); */
         makeNucleiMask(resized, globalBin);
        resized.setDimensions(1, resized.getNSlices(), 1);
        resized.setCalibration(cal);
        
       return resized; 
      
    }
    
     /** Look for all nuclei
     Do z slice by slice stardist 
     * return nuclei population
     */
    public Objects3DPopulation stardistNucleiPop(ImagePlus imgNuc){
        double factor = 300.0/imgNuc.getWidth();
        ImagePlus resized = imgNuc.resize((int)(imgNuc.getWidth()*factor), (int)(imgNuc.getHeight()*factor), "bilinear");
        resized.setCalibration(cal);
        IJ.run(resized, "Remove Outliers", "block_radius_x=2 block_radius_y=2 standard_deviations=1 stack");
        StarDist2D star = new StarDist2D();
        star.loadInput(resized);
        star.run();
    
        // try to get rid of extreme slices without nuclei
        ImagePlus globalBin = new Duplicator().run(resized);
        if ( globalBin.getNFrames()>1) globalBin.setDimensions(1, globalBin.getNFrames(), 1);
        //globalBin.setT(globalBin.getNFrames()/2);
        IJ.setAutoThreshold(globalBin, "Otsu dark stack");
        Prefs.blackBackground = false;
        IJ.run(globalBin, "Convert to Mask", "method=Otsu background=Dark stack");
        
        globalBin.setSlice(globalBin.getNSlices()/2);
        ImageStatistics stat = globalBin.getStatistics();
        // threshold didn't work well, very few pixels found
        if ( stat.mean < 10 ) {
            closeImages(globalBin);
            globalBin = new Duplicator().run(resized);
            if ( globalBin.getNFrames()>1) globalBin.setDimensions(1, globalBin.getNFrames(), 1);
            IJ.setAutoThreshold(globalBin, "Default dark stack");
            Prefs.blackBackground = false;
            IJ.run(globalBin, "Convert to Mask", "method=Default background=Dark stack");
        }
        
        //globalBin.show();
        //new WaitForUserDialog("t").show();
    
        ImagePlus back = resized.resize(imgNuc.getWidth(), imgNuc.getHeight(), "bilinear");
        closeImages(resized);
        closeImages(imgNuc);
        
        makeNucleiResizeMask(back, 1.0/factor, globalBin);
        if (back.getNFrames()>1) back.setDimensions(1, back.getNFrames(), 1);
        back.setCalibration(cal);
        // back.show();
       // new WaitForUserDialog("test").show(); 
        
       Objects3DPopulation nucPop = new Objects3DPopulation(getPopFromImage(back).getObjectsWithinVolume(minNuc, maxNuc, true));
       closeImages(back);
       return(nucPop);
    }
      
    /** Look for one nuclei in a Roi
     Do z slice by slice stardist 
     */
    public Object3D stardistRoi(ImagePlus imgNuc){
        double factor = 40.0/imgNuc.getWidth();
        ImagePlus resized = imgNuc.resize((int)(imgNuc.getWidth()*factor), (int)(imgNuc.getHeight()*factor), "bilinear");
        StarDist2D star = new StarDist2D();
        star.loadInput(resized);
        star.run();
        
        // try to get rid of extreme slices without nuclei
        ImagePlus globalBin = new Duplicator().run(imgNuc);
        IJ.setAutoThreshold(globalBin, "Otsu dark stack");
        Prefs.blackBackground = false;
        IJ.run(globalBin, "Convert to Mask", "method=Otsu background=Dark stack");
    
        ImagePlus back = resized.resize(imgNuc.getWidth(), imgNuc.getHeight(), "bilinear");
        closeImages(resized);
        closeImages(imgNuc);
        
        makeNucleiResizeMask(back, 1.0/factor, globalBin);
        back.setDimensions(1, back.getNSlices(), 1);
        back.setCalibration(cal);
       
        Objects3DPopulation nucPop = new Objects3DPopulation(getPopFromImage(back).getObjectsWithinVolume(minNuc, maxNuc, true));
        // Find bigger object
        int ind = 0;    
        if (nucPop.getNbObjects()>1){
            double maxv = 0;
            for (int k=0; k<nucPop.getNbObjects(); k++) {
                Object3D obj = nucPop.getObject(k);
                double vol = obj.getVolumePixels();
                if (vol>maxv){
                    ind = k;
                    maxv = vol;
                }
            }
            
        }
        if ( nucPop.getNbObjects()==0 ) {
            back.show();
            ind = 0;
            return null;
        }
        Object3D nucObj = nucPop.getObject(ind);
        closeImages(back);
        return(nucObj);
      
    }
    
    
   /** Looks for the closest nuclei at each time (track it). Loose it if too far from given distance */ 
    public Objects3DPopulation trackNucleus( Object3D obj, ArrayList<Objects3DPopulation> pop) {
        Objects3DPopulation nucl = new Objects3DPopulation();    
        for (int i=0; i<pop.size(); i++) {
            Object3D closest = (pop.get(i)).closestCenter(obj.getCenterAsPoint());
            // threshold distance to loose the nuclei (not aligned image so can move)
            if (obj.distCenterUnit(closest) > 6) {
                return nucl;
            }
            // within distance, continue
            nucl.addObject(closest);
            obj = closest;
        }
        return nucl;
    }
    
    public int[] getBoundingBoxXY(Objects3DPopulation pop) {
        int[] bb = new int[]{10000,0,100000,0};
        for (int k=0; k<pop.getNbObjects(); k++) {
                Object3D obj = pop.getObject(k);
                int[] bbobj = obj.getBoundingBox();
                if ( bbobj[0] < bb[0] ) bb[0] = bbobj[0];
                if ( bbobj[2] < bb[2] ) bb[2] = bbobj[2];
                if ( bbobj[1] > bb[1] ) bb[1] = bbobj[1];
                if ( bbobj[3] > bb[3] ) bb[3] = bbobj[3];
        }
        return bb;
    }

    public void translateToRoi( Objects3DPopulation pop, int[] roi) {
        //Objects3DPopulation res = new Objects3DPopulation();
        for (int k=0; k<pop.getNbObjects(); k++) {
                Object3D obj = pop.getObject(k);
                obj.translate(-roi[0], -roi[2], 0);
               // res.addObject(obj);
        }
        //return res;
    }
    
     public void translateRoiBack( Objects3DPopulation pop, Roi roi) {
        for (int k=0; k<pop.getNbObjects(); k++) {
                Object3D obj = pop.getObject(k);
                obj.translate(roi.getBounds().getMinX(), roi.getBounds().getMinY(), 0);
        }
      }
     
     public void translateRoiBack( ArrayList<Objects3DPopulation> array, Roi roi) {
        for (int i=0; i<array.size(); i++){
            Objects3DPopulation pop = array.get(i);
            for (int k=0; k<pop.getNbObjects(); k++) {
                Object3D obj = pop.getObject(k);        
                obj.translate(roi.getBounds().getMinX(), roi.getBounds().getMinY(), 0);
            }
        }
    }
}