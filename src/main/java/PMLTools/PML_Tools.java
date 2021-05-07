package PMLTools;


import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.GaussianBlur3D;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import static java.lang.Double.parseDouble;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.ThreadLocalRandom;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Object3D_IJUtils;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.FastFilters3D;
import mcib3d.image3d.regionGrowing.Watershed3D;
import mcib3d.spatial.analysis.SpatialStatistics;
import mcib3d.spatial.descriptors.F_Function;
import mcib3d.spatial.descriptors.SpatialDescriptor;
import mcib3d.spatial.sampler.SpatialModel;
import mcib3d.spatial.sampler.SpatialRandomHardCore;
import mcib3d.utils.ArrayUtil;
import mcib3d.utils.CDFTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
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
   
    // min nucleus volume in pixels^3
    public static double minNuc = 10;
    // max nucleus volume in pixels^3
    public static double maxNuc = 500;
    
    // min volume in pixels^3 for dots
    public static double minPML = 0.1;
    // max volume in picxels^3 for dots
    public static double maxPML = 5;
    private static Calibration cal = new Calibration(); 
    
    // dots dog paramaters
    private static int sig1 = 2;
    private static int sig2 = 4;
    // dots threshold method
    private static String thMet = "Moments";
    
    public static CLIJ2 clij2 = CLIJ2.getInstance();
    
     /**
     * check  installed modules
     * @return 
     */
    public static boolean checkInstalledModules() {
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
    public static ArrayList findImages(String imagesFolder, String imageExt) {
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
    public static List<String> findChannels (String imageName) throws DependencyException, ServiceException, FormatException, IOException {
        List<String> channels = new ArrayList<>();
        // create OME-XML metadata store of the latest schema version
        ServiceFactory factory;
        factory = new ServiceFactory();
        OMEXMLService service = factory.getInstance(OMEXMLService.class);
        IMetadata meta = service.createOMEXMLMetadata();
        ImageProcessorReader reader = new ImageProcessorReader();
        reader.setMetadataStore(meta);
        reader.setId(imageName);
        int chs = reader.getSizeC();
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                String channelsID = meta.getImageName(0);
                channels = Arrays.asList(channelsID.replace("_", "-").split("/"));
                break;
            case "lif" :
                String[] ch = new String[chs];
                if (chs > 1) {
                    for (int n = 0; n < chs; n++) 
                        if (meta.getChannelExcitationWavelength(0, n) == null)
                            channels.add(Integer.toString(n));
                        else 
                            channels.add(meta.getChannelExcitationWavelength(0, n).value().toString());
                }
                break;
            default :
                chs = reader.getSizeC();
                for (int n = 0; n < chs; n++)
                    channels.add(Integer.toString(n));
        }
        return(channels);         
    }
    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public static Calibration findImageCalib(IMetadata meta) {
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
    public static void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    public static Objects3DPopulation getPopFromImage(ImagePlus img) {
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
    public static ClearCLBuffer medianFilter(ClearCLBuffer imgCL, double sizeX, double sizeY, double sizeZ) {
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
    public static ClearCLBuffer DOG(ClearCLBuffer imgCL, double size1, double size2) {
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, size1, size1, size1, size2, size2, size2);
        clij2.release(imgCL);
        return(imgCLDOG);
    }
    
    /**
     * Fill hole
     * USING CLIJ2
     */
    private static void fillHole(ClearCLBuffer imgCL) {
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
    public static ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed, boolean fill) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        if (fill)
            fillHole(imgCLBin);
        return(imgCLBin);
    }
    
    /**
     * return objects population in an binary image
     * Using CLIJ2
     * @param imgCL
     * @return pop
     */

    private static Objects3DPopulation getPopFromClearBuffer(ClearCLBuffer imgCL, Calibration cal) {
        ClearCLBuffer output = clij2.create(imgCL);
        clij2.connectedComponentsLabelingBox(imgCL, output);
        ImagePlus imgLab  = clij2.pull(output);
        imgLab.setCalibration(cal); 
        ImageInt labels = new ImageLabeller().getLabels(ImageHandler.wrap(imgLab));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        clij2.release(output);
        imgLab.close();
        return pop;
    }  
    
    
    /**
     * Dialog 
     * A finir paramètres de taille etc ....
     */
    public static String dialog() {
        String dir = "";
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addDirectoryField("Choose Directory Containing Image Files : ", "");
        gd.addNumericField("Min PML size (µm3): ", minPML, 3);
        gd.addNumericField("Max PML size (µm3): ", maxPML, 3);
        gd.showDialog();
        dir = gd.getNextString();
        minPML = gd.getNextNumber();
        maxPML = gd.getNextNumber();
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
    public static ImagePlus removeOutliers(ImagePlus img, int radX, int radY, float factor) {
        
        for (int i = 0; i < img.getNSlices(); i++) {
            img.setSlice(i);
            ImageProcessor ip = img.getProcessor();
            RemoveOutliers removeOut = new RemoveOutliers(ip.convertToFloatProcessor());
            removeOut.removeOutliers(radX, radY, factor);
        }
        return(img);
    }
    
    
/**
     * Nucleus segmentation 2
     * @param imgNuc
     * @param time
     * @return 
     */
    public static Object3D findnucleus(ImagePlus imgNuc, int time) {
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
        Objects3DPopulation nucPop = new Objects3DPopulation(getPopFromImage(imgStack).getObjectsWithinVolume(minNuc, maxNuc, true));
        nucPop.removeObjectsTouchingBorders(imgStack, false);
        Object3D nucObj = nucPop.getObject(0);
        nucObj.setValue(time);
        closeImages(imgStack);
        return(nucObj);
    }
    
    
    // Threshold images and fill holes
    public static void threshold(ImagePlus img, AutoThresholder.Method thMed, boolean fill) {
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
     * Fin image drift correction
     * @param img
     * @return 
     */
    public static double[] correctDrift(ImagePlus img) {
       double[] mat = new double[3];
        
       
       return(mat); 
    }
    
    /** 
     * Find dots
     * @param img channel
     * @return dots population
     */
    public static Objects3DPopulation findDots(ImagePlus img, String name) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLMed = clij2.create(imgCL);
        clij2.mean3DBox(imgCL, imgCLMed, 1, 1, 1);
        clij2.copy(imgCLMed, imgCL);
        clij2.release(imgCLMed);
        ClearCLBuffer imgCLDOG = DOG(imgCLMed, sig1, sig2);
        clij2.release(imgCL);
        ImagePlus imgBin = clij2.pull(threshold(imgCLDOG, thMet, false));
        clij2.release(imgCLDOG);
        imgBin.setCalibration(img.getCalibration());
        Objects3DPopulation pmlPop = new Objects3DPopulation(getPopFromImage(imgBin).getObjectsWithinVolume(minPML, maxPML, true));
//        FileSaver fociMaskFile = new FileSaver(imgBin);
//        fociMaskFile.saveAsTiff(name);
        closeImages(imgBin);
        return(pmlPop);
    } 
    
    
    public static Objects3DPopulation getPopFromImage(ImagePlus img, Calibration cal) {
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
     * @param integrated
     */
    public static void pmlDiffus(Objects3DPopulation pmlPop, Object3D nucObj, ImagePlus imgPML, boolean integrated) {
        ImageHandler imhDotsDiffuse = ImageHandler.wrap(imgPML.duplicate());
        double pmlIntDiffuse ;
        int pmlDots = pmlPop.getNbObjects();
        double volPMLDilated = 0;
        float dilate = 1.5f;
        for (int p = 0; p < pmlDots; p++) {
            Object3D pmlObj = pmlPop.getObject(p);
            // dilate 
            Object3DVoxels pmlDilatedObj = pmlObj.getDilatedObject(dilate, dilate, dilate);
            pmlDilatedObj.draw(imhDotsDiffuse, 0);
            volPMLDilated += pmlDilatedObj.getVolumeUnit();
        }
        double nucVolume = nucObj.getVolumeUnit();
        if (integrated)
            pmlIntDiffuse = (nucObj.getIntegratedDensity(imhDotsDiffuse)/nucVolume) * (nucVolume - volPMLDilated); 
        else {
            ArrayUtil pixelInt = nucObj.listValues​(imhDotsDiffuse);
            pmlIntDiffuse = pixelInt.getMean();
        }
        // put in comment nucleus object diffus intensity 
        nucObj.setComment(Double.toString(pmlIntDiffuse));
        imhDotsDiffuse.closeImagePlus();

    }
    
    /**
     * Create diffuse PML image
     * Fill zero in pml dots
     * @param pmlPop
     * @param nucObj
     * @param imgDotsOrg
     * @param outDirResults
     * @param rootName
     * @param index
     * @param seriesName
     */
    public static void saveDiffusImage(Objects3DPopulation pmlPop, Object3D nucObj, ImagePlus imgDotsOrg, String outDirResults,
            String rootName, String seriesName, int index) {
        ImageHandler imhDotsDiffuse = ImageHandler.wrap(imgDotsOrg.duplicate());
       
        float dilate = 1.5f;
        for (int p = 0; p < pmlPop.getNbObjects(); p++) {
            Object3D pmlObj = pmlPop.getObject(p);
            // dilate 
            Object3DVoxels pmlDilatedObj = pmlObj.getDilatedObject(dilate, dilate, dilate);
            pmlDilatedObj.draw(imhDotsDiffuse, 0);
        }
        ImagePlus imgColor = imhDotsDiffuse.getImagePlus();
        IJ.run(imgColor, "RGB Color", "");
        drawCountours(nucObj, imgColor, Color.white);
        // Save diffus
        FileSaver imgDiffus = new FileSaver(imgColor);
        imgDiffus.saveAsTiff(outDirResults + rootName + "_" + seriesName + "-Nuc"+index+"_PMLDiffus.tif");
        closeImages(imgColor);
        imhDotsDiffuse.closeImagePlus();
    }
    
    
    
    /*
    * Find PML inside nucleus
    * tag nucleus name as nucleus index
    * tag nucleus value as number of pml inside
    */
    public static Objects3DPopulation coloc(Object3D nucObj, Objects3DPopulation pmlPop) {
        Objects3DPopulation pmlInNuc = new Objects3DPopulation();
        int pml = 0;
        for (int p = 0; p < pmlPop.getNbObjects(); p++) {
            Object3D pmlObj = pmlPop.getObject(p);
            if (pmlObj.hasOneVoxelColoc(nucObj)) {
                pmlInNuc.addObject(pmlObj);
                pml++;
            }
        }
        nucObj.setValue(pml);
        return(pmlInNuc);
    }
    
    
    // write object labels
    public static void labelsObject (Object3D obj, ImagePlus img, int number, int color) {
        Font tagFont = new Font("SansSerif", Font.PLAIN, 32);
        int[] box = obj.getBoundingBox();
        int z = (int)obj.getCenterZ();
        int x = box[0] - 2;
        int y = box[2] - 2;
        img.setSlice(z+1);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(tagFont);
        ip.setColor(color);
        ip.drawString(Integer.toString(number), x, y);
        img.updateAndDraw();    
    }
    
    // tag object number with random color
    public static void tagsObject(ImageHandler imh, Object3D obj) {        
        int col = ThreadLocalRandom.current().nextInt(2, 255 + 1);
        obj.draw(imh, col);  
    }
    
   /**
     * ramdom color nucleus population
     */
    public static ImagePlus randomColorPop (Objects3DPopulation cellsPop,  ImageHandler img, boolean label) {
        //create image objects population
        img.set332RGBLut();
        img.setCalibration(img.getCalibration());
        for (int i = 0; i < cellsPop.getNbObjects(); i++) {
            Object3D obj = cellsPop.getObject(i);
            int col = ThreadLocalRandom.current().nextInt(2, 255 + 1);
            obj.draw(img, col);
            if (label)
               labelsObject(obj, img.getImagePlus(), (i+1), col); 
        } 
        return(img.getImagePlus());
    } 

    
    public static Point3D[] createEvaluationPoints(int numPoints, Objects3DPopulation population, Object3D mask) {
        Point3D[] evaluationPoints = new Point3D[numPoints];
        population.setMask(mask);
        for (int i = 0; i < numPoints; ++i) {
            evaluationPoints[i] = population.getRandomPointInMask();
        }
        return evaluationPoints;
    }
   

    /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public static void clearOutSide(ImagePlus img, Roi roi) {
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(roi);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(roi);
        }
        img.updateAndDraw();
    }
    
    public static Plot createPlot(ArrayUtil xEvals, ArrayUtil[] sampleDistances, ArrayUtil observedDistances, ArrayUtil observedCD, ArrayUtil averageCD, String function) {
     
        Color ColorAVG = Color.blue;
        Color ColorENV = Color.gray;
        Color ColorOBS = Color.red;
        double plotMaxX = observedDistances.getMaximum();
        double plotMaxY = observedCD.getMaximum();
        int nbBins = 100;
        double env = 0.5;
        // low env
        double max = xEvals.getMaximum();
        ArrayUtil xEval0 = new ArrayUtil(nbBins);
        for (int i = 0; i < nbBins; i++) {
            xEval0.addValue(i, ((double) i) * max / ((double) nbBins));
        }
        // get the values
        ArrayUtil samplesPc5 = CDFTools.cdfPercentage(sampleDistances, xEval0, env / 2.0);
        ArrayUtil samplesPc95 = CDFTools.cdfPercentage(sampleDistances, xEval0, 1.0 - env / 2.0);
        // get the limits
        if (xEval0.getMaximum() > plotMaxX) {
            plotMaxX = xEval0.getMaximum();
        }
        if (samplesPc5.getMaximum() > plotMaxY) {
            plotMaxY = samplesPc5.getMaximum();
        }
        if (samplesPc95.getMaximum() > plotMaxY) {
            plotMaxY = samplesPc95.getMaximum();
        }
        if (xEvals.getMaximum() > plotMaxX) {
            plotMaxX = xEvals.getMaximum();
        }
        if (averageCD.getMaximum() > plotMaxY) {
            plotMaxY = averageCD.getMaximum();
        }
        if (observedCD.getMaximum() > plotMaxY) {
            plotMaxY = observedCD.getMaximum();
        }
        if (observedDistances.getMaximum() > plotMaxX) {
            plotMaxX = observedDistances.getMaximum();
        }
        // create the plot
        Plot plot = new Plot(function + "-function", "distance", "cumulated frequency");
        plot.setLimits(0, plotMaxX, 0, plotMaxY);

        // envelope  for e.g 10 % at 5 and 95 %
        plot.setColor(ColorENV);
        plot.addPoints(xEval0.getArray(), samplesPc5.getArray(), Plot.LINE);

        // envelope  for e.g 10 % at 5 and 95 %
        plot.setColor(ColorENV);
        plot.addPoints(xEval0.getArray(), samplesPc95.getArray(), Plot.LINE);

        // average
        plot.setColor(ColorAVG);
        plot.addPoints(xEvals.getArray(), averageCD.getArray(), Plot.LINE);

        // observed
        plot.setColor(ColorOBS);
        plot.addPoints(observedDistances.getArray(), observedCD.getArray(), Plot.LINE);

        return plot;
    }
    
    
     
    /**
    * For each nucleus compute F function
     * @param pop
     * @param mask
     * @param nuc
     * @param imgName
     * @param outDirResults
     * @return F SDI
    **/ 
    public static double processF (Objects3DPopulation pop, Object3D mask, String imgName, int nuc, String outDirResults) {

        // define spatial descriptor, model
        SpatialDescriptor spatialDesc = new F_Function(pop.getNbObjects(), mask);
        pop.setMask(mask);
        SpatialModel spatialModel = new SpatialRandomHardCore(pop.getNbObjects(), 0.8, mask);
        SpatialStatistics spatialStatistics = new SpatialStatistics(spatialDesc, spatialModel, 50, pop);
        spatialStatistics.setEnvelope(0.25);
        spatialStatistics.setVerbose(false);
//        Plot fPlot = spatialStatistics.getPlot();
//        fPlot.draw();
//        fPlot.addLabel(0.1, 0.1, "p = " + String.valueOf(spatialStatistics.getSdi()));
//        ImagePlus imgPlot = fPlot.getImagePlus();
//        FileSaver plotSave = new FileSaver(imgPlot);
//        plotSave.saveAsTiff(outDirResults + imgName + "_Fplot_" + nuc + ".tif");
//        flush_close(imgPlot);
        System.out.println("Nucleus" + nuc + " Sdi = " + spatialStatistics.getSdi());
        return(spatialStatistics.getSdi());
    }
    
    
    public static ImagePlus WatershedSplit(ImagePlus binaryMask, float rad) {
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
    * Compute nucleus and pml results
    * @param nucObj nucleus
    * @param pmlPop pml population
    * @param imgPML read pml intensity
    * @param imgName image file
     * @param results buffer
    **/
    public static void computeNucParameters(Object3D nucObj, Objects3DPopulation pmlPop, ImagePlus imgPML, String imgName, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing nucleus parameters ....");
        ImageHandler imhPML = ImageHandler.wrap(imgPML);
        // measure nucleus volume
        // measure pml integrated intensity and volume
        double  nucVolume = nucObj.getVolumeUnit();
        for (int p = 0; p < pmlPop.getNbObjects(); p++) {
            Object3D pmlObj = pmlPop.getObject(p);
            double pmlIntensity = pmlObj.getIntegratedDensity(imhPML);
            double pmlVolume = pmlObj.getVolumeUnit();
            double minDistBorder = pmlObj.distCenterBorderUnit(nucObj);
            results.write(imgName+"\t"+nucObj.getName()+"\t"+nucVolume+"\t"+(p+1)+"\t"+pmlIntensity+"\t"+pmlVolume+"\t"+minDistBorder+"\n");
            results.flush();
        }
    }
    
    
    /**
    * Compute global nucleus and pml parameters for fixed cells
    * @param nucObj nucleus
     * @param nucIndex
    * @param pmlPop pml population
    * @param imgPML read pml intensity
    * @param imgName image file
    * @param outDirResults results file
     * @param results buffer
     * @throws java.io.IOException
    **/
    public static void computeNucParameters2(Object3D nucObj, int nucIndex, Objects3DPopulation pmlPop, ImagePlus imgPML, 
            String imgName, String outDirResults, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing nucleus parameters ....");
        ImageHandler imhPML = ImageHandler.wrap(imgPML);
        // measure nucleus volume and integrated intensity in PML diffuse
        // measure pml integrated intensity and volume
        DescriptiveStatistics pmlIntensity = new DescriptiveStatistics();
        DescriptiveStatistics pmlVolume = new DescriptiveStatistics();
        DescriptiveStatistics minDistBorder = new DescriptiveStatistics();
        double minDistCenterMean = Double.NaN;
        double minDistCenterSD = Double.NaN;
        double sdiF = Double.NaN;
        int pmlNuc = pmlPop.getNbObjects();
        double nucVolume = nucObj.getVolumeUnit();
        String nucIntDiffuse = nucObj.getComment();
        double nucShericity = nucObj.getSphericity(true);
        for (int p = 0; p < pmlNuc; p++) {
            Object3D pmlObj = pmlPop.getObject(p);
            pmlIntensity.addValue(pmlObj.getIntegratedDensity(imhPML));
            pmlVolume.addValue(pmlObj.getVolumeUnit());
            minDistBorder.addValue(pmlObj.distCenterBorderUnit(nucObj));
        }
        if (pmlPop.getNbObjects() > 4) {
            sdiF = processF(pmlPop, nucObj, imgName, nucIndex, outDirResults);
        }
        if (pmlPop.getNbObjects() > 2) {
            minDistCenterMean = pmlPop.distancesAllClosestCenter().getMean(); 
            minDistCenterSD = pmlPop.distancesAllClosestCenter().getStdDev();
        }
        // compute statistics
        double pmlIntMean = pmlIntensity.getMean();
        double pmlIntSD = pmlIntensity.getStandardDeviation();
        double pmlIntMin = pmlIntensity.getMin();
        double pmlIntMax = pmlIntensity.getMax();
        double pmlVolumeMean = pmlVolume.getMean();
        double pmlVolumeSD = pmlVolume.getStandardDeviation();
        double pmlVolumeMin = pmlVolume.getMin();
        double pmlVolumeMax = pmlVolume.getMax();
        double pmlVolumeSum = pmlVolume.getSum();
        double minDistBorderMean = minDistBorder.getMean();
        double minDistBorderSD = minDistBorder.getStandardDeviation();

        results.write(imgName+"\t"+nucIndex+"\t"+nucVolume+"\t"+nucShericity+"\t"+pmlNuc+"\t"+nucIntDiffuse+"\t"+pmlIntMean+"\t"+
                pmlIntSD+"\t"+pmlIntMin+"\t"+pmlIntMax+"\t"+pmlVolumeMean+"\t"+pmlVolumeSD+"\t"+pmlVolumeMin+"\t"+pmlVolumeMax+"\t"+pmlVolumeSum+"\t"+
                minDistCenterMean+"\t"+minDistCenterSD+"\t"+minDistBorderMean+"\t"+minDistBorderSD+"\t"+sdiF+"\n");
        results.flush();
    }
    
    
    
     /**
    * Compute individual nucleus/pml parameters for live cells at time t
    * @param nucObj nucleus object
    * @param pmlPop pml population
    * @param imgPML read pml intensity
    * @param imgPMLDif read PML diffuse intensity
     * @param time
    * @param results buffer
    **/
    public static void computeNucParameters(Object3D nucObj, Objects3DPopulation pmlPop, ImagePlus imgPML, ImageHandler imgPMLDif, int time, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing nucleus parameters ....");

        // find dots inside nucleus
        // measure nucleus volume and integrated intensity in PML diffuse
        // measure pml integrated intensity and volume
            DescriptiveStatistics pmlIntensity = new DescriptiveStatistics();
            DescriptiveStatistics pmlVolume = new DescriptiveStatistics();
            int pmlinNuc = pmlPop.getNbObjects();
            ImageHandler imh = ImageHandler.wrap(imgPML);
            double nucVolume = nucObj.getVolumeUnit();
            double nucIntDiffuse = nucObj.getIntegratedDensity(imgPMLDif);
            double minDistCenterMean = pmlPop.distancesAllClosestCenter().getMean(); 
            double minDistCenterSD = pmlPop.distancesAllClosestCenter().getStdDev();
            int pmlIndex= 0;
            for (int p = 0; p < pmlPop.getNbObjects(); p++) {
                pmlIndex++;
                Object3D pml = pmlPop.getObject(p);
                pml.setName(String.valueOf(pmlIndex));
                pmlIntensity.addValue(pml.getIntegratedDensity(imh));
                pmlVolume.addValue(pml.getVolumeUnit());
            }
            // compute statistics
            double pmlIntMean = pmlIntensity.getMean();
            double pmlIntSD = pmlIntensity.getStandardDeviation();
            double pmlIntMin = pmlIntensity.getMin();
            double pmlIntMax= pmlIntensity.getMax();
            double pmlVolumeMean = pmlVolume.getMean();
            double pmlVolumeSD = pmlVolume.getStandardDeviation();
            double pmlVolumeMin = pmlVolume.getMin();
            double pmlVolumeMax = pmlVolume.getMax();
            double pmlVolumeSum = pmlVolume.getSum();
            
            results.write(time+"\t"+nucVolume+"\t"+pmlinNuc+"\t"+nucIntDiffuse+"\t"+pmlIntMean+"\t"+pmlIntSD+"\t"+pmlIntMin+"\t"+pmlIntMax+"\t"+
                    pmlVolumeMean+"\t"+pmlVolumeSD+"\t"+pmlVolumeMin+"\t"+pmlVolumeMax+"\t"+pmlVolumeSum+ "\t"+minDistCenterMean+"\t"+minDistCenterSD+"\n");
            results.flush();
    }
    
    /* create xml file for fixed cells
    * Write to xml file nucleus and pml parameters
    *   @param xmlFile file to write
    *   @param nucPop nucleus Population
    *   @param pmlPop pml population
    *   @param ims plm Diffuse image
    */

    public static void writeXml(String xmlFile, Objects3DPopulation nucPop, Objects3DPopulation pmlPop, ImageHandler imgPml, ImageHandler imgDiffus) throws ParserConfigurationException, TransformerConfigurationException, TransformerException {
        // create new XML doc
        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
        Document doc = docBuilder.newDocument();
        // XML root = nucleus
        Element rootElement = doc.createElement("Nucleus");
        doc.appendChild(rootElement);
        for (int n = 0; n < nucPop.getNbObjects(); n++) {
            Object3D nucObj = nucPop.getObject(n);
            // Add nucleus parameters
            Element nuc_param = doc.createElement("nucleus");
            rootElement.appendChild(nuc_param);
            // id
            Attr nucleus_id = doc.createAttribute("id");
            nucleus_id.setValue(nucObj.getName());
            nuc_param.setAttributeNode(nucleus_id);
            // volume
            Attr nucleus_vol = doc.createAttribute("volume");
            nucleus_vol.setValue(Double.toString(nucObj.getVolumeUnit()));
            nuc_param.setAttributeNode(nucleus_vol);
            // pml number
            Attr pml_nb = doc.createAttribute("pml_number");
            pml_nb.setValue(Integer.toString(nucObj.getValue()));
            nuc_param.setAttributeNode(pml_nb);
            // nucleus integrated intensity in diffuse pml image
            Attr diff_Int = doc.createAttribute("Diffuse_intensity");
            diff_Int.setValue(Double.toString(nucObj.getIntegratedDensity(imgDiffus)));
            nuc_param.setAttributeNode(diff_Int);
            for (int p = 0; p < pmlPop.getNbObjects(); p++) {
                Object3D pmlObj =  pmlPop.getObject(p);
                // write pml parameters for nucleus index
                if(pmlObj.getValue() == Integer.parseInt(nucObj.getName())) {
                    // Add pml parameters
                    Element pml_param = doc.createElement("pml");
                    nuc_param.appendChild(pml_param);
                    // id
                    Attr pml_id = doc.createAttribute("id");
                    pml_id.setValue(pmlObj.getComment());
                    pml_param.setAttributeNode(pml_id);
                    // Volume
                    Attr pml_vol = doc.createAttribute("volume");
                    pml_vol.setValue(Double.toString(pmlObj.getVolumeUnit()));
                    pml_param.setAttributeNode(pml_vol);
                    // Intensity
                    Attr pml_intensity = doc.createAttribute("intensity");
                    pml_intensity.setValue(Double.toString(pmlObj.getIntegratedDensity(imgPml)));
                    pml_param.setAttributeNode(pml_intensity);
                    // pml centroids           
                    Attr pml_XCentroid = doc.createAttribute("x_centroid");
                    pml_XCentroid.setValue(Double.toString(pmlObj.getCenterX()));
                    pml_param.setAttributeNode(pml_XCentroid);
                    Attr pml_YCentroid = doc.createAttribute("y_centroid");
                    pml_YCentroid.setValue(Double.toString(pmlObj.getCenterY()));
                    pml_param.setAttributeNode(pml_YCentroid);
                    Attr pml_ZCentroid = doc.createAttribute("z_centroid");
                    pml_ZCentroid.setValue(Double.toString(pmlObj.getCenterZ()));
                    pml_param.setAttributeNode(pml_ZCentroid);
                }
            }
        }
        // write the content into xml file
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        Transformer transformer = transformerFactory.newTransformer();
        DOMSource source = new DOMSource(doc);
        StreamResult result = new StreamResult(new File(xmlFile));
        transformer.transform(source, result);
    }
    
    /* create xml file for live cells
    * Write to xml file nucleus and pml parameters
    *   @param xmlFile file to write
    *   @param nucPop nucleus Population
    *   @param pmlPop pml population
    *   @param ims plm Diffuse image
    */

    public static void writeXml(String xmlFile, Object3D nucObj, Objects3DPopulation pmlPop, ImageHandler imgPml, ImageHandler imgDiffus, int t) throws ParserConfigurationException, TransformerConfigurationException,
            TransformerException, SAXException, IOException {
        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
        Document doc;
        Element rootElement = null;
        // if t = 1 create new XML doc
        // and root as nucleus
        if (t == 1) {
            doc = docBuilder.newDocument();
            // XML root = nucleus
            rootElement = doc.createElement("Nucleus");
            doc.appendChild(rootElement);
            
        }
        // else parse xml
        else {
            doc = docBuilder.parse(xmlFile);
            rootElement = doc.getDocumentElement();
        }
        // Add time
        Element time = doc.createElement("time");
        rootElement.appendChild(time);
        // Add time index
        Attr t_index = doc.createAttribute("index");
        t_index.setValue(Integer.toString(t));
        time.setAttributeNode(t_index);
        // Add nucleus
        Element nucleus = doc.createElement("nucleus");
        time.appendChild(nucleus);
        // nucleus id
        Attr nucleus_id = doc.createAttribute("id");
        nucleus_id.setValue(nucObj.getName());
        nucleus.setAttributeNode(nucleus_id);
        // volume
        Attr nucleus_vol = doc.createAttribute("volume");
        nucleus_vol.setValue(Double.toString(nucObj.getVolumeUnit()));
        nucleus.setAttributeNode(nucleus_vol);
        // pml number
        Attr pml_nb = doc.createAttribute("pml_number");
        pml_nb.setValue(Integer.toString(nucObj.getValue()));
        nucleus.setAttributeNode(pml_nb);
        // nucleus integrated intensity in diffuse pml image
        Attr intensity = doc.createAttribute("diffuse_intensity");
        intensity.setValue(Double.toString(nucObj.getIntegratedDensity(imgDiffus)));
        nucleus.setAttributeNode(intensity);
        // add plm
        int pmlIndex = 0;
        for (int p = 0; p < pmlPop.getNbObjects(); p++) {
            pmlIndex++;
            Object3D pmlObj =  pmlPop.getObject(p);
            // Add pml parameters
            Element pml_param = doc.createElement("pml");
            nucleus.appendChild(pml_param);
            // id
            Attr pml_id = doc.createAttribute("id");
            pml_id.setValue(Integer.toString(pmlIndex));
            pml_param.setAttributeNode(pml_id);
            // Volume
            Attr pml_vol = doc.createAttribute("volume");
            pml_vol.setValue(Double.toString(pmlObj.getVolumeUnit()));
            pml_param.setAttributeNode(pml_vol);
            // Intensity
            Attr pml_intensity = doc.createAttribute("intensity");
            pml_intensity.setValue(Double.toString(pmlObj.getIntegratedDensity(imgPml)));
            pml_param.setAttributeNode(pml_intensity);
            // pml centroids           
            Attr pml_XCentroid = doc.createAttribute("x_centroid");
            pml_XCentroid.setValue(Double.toString(pmlObj.getCenterX()));
            pml_param.setAttributeNode(pml_XCentroid);
            Attr pml_YCentroid = doc.createAttribute("y_centroid");
            pml_YCentroid.setValue(Double.toString(pmlObj.getCenterY()));
            pml_param.setAttributeNode(pml_YCentroid);
            Attr pml_ZCentroid = doc.createAttribute("z_centroid");
            pml_ZCentroid.setValue(Double.toString(pmlObj.getCenterZ()));
            pml_param.setAttributeNode(pml_ZCentroid);
        }

        // write the content into xml file
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        Transformer transformer = transformerFactory.newTransformer();
        DOMSource source = new DOMSource(doc);
        StreamResult result = new StreamResult(new File(xmlFile));
        transformer.transform(source, result);
    }

      /*
    Draw countours of objects
    */
    public static void drawCountours(Object3D obj, ImagePlus img, Color col) {
        ImagePlus imgMask = IJ.createImage("mask", img.getWidth(), img.getHeight(), img.getNSlices(), 8);
        for (int z = obj.getZmin(); z < obj.getZmax(); z++) {
            imgMask.setZ(z+1);
            ImageProcessor ip = imgMask.getProcessor();
            ByteProcessor bp = new ByteProcessor(ip, true);
            Object3D_IJUtils.draw(obj, bp, z, 255);
            ImagePlus maskPlus = new ImagePlus("mask " + z, bp);
            maskPlus.getProcessor().setAutoThreshold(AutoThresholder.Method.Default, true);
            IJ.run(maskPlus, "Create Selection", "");
            Roi roi = maskPlus.getRoi();
            img.setZ(z+1);
            img.getProcessor().setColor(col);
            img.getProcessor().drawRoi(roi);
            img.updateAndDraw();
        }   
        closeImages(imgMask);
    }
    
    
}
