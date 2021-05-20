package PML;

/*
 * Find dots (PML) in nucleus
 * Measure integrated intensity, nb of dots per nucleus 
 * Author Philippe Mailly
 */



import ij.*;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import loci.common.services.ServiceFactory;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Object3D;
import ij.measure.Calibration;
import ij.plugin.frame.RoiManager;
import java.awt.Rectangle;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.Region;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


public class PML_LiveCells implements PlugIn {

    private final boolean canceled = false;
    private String imageDir = "";
    private String outDirResults = "";
    private BufferedWriter outPutResults;
    private Calibration cal = new Calibration();
    
    private PML_Tools pml = new PML_Tools();
    
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            if (!pml.checkInstalledModules()) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Images folder");
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            String fileExt = pml.findImageType(inDir);
            ArrayList<String> imageFiles = pml.findImages(imageDir, fileExt);
            if (imageFiles == null) {
                return;
            }
            
            // create output folder
            outDirResults = imageDir + "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }

            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            int nucIndex = 0;
            // Find channel names , calibration
            reader.setId(imageFiles.get(0));
            int[] channelIndex = new int[reader.getSizeC()];
            cal = pml.findImageCalib(meta);
            String[] chsName = pml.findChannels(imageFiles.get(0), meta, reader);
            channelIndex = pml.dialog(chsName);
            if (channelIndex == null)
                return;
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                int series = reader.getSeries();
                int width = meta.getPixelsSizeX(series).getValue();
                int height = meta.getPixelsSizeY(series).getValue();
                reader.setSeries(series);
                boolean roiFile = true;
                // Test if ROI file exist
                String roi_file  = (new File(imageDir+rootName+".zip").exists()) ? imageDir+rootName+".zip" :  ((new File(imageDir+rootName+".roi").exists()) ? imageDir+rootName+".roi" : null);
                if (roi_file == null) {
                    IJ.showStatus("No ROI file found !") ;
                    roiFile = false;
                }
                List<Roi> rois = new ArrayList<>();
                if (roiFile) {
                    // find rois
                    RoiManager rm = new RoiManager(false);
                    rm.runCommand("Open", roi_file);
                    rois = Arrays.asList(rm.getRoisAsArray());
                }
                // define roi as all image
                else {
                    rois.add(new Roi(0, 0, width, height));
                }
                
                // For each roi open cropped image
                for (Roi roi : rois) {
                    nucIndex++;
                    
                    // Write headers results for results file{
                    FileWriter fileResults = new FileWriter(outDirResults + rootName + "_Nucleus_" + nucIndex +"_results.xls", false);
                    outPutResults = new BufferedWriter(fileResults);
                    outPutResults.write("Time\tNucleus Volume\tPML dot number\tNucleus Diffuse IntDensity\tPML Mean dots IntDensity\tPML dots Mean Volume"
                            + "\tPML dots STD IntDensity\tPML dot STD Volume\tPML Sum Vol\tPML dot Mean center-center distance\tPML dot SD center-center distance\n");
                    outPutResults.flush();
                    
                    Rectangle rectRoi = roi.getBounds();
                    ImporterOptions options = new ImporterOptions();
                    options.setId(f);
                    options.setCrop(true);
                    options.setCropRegion(0, new Region(rectRoi.x, rectRoi.y, rectRoi.width, rectRoi.height));
                    // open Dapi channel
                    
                    options.setCBegin(0, channelIndex[0]);
                    options.setCEnd(0, channelIndex[0]);
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setQuiet(true);
                    
                    Objects3DPopulation nucPop = new Objects3DPopulation();
                    ArrayList<Objects3DPopulation> pmlPopList = new ArrayList<>();
                    
                    // Stack registration
                    ImagePlus imgNuc = BF.openImagePlus(options)[0];
                    ArrayList<Transformer> trans = new StackReg_Plus().stackRegister(pml.stackProj(imgNuc));
                    ArrayList<Double> pmlDiffusInt = new ArrayList<>();
                    ArrayList<DescriptiveStatistics> pmlInt = new ArrayList<>();
                    pml.closeImages(imgNuc);
                    //int time = reader.getSizeT();
                    int time = 10;
                    // for each time find nucleus, plml
                    ImagePlus[] imgDiffusArray = new ImagePlus[time];
                    for (int t = 0; t < time; t++) {
                        IJ.showStatus("Reading time " + t + "/"+time);
                        options.setTBegin(0, t);
                        options.setTEnd(0, t);
                        options.setCBegin(0, channelIndex[0]);
                        options.setCEnd(0, channelIndex[0]);
                        // Open nucleus channel
                        imgNuc = BF.openImagePlus(options)[0];
                        // apply image drift correction
                        if (t>0) trans.get(t-1).doTransformation(imgNuc);
                        // find nuc object
                        Object3D nucObj = pml.findnucleus(imgNuc);
                        nucPop.addObject(nucObj);
                        // Open pml channel
                        options.setCBegin(0, channelIndex[1]);
                        options.setCEnd(0, channelIndex[1]);
                        ImagePlus imgPML = BF.openImagePlus(options)[0];
                        Objects3DPopulation pmlPop = new Objects3DPopulation();
                        // apply image drift correction
                        if ( t>0) trans.get(t-1).doTransformation(imgPML);
                        if (pml.trackMate_Detector_Method.equals("DoG"))
                            pmlPop = pml.findDotsDoG(imgPML, nucObj);
                        else
                            pmlPop = pml.findDotsLoG(imgPML, nucObj);
                        //System.out.println(pmlPop.getNbObjects()+" pml found");
                        pmlPopList.add(pmlPop);
                        pmlDiffusInt.add(pml.pmlDiffus(pmlPop, nucObj, imgPML));
                        pmlInt.add(pml.getPMLIntensity(pmlPop, imgPML));
                        imgDiffusArray[t] = imgPML;
                    }

                    // Save diffus image
                    pml.saveDiffusImage(pmlPopList, imgDiffusArray, outDirResults+rootName+"_Diffuse-"+nucIndex+".tif");
                    // Save images objects
                     pml.saveImageObjects(pmlPopList, nucPop, imgDiffusArray, outDirResults+rootName+"_Objects-"+nucIndex+".tif");
                     ImagePlus dotBin = pml.saveImagePMLs(nucPop, imgDiffusArray, outDirResults+rootName+"_PMLs-"+nucIndex+".tif");
                    
                     
                    
                    double meanVol = 0.0;
                    // find parameters
                    for (int i = 0; i < nucPop.getNbObjects(); i++) {
                        Object3D nucObj = nucPop.getObject(i);
                        double nucVol = nucObj.getVolumeUnit();
                        Objects3DPopulation pmlPop = pmlPopList.get(i);
                        int pmlDots = pmlPop.getNbObjects();
                        double pmlVolMean = pml.getPMLVolume(pmlPop).getMean();
                        meanVol += pmlVolMean;        
                       double pmlVolStd = pml.getPMLVolume(pmlPop).getStandardDeviation();
                        double pmlVolTotal = pml.getPMLVolume(pmlPop).getSum();
                        double minDistCenterMean = pmlPop.distancesAllClosestCenter().getMean(); 
                        double minDistCenterSD = pmlPop.distancesAllClosestCenter().getStdDev();
                        outPutResults.write(i+"\t"+nucVol+"\t"+pmlPop.getNbObjects()+"\t"+pmlDiffusInt.get(i)+"\t"+pmlInt.get(i).getMean()+"\t"
                                +pmlInt.get(i).getStandardDeviation()+"\t"+pmlVolMean+"\t"+pmlVolStd+"\t"+pmlVolTotal+"\t"+minDistCenterMean+"\t"+minDistCenterSD+"\n");
                        outPutResults.flush();
                    }                    
                    // Do Tracking
                    IJ.showStatus("Track PMLs");
                    TrackMater track = new TrackMater();
                    track.run(dotBin, outDirResults+rootName+"nuc_"+nucIndex+"_trackmateSaved.xml", outDirResults+rootName+"nuc_"+nucIndex+"_trackmateExport.xml", outDirResults+rootName+"nuc_"+nucIndex+"_trackmateSpotsStats.csv", outDirResults, rootName+"_PMLs-"+nucIndex+".tif");
                    }
            }
            //}
            IJ.showStatus("Process done"); 
        } catch (DependencyException ex) {
            Logger.getLogger(PML_LiveCells.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ServiceException ex) {
            Logger.getLogger(PML_LiveCells.class.getName()).log(Level.SEVERE, null, ex);
        } catch (FormatException ex) {
            Logger.getLogger(PML_LiveCells.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(PML_LiveCells.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}