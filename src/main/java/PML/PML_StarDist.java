package PML;

/*
 * Find dots (PML) in nucleus
 * Measure integrated intensity, nb of dots per nucleus 
 * Use StarDist to find nucleus 
 * Author Orion
 */



import ij.*;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
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
import ij.plugin.Concatenator;
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


public class PML_StarDist implements PlugIn {

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
            cal = pml.findImageCalib(meta);
            String[] chsName = pml.findChannels(imageFiles.get(0), meta, reader);
            //System.out.println(chsName[0]);
            int[] channelIndex = pml.dialog(chsName);
            cal = pml.getCalib();
            if (channelIndex == null)
                return;
            for (String f : imageFiles) {
               String rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                int series = reader.getSeries();
                int width = meta.getPixelsSizeX(series).getValue();
                int height = meta.getPixelsSizeY(series).getValue();
                reader.setSeries(series);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setCrop(false);
                // open Dapi channel
                options.setCBegin(0, channelIndex[0]);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setCEnd(0, channelIndex[0]);
                options.setQuiet(true);
                
                int time = reader.getSizeT();
                time = 5;
                ImagePlus[] nucleiMaskArray = new ImagePlus[time];
                ArrayList<Objects3DPopulation> nucPops = new ArrayList<>();
                
                // Look for nuclei
                 for (int t = 0; t < time; t++) {
                        IJ.showStatus("Reading time " + t + "/"+time);
                        options.setTBegin(0, t);
                        options.setTEnd(0, t);
                        options.setCBegin(0, channelIndex[0]);
                        options.setCEnd(0, channelIndex[0]);
                        // Open nucleus channel
                        ImagePlus imgNuc = BF.openImagePlus(options)[0];
                        imgNuc.setCalibration(cal);          
                       // nucleiMaskArray[t] = pml.stardistNuclei(imgNuc);
                       Objects3DPopulation nucPop = pml.stardistNucleiPop(imgNuc);
                       nucPops.add(nucPop);
                 }
                
                 // "Track" nuclei: sorted by time to sorted by nuclei index
                 ArrayList<Objects3DPopulation> nuclei = new ArrayList<>();
                 Objects3DPopulation init = nucPops.get(0); // population at t=0
                 // find one nuclei in pop
                 for (int i = 0; i < init.getNbObjects(); i++) {
                        Object3D nucObj = init.getObject(i);
                        
                }
                //ImagePlus nucleiImg = new Concatenator().concatenate(nucleiMaskArray, false);
                //nucleiImg.show();
                
                
                //new WaitForUserDialog("show").show();
                    
                 //TrackMater track = new TrackMater();
                 // nuclei radius
                 //track.run(nucleiImg, outDirResults, rootName+"_Nuclei.tif", 10, 1, pml.trackMate_Detector_Method, 3, 2, false, false);
                 //track.getTracks();
                
                /** 
                    // Write headers results for results file{
                    FileWriter fileResults = new FileWriter(outDirResults + rootName + "_Nucleus_" + nucIndex +"_results.xls", false);
                    outPutResults = new BufferedWriter(fileResults);
                    outPutResults.write("Time\tNucleus Volume\tPML dot number\tNucleus Diffuse IntDensity\tPML Mean dots IntDensity\tPML dots Mean Volume"
                            + "\tPML dots STD IntDensity\tPML dot STD Volume\tPML Sum Vol\tPML dot Mean center-center distance\tPML dot SD center-center distance\n");
                    outPutResults.flush();
                    
                    Objects3DPopulation nucPop = new Objects3DPopulation();
                    ArrayList<Objects3DPopulation> pmlPopList = new ArrayList<>();
                    
                    ArrayList<Transformer> trans = new StackReg_Plus().stackRegister(pml.stackProj(imgNuc));
                    ArrayList<Double> pmlDiffusInt = new ArrayList<>();
                    ArrayList<DescriptiveStatistics> pmlInt = new ArrayList<>();
                    pml.closeImages(imgNuc);
                    //int time = reader.getSizeT();
                    int time = reader.getSizeT();
                    //time = 4;
                    // for each time find nucleus, plml
                    * 
                    ImagePlus[] imgDiffusArray = new ImagePlus[time];
                    * 
                        // apply image drift correction
                        Object3D nucObj = pml.findnucleus(imgNuc, trans , t);
                        if (nucObj == null) {
                            IJ.log("No nucleus found in Roi "+nucIndex+" at time "+(t+1)+"\n Skipping it \n");
                            break;
                        }
                        nucPop.addObject(nucObj);
                        // Open pml channel
                        options.setCBegin(0, channelIndex[1]);
                        options.setCEnd(0, channelIndex[1]);
                        ImagePlus imgPML = BF.openImagePlus(options)[0];
                        imgPML.setCalibration(cal);
                        Objects3DPopulation pmlPop = new Objects3DPopulation();
                        // apply image drift correction
                        if (t > 0)
                            trans.get(t - 1).doTransformation(imgPML);
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
                    
                    // Save images objects
                    pml.saveImageObjects(pmlPopList, nucPop, imgDiffusArray, outDirResults+rootName+"_Objects-"+nucIndex+".tif");
                    ImagePlus dotBin = pml.saveImagePMLs(pmlPopList, imgDiffusArray, outDirResults+rootName+"_PMLs-"+nucIndex+".tif");
                    // Save diffus image
                    pml.saveDiffusImage(pmlPopList, imgDiffusArray, outDirResults+rootName+"_Diffuse-"+nucIndex+".tif");
                        
                     
                    
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
                        //double minDistCenterMean = pmlPop.distancesAllClosestCenter().getMean(); 
                        //double minDistCenterSD = pmlPop.distancesAllClosestCenter().getStdDev();
                        outPutResults.write(i+"\t"+nucVol+"\t"+pmlPop.getNbObjects()+"\t"+pmlDiffusInt.get(i)+"\t"+pmlInt.get(i).getMean()+"\t"
                                +pmlInt.get(i).getStandardDeviation()+"\t"+pmlVolMean+"\t"+pmlVolStd+"\t"+pmlVolTotal+"\n"); //+minDistCenterMean+"\t"+minDistCenterSD+"\n");
                        outPutResults.flush();
                    }                    
                    // Do Tracking
                    IJ.showStatus("Track PMLs");
                    TrackMater track = new TrackMater();
                    String resName = outDirResults+rootName+"nuc_"+nucIndex;
                    track.run(dotBin, resName+"_trackmateSaved.xml", resName+"_trackmateExport.xml", resName+"_trackmateSpotsStats.csv", resName+"_trackmateTrackStats.csv", outDirResults, rootName+"_PMLs-"+nucIndex+".tif");
                    }*/
           }
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