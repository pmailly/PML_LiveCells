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
import ij.io.FileSaver;
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
            
            // use bioformat to open  images or not
            boolean bioformat = false;
            // Find channel names , calibration
            reader.setId(imageFiles.get(0));
            cal = pml.findImageCalib(meta);
            String[] chsName = pml.findChannels(imageFiles.get(0), meta, reader, bioformat);
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
                if (bioformat){
                    options.setId(f);
                    options.setCrop(false);
                    // open Dapi channel
                    options.setCBegin(0, channelIndex[0]);
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setCEnd(0, channelIndex[0]);
                    options.setQuiet(true);
                }
                
                int time = reader.getSizeT();
                time = 5;
                ArrayList<Objects3DPopulation> nucPops = new ArrayList<>();
                
                ImagePlus[] imgWholeArray = new ImagePlus[time];
                int nz = 0;
                
                // Look for nuclei
                 for (int t = 0; t < time; t++) {
                        IJ.showStatus("Reading time " + t + "/"+time);
                        ImagePlus imgNuc;
                        if (bioformat){
                            options.setTBegin(0, t);
                            options.setTEnd(0, t);
                            options.setCBegin(0, channelIndex[0]);
                            options.setCEnd(0, channelIndex[0]);
                           //Open nucleus channel
                            imgNuc = BF.openImagePlus(options)[0];
                        } else {
                            // faster ??
                            imgNuc = IJ.openImage(inDir+"/"+rootName+"_"+chsName[channelIndex[0]]+"_t"+(t+1)+".TIF");    
                        }
                       imgNuc.setCalibration(cal);          
                       imgWholeArray[t] = imgNuc.duplicate();
                       IJ.run(imgWholeArray[t], "8-bit", ""); // for smaller file size, only used for drawing
                       IJ.run(imgWholeArray[t], "Multiply...", "value=0.25 stack"); // put to low intensity to not interfere with drawings
                       nz = imgNuc.getNSlices();
                       Objects3DPopulation nucPop = pml.stardistNucleiPop(imgNuc);
                       nucPops.add(nucPop);    
                       pml.closeImages(imgNuc);
                 }
                
                 // "Track" nuclei: sorted by time to sorted by nuclei index
                 ArrayList<Objects3DPopulation> nuclei = new ArrayList<>();
                 Objects3DPopulation init = nucPops.get(0); // population at t=0
                 System.out.println("Found init "+init.getNbObjects()+" nuclei");
                 // find one nuclei in pop
                 for (int i = 0; i < init.getNbObjects(); i++) {
                        Object3D nucObj = init.getObject(i);
                        Objects3DPopulation one = pml.trackNucleus(nucObj, nucPops);
                        if ( one.getNbObjects() >= Math.min(time*0.25,20) ) {
                            nuclei.add(one);
                        }
                }
                System.out.println("Keep "+nuclei.size()+" nuclei");
                 
                 
                 // Work on each nuclei: crop, get pml, track...
                 for (int nucIndex=0; nucIndex<nuclei.size(); nucIndex++) {
                     
                     Objects3DPopulation nuc = nuclei.get(nucIndex);
                     int[] roilim = pml.getBoundingBoxXY(nuc);
                     // Crop a little larger than found nuclei
                     int extend = 5;
                     if ((roilim[0]-extend)>0) roilim[0] -= extend;
                     else roilim[0] = 0;
                     roilim[1] += extend;
                     if ((roilim[2]-extend)>0) roilim[2] -= extend;
                     else roilim[2] = 0;
                     roilim[3] += extend;
                     // open PML cropped
                     Roi croproi = new Roi(roilim[0], roilim[2], roilim[1]-roilim[0], roilim[3]-roilim[2]);
                     // move nucleus coordinates to ROI coordinates
                     pml.translateToRoi(nuc, roilim);
                          
                      // nb of times the nuclei is found, not necessarily until the end
                     int nuctime = nuc.getNbObjects();
                     ImagePlus[] imgDiffusArray = new ImagePlus[nuctime];
                     ArrayList<Objects3DPopulation> pmlPopList = new ArrayList<>();
                     ArrayList<Double> pmlDiffusInt = new ArrayList<>();
                     ArrayList<DescriptiveStatistics> pmlInt = new ArrayList<>();
                    
                     // Draw nucleus hyperstack (z, time)
                    ImagePlus unalignednucleus = pml.drawNucleusStack(nuc, croproi, time, nz);
                    // Get transformations to do to align stack
                    ArrayList<Transformer> trans = new StackReg_Plus().stackRegister(pml.stackProj(unalignednucleus));
                    ImagePlus alignednucleus = pml.alignStack(trans, unalignednucleus, true);
                    pml.alignPop(alignednucleus, nuc);
                    
                    
                    for (int t = 0; t < nuctime; t++) { 
                        ImagePlus imgPML;
                          // Open pml channel
                        if (bioformat){ 
                            options.setCBegin(0, channelIndex[1]);
                            options.setCEnd(0, channelIndex[1]);
                            options.setTBegin(0, t);
                            options.setTEnd(0, t);
                            imgPML = BF.openImagePlus(options)[0];
                        } else {
                            imgPML = IJ.openImage(inDir+"/"+rootName+"_"+chsName[channelIndex[1]]+"_t"+(t+1)+".TIF");
                        }
                            imgPML.setCalibration(cal);
                              // crop to ROI
                            imgPML.setRoi(croproi);
                            IJ.run(imgPML, "Crop", ""); 
                            IJ.run(imgPML, "Select None", ""); 
                            
                            imgPML.setCalibration(cal);
                            Objects3DPopulation pmlPop = new Objects3DPopulation();
                            Transformer curtrans = null;
                            if (t>0) curtrans = trans.get(t-1);
                           
                            if (pml.trackMate_Detector_Method.equals("DoG"))   
                                pmlPop = pml.findDotsAlign(imgPML, alignednucleus, curtrans, 1);
                            else
                                pmlPop = pml.findDotsAlign(imgPML, alignednucleus, curtrans, 0);
                           
                            pmlPopList.add(pmlPop);
                            if (t>0) (trans.get(t-1)).doTransformation(imgPML);
                            
                            pmlDiffusInt.add(pml.pmlDiffus(pmlPop, nuc.getObject(t), imgPML));
                            pmlInt.add(pml.getPMLIntensity(pmlPop, imgPML));
                            imgDiffusArray[t] = imgPML.duplicate();
                            pml.closeImages(imgPML);    
                    }
                    
                     // Save images objects
                    // Align populations
                    ImagePlus dotBin = pml.drawOneNucleiWithPML(pmlPopList, alignednucleus, outDirResults+rootName, nucIndex);
                    //pml.saveDiffusImage(pmlPopList, imgDiffusArray, outDirResults+rootName+"_Diffuse-"+nucIndex+".tif");
 
                    // Write headers results for results file{
                    FileWriter fileResults = new FileWriter(outDirResults + rootName + "_Nucleus_" + nucIndex +"_results.xls", false);
                    outPutResults = new BufferedWriter(fileResults);
                    outPutResults.write("Time\tNucleus Volume\tPML dot number\tNucleus Diffuse IntDensity\tPML Mean dots IntDensity\tPML dots Mean Volume"
                            + "\tPML dots STD IntDensity\tPML dot STD Volume\tPML Sum Vol\n");
                    outPutResults.flush();
                    
                    double meanVol = 0.0;
                    // find parameters
                    for (int i = 0; i < nuctime; i++) {
                        Object3D nucObj = nuc.getObject(i);
                        double nucVol = nucObj.getVolumeUnit();
                        Objects3DPopulation pmlPop = pmlPopList.get(i);
                        int pmlDots = pmlPop.getNbObjects();
                        double pmlVolMean = pml.getPMLVolume(pmlPop).getMean();
                        meanVol += pmlVolMean;        
                        double pmlVolStd = pml.getPMLVolume(pmlPop).getStandardDeviation();
                        double pmlVolTotal = pml.getPMLVolume(pmlPop).getSum();
                        outPutResults.write((i+1)+"\t"+nucVol+"\t"+pmlPop.getNbObjects()+"\t"+pmlDiffusInt.get(i)+"\t"+pmlInt.get(i).getMean()+"\t"
                                +pmlInt.get(i).getStandardDeviation()+"\t"+pmlVolMean+"\t"+pmlVolStd+"\t"+pmlVolTotal+"\n"); 
                        outPutResults.flush();
                    }                    
                    // Do Tracking
                    IJ.showStatus("Track PMLs");
                    TrackMater track = new TrackMater();
                    String resName = outDirResults+rootName+"nuc_"+nucIndex;
                    //track.run(dotBin, outDirResults, rootName+"_PMLs-"+nucIndex+".tif", pml.radius, pml.threshold, pml.trackMate_Detector_Method, 1.0, pml.merging_dist, true, true);
                    
                    boolean success = track.trackmateObjects(dotBin, pmlPopList, outDirResults, rootName+"_NucleusPMLs-"+nucIndex+".tif", pml.radius, 1.0, pml.merging_dist);
                    if (success) track.saveResults(resName+"_trackmateSaved.xml", resName+"_trackmateExport.xml", resName+"_trackmateSpotsStats.csv", resName+"_trackmateTrackStats.csv", outDirResults, rootName+"_PMLs-"+nucIndex+".tif");           
   
                    // Attention, This retranslate the populations to whole image position
                    pml.drawOnWholeImage(pmlPopList, nuc, imgWholeArray, croproi, nucIndex);
                 }
                 
                pml.saveWholeImage(imgWholeArray, outDirResults+rootName+"_Objects.tif");
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
