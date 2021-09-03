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
import ij.plugin.SubHyperstackMaker;
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
                time = 4;
                ArrayList<Objects3DPopulation> nucPops = new ArrayList<>();
                if (pml.verbose) IJ.log("\\Clear");
                
                // to save full image
                ImagePlus[] imgWholeArray = null;
                if (pml.saveWhole) imgWholeArray = new ImagePlus[time];
                double resizeWhole = 300;  // diminue la taille pour pas enregistrer une image trop big
                double factorResize = 0;       
                int nz = 0;
                         
                // Look for nuclei
                for (int t = 0; t < time; t++) {
                        IJ.showStatus("Reading time " + t + "/"+time);
                        if (pml.verbose) IJ.log("Doing time " + t + "/"+time);
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
                       if (pml.saveWhole){
                        ImagePlus whole  = imgNuc.duplicate();
                        IJ.run(whole, "8-bit", ""); // for smaller file size, only used for drawing
                        IJ.run(whole, "Multiply...", "value=0 stack"); // empty image
                        factorResize = resizeWhole/whole.getWidth();
                        ImagePlus resized = whole.resize((int)resizeWhole, (int)(whole.getHeight()*factorResize), "bilinear");
                        imgWholeArray[t] = resized;
                        pml.closeImages(whole);
                       }
                       
                       nz = imgNuc.getNSlices();
                       Objects3DPopulation nucPop = pml.stardistNucleiPop(imgNuc);
                       nucPops.add(nucPop);  
                       if (pml.verbose) IJ.log("Found "+nucPop.getNbObjects()+" nuclei");
                       pml.closeImages(imgNuc);
                 }
                
                // Do each nuclei that is kept long enough
                 // "Track" nuclei: sorted by time to sorted by nuclei index
                 Objects3DPopulation init = nucPops.get(0); // population at t=0
                 if (pml.verbose) IJ.log("Found init "+init.getNbObjects()+" nuclei");
                 
                 int nucIndex = 0;                 
                 for (int i = 0; i < init.getNbObjects(); i++) {
                        Object3D nucObj = init.getObject(i);
                        Objects3DPopulation nuc = pml.trackNucleus(nucObj, nucPops);
                        
                        // nb of times the nuclei is found, not necessarily until the end
                        int nuctime = nuc.getNbObjects();
                            
                        // Keep the nuclei and work on it
                        if ( nuc!=null && nuctime >= Math.min(time*0.25,20) ) {
                            
                            // Work on each nuclei: crop, get pml, track...
                            if (pml.verbose) IJ.log("Working on nuclei "+nucIndex);
                            int[] roilim = pml.getBoundingBoxXY(nuc);
                            // Crop a little larger than found nuclei
                            int extend = 10;
                            if ((roilim[0]-extend)>=0) roilim[0] -= extend;
                            else roilim[0] = 0;
                            if ((roilim[1]+extend)<width) roilim[1] += extend;
                            else roilim[1] = (width-1);
                            if ((roilim[2]-extend)>=0) roilim[2] -= extend;
                            else roilim[2] = 0;
                            if ((roilim[3]+extend)<height) roilim[3] += extend;
                            else roilim[3] = (height-1);
                             // open PML cropped
                             Roi croproi = new Roi(roilim[0], roilim[2], roilim[1]-roilim[0], roilim[3]-roilim[2]);
                             // move nucleus coordinates to ROI coordinates
                             pml.translateToRoi(nuc, roilim);

                            ArrayList<Objects3DPopulation> pmlPopList = new ArrayList<>();
                            ImagePlus[] dotBins = new ImagePlus[nuctime];    

                            // Draw nucleus hyperstack (z, time)
                            ImagePlus unalignednucleus = pml.drawNucleusStack(nuc, croproi, time, nz);
                            // Get transformations to do to align stack
                            ArrayList<Transformer> trans = new StackReg_Plus().stackRegister(pml.stackProj(unalignednucleus), nucIndex);
                            ImagePlus alignednucleus = pml.alignStack(trans, unalignednucleus, true, nucIndex);
                            pml.closeImages(unalignednucleus);
                            pml.alignPop(alignednucleus, nuc);

                            // Write headers results for results file{
                            FileWriter fileResults = new FileWriter(outDirResults + rootName + "_Nucleus_" + nucIndex +"_results.xls", false);
                            outPutResults = new BufferedWriter(fileResults);
                            outPutResults.write("Time\tNucleus Volume\tPML dot number\tPML Mean dots IntDensity\tPML dots Mean Volume"
                                    + "\tPML dots STD IntDensity\tPML dot STD Volume\tPML Sum Vol\n");
                            outPutResults.flush();

                            SubHyperstackMaker sub = new SubHyperstackMaker();
                            ImagePlus[] pmlImgs = null;
                            if (pml.savePMLImg) pmlImgs = new ImagePlus[nuctime];
                            for (int t = 0; t < nuctime; t++) { 
                                 if (pml.verbose) IJ.log("Working on nuclei "+nucIndex+" time "+t+"∕"+nuctime);

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

                                    ImagePlus anucleus = sub.makeSubhyperstack(alignednucleus, "1-1", "1-"+alignednucleus.getNSlices(), (t+1)+"-"+(t+1));
                                    if (pml.trackMate_Detector_Method.equals("DoG"))   
                                        pmlPop = pml.findDotsAlign(imgPML, anucleus, curtrans, 1, nucIndex);
                                    if (pml.trackMate_Detector_Method.equals("LoG"))
                                        pmlPop = pml.findDotsAlign(imgPML, anucleus, curtrans, 0, nucIndex);
                                    if (pml.trackMate_Detector_Method.equals("StarDist"))
                                        pmlPop = pml.findDotsStarDist(imgPML, anucleus, curtrans, nucIndex);

                                     // draw current time point
                                     dotBins[t] = pml.drawOneNucleiWithPMLOneTime(pmlPop, anucleus); 
                                    pmlPopList.add(pmlPop);
                                    if (t>0) (trans.get(t-1)).doTransformation(imgPML, false, nucIndex);

                                    int pmlNbDots = pmlPop.getNbObjects();
                                    double[] pmlVols = pml.getPMLVolumes(pmlPop);
                                    double[] pmlInt = pml.getPMLIntensity(pmlPop, imgPML);
                                    if ( pml.savePMLImg) pmlImgs[t] = imgPML.duplicate();
                                    pml.closeImages(imgPML);
                                    double nucleusVolume = (nuc.getObject(t)).getVolumeUnit();
                                    outPutResults.write((t+1)+"\t"+nucleusVolume+"\t"+pmlNbDots+"\t"+pmlInt[0]+"\t"
                                        +pmlInt[1]+"\t"+pmlVols[0]+"\t"+pmlVols[1]+"\t"+pmlVols[2]+"\n"); 
                                    outPutResults.flush();
                            }

                             if (pml.verbose) IJ.log("Write image nuclei "+nucIndex);
                            // Save images objects
                            ImagePlus dotBin = new Concatenator().concatenate(dotBins, false);    
                            FileSaver nucPML = new FileSaver(dotBin);
                            nucPML.saveAsTiff(outDirResults+rootName+"_NucleusPMLs-"+nucIndex+".tif");
                            if (pml.savePMLImg){
                                ImagePlus pmls = new Concatenator().concatenate(pmlImgs, false);    
                                FileSaver imgPML = new FileSaver(pmls);
                                imgPML.saveAsTiff(outDirResults+rootName+"_PMLs-"+nucIndex+".tif");     
                            }

                            // Do Tracking
                            IJ.showStatus("Track PMLs");
                            TrackMater track = new TrackMater();
                            String resName = outDirResults+rootName+"nuc_"+nucIndex;

                            if (pml.verbose) IJ.log("Track PMLs for nuclei "+nucIndex);               
                            boolean success = track.trackmateObjects(dotBin, pmlPopList, outDirResults, rootName+"_NucleusPMLs-"+nucIndex+".tif", pml.radius, 1.0, pml.merging_dist);
                            if (success) track.saveResults(resName+"_trackmateSaved.xml", resName+"_trackmateExport.xml", resName+"_trackmateSpotsStats.csv", resName+"_trackmateTrackStats.csv", outDirResults, rootName+"_PMLs-"+nucIndex+".tif");           

                            if ( pml.saveWhole){
                                // Attention, This retranslate the populations to whole image position
                                pml.drawOnWholeImage(pmlPopList, nuc, imgWholeArray, croproi, nucIndex, factorResize);
                            }
                            nucIndex++;
                         }
                 }
                    if (pml.saveWhole) pml.saveWholeImage(imgWholeArray, outDirResults+rootName+"_ObjectsLowRes.tif");
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
