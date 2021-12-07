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
import java.util.concurrent.atomic.AtomicInteger;
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
import mcib3d.utils.ThreadUtil;

public class PML_StarDist_Parallel implements PlugIn {

    private final boolean canceled = false;
    private String imageDir = "";
    private String outDirResults = "";
    private Calibration cal = new Calibration();
    
    private PML_Tools pml = new PML_Tools();
    // use bioformat to open  images or not
    private boolean bioformat = false;
    private int width;
    private int height;
    private String rootName;
    private File inDir;
    private int[] channelIndex;
    private String f;
    private String[] chsName;
    
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
            inDir = new File(imageDir);
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
            
            
            // Find channel names , calibration
            reader.setId(imageFiles.get(0));
            cal = pml.findImageCalib(meta);
            chsName = pml.findChannels(imageFiles.get(0), meta, reader, bioformat);
            //System.out.println(chsName[0]);
            channelIndex = pml.dialog(chsName);
            cal = pml.getCalib();
            if (channelIndex == null)
                return;
            for (String cf : imageFiles) {
                f = cf;
                rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                int series = reader.getSeries();
                width = meta.getPixelsSizeX(series).getValue();
                height = meta.getPixelsSizeY(series).getValue();
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
                
                final int time = reader.getSizeT();
                //final int time = 5;
                Objects3DPopulation[] nucPops = new Objects3DPopulation[time];
                if (pml.verbose) IJ.log("\\Clear");
                
                int nz = reader.getSizeZ();
                         
                // try parallelize
                final AtomicInteger ai = new AtomicInteger(0);
                int ncpu = (int) Math.ceil(ThreadUtil.getNbCpus()*0.9);
                Thread[] threads = ThreadUtil.createThreadArray(ncpu);
                final int neach = (int) Math.ceil((double)time/(double)ncpu);
                // Look for nuclei
                for (int iThread=0; iThread<threads.length; iThread++) {
                    threads[iThread] = new Thread(){
                        public void run(){
                            for (int k=ai.getAndIncrement(); k<ncpu; k=ai.getAndIncrement()) {
                                for (int t = neach*k; ((t<(neach*(k+1))&&(t<time))); t++) {
                                    
                                   try{ 
                                    //IJ.showStatus("Reading time " + t + "/"+time);
                                    //System.out.println("Doing parallel "+t);
                                    if (pml.verbose) IJ.log("Doing time " + t + "/"+time);
                                    ImagePlus imgNuc;
                                    if (bioformat){
                                        ImporterOptions newopt = new ImporterOptions();
                                        newopt.setId(f);
                                        newopt.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                                        newopt.setQuiet(true);
                                        newopt.setCBegin(0, channelIndex[0]);
                                        newopt.setCEnd(0, channelIndex[0]);
                                        newopt.setTBegin(0, t);
                                        newopt.setTEnd(0, t);
                                        //Open nucleus channel
                                        imgNuc = BF.openImagePlus(newopt)[0];
                                    } else {
                                        String imageName = inDir+File.separator+rootName+"_"+chsName[channelIndex[0]]+"_t"+(t+1)+".TIF";
                                        if (pml.multiPos)
                                            imageName = inDir+File.separator+rootName+"_"+chsName[channelIndex[0]]+"_s1_t"+(t+1)+".TIF";
                                        imgNuc = IJ.openImage(imageName);    
                                    }
                                   imgNuc.setCalibration(cal);          
                                  
                                   //nz = imgNuc.getNSlices();
                                   //Thread.sleep((long)(Math.random() * 10000));
                                   Objects3DPopulation nucPop = pml.stardistNucleiPop(imgNuc);
                                   nucPops[t] = nucPop;  
                                   if (pml.verbose) IJ.log("Found "+nucPop.getNbObjects()+" nuclei");
                                   pml.closeImages(imgNuc);
                                    }
                                   catch (Exception e) {Logger.getLogger(PML_LiveCells.class.getName()).log(Level.SEVERE, null, e);}
                                    }
                              }
                        }
                     };
              }
              ThreadUtil.startAndJoin(threads);
                 
                
                // Do each nuclei that is kept long enough
                 // "Track" nuclei: sorted by time to sorted by nuclei index
                 Objects3DPopulation init = nucPops[0]; // population at t=0
                  if (pml.verbose) IJ.log(" ********************************************** \n");
                   if (pml.verbose) IJ.log(" ********************************************** \n");
                    if (pml.verbose) IJ.log(" ********************************************** \n");
                     if (pml.verbose) IJ.log(" ********** Nucleus segmentation done ************ \n");
                 
                 if (pml.verbose) IJ.log("Found init "+init.getNbObjects()+" nuclei");
                ArrayList<Objects3DPopulation> nuclei = new ArrayList<Objects3DPopulation>();
                 for (int i = 0; i < init.getNbObjects(); i++) {
                        Object3D nucObj = init.getObject(i);
                        Objects3DPopulation nuc = trackNucleus(nucObj, nucPops);
                        
                        // nb of times the nuclei is found, not necessarily until the end
                        int nuctime = nuc.getNbObjects();
                            
                        // Keep the nuclei and work on it
                        if ( nuc!=null && nuctime >= Math.min(time*0.25,20) ) {
                               nuclei.add(nuc);
                        }
                 }
                 
                 int nnuclei = nuclei.size();
                 if (pml.verbose) IJ.log("Keep "+nnuclei+" nuclei to analyze \n");
                 if (pml.verbose) IJ.log(" ********************************************** \n");
                 if (pml.verbose) IJ.log(" ********************************************** \n");
                 if (pml.verbose) IJ.log(" ********************************************** \n");                
// new WaitForUserDialog("nuclei ready "+nnuclei).show();
                // try parallelize
                final AtomicInteger rai = new AtomicInteger(0);
                final int neac = (int) Math.ceil((double)nnuclei/(double)ncpu);
                // Look for nuclei
                for (int iThread=0; iThread<threads.length; iThread++) {
                    threads[iThread] = new Thread(){
                        public void run(){
                            for (int k=rai.getAndIncrement(); k<ncpu; k=rai.getAndIncrement()) {
                                for (int nucIndex = neac*k; ((nucIndex<(neac*(k+1))&&(nucIndex<nnuclei))); nucIndex++) {
                                    Objects3DPopulation nuc = nuclei.get(nucIndex);
                                    doOneNucleus(nucIndex, nuc, time, nz);
                                }
                            }
                        }
                    };
                           
                 }
                ThreadUtil.startAndJoin(threads);
              }
            IJ.showStatus("Process done");
            if (pml.verbose) IJ.log("Process done !");
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
    
       /** Looks for the closest nuclei at each time (track it). Loose it if too far from given distance, do it by association instead ? */ 
    public Objects3DPopulation trackNucleus( Object3D obj, Objects3DPopulation[] pop) {
        Objects3DPopulation nucl = new Objects3DPopulation();
        if (pop.length>1){
        for (int i=0; i<pop.length; i++) {
            if (pop[i]==null) return nucl;
            Object3D closest = (pop[i]).closestCenter(obj.getCenterAsPoint());
            // threshold distance to loose the nuclei (not aligned image so can move)
            if (obj.distCenterUnit(closest) > 4) {
                return nucl;
            }
            // within distance, continue
            nucl.addObject(closest);
            obj = closest;
        }
        return nucl;
        }
        return null;
    }
    
    public void doOneNucleus(int nucIndex, Objects3DPopulation nuc, int time, int nz){
          try{             

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

                    int nuctime = nuc.getNbObjects();
                Objects3DPopulation[] pmlPopList = new Objects3DPopulation[nuctime];
            
                 if (pml.verbose) IJ.log("Working on nuclei "+nucIndex+" timefinal "+nuctime);
                ImagePlus[] dotBins = new ImagePlus[nuctime];  
                
                // Draw nucleus hyperstack (z, time)
                ImagePlus unalignednucleus = pml.drawNucleusStack(nuc, croproi, nuctime, nz);

                // Get transformations to do to align stack
                ArrayList<Transformer> trans = new StackReg_Plus().stackRegister(pml.stackProj(unalignednucleus), nucIndex);

                ImagePlus alignednucleus = pml.alignStack(trans, unalignednucleus, true, nucIndex);
                pml.closeImages(unalignednucleus);
                pml.alignPop(alignednucleus, nuc);

                // Write headers results for results file{
                FileWriter fileResults = new FileWriter(outDirResults + rootName + "_Nucleus_" + nucIndex +"_results.xls", false);
                BufferedWriter outPutResults = new BufferedWriter(fileResults);
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
                        ImporterOptions newopt = new ImporterOptions();
                        newopt.setId(f);
                        newopt.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                        newopt.setQuiet(true);
                        newopt.setCBegin(0, channelIndex[1]);
                        newopt.setCEnd(0, channelIndex[1]);
                        newopt.setTBegin(0, t);
                        newopt.setTEnd(0, t);
                        imgPML = BF.openImagePlus(newopt)[0];
                    } else {
                        String imageName = inDir+File.separator+rootName+"_"+chsName[channelIndex[1]]+"_t"+(t+1)+".TIF";
                        if (pml.multiPos)
                            imageName = inDir+File.separator+rootName+"_"+chsName[channelIndex[1]]+"_s1_t"+(t+1)+".TIF";
                        imgPML = IJ.openImage(imageName);
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
                        pmlPop = pml.findDotsStarDist(imgPML, anucleus, curtrans, nucIndex, false);

                     // draw current time point
                     dotBins[t] = pml.drawOneNucleiWithPMLOneTime(pmlPop, anucleus); 
                    pmlPopList[t] = pmlPop;
                    if (t>0) (trans.get(t-1)).doTransformation(imgPML, false, nucIndex);

                    int pmlNbDots = pmlPop.getNbObjects();
                    String pmlVols = pml.getPMLVolumesAsString(pmlPop);
                    String pmlInt = pml.getPMLIntensity(pmlPop, imgPML);
                    if ( pml.savePMLImg) pmlImgs[t] = imgPML.duplicate();
                    pml.closeImages(imgPML);
                    double nucleusVolume = (nuc.getObject(t)).getVolumeUnit();
                    outPutResults.write((t+1)+"\t"+nucleusVolume+"\t"+pmlNbDots+"\t"+pmlInt+"\t"+pmlVols+"\n"); 
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
                        pml.closeImages(pmls);
                    }

                    // Do Tracking
                    IJ.showStatus("Track PMLs");
                    TrackMater track = new TrackMater();
                    String resName = outDirResults+rootName+"nuc_"+nucIndex;

                    if (pml.verbose) IJ.log("Track PMLs for nuclei "+nucIndex);               
                    boolean success = track.trackmateObjects(dotBin, pmlPopList, outDirResults, rootName+"_NucleusPMLs-"+nucIndex+".tif", pml.radius, 1.0, pml.merging_dist);
                    if (success) track.saveResults(resName+"_trackmateSaved.xml", resName+"_trackmateExport.xml", resName+"_trackmateSpotsStats.csv", resName+"_trackmateTrackStats.csv", outDirResults, rootName+"_PMLs-"+nucIndex+".tif");           
                    if (pml.verbose) IJ.log("Tracking PMLs for nuclei "+nucIndex+" finised: "+success);
                    pml.closeImages(dotBin);

                }
                catch (Exception e) {Logger.getLogger(PML_LiveCells.class.getName()).log(Level.SEVERE, null, e);}
             }
}
