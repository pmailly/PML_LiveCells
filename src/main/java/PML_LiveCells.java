/*
 * Find dots (PML) in nucleus
 * Measure integrated intensity, nb of dots per nucleus 
 * Author Philippe Mailly
 */




import PMLTools.PML_Tools;
import static PMLTools.PML_Tools.closeImages;
import static PMLTools.PML_Tools.stackRegister;
import static PMLTools.PML_Tools.dialog;
import static PMLTools.PML_Tools.driftCorrection;
import static PMLTools.PML_Tools.findDots;
import static PMLTools.PML_Tools.findImageCalib;
import static PMLTools.PML_Tools.findImages;
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
import mcib3d.geom.Object3DVoxels;
import ij.gui.PolygonRoi;
import ij.measure.Calibration;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.Region;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import org.apache.commons.io.FilenameUtils;


public class PML_LiveCells implements PlugIn {

    private final boolean canceled = false;
    private String imageDir = "";
    public static String outDirResults = "";
    private BufferedWriter outPutResults;
    private Calibration cal = new Calibration();

 
  
    /**
     * Find convex hull of PML dots
     * @param img
     * @return nucleusObj
     */
    private Object3D find_nucleus(ImagePlus img) {
        MaximumFinder max = new MaximumFinder();
        for (int z = 1; z <= img.getNSlices(); z++ ) {
            ImageProcessor ip = img.getStack().getProcessor(z);
            ip.blurGaussian(2);
            Polygon polyMax = max.getMaxima(ip, 10, true);
            if (polyMax.npoints > 1) {
                Polygon convexhull = new PolygonRoi(polyMax, Roi.POLYGON).getConvexHull();
                ip.setRoi(convexhull);
                ip.setColor(Color.black);
                ip.fillPolygon(convexhull);
                for (int i =0; i < 10; i++)
                    ip.dilate();
            }
        }
        img.updateAndDraw();
        Object3DVoxels nucObj = new Object3DVoxels(img, 255);
        return(nucObj);
    }
           
    
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
            imageDir = dialog();
            if (imageDir == null) {
                return;
            }
            String fileExt = "nd";
            File inDir = new File(imageDir);
            ArrayList<String> imageFiles = findImages(imageDir, fileExt);
            if (imageFiles == null) {
                return;
            }
            // create output folder
            outDirResults = inDir + File.separator+ "Results"+ File.separator;
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
            int series = 0;
            int imageNum = 0, nucIndex = 0;
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                reader.setSeries(series);
                if (imageNum == 0) {
                    cal = findImageCalib(meta);
                }
                
                // Find ROI file
                String roi_file = imageDir+rootName+".zip";
                if (!new File(roi_file).exists()) {
                    IJ.showStatus("No ROI file found !") ;
                    return;
                }
                else {
                    roi_file = imageDir+rootName+".roi";
                    if (!new File(roi_file).exists()) {
                        IJ.showStatus("No ROI file found !") ;
                        return;
                    }
                }
                // Write headers results for results file
                if (imageNum == 0) {
                    FileWriter fileResults = new FileWriter(outDirResults + rootName + "_Nucleus_" + nucIndex +"_results.xls", false);
                    outPutResults = new BufferedWriter(fileResults);
                    outPutResults.write("Time\tNucleus Volume\tPML dot number\tNucleus Diffuse IntDensity\tPML Mean dot IntDensity\tPML dot SD IntDensity"
                            + "\tPML dot Min IntDensity\tPML dot Max IntDensity\tPML dot Mean Volume\tPML dot SD Volume\tPML Min Vol\tPML Max Vol\tPML Sum Vol"
                            + "\tPML dot Mean center-center distance\tPML dot SD center-center distance\n");
                    outPutResults.flush();
                }
                imageNum++;
                // find rois
                RoiManager rm = new RoiManager(false);
                rm.runCommand("Open", roi_file);
                Roi[] rois = rm.getRoisAsArray();
                
                // For each roi open cropped image
                for (Roi roi : rois) {
                    nucIndex++;

                    Rectangle rectRoi = roi.getBounds();
                    ImporterOptions options = new ImporterOptions();
                    options.setId(f);
                    options.setCrop(true);
                    options.setCropRegion(0, new Region(rectRoi.x, rectRoi.y, rectRoi.width, rectRoi.height));
                    options.setCBegin(0, 0);
                    options.setCEnd(0, 1);
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setQuiet(true);
                    
                    Objects3DPopulation nucPop = new Objects3DPopulation();
                    ArrayList<Objects3DPopulation> plmPop = new ArrayList<>();
                    
                    // Stack registration
                    ImagePlus imgNuc = BF.openImagePlus(options)[0];
                    ArrayList<double[]> matReg = stackRegister(imgNuc);
                    closeImages(imgNuc);
                    // for each time find nucleus, plml
                    for (int t = 0; t < reader.getSizeT(); t++) {
                        options.setTBegin(0, t);
                        options.setTEnd(0, t);
                        // Open nucleus channel
                        imgNuc = BF.openImagePlus(options)[0];
                        // apply image drift correction
                        driftCorrection(imgNuc, matReg);
                        // find nuc object
                        nucPop.addObject(PML_Tools.findnucleus(imgNuc));
                        // Open pml channel
                        ImagePlus imgPML = BF.openImagePlus(options)[1];
                        // apply image drift correction
                        driftCorrection(imgPML, matReg);
                        plmPop.add(findDots(imgPML, "fileName"));
                    }
                }
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