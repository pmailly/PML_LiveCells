package PML;


import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.measure.*;

/** Select focused slices from a Z stack. Based on the autofocus algorithm "Normalized Variance" (Groen et al., 1985; Yeo et
 * al., 1993). However, the images are first treated by a sobel edge filter. This provided a better result for fluorescent bead images.
 * Code modified from the "Select Frames With Best Edges" plugin from Jennifer West (http://www.umanitoba.ca/faculties/science/astronomy/jwest/plugins.html)
 * First version 2009-9-27
 * Second version 2010-11-27
 * Third version 2011-2-15
 * Forth version 2011-3-2
 * Fifth version 2020-2-17 output to result table
 * By TSENG Qingzong; qztseng at gmail.com
 */
 
public class Find_focused_slices implements PlugInFilter, Measurements {

    boolean abort = false;
    double percent = 30;
    double vThr = 1.00;
    boolean consecutive = true;
    boolean edge = false;


    public void run(ImagePlus imp) {
        ImageStack stack = imp.getStack();
        int width = imp.getWidth();
        int height = imp.getHeight();
        String name = imp.getTitle();
        ImageStack stack2 = new ImageStack(width, height, imp.getProcessor().getColorModel());

        int size = stack.getSize();
        if (size == 1){
        	IJ.error("Stack required.");
            return;
        }

        double vMax = 0;
        double[] varA = new double[size];
        
        for (int slice = 1; slice <= size; slice++) {
            imp.setSlice(slice);
            ImageProcessor ip = imp.getProcessor();
            varA[slice - 1] = calVar(ip);
            if (varA[slice - 1] > vMax) {
                vMax = varA[slice - 1];
            }

        }
        if (vMax < vThr) {
            IJ.error("All slices are below the variance threshold value");
            return;
        }
        for (int slice = 1; slice <= size; slice++) {
            imp.setSlice(slice);
            ImageProcessor ip = imp.getProcessor();
            ip.resetRoi();
            if (varA[slice - 1] / vMax >= percent / 100 && varA[slice - 1] > vThr) {
                ip = ip.crop();
                stack2.addSlice("",ip,slice-1);                
            }
            else {
                ip.setColor(Color.black);
                ip.fill();
                stack2.addSlice("",ip,slice-1);
            }
        }
		
        ImagePlus focusstack = imp.createImagePlus();
        focusstack.setStack("Focused slices of " + name + "_" + percent + "%", stack2);
        focusstack.setCalibration(imp.getCalibration());
        //focusstack.show();
        //new WaitForUserDialog("focused").show();
    }

    double calVar(ImageProcessor ip) {

        double variance = 0;
        int W = ip.getWidth();
        int H = ip.getHeight();

        Rectangle r = ip.getRoi();
        if (r == null) {
            r.x = 0;
            r.y = 0;
            r.height = H;
            r.width = W;
        }
        ImageProcessor edged = ip.duplicate();
        if (edge) edged.findEdges();
        double mean = ImageStatistics.getStatistics(edged, MEAN, null).mean;
        double a = 0;
        for (int y = r.y; y < (r.y + r.height); y++) {
            for (int x = r.x; x < (r.x + r.width); x++) {
                a += Math.pow(edged.getPixel(x, y) - mean, 2);
            }
        }
        variance = (1 / (W * H * mean)) * a;
        return variance;

    }

    @Override
    public void run(ImageProcessor arg0) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

    @Override
    public int setup(String arg0, ImagePlus arg1) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

}
