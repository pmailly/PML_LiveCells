/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package PML;

import ij.ImagePlus;
import ij.gui.WaitForUserDialog;
import ij.plugin.ContrastEnhancer;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Rectangle;

/**
 *
 * @author Gaelle Letort, Collège de France
 */
public class DoG extends GaussianBlur {
  private ImageProcessor ip1;
  private ImageProcessor ip2;
  private PlugInFilterRunner pfr;
  
  
  public static void stackDOG(ImagePlus imp, double sigma1, double sigma2) {
      //imp.show();
      for (int i=1; i<=imp.getNSlices(); i++) {
          imp.setSlice(i);
          ImageProcessor ip = imp.getProcessor();
          run(ip, sigma1, sigma2);
      }
     // new WaitForUserDialog("test").show();
  } 
  
  /**
   * Perform a Difference of Gaussians (filteredImage2 - filteredImage1) on the image. Sigma1 should
   * be greater than sigma2.
   *
   * @param ip the image
   * @param sigma1 the sigma 1
   * @param sigma2 the sigma 2
   */
  public static void run(ImageProcessor ip, double sigma1, double sigma2) {
    final ImageProcessor ip1 = (sigma1 > 0) ? duplicateProcessor(ip) : ip;
    final ImageProcessor ip2 = (sigma2 > 0) ? duplicateProcessor(ip) : ip;
    final DoG filter = new DoG();
    //filter.noProgress = true;
    filter.showProgress(false);
    filter.blurGaussian(ip1, sigma1);
    filter.blurGaussian(ip2, sigma2);
    differenceOfGaussians(ip, ip1, ip2);
  }

  /**
   * Subtract one image from the other (ip2 - ip1) and store in the result processor.
   *
   * @param resultIp the result ip
   * @param ip1 the image 1
   * @param ip2 the image 2
   */
  private static void differenceOfGaussians(ImageProcessor resultIp, ImageProcessor ip1,
      ImageProcessor ip2) {
    // Reuse the processor space
    FloatProcessor fp1 = null;
    FloatProcessor fp2 = null;
    FloatProcessor fp3 = null;

    final Rectangle roi = resultIp.getRoi();
    final int yTo = roi.y + roi.height;

    for (int i = 0; i < resultIp.getNChannels(); i++) {
      fp1 = ip1.toFloat(i, fp1);
      fp2 = ip2.toFloat(i, fp2);
      final float[] ff1 = (float[]) fp1.getPixels();
      final float[] ff2 = (float[]) fp2.getPixels();

      // If an ROI is present start with the original image
      if (roi.height != resultIp.getHeight() || roi.width != resultIp.getWidth()) {
        fp3 = resultIp.toFloat(i, fp3);
        final float[] ff3 = (float[]) fp3.getPixels();
        // Copy within the ROI
        for (int y = roi.y; y < yTo; y++) {
          int index = y * resultIp.getWidth() + roi.x;
          for (int x = 0; x < roi.width; x++, index++) {
            ff3[index] = ff2[index] - ff1[index];
          }
        }
      } else {
        fp3 = new FloatProcessor(fp1.getWidth(), fp2.getHeight());
        final float[] ff3 = (float[]) fp3.getPixels();
        for (int j = ff1.length; j-- > 0;) {
          ff3[j] = ff2[j] - ff1[j];
        }
      }

      if (Thread.currentThread().isInterrupted()) {
        return; // interruption for new parameters during preview?
      }
      resultIp.setPixels(i, fp3);
    }
  }

  /**
   * Perform a Gaussian blur on the image processor.
   *
   * @param ip the image
   * @param sigma The Gaussian width
   */
  @Override
  public void blurGaussian(ImageProcessor ip, double sigma) {
    final double accuracy =
        (ip instanceof ByteProcessor || ip instanceof ColorProcessor) ? 0.002 : 0.0002;
    blurGaussian(ip, sigma, sigma, accuracy);
  }

  private static ImageProcessor duplicateProcessor(ImageProcessor ip) {
    final ImageProcessor duplicateIp = ip.duplicate();
    if (ip.getRoi().height != ip.getHeight() || ip.getRoi().width != ip.getWidth()) {
      duplicateIp.snapshot();
      duplicateIp.setRoi(ip.getRoi());
      duplicateIp.setMask(ip.getMask());
    }
    return duplicateIp;
  }
}
