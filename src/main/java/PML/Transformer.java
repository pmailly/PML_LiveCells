package PML;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.process.ImageConverter;
import ij.process.ImageStatistics;
import ij.process.ShortProcessor;
import java.lang.reflect.Method;

/**
 * Transformation matrix from StackRegPlus (globalTransform)
 * Anchor points from StackRegPlus
 * Apply TurboReg transformation
 * 
 */
public class Transformer {
    protected double[][] sourcePts;

    private int width;
    private int height;
    
    /** \brief initialise */
    public Transformer()
    {
        sourcePts = new double[3][2];
        for ( int i=0; i < 3; i ++)
        {
            for ( int j=0; j<2; j++ )
                sourcePts[i][j] = 0;
        }
    }
    
    public void setSize(int w, int h)
    {
        width = w;
        height = h;
    }
    

        /** \brief Fill matrix  of source points */
    public void setSource(double[][] mat)
    {
        for ( int i=0; i < 3; i ++)
        {
            for ( int j=0; j<2; j++ )
                sourcePts[i][j] = mat[i][j];
        }
    }
    
    /** \brief Apply the saved transformation with TurboReg 
     * Modify directly the input imageplus
     */
    public void doTransformation(ImagePlus imp)
    {
        try 
        {
            Object turboReg = null;
            IJ.run(imp, "16-bit","");
            //ImageStatistics stat = imp.getStatistics();
            for ( int z=1; z<=imp.getNSlices(); z++)
            {
                imp.setSlice(z);
                ImagePlus source = new ImagePlus("StackRegSource", new ShortProcessor(
					width, height, (short[])imp.getProcessor().getPixels(),
					imp.getProcessor().getColorModel()));
                final FileSaver sourceFile = new FileSaver(source);
                final String sourcePathAndFileName = IJ.getDirectory("temp") + source.getTitle();
                sourceFile.saveAsTiff(sourcePathAndFileName);
                turboReg = IJ.runPlugIn("TurboReg_", "-transform"
							+ " -file " + sourcePathAndFileName
							+ " " + width
							+ " " + height
							+ " -rigidBody"
							+ " " + sourcePts[0][0]
							+ " " + sourcePts[0][1]
							+ " " + (width / 2)
							+ " " + (height / 2)
							+ " " + sourcePts[1][0]
							+ " " + sourcePts[1][1]
							+ " " + (width / 2)
							+ " " + (height / 4)
							+ " " + sourcePts[2][0]
							+ " " + sourcePts[2][1]
							+ " " + (width / 2)
							+ " " + ((3 * height) / 4)
							+ " -hideOutput"
						);
                if (turboReg == null) throw(new ClassNotFoundException());
	
                Method method = turboReg.getClass().getMethod("getTransformedImage",(Class[])null);
                ImagePlus transformedSource =(ImagePlus) method.invoke(turboReg);
                transformedSource.getStack().deleteLastSlice();
                transformedSource.getProcessor().setMinAndMax(0.0, 65535.0);
                final ImageConverter converter = new ImageConverter(transformedSource);
                converter.convertToGray16();
                imp.setProcessor(null, transformedSource.getProcessor());
               }        
             //double newval = stat.mean - 1*stat.stdDev;    
             //IJ.run(imp, "Macro...", "code=v=v+(v==0)*"+newval+" stack");
        }
        catch (Exception e)
        {
            IJ.error("Error in alignement with TurboReg "+e.toString());
            return;
        }
    }
    
    
    
    
}