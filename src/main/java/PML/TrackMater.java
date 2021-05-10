/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package PML;

import fiji.plugin.trackmate.Logger;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.Settings;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.TrackMatePlugIn_;
import fiji.plugin.trackmate.action.ExportAllSpotsStatsAction;
import fiji.plugin.trackmate.action.ExportStatsToIJAction;
import fiji.plugin.trackmate.action.ExportTracksToXML;
import fiji.plugin.trackmate.detection.DogDetectorFactory;
import fiji.plugin.trackmate.features.edges.EdgeAnalyzer;
import fiji.plugin.trackmate.features.spot.SpotAnalyzerFactory;
import fiji.plugin.trackmate.features.track.TrackAnalyzer;
import fiji.plugin.trackmate.io.TmXmlWriter;
import fiji.plugin.trackmate.providers.EdgeAnalyzerProvider;
import fiji.plugin.trackmate.providers.SpotAnalyzerProvider;
import fiji.plugin.trackmate.providers.TrackAnalyzerProvider;
import fiji.plugin.trackmate.tracking.sparselap.SimpleSparseLAPTrackerFactory;
import fiji.plugin.trackmate.util.LogRecorder;
import fiji.plugin.trackmate.util.TMUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.measure.ResultsTable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import net.imglib2.util.ValuePair;


public class TrackMater extends TrackMatePlugIn_ {
 
    private double radius;
    private double threshold;
    private boolean subpixel = true;
    private boolean median = true;
    private Logger logger = new LogRecorder( Logger.DEFAULT_LOGGER );
    
    public void setDetectorParameters(double rad, double thres)
    {
        radius = rad;
        threshold = thres;
    }
    
    public void run(ImagePlus imp, String savefile, String exportfile, String statfile)
    {
       //Initialisation TrackMate.
        logger = new LogRecorder( Logger.IJ_LOGGER );
        settings = createSettings( imp );
	model = createModel();
	model.setLogger( logger );
	trackmate = createTrackMate();
    
    
       // Configure default settings.
        // Default detector.
        settings.detectorFactory = new DogDetectorFactory< >();
        settings.detectorSettings = settings.detectorFactory.getDefaultSettings();

        // Default tracker.
        settings.trackerFactory = new SimpleSparseLAPTrackerFactory();
        settings.trackerSettings = settings.trackerFactory.getDefaultSettings();
        
        // Set-up detector
        settings.detectorSettings.put( "radius", radius );
        settings.detectorSettings.put( "threshold", threshold );
        settings.detectorSettings.put( "subpixel", subpixel );
        settings.detectorSettings.put( "median", median );
        settings.detectorSettings.put( "channel", 1 );
        
        // Set-up tracker
        //settings.trackerSettings.put("max_distance", 1);
        settings.trackerSettings.replace("max_frame_gap", 0);
        settings.trackerSettings.replace("max_gap_distance", 1);
     
        trackmate.getSettings().trackerSettings.replace("max_distance", 1);
        Set keys = settings.trackerSettings.keySet();
        System.out.println(keys.toString());
        
        // Run trackMate with the settings
        final String welcomeMessage = TrackMate.PLUGIN_NAME_STR + " v" + TrackMate.PLUGIN_NAME_VERSION + " started on:\n" + TMUtils.getCurrentTimeString() + '\n';
	logger.log( welcomeMessage );
        if ( !trackmate.checkInput() || !trackmate.process() )
        {
                logger.error( "Error while performing tracking:\n" + trackmate.getErrorMessage() );
                return;
        }
        
        // Save resultats
        final String save_path_str = savefile;
        final File save_path = new File( save_path_str );
        final TmXmlWriter writer = new TmXmlWriter( save_path, logger );

        writer.appendLog( logger.toString() );
        writer.appendModel( trackmate.getModel() );
        writer.appendSettings( trackmate.getSettings() );
        try
        {
                writer.writeToFile();
                logger.log( "Data saved to: " + save_path.toString() + '\n' );
        }
        catch ( final FileNotFoundException e )
        {
                logger.error( "When saving to " + save_path + ", file not found:\n" + e.getMessage() + '\n' );
                return;
        }
        catch ( final IOException e )
        {
                logger.error( "When saving to " + save_path + ", Input/Output error:\n" + e.getMessage() + '\n' );
                return;
        }
        
        // Export results to XML
        final String export_path_str = exportfile;
        final File export_path = new File( export_path_str );
        try
        {
                ExportTracksToXML.export( model, settings, export_path );
                logger.log( "Data exported to: " + export_path.toString() + '\n' );
        }
        catch ( final FileNotFoundException e )
        {
                logger.error( "When exporting to " + export_path + ", file not found:\n" + e.getMessage() + '\n' );
                return;
        }
        catch ( final IOException e )
        {
                logger.error( "When exporting to " + export_path + ", Input/Output error:\n" + e.getMessage() + '\n' );
                return;
        }
        
         // Export statistics file ?        
        final SelectionModel selectionModel = new SelectionModel( model );
        ExportStatsToIJAction stat = new ExportStatsToIJAction(selectionModel);
        stat.execute(trackmate);
        ResultsTable statTable = stat.getSpotTable();
        statTable.save(statfile);
        logger.log( "Stat data saved in: " + statfile + '\n' );
        

    }
    
	@Override
	protected Settings createSettings( final ImagePlus imp )
	{
		final Settings s = super.createSettings( imp );

		s.clearSpotAnalyzerFactories();
		final SpotAnalyzerProvider spotAnalyzerProvider = new SpotAnalyzerProvider();
		final List< String > spotAnalyzerKeys = spotAnalyzerProvider.getKeys();
		for ( final String key : spotAnalyzerKeys )
		{
			final SpotAnalyzerFactory< ? > spotFeatureAnalyzer = spotAnalyzerProvider.getFactory( key );
			s.addSpotAnalyzerFactory( spotFeatureAnalyzer );
		}

		s.clearEdgeAnalyzers();
		final EdgeAnalyzerProvider edgeAnalyzerProvider = new EdgeAnalyzerProvider();
		final List< String > edgeAnalyzerKeys = edgeAnalyzerProvider.getKeys();
		for ( final String key : edgeAnalyzerKeys )
		{
			final EdgeAnalyzer edgeAnalyzer = edgeAnalyzerProvider.getFactory( key );
			s.addEdgeAnalyzer( edgeAnalyzer );
		}

		s.clearTrackAnalyzers();
		final TrackAnalyzerProvider trackAnalyzerProvider = new TrackAnalyzerProvider();
		final List< String > trackAnalyzerKeys = trackAnalyzerProvider.getKeys();
		for ( final String key : trackAnalyzerKeys )
		{
			final TrackAnalyzer trackAnalyzer = trackAnalyzerProvider.getFactory( key );
			s.addTrackAnalyzer( trackAnalyzer );
		}

		return s;
	}
    
}
