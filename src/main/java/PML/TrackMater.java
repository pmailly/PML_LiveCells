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
import fiji.plugin.trackmate.action.ExportTracksToXML;
import fiji.plugin.trackmate.detection.DogDetectorFactory;
import fiji.plugin.trackmate.detection.LogDetectorFactory;
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
import java.util.Iterator;
import java.util.List;
import java.util.Set;




public class TrackMater extends TrackMatePlugIn_ {
    
    //private PML_Tools pml = new PML_Tools();
   
    private Logger logger = new LogRecorder( Logger.DEFAULT_LOGGER );
    private String welcomeMessage;
    
    
    public void run(ImagePlus imp, String path, String imgname, double radius, double threshold, String detect, double link, double merging_dist, boolean subpixel, boolean median)
    {
       //Initialisation TrackMate.
        logger = new LogRecorder( Logger.VOID_LOGGER );
        settings = createSettings( imp );
        settings.imageFileName = imgname;
        settings.imageFolder = path;
	model = createModel();
	model.setLogger( logger );
	trackmate = createTrackMate();
    
    
       // Configure default settings.
        // Default detector.
        if (detect.equals("LoG"))
            settings.detectorFactory = new LogDetectorFactory();
        else
            settings.detectorFactory = new DogDetectorFactory();
        
        settings.detectorSettings = settings.detectorFactory.getDefaultSettings();

        // Default tracker.
        settings.trackerFactory = new SimpleSparseLAPTrackerFactory();
        settings.trackerSettings = settings.trackerFactory.getDefaultSettings();
        
        // Set-up detector
        settings.detectorSettings.put( "RADIUS", radius );
        settings.detectorSettings.put( "THRESHOLD", threshold );
        settings.detectorSettings.put( "DO_SUBPIXEL_LOCALIZATION", subpixel );
        settings.detectorSettings.put( "DO_MEDIAN_FILTERING", median );
        settings.detectorSettings.put( "TARGET_CHANNEL", 0 );
        
        // Set-up tracker
        settings.trackerSettings.put("LINKING_MAX_DISTANCE", link);
        settings.trackerSettings.put("MAX_FRAME_GAP", 0);
        settings.trackerSettings.put("GAP_CLOSING_MAX_DISTANCE", 1.0);
        
        settings.trackerSettings.put("ALLOW_TRACK_MERGING", true);
        settings.trackerSettings.put("MERGING_MAX_DISTANCE", merging_dist);
        settings.trackerSettings.put("ALLOW_TRACK_SPLITTING", true);
        settings.trackerSettings.put("SPLITTING_MAX_DISTANCE", merging_dist);
        // Run trackMate with the settings
        welcomeMessage = TrackMate.PLUGIN_NAME_STR + " v" + TrackMate.PLUGIN_NAME_VERSION + " started on:\n" + TMUtils.getCurrentTimeString() + '\n';
        if ( !trackmate.checkInput() || !trackmate.process() )
        {
                IJ.error( "Error while performing tracking:\n" + trackmate.getErrorMessage() );
                return;
        }
    }
    
    public void saveResults(String savefile, String exportfile, String statfile, String trackfile, String path, String imgname) {
        // Save resultats
        final String save_path_str = savefile;
        final File save_path = new File( save_path_str );
        final TmXmlWriter writer = new TmXmlWriter( save_path, logger );

        writer.appendLog(welcomeMessage);
        writer.appendModel( trackmate.getModel() );
        writer.appendSettings( trackmate.getSettings() );
        try
        {
                writer.writeToFile();
        }
        catch ( final Exception e )
        {
                IJ.error( "When saving to " + save_path + ", file not found:\n" + e.getMessage() + '\n' );
                return;
        }
        
        // Export results to XML
//        final String export_path_str = exportfile;
//        final File export_path = new File( export_path_str );
//        try
//        {
//                ExportTracksToXML.export( model, settings, export_path );
//        }
//        catch ( final Exception e )
//        {
//                IJ.error( "When exporting to " + export_path + ", file not found:\n" + e.getMessage() + '\n' );
//                return;
//        }
      
        
         // Export statistics file ?        
        final SelectionModel selectionModel = new SelectionModel( model );
        ExportStats stat = new ExportStats(selectionModel);
        stat.execute(trackmate);
        ResultsTable statTable = stat.getSpotTable();
        statTable.save(statfile);
        
        // Export fusion/splitting stats
        ResultsTable trackTable = stat.getTrackTable();
        trackTable.save(trackfile);
        
        
        IJ.showStatus("Tracking done "+imgname);

    }
    
    public void getTracks(){
        Set<Integer> trackIDs = model.getTrackModel().trackIDs(true); 
        final Iterator<Integer> it = trackIDs.iterator();
	while (it.hasNext()) {
            final Integer id = it.next();
            System.out.println(id);
	}

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
