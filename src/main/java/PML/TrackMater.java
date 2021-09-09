/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package PML;

import fiji.plugin.trackmate.Logger;
import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.Settings;
import fiji.plugin.trackmate.Spot;
import fiji.plugin.trackmate.SpotCollection;
import fiji.plugin.trackmate.TrackMatePlugIn_;
import fiji.plugin.trackmate.detection.DogDetectorFactory;
import fiji.plugin.trackmate.detection.LogDetectorFactory;
import fiji.plugin.trackmate.detection.ManualDetectorFactory;
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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Object3D;
import mcib3d.geom.Vector3D;



public class TrackMater extends TrackMatePlugIn_ {
    
    private Logger logger = new LogRecorder( Logger.DEFAULT_LOGGER );
    private String welcomeMessage;
    private int nSpots = 0;
    
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
    
    public boolean trackmateObjects(ImagePlus imp, ArrayList<Objects3DPopulation> pmlPopList, String path, String imgname, double radius, double link, double merging_dist)
    {
       //Initialisation TrackMate.
        logger = new LogRecorder( Logger.VOID_LOGGER );
        //logger = new LogRecorder( Logger.IJ_LOGGER );
        
        final String spaceUnit = imp.getCalibration().getUnit();
        final String timeUnit = imp.getCalibration().getTimeUnit();
	final double frameInterval = imp.getCalibration().frameInterval;

	getModel(pmlPopList, frameInterval, spaceUnit, timeUnit );
        if (nSpots ==0) return false; // no spots found
        settings = createSettings( imp );
        settings.imageFileName = imgname;
        settings.imageFolder = path;
	//fillSettings();
        // Default tracker.
        settings.detectorFactory = new ManualDetectorFactory<>();
	settings.detectorSettings = settings.detectorFactory.getDefaultSettings();

        settings.trackerFactory = new SimpleSparseLAPTrackerFactory();
        settings.trackerSettings = settings.trackerFactory.getDefaultSettings();
        
        // Set-up tracker
        settings.trackerSettings.put("LINKING_MAX_DISTANCE", link);
        settings.trackerSettings.put("MAX_FRAME_GAP", 0);
        settings.trackerSettings.put("GAP_CLOSING_MAX_DISTANCE", 1.0);
        
        settings.trackerSettings.put("ALLOW_TRACK_MERGING", true);
        settings.trackerSettings.put("MERGING_MAX_DISTANCE", merging_dist);
        settings.trackerSettings.put("ALLOW_TRACK_SPLITTING", true);
        settings.trackerSettings.put("SPLITTING_MAX_DISTANCE", merging_dist);

	logger.log( "Computing features.\n" );
	 // Run trackMate with the settings
        trackmate = createTrackMate();
        if ( !runTracking()  ) { 
            IJ.error( "Error while performing tracking:\n" + trackmate.getErrorMessage() );
            return false;
        }
        return true;
    }
    
    public boolean runTracking(){
        if ( !trackmate.computeSpotFeatures( true ) ) { return false; }
        if ( !trackmate.execSpotFiltering( true ) ) { return false; }
        if ( !trackmate.execTracking() ) { return false; }
        if ( !trackmate.computeTrackFeatures( true ) ) { return false; }
        if ( !trackmate.execTrackFiltering( true ) ) { return false; }
        if ( !trackmate.computeEdgeFeatures( true ) ) { return false; }
        return true;
    }
    
  /** public void fillSettings()
    {
	//settings.detectorFactory = new ManualDetectorFactory<>();
	//settings.detectorSettings = settings.detectorFactory.getDefaultSettings();

        // Spot features.
        settings.addSpotAnalyzerFactory( new ManualSpotColorAnalyzerFactory<>() );

        // Edge features.
        settings.addEdgeAnalyzer( (EdgeAnalyzer) new EdgeTargetAnalyzer() );
        settings.addEdgeAnalyzer( new EdgeTimeLocationAnalyzer() );
        settings.addEdgeAnalyzer( new EdgeVelocityAnalyzer() );
        settings.addEdgeAnalyzer( (EdgeAnalyzer) new ManualEdgeColorAnalyzer() );

        // Track features.
        settings.addTrackAnalyzer( (TrackAnalyzer) new TrackDurationAnalyzer() );
        settings.addTrackAnalyzer((TrackAnalyzer)  new TrackIndexAnalyzer() );
        settings.addTrackAnalyzer( (TrackAnalyzer) new TrackLocationAnalyzer() );
        settings.addTrackAnalyzer( (TrackAnalyzer) new TrackSpeedStatisticsAnalyzer() );
        settings.addTrackAnalyzer( (TrackAnalyzer) new TrackSpotQualityFeatureAnalyzer() );

        logger.log( "Added the following features to be computed:\n" + settings.toStringFeatureAnalyzersInfo() );
    }*/
    
    public void getModel(ArrayList<Objects3DPopulation> pmlPopList, final double frameInterval, final String spaceUnit, final String timeUnit )
    {
		Map< Integer, Set< Spot > > spots = new HashMap<>();
		Map< Integer, List< Spot > > tracks = new HashMap<>();
                
                // Iterate over objects.
		logger.log( String.format( "Parsing records.\n" ) );
		nSpots = 0;
                double q = 1.;
                int t = 0;
                int sizepop = pmlPopList.size();
                // Read each time point population
                while (!pmlPopList.isEmpty())
		//for ( Objects3DPopulation pop: pmlPopList )
		{
                    Objects3DPopulation pop = pmlPopList.get(0);
                    pmlPopList.remove(0);
                    logger.setProgress( ( double ) nSpots / sizepop );

                    // list of objects at time t
                    Set< Spot > list = new HashSet<>();
                    spots.put( Integer.valueOf( t ), list );

                    for ( int i=0; i < pop.getNbObjects(); i++ ) {
                        nSpots++;
                        Object3D obj = pop.getObject(i);
                        double r = Math.pow(4.0/3.0/Math.PI*obj.getVolumeUnit(), 0.33);	
                        Vector3D cent = obj.getCenterUnit();
                        Spot spot = new Spot( cent.x, cent.y, cent.z, r, q, null );
                        spot.putFeature( Spot.FRAME, ( double ) t );
                        spot.putFeature( Spot.POSITION_T, frameInterval * t );
                        list.add( spot );
                    }
                    t++;
                    pop = null;
                }	
		
		logger.log( String.format( "Parsing done. Loaded %d objects.\n", nSpots ) );
		if ( nSpots== 0 ) return;
               
		//Generate a Model object. 
		SpotCollection sc = SpotCollection.fromMap( spots );
		sc.setVisible( true );
		logger.log( String.format( "Found %d spots.\n", sc.getNSpots( true ) ) );

		NavigableSet< Integer > frames = sc.keySet();
		for ( final Integer frame : frames )
			logger.log( String.format( "- frame %4d, n spots = %d\n", frame, sc.getNSpots( frame, true ) ) );

		model = new Model();
		model.setPhysicalUnits( spaceUnit, timeUnit );
		model.setLogger( logger );
		model.setSpots( sc, false );
		logger.setProgress( 0. );
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
