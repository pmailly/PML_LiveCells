/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package PML;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import org.jgrapht.graph.DefaultWeightedEdge;
import fiji.plugin.trackmate.FeatureModel;
import fiji.plugin.trackmate.Spot;
import fiji.plugin.trackmate.action.AbstractTMAction;
import fiji.plugin.trackmate.action.TrackMateAction;
import fiji.plugin.trackmate.action.TrackMateActionFactory;
import fiji.plugin.trackmate.features.edges.EdgeTimeLocationAnalyzer;
import fiji.plugin.trackmate.util.ModelTools;
import ij.measure.ResultsTable;

import javax.swing.ImageIcon;


import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.gui.displaysettings.DisplaySettings;
import java.awt.Frame;


public class ExportStats extends AbstractTMAction
{

	public static final ImageIcon ICON = new ImageIcon();

	public static final String NAME = "Export statistics to tables";

	public static final String KEY = "EXPORT_STATS_TO_IJ";

	public static final String INFO_TEXT = "<html>"
			+ "Compute and export all statistics to 3 ImageJ results table. "
			+ "Statistisc are separated in features computed for: "
			+ "<ol> "
			+ "	<li> spots in filtered tracks; "
			+ "	<li> links between those spots; "
			+ "	<li> filtered tracks. "
			+ "</ol> "
			+ "For tracks and links, they are recalculated prior to exporting. Note "
			+ "that spots and links that are not in a filtered tracks are not part "
			+ "of this export."
			+ "</html>";

	private static final String SPOT_TABLE_NAME = "Spots in tracks statistics";

	private static final String EDGE_TABLE_NAME = "Links in tracks statistics";

	private static final String TRACK_TABLE_NAME = "Track statistics";

	private static final String ID_COLUMN = "ID";

	private static final String TRACK_ID_COLUMN = "TRACK_ID";

	private ResultsTable spotTable;

	private ResultsTable edgeTable;

	private ResultsTable trackTable;

	private final SelectionModel selectionModel;

        public ExportStats(){this.selectionModel = null;}

	
        @Override
	public void execute( final TrackMate trackmate, final SelectionModel selectionModel, final DisplaySettings displaySettings, final Frame parent  )
	{
		logger.log( "Exporting statistics.\n" );

		// Model
		final Model model = trackmate.getModel();
		final FeatureModel fm = model.getFeatureModel();

		// Export spots
		logger.log( "  - Exporting spot statistics..." );
		final Set< Integer > trackIDs = model.getTrackModel().trackIDs( true );
		final Collection< String > spotFeatures = trackmate.getModel().getFeatureModel().getSpotFeatures();

		this.spotTable = new ResultsTable();

		// Parse spots to insert values as objects
		for ( final Integer trackID : trackIDs )
		{
			final Set< Spot > track = model.getTrackModel().trackSpots( trackID );
			// Sort by frame
			final List< Spot > sortedTrack = new ArrayList<>( track );
			Collections.sort( sortedTrack, Spot.frameComparator );

			for ( final Spot spot : sortedTrack )
			{
				spotTable.incrementCounter();
				spotTable.addLabel( spot.getName() );
				spotTable.addValue( ID_COLUMN, "" + spot.ID() );
				spotTable.addValue( "TRACK_ID", "" + trackID.intValue() );
				for ( final String feature : spotFeatures )
				{
					final Double val = spot.getFeature( feature );
					if ( null == val )
					{
						spotTable.addValue( feature, "None" );
					}
					else
					{
						if ( fm.getSpotFeatureIsInt().get( feature ).booleanValue() )
						{
							spotTable.addValue( feature, "" + val.intValue() );
						}
						else
						{
							spotTable.addValue( feature, val.doubleValue() );
						}
					}
				}
			}
		}
		logger.log( " Done.\n" );

		// Export edges
		logger.log( "  - Exporting links statistics..." );
		// Yield available edge feature
		final Collection< String > edgeFeatures = fm.getEdgeFeatures();

		this.edgeTable = new ResultsTable();

		// Sort by track
		for ( final Integer trackID : trackIDs )
		{
			// Comparators
			final Comparator< DefaultWeightedEdge > edgeTimeComparator = ModelTools.featureEdgeComparator( EdgeTimeLocationAnalyzer.TIME, fm );
			final Comparator< DefaultWeightedEdge > edgeSourceSpotTimeComparator = new EdgeSourceSpotFrameComparator( model );

			final Set< DefaultWeightedEdge > track = model.getTrackModel().trackEdges( trackID );
			final List< DefaultWeightedEdge > sortedTrack = new ArrayList<>( track );

			/*
			 * Sort them by frame, if the EdgeTimeLocationAnalyzer feature is
			 * declared.
			 */

			if ( model.getFeatureModel().getEdgeFeatures().contains( EdgeTimeLocationAnalyzer.KEY ) )
				Collections.sort( sortedTrack, edgeTimeComparator );
			else
				Collections.sort( sortedTrack, edgeSourceSpotTimeComparator );

			for ( final DefaultWeightedEdge edge : sortedTrack )
			{
				edgeTable.incrementCounter();
				edgeTable.addLabel( edge.toString() );
				edgeTable.addValue( TRACK_ID_COLUMN, "" + trackID.intValue() );
				for ( final String feature : edgeFeatures )
				{
					final Object o = fm.getEdgeFeature( edge, feature );
					if ( o instanceof String )
					{
						continue;
					}
					final Number d = ( Number ) o;
					if ( d == null )
					{
						edgeTable.addValue( feature, "None" );
					}
					else
					{
						if ( fm.getEdgeFeatureIsInt().get( feature ).booleanValue() )
						{
							edgeTable.addValue( feature, "" + d.intValue() );
						}
						else
						{
							edgeTable.addValue( feature, d.doubleValue() );
						}

					}
				}

			}
		}
		logger.log( " Done.\n" );

		// Export tracks
		logger.log( "  - Exporting tracks statistics..." );
		// Yield available edge feature
		final Collection< String > trackFeatures = fm.getTrackFeatures();

		this.trackTable = new ResultsTable();

		// Sort by track
		for ( final Integer trackID : trackIDs )
		{
			trackTable.incrementCounter();
			trackTable.addLabel( model.getTrackModel().name( trackID ) );
			trackTable.addValue( TRACK_ID_COLUMN, "" + trackID.intValue() );
			for ( final String feature : trackFeatures )
			{
				final Double val = fm.getTrackFeature( trackID, feature );
				if ( null == val )
				{
					trackTable.addValue( feature, "None" );
				}
				else
				{
					if ( fm.getTrackFeatureIsInt().get( feature ).booleanValue() )
					{
						trackTable.addValue( feature, "" + val.intValue() );
					}
					else
					{
						trackTable.addValue( feature, val.doubleValue() );
					}
				}
			}
		}
		logger.log( " Done.\n" );

		// Show tables
		//spotTable.show( SPOT_TABLE_NAME );
		//edgeTable.show( EDGE_TABLE_NAME );
		//trackTable.show( TRACK_TABLE_NAME );

		
	}

	/**
	 * Returns the results table containing the spot statistics, or
	 * <code>null</code> if the {@link #execute(TrackMate)} method has not been
	 * called.
	 *
	 * @return the results table containing the spot statistics.
	 */
	public ResultsTable getSpotTable()
	{
		return spotTable;
	}

	/**
	 * Returns the results table containing the edge statistics, or
	 * <code>null</code> if the {@link #execute(TrackMate)} method has not been
	 * called.
	 *
	 * @return the results table containing the edge statistics.
	 */
	public ResultsTable getEdgeTable()
	{
		return edgeTable;
	}

	/**
	 * Returns the results table containing the track statistics, or
	 * <code>null</code> if the {@link #execute(TrackMate)} method has not been
	 * called.
	 *
	 * @return the results table containing the track statistics.
	 */
	public ResultsTable getTrackTable()
	{
		return trackTable;
	}

	// Invisible because called on the view config panel.
	//@Plugin( type = TrackMateActionFactory.class, visible = false )
	public static class Factory implements TrackMateActionFactory
	{

		@Override
		public String getInfoText()
		{
			return INFO_TEXT;
		}

		@Override
		public String getKey()
		{
			return KEY;
		}

		@Override
		public TrackMateAction create()
		{
			return new ExportStats();
		}

		@Override
		public ImageIcon getIcon()
		{
			return ICON;
		}

		@Override
		public String getName()
		{
			return NAME;
		}
	}

	private static final class EdgeSourceSpotFrameComparator implements Comparator< DefaultWeightedEdge >
	{

		private final Model model;

		public EdgeSourceSpotFrameComparator( final Model model )
		{
			this.model = model;
		}

		@Override
		public int compare( final DefaultWeightedEdge e1, final DefaultWeightedEdge e2 )
		{
			final double t1 = model.getTrackModel().getEdgeSource( e1 ).getFeature( Spot.FRAME ).doubleValue();
			final double t2 = model.getTrackModel().getEdgeSource( e2 ).getFeature( Spot.FRAME ).doubleValue();
			if ( t1 < t2 ) { return -1; }
			if ( t1 > t2 ) { return 1; }
			return 0;
		}

	}

}