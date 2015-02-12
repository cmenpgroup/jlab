package infn.bed.bedview;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.List;

import cnuphys.bCNU.attributes.AttributeType;
import cnuphys.bCNU.drawable.DrawableAdapter;
import cnuphys.bCNU.drawable.IDrawable;
import cnuphys.bCNU.graphics.GraphicsUtilities;
import cnuphys.bCNU.graphics.container.IContainer;
import cnuphys.bCNU.graphics.style.Styled;
import cnuphys.bCNU.graphics.toolbar.BaseToolBar;
import cnuphys.bCNU.graphics.world.WorldGraphicsUtilities;
import cnuphys.bCNU.layer.LogicalLayer;
import cnuphys.bCNU.util.X11Colors;
import infn.bed.geometry.GeoConstants;
import infn.bed.component.ControlPanel;
import infn.bed.item.FullSideViewBar;
import infn.bed.item.FullSideViewVeto;

/**
 * This class handles the drawing of the full side view of the bars and the
 * vetoes. It orders the bars in a 3x3 grid and numbers them going from
 * left-to-right, top-to-bottom. It orders the vetoes going clockwise starting
 * from the upstream internal veto and fully around the 3x3 grid then to the
 * left veto followed by the right veto, then repeated for the external vetoes:
 * 
 * 1. Internal Upstream 
 * 2. Internal Top 
 * 3. Internal Downstream 
 * 4. Internal Bottom 
 * 5. Internal Left 
 * 6. Internal Right 
 * 7. External Upstream 
 * 8. External Top-Upstream 
 * 9. External Top-Downstream 
 * 10. External Downstream 
 * 11. External Bottom-Downstream 
 * 12. External Bottom-Upstream 
 * 13. External Left 
 * 14. External Right
 * 
 * This view is used to highlight hits in the vetoes and the bars as well as
 * display information for energy and time of hits.
 * 
 * @author Andy Beiter
 * 
 */
@SuppressWarnings("serial")
public class FullSideView extends BedView {

	/**
	 * A bar rectangle for each bar.
	 */
	private Rectangle2D.Double _barWorldRects[];

	/**
	 * A veto rectangle for each veto.
	 */
	private Rectangle2D.Double _vetoWorldRects[];

	/**
	 * Used for drawing and customizing the bar (and veto) rectangles.
	 */
	private Styled _barStyle;

	/**
	 * Used for the before draw for rectangles (not very important).
	 */
	private IDrawable _beforeDraw;

	/**
	 * Array of bar instances. These handle displaying hits and info.
	 */
	private FullSideViewBar _superLayerBars[];

	/**
	 * Array of veto instances. These handle displaying hits and info.
	 */
	private FullSideViewVeto _superLayerVetoes[];

	/**
	 * The grid used for positioning the rectangles. The shape is a 3x3 grid.
	 */
	private static Rectangle2D.Double _defaultWorldRectangle = new Rectangle2D.Double(
			0.0, 0.0, 3.0, 3.0);

	/**
	 * Constructor used to set up parameters and to handle the drawing and
	 * adding of bar/veto rectangles and objects.
	 * 
	 * @param keyVals
	 *            used in the super class (BedView) to set up parameters of this
	 *            view
	 */
	private FullSideView(Object... keyVals) {
		super(keyVals);

		setbarWorldRects();
		setBeforeDraw();
		setAfterDraw();
		addItems();
	}

	/**
	 * Method used to create an instance of this class
	 * 
	 * @return the new instance of the view
	 */
	public static FullSideView createFullSideView() {
		FullSideView view = null;

		// set to a fraction of screen
		Dimension d = GraphicsUtilities.screenFraction(0.5);

		// create the view
		view = new FullSideView(
				AttributeType.WORLDSYSTEM,
				_defaultWorldRectangle,
				AttributeType.WIDTH,
				d.width, // container width, not total view width
				AttributeType.HEIGHT,
				d.height, // container height, not total view width
				AttributeType.TOOLBAR, true, AttributeType.TOOLBARBITS,
				BaseToolBar.NODRAWING & ~BaseToolBar.RANGEBUTTON
						& ~BaseToolBar.TEXTFIELD
						& ~BaseToolBar.CONTROLPANELBUTTON
						& ~BaseToolBar.TEXTBUTTON & ~BaseToolBar.DELETEBUTTON,
				AttributeType.VISIBLE, true, AttributeType.HEADSUP, false,
				AttributeType.TITLE, "Full Side View",
				AttributeType.STANDARDVIEWDECORATIONS, true);

		view._controlPanel = new ControlPanel(view, ControlPanel.FEEDBACK, 0);

		view.add(view._controlPanel, BorderLayout.EAST);
		view.pack();
		return view;
	}

	/**
	 * This method creates the shapes for each of the bar rectangles and veto
	 * rectangles.
	 */
	private void setbarWorldRects() {

		_barWorldRects = new Rectangle2D.Double[9];
		_vetoWorldRects = new Rectangle2D.Double[14];

		Rectangle2D.Double worldRect = _defaultWorldRectangle;

		double gap = worldRect.width / 48;
		double boxWidth = worldRect.width / 12 - 2 * gap;
		double boxHeight = worldRect.height / 12 - 2 * gap;

		double left = worldRect.getMinX();
		double right = worldRect.getMaxX();

		double barLeft = 1.5 - (1.5 * boxWidth);
		double barBot = 1.5 - boxWidth;

		double x13 = barLeft + (boxWidth);
		double x23 = x13 + (boxWidth);
		double y13 = barBot + (boxWidth);
		double y23 = y13 + (boxWidth);

		_barWorldRects[0] = new Rectangle2D.Double(barLeft, y23, boxWidth,
				boxHeight);
		_barWorldRects[1] = new Rectangle2D.Double(x13, y23, boxWidth,
				boxHeight);
		_barWorldRects[2] = new Rectangle2D.Double(x23, y23, boxWidth,
				boxHeight);

		_barWorldRects[3] = new Rectangle2D.Double(barLeft, y13, boxWidth,
				boxHeight);
		_barWorldRects[4] = new Rectangle2D.Double(x13, y13, boxWidth,
				boxHeight);
		_barWorldRects[5] = new Rectangle2D.Double(x23, y13, boxWidth,
				boxHeight);

		_barWorldRects[6] = new Rectangle2D.Double(barLeft, barBot, boxWidth,
				boxHeight);
		_barWorldRects[7] = new Rectangle2D.Double(x13, barBot, boxWidth,
				boxHeight);
		_barWorldRects[8] = new Rectangle2D.Double(x23, barBot, boxWidth,
				boxHeight);

		// ------ VETOES ------
		// organized by location for ease of debugging
		// for numbering, see class documentation

		// upstream
		_vetoWorldRects[0] = new Rectangle2D.Double(barLeft - 2 * gap, barBot,
				gap, 3 * boxHeight);// internal
		_vetoWorldRects[6] = new Rectangle2D.Double(barLeft - 4 * gap, barBot
				- 2 * gap, gap, 5 * boxHeight);// external

		// downstream
		_vetoWorldRects[2] = new Rectangle2D.Double(barLeft + 3 * boxWidth
				+ gap, barBot, gap, 3 * boxHeight);// internal
		_vetoWorldRects[9] = new Rectangle2D.Double(barLeft + 3 * boxWidth + 3
				* gap, barBot - 2 * gap, gap, 5 * boxHeight);// external

		// left
		_vetoWorldRects[4] = new Rectangle2D.Double(left + (3 * boxWidth) + 3
				* gap, barBot - 0.5 * gap, 3 * boxWidth + gap, 3 * boxHeight
				+ gap); // internal
		_vetoWorldRects[12] = new Rectangle2D.Double(left + gap, barBot - 0.5
				* gap, 3 * boxWidth + gap, 3 * boxHeight + gap); // external

		// right
		_vetoWorldRects[5] = new Rectangle2D.Double(right - 6 * boxWidth - 4
				* gap, barBot - 0.5 * gap, 3 * boxWidth + gap, 3 * boxHeight
				+ gap); // internal
		_vetoWorldRects[13] = new Rectangle2D.Double(right - 3 * boxWidth - 2
				* gap, barBot - 0.5 * gap, 3 * boxWidth + gap, 3 * boxHeight
				+ gap); // external

		// bot
		_vetoWorldRects[3] = new Rectangle2D.Double(barLeft, barBot - 2 * gap,
				3 * boxWidth, gap);// internal
		_vetoWorldRects[11] = new Rectangle2D.Double(barLeft - 2 * gap, barBot
				- 4 * gap, 2 * boxWidth + gap, gap);// external - upstream
		_vetoWorldRects[10] = new Rectangle2D.Double(barLeft + 1.5 * boxWidth,
				barBot - 4 * gap, 2 * boxWidth + gap, gap);// external -
															// downstream

		// top
		_vetoWorldRects[1] = new Rectangle2D.Double(barLeft, barBot + 3
				* boxHeight + gap, 3 * boxWidth, gap);// internal
		_vetoWorldRects[7] = new Rectangle2D.Double(barLeft - 2 * gap, barBot
				+ 3 * boxHeight + 3 * gap, 2 * boxWidth + gap, gap);// external
																	// -
																	// upstream
		_vetoWorldRects[8] = new Rectangle2D.Double(barLeft + 1.5 * boxWidth,
				barBot + 3 * boxHeight + 3 * gap, 2 * boxWidth + gap, gap);// external
																			// -
																			// downstream

	}

	/**
	 * Draws the rectangle backgrounds (probably not necessary)
	 */
	private void setBeforeDraw() {
		// style for bar rects
		_barStyle = new Styled(X11Colors.getX11Color("Dark Blue"));
		_barStyle.setLineColor(Color.black);

		// use a before-drawer to bar dividers and labels
		_beforeDraw = new DrawableAdapter() {

			@Override
			public void draw(Graphics g, IContainer container) {
				for (int bar = 0; bar < GeoConstants.NUM_BAR; bar++) {
					WorldGraphicsUtilities.drawWorldRectangle(g, container,
							_barWorldRects[bar], _barStyle);
				}
				for (int veto = 0; veto < GeoConstants.NUM_VETOES; veto++) {
					WorldGraphicsUtilities.drawWorldRectangle(g, container,
							_vetoWorldRects[veto], _barStyle);
				}
			}

		};

		getContainer().setBeforeDraw(_beforeDraw);
	}

	/**
	 * Can be used to draw things based on the final layout, but unused at the
	 * current moment
	 */
	private void setAfterDraw() {
		IDrawable _afterDraw = new DrawableAdapter() {

			@Override
			public void draw(Graphics g, IContainer container) {
			}

		};
		getContainer().setAfterDraw(_afterDraw);
	}

	/**
	 * Creates the bar and veto instances that will handle and display hits
	 */
	private void addItems() {
		LogicalLayer detectorLayer = getContainer().getLogicalLayer(
				_detectorLayerName);

		_superLayerBars = new FullSideViewBar[GeoConstants.NUM_BAR];

		_superLayerVetoes = new FullSideViewVeto[GeoConstants.NUM_VETOES];

		for (int bar = 0; bar < GeoConstants.NUM_BAR; bar++) {
			_superLayerBars[bar] = new FullSideViewBar(detectorLayer, this,
					_barWorldRects[bar], bar);
		}
		for (int veto = 0; veto < GeoConstants.NUM_VETOES; veto++) {
			_superLayerVetoes[veto] = new FullSideViewVeto(detectorLayer, this,
					_vetoWorldRects[veto], veto);
		}
	}

	/**
	 * Some view specific feedback. Should always call super.getFeedbackStrings
	 * first.
	 * 
	 * Currently does not do anything but use the super call.
	 * 
	 * @param container
	 *            the base container for the view.
	 * @param screenPoint
	 *            the pixel point
	 * @param worldPoint
	 *            the corresponding world location.
	 */
	@Override
	public void getFeedbackStrings(IContainer container, Point screenPoint,
			Point2D.Double worldPoint, List<String> feedbackStrings) {
		super.getFeedbackStrings(container, screenPoint, worldPoint,
				feedbackStrings);
	}

	/**
	 * Gets which bar or veto the point is contained in.
	 * 
	 * @param container
	 *            the base container for the view.
	 * @param screenPoint
	 *            the pixel point
	 * @param worldPoint
	 *            the corresponding world location.
	 * @return the bar (1-9), veto(1-14), or -1 for none.
	 */
	@Override
	public int getSector(IContainer container, Point screenPoint,
			Point2D.Double worldPoint) {
		for (int bar = 0; bar < GeoConstants.NUM_BAR; bar++) {
			if (_barWorldRects[bar].contains(worldPoint)) {
				return bar + 1; // convert to 1-based index
			}
		}
		for (int veto = 0; veto < GeoConstants.NUM_VETOES; veto++) {
			if (_vetoWorldRects[veto].contains(worldPoint)) {
				return veto + 1;
			}
		}
		return -1;
	}
	
	public FullSideViewBar[] getBars() {
		return _superLayerBars;
	}
	
	public FullSideViewVeto[] getVetoes() {
		return _superLayerVetoes;
	}
}
