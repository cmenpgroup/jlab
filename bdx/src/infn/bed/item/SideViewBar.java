package infn.bed.item;

import infn.bed.bedview.BarSideView;
import infn.bed.event.ChargeTimeData;
import infn.bed.event.EventManager;
import infn.bed.frame.Bed;
import infn.bed.bedview.BedView;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import cnuphys.bCNU.event.EventControl;
import cnuphys.bCNU.format.DoubleFormat;
import cnuphys.bCNU.graphics.container.IContainer;
import cnuphys.bCNU.graphics.world.WorldGraphicsUtilities;
import cnuphys.bCNU.item.RectangleItem;
import cnuphys.bCNU.layer.LogicalLayer;
import cnuphys.bCNU.log.Log;
import cnuphys.bCNU.util.FileUtilities;
import cnuphys.bCNU.util.Fonts;

/**
 * This class handles converting charge-time information to energy-time
 * information. It then draws rectangles and colors them if they are hit. It
 * displays charge-time information for each bar in the info panel.
 * 
 * @author Andy Beiter
 * 
 */
public class SideViewBar extends RectangleItem {

	/**
	 * Font for label text
	 */
	private static final Font labelFont = Fonts.commonFont(Font.PLAIN, 11);

	/**
	 * The number of the bar
	 */
	private int _bar;

	/**
	 * The array of the hit sectors (detectors)
	 */
	private int hitSectors[];

	/**
	 * The array of the hit layers (columns)
	 */
	private int hitLayers[];

	/**
	 * The array of the hit paddles (rows)
	 */
	private int hitPaddles[];

	/**
	 * The array of the charges from the left PMTs
	 */
	private int chargeLeft[];

	/**
	 * The array of the charges from the right PMTs
	 */
	private int chargeRight[];

	/**
	 * The array of the times from the left PMTs
	 */
	private int timeLeft[];

	/**
	 * The array of the times from the right PMTs
	 */
	private int timeRight[];

	/**
	 * The array of energies from each hit
	 */
	private double totalE[];

	/**
	 * The array of times from each hit
	 */
	private double totalT[];

	/**
	 * The effective speed of light in the bar
	 */
	private double v_eff;

	/**
	 * The charge-to-energy conversion factor for the left PMT for the bar
	 */
	private double A_left;

	/**
	 * The charge-to-energy conversion factor for the right PMT for the bar
	 */
	private double A_right;

	/**
	 * The attenuation length
	 */
	private double lambda;

	/**
	 * Time delay for the left PMT
	 */
	private double delta_left;

	/**
	 * Time delay for the right PMT
	 */
	private double delta_right;

	/**
	 * The tdc-to-time conversion factor for the left PMT for the bar
	 */
	private double tdcConvLeft;

	/**
	 * The tdc-to-time conversion factor for the right PMT for the bar
	 */
	private double tdcConvRight;

	/**
	 * The length of the bar
	 */
	private double length;

	/**
	 * The view this bar is in
	 */
	private BarSideView _view;

	/**
	 * Color for hit cells
	 */
	private static final Color defaultHitCellFill = Color.red;

	/**
	 * The rectangle the bar is drawn in
	 */
	private Rectangle2D.Double _worldRectangle;

	/**
	 * Constructor for the bar used in the side view
	 * 
	 * @param layer
	 *            the Layer this item is on.
	 * @param view
	 *            the BarSideView parent
	 * @param worldRectangle
	 *            the rectangle the bar is in
	 * @param bar
	 *            the bar (0-8)
	 */
	public SideViewBar(LogicalLayer layer, BarSideView view,
			Rectangle2D.Double worldRectangle, int bar) {
		super(layer, worldRectangle);
		_worldRectangle = worldRectangle;
		_view = view;

		_style.setFillColor(Color.white);
		_style.setLineColor(Color.black);
		_bar = bar + 1; // convert to 1-based

		_name = "Bar: " + _bar;
		getConstants();
	}

	/**
	 * Gets the calibration constants from a file and stores them
	 * 
	 * TODO change file for real constants
	 */
	private void getConstants() {
		File file = FileUtilities.findFile(Bed.dataPath,
				"calibrationConstantsSimulation.dat");

		if ((file != null) && file.exists()) {
			Log.getInstance().info(
					"Calibration constants file found: " + file.getPath());
			try {
				FileReader fileReader = new FileReader(file);
				BufferedReader bufferedReader = new BufferedReader(fileReader);
				boolean notFound = true;
				while (notFound) {
					String s = bufferedReader.readLine();
					if (s == null) {
						break;
					} else {
						if (!s.startsWith("#") && (s.length() > 0)) {
							String tokens[] = FileUtilities.tokens(s);
							// if this line is me
							if (Integer.parseInt(tokens[0]) == _bar) {
								notFound = false;
								v_eff = Double.parseDouble(tokens[1]);
								A_left = Double.parseDouble(tokens[2]);
								A_right = Double.parseDouble(tokens[3]);
								lambda = Double.parseDouble(tokens[4]);
								delta_left = Double.parseDouble(tokens[5]);
								delta_right = Double.parseDouble(tokens[6]);
								tdcConvLeft = Double.parseDouble(tokens[7]);
								tdcConvRight = Double.parseDouble(tokens[8]);
								length = Double.parseDouble(tokens[9]);
							}
						}
					}
				}
				bufferedReader.close();
			} catch (FileNotFoundException e) {
				Log.getInstance().exception(e);
				e.printStackTrace();
			} catch (IOException e) {
				Log.getInstance().exception(e);
				e.printStackTrace();
			}
			Log.getInstance().info(
					"Successfully read in calibration constants.");
		} else {
			Log.getInstance().warning(
					"Calibration constants file not found at: "
							+ ((file == null) ? "???" : file.getPath()));
		}
	}

	/**
	 * Custom drawer for the bar.
	 * 
	 * @param g
	 *            the graphics context.
	 * @param container
	 *            the graphical container being rendered.
	 */
	@Override
	public void drawItem(Graphics g, IContainer container) {

		if (EventControl.getInstance().isAccumulating()) {
			return;
		}

		super.drawItem(g, container); // draws rectangular shell

		g.setFont(labelFont);
		g.setColor(Color.yellow);

		if (_view.getMode() == BedView.Mode.SINGLE_EVENT) {
			singleEventDrawItem(g, container);
		} else {
			accumulatedDrawItem(g, container);
		}

		g.setColor(_style.getLineColor());
		g.drawPolygon(_lastDrawnPolygon);
	}

	/**
	 * Draw in single event mode
	 * 
	 * @param g
	 *            the graphics context
	 * @param container
	 *            the rendering container
	 */
	private void singleEventDrawItem(Graphics g, IContainer container) {

		// draw default, blank rectangle
		WorldGraphicsUtilities.drawWorldRectangle(g, container,
				_worldRectangle, Color.white, _style.getLineColor());

		// get the data and make sure it's not null
		ChargeTimeData ctData = EventManager.getInstance().getChargeTimeData();
		if (ctData != null) {
			hitSectors = ctData.getSector();
			hitLayers = ctData.getLayer();
			hitPaddles = ctData.getPaddle();
			chargeLeft = ctData.getChargeLeft();
			timeLeft = ctData.getTimeLeft();
			chargeRight = ctData.getChargeRight();
			timeRight = ctData.getTimeRight();

			// if we have hits
			if (chargeLeft != null && chargeRight != null) {

				// convert to energy
				chargeToEnergy();

				for (int i = 0; i < totalE.length; i++) {

					// if the hits are in this bar
					if (inThisBar(hitSectors[i], hitLayers[i], hitPaddles[i])) {

						// if the energy is above 0 (extra check)
						if (totalE[i] > 0) {

							// draw red rectangle
							_style.setFillColor(defaultHitCellFill);
							WorldGraphicsUtilities.drawWorldRectangle(g,
									container, _worldRectangle,
									_style.getFillColor(),
									_style.getLineColor());
						}
					}
				}
			}
		}
	}

	/**
	 * Converts the charge-time information to energy-time information
	 */
	private void chargeToEnergy() {
		double t_l[] = new double[timeLeft.length];
		double t_r[] = new double[timeRight.length];
		// change tdc to time
		for (int i = 0; i < timeLeft.length; i++) {
			t_l[i] = (timeLeft[i] * 1.0 / tdcConvLeft) - delta_left;
		}
		for (int i = 0; i < timeRight.length; i++) {
			t_r[i] = (timeRight[i] * 1.0 / tdcConvRight) - delta_right;
		}
		// get the position of the hit from the left side
		double posFromLeft[] = new double[chargeLeft.length];
		for (int i = 0; i < posFromLeft.length; i++) {
			posFromLeft[i] = (v_eff * (t_l[i] - t_r[i]) + length) / 2.0;
		}
		// change charge to energy
		double e_l[] = new double[chargeLeft.length];
		double e_r[] = new double[chargeRight.length];
		for (int i = 0; i < e_l.length; i++) {
			e_l[i] = chargeLeft[i] * A_left;
		}
		for (int i = 0; i < e_r.length; i++) {
			e_r[i] = chargeRight[i] * A_right;
		}
		totalE = new double[e_l.length];
		totalT = new double[t_l.length];
		for (int i = 0; i < totalE.length; i++) {
			double e_l_prime = e_l[i] * Math.exp(posFromLeft[i] / lambda);
			double e_r_prime = e_r[i]
					* Math.exp((length - posFromLeft[i]) / lambda);
			totalE[i] = (e_l_prime + e_r_prime) / 2;
			totalT[i] = (t_l[i] + t_r[i] - (length / v_eff)) / 2.0;
		}
	}

	/**
	 * Checks if a given sector/layer/paddle is this bar
	 * 
	 * @param sector
	 *            Which detector (currently not used)
	 * @param layer
	 *            Which column
	 * @param paddle
	 *            Which row
	 * @return True if the s/l/p is this bar, false if not
	 */
	private boolean inThisBar(int sector, int layer, int paddle) {
		switch (layer) {
		case 0:
			switch (paddle) {
			case 0:
				return (7 == _bar);
			case 1:
				return (4 == _bar);
			case 2:
				return (1 == _bar);
			}
			break;
		case 1:
			switch (paddle) {
			case 0:
				return (8 == _bar);
			case 1:
				return (5 == _bar);
			case 2:
				return (2 == _bar);
			}
			break;
		case 2:
			switch (paddle) {
			case 0:
				return (9 == _bar);
			case 1:
				return (6 == _bar);
			case 2:
				return (3 == _bar);
			}
		}
		return false;
	}

	/**
	 * Draw hits in accumulated mode
	 * 
	 * Currently unused
	 * 
	 * @param g
	 *            the graphics context
	 * @param container
	 *            the rendering container
	 */
	private void accumulatedDrawItem(Graphics g, IContainer container) {
	}

	/**
	 * Add any appropriate feedback strings for the heads-up display or feedback
	 * panel.
	 * 
	 * @param container
	 *            the Base container.
	 * @param screenPoint
	 *            the mouse location.
	 * @param worldPoint
	 *            the corresponding world point.
	 * @param feedbackStrings
	 *            the List of feedback strings to add to.
	 */
	@Override
	public void getFeedbackStrings(IContainer container, Point screenPoint,
			Point2D.Double worldPoint, List<String> feedbackStrings) {
		if (_worldRectangle.contains(worldPoint)) {

			double x = 0;
			double y = worldPoint.y;
			double z = 3 - worldPoint.x;
			z *= 10;
			y *= 10;

			String rtp = "approx xyz " + DoubleFormat.doubleFormat(x, 1)
					+ "cm, " + DoubleFormat.doubleFormat(y, 1) + "cm, "
					+ DoubleFormat.doubleFormat(z, 1) + "cm";
			feedbackStrings.add(rtp);

			if (_view.getMode() == BedView.Mode.SINGLE_EVENT) {
				singleEventFeedbackStrings(feedbackStrings);
			} else {
				accumulatedFeedbackStrings(feedbackStrings);
			}

		}
	}

	/**
	 * Get the feedback strings for single event mode. Displays charge left,
	 * charge right, time left, and time right for each hit.
	 * 
	 * @param feedbackStrings
	 *            The list of feedback strings
	 */
	private void singleEventFeedbackStrings(List<String> feedbackStrings) {
		if (EventManager.getInstance().getChargeTimeData() != null) {
			if (chargeLeft != null && chargeRight != null) {
				String energyStr = "";
				for (int i = 0; i < totalE.length; i++) {
					if (inThisBar(hitSectors[i], hitLayers[i], hitPaddles[i])) {
						if (totalE[i] > 0) {
							energyStr += "$orange$" + "Left PMT Charge:  "
									+ chargeLeft[i] + "\nLeft PMT Time:  "
									+ timeLeft[i] + "\nRight PMT Charge:  "
									+ chargeRight[i] + "\nRight PMT Time:  "
									+ timeRight[i];
						}
					}
				}
				feedbackStrings.add(energyStr);
			}
		}
	}

	/**
	 * Get the feedback strings for accumulated mode
	 * 
	 * Currently unused
	 * 
	 * @param feedbackStrings
	 *            The list of feedback strings
	 */
	private void accumulatedFeedbackStrings(List<String> feedbackStrings) {
	}
}
