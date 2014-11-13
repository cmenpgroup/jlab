package infn.bed.event;

import infn.bed.frame.Bed;
import infn.bed.geometry.GeoConstants;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Vector;

import org.jlab.coda.jevio.IEvioStructure;

import cnuphys.bCNU.log.Log;
import cnuphys.bCNU.util.FileUtilities;
import cnuphys.lund.LundId;

/**
 * This class handles the reading in of a charge-time data file. It also
 * converts the full waveform data to charge-time data in a special constructor.
 * This allows the view classes to base their code off this class.
 * 
 * @author Andy Beiter
 * 
 */
public class ChargeTimeData implements ILoad {

	/**
	 * An array of the hit sectors (which detectors)
	 */
	private int sectorArr[];

	/**
	 * An array of the hit layers (columns in the detector)
	 */
	private int layerArr[];

	/**
	 * An array of the hit paddles (rows in the detectors)
	 */
	private int paddleArr[];

	/**
	 * An array of the charges from hits in the left PMTs
	 */
	private int chargeLeft[];

	/**
	 * An array of the charges from hits in the right PMTs
	 */
	private int chargeRight[];

	/**
	 * An array of the times from hits in the left PMTs
	 */
	private int timeLeft[];

	/**
	 * An array of the times from hits in the right PMTs
	 */
	private int timeRight[];

	/**
	 * An array of the hit veto sectors (which detector)
	 */
	private int vetoSector[];

	/**
	 * An array indicating if the hit veto was internal or external
	 */
	private int vetoIntOrExt[];

	/**
	 * An array representing which of the vetoes was hit. For numbering, see
	 * FullSideView.java
	 */
	private int vetoChannel[];

	/**
	 * An array of the charges from hits in the vetoes.
	 */
	private int vetoCharge1[];

	/**
	 * An array of the charges from hits in the vetoes if the vetoes use two
	 * PMTs. It should only contain data for the top and bottom external vetoes.
	 */
	private int vetoCharge2[];

	/**
	 * An array of the times from hits in the vetoes.
	 */
	private int vetoTime1[];

	/**
	 * An array of the times from hits in the vetoes if the vetoes use two PMTs.
	 * It should only contain data for the top and bottom external vetoes.
	 */
	private int vetoTime2[];

	/**
	 * Constructor used for charge-time files
	 */
	public ChargeTimeData() {
	}

	/**
	 * Constructor used for full waveshape files. Converts the waveshape to
	 * charge-time information.
	 * 
	 * @param samples
	 *            The waveshape from the PMTs
	 */
	public ChargeTimeData(ArrayList<ArrayList<Short>> samples) {
		ArrayList<Double> leftCharges = new ArrayList<Double>();
		ArrayList<Double> leftTimes = new ArrayList<Double>();
		ArrayList<Double> rightCharges = new ArrayList<Double>();
		ArrayList<Double> rightTimes = new ArrayList<Double>();
		ArrayList<Integer> sectors = new ArrayList<Integer>();
		ArrayList<Integer> layers = new ArrayList<Integer>();
		ArrayList<Integer> paddles = new ArrayList<Integer>();

		// TODO this can be cleaned up, I think
		// relies on detectorTable.dat which converts from the index in samples
		// to sector/layer/paddle or sector/intOrExt/channel
		for (int i = 0; i < samples.size(); i++) {

			// TODO convert to charge for a veto
			// same algorithm except for double i increment?
			if (i < (GeoConstants.NUM_BAR * 2)) { // if a bar
				ArrayList<Short> leftSamples = samples.get(i); // assumes left
																// info is first
				ArrayList<Short> rightSamples = samples.get(i + 1);

				// for our timing algorithm, see fadc_time.pdf in docs

				// declare time conversion variables
				File file = FileUtilities.findFile(Bed.dataPath,
						"detectorTable.dat"); // translation table as mentioned
												// above
				int index = -1, sector = -1, layer = -1, paddle = -1;
				if ((file != null) && file.exists()) {
					Log.getInstance().info(
							"Calibration constants file found: "
									+ file.getPath());
					try {
						FileReader fileReader = new FileReader(file);
						BufferedReader bufferedReader = new BufferedReader(
								fileReader);
						boolean notFound = true;
						while (notFound) {
							String s = bufferedReader.readLine();
							if (s == null) {
								break;
							} else {
								if (!s.startsWith("#") && (s.length() > 0)) {
									String tokens[] = FileUtilities.tokens(s);
									index = Integer.parseInt(tokens[0]);
									sector = Integer.parseInt(tokens[1]);
									layer = Integer.parseInt(tokens[2]);
									paddle = Integer.parseInt(tokens[3]);
									if (index == i) {
										notFound = false; // we found our info
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
				} else {
					Log.getInstance()
							.warning(
									"Translation table file not found at: "
											+ ((file == null) ? "???" : file
													.getPath()));
					System.out.println("Failed to read in translation table!");
				}
				/*
				 * TODO The current charge algorithm is to just sum the ADC
				 * values divided by the resistance multiplied by the time which
				 * should be 4x the index. This is subject to change.
				 */
				int leftHits = convertHits(leftSamples, leftCharges, leftTimes);
				for(int hit = 0; hit < leftHits; hit++) {
					sectors.add(sector);
					layers.add(layer);
					paddles.add(paddle);
				}
				int rightHits = convertHits(rightSamples, rightCharges, rightTimes);
				for(int hit = 0; hit < rightHits; hit++) {
					sectors.add(sector);
					layers.add(layer);
					paddles.add(paddle);
				}
				i++; // must increment twice since we're doing two bars
			}
		}

		// set arrays to final charge time data
		sectorArr = getIntArray(sectors);
		layerArr = getIntArray(layers);
		paddleArr = getIntArray(paddles);
		chargeLeft = getIntArrayDoubleCast(leftCharges);
		chargeRight = getIntArrayDoubleCast(rightCharges);
		timeLeft = getIntArrayDoubleCast(leftTimes);
		timeRight = getIntArrayDoubleCast(rightTimes);
	}
	
	private int convertHits(ArrayList<Short> samples, ArrayList<Double> charges, ArrayList<Double> times) {
		int numHits = 0;
		double a_L = 0, b_L = 0, charge = 0, time = 0, threshold = 0;
		double fADCResistance = 50.0; // in ohms, resistance of fADC
		boolean collectingPulse = false; // if we have a hit
		threshold = getThreshold(); // gets threshold to look for hits
									// above
		for (int j = 1; j < (samples.size() - 1); j++) {
			if (samples.get(j) > threshold
					&& samples.get(j - 1) < threshold) {
				a_L = samples.get(j + 1) - samples.get(j - 1)
						* 1.0 / 4.0;
				b_L = samples.get(j + 1) - a_L * (j - 1) * 4;
				charge += (samples.get(j) / fADCResistance)
						* (j - 1) * 4;
				collectingPulse = true;
			} else if ((samples.get(j + 1) < samples.get(j))
					&& (samples.get(j - 1) < samples.get(j))
					&& (samples.get(j) > threshold)) {
				time = samples.get(j) / 2.0;
				time -= b_L;
				time /= a_L;
				charge += (samples.get(j) / fADCResistance)
						* (j - 1) * 4;
			} else if ((samples.get(j) > threshold)
					&& (samples.get(j + 1) < threshold)) { // end of
																// hit
				charge += (samples.get(j) / fADCResistance)
						* (j - 1) * 4;
				// add and reset
				charges.add(charge);
				times.add(time);
				numHits ++;
				a_L = 0;
				b_L = 0;
				time = 0;
				charge = 0;
				collectingPulse = false;
			} else if (collectingPulse) { // collecting, but not the
											// extremes
				charge += (samples.get(j) / fADCResistance)
						* (j - 1) * 4;
			}
		}
		return numHits;
	}

	/**
	 * Returns the threshold for the ADC. Should be in ADC channel units
	 * (uncalibrated)
	 * 
	 * @return The ADC channel above which hits are detected.
	 */
	private double getThreshold() {
		return 300;	//TODO change threshold for final product
	}

	/**
	 * Converts ArrayList of Integers to an int array
	 * 
	 * @param ints
	 *            ArrayList of Integers
	 * @return primitive array of ints
	 */
	private int[] getIntArray(ArrayList<Integer> ints) {
		int primArr[] = new int[ints.size()];
		for (int i = 0; i < primArr.length; i++) {
			primArr[i] = ints.get(i);
		}
		return primArr;
	}

	/**
	 * Converts ArrayList of Doubles to an int array
	 * 
	 * @param doubles
	 *            ArrayList of Doubles
	 * @return primitive array of ints
	 */
	private int[] getIntArrayDoubleCast(ArrayList<Double> doubles) {
		int primArr[] = new int[doubles.size()];
		for (int i = 0; i < primArr.length; i++) {
			// must do two steps
			double temp = doubles.get(i);
			primArr[i] = (int) temp;
		}
		return primArr;
	}

	/**
	 * Loads in the charge-time data from a charge-time file. Tag 102 is for
	 * bars and tag 202 is for vetoes. Numbers 1-7 gets the sector,
	 * layer/intOrExt, paddle/channel, chargeLeft/1, chargeRight/2, timeLeft/1,
	 * and timeRight/2 arrays, respectively.
	 * 
	 * @param structure
	 *            The data
	 * @param tag
	 *            The tag used to identify the data
	 * @param num
	 *            The number used to identify the data
	 */
	@Override
	public void load(IEvioStructure structure, int tag, int num) {
		try {
			if (tag == 102) {
				switch (num) {
				case 1:
					sectorArr = structure.getIntData();
					break;
				case 2:
					layerArr = structure.getIntData();
					break;
				case 3:
					paddleArr = structure.getIntData();
					break;
				case 4:
					chargeLeft = structure.getIntData();
					break;
				case 5:
					chargeRight = structure.getIntData();
					break;
				case 6:
					timeLeft = structure.getIntData();
					break;
				case 7:
					timeRight = structure.getIntData();
					break;
				}
			} else if (tag == 202) {
				switch (num) {
				case 1:
					vetoSector = structure.getIntData();
					break;
				case 2:
					vetoIntOrExt = structure.getIntData();
					break;
				case 3:
					vetoChannel = structure.getIntData();
					break;
				case 4:
					vetoCharge1 = structure.getIntData();
					break;
				case 5:
					vetoCharge2 = structure.getIntData();
					break;
				case 6:
					vetoTime1 = structure.getIntData();
					break;
				case 7:
					vetoTime2 = structure.getIntData();
					break;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Currently unused.
	 * 
	 * @see infn.bed.event.ILoad#uniqueLundIds()
	 */
	@Override
	public Vector<LundId> uniqueLundIds() {
		return null;
	}

	/**
	 * Gets the hit sector array
	 * 
	 * @return The array of hit sectors
	 */
	public int[] getSector() {
		return sectorArr;
	}

	/**
	 * Gets the hit layer array
	 * 
	 * @return The array of hit layers
	 */
	public int[] getLayer() {
		return layerArr;
	}

	/**
	 * Gets the hit sector paddle
	 * 
	 * @return The array of hit paddles
	 */
	public int[] getPaddle() {
		return paddleArr;
	}

	/**
	 * Gets the charges from the left PMT
	 * 
	 * @return The array of charges from the left PMT
	 */
	public int[] getChargeLeft() {
		return chargeLeft;
	}

	/**
	 * Gets the charges from the right PMT
	 * 
	 * @return The array of charges from the right PMT
	 */
	public int[] getChargeRight() {
		return chargeRight;
	}

	/**
	 * Gets the times from the left PMT
	 * 
	 * @return The array of times from the left PMT
	 */
	public int[] getTimeLeft() {
		return timeLeft;
	}

	/**
	 * Gets the times from the right PMT
	 * 
	 * @return The array of times from the right PMT
	 */
	public int[] getTimeRight() {
		return timeRight;
	}

	/**
	 * Gets the hit veto sector array
	 * 
	 * @return The array of hit veto sectors
	 */
	public int[] getVetoSector() {
		return vetoSector;
	}

	/**
	 * Gets the hit intOrExt array
	 * 
	 * @return The array of hit intOrExts
	 */
	public int[] getVetoIntOrExt() {
		return vetoIntOrExt;
	}

	/**
	 * Gets the hit veto channel array
	 * 
	 * @return The array of hit veto channels
	 */
	public int[] getVetoChannel() {
		return vetoChannel;
	}

	/**
	 * Gets the veto charges from the first PMT
	 * 
	 * @return The array of veto charges from the first PMT
	 */
	public int[] getVetoCharge1() {
		return vetoCharge1;
	}

	/**
	 * Gets the veto charges from the second PMT
	 * 
	 * @return The array of veto charges from the second PMT
	 */
	public int[] getVetoCharge2() {
		return vetoCharge2;
	}

	/**
	 * Gets the veto times from the first PMT
	 * 
	 * @return The array of veto times from the first PMT
	 */
	public int[] getVetoTime1() {
		return vetoTime1;
	}

	/**
	 * Gets the veto times from the second PMT
	 * 
	 * @return The array of veto times from the second PMT
	 */
	public int[] getVetoTime2() {
		return vetoTime2;
	}

}
