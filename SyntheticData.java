package synthesis;
import org.apache.commons.math3.distribution.BetaDistribution;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.stream.Collectors;
import static java.lang.Math.abs;


public class SyntheticData {

	static String[] filenames = {"parameters_int_gauss_pca_0413.csv"}; // parameters_a_rt_2

	static class Spectra implements Comparable<Spectra> {
		double rt;
		int rtIndex;
		ArrayList<double[]> mzIntensity; // [mz, intensity] pairs

		Spectra(double rt, int rtIndex, ArrayList<double[]> mzIntensity) {
			this.rt = rt;
			this.rtIndex = rtIndex;
			this.mzIntensity = mzIntensity;
		}

		// returns true if this spectra contains mz looking for
		public int containsMz(double mz) {
			int i = 0;
			for (double[] pair : this.mzIntensity){
				if (pair.length > 0 && doubleEqual(mz, pair[0]))
					return i;
				i++;
			}
			return -1;
		}

		@Override
		public int compareTo(Spectra o) {
			return Double.compare(this.rt, o.rt);
		}
	}

	static class Peak {
		double rtStart;
		int rtStartIndex;
		double rtEnd;
		int rtEndIndex;
		double mzMin;
		double mzMax;
		ArrayList<Spectra> spectras;

		public Peak(double rtStart, int rtStartIndex,double rtEnd, int rtEndIndex, double mzMin, double mzMax, ArrayList<Spectra> spectras) {
			this.rtStart = rtStart;
			this.rtStartIndex = rtStartIndex;
			this.rtEnd = rtEnd;
			this.rtEndIndex = rtEndIndex;
			this.mzMin = mzMin;
			this.mzMax = mzMax;
			this.spectras = spectras;
		}

		public Peak copy() {
			ArrayList<Spectra> newSpect = new ArrayList<>();
			for (Spectra s : this.spectras){
				if (s.mzIntensity.size() > 0){
					ArrayList<double[]> mzints = new ArrayList<>();
					for (double[] doubles : s.mzIntensity){
						double[] mzint = new double[]{doubles[0], doubles[1]};
						mzints.add(mzint);
					}
					Spectra spect = new Spectra(s.rt, s.rtIndex, mzints);
					newSpect.add(spect);
				}
			}
			return new Peak(this.rtStart, this.rtStartIndex, this.rtEnd, this.rtEndIndex, this.mzMin, this.mzMax, newSpect);
		}

		// EFFECT: sets rt of each peak in the feature to rt at current idex + offset and
		//		   update index
		// Note: peak fields never updated, but will not influence outcome.
		public void shiftRt(int offset, List<double[]> rtData) {
			if (this.spectras != null) {
				for (Spectra spectra : this.spectras) {
					if (spectra.rtIndex + offset >= 0 && spectra.rtIndex + offset < rtData.size()) {
						spectra.rtIndex += offset;
						spectra.rt = rtData.get(spectra.rtIndex)[1];
					}
				}
			}
		}

		// EFFECT: changes width of feature from rt to (rt-mid)*proportion where mid = (max rt -min rt)/2
		public void scalePeakWidth(double proportion, List<double[]> rtData) {
			double minRt = this.spectras.get(0).rt;
			double maxRt = this.spectras.get(spectras.size()-1).rt;
			double mid = (minRt + maxRt)/2;
			// if scaled peak is outside of boundary of rtData, skip
			if ((minRt - mid)*proportion + mid < rtData.get(0)[1] ||
					(maxRt - mid)*proportion + mid > rtData.get(rtData.size()-1)[1]){
				return;
			}
			// for each peak in metabolic feature, change its rt value and rt index (index of rt value in raw data)
			for (Spectra spect : this.spectras) {
				double oldRt = spect.rt;
				double newRt = (spect.rt - mid)*proportion + mid;
				int oldidx = spect.rtIndex;
				int newidx = spect.rtIndex;
				double mindist = Double.POSITIVE_INFINITY;
				// find new index and new rt
				if (newRt < oldRt) {
					for (int i = oldidx; i >= 0; i--) {
						double dist = abs(rtData.get(i)[1]-newRt);
						if (mindist < dist){
							break;
						} else {
							mindist = dist;
							newidx = i;
						}
					}
				} else {
					for (int i = oldidx; i < rtData.size(); i++) {
						double dist = abs(rtData.get(i)[1]-newRt);
						if (mindist < dist){
							break;
						} else {
							mindist = dist;
							newidx = i;
						}
					}
				}
				// set new values
				spect.rtIndex = newidx;
				spect.rt = rtData.get(newidx)[1];
			}
		}

		// EFFECT: changes the intensity of each peak in metabolomic feature from intensity to intensity * proportion
		public void scaleIntensity(double proportion) {
			for (Spectra spect : this.spectras) {
				for (double[] arr : spect.mzIntensity) {
					arr[1] = arr[1]*proportion;
				}
			}
		}
	}

	public static void main(String[] args) throws Exception {

		for (String file : filenames){


			List<List<String>> user_input = readParameters(file);
			List<String> samples = user_input.get(0); // "CD-9OS5Y",
			List<String> copies = user_input.get(1);

			List<String> rt_model = user_input.get(2);
			List<String> int_model = user_input.get(3);
			List<String> pw_model = user_input.get(4);

			double massTol = 0.01;
			int num_models = 4;

			for (int x = 0; x < samples.size(); x ++) {


				// Read input files
				List<double[]> intensitiesData = readInput(samples.get(x) + "_intensity.csv");
				List<double[]> mzsData = readInput(samples.get(x) + "_mz.csv");
				List<double[]> regions = readInput(samples.get(x) + "_features.csv");
				List<double[]> rtData = readInput(samples.get(x) + "_rt.csv");

				// mzsData arrays same length (pad with 0)
				boolean[][] read = new boolean[mzsData.size()][mzsData.get(0).length];
				int[][] readed = new int[mzsData.size()][mzsData.get(0).length];

				// Construct peaks from input files
				LinkedList<Peak> peakList = new LinkedList<>();
				for (int i = 0; i < regions.size(); i++) {
					double rtStart = regions.get(i)[4];
					double rtEnd = regions.get(i)[5];
					double mzMin = regions.get(i)[1];
					double mzMax = regions.get(i)[2];

					int rowStartIndex = getRtRowIndex(rtData, rtStart);
					int rowEndIndex = getRtRowIndex(rtData, rtEnd);
//					if (abs(mzMin - 283.08) < 0.01){
//						int akc = 0;
//					}
					Peak grid = getPeak(rtData, mzsData, intensitiesData, rowStartIndex, rowEndIndex, mzMin, mzMax, read, readed);
					peakList.add(grid);
				}

				boolean rt_random = rt_model.get(rt_model.size()-1).equals("b");
				boolean int_random = int_model.get(int_model.size()-1).equals("b");
				boolean pw_random = pw_model.get(pw_model.size()-1).equals("b");
				int[] rt_feature_assignment = new int[0];
				int[] int_feature_assignment = new int[0];
				int[] pw_feature_assignment = new int[0];

				if (rt_random) {
					rt_feature_assignment = getFeatureAssignment(peakList.size(), num_models, Double.parseDouble(rt_model.get(rt_model.size()-2))); // TODO: coupling: assumes na index 5
				}
				if (pw_random) {
					pw_feature_assignment = getFeatureAssignment(peakList.size(), num_models, Double.parseDouble(pw_model.get(pw_model.size()-2)));
				}
				if (int_random) {
					int_feature_assignment = getFeatureAssignment(peakList.size(), num_models, Double.parseDouble(int_model.get(int_model.size()-2)));
				}

				// todo:
				ArrayList<ArrayList<String>> sample_factors = new ArrayList<>();
				for (int y = 1; y <= Integer.parseInt(copies.get(x)); y++){

					double offset = 0;
					int newIdx = 0;
					double proportion = 0;

					ArrayList<String> histo = new ArrayList<>();
					ArrayList<Spectra> spectraNew = new ArrayList<>();
					for (int i = 0; i <peakList.size(); i++) {
						if (peakList.get(i).spectras.isEmpty())
							continue;
						Peak newPeak = peakList.get(i).copy();

						int model = -1;
						if (rt_random) {
							model = rt_feature_assignment[i];
							List<String> rand_model = new ArrayList<>();
							rand_model.add(0, String.valueOf(model));
							rand_model.add(1, rt_model.get(model)); 	// rand_model stores parameters
							offset = getOffset(model, rand_model, y, Integer.parseInt(copies.get(x))); 			// note: if random, user need to put place holder 0
						} else {
							model = Integer.parseInt(rt_model.get(0));
							offset = getOffset(model, rt_model, y, Integer.parseInt(copies.get(x))); 			// rt_model stores parameters
						}
						if (offset == Integer.MIN_VALUE)
							throw new Exception("Invalid model number supplied. ");
						double oldRT = newPeak.spectras.get(0).rt;
						int oldIdx = newPeak.spectras.get(0).rtIndex;
						newIdx = findNearestRTIdx(oldRT + offset, oldIdx, rtData);
						newPeak.shiftRt(newIdx - oldIdx, rtData);


						if (pw_random){
							model = pw_feature_assignment[i];
							List<String> rand_model = new ArrayList<>();
							rand_model.add(0, String.valueOf(model));
							rand_model.add(1, pw_model.get(model));
							proportion = getProportion(model, rand_model, y, Integer.parseInt(copies.get(x)));
						} else {
							model = Integer.parseInt(pw_model.get(0));
							proportion = getProportion(model, pw_model, y, Integer.parseInt(copies.get(x)));
						}
						if (proportion == Integer.MIN_VALUE)
							throw new Exception("Invalid model number supplied. ");
						newPeak.scalePeakWidth(proportion, rtData); // stepped <-- is correct

						if (int_random){
							model = int_feature_assignment[i];
							List<String> rand_model = new ArrayList<>();
							rand_model.add(0, String.valueOf(model));
							rand_model.add(1, int_model.get(model));
							proportion = getProportion(model, rand_model, y, Integer.parseInt(copies.get(x)));
						} else {
							model = Integer.parseInt(int_model.get(0));
							proportion = getProportion(model, int_model, y, Integer.parseInt(copies.get(x)));
						}
						if (proportion == Integer.MIN_VALUE)
							throw new Exception("Invalid model number supplied. ");
						newPeak.scaleIntensity(proportion);
						// Note: may be good to refactor
						histo.add(String.valueOf(proportion));

						spectraNew.addAll(newPeak.spectras);
					}

					// todo:
//					String[] myArrays = histo.toArray(new String[0]);
					sample_factors.add(histo);

					// Merge background into spectraNew
//					ArrayList<Spectra> spectraNewTwo = spectraNew;


//					if (y == 1 || y == 2) {
//						Double[] myArray = histo.toArray(new Double[0]);
//						int nRanges = 100;
//						int[] buckets = new int[nRanges];
//						double max = Collections.max(Arrays.asList(myArray));
//						double min = Collections.min(Arrays.asList(myArray));
//						double sizeOfRange = (max - min) / (nRanges - 1);
//
//						for (double elem : myArray) {
//							for (int i = 0; i < nRanges; i++) {
//								if ((elem >= sizeOfRange * i) && (elem < sizeOfRange * (i + 1)))
//									buckets[i]++;
//							}
//						}
//						for (int i = 0; i < nRanges; i++){
//							System.out.println(sizeOfRange * i + " - " + sizeOfRange * (i + 1) + ": " + buckets[i]);
//						}
//					}

					sortSpectras(spectraNew, massTol);
					mergeBackground(intensitiesData, mzsData, rtData, read, spectraNew);


					// sort spectraNew
					sortSpectras(spectraNew, massTol);

					// Convert spectras to txt and output
					writeRT(spectraNew, samples.get(x), y);
					writeMzIntensity(spectraNew, samples.get(x), y);
				}

				// todo:
				String[][] stringArray = sample_factors.stream().map(u -> u.toArray(new String[0])).toArray(String[][]::new);

				List<String> list = Arrays.stream(stringArray).map(line -> String.join(",", line)).collect(Collectors.toList());
				try {
					Files.write(Paths.get(System.getProperty("user.dir"),"out.csv"), list, StandardOpenOption.CREATE);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}

	private static int findNearestRTIdx(double new_rt, int old_idx, List<double[]> rtData){
		double newRt = new_rt;
		double oldRt = rtData.get(old_idx)[1];
		int oldidx = old_idx;
		int newidx = old_idx;
		double mindist = Double.POSITIVE_INFINITY;

		if (newRt < rtData.get(0)[1] || newRt > rtData.get(rtData.size()-1)[1]){
			return old_idx;
		}

		if (newRt < oldRt) {
			for (int i = oldidx; i >= 0; i--) {
				double dist = abs(rtData.get(i)[1]-newRt);
				if (mindist < dist){
					break;
				} else {
					mindist = dist;
					newidx = i;
				}
			}
		} else {
			for (int i = oldidx; i < rtData.size(); i++) {
				double dist = abs(rtData.get(i)[1]-newRt);
				if (mindist < dist){
					break;
				} else {
					mindist = dist;
					newidx = i;
				}
			}
		}
		return newidx;
	}

	// model: 		user model choice
	// rt_model:	user defined model parameters
	// curr:		curr simulation number
	// total: 		total number of simulations needed for this sample
	// effect:		returns offset for different models based on user's input parameters
	//				default mean and stdev for linear is 0.1 and 0.1
	private static double getProportion(int model, List<String> user_model, int curr, int total) throws Exception {
		switch(model){
			case 1: // constant
				return Double.parseDouble(user_model.get(1));
			case 2: // linear
				String line = user_model.get(1).replaceAll("\\(", "").replaceAll("\\)", "");
				String[] linear_params = line.split("/");
				double slope = Double.parseDouble(linear_params[0]);
				double min = Double.parseDouble(linear_params[1]);
				double max = Double.parseDouble(linear_params[2]);
				if (min < 0)
					throw new Exception("Invalid min input for linear: min can't be negative");
				if (max - min < 0)
					throw new Exception("Invalid range for linear model: max can't be less than min");
				Random rand = new Random();
				double default_mean = 0.1;
				double default_stev = 0.1;// todo
				return min + (curr+1.0-1.0)* (max-min)/total + rand.nextGaussian()*default_stev + default_mean;
			case 3: // gaussian
				line = user_model.get(1).replaceAll("\\(", "").replaceAll("\\)", "");
				String[] gauss_params = line.split("/");
				rand = new Random();
				double mean = Double.parseDouble(gauss_params[0]);
				double stev = Double.parseDouble(gauss_params[1]);
				return abs(rand.nextGaussian()*stev + mean); // assumption!
			case 4: // beta
				line = user_model.get(1).replaceAll("\\(", "").replaceAll("\\)", "");
				String[] beta_params = line.split("/");
				double alpha = Double.parseDouble(beta_params[0]);
				double beta = Double.parseDouble(beta_params[1]);
				double range = Double.parseDouble(beta_params[2]);
				BetaDistribution betaDistribution = new BetaDistribution(alpha, beta);
				return betaDistribution.sample()*range;
			case 5: // no change
				return 1;
		}
		return Double.MIN_VALUE;
	}

	// size: 	how many peaks there are in the sample (# model = # peaks)
	// bound:	number of available models to draw from (not including NA, NA dealt separately)
	// theta: 	% NA
	// returns:	assignment of random models to each peak as int array
	private static int[] getFeatureAssignment(int size, int bound, double theta) {
		Random random = new Random();
		int[] assignment = new int[size];
		for (int i = 0; i < size; i++) {
			if (random.nextDouble() < theta) {
				assignment[i] = 5; // no change is case 5
			} else {
				assignment[i] = random.nextInt(bound) + 1; // 0-3 --> 1-4
			}
		}
		return assignment;
	}

	// model: 		user model choice
	// rt_model:	user defined model parameters
	// curr: 		number of current simulation
	// t: 			total number of simulations for this sample
	// effect:		returns time in seconds for different models based on user's input parameters
	// 				default mean and stdev for linear is 1 and 2, respectively
	private static double getOffset(int model, List<String> rt_model, int curr, int t) {
		switch (model) {
			case 1: // constant
				double seconds = Double.parseDouble(rt_model.get(1));
				return seconds;
			case 2: // linear
				String line = rt_model.get(1).replaceAll("\\(", "").replaceAll("\\)", "");
				String[] linear_params = line.split("/");
				double slope = Double.parseDouble(linear_params[0]);
				double min = Double.parseDouble(linear_params[1]);
				double max = Double.parseDouble(linear_params[2]);
				Random rand = new Random();
				double default_mean = 5;
				double default_stev = 2;
				return (curr+1.0-1.0)*(max-min)/t + rand.nextGaussian()*default_stev + default_mean;
//				double x = rand.nextDouble()*(max-min)+min
			case 3: // Gaussian
				line = rt_model.get(1).replaceAll("\\(", "").replaceAll("\\)", "");
				String[] gauss_params = line.split("/");
				rand = new Random();
				double mean = Double.parseDouble(gauss_params[0]);
				double stev = Double.parseDouble(gauss_params[1]);
				return rand.nextGaussian()*stev + mean;
			case 4: // Beta
				line = rt_model.get(1).replaceAll("\\(", "").replaceAll("\\)", "");
				String[] beta_params = line.split("/");
				double alpha = Double.parseDouble(beta_params[0]);
				double beta = Double.parseDouble(beta_params[1]);
				double range = Double.parseDouble(beta_params[2]);
				BetaDistribution betaDistribution = new BetaDistribution(alpha, beta);
				double ratio = betaDistribution.sample();
				return ratio*range;
			case 5: // no change
				return 0;
		}
		return Integer.MIN_VALUE;
	}

	// EFFECT: join all spectra while merging peaks with rt and mass tolerance
	private static void sortSpectras(ArrayList<Spectra> spectraNew, double massTol) {
		Collections.sort(spectraNew);
		for (int i = 0; i < spectraNew.size(); i++) {
			boolean ifSort = false;
			while (i + 1 < spectraNew.size() && doubleEqual(spectraNew.get(i).rt, spectraNew.get(i + 1).rt)) {
				spectraNew.get(i).mzIntensity.addAll(spectraNew.get(i + 1).mzIntensity);
				spectraNew.remove(i + 1);
				ifSort = true;
			}
			if (ifSort) {
				Collections.sort(spectraNew.get(i).mzIntensity, new Comparator<double[]>() {
					public int compare(double[] arr1, double[] arr2) {
						if (arr1[0] > arr2[0]) {
							return 1;
						} else if (arr1[0] < arr2[0]) {
							return -1;
						}
						return 0;
					}
				});
			}
			for (int j = 0; j < spectraNew.get(i).mzIntensity.size(); j++){
				double currMz = spectraNew.get(i).mzIntensity.get(j)[0];
				while (j+1 < spectraNew.get(i).mzIntensity.size() && withinMzTol(currMz, spectraNew.get(i).mzIntensity.get(j+1)[0], massTol)) {
					spectraNew.get(i).mzIntensity.remove(j+1);
				}
			}
		}
	}

	private static void mergeBackground(List<double[]> intensitiesData, List<double[]> mzsData, List<double[]> rtData, boolean[][] read, ArrayList<Spectra> spectraNew) {
		for (int i = 0; i < read.length; i++) {
			ArrayList<double[]> mzInt = new ArrayList<>();
			for (int j = 0; j < read[i].length; j++) {
				if (!read[i][j] && mzsData.get(i)[j] != 0) {
					mzInt.add(new double[]{mzsData.get(i)[j], intensitiesData.get(i)[j]});
				}
			}
			if (mzInt.size() > 0) {
				Spectra spect = new Spectra(rtData.get(i)[1], i, mzInt);
				spectraNew.add(spect);
			}
		}
	}

	// Requires: peak must be found
	static Peak getPeak(List<double[]> rtData, List<double[]> mzData, List<double[]> intensityData,
						int rowStartIndex, int rowEndIndex, double mzMin1, double mzMax1, boolean[][] read, int[][] readded) {
		double rtStart = rtData.get(rowStartIndex)[1];
		double rtEnd = rtData.get(rowEndIndex)[1];
		double mzMin = mzMin1;
		double mzMax = mzMax1;
		ArrayList<ArrayList<Double>> mz = new ArrayList<>();
		ArrayList<ArrayList<Double>> intensity = new ArrayList<>();
		//double[] rt = new double[rowEndIndex - rowStartIndex + 1];
		for (int i = rowStartIndex; i <= rowEndIndex; i++) {
			mz.add(new ArrayList<Double>());
			intensity.add(new ArrayList<Double>());
		}

		// get mz range
		for (int i = rowStartIndex; i <= rowEndIndex; i++) {
			int from = 0;
			for (int j = 0; j < mzData.get(i).length; j++) {
				if (mzData.get(i)[j] >= mzMin1) {
					from = j;
					break;
				}
			}
			int to = from;
			for (int j = from; j < mzData.get(i).length; j++) {
				if (mzData.get(i)[j] > mzMax) {
					to = j - 1;
					break;
				}
			}

			for (int j = from; j <= to; j++) {
				read[i][j] = true;
				mz.get(i - rowStartIndex).add(mzData.get(i)[j]);
				intensity.get(i - rowStartIndex).add(intensityData.get(i)[j]);
				readded[i][j]++;
			}
		}

		ArrayList<Spectra> rtMzIntensity = new ArrayList<>();
		for (int i = 0; i < mz.size(); i++) {
			ArrayList<Double> tempMz = mz.get(i);
			ArrayList<Double> tempIntensity = intensity.get(i);
			ArrayList<double[]> tempList = new ArrayList<>();
			for (int j = 0; j < tempMz.size(); j++) {
				double[] temp = new double[] { tempMz.get(j), tempIntensity.get(j) };
				tempList.add(temp);
			}
			Spectra temp = new Spectra(rtData.get(i+rowStartIndex)[1], i+rowStartIndex, tempList);
			rtMzIntensity.add(temp);
		}

		Peak result = new Peak(rtStart, rowStartIndex, rtEnd, rowEndIndex, mzMin, mzMax, rtMzIntensity);

		return result;

	}

	static int getRtRowIndex(List<double[]> rt, double key) {
		int result = -1;
		int low = 0;
		int high = rt.size() - 1;
		while (low <= high) {
			int mid = (low + high) / 2;
			if (key == rt.get(low)[1]) {
				result = (int) rt.get(low)[0] - 1;
				break;
			} else if (key == rt.get(mid)[1]) {
				result = (int) rt.get(mid)[0] - 1;
				break;
			} else if (key < rt.get(mid)[1])
				high = mid - 1;
			else
				low = mid + 1;
		}
		return result;
	}

	static boolean doubleEqual(double d1, double d2) {
		return abs(d1 - d2) < 1E-14;
	}

	static boolean withinMzTol(double d1, double d2, double tol) {
		return abs(d1 - d2) < tol;
	}

//	static List<List<String>>[] readParametersHelper(String[] fileNames) throws IOException {
//		List<List<String>>[] val = new List[fileNames.length];
//		for (String fileName: fileNames) {
//			List<List<String>> param = readParameters(fileName);
//			val.
//		}
//		return val;
//	}

	static List<List<String>> readParameters(String fileName) throws IOException {
		List<List<String>> list = new ArrayList<>();
		BufferedReader reader = new BufferedReader(new FileReader(fileName));

		try{
			int row = 0;
			String line;
			while((line = reader.readLine()) != null){
				if (row > 7)
					break; // only reads first 7 rows
				line = line.replaceAll("\"", "");
				String[] values = line.split("\t|\\,");
				// this adds the currently parsed line to the 2-dimensional string array
				List<String> genList = new ArrayList<>(Arrays.asList(values));

				if (!genList.isEmpty()) {
					genList.remove(0);
				}
				// this removes the first column

				if (row != 2 && row != 3)
					list.add(genList);
				row ++;
			}

		}catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}

		return list;
	}

	static List<double[]> readInput(String fileName) throws IOException {

		File file= new File(fileName);
		List<double[]> list = new ArrayList<>();
		// this gives you a 2-dimensional array of strings
		Scanner inputStream;

		try{
			inputStream = new Scanner(file);
			inputStream.next();

			while(inputStream.hasNext()){
				String line= inputStream.next();
				String[] values = line.split(",");
				// this adds the currently parsed line to the 2-dimensional string array
				List<String> genList = new ArrayList<>(Arrays.asList(values));

				if (!fileName.contains("_rt.csv")) {
					genList.remove(0);
				}
				double[] rslt = new double[genList.size()];
				rslt[0] = Double.parseDouble(genList.get(0).replaceAll("^\"|\"$", ""));
				for (int i =1; i < genList.size(); i++) {
					if (genList.get(i).equals("NA")) break;
					rslt[i] = Double.parseDouble(genList.get(i).trim());
				}
				list.add(rslt);
			}

			inputStream.close();
		}catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		return list;
	}

	static void writeRT(ArrayList<Spectra> spectraNew, String filename, int num) throws IOException {
		String[] rt = new String[spectraNew.size()];
		FileWriter writer = new FileWriter(filename +"_" + num + "_rt.csv");

		for (int i = 0; i < spectraNew.size(); i++){
			rt[i] = spectraNew.get(i).rt + "";
		}
		String csv = String.join(",", rt);
		writer.write(csv);
		writer.close();
	}

	static void writeMzIntensity(ArrayList<Spectra> spectraNew, String filename, int num) throws IOException {
		double[][] mz = new double[spectraNew.size()][];
		double[][] intensity = new double[spectraNew.size()][];

		for (int i = 0; i < spectraNew.size(); i++){
			mz[i] = new double[spectraNew.get(i).mzIntensity.size()];
			intensity[i] = new double[spectraNew.get(i).mzIntensity.size()];
			for (int j = 0; j < spectraNew.get(i).mzIntensity.size(); j++) {
				mz[i][j] = spectraNew.get(i).mzIntensity.get(j)[0];
				intensity[i][j] = spectraNew.get(i).mzIntensity.get(j)[1];
			}
		}

		List<String> mzList = Arrays.stream(mz)
				.map(line -> Arrays.stream(line)
						.mapToObj(String::valueOf)
						.collect(Collectors.joining(",")))
				.collect(Collectors.toList());
		try {
			Files.write(Paths.get(System.getProperty("user.dir"),filename +"_" + num + "_mz.csv"), mzList);
		} catch (IOException e) {
			e.printStackTrace();
		}
		List<String> intensityList = Arrays.stream(intensity)
				.map(line -> Arrays.stream(line)
						.mapToObj(String::valueOf)
						.collect(Collectors.joining(",")))
				.collect(Collectors.toList());
		try {
			Files.write(Paths.get(System.getProperty("user.dir"),filename +"_" + num + "_intensity.csv"), intensityList); // StandardOpenOption.CREATE
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
