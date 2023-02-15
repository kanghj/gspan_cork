package smu.hongjin;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.mskcc.cbio.portal.stats.FisherExact;

import io.github.tonyzzx.gspan.gSpan;
import io.github.tonyzzx.gspan.gSpan.GRAPH_LABEL;

public class CountingUtils {

//	public static double initialFeatureScore(int A_S0, int A_S1, int B_S0, int B_S1, int U_S0, int U_S1, double AWeight,
//			double BWeight, double UWeight) {
//		// first component: from original CORK paper
//		// ("Near-optimal supervised feature selection among frequent subgraphs" by
//		// Thoma et al.)
//		// SIAM International Conference on Data Mining. 2009
//		// SIAM is a rank A conference.
//		double correspondanceBetweenLabels = -1 * (AWeight * A_S0 * BWeight * B_S0 + AWeight * A_S1 * BWeight * B_S1);
//
//		// second component: incentize skewedness. The more the unlabeled data shifts to
//		// the majority class, the better
//		// i.e. incentize correspondense between U and A
//		double skewedness = UWeight * U_S1 * AWeight * A_S1;
//		
//		// some U and many of B
//		skewedness += UWeight * U_S0 * BWeight * B_S1;
//		
//		// third component: penalty for unlabeled data that is still unlabeled;
////		double lackOfLabels = -1 * UWeight * U_S0 * UWeight * U_S0;		// 
//
//		return correspondanceBetweenLabels; // + skewedness;
//		
//		// a feature is good w.r.t. to the misuse distribution if
//		// 		1. unlabeled distribution is skewed
//		// 		2.  
//	}

	public static double initialFeatureScore(int A_S0, int A_S1, int B_S0, int B_S1, int U_S0, int U_S1, int A_N,
			int B_N, double AWeight, double BWeight,  long ID) {

//		System.out.println("\t==getting score==");
//		System.out.print("\tfirst component=" + Math.abs(AWeight * A_S1 - BWeight * B_S1));
		
//		if ()
		if (Double.isInfinite(BWeight) || Double.isInfinite(AWeight)) {
			
			
			return Math.max(A_S1, B_S1);
		}

		double expectedRatio = AWeight / (AWeight + BWeight);
		LoggingUtils.logOnce("\t\t\t, expectedRatio=" + expectedRatio);

		double score = Math.abs(AWeight * A_S1 - BWeight * B_S1) // difference bet. percentages [0..100]
//				- Math.abs(gSpan.skewnessImportance * expectedRatio - UWeight * U_S1) // [0..skewnessImportance/2]
		;
//		if (Math.abs(AWeight * A_S1 - BWeight * B_S1) > 0 && score <= 0) {
//			gSpan.wouldNotBePrunedWithoutSemiSupervisedFilters += 1;
//			System.out.println("\t\t" + ID + " would be accepted if not for the distribution penalty!");
//			System.out.println("\t\tA_S1=" + A_S1 + ",B_S1=" + B_S1);
//		}

//		long[][] currentCounts = { { A_S1, A_S0 }, { B_S1, B_S0 } };
//		ChiSquareTest test = new ChiSquareTest();
//		
//		System.out.println("\t==debug==");
//		System.out.println("\t\tA_S1=" + A_S1 + ",B_S1=" + B_S1);
//		System.out.println("\t\tA_S0=" + A_S0 + ",B_S0=" + B_S0);
//		
//
//		double score = test.chiSquare(currentCounts);
//		System.out.println("\t\tvalue=" + score);
		return score;

	}


	

//	public static double upperBound(double q_s, int A_S0, int A_S1, int B_S0, int B_S1, int U_S0, int U_S1,  double AWeight,
//			double BWeight, double UWeight) {
//		// first component upper bound: from original CORK paper
//		// "Near-optimal supervised feature selection among frequent subgraphs" by Thoma
//		// et al.
//		double maxCorrespondanceIncrease = Math.max(
//				Math.max(
//						AWeight * A_S1 * (BWeight * (B_S1 - B_S0)), BWeight * B_S1 * (AWeight * (A_S1 - A_S0))), 0);
//
//		// second component upper bound: 
//		// both U_S1 and A_S1 cannot increase
//		// but the penalty for U and B can decrease.
//		// U_S0 can increase, although B_S1 cannot. Thus, the best case is that all of U_S1 becomes U_S0
//		double maxSkewIncrease = 0 ;// UWeight * (U_S1) * BWeight * B_S1;
//				
//
//		return q_s + maxCorrespondanceIncrease + maxSkewIncrease;
//	}

	public enum UpperBoundReturnType {
		BAD, GOOD;
	}

	public static UpperBoundReturnType upperBound(double current, int A_S0, int A_S1, int B_S0, int B_S1, int U_S0,
			int U_S1, double AWeight, double BWeight, boolean isDebug) {

		if (A_S1 == 0 && B_S1 == 0) {
			return UpperBoundReturnType.BAD;
		}
		
		
		ChiSquareTest test = new ChiSquareTest();

		long[][] currentCounts = { { A_S1, A_S0 }, { B_S1, B_S0 } };
		double currentPValue = test.chiSquareTest(currentCounts);

		long[][] bestCounts1 = { { A_S1, A_S0 }, { 0, B_S0 + B_S1 } };
		double bestPValue1 = test.chiSquareTest(bestCounts1);
		if (Double.isNaN(bestPValue1)) { // can reach Nan if aA_S1 was already 0
			bestPValue1 = currentPValue;
		}

		long[][] bestCounts2 = { { 0, A_S0 + A_S1 }, { B_S1, B_S0 } };
		double bestPValue2 = test.chiSquareTest(bestCounts2);
		if (Double.isNaN(bestPValue2)) {// can reach Nan if B_S1 was already 0
			bestPValue2 = currentPValue;
		}

		if (bestPValue1 > 0.1 && bestPValue2 > 0.1) { // best case still bad

			return UpperBoundReturnType.BAD;
			
		} else if (fuzzyEquals(Math.min(bestPValue1, bestPValue2), currentPValue, 0.0005)) {


			return UpperBoundReturnType.BAD;
			
		} else {
		
			return UpperBoundReturnType.GOOD;
		

		}
	}

	public static boolean fuzzyEquals(double a, double b, double tolerance) {
		return Math.copySign(a - b, 1.0) <= tolerance
				// copySign(x, 1.0) is a branch-free version of abs(x), but with different NaN
				// semantics
				|| (a == b) // needed to ensure that infinities equal themselves
				|| (Double.isNaN(a) && Double.isNaN(b));
	}

	public static int minimumCountForSignificanceMinority(int classACounts, int classBCounts) {
		ChiSquareTest test = new ChiSquareTest();

		int majorityCount = classACounts > classBCounts ? classACounts : classBCounts;
		int minorityCount = classACounts > classBCounts ? classBCounts : classACounts;

		for (int i = 3; i < minorityCount; i++) {
			long[][] currentCounts = { { 0, majorityCount }, { i, minorityCount - i } };
			double currentPValue = test.chiSquareTest(currentCounts);
			if (currentPValue < 0.05) {
				LoggingUtils.logOnce("MINIMUM COUNT FOR signifance in the minority class=" + i);
				return i;
			}
		}

		System.out.println("WARN: features of minority class never reach significance. Please double-check the counts.");
		return 3;
	}

	// actually we don't have to do this for all subgraphs.
	public static Map<Long, Double> findClosestLabelledPointForKUnLabelled2(int k, Set<Long> subgraphIds,
			Map<Long, Set<Long>> misuseSubgraphCoverage, Map<Long, Set<Long>> correctUseSubgraphCoverage,
			Map<Long, Set<Long>> unlabeledGraphsCoverage, Set<Long> allMisuses, Set<Long> allCorrect,
			Set<Long> allUnlabeled) {

		Map<Long, Map<Long, Double>> unlabeledToLabeledDistance = new HashMap<>();
		// OR:
		// for all unlabeled points
		for (long unlabeled : allUnlabeled) {
			if (!unlabeledToLabeledDistance.containsKey(unlabeled)) {
				unlabeledToLabeledDistance.put(unlabeled, new HashMap<>());
			}
			Map<Long, Double> map = unlabeledToLabeledDistance.get(unlabeled);

			iterateLabeledAndCountDistance(subgraphIds, misuseSubgraphCoverage, unlabeledGraphsCoverage, allMisuses,
					unlabeled, map);
			iterateLabeledAndCountDistance(subgraphIds, correctUseSubgraphCoverage, unlabeledGraphsCoverage, allCorrect,
					unlabeled, map);
		}
		// now,
		// unlabeledGraphToLabeledGraphDistance contains a mapping of all unlabeled ->
		// labeled -> value proportional to distance

		List<Map.Entry<Long, Double>> shortest = new ArrayList<>();
		for (Entry<Long, Map<Long, Double>> entry : unlabeledToLabeledDistance.entrySet()) {
			Long unlabelledId = entry.getKey();

			double minDistance = 9999; // for this unlabeled point, the smallest distance to a labeled point
			for (Entry<Long, Double> value : entry.getValue().entrySet()) {
				if (value.getValue() < minDistance) {
					minDistance = value.getValue();
				}
			}

			if (minDistance > shortest.get(k - 1).getValue())
				continue;

			// insert into sorted list
			for (int i = 0; i < k; i++) {
				if (minDistance <= shortest.get(i).getValue()) {
					shortest.add(i, new AbstractMap.SimpleEntry<Long, Double>(unlabelledId, minDistance));

					shortest.remove(shortest.size() - 1);
					break;
				}
			}
		}

		assert shortest.size() == k;
		Map<Long, Double> result = new HashMap<>();
		for (Entry<Long, Double> entry : shortest) {
			result.put(entry.getKey(), Math.sqrt(entry.getValue()));
		}

		return result;
	}

	private static void iterateLabeledAndCountDistance(Set<Long> subgraphIds // the features
			, Map<Long, Set<Long>> labeledSubgraphCoverage // subgraphId -> which graph IDs contain the subgraph
			, Map<Long, Set<Long>> unlabeledGraphsCoverage // subgraphId -> which graph IDs contain the subgraph
			, Set<Long> labeledGraphIds // labeled graph IDS
			, long unlabeled, Map<Long, Double> map) {

		for (long labeled : labeledGraphIds) {
			if (!map.containsKey(labeled))
				map.put(labeled, 0.0);

			for (long subgraphId : subgraphIds) {
				boolean isSubgraphInUnlabelled = unlabeledGraphsCoverage.get(subgraphId).contains(unlabeled);
				boolean isSubgraphInLabelled = labeledSubgraphCoverage.get(subgraphId).contains(labeled);
				if (isSubgraphInUnlabelled != isSubgraphInLabelled) {
					// increase distance,
					map.put(labeled, map.get(labeled) + 1.0);
				}
			}
		}
	}

	public static Map<Long, Double> findClosestLabelledPointForKUnLabelled(int k, Set<Long> subgraphIds,
			Map<Long, Set<Long>> misuseSubgraphCoverage, Map<Long, Set<Long>> correctUseSubgraphCoverage,
			Map<Long, Set<Long>> unlabeledGraphsCoverage) {

		Map<Long, Map<Long, Double>> unlabeledGraphToLabeledGraphDistance = new HashMap<>();
		for (long subgraphId : subgraphIds) {
			// iterate all unlabeledGraphs
			for (long unlabeledGraph : unlabeledGraphsCoverage.get(subgraphId)) {
				// in each unlabeled graph
				// find distance to all misuseGraphs and correctUsageGraphs

				if (!unlabeledGraphToLabeledGraphDistance.containsKey(unlabeledGraph)) {
					unlabeledGraphToLabeledGraphDistance.put(unlabeledGraph, new HashMap<>());
				}
				for (long labeledGraph : misuseSubgraphCoverage.get(subgraphId)) {
					Map<Long, Double> map = unlabeledGraphToLabeledGraphDistance.get(unlabeledGraph);
					if (!map.containsKey(labeledGraph)) {
						map.put(labeledGraph, 0.0);
					}
					map.put(labeledGraph, map.get(labeledGraph) + 1.0); // TODO not correct
				}

				for (long labeledGraph : correctUseSubgraphCoverage.get(subgraphId)) {
					Map<Long, Double> map = unlabeledGraphToLabeledGraphDistance.get(unlabeledGraph);
					if (!map.containsKey(labeledGraph)) {
						map.put(labeledGraph, 0.0);
					}
					map.put(labeledGraph, map.get(labeledGraph) + 1.0); // TODO not correct
				}
			}
		}

		List<Map.Entry<Long, Double>> shortest = new ArrayList<>();
		for (Entry<Long, Map<Long, Double>> entry : unlabeledGraphToLabeledGraphDistance.entrySet()) {
			Long unlabelledId = entry.getKey();
			double minDistance = 9999;
			for (Entry<Long, Double> value : entry.getValue().entrySet()) {
				if (value.getValue() < minDistance) {
					minDistance = value.getValue();
				}
			}

			if (minDistance > shortest.get(k - 1).getValue())
				continue;

			// insert into sorted list
			for (int i = 0; i < k; i++) {
				if (minDistance <= shortest.get(i).getValue()) {
					shortest.add(i, new AbstractMap.SimpleEntry<Long, Double>(unlabelledId, minDistance));

					shortest.remove(shortest.size() - 1);
					break;
				}
			}
		}

		assert shortest.size() == k;
		Map<Long, Double> result = new HashMap<>();
		for (Entry<Long, Double> entry : shortest) {
			result.put(entry.getKey(), entry.getValue());
		}

		return result;
	}

	public static void writeUnlabelledGraphFeatures(gSpan gSpan, Map<Long, Set<Integer>> labeledCoverage, 
			Map<Long, Set<Integer>> unlabeledCoverage,
			BufferedWriter writer) throws IOException {
		System.out.println("\tConsolidating and writing graph and their subgraph features");

		if (labeledCoverage.size() != gSpan.selectedSubgraphFeatures.size()) {
			
			
			throw new RuntimeException("wrong size!");
		}

		List<Integer> graphs = new ArrayList<>();
		for (Entry<Long, Set<Integer>> entry : unlabeledCoverage.entrySet()) {
			graphs.addAll(entry.getValue());
		}

		List<Long> features = labeledCoverage.keySet().stream().sorted().collect(Collectors.toList());
		writer.write("graph_id,is_correct");
		for (Long feature : features) {
			writer.write(",feature_" + feature);
		}
		writer.write("\n");

		System.out.println("\tWriting to unlabelled feature vector file");

		// <graph id>, is_correct, feature_1, feature_2, feature_3, ... \n
		for (Integer graph : graphs) {

			for (int graphNum = 0; graphNum < 1; graphNum++) {
				int timesToRepeat = 1;
				for (int repeated = 0; repeated < timesToRepeat; repeated++) {

					writer.write(graph + "_" + graphNum + "_" + repeated + ",");

					writer.write("?");

					for (Long feature : features) {
						if (unlabeledCoverage.get(feature).contains(graph)) {
							writer.write(",1");
						} else {
							writer.write(",0");
						}
					}
					writer.write("\n");
				}
			}
		}

	}

	public static void writeGraphFeatures(gSpan gSpan, Map<Long, Set<Integer>> coverage, Map<Long, Double> sortedSubgraphFeatures, BufferedWriter writer)
			throws IOException {
		System.out.println("\tConsolidating and writing graph and their subgraph features");

		if (coverage.size() != gSpan.selectedSubgraphFeatures.size()) {
			throw new RuntimeException("wrong size!");
		}

		List<Integer> graphs = new ArrayList<>();

		graphs.addAll(gSpan.correctUses);
		graphs.addAll(gSpan.misuses);
		
		System.out.println("correct sz " + gSpan.correctUses.size());
		System.out.println("misuses sz " + gSpan.misuses.size());
		
		System.out.println("graphs sz " + graphs.size()); // graphs are different!!

//		List<Long> features = coverage.keySet().stream().sorted().collect(Collectors.toList());
		
		List<Long> features = new ArrayList<>(sortedSubgraphFeatures.keySet());
		
		writer.write("graph_id,is_correct");
		for (Long feature : features) {
			if (!coverage.containsKey(feature)) {
				continue;
			}
			writer.write(",feature_" + feature);
		}
		writer.write("\n");

		// amplify minority class
		int numCorrect = 0;
		int numMisuses = 0;
		for (Integer graph : graphs) {
			if (gSpan.misuses.contains(graph)) {
				numMisuses += gSpan.quantities.get(graph);
			} else if (gSpan.correctUses.contains(graph)) {
				numCorrect += gSpan.quantities.get(graph);
			}
		}
		
		int repeatMisuses;
		int repeatCorrect;
		if (numCorrect > numMisuses) {
			repeatCorrect = 1;
			repeatMisuses = Math.toIntExact(Math.round(Math.ceil((float) numCorrect / Math.max(1, numMisuses) )));
//			repeatMisuses = 1;
		} else {
			repeatMisuses = 1;
			repeatCorrect = Math.toIntExact(Math.round(Math.floor((float) numMisuses / Math.max(1, numCorrect))));
//			repeatCorrect = 1;
		}

		System.out.println("\tWriting to feature vector file");
		System.out.println("\tRepeating; repeatCorrect= " + repeatCorrect + " , and repeatMisuses=" + repeatMisuses);

		int debugUniqGraphCountMisuses = 0;
		int debugUniqGraphCountCorrects = 0;
		int debugGraphCountMisuses = 0;
		int debugGraphCountCorrects = 0;
		
		// ok, so we know misuses and correctUses remain constant.
		System.out.println("gSpan.misuses sz" + gSpan.misuses.size());
		System.out.println("gSpan.correctUses sz" + gSpan.correctUses.size());
		
		// <graph id>, is_correct, feature_1, feature_2, feature_3, ... \n
		for (Integer graph : graphs) {
			if (gSpan.misuses.contains(graph)) {
				debugUniqGraphCountMisuses += 1;
				debugGraphCountMisuses += gSpan.quantities.get(graph);
			} else if (gSpan.correctUses.contains(graph)) {
				debugUniqGraphCountCorrects += 1;
				debugGraphCountCorrects += gSpan.quantities.get(graph);
			}

			for (int graphNum = 0; graphNum < gSpan.quantities.get(graph); graphNum++) {
				boolean isCorrect = gSpan.correctUses.contains(graph);
				boolean isMisuse = gSpan.misuses.contains(graph);

				int timesToRepeat = isCorrect ? repeatCorrect : repeatMisuses;
				for (int repeated = 0; repeated < timesToRepeat; repeated++) {

					writer.write(graph + "_" + graphNum + "_" + repeated + ",");

					if (isCorrect) {
						writer.write("1");
					} else if (isMisuse) {
						writer.write("0");
					} else {
						throw new RuntimeException("graph label is incorrect somehow " + graph);
					}

					for (Long feature : features) {
						if (!coverage.containsKey(feature)) {
							continue;
						}
						if (coverage.get(feature).contains(graph)) {
							writer.write(",1");
						} else {
							writer.write(",0");
						}
					}
					writer.write("\n");
				}
			}
		}

		System.out.println("debugUniqGraphCountMisuses : " + debugUniqGraphCountMisuses);
		System.out.println("debugUniqGraphCountCorrects : " + debugUniqGraphCountCorrects);
		System.out.println("debugGraphCountMisuses : " + debugGraphCountMisuses);
		System.out.println("debugGraphCountCorrects : " + debugGraphCountCorrects);
		System.out.println("\tCompleted consolidation and writing");
	}

}
