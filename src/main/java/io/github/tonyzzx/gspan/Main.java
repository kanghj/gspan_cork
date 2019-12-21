package io.github.tonyzzx.gspan;

import org.apache.commons.cli.*;

import io.github.tonyzzx.gspan.model.Graph;
import smu.hongjin.CountingUtils;
import smu.hongjin.ExcludeBoringNodes;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import java.util.Scanner;
import java.util.Set;

public class Main {

	public static void main(String[] args) throws IOException {
		Arguments arguments = Arguments.getInstance(args);
		System.out.println("args: " + Arrays.toString(args));

		File inFile = new File(arguments.inFilePath);
		
		//
		ExcludeBoringNodes.pathToVertMap = arguments.boringNodePath;
		Graph.boringNodes = ExcludeBoringNodes.IdsOfBoringNodes();
		//
		
		
		File outFile = new File(arguments.outFilePath);
		try (FileReader reader = new FileReader(inFile)) {
			try (FileWriter writer = new FileWriter(outFile)) {
				gSpan gSpan = new gSpan();
				System.out.println("gSpan is mining...");
				gSpan.run(reader, writer, arguments.minSup, arguments.maxNodeNum, arguments.minNodeNum);

				System.out.println("It's done! The result is in  " + arguments.outFilePath + " .");

				postprocess(arguments, gSpan);
			}
		}

	}

	private static void postprocess(Arguments arguments, gSpan gSpan) throws IOException {
		
		List<Integer> graphs = new ArrayList<>();
		graphs.addAll(gSpan.correctUses);
		graphs.addAll(gSpan.misuses);
		
		Map<Integer, Set<Long>> graphsToFeaturesContained = new HashMap<>();
		
		for (int graph : graphs) {
			graphsToFeaturesContained.putIfAbsent(graph, new HashSet<>());
		}

		Map<Long, Double> sortedSubgraphFeatures = gSpan.selectedSubgraphFeatures.entrySet().stream()
				.sorted(Entry.comparingByValue(Comparator.reverseOrder()))
//				.sorted(Entry.comparingByValue())
				.collect(Collectors.toMap(Entry::getKey, Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new));


		Set<Long> weakSubgraphFeatures = new HashSet<>();
		
		// the worst CORK value is when all the misuses and correct uses cannot be disambiguated.
		int previousCorkScore = gSpan.totalCorrectUses * gSpan.totalMisuses;// Integer.MAX_VALUE;
		System.out.println("worst score : " + previousCorkScore);
		for (Map.Entry<Long, Double> entry : sortedSubgraphFeatures.entrySet()) {
			Long featureId = entry.getKey();

			Set<Integer> subgraphCovers = new HashSet<>(gSpan.coverage.get(featureId));
			
			for (Integer graph : subgraphCovers) {
				graphsToFeaturesContained.get(graph).add(featureId);
			}
			
			// see if the number of graphs that cannot be disambiguated has decreased?
			// if so, the subgraph is not weak
			// otherwise, drop this subgraph
			
			int corkScore = 0;
			
			
			for (int i = 0; i < graphs.size(); i++) {
				int graph1 = graphs.get(i);
				Set<Long> features1 = graphsToFeaturesContained.get(graph1);
				boolean isMisuse1 = gSpan.misuses.contains(graph1);

				for (int j = i + 1; j < graphs.size(); j++) {
					int graph2 = graphs.get(j);
					
					boolean isMisuse2 = gSpan.misuses.contains(graph2);
					
					if (isMisuse1 == isMisuse2) {
						continue;
					}
					
					Set<Long> intersection = new HashSet<>(graphsToFeaturesContained.get(graph2));
					intersection.retainAll(features1);
					
					if (intersection.size() == features1.size()) { // cannot disambiguate
						corkScore += 1;
					}
				}
			}
			
			boolean isWeak = true;
			if (corkScore < previousCorkScore - 1) {
				isWeak = false;
			}
			previousCorkScore = corkScore;

			System.out.println("feature " + featureId + " - CORK score now=" + corkScore);
			if (isWeak) {
				System.out.println("\t" + featureId + " is getting dropped");
				weakSubgraphFeatures.add(featureId);
				
				for (Integer graph : subgraphCovers) {
					graphsToFeaturesContained.get(graph).remove(featureId); // shouldn't matter, but let's remove anyway
				}
			}
		}

		// remove the lame features
		for (long weakFeature : weakSubgraphFeatures) {
			gSpan.selectedSubgraphFeatures.remove(weakFeature);
			gSpan.coverage.remove(weakFeature);
		}

		try (BufferedWriter selectedSubGraphWriter = new BufferedWriter(
				new FileWriter(arguments.outFilePath + "_best_subgraphs.txt"))) {
			for (Map.Entry<Long, Double> subgraphFeature : sortedSubgraphFeatures.entrySet()) {
				if (!gSpan.coverage.containsKey(subgraphFeature.getKey())) {
					if (!weakSubgraphFeatures.contains(subgraphFeature.getKey())) {
						throw new RuntimeException("huhhh");
					}
					continue; // we have removed this weak feature
				}
				selectedSubGraphWriter.write(subgraphFeature.getKey() + "," + subgraphFeature.getValue());
				selectedSubGraphWriter.write("\n");
			}
		}

		// print subgraph features
		System.out.println(
				"The identified discriminative subgraphs are in  " + arguments.outFilePath + "_best_subgraphs.txt");
		try (BufferedWriter featuresWriter = new BufferedWriter(
				new FileWriter(arguments.outFilePath + "_features.txt"))) {
			CountingUtils.writeGraphFeatures(gSpan, gSpan.coverage, sortedSubgraphFeatures, featuresWriter);
		}
		System.out.println("The feature vectors of labeled graphs are in " + arguments.outFilePath + "_features.txt");

//		try (BufferedWriter featuresWriter = new BufferedWriter(
//				new FileWriter(arguments.outFilePath + "_unlabelled_features.txt"))) {
//			CountingUtils.writeUnlabelledGraphFeatures(gSpan, gSpan.coverage, gSpan.unlabeledCoverage, featuresWriter);
//		}
//		System.out.println(
//				"The feature vectors of unlabeled graphs are in " + arguments.outFilePath + "_unlabelled_features.txt");

		// find new examples to label
		System.out.println("Computing which unlabeled graphs were not covered");
		System.out.println("\tand which labeled graphs were not covered");

		for (Long feature : gSpan.selectedSubgraphFeatures.keySet()) {
			Set<Integer> coveredGraphs = gSpan.unlabeledCoverage.get(feature);
			gSpan.uncoveredUnlabeledGraphs.removeAll(coveredGraphs);

		}	
	
		writeClingo(arguments, gSpan);

		System.out.println("The IDs of the uncovered methods have been written to " + arguments.outFilePath
				+ "_interesting_unlabeled.txt");
		System.out.println("The prefixes of the files we care about are " + arguments.outFilePath);
	}

	private static void writeClingo(Arguments arguments, gSpan gSpan) throws IOException {
		System.out.println("writing clingo to " + arguments.outFilePath + "_interesting_subgraphs_coverage.txt");
		try (BufferedWriter selectedSubGraphWriter = new BufferedWriter(
				new FileWriter(arguments.outFilePath + "_interesting_subgraphs_coverage.txt"))) {
			
			selectedSubGraphWriter.write("#const n=" + gSpan.minimumToLabel + ".\n");
			selectedSubGraphWriter.write("#const s=" + (gSpan.totalCorrectUses + gSpan.totalMisuses + gSpan.totalUnlabeled) / 200 + ".\n");
			
			Set<Integer> graphs = new HashSet<>();
			
			
			for (Set<Integer> graphs1 : gSpan.unlabeledCoverage.values()) {
				graphs.addAll(graphs1);
			}
			graphs.retainAll(gSpan.uncoveredUnlabeledGraphs);
			
			for (int graph : graphs) {
				selectedSubGraphWriter.write("graph(" + graph + ")" + ".\n");
			}
		
			
			boolean pruneMore = gSpan.interestingSubgraphs.size() > 50000;
			
			Map<Long, Long> countsOfSubgraphs = new HashMap<>();
			for (long interestingSubgraph : gSpan.interestingSubgraphs ) {
				
				long size = gSpan.unlabeledCoverage.get(interestingSubgraph)
						.stream()
						.filter(graph -> graphs.contains(graph))
						.count();
				if (size < gSpan.minimumToLabel) {
					continue;
				}
				// if too many subgraphs, we prune more subgraphs..
				if (pruneMore && size < gSpan.minimumToLabel * 3) {
					continue;
				}
				countsOfSubgraphs.put(interestingSubgraph, size);			
			}
			for (long subgraph: gSpan.unlabeledCoverage.keySet()) {
				if (!countsOfSubgraphs.containsKey(subgraph)) {
					continue;
				}
				selectedSubGraphWriter.write("subgraph(" + subgraph + ")" + ".\n");
			}
			List<Long> top = countsOfSubgraphs.entrySet().stream()
				.sorted(Entry.comparingByValue(Comparator.reverseOrder()))
//				.limit(countsOfSubgraphs.entrySet().size() / 2) 	// 
				.map(entry -> entry.getKey())
				.collect(Collectors.toList());
			
			
			for (long interestingSubgraph : top) {
				for (int graph : gSpan.unlabeledCoverage.get(interestingSubgraph)) {
					if (graphs.contains(graph )) {
						selectedSubGraphWriter.write("covers(" + graph + "," + interestingSubgraph + ")" + ".\n");
					}
				}
				
			}
			
			selectedSubGraphWriter.write("{ selected(G) : graph(G) } <= 1 :- graph(G).\n");


			selectedSubGraphWriter.write("sufficiently_covered(SG) :- #count { G : selected(G), covers(G, SG), graph(G) } >= n, subgraph(SG).\n");

			selectedSubGraphWriter.write(":- { selected(G)} > s .\n");

			selectedSubGraphWriter.write("covered_count(X) :- X = #count {SG : sufficiently_covered(SG), subgraph(SG)}.\n");
//			selectedSubGraphWriter.write("num_selected(Y) :- Y =  #count{G : selected(G)}.\n");

//			selectedSubGraphWriter.write("#minimize {Y@1 : num_selected(Y)}.\n");
			selectedSubGraphWriter.write("#maximize {X : covered_count(X) }.\n");

			selectedSubGraphWriter.write("#show selected/1.\n");
			selectedSubGraphWriter.write("#show covered_count/1.\n");
//			selectedSubGraphWriter.write("#show num_selected/1.\n");
			
		}
	}

	private static class Arguments {
		private static Arguments arguments;

		private String[] args;

		String inFilePath;
		
		String boringNodePath;
		
		long minSup;
		long minNodeNum = 0;
		long maxNodeNum = Long.MAX_VALUE;
		String outFilePath;

		private Arguments(String[] args) {
			this.args = args;
		}

		static Arguments getInstance(String[] args) {
			arguments = new Arguments(args);
			if (args.length > 0) {
				arguments.initFromCmd();
			} else {
				arguments.initFromRun();
			}
			return arguments;
		}

		/***
		 * User inputs args.
		 */
		private void initFromCmd() {
			Options options = new Options();
			options.addRequiredOption("d", "data", true, "(Required) File path of data set");
			options.addRequiredOption("s", "sup", true, "(Required) Minimum support");
			options.addOption("b", "boring", true, "Boring Nodes path");
			options.addOption("i", "min-node", true, "Minimum number of nodes for each sub-graph");
			options.addOption("a", "max-node", true, "Maximum number of nodes for each sub-graph");
			options.addOption("r", "result", true, "File path of result");
			options.addOption("h", "help", false, "Help");

			CommandLineParser parser = new DefaultParser();
			HelpFormatter formatter = new HelpFormatter();
			CommandLine cmd = null;
			try {
				cmd = parser.parse(options, args);
				if (cmd.hasOption("h")) {
					formatter.printHelp("gSpan", options);
					System.exit(0);
				}
			} catch (ParseException e) {
				formatter.printHelp("gSpan", options);
				System.exit(1);
			}

			inFilePath = cmd.getOptionValue("d");
			boringNodePath = cmd.getOptionValue("b");
			minSup = Long.parseLong(cmd.getOptionValue("s"));
			minNodeNum = Long.parseLong(cmd.getOptionValue("i", "0"));
			maxNodeNum = Long.parseLong(cmd.getOptionValue("a", String.valueOf(Long.MAX_VALUE)));
			outFilePath = cmd.getOptionValue("r", inFilePath.replace(".txt", "") + "_result");
		}

		/***
		 * User runs it directly.
		 */
		private void initFromRun() {
			try (Scanner sc = new Scanner(System.in)) {
				System.out.println("Please input the file path of data set: ");
				inFilePath = sc.nextLine();
				System.out.println("Please set the minimum support: ");
				minSup = sc.nextLong();
				outFilePath = inFilePath + "_result";

				maxNodeNum = 6;
			}
		}
	}
}
