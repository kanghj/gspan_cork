package smu.hongjin;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class ExcludeBoringNodes {

	public static String pathToVertMap;
	
	public static List<Integer> IdsOfBoringNodes() {
		List<String> lines;
		try {
			lines = Files.readAllLines(Paths.get(pathToVertMap));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException(e);
		}
		List<Integer> result = new ArrayList<>();
		
		for (String line : lines) {
			String[] splitted = line.split(",");
			if (splitted[0].equals("Object")) {
				result.add(Integer.parseInt(splitted[splitted.length - 1]));
			}
			if (splitted[0].equals("String")) {
				result.add(Integer.parseInt(splitted[splitted.length - 1]));
			}
			if (splitted[0].equals("Integer")) {
				result.add(Integer.parseInt(splitted[splitted.length - 1]));
			}
			if (splitted[0].equals("int")) {
				result.add(Integer.parseInt(splitted[splitted.length - 1]));
			}
			
		}
		
		System.out.println("Boring nodes: ");
		System.out.println(result);
		return result;
	}
	
}
