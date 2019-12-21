package smu.hongjin;

import java.time.Duration;
import java.time.Instant;
import java.time.temporal.Temporal;
import java.util.HashSet;
import java.util.Set;

public class LoggingUtils {

	public static Set<String> logged = new HashSet<>();
	private static Temporal start;

	public static void logOnce(String log) {

		if (!logged.contains(log)) {
			System.out.println(log);
			logged.add(log);
		}
	}

	public static void logTimingStatistics() {
		if (start == null) {
			start = Instant.now();
			System.out.println("Start tracking time at " + start);
			return;
		}

		Instant currentTime = Instant.now();
		long timeElapsed = Duration.between(start, currentTime).toMillis();
		long minutes = (timeElapsed / 1000) / 60;
		long seconds = (timeElapsed / 1000) % 60;
		long ms = (timeElapsed % 1000);
		System.out.println("\tTotal elapsed time: " + minutes + " minutes " + seconds + " seconds " + ms + "ms");
	}

}
