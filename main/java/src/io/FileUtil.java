package io;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;

public class FileUtil {
	/**
	 * Writes the given lines into a file with the given filename
	 * 
	 * @throws IOException
	 */
	public static void writeFile(String fileName, String content) throws IOException {
		File f = new File(fileName);
		FileWriter fw = new FileWriter(f);
		fw.append(content);
		fw.close();
	}

	/**
	 * Writes the given lines into a file with the given filename
	 * 
	 * @throws IOException
	 */
	public static void writeFile(String fileName, String[] lines, boolean append) throws IOException {
		File f = new File(fileName);
		FileWriter fw = new FileWriter(f, append);
		for (String line : lines) {
			fw.append(line + "\n");
		}
		fw.close();
	}

	/**
	 * Stores the given properties in a file with the given filename
	 * 
	 * @throws IOException
	 */
	public static void saveProperties(String fileName, Properties p) throws IOException {
		File f = new File(fileName);
		FileWriter fw = new FileWriter(f);
		p.store(fw, "");
		fw.close();
	}
}
