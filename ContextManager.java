package scyther;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.RandomAccessFile;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;

/**
 * Manager class for building a term-context matrix. Implements algorithms for calculating cosine similarities and positive point wise mutual information.
 * @author Jacob Burton
 *
 */
public class ContextManager {
	
	// Constants
	public static final int WINDOW_SIZE = 4;
	
	// Class variables
	private float[][] TermContextMatrix_PPMI;
	private float[][] TermContextMatrix_FREQ;
	
	private int V;
	private float matrixSum = 0f;
	private float[] rowSums;
	private float[] colSums;
	
	private HashMap<String, Integer> wordIndex;
	
	/**
	 * Constructs a term-context matrix from a corpus file
	 * This method has a time complexity of O(D^2).
	 * This method has a space complexity of O(V^2).
	 * D represents the number of documents in the corpus and V represents the number of unique tokens across all documents. 
	 * @param D A file object that points to a corpus
	 */
	public void BuildTermContextMatrix(File D, boolean writeToFile) {
		
		System.out.println("Building term-context matrix...");
		System.out.println("Memory used: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / Math.pow(1000, 3) + " GB");
		
		// If necessary, filter corpus to english words only
		File filteredCorpus = new File("en-context.txt");
		if(!filteredCorpus.isFile()) FilterCorpus(D);
		
		// Build word index
		wordIndex = new HashMap<>();
		
		// Initialize variables
		float[][] C = null;
		int matrixSize = 0;
		
		String[] window = null;
		Set<String> uniqueWindow;
		
		int termIndex = 0;
		int contextIndex = 0;

		try {
		
			// Read corpus file
			FileReader fr1 = new FileReader(filteredCorpus);
			BufferedReader br1 = new BufferedReader(fr1);
		
			// Populate word index
			String currentToken1;
			int index = 0;
			while((currentToken1 = br1.readLine()) != null) {
				if(!wordIndex.containsKey(currentToken1)) {
					wordIndex.put(currentToken1, index++);
				}
			}
			br1.close();
			
			System.out.println("\nInitializing matrix...");
			
			matrixSize = wordIndex.size();
			System.out.println("V=" + matrixSize);
			C = new float[matrixSize][matrixSize];
			window = new String[(WINDOW_SIZE * 2) + 1];
			
			System.out.println("Matrix initialized.");
			System.out.println("Memory used: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / Math.pow(1000, 3) + " GB");
			
			FileReader fr2 = new FileReader(filteredCorpus);
			BufferedReader br2 = new BufferedReader(fr2);
			
			// Read from corpus, get context words, and update term-context matrix accordingly
			System.out.println("\nUpdating term-context matrix...");
			
			String currentToken2;
			String term;
			boolean doneFindingContextWords = false;
			while((currentToken2 = br2.readLine()) != null) {
				
				// If middle spot is empty, put token there
				if(window[WINDOW_SIZE] == null) {
					
					window[WINDOW_SIZE] = currentToken2;
				
				// Find the next empty spot to put token
				} else {
					for(int i = WINDOW_SIZE + 1; i < window.length; i++) {
						if(window[i] == null) {
							window[i] = currentToken2;

							// If this is the last spot in the array
							if(i == window.length - 1) doneFindingContextWords = true; 

							break;							
						} 
					}
				}
				
				// Update term-context matrix from context words
				if(doneFindingContextWords) {

					term = window[WINDOW_SIZE];
					uniqueWindow = new HashSet<>(Arrays.asList(window));
					for(String contextWord : uniqueWindow) {
						if(contextWord != null && !contextWord.equals(term)) {
							
							try {
								
								// Get index for term and context word
								termIndex = GetIndex(term);
								contextIndex = GetIndex(contextWord);
								
								// Get current value in matrix
								float currentFrequency = C[termIndex][contextIndex];
								
								// Update matrix
								C[termIndex][contextIndex] = currentFrequency + 1;
						
							} catch (Exception e) {
								e.printStackTrace();
							}		
						}
					}
					doneFindingContextWords = false;
					window = shiftLeft(window);
				}
			}
			
			br2.close();
			
			// Write word index to file
			File wordIndexFile = new File("wordIndex.txt");
	        FileOutputStream fos = new FileOutputStream(wordIndexFile);
	        ObjectOutputStream oos = new ObjectOutputStream(fos);
	        oos.writeObject(wordIndex);
	        oos.close();
			
			System.out.println("Term-context matrix updated sucessfully.");

		} catch(IOException e) {
			e.printStackTrace();
		} catch(OutOfMemoryError e) {
			System.out.println("Memory required: " + (Math.pow(matrixSize, 2) * 2) / Math.pow(1000, 2) + " MB");
			e.printStackTrace();
		}
		
		// Write matrix to file, if necessary
		if(writeToFile) {
			
			this.V = matrixSize;
			
			// Store matrix sums
			this.rowSums = new float[this.V];
			this.colSums = new float[this.V];

			System.out.println(C.length);
			
			WriteMatrixToDisk(C);

			// Store matrix sums
			float currentNum = 0f;
			for(int i = 0; i < this.V; i++) {
				for(int j = 0; j < this.V; j++) {
					currentNum = C[i][j];
					if(currentNum != 0) {
						this.rowSums[i] = currentNum;
						this.colSums[j] += currentNum;
						this.matrixSum += currentNum;
					}  
				}
			}
			
			// Release all memory used by frequency term-context matrix
			C = null;
			this.TermContextMatrix_FREQ = null;
			System.gc();
			//this.TermContextMatrix_PPMI = null;
		} else {
			this.TermContextMatrix_FREQ = C;
		}
	}
	
	/**
	 * Filters the corpus if necessary to contain only English words
	 * @param D
	 */
	private void FilterCorpus(File D) {
		
		System.out.println("\nFiltering corpus...");
		
		try {

			// Load Google's dictionary
			Set<String> dict = new HashSet<>();
			
			FileReader dictFr = new FileReader("enable1.txt");
			BufferedReader dictBr = new BufferedReader(dictFr);
			String word;
			while((word = dictBr.readLine()) != null) {
				dict.add(word.trim());
			}
			dictBr.close();
			
			// Filter corpus
			FileWriter fw = new FileWriter("en-context.txt");
			BufferedWriter bw = new BufferedWriter(fw);
			
			FileReader fr = new FileReader(D);
			BufferedReader br = new BufferedReader(fr);
			String line;			
			while((line = br.readLine()) != null) {
				if(dict.contains(line.trim())) bw.write(line + "\n");
			}
			bw.close();
			br.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("Corpus filtered successfully.");
	}
	
	/**
	 * Shifts the entire window array to the left by 1 index
	 * @param window
	 * @return
	 */
	private static String[] shiftLeft(String[] window) {
		String[] newWindow = new String[(WINDOW_SIZE * 2) + 1];
		
		// Loop through old window and populate new
		for(int i = 1; i < window.length; i++) {
			newWindow[i-1] = window[i];
		}

		return newWindow;
	}
	
	/**
	 * Calculates the corresponding PPMI value between two words in the matrix
	 * @param C The context matrix
	 * @param w1 Word 1
	 * @param w2 Word 2
	 * @return
	 */
	public float CalculatePPMIOfTwoWords(int w1, int w2) {
		
		float[][] C = this.TermContextMatrix_FREQ;
		
		try {

			// Calculate matrix sum
			float sum = 0f;
			float sumDelta = 0f;
			for(int i = 0; i < C.length; i++) {
				for(int j = 0; j < C[i].length; j++) {
					sumDelta = C[i][j];
					sum += sumDelta;
				}
			}
					
			float intersection = C[w1][w2] / sum;
			
			float row = 0f;
			float rowDelta = 0f;
			for(int i = 0; i < C.length; i++) {
				rowDelta = C[w1][i];
				row += rowDelta;
			}
			float rowProbability = row / sum;
			
			float col = 0f;
			float colDelta = 0f;
			for(int i = 0; i < C[w1].length; i++) {
				colDelta = C[i][w2];
				col += colDelta;
			}
			float colProbability = (float) (Math.pow(col, 0.75) / Math.pow(sum, 0.75));
			
			float logValue = intersection / (rowProbability * colProbability);
			float maxValue = (float) (Math.log10(logValue) / Math.log10(2));
			
			float PPMI = (float) (Math.max(maxValue, 0));
			
			return PPMI;
			
		} catch(Exception e) {
			e.printStackTrace();
		}
		
		return 0f;
	}
	
	/**
	 * Calculates the PPMI value for an entire term-context matrix and updates each cell accordingly
	 * @param C The matrix
	 */
	public void CalculatePPMI(boolean readFromFile, boolean writeToFile) {
		
		long startTime = System.currentTimeMillis();

		System.out.println("\nConverting matrix to use PPMI values...");
		System.out.println("Memory used: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / Math.pow(1000, 3) + " GB");

		float[][] C = null;
		float[][] newC = null;

		// Variables
		float matrixSum = 0f;
		float matrixDelta = 0f;
		
		float cellValue = 0f;
		float cellProbability = 0f;
		
		float rowSum = 0f;
		float rowDelta = 0f;
		float rowProbability = 0f;
		
		float colSum = 0f;
		float colDelta = 0f;
		float colProbability = 0f;
		
		float logValue = 0f;
		float maxValue = 0f;
		float PPMI = 0f;

		// Read from file if necessary
		if(readFromFile) {
			
			newC = new float[this.V][this.V];
			System.out.println("Memory used: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / Math.pow(1000, 3) + " GB");
					
			try {

				// Read file
				RandomAccessFile file = new RandomAccessFile("context-matrix.raf", "r");
				
				String line;
				String[] values;
				int rowNum = 0;
				int colNum = 0;
				float currentNum = 0f;
				int count = 0;
				while((line = file.readLine()) != null) {
					
					// Check row num
					if(rowNum >= this.V) break;
					
					values = line.split(" ");
					
					for(String value : values) {

						count++;
						
						currentNum = Float.parseFloat(value);

						// Calculate PPMI and store value in newC
						if(currentNum != 0 && currentNum < 16) {
							cellProbability = currentNum / this.matrixSum;
							rowProbability = this.rowSums[rowNum] / this.matrixSum;
							colProbability = this.colSums[colNum] / this.matrixSum;
							
							logValue = cellProbability / (rowProbability * colProbability);
							maxValue = (float) (Math.log10(logValue) / Math.log10(2));
							
							PPMI = Math.max(maxValue, 0);
							
							newC[rowNum][colNum] = PPMI;
						}
						
						if(count % 1000000 == 0)
							System.out.println("Progress: " + (count / 1000000) + " / " + (((long)newC.length * (long)newC[0].length) / 1000000));

						colNum++;
					}
					
					colNum = 0;
					rowNum++;
				}
				
				// Close file
				file.close();
				
				// Store new matrix on instance object
				this.TermContextMatrix_PPMI = newC;
				if(this.TermContextMatrix_PPMI == null) System.out.println("PPMI term-context matrix is null.");

			} catch(IOException e0) {
				e0.printStackTrace();
			}
			
		} else {
			C = this.TermContextMatrix_FREQ;
			newC = new float[C.length][C[0].length];

			// Calculate matrix sum
			for(int i = 0; i < C.length; i++) {
				for(int j = 0; j < C[i].length; j++) {
					matrixDelta = C[i][j];
					matrixSum += matrixDelta;
				}
			}
			
			// For each cell
			for(int u = 0; u < C.length; u++) {
				
				// Calculate row sum
				for(int r = 0; r < C[u].length; r++) {
					rowDelta = C[u][r];
					rowSum += rowDelta;
				}
				rowProbability = rowSum / matrixSum;
				
				for(int v = 0; v < C[u].length; v++) {
					
					cellValue = C[u][v];
					
					if(u != v && cellValue != 0) {

						// Calculate cell probability
						cellProbability = cellValue / matrixSum;
						
						// Calculate col sum
						for(int c = 0; c < C.length; c++) {
							colDelta = C[c][v];
							colSum += colDelta;
						}
						colProbability = (float) (Math.pow(colSum, 0.75) / Math.pow(matrixSum, 0.75));

						logValue = cellProbability / (rowProbability * colProbability);
						maxValue = (float) (Math.log10(logValue) / Math.log10(2));
						
						PPMI = Math.max(maxValue, 0);
						
						newC[u][v] = PPMI;
					}
				}
			}

			this.TermContextMatrix_PPMI = newC;			
		}

		System.out.println("Matrix converted.");

		long endTime = System.currentTimeMillis();
		
		System.out.println("Time taken: " + (endTime - startTime) + " ms");
		
		if(writeToFile) {
			
			try {
				
				WriteMatrixToDisk(this.TermContextMatrix_PPMI);
			
				// Save any settings pertaining to the matrix for later use
				BufferedWriter settingsBw = new BufferedWriter(new FileWriter("settings.txt"));
				settingsBw.write("V " + this.V + "\n");
				settingsBw.close();
				
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}
	
	/**
	 * Write the specified matrix to a random access file
	 * @param C
	 */
	private void WriteMatrixToDisk(float[][] C) {
		
		System.out.println("Writing matrix to disk...");

		try {
		
			BufferedWriter bw = new BufferedWriter(new FileWriter("context-matrix.raf"));
			int row = 0;
			int col = -1;
			for(long i = 0; i < (long)C.length * (long)C[0].length; i++) {
				col++;
				if(col % (C[0].length - 1) == 0 && col != 0) {
					bw.write(String.format("%07.4f %n", C[row][col]));
					row++;
					col = -1;
				} else {
					bw.write(String.format("%07.4f ", C[row][col]));
				}
			}
			bw.close();
			
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("Finished writing matrix to disk.");
	}
	
	/**
	 * Prints the entire matrix for debugging purposes
	 * @param C The matrix
	 */
	public void PrintMatrix() {
		
		float[][] C = this.TermContextMatrix_PPMI;
		
		System.out.println("\nMatrix:");
		for(int i = 0; i < C.length; i++) {
			for(int j = 0; j < C[i].length; j++) {
				if(j >= C[i].length - 1) {
					System.out.printf(" %.1f%n", C[i][j]);
				} else {
					System.out.printf(" %.1f ", C[i][j]);
				}
			}
		}
		System.out.println("End of matrix.");
		
	}
	
	/**
	 * Calculates the cosine similarity between two words in the matrix
	 * This method has a time complexity of O(V).
	 * This method has a space complexity of O(V).
	 * V represents the dimension of the matrix.
	 * @param u Index for word 1.
	 * @param v Index for word 2.
	 * @param fromFile Should we read the matrix from a file?
	 * @param V Length of the matrix.
	 * @return similarity Cosine similarity float.
	 */
	public float CalculateSimilarity(int u, int v, boolean fromFile, int V) {
		
		float[] uRow = new float[V], vRow = new float[V];
		
		System.out.println("u="+u+", v="+v+", V="+V);

		if(fromFile) {
			try {

				RandomAccessFile matrix = new RandomAccessFile("context-matrix.raf", "r");

				// Get uRow
				matrix.seek((long)u * (((long)V * 8) + 2));
				String uRowLine = matrix.readLine();
//				System.out.println(uRowLine);
				String[] uRowString = uRowLine.trim().split(" ");
				for(int i = 0; i < uRowString.length; i++) 
					uRow[i] = Float.parseFloat(uRowString[i]);
				
				// Get vRow
				matrix.seek((long)v * (((long)V * 8) + 2));
				String vRowLine = matrix.readLine();
//				System.out.println(vRowLine);
				String[] vRowString = vRowLine.trim().split(" ");
				for(int i = 0; i < vRowString.length; i++) 
					vRow[i] = Float.parseFloat(vRowString[i]);
				
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
	
			if(this.TermContextMatrix_PPMI == null) 
				System.out.println("PPMI term-context matrix is null.");
			
			float[][] C = this.TermContextMatrix_PPMI;
			if(u != v && u < C.length && v < C.length) {
				uRow = C[u];
				vRow = C[v];
			} else {
				return 0f;
			}
		}

		
		float similarity = 0f;
		
		try {


			float numerator = 0f;
			for(int i = 0; i < V; i++) {
				numerator += uRow[i] * vRow[i];
			}
			
			float uDenom = 0f;
			for(int i = 0; i < uRow.length; i++) {
				uDenom += Math.pow(uRow[i], 2);
			}
			float vDenom = 0f;
			for(int i = 0; i < vRow.length; i++) {
				vDenom += Math.pow(vRow[i], 2);
			}
			
			float denominator = (float) (Math.sqrt(uDenom) * Math.sqrt(vDenom));
			
			similarity = numerator / denominator;

		
		} catch(Exception e) {
			e.printStackTrace();
		}
				
		return similarity;
	}
	
	/**
	 * Calculates cosine similarity from two predefined vectors each representing a single token.
	 * @param Vector for word 1.
	 * @param Vector for word 2.
	 * @return Cosine similarity float.
	 */
	public float CalculateSimilarityByVector(float[] v1, float[] v2) {
		
		float similarity = 0f;
		
		float numerator = 0f;
		for(int i = 0; i < v1.length; i++) {
			numerator += v1[i] * v2[i];
		}
		
		float v1Denom = 0f;
		for(int i = 0; i < v1.length; i++) {
			v1Denom += Math.pow(v1[i], 2);
		}
		float v2Denom = 0f;
		for(int i = 0; i < v2.length; i++) {
			v2Denom += Math.pow(v2[i], 2);
		}
		
		float denominator = (float) (Math.sqrt(v1Denom) * Math.sqrt(v2Denom));
		
		similarity = numerator / denominator;
		
		return similarity;
	}
	
	/**
	 * Gets the top K context words for a given term ranked by cosine similarity
	 * This method has a time complexity of O(V + K).
	 * This method has a space complexity of O(V^2 + K).
	 * V represents the dimension of the matrix and K represents the number of context words to return
	 * @param u Index for the word to get context words for.
	 * @param k Number of context words to return.
	 * @param V Size of the matrix.
	 */
	public String[] GetContext(int u, int k, int V) {
		
		// Check for invalid indexes
		if(u >= V) {
			return null;
		}
		
		// Load word index if null
		if(this.wordIndex == null) {
			LoadWordIndexFromFile();
		}

		float[][] C = this.TermContextMatrix_PPMI;
		float[] uRow = new float[V];
		
		if(C == null) {
			
			try {
				
				RandomAccessFile matrix = new RandomAccessFile("context-matrix.raf", "r");
				matrix.seek((long)u * (((long)V * 8) + 2));
				String uRowString = matrix.readLine();
				String[] uRowValues = uRowString.trim().split(" ");
				for(int i = 0; i < uRowValues.length; i++)
					uRow[i] = Float.parseFloat(uRowValues[i]);
				
			} catch(IOException e) {
				e.printStackTrace();
			}			
		}		
		
		Queue<ContextWordEntry> entriesQueue = new PriorityQueue<>();
		
		// Loop through all context words for index u and calculate similarity
		float similarity = 0f;
		for(int i = 0; i < uRow.length; i++) {
			if(i != u) {
				similarity = CalculateSimilarity(u, i, true, this.V);
				entriesQueue.offer(new ContextWordEntry(i, similarity));
			}
		}
		
		// Get top K context word entries
		ContextWordEntry[] entries = new ContextWordEntry[k];
		for(int i = 0; i < k; i++) {
			entries[i] = entriesQueue.poll();
		}
		
		// Loop through each entry and get its corresponding token
		Set<Map.Entry<String, Integer>> entrySet = wordIndex.entrySet();
		Iterator<Map.Entry<String, Integer>> it;
		
		String[] contextWords = new String[k];
		ContextWordEntry entry;
		Map.Entry<String, Integer> c;
		for(int i = 0; i < entries.length; i++) {
			if(entries[i] != null) {
				entry = entries[i];
				it = entrySet.iterator();
				
				while(it.hasNext()) {
					c = it.next();
					if(c.getValue() == entry.getIndex()) {
						contextWords[i] = c.getKey();
						break;
					}					
				}
			}
		}
		
		return contextWords;
	}
	
	private void LoadWordIndexFromFile() {
		
		File indexFile = new File("wordIndex.txt");

        try {

        	// Read word index from file
		    FileInputStream fis = new FileInputStream(indexFile);
		    ObjectInputStream ois = new ObjectInputStream(fis);
		    this.wordIndex = (HashMap<String, Integer>) ois.readObject();
		    ois.close();

        } catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Gets the corresponding index in the matrix for the specified word
	 * @param Word used to search the word index.
	 */
	public int GetIndex(String word) {
		String s = word.toLowerCase();
		
		if(this.wordIndex == null) {
			LoadWordIndexFromFile();
		}
    	return this.wordIndex.get(s);
	}
	
	/**
	 * Gets the dimension V for the matrix.
	 */
	public int GetSize() {
		if(this.V == 0) {
			
			// Read from settings file
			try {
				BufferedReader br = new BufferedReader(new FileReader("settings.txt"));
				String line;
				if((line = br.readLine()) != null) {
					this.V = Integer.parseInt(line.split(" ")[1]);
				}
			} catch(IOException e) {
				e.printStackTrace();
			}
			
		}
		return this.V;
	}
	
	/**
	 * Gets the context matrix.
	 */
	public float[][] GetTermContextMatrix() {
		return this.TermContextMatrix_PPMI;
	}
	
}