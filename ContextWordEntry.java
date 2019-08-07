package scyther;
/**
 * A context word for a predetermined term
 * @param index The current index of the context word in the word index
 * @param similarity The similarity between the context word and its relating term
 * @author Jacob Burton
 *
 */
class ContextWordEntry implements Comparable<ContextWordEntry> {

	private int index;
	private float similarity;
	
	public ContextWordEntry(int index, float similarity) {
		this.index = index;
		this.similarity = similarity;
	}
	
	@Override
	public int compareTo(ContextWordEntry otherEntry) {
		return -Float.compare(this.similarity, otherEntry.similarity);
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}

	public float getSimilarity() {
		return similarity;
	}

	public void setSimilarity(float similarity) {
		this.similarity = similarity;
	}
	
}