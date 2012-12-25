package tester;

import weka.core.SparseInstance;

public class Subject {
	boolean adopted;
	boolean treated;
	SparseInstance instance;
	double propensity;
	public Subject(boolean adopted, boolean treated, SparseInstance instance) {
		this.adopted = adopted;
		this.treated = treated;
		this.instance = instance;
	}
}
