package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.AssemblyComplexity;
import org.broadinstitute.hellbender.tools.walkers.annotator.FeaturizedReadSets;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReferenceBases;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.eclipse.jetty.util.ArrayQueue;

import java.io.File;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Queue;
import java.util.stream.IntStream;

public class Mutect3DatasetEngine {

    private enum VariantType {
        SNV, INSERTION, DELETION
    }

    // number of features for each vectorized read
    private static final int FEATURES_PER_READ = FeaturizedReadSets.FEATURES_PER_READ;

    // number of additional variant features for assembly complexity, reference context
    private static final int NUM_EXTRA_FEATURES = 9;

    // threshold of negative log-10 population allele frequency to consider something an artifact for the purposes of training data
    // we want to be really sure we don't get germline variants
    // TODO: is this really necessary?
    private static final double RARE_POPAF_THRESHOLD = 5.9;

    // very cautious threshold of negative log-10 population allele frequency to consider something germline for training data.
    // There are so many germline variants we can be wasteful!
    private static final double COMMON_POPAF_THRESHOLD = 1;

    // TODO: is this necessary?
    private static final int MIN_REF = 5;

    private final File datasetFile;

    private final int maxRefCount;

    // number of nonartifact data to keep for each artifact datum
    private final int nonArtifactPerArtifact;

    private final FeaturizedReadSets readFeaturizer;

    // simple method to balance data: for each k-alt-read artifact there are
    // nonArtifactPerArtifact (downsampled) k-alt-read non-artifacts.
    private final EnumMap<VariantType, ArrayQueue<Integer>> unmatchedCounts;









    public Mutect3DatasetEngine(final File datasetFile, final int maxRefCount, final int nonArtifactPerArtifact) {
        this.datasetFile = datasetFile;
        this.maxRefCount = maxRefCount;
        this.nonArtifactPerArtifact = nonArtifactPerArtifact;
        readFeaturizer = new FeaturizedReadSets(maxRefCount);

        unmatchedCounts = new EnumMap<VariantType, ArrayQueue<Integer>>(VariantType.class);
        for (final VariantType type : VariantType.values()) {
            unmatchedCounts.put(type, new ArrayQueue<>());
        }
    }

    // add one datum per alt allele
    public void addData(final ReferenceContext ref, final VariantContext vc, final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                         final AlleleLikelihoods<Fragment, Haplotype> logFragmentLikelihoods) {
        final String refBases = ReferenceBases.annotate(ref, vc);
        final String refAllele = vc.getReference().getBaseString();
        final String contig = vc.getContig();
        final int position = vc.getStart();

        // haplotype equivalence counts, haplotype complexity, haplotype dominance
        final Triple<int[], int[], double[]> assemblyComplexity = AssemblyComplexity.annotate(vc, logFragmentLikelihoods);

        //TODO: need this for tumor and normal sample
        final List<List<List<Integer>>> readVectorsByAllele =  readFeaturizer.getReadVectors(vc, sample, likelihoods, logFragmentLikelihoods);



        for (int n = 0; n < vc.getNAlleles() - 1; n++) {


            final Allele altAllele = vc.getAlternateAllele(n);
            final String altBases = altAllele.getBaseString();
            final int diff = altAllele.length() - refAllele.length();
            final VariantType type = diff == 0 ? VariantType.SNV : ( diff > 0 ? VariantType.INSERTION : VariantType.DELETION);



            // TODO: write contig position ref alt
            // TODO: write ref bases

            final List<Double> variantFeatureVector = variantFeatures(n, assemblyComplexity, refBases);
            // TODO: write variant vector

            final List<List<Integer>> refReadVectors = readVectorsByAllele.get(0);
            final List<List<Integer>> altReadVectors = readVectorsByAllele.get(n + 1);
            // TODO: write read vectors




            }

        

    }

    private List<Double> variantFeatures(final int altAlleleIndex, Triple<int[], int[], double[]> assemblyComplexity, final String refBases) {
        final int[] haplotypeEquivalenceCounts = assemblyComplexity.getLeft();
        final int haplotypeComplexity = assemblyComplexity.getMiddle()[altAlleleIndex];
        final double haplotypeDominance = assemblyComplexity.getRight()[altAlleleIndex];

        final List<Double> result = new ArrayList<>(NUM_EXTRA_FEATURES);

        // take integer haplotype equivalence counts (already in order from greatest to least from Mutect)
        // and calculate the fractional share of the 2nd and 3rd, or 0 if none exist
        final double total = MathUtils.sum(haplotypeEquivalenceCounts);
        result.add(haplotypeEquivalenceCounts.length < 2 ? 0.0 : haplotypeEquivalenceCounts[1] / total);
        result.add(haplotypeEquivalenceCounts.length < 3 ? 0.0 : haplotypeEquivalenceCounts[2] / total);
        result.add((double) haplotypeComplexity);
        result.add(haplotypeDominance);

        IntStream.range(1, 6).forEach(repeatLength -> result.add((double) countRepeats(refBases.getBytes(), repeatLength)));

        Utils.validate(result.size() == NUM_EXTRA_FEATURES, "produced a variant feature vector of wrong size");
        return result;
    }

    // count how many repeats of length k surround the middle base
    // example: countRepeats(GACTACTACTG,3) = 3
    private int countRepeats(final byte[] refBases, final int k) {
        final int N = refBases.length;
        final int n = (N - 1) / 2;
        Utils.validateArg(k <= n, "Too few ref bases for given repeat length");

        // extend a repeat forward, to front and then to back(both exclusive)
        // note that this only extends backward if they match the bases going forward
        // that is AAAAAGTGTCC(first G in the middle) will get extended forward through
        // both GTs but won 't be extended back
        int front = n + k;
        while (front < N && refBases[front] == refBases[front - k]) {
            front++;
        }
        int back = n - 1;
        while (back >= 0 && refBases[back] == refBases[back + k]){
            back--;
        }
        final int forwardRepeats = (front - back - 1) / k;

        // same idea but extending backwards first (now back is exclusive and front is inclusive)
        back = n - k;
        while (back >= 0 && refBases[back] == refBases[back + k]) {
            back--;
        }
        front = n + 1;
        while (front < N && refBases[front] == refBases[front - k]) {
            front++;
        }
        final int backwardRepeats = (front - back - 1) / k;

        return FastMath.max(forwardRepeats, backwardRepeats);
    }



    public void shutdown() {
        // close the writer
    }
}
