package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
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
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;
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

    private final int maxRefCount;

    // TODO: is this necessary?
    private static final int MIN_REF = 5;

    private final PrintWriter printWriter;

    // number of nonartifact data to keep for each artifact datum
    private final int nonArtifactPerArtifact;

    private final Set<String> normalSamples;

    // simple method to balance data: for each k-alt-read artifact there are
    // nonArtifactPerArtifact (downsampled) k-alt-read non-artifacts.
    private final EnumMap<VariantType, ArrayQueue<Integer>> unmatchedCounts;


    public Mutect3DatasetEngine(final File datasetFile, final int maxRefCount, final int nonArtifactPerArtifact, final Set<String> normalSamples) {
        try {
            printWriter = new PrintWriter(new FileWriter(Utils.nonNull(datasetFile)));
        } catch (IOException ex) {
            throw new UserException.BadInput("Could not create dataset file writer");
        }

        this.normalSamples = Utils.nonNull(normalSamples);
        this.nonArtifactPerArtifact = nonArtifactPerArtifact;
        this.maxRefCount = maxRefCount;

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

        final List<String> tumorSamples = likelihoods.samples().stream().filter(sample -> !normalSamples.contains(sample)).collect(Collectors.toList());
        final List<List<List<Integer>>> normalReadVectorsByAllele =  FeaturizedReadSets.getReadVectors(vc, normalSamples, likelihoods, logFragmentLikelihoods, maxRefCount);
        final List<List<List<Integer>>> tumorReadVectorsByAllele =  FeaturizedReadSets.getReadVectors(vc, tumorSamples, likelihoods, logFragmentLikelihoods, maxRefCount);

        // ref reads have already been downsampled by the read featurizer
        final List<List<Integer>> tumorRefReads = tumorReadVectorsByAllele.get(0);
        final List<List<Integer>> normalRefReads = normalReadVectorsByAllele.get(0);

        for (int n = 0; n < vc.getNAlleles() - 1; n++) {
            final String altAllele = vc.getAlternateAllele(n).getBaseString();
            final int diff = altAllele.length() - refAllele.length();
            final VariantType type = diff == 0 ? VariantType.SNV : ( diff > 0 ? VariantType.INSERTION : VariantType.DELETION);


            // format is (all vectors are one vector per line, single-spaced)
            // CONTIG:POSITION,REF->ALT
            // REFERENCE CONTEXT BASES
            // variant feature vector
            // tumor_ref_count tumor_alt_count normal_ref_count normal_alt_count
            // tumor ref reads, one per line
            // tumor alt reads
            // normal ref reads
            // normal alt reads
            printWriter.printf("%s:%d,%s->%s", contig, position, refAllele, altAllele);
            printWriter.print(refBases);

            // print variant feature vector with 3 decimal places, single-spaced
            final List<Double> variantFeatureVector = variantFeatures(n, assemblyComplexity, refBases);
            printWriter.print(numberString(variantFeatureVector, "%.3f", " "));

            final List<List<Integer>> tumorAltReads = tumorReadVectorsByAllele.get(n+1);
            final List<List<Integer>> normalAltReads = normalReadVectorsByAllele.get(n+1);
            printWriter.printf("%d %d %d %d", tumorRefReads.size(), tumorAltReads.size(), normalRefReads.size(), normalAltReads.size());
            tumorRefReads.forEach(r -> printWriter.print(numberString(r)));
            tumorAltReads.forEach(r -> printWriter.print(numberString(r)));
            normalRefReads.forEach(r -> printWriter.print(numberString(r)));
            normalAltReads.forEach(r -> printWriter.print(numberString(r)));
            }
    }

    private String numberString(final List<? extends Number> numbers) {
        return numberString(numbers, "%.3f", " ");
    }

    private String numberString(final List<? extends Number> numbers, final String formatString, final String separator) {
        return numbers.stream().map(x -> String.format(formatString, x)).collect(Collectors.joining(separator));
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
