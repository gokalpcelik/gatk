package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination;
import org.broadinstitute.hellbender.tools.walkers.contamination.PileupSummary;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * <h3>Usage example</h3>
 *
 * <pre>
 * gatk GetNormalArtifactData \
 *   -I tumor.bam \
 *   -I normal.bam \
 *   -normal normal_sample \
 *   -L intervals.list \
 *   -O normal-artifact.table
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Collects data for training normal artifact filter",
        oneLineSummary = "Collects data for training normal artifact filter",
        programGroup = CoverageAnalysisProgramGroup.class)
@DocumentedFeature
public class GetNormalArtifactData extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table", optional=false)
    private File outputTable;

    private List<NormalArtifactRecord> data = new ArrayList<>();

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return false;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public boolean requiresFeatures() {
        return false;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public void onTraversalStart() {
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {


        data.add
    }

    @Override
    public Object onTraversalSuccess() {
        NormalArtifactRecord.writeToFile(data, outputTable);
        return "SUCCESS";
    }
}
