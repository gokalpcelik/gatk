package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary="Combine output files from GetNormalArtifactData in the order defined by a sequence dictionary",
        oneLineSummary = "Combine output files from GetNormalArtifactData in the order defined by a sequence dictionary",
        programGroup = CoverageAnalysisProgramGroup.class
)
public class GatherNormalArtifactData extends CommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "an output of GetNormalArtifacatData")
    final List<File> input = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "output")
    final File output = null;

    @Override
    protected Object doWork() {
        final List<NormalArtifactRecord> data = input.stream()
                .flatMap(file -> NormalArtifactRecord.readFromFile(file).stream()).collect(Collectors.toList());
        NormalArtifactRecord.writeToFile(data, output);
        return "SUCCESS";
    }
}
