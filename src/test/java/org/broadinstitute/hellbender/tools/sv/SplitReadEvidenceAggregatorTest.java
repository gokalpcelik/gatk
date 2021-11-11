package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class SplitReadEvidenceAggregatorTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SAMSequenceDictionaryExtractor.extractDictionary(new File(GATKBaseTest.FULL_HG19_DICT).toPath());
    private static final String TEST_EVIDENCE = toolsTestDir + "/walkers/sv/pesr/NA12878.sr.txt.gz";

    @Test
    public void testGetEvidenceQueryInterval() {
        final FeatureDataSource<SplitReadEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final SplitReadEvidenceAggregator startAggregator = new SplitReadEvidenceAggregator(source, DICTIONARY, 0, true);

        final List<SplitReadEvidence> test1 = startAggregator.collectEvidence(SVTestUtils.newCallRecordWithCoordinates("test1", "20", 10000325, "20", 20000000));
        Assert.assertEquals(test1.size(), 1);

    }

}