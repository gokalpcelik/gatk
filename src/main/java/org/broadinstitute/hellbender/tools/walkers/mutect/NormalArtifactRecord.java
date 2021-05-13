package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

public class NormalArtifactRecord {
    private int normalAltCount;
    private int normalRefCount;
    private int tumorAltCount;
    private int tumorRefCount;
    private double downsampling;
    private String type;

    public int getNormalAltCount() {
        return normalAltCount;
    }

    public int getNormalRefCount() {
        return normalRefCount;
    }

    public int getTumorAltCount() {
        return tumorAltCount;
    }

    public int getTumorRefCount() {
        return tumorRefCount;
    }

    public double getDownsampling() {
        return downsampling;
    }

    public String getType() {
        return type;
    }

    public NormalArtifactRecord(final int normalAltCount, final int normalRefCount, final int tumorAltCount,
                                final int tumorRefCount, final double downsampling, final String type) {
        this.normalAltCount = normalAltCount;
        this.normalRefCount = normalRefCount;
        this.tumorAltCount = tumorAltCount;
        this.tumorRefCount = tumorRefCount;
        this.downsampling = downsampling;
        this.type = type;
    }


    //----- The following two public static methods read and write contamination files
    public static void writeToFile(final List<NormalArtifactRecord> records, final File outputTable) {
        try ( NormalArtifactWriter writer = new NormalArtifactWriter(IOUtils.fileToPath(outputTable)) ) {
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    public static List<NormalArtifactRecord> readFromFile(final File tableFile) {
        try( NormalArtifactReader reader = new NormalArtifactReader(IOUtils.fileToPath(tableFile)) ) {
            return reader.toList();
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }

    //-------- The following methods are boilerplate for reading and writing contamination tables
    public static class NormalArtifactWriter extends TableWriter<NormalArtifactRecord> {
        public NormalArtifactWriter(final Path output) throws IOException {
            super(output, NormalArtifactColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final NormalArtifactRecord record, final DataLine dataLine) {
            dataLine.set(NormalArtifactColumn.NORMAL_ALT_COUNT.toString(), record.getNormalAltCount())
                    .set(NormalArtifactColumn.NORMAL_REF_COUNT.toString(), record.getNormalRefCount())
                    .set(NormalArtifactColumn.TUMOR_ALT_COUNT.toString(), record.getTumorAltCount())
                    .set(NormalArtifactColumn.TUMOR_REF_COUNT.toString(), record.getTumorRefCount())
                    .set(NormalArtifactColumn.DOWNSAMPLING.toString(), record.getDownsampling())
                    .set(NormalArtifactColumn.TYPE.toString(), record.getType());
        }
    }

    private static class NormalArtifactReader extends TableReader<NormalArtifactRecord> {
        public NormalArtifactReader(final Path path) throws IOException {
            super(path);
        }

        @Override
        protected NormalArtifactRecord createRecord(final DataLine dataLine) {
            final int normalAltCount = dataLine.getInt(NormalArtifactColumn.NORMAL_ALT_COUNT);
            final int normalRefCount = dataLine.getInt(NormalArtifactColumn.NORMAL_REF_COUNT);
            final int tumorAltCount = dataLine.getInt(NormalArtifactColumn.TUMOR_ALT_COUNT);
            final int tumorRefCount = dataLine.getInt(NormalArtifactColumn.TUMOR_REF_COUNT);
            final double downsampling = dataLine.getDouble(NormalArtifactColumn.DOWNSAMPLING);
            final String type = dataLine.get(NormalArtifactColumn.TYPE);
            return new NormalArtifactRecord(normalAltCount, normalRefCount, tumorAltCount, tumorRefCount, downsampling, type);
        }
    }

    private enum NormalArtifactColumn {
        NORMAL_ALT_COUNT("normal_alt"),
        NORMAL_REF_COUNT("normal_ref"),
        TUMOR_ALT_COUNT("tumor_alt"),
        TUMOR_REF_COUNT("tumor_ref"),
        DOWNSAMPLING("downsampling"),
        TYPE("type");

        private final String columnName;

        NormalArtifactColumn(final String columnName) {
            this.columnName = Utils.nonNull(columnName);
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection(NORMAL_ALT_COUNT, NORMAL_REF_COUNT,
                TUMOR_ALT_COUNT, TUMOR_REF_COUNT, DOWNSAMPLING, TYPE);
    }
}
